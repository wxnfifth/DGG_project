#include "stdafx.h"
#include <windows.h>
#include "Shlwapi.h"
#include "wxnTime.h"
#include "svg_precompute.h"
#include "svg_precompute_debug.h"
#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"

#include "MMP\geodesic_algorithm_vg_mmp.h"
#include "JIAJUN\dgg_pruning.h"
#include "svg_definition.h"
#include "wxn\wxnMath.h"
#include "wxn\wxn_path_helper.h"
#include "wxn\wxn_dijstra.h"
#include "YXMetric\YXPathTracer.h"
#include "wxn\wxn_path_helper.h"
#include <thread>
#include <random>


template <class T>
void dijkstra_pruning_disk_graph_debug(const vector<vector<int>>&  graph_neighbor,
	const vector<vector<T>>& graph_neighbor_dis,
	vector<bool>& current_graph_neighbor_deleted,
	int src, double eps_vg, vector<T>& dis, vector<bool>& mark)
{
	const double max_error = 1e-5;
	struct QueueNode{
		T dis;
		int node_index;
		QueueNode(){}
		QueueNode(int _node_index, double _dis) {
			dis = _dis;
			node_index = _node_index;
		}
		bool operator<(const QueueNode& other)const{
			return dis > other.dis;
		}
	};

	fill(dis.begin(), dis.end(), JiajunMaxDist);
	fill(mark.begin(), mark.end(), false);

	priority_queue<QueueNode> que;
	dis[src] = 0;
	mark[src] = true;
	T max_dis_radius = 0;
	map<int, int> node_map;
	for (int i = 0; i < graph_neighbor[src].size(); ++i) {
		int v = graph_neighbor[src][i];
		T d = graph_neighbor_dis[src][i];
		max_dis_radius = max(max_dis_radius, d);
		dis[v] = d;
		que.push(QueueNode(v, d));
		node_map[v] = i;
	}
	int cnt = 0;
	while (!que.empty()) {
		QueueNode u = que.top();
		cnt++;
		que.pop();
		if (mark[u.node_index]) continue;
		mark[u.node_index] = true;
		if (u.dis > max_dis_radius) {
			break;
		}
		bool found_flag = false;
		for (int i = 0; i < graph_neighbor[u.node_index].size(); ++i) {
			int v = graph_neighbor[u.node_index][i];
			T d = graph_neighbor_dis[u.node_index][i];
			if (u.node_index == v) {
				continue;
			}
			if (u.dis + d < dis[v] * (1 + eps_vg)) {
				if (src == 0) {
					//printf("v %d\n", v);
				}
				QueueNode b;
				b.node_index = v;
				b.dis = min(u.dis + d, dis[v]);
				dis[v] = b.dis;
				current_graph_neighbor_deleted[node_map[v]] = true;
				que.push(b);
			}
		}
	}
	if (src == 0) {
		printf("cnt %d\n", cnt);
	}
	//printf("cnt %d\n", cnt);
	if (false) {
		for (int v : graph_neighbor[src]) {
			dis[v] = JiajunMaxDist;
			mark[v] = false;
		}
		dis[src] = JiajunMaxDist;
		mark[src] = false;
	}
}


template<class T>
void dijkstraPruningThreadDebug(int thread_id, int thread_num, int node_number,
	const vector<vector<int>>& graph_neighbor,
	const vector<vector<T>>& graph_neighbor_dis,
	vector<vector<bool>>& graph_neighbor_deleted, double eps_vg, CRichModel& model)
{
	vector<T> dis(node_number);
	vector<bool> mark(node_number);
	fill(dis.begin(), dis.end(), JiajunMaxDist);
	fill(mark.begin(), mark.end(), false);
	int part_size = node_number / thread_num;
	int begin = thread_id * part_size;
	int end;
	if (thread_id != thread_num - 1) {
		end = (thread_id + 1) * part_size - 1;
	}
	else {
		end = node_number - 1;
	}
	ElapasedTime time_once;
	double past_time = 0;
	for (int i = begin; i <= end; ++i) {
		if (time_once.getTime() - past_time > 5) {
			past_time = time_once.getTime();
			char buf[128];
			double percent = (double)(i - begin)  * 100. / double(end - begin + 1);
			double current_time = time_once.getTime();
			double remain_time = current_time / percent * (100 - percent);
			printf("Computed %.0lf percent, time %lf, estimate_remain_time %lf\n",
				percent, current_time, remain_time);
		}

		dijkstra_pruning_induced_graph_debug<T>(graph_neighbor,
			graph_neighbor_dis, graph_neighbor_deleted[i], i,
			eps_vg, dis, mark);
		if (i == 0) {//debug here!!
			printf("start debug!\n");
			vector<int> dests;
			int j = 0;
			for (auto flag_deleted : graph_neighbor_deleted[i]) {
				if (!flag_deleted) {
					dests.push_back(graph_neighbor[i][j]);
				}
				++j;
			}

			CylinderPath cylinder_path(0.01);
			cylinder_path.cntGeodesicPaths(model, i, dests);
			printf("after pruning neigh sz %d\n", dests.size());

		}
	}
}

template<class T>
void readInputFile(const string& svg_file_name,
	vector<vector<int>>& graph_neighbor,
	vector<vector<T>>& graph_neighbor_dis,
	vector<vector<bool>>& graph_neighbor_deleted,
	int& node_number)
{

	std::ifstream input_file(svg_file_name, std::ios::in | std::ios::binary);
	HeadOfSVG head_of_svg;
	input_file.read((char*)&head_of_svg, sizeof(head_of_svg));
	head_of_svg.print();
	node_number = head_of_svg.num_of_vertex;
	graph_neighbor.reserve(node_number);
	graph_neighbor.resize(node_number);
	graph_neighbor_dis.reserve(node_number);
	graph_neighbor_dis.resize(node_number);
	graph_neighbor_deleted.reserve(node_number);
	graph_neighbor_deleted.resize(node_number);

	for (int i = 0; i < head_of_svg.num_of_vertex; ++i) {
		BodyHeadOfSVG body_head;
		input_file.read((char*)&body_head, sizeof(body_head));
		std::vector<BodyPartOfSVGWithAngle> body_parts;
		for (int j = 0; j < body_head.neighbor_num; ++j) {
			BodyPartOfSVGWithAngle body_part;
			input_file.read((char*)&body_part, sizeof(body_part));
			body_parts.push_back(body_part);
		}
		int u = body_head.source_index;
		int number_of_neighbor = body_parts.size();
		graph_neighbor[u].reserve(number_of_neighbor);
		graph_neighbor_dis[u].reserve(number_of_neighbor);
		graph_neighbor_deleted[u].reserve(number_of_neighbor);
		for (auto body : body_parts) {
			graph_neighbor[u].push_back(body.dest_index);
			graph_neighbor_dis[u].push_back(body.dest_dis);
			graph_neighbor_deleted[u].push_back(false);
		}
		if (i > 0 && i % (head_of_svg.num_of_vertex / 10) == 0){
			std::cerr << "read " << i * 100 / head_of_svg.num_of_vertex << " percent \n";
		}
	}
	input_file.close();
}

template<class T>
void wxn_pruning_debug(const string& svg_file_name, double eps_vg, string& test_output_filename, int thread_num, double& prune_time, CRichModel& model)
{
	printf("debug!\n");
	vector<vector<int>> graph_neighbor;
	vector<vector<T>> graph_neighbor_dis;
	vector<vector<bool>> graph_neighbor_deleted;
	int node_number;

	readInputFile(svg_file_name, graph_neighbor, graph_neighbor_dis, graph_neighbor_deleted, node_number);

	ElapasedTime t;
	{
		int part_size = node_number / thread_num;
		vector<std::thread> tt(thread_num);
		for (int thread_id = 0; thread_id < thread_num; ++thread_id) {
			tt[thread_id] = std::thread(&dijkstraPruningThreadDebug<T>, thread_id,
				thread_num, node_number, ref(graph_neighbor),
				ref(graph_neighbor_dis), ref(graph_neighbor_deleted),
				eps_vg, ref(model));

		}
		for (int i = 0; i < thread_num; ++i) {
			tt[i].join();
		}

	}
	t.printTime();
	prune_time = t.getTime();
	std::ifstream input_file(svg_file_name, std::ios::in | std::ios::binary);
	HeadOfSVG head_of_svg;
	input_file.read((char*)&head_of_svg, sizeof(head_of_svg));
	head_of_svg.print();

	ofstream output_file(test_output_filename, std::ios::out | std::ios::binary);
	output_file.write((char*)&head_of_svg, sizeof(head_of_svg));

	int cnt = 0;
	for (int i = 0; i < head_of_svg.num_of_vertex; ++i) {
		BodyHeadOfSVG body_head;
		input_file.read((char*)&body_head, sizeof(body_head));
		vector<int> new_graph_neighbor;
		vector<T> new_graph_neighbor_dis;
		vector<int> origin2current(body_head.neighbor_num, -1);
		for (int j = 0; j < body_head.neighbor_num; ++j) {
			if (!graph_neighbor_deleted[i][j]) {
				cnt++;
				origin2current[j] = new_graph_neighbor.size();
				new_graph_neighbor.push_back(graph_neighbor[i][j]);
				new_graph_neighbor_dis.push_back(graph_neighbor_dis[i][j]);
			}
		}
		int start_pos = 0;
		for (int j = 0; j < origin2current.size(); ++j) {
			if (origin2current[j] == -1) {
				origin2current[j] = start_pos;
			}
			else{
				start_pos = origin2current[j];
			}
		}
		vector<BodyPartOfSVGWithAngle> body_parts;
		for (int j = 0; j < body_head.neighbor_num; ++j){
			BodyPartOfSVGWithAngle b;
			input_file.read((char*)&b, sizeof(b));
			body_parts.push_back(b);
		}
		body_head.neighbor_num = new_graph_neighbor.size();
		output_file.write((char*)&body_head, sizeof(body_head));
		for (int j = 0; j < body_parts.size(); ++j) {
			if (!graph_neighbor_deleted[i][j]) {
				BodyPartOfSVGWithAngle b = body_parts[j];
				if (b.begin_pos == -1) {
					b.begin_pos = -1;
					b.end_pos = -1;
				}
				else {
					b.begin_pos = origin2current[b.begin_pos];
					b.end_pos = origin2current[b.end_pos];
				}
				output_file.write((char*)&b, sizeof(b));
			}
		}
	}
	printf("average neigh %lf\n", (double)cnt / node_number);
	input_file.close();
	output_file.close();



}


template <class T>
void dijkstra_pruning_induced_graph_debug(const vector<vector<int>>&  graph_neighbor,
	const vector<vector<T>>& graph_neighbor_dis,
	vector<bool>& current_graph_neighbor_deleted,
	int src, double eps_vg, vector<T>& dis, vector<bool>& mark)
{
	//int node_number = graph_neighbor.size();
	//vector<T> dis(node_number);
	//vector<bool> mark(node_number);
	//fill(dis.begin(), dis.end(), JiajunMaxDist);
	//fill(mark.begin(), mark.end(), false);
	const double max_error = 1e-5;
	struct QueueNode{
		T dis;
		int node_index;
		QueueNode(){}
		QueueNode(int _node_index, double _dis) {
			dis = _dis;
			node_index = _node_index;
		}
		bool operator<(const QueueNode& other)const{
			return dis > other.dis;
		}
	};
	priority_queue<QueueNode> que;
	dis[src] = 0;
	mark[src] = true;
	//printf("line 1834\n");
	map<int, int> node_map;
	//printf("sz %d\n", graph_neighbor[src].size());
	for (int i = 0; i < graph_neighbor[src].size(); ++i) {
		int v = graph_neighbor[src][i];
		//printf("src %d i %d v %d ", src, i, v);
		T d = graph_neighbor_dis[src][i];
		//printf("d %lf\n", d);
		dis[v] = d;
		que.push(QueueNode(v, d));
		//printf("line 1844\n");
		node_map[v] = i;
		//printf("line 1845\n");
	}
	//printf("line 1842\n");
	int cnt = 0;
	while (!que.empty()) {
		QueueNode u = que.top();
		cnt++;
		que.pop();
		if (mark[u.node_index]) continue;
		if (src == 0) {
			//if (cnt < 4) {
			//printf("u %d ", u.node_index);
			//}
		}
		mark[u.node_index] = true;
		bool found_flag = false;
		for (int i = 0; i < graph_neighbor[u.node_index].size(); ++i) {
			int v = graph_neighbor[u.node_index][i];
			T d = graph_neighbor_dis[u.node_index][i];
			//if (src == 0 && v == 205) {
			//	printf("u %d v %d d %lf\n", u.node_index, v, d);
			//}
			if (fabs(dis[v] - JiajunMaxDist) < max_error || u.node_index == v) {
				continue;
			}
			if (u.dis + d < dis[v] * (1 + eps_vg)) {
				if (src == 0) {
					//printf("v %d\n", v);
				}
				QueueNode b;
				b.node_index = v;
				b.dis = min(u.dis + d, dis[v]);
				dis[v] = b.dis;
				current_graph_neighbor_deleted[node_map[v]] = true;
				que.push(b);
			}
		}
	}
	if (src == 0) {
		printf("cnt %d\n", cnt);
	}
	//printf("cnt %d\n", cnt);
	for (int v : graph_neighbor[src]) {
		dis[v] = JiajunMaxDist;
		mark[v] = false;
	}
	dis[src] = JiajunMaxDist;
	mark[src] = false;
}


void svg_precompute_ich_multithread_debug(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta, int thread_num)
{
	double ich_multi_time = 0;
	svg_precompute_ich_multithread_before_pruning(input_obj_name, eps_vg, svg_file_name, const_for_theta, thread_num, ich_multi_time);
	string wxn_output_filename(svg_file_name.substr(0, svg_file_name.length() - 7) + "_pruning.binary");
	double prune_time;
	ElapasedTime wxn_pruning_time;
	CRichModel model(input_obj_name);
	model.Preprocess();
	wxn_pruning_debug<double>(svg_file_name, eps_vg, wxn_output_filename, thread_num, prune_time, model);
	wxn_pruning_time.printTime("wxn pruning time");
	fprintf(stderr, "prunning time %lf\n", prune_time);
	fprintf(stderr, "total_time_and_pruning %lf\n", ich_multi_time + prune_time);
}

