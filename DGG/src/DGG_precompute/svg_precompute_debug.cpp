#include "stdafx.h"
#include <windows.h>
#include "Shlwapi.h"
#include "svg_precompute_debug.h"
#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"

#include "MMP\geodesic_algorithm_vg_mmp.h"
#include "JIAJUN\dgg_pruning.h"
#include "svg_definition.h"
#include "wxn\wxnTime.h"
#include "wxn\wxnBuffer.h"
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
		ElapasedTime one_itr_time;
		
		double percent = (double)(i - begin)  * 100. / double(end - begin + 1);
		time_once.printEstimateTime(5, percent);

		dijkstra_pruning_induced_graph_debug<T>(graph_neighbor,
			graph_neighbor_dis, graph_neighbor_deleted[i], i,
			eps_vg, dis, mark);
		if (i == 0) {//debug here!!
			printf("start debug!\n");
			//find if 2170 is in graph_neighbor
			bool flag_found_2170 = false;
			double v0_to_v2170;
			{
				int cnt = 0;
				for (auto n : graph_neighbor[i]) {
					if (n == 2170) {
						flag_found_2170 = true;
					}
					if (n == 2170) {
						printf("v_dest %d, deleted %d, neighbor_dis %.10lf\n", n, graph_neighbor_deleted[i][cnt] ? 1 : 0, graph_neighbor_dis[i][cnt]);
						v0_to_v2170 = graph_neighbor_dis[i][cnt];
					}
					if (n == 2812) {
						printf("v_dest %d, deleted %d, neighbor_dis %.10lf\n", n, graph_neighbor_deleted[i][cnt] ? 1 : 0, graph_neighbor_dis[i][cnt]);
						
					}
					cnt++;
				}
			}
			double v2170_to_v2812;
			{
				int cnt = 0;
				for (auto n : graph_neighbor[2170]) {
					if (n == 2812) {
						printf("2170 to 2812, dis %.10lf\n", graph_neighbor_dis[2170][cnt]);
						v2170_to_v2812 = graph_neighbor_dis[2170][cnt];
						printf("0 to 2170 to 2812 sum %.10lf\n", v0_to_v2170 + v2170_to_v2812);
					}
					cnt++;
				}
			}

			if (flag_found_2170) {
				printf("found 2170 true!\n");
			}
			else{
				printf("found 2170 false\n");
			}
			vector<int> dests;
			int j = 0;
			for (auto flag_deleted : graph_neighbor_deleted[i]) {
				if (!flag_deleted) {
					dests.push_back(graph_neighbor[i][j]);
				}
				++j;
			}
			//CylinderPath cylinder_path(0.01);
			//cylinder_path.cntGeodesicPaths(model, i, dests);
			
			printf("after pruning neigh sz %d\n", dests.size());

		}
		if (i <= 10) {
			int cnt = 0;
			for (auto flag_deleted : graph_neighbor_deleted[i]) {
				if (!flag_deleted) {
					cnt++;
				}
			}
			printf("source %d , after pruning neigh sz %d, time %lf\n", i, cnt, time_once.getTime() - past_time);
			past_time = time_once.getTime();
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
	printf("pruning debug!\n");
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
			//printf("u %d ", u.node_index);c
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

void getDisAndAngles_debug(const CRichModel& model, int source, double eps_vg, vector<int>& dests, vector<double>& angles, vector<double>& dis)
{
	CICHWithFurtherPriorityQueue alg(model, vector < int > {source});
	set<int> fixed_dests;
	alg.ExecuteLocally_DGG(eps_vg, fixed_dests);
	//dests.assign(fixed_dests.begin(), fixed_dests.end());
	dests.clear();
	dests.reserve(fixed_dests.size());
	dis.reserve(fixed_dests.size());
	for (auto d : fixed_dests) {
		map<int, CICHWithFurtherPriorityQueue::InfoAtVertex>::const_iterator it = alg.m_InfoAtVertices.find(d);
		int indexOfParent = it->second.indexOfParent;
		int parent = 0;
		if (it->second.fParentIsPseudoSource) {
			parent = indexOfParent;
		}
		else {
			parent = it->second.indexOfRootVertOfParent;
		}
		if (source == 0 && d == 2812) {
			printf("source %d dest %d parent %d dis %.10lf\n", source, d, parent, it->second.disUptodate);
			vector<IntersectionWithPath> resultingPath;
			alg.FindSourceVertex(d, resultingPath);
			CylinderPath path(0.00001);
			vector<CPoint3D> path_points;
			for (auto& p : resultingPath) {
				path_points.push_back(p.GetPosition(model));
			}
			path.addLines(path_points);
			path.write_to_file("path_0_to_2812.obj");
			path.write_path_points_to_file(path_points, resultingPath, "path_0_to_2812_dgg_points.txt");
		}
		if (parent == source) {
			dests.push_back(d);
			dis.push_back(it->second.disUptodate);
		}
	}
	angles.reserve(dests.size());
	for (auto d : dests) {
		bool isVert;
		int id;
		CPoint3D t = alg.BackTraceDirectionOnly(d, isVert, id);
		auto& neighs = model.Neigh(source);
		vector<double> sum_angle(neighs.size() + 1);
		sum_angle[0] = 0;
		for (int j = 1; j <= neighs.size(); ++j) {
			sum_angle[j] = sum_angle[j - 1] + neighs[j - 1].second;
		}
		double angle = 0;
		if (isVert) {
			bool flag_found = false;
			for (int j = 0; j < neighs.size(); ++j) {
				auto& neigh = neighs[j];
				if (id == model.Edge(neigh.first).indexOfRightVert) {
					//printf("yes\n");
					flag_found = true;
					angle = sum_angle[j];
					break;
				}
			}
			if (!flag_found) {
				//printf("not found vert!\n");
				angle = -1;
			}
		}
		else { // is edge
			int v0 = model.Edge(id).indexOfLeftVert;
			int v1 = model.Edge(id).indexOfRightVert;
			bool flag = false;
			//CPoint3D p_cpoint3d(p.x(), p.y(), p.z());
			for (int j = 0; j < neighs.size(); ++j) {
				auto& neigh = neighs[j];
				if (v0 == model.Edge(neigh.first).indexOfRightVert) {
					int jminus1 = (j - 1 + neighs.size()) % neighs.size();
					int vjminus1 = model.Edge(neighs[jminus1].first).indexOfRightVert;
					int jplus1 = (j + 1) % neighs.size();
					int vjplus1 = model.Edge(neighs[jplus1].first).indexOfRightVert;
					//printf("v1 %d j -1 %d j + 1 %d\n" , v1, vjminus1, vjplus1); 

					if (v1 == vjminus1) {//v1 first
						double l = model.Edge(neighs[jminus1].first).length;
						double r = (model.Vert(source) - t).Len();
						double b = (model.Vert(vjminus1) - t).Len();
						angle = sum_angle[jminus1] + acos((l * l + r * r - b * b) / (2 * l * r));
					}
					else if (v1 == vjplus1) {//v0 first
						double l = model.Edge(neighs[j].first).length;
						double r = (model.Vert(source) - t).Len();
						double b = (model.Vert(v0) - t).Len();
						angle = sum_angle[j] + acos((l * l + r * r - b * b) / (2 * l * r));
					}
					else{
						//fprintf(stderr, "error line 1081\n");
					}
					flag = true;
					break;
				}
			}
			if (!flag) {
				//fprintf(stderr, "flag %d\n", flag);
				//printf("source %d\n", source);
				//printf("v0 %d v1 %d\n", v0, v1);
				//printBallToObj(vector < CPoint3D > {model.Vert(source)},"bunny_nf10k_source.obj",0.01);
				//printBallToObj(vector < CPoint3D > {model.Vert(v0),model.Vert(v1)}, "bunny_nf10k_edge_point.obj", 0.01);
				////printf("")
				//	for (int j = 0; j < neighs.size(); ++j) {
				//		auto& neigh = neighs[j];
				//		printf("right %d\n", model.Edge(neigh.first).indexOfRightVert);
				//	}
				angle = -1;
			}
		}
		angles.push_back(angle);
	}
}

void getFanOutput_debug(const vector<int>& dests, const vector<double>& angles,
	const vector<double>& dis, const CRichModel& model,
	double theta, int source,
	vector<BodyPartOfSVGWithAngle>& body_parts_with_angle)
{

	body_parts_with_angle.reserve(dests.size());
	for (int i = 0; i < dests.size(); ++i) {
		if (angles[i] >= 0) {
			BodyPartOfSVGWithAngle b_with_angle(dests[i], dis[i], angles[i], 0, 0);
			body_parts_with_angle.push_back(b_with_angle);
		}
	}
	sort(body_parts_with_angle.begin(), body_parts_with_angle.end());
	float angle_sum = model.AngleSum(source);
	vector<double> tmp_angles(body_parts_with_angle.size() * 2);
	for (int i = 0; i < body_parts_with_angle.size(); ++i) {
		tmp_angles[i] = body_parts_with_angle[i].angle;
	}
	for (int i = body_parts_with_angle.size(); i < tmp_angles.size(); ++i) {
		tmp_angles[i] = body_parts_with_angle[i - body_parts_with_angle.size()].angle + angle_sum;
	}
	for (int i = 0; i < body_parts_with_angle.size(); ++i) {//assume i is father
		float father_angle = body_parts_with_angle[i].angle;
		//based on father_angle as 0
		float start_angle = M_PI - theta + father_angle;
		float end_angle = angle_sum - (M_PI - theta) + father_angle;
		if (start_angle > end_angle) {
			body_parts_with_angle[i].begin_pos = -1;
			body_parts_with_angle[i].end_pos = -1;
			continue;
		}

		int start_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), start_angle) - tmp_angles.begin();
		if (start_pos > 0) start_pos--;
		int end_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), end_angle) - tmp_angles.begin();
		//printf("start_angle %lf start_pos_angle %lf end_angle  %lf end_pos_angle %lf\n" , start_angle, tmp_angles[start_pos], end_angle, tmp_angles[end_pos]);
		if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
		if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
		body_parts_with_angle[i].begin_pos = start_pos;
		body_parts_with_angle[i].end_pos = end_pos;
	}
}

void ichPropogateHead_debug(const HeadOfSVG& head, const string& part_svg_filename, double eps_vg, double theta, const CRichModel& model)
{
	WxnBuffer wxn_buffer;
	wxn_buffer.open(part_svg_filename);
	wxn_buffer.addStruct(&head, sizeof(head));
	ElapasedTime time_once;

	double past_time;
	double average_degree = 0;
	for (int source = head.begin_vertex_index; source <= head.end_vertex_index; ++source) {
		if (time_once.getTime() - past_time > 5) {
			past_time = time_once.getTime();
			char buf[128];
			sprintf(buf, "Computed %.0lf percent", (double)(source - head.begin_vertex_index) * 100. / (head.end_vertex_index - head.begin_vertex_index));
			time_once.printTime(buf);
		}

		vector<int> dests;
		vector<double> angles;
		vector<double> dis;
		getDisAndAngles_debug(model, source, eps_vg, dests, angles, dis);

		vector<BodyPartOfSVGWithAngle> body_parts_with_angle;
		getFanOutput_debug(dests, angles, dis, model, theta, source, body_parts_with_angle);

		BodyHeadOfSVG body_header(source, body_parts_with_angle.size());
		//output_file.write((char*)&body_header, sizeof(body_header));
		wxn_buffer.addStruct((void*)&body_header, sizeof(body_header));
		for (auto& b : body_parts_with_angle) {
			//output_file.write((char*)&b, sizeof(b));
			wxn_buffer.addStruct((void*)&b, sizeof(b));
		}
		average_degree += body_parts_with_angle.size();
	}
	wxn_buffer.close();

}

void combinePartPrecomputeFiles_debug(const vector<string>& part_filenames, const string& svg_filename, int num_of_vertex, int thread_num)
{
	HeadOfSVG head;
	head.begin_vertex_index = 0; head.end_vertex_index = num_of_vertex - 1;
	head.num_of_vertex = num_of_vertex;
	ofstream output_file(svg_filename.c_str(), ios::out | ios::binary);
	output_file.write((char*)&head, sizeof(head));

	for (int thread_id = 0; thread_id < thread_num; ++thread_id) {
		const string& part_svg_filename = part_filenames[thread_id];
		std::ifstream input_file(part_svg_filename, std::ios::in | std::ios::binary);
		HeadOfSVG head_of_svg;
		input_file.read((char*)&head_of_svg, sizeof(head_of_svg));
		//		printf("head %d %d vert_sz %d\n" , head_of_svg.begin_vertex_index , head_of_svg.end_vertex_index, head_of_svg.num_of_vertex);

		for (int i = head_of_svg.begin_vertex_index; i <= head_of_svg.end_vertex_index; ++i) {
			//printf("i%d\n" , i);
			BodyHeadOfSVG body_head;
			input_file.read((char*)&body_head, sizeof(body_head));
			//printf("readed head\n");
			vector<BodyPartOfSVGWithAngle> body_parts;
			//printf("neigh %d\n" , body_head.neighbor_num);
			body_parts.reserve(body_head.neighbor_num);
			//printf("neigh %d\n" , body_head.neighbor_num);
			for (int j = 0; j < body_head.neighbor_num; ++j) {
				//printf("j%d\n" , j);
				BodyPartOfSVGWithAngle body_part;
				input_file.read((char*)&body_part, sizeof(body_part));
				body_parts.push_back(body_part);
			}
			output_file.write((char*)&body_head, sizeof(body_head));
			for (auto& b : body_parts) {
				output_file.write((char*)&b, sizeof(b));
			}
		}
		input_file.close();
	}
	output_file.close();
}


void svg_precompute_ich_multithread_before_pruning_debug(const string& input_obj_name, double eps_vg,
	string& svg_file_name, double const_for_theta,
	int thread_num, double& ich_multi_time)
{
	char buf[1024];
	sprintf(buf, "%s_DGGICH%.10lf_c%.0lf.binary", input_obj_name.substr(0, input_obj_name.length() - 4).c_str(), eps_vg, const_for_theta);
	svg_file_name = string(buf);
	printf("svg filename is %s\n", svg_file_name.c_str());
	int svg_file_exists = PathFileExists(svg_file_name.c_str());
	if (svg_file_exists == 1)
	{
		return;
	}

	ElapasedTime total_t;
	double theta = asin(sqrt(eps_vg));
	theta *= const_for_theta;
	fprintf(stderr, "******** eps %lf const %lf theta %lf du\n", eps_vg, const_for_theta, theta / M_PI * 180.0);
	CRichModel model(input_obj_name);
	model.Preprocess();

	vector<HeadOfSVG> heads;
	vector<string> svg_part_file_names;
	ElapasedTime time_multi;
	int part_size = model.GetNumOfVerts() / thread_num;
	for (int i = 0; i < thread_num; ++i) {
		HeadOfSVG head;
		head.num_of_vertex = model.GetNumOfVerts();
		head.begin_vertex_index = i * part_size;
		if (i != thread_num - 1) {
			head.end_vertex_index = (i + 1)*part_size - 1;
		}
		else{
			head.end_vertex_index = model.GetNumOfVerts() - 1;
		}
		heads.push_back(head);
		printf("head %d %d vert_sz %d\n", head.begin_vertex_index, head.end_vertex_index, head.num_of_vertex);
		char buf[1024];
		sprintf(buf, "%s_DGGICH%.10lf_c%.0lf_part%d.binary", input_obj_name.substr(0, input_obj_name.length() - 4).c_str(), eps_vg, const_for_theta, i);
		svg_part_file_names.push_back((string)buf);
	}

	vector<std::thread> tt(thread_num);
	for (int i = 0; i < thread_num; ++i) {
		tt[i] = std::thread(&ichPropogateHead_debug, heads[i], svg_part_file_names[i], eps_vg, theta, std::ref(model));
	}
	for (int i = 0; i < thread_num; ++i) {
		tt[i].join();
	}
	ich_multi_time = time_multi.getTime();
	time_multi.printTime("ich_multi_time");

	ElapasedTime combine_time;
	combinePartPrecomputeFiles_debug(svg_part_file_names, svg_file_name, model.GetNumOfVerts(), thread_num);
	combine_time.printTime("combine time");

}


void svg_precompute_ich_multithread_debug(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta, int thread_num)
{
	double ich_multi_time = 0;
	svg_precompute_ich_multithread_before_pruning_debug(input_obj_name, eps_vg, svg_file_name, const_for_theta, thread_num, ich_multi_time);
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

