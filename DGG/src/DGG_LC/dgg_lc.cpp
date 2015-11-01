
#include <string>
#include <ctime>
#include <random>
#include <windows.h>
#include <psapi.h>
#include <stdio.h>   
#include <tchar.h>
#include <cmath>
#include <cfloat>
#include "wxn\wxnTime.h"
#include "wxn\wxn_dijstra.h"
#include "svg_definition.h"
#include "ich\ICHWithFurtherPriorityQueue.h"

using namespace std;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void output_error_distribute(double max_dis, const vector<pair<double,double>>& errors_list, double& all_average_error)
{
  all_average_error = 0;
  double longest_dis = max_dis;
  int interval_num = 100;
  vector<double> average_error_sep(interval_num+1, 0.0);
  vector<double> average_error_cnt(interval_num+1, 0);

  for (int i = 0; i < errors_list.size(); ++i) {
    double dis = errors_list[i].first;
    double error = errors_list[i].second;
    double percent =  dis / longest_dis;
    int pos = percent * interval_num;
    average_error_sep[pos] += error;
    average_error_cnt[pos]++;
    all_average_error += errors_list[i].second;
  }

			  if(false) {
    for (int i = 0; i < average_error_sep.size(); ++i) {
      if (average_error_cnt[i] != 0) {
         fprintf(stderr,"dis percent %d to %d: average error %.10lf\n" , i , i + 1 , average_error_sep[i] / average_error_cnt[i]);
      }else{
         fprintf(stderr,"dis percent %d to %d: average error non\n" , i , i + 1 );
      }
    }
  }
}

template<class T>
void check_errors(const string& obj_file_name, const string& svg_file_name, const string& method)
{
	CRichModel rich_model(obj_file_name);
	rich_model.Preprocess();

	SparseGraph<T>* s_graph = NULL;

	if (method == "dij") {
		s_graph = new Dijstra_vector<T>();
		// fprintf(stderr,"size before %d\n" , sizeof(*s_graph));
		s_graph->read_svg_file_binary((string)svg_file_name);
		// fprintf(stderr,"size after %d\n" , sizeof(*s_graph));
	}
	else if (method == "lll") {
		s_graph = new LC_LLL<T>();
		s_graph->read_svg_file_with_angle((string)svg_file_name);//load the precomputed infomation 
	}
	else if (method == "fim") {
		s_graph = new LC_FIM<T>();//temp modify
		s_graph->read_svg_file_binary((string)svg_file_name);//load the precomputed infomation 
	}
	else if (method == "hy") {
		s_graph = new LC_HY<T>();
		s_graph->read_svg_file_with_angle((string)svg_file_name);
		dynamic_cast<LC_HY<T>*>(s_graph)->setModel(rich_model);
		//s_graph->model_ptr = &rich_model;
	}
	else if (method == "dggdij") {
		s_graph = new Dijstra_vector<T>();
		s_graph->read_svg_file_with_angle((string)svg_file_name);
	}
	else{
		fprintf(stderr, "invalid choice\n");
		exit(1);
	}

	vector<pair<double, double>> erros_list;//<dis,error>
	std::mt19937 rng;
	std::uniform_int_distribution<int> uint_dist(0, s_graph->NodeNum() - 1);
	double average_time = 0;
	int iteration_times = 10;
	double max_dis = 0;
	for (int itr = 0; itr < iteration_times; ++itr)
	{
		int source_vert = uint_dist(rng);//index of source vertex
		ElapasedTime t;
		s_graph->findShortestDistance(source_vert);  //find the geodesic distance from a single source to all vertex of a mesh
		average_time += t.getTime();

		t.start();
		vector<int> sources;
		sources.push_back(source_vert);
		CICHWithFurtherPriorityQueue ich_algorithm(rich_model, sources);
		ich_algorithm.Execute();
		vector<double> correct_dis(s_graph->NodeNum());
		for (int i = 0; i < correct_dis.size(); ++i) {
			double d = ich_algorithm.m_InfoAtVertices[i].disUptodate;
			if (d == d) {
				correct_dis[i] = d;
			} else {
				correct_dis[i] = 0;
			}
		}

		max_dis = std::max(max_dis, *std::max_element(correct_dis.begin(), correct_dis.end()));

		int cnt_error_dis = 0;
		double average_error = 0;
		int cnt_finite = 0;
		for (int i = 0; i < correct_dis.size(); ++i) {
			if (fabs(correct_dis[i]) < 1e-7) continue;
			double dis = s_graph->distanceToSource(i);
			if (!_finite(correct_dis[i]) || !_finite(dis)) continue;
			double error = fabs(dis - correct_dis[i]) / correct_dis[i];
			if (error > 1) {
				cnt_error_dis++;
				continue;
			}
			cnt_finite++;
			average_error += error;
			erros_list.push_back(make_pair(correct_dis[i], error));
		}
	}

	double all_average_error;
	output_error_distribute(max_dis, erros_list, all_average_error);

	fprintf(stderr, "total_average_running_time %lf\n", average_time / iteration_times);
	fprintf(stderr, "total_average_error %.10lf\n", all_average_error / erros_list.size());

	if (method == "hy") {
		int cnt = 0;

		auto& lst = split(svg_file_name, '_');
		string eps_str = lst[lst.size() - 3].substr(2);
		fprintf(stderr, "eps %s\n", eps_str.c_str());
		double eps = atof(eps_str.c_str());
		double max_error = 0;
		for (auto& e : erros_list) {
			if (e.second > eps) {
				cnt++;
			}
			max_error = max(e.second, eps);
		}
		fprintf(stderr, "max_error %lf\n", max_error);
	}

	HANDLE current_process = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;

	if (GetProcessMemoryInfo(current_process, &pmc, sizeof(pmc)))
	{
		fprintf(stderr, "PeakWorkingSetSize: %d M\n", pmc.PeakWorkingSetSize / 1024 / 1024);
	}

	delete s_graph;



}


enum FloatType
{
	DGG_FLOAT, DGG_DOUBLE
};

int main(int argc, char** argv)
{
  bool test_performance = false;
  string obj_file_name;
  string svg_file_name;
  FloatType float_type{ DGG_FLOAT };

  if (argc < 4) {
     fprintf(stderr,"parameter insufficient !\n usage SVG_LC.exe [obj_file_name] [svg_file_name] [dij  or lll or fim or hy] [flt or dbl]\n");
    return 1;
  } else if (argc == 4) {
	  float_type = DGG_FLOAT;
  } else if (argc == 5) {
	  if (string(argv[4]) == "flt") {
		  float_type = DGG_FLOAT;
	  } else if (string(argv[4]) == "dbl") {
		  float_type = DGG_DOUBLE;
	  } else {
		  printf("type error!\n");
		  exit(1);
	  }
  } 

  obj_file_name = argv[1];
  svg_file_name = argv[2];
  string method = argv[3];
  if (float_type == DGG_FLOAT) {
	  check_errors<float>(obj_file_name, svg_file_name, method);
  } else if (float_type == DGG_DOUBLE) {
	  printf("dbl\n");
	  check_errors<double>(obj_file_name, svg_file_name, method);
  }



  return 0;
}