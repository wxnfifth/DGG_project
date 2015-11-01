#include "stdafx.h"
#include <windows.h>
#include "Shlwapi.h"
#include "wxnTime.h"
#include "svg_precompute.h"
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

bool flag_first_output = true;
struct WxnBuffer{
	char* buf;
	int len;
	int capacity;
	ofstream output_file;
	WxnBuffer() {
		buf = nullptr;
		len = 0;
		capacity = 0;
	}
	WxnBuffer(const string& output_filename) {
		output_file.open(output_filename.c_str(), ios::out | ios::binary);
		capacity = 16 * 1024 * 1024;
		buf = new char[capacity];
		len = 0;
	}
	void open(const string& output_filename) {
		output_file.open(output_filename.c_str(), ios::out | ios::binary);
		capacity = 16 * 1024 * 1024;
		buf = new char[capacity];
		len = 0;
	}

	void addStruct(const void* const ptr, int struct_size)
	{
		//fwrite(ptr, struct_size, 1, file);
		//return;
		if (len + struct_size > capacity) {
			//fwrite(buf, sizeof(char), len, file);
			output_file.write(buf, len);
			len = 0;
		}
		if (struct_size > capacity) {
			fprintf(stderr, "str too large!");
			exit(1);
		}
		memcpy((void*)(buf + len), ptr, struct_size);
		len += struct_size;
	}

	void close() {
		if (len != 0) {
			output_file.write(buf, len);
			len = 0;
		}
		output_file.close();
	}

	~WxnBuffer() {
		if (buf != nullptr) {
			delete[] buf;
		}
	}
};

string get_DGG_filename(const string& input_obj_name, const string& method_name, double eps_vg, double const_for_theta)
{
	char buf[1024];
	sprintf(buf, "%s_%s%.10lf_c%.0lf.binary", input_obj_name.substr(0, input_obj_name.length() - 4).c_str(), method_name, eps_vg, const_for_theta);
	//svg_file_name = string(buf);
	return string(buf);
}

void generate_output(const string& output_file_name,stringstream& output_str,bool force_output = false)
{
    if( force_output || output_str.tellp() > 50 * 1024 * 1024 ){
        FILE* fout = NULL;
        if(flag_first_output){
            fout = fopen(output_file_name.c_str(),"w");
            flag_first_output = false;
        }else{
            fout = fopen(output_file_name.c_str(),"a");
        }
        string& str = output_str.str();
        ElapasedTime time;
        fwrite(str.c_str(),1,str.length(),fout);
        fclose(fout);
        cerr << "write to file time " << time.getTime() << "\n";
        output_str.str(std::string());
        output_str.clear();
    }
}

void svg_precompute_fix_neighbor(const string& input_obj_name, const int fixed_k, string& svg_file_name) {
  //Step 1. Initialize models
  CRichModel model(input_obj_name);
  model.Preprocess();

  int begin_vertex_index = 0;

  int end_vertex_index = model.GetNumOfVerts() - 1;

  svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
    + "_SVG_k" + to_string(fixed_k) +  ".binary";

  ofstream output_file (svg_file_name.c_str() , ios::out | ios::binary);
  int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
  HeadOfSVG head_of_svg(begin_vertex_index , end_vertex_index , num_of_vertex );     
  output_file.write((char*)&head_of_svg , sizeof(head_of_svg));

  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //#pragma omp parallel for
  for (int tmp_source = begin_vertex_index;tmp_source <= end_vertex_index;++tmp_source) {

    if (time_once.getTime() -  past_time > 5 ) {
      past_time = time_once.getTime();
      char buf[128];
      sprintf(buf, "Computed %.0lf percent", (double) tmp_source  * 100. / (end_vertex_index - begin_vertex_index));
      time_once.printTime(buf );
    }
    ElapasedTime dis_time;
    int source_index = tmp_source;
    double farthest = 0.0;
    //Step 2: Construct a group of source points;
    vector<int> sources;
    sources.push_back(source_index);
    //Step 3: Construct a new algorithm object
    CICHWithFurtherPriorityQueue alg(model, sources);
    //Step 4: Locally propagate wavefronts and stop at the prescribed geodesic distance.
    set<int> fixedDests;
    //The first parameter is the distance threshold, 
    //and the second is to return those vertices where the geodesic distance makes sense.
    double max_radius = 1e10;
    alg.ExecuteLocally_SVG(max_radius, fixedDests,fixed_k);
    dis_time_total += dis_time.getTime();
    //printf("Totally collected: %d\n", fixedDests.size());
    set<int>::iterator itr;
    int cnt = 0;
    vector<pair<int,double>> dests;

    struct node {
      int id;
      double dis;
      int operator<(const node & other) const{
        return dis < other.dis;
      }
    };
    vector<node> covered_points;
    covered_points.resize(fixedDests.size());
    int _index = 0;
    for (itr = fixedDests.begin(); itr != fixedDests.end(); ++itr) {
      int v = *itr;
      map<int, CICHWithFurtherPriorityQueue::InfoAtVertex>::const_iterator it = alg.m_InfoAtVertices.find(v);
      int indexOfParent = it->second.indexOfParent;
      int parent = 0;
      covered_points[_index].id = v;
      covered_points[_index].dis = it->second.disUptodate;
      _index ++;
      if (farthest < it->second.disUptodate){
        farthest = it->second.disUptodate;
      }
      if (it->second.fParentIsPseudoSource) {
        parent = indexOfParent;
      } else {
        parent = it->second.indexOfRootVertOfParent;
      }
      //if (parent == sources[0]) {
      if (true) {
        dests.push_back(pair<int,double>(v , alg.m_InfoAtVertices[v].disUptodate));
        ++cnt;
      }
    }
    std::sort(covered_points.begin(), covered_points.end());
    std::map<int, int> mp;
    for(int i = 0; i < covered_points.size(); ++i){
      mp[covered_points[i].id] = i;
    }

    BodyHeadOfSVG body_header(source_index , dests.size());
    output_file.write((char*)&body_header , sizeof(body_header));

    vector<BodyPartOfSVGWithK> body_parts(dests.size());
    for(int i = 0; i < dests.size(); ++i) {
      BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , mp[dests[i].first]);
      body_parts[i] = body_part;
    }
    sort(body_parts.begin() , body_parts.end());
    for (int i = 0; i < body_parts.size(); ++i) {
      BodyPartOfSVG b = body_parts[i];
      output_file.write( (char*)&b , sizeof(b));
    }
  }

  time_once.printTime("time past ");
  output_file.close();

}


void svg_precompute_fix_neighbor_debug(const string& input_obj_name, const int fixed_k) {
	//Step 1. Initialize models
	CRichModel model(input_obj_name);
	model.Preprocess();
	ElapasedTime time_once;

	double dis_time_total(0);
	double past_time(0);
	//#pragma omp parallel for
	//for (int tmp_source = begin_vertex_index;tmp_source <= end_vertex_index;++tmp_source) {
	for (int tmp_source = 0; tmp_source < 10; ++tmp_source) {
		if (time_once.getTime() -  past_time > 5 ) {
			past_time = time_once.getTime();
			char buf[128];
			time_once.printTime(buf);
		}
		ElapasedTime dis_time;
		int source_index = tmp_source;
		double farthest = 0.0;
		//Step 2: Construct a group of source points;
		vector<int> sources;
		sources.push_back(source_index);
		//Step 3: Construct a new algorithm object
		CICHWithFurtherPriorityQueue alg(model, sources);
		//Step 4: Locally propagate wavefronts and stop at the prescribed geodesic distance.
		set<int> fixedDests;
		//The first parameter is the distance threshold, 
		//and the second is to return those vertices where the geodesic distance makes sense.
		double max_radius = 1e10;
		alg.ExecuteLocally_SVG(max_radius, fixedDests,fixed_k);
		dis_time_total += dis_time.getTime();
		//printf("Totally collected: %d\n", fixedDests.size());
		set<int>::iterator itr;
		int cnt = 0;
		vector<pair<int,double>> dests;

		struct node {
			int id;
			double dis;
			int operator<(const node & other) const{
				return dis < other.dis;
			}
		};
		vector<node> covered_points;
		covered_points.resize(fixedDests.size());
		int _index = 0;
		for (itr = fixedDests.begin(); itr != fixedDests.end(); ++itr) {
			int v = *itr;
			map<int, CICHWithFurtherPriorityQueue::InfoAtVertex>::const_iterator it = alg.m_InfoAtVertices.find(v);
			int indexOfParent = it->second.indexOfParent;
			int parent = 0;
			covered_points[_index].id = v;
			covered_points[_index].dis = it->second.disUptodate;
			_index ++;
			if (farthest < it->second.disUptodate){
				farthest = it->second.disUptodate;
			}
			if (it->second.fParentIsPseudoSource) {
				parent = indexOfParent;
			} else {
				parent = it->second.indexOfRootVertOfParent;
			}
			//if (parent == sources[0]) {
			if (true) {
				dests.push_back(pair<int,double>(v , alg.m_InfoAtVertices[v].disUptodate));
				++cnt;
			}
		}
		std::sort(covered_points.begin(), covered_points.end());
		std::map<int, int> mp;
		for(int i = 0; i < covered_points.size(); ++i){
			mp[covered_points[i].id] = i;
		}

		vector<BodyPartOfSVGWithK> body_parts(dests.size());
		for(int i = 0; i < dests.size(); ++i) {
			BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , mp[dests[i].first]);
			body_parts[i] = body_part;
		}
		sort(body_parts.begin() , body_parts.end());
		printf("source %d neigh_size %d\n", tmp_source, body_parts.size());
	}

	time_once.printTime("time past ");

}


#if 0
void svg_precompute_jiajun_output(const string& input_obj_name, double eps_vg, string& svg_file_name)
{
  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;

  clock_t start = clock();
  bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(),points,faces, realIndex, originalVertNum);
  if(!success)
  {
     fprintf(stderr, "something is wrong with the input file" );
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  fprintf(stderr,  "loading model took %lf seconds" , (double)(end - start) / (double)CLOCKS_PER_SEC );

  svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
    + "_VG" + to_string(eps_vg) +  ".binary";

  int begin_vertex_index = 0;

  int end_vertex_index = points.size() / 3 - 1;

	char buf[1024];
  string prefix = input_obj_name.substr(0,input_obj_name.length()-4);
  sprintf(buf, "%s_K%.6lf_degree.txt", prefix.c_str(), eps_vg);
	FILE * fd = fopen(buf, "w");
  fprintf(fd,"%d\n" , mesh.vertices().size());
	sprintf(buf, "%s_K%.6lf.edge", prefix.c_str(), eps_vg);
	FILE * fe = fopen(buf, "wb");
	int total_edge = 0;
  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //end_vertex_index = 1;

  for (int tmp_source = begin_vertex_index;tmp_source <= end_vertex_index;++tmp_source) {

    if (time_once.getTime() -  past_time > 5 ) {
      past_time = time_once.getTime();
      char buf[128];
      sprintf(buf, "Computed %.0lf percent", (double) tmp_source  * 100. / (end_vertex_index - begin_vertex_index));
      time_once.printTime(buf );
    }
    ElapasedTime dis_time;
    int source_index = tmp_source;
    vector<int> srcs;
    srcs.push_back(source_index);
    vector<geodesic::SurfacePoint> sources;
    for (unsigned i = 0; i < srcs.size(); ++i)
    {
      srcs[i] %= originalVertNum;
      srcs[i] = realIndex[srcs[i]];
      sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
    }

    geodesic::GeodesicAlgorithmBase *algorithm;

    algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh); 

    double step = 1.0;
    const double eta = 100;

    algorithm->step = step;
    algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
    map<int,double> fixedDests;
    algorithm->propagate_vg(sources, eps_vg, fixedDests);



    vector<pair<int,double>> dests;

    struct node {
      int id;
      double dis;
      int operator<(const node & other) const{
        return dis < other.dis;
      }
    };

    //printf("tmp_source %d dest size %d\n" , tmp_source, fixedDests.size());
    vector<node> covered_points;
    covered_points.resize(fixedDests.size());
    int _index = 0;
    for (auto& d:fixedDests) {
      geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
      double dis;
      //printf("line 274\n");
      algorithm->best_source(dest_p, dis);
      //printf("line 276\n");
      covered_points[_index].id = d.first;
      covered_points[_index].dis = d.second;
      _index ++;
      //if (fabs(dis - d.second) > 1e-6) {
      //  fprintf(stderr,"dis %lf origin %lf\n", dis, d.second);
      //}
      dests.push_back(make_pair(d.first,d.second));
    }
    delete algorithm;
    algorithm = NULL;

    std::sort(covered_points.begin(), covered_points.end());
    std::map<int, int> mp;
    for(int i = 0; i < covered_points.size(); ++i){
      mp[covered_points[i].id] = i;
    }

    fprintf(fd, "%d\n", dests.size());

    //BodyHeadOfSVG body_header(source_index , dests.size());
    //output_file.write((char*)&body_header , sizeof(body_header));

    vector<BodyPartOfSVGWithK> body_parts(dests.size());
    for(int i = 0; i < dests.size(); ++i) {
      BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , mp[dests[i].first]);
      //BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , i);
      body_parts[i] = body_part;
    }
    sort(body_parts.begin() , body_parts.end());
    for (int i = 0; i < body_parts.size(); ++i) {
      //BodyPartOfSVG b = body_parts[i];
      //output_file.write( (char*)&b , sizeof(b));
      SVGEdge e;
      e.v = body_parts[i].dest_index;
      e.dis = body_parts[i].dest_dis;
      fwrite(&e, sizeof(e), 1, fe); //将edge写入文件
    }

  }



  time_once.printTime("time past ");

  fclose(fd);
  fclose(fe);


}
#endif

void svg_precompute_jiajun(const string& input_obj_name, double eps_vg, string& svg_file_name)
{
  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;

  clock_t start = clock();
  bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(),points,faces, realIndex, originalVertNum);
  if(!success)
  {
     fprintf(stderr, "something is wrong with the input file" );
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  fprintf(stderr,  "loading model took %lf seconds" , (double)(end - start) / (double)CLOCKS_PER_SEC );

  svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
    + "_VG" + to_string(eps_vg) +  ".binary";

  int begin_vertex_index = 0;

  int end_vertex_index = points.size() / 3 - 1;

  ofstream output_file (svg_file_name.c_str() , ios::out | ios::binary);
  int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
  HeadOfSVG head_of_svg(begin_vertex_index , end_vertex_index , num_of_vertex );     
  output_file.write((char*)&head_of_svg , sizeof(head_of_svg));
  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //end_vertex_index = 0;

  for (int tmp_source = begin_vertex_index;tmp_source <= end_vertex_index;++tmp_source) {

    if (time_once.getTime() -  past_time > 5 ) {
      past_time = time_once.getTime();
      char buf[128];
      sprintf(buf, "Computed %.0lf percent", (double) tmp_source  * 100. / (end_vertex_index - begin_vertex_index));
      time_once.printTime(buf );
    }
    ElapasedTime dis_time;
    int source_index = tmp_source;
    vector<int> srcs;
    srcs.push_back(source_index);
    vector<geodesic::SurfacePoint> sources;
    for (unsigned i = 0; i < srcs.size(); ++i)
    {
      srcs[i] %= originalVertNum;
      srcs[i] = realIndex[srcs[i]];
      sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
    }

    geodesic::GeodesicAlgorithmBase *algorithm;

    algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh); 

    double step = 1.0;
    const double eta = 100;

    algorithm->step = step;
    algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
    map<int,double> fixedDests;
    algorithm->propagate_vg(sources, eps_vg, fixedDests);



    vector<pair<int,double>> dests;

    struct node {
      int id;
      double dis;
      int operator<(const node & other) const{
        return dis < other.dis;
      }
    };

    //printf("tmp_source %d dest size %d\n" , tmp_source, fixedDests.size());
    vector<node> covered_points;
    covered_points.resize(fixedDests.size());
    int _index = 0;
    for (auto& d:fixedDests) {
      geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
      double dis;
      //printf("line 274\n");
      algorithm->best_source(dest_p, dis);
      //printf("line 276\n");
      covered_points[_index].id = d.first;
      covered_points[_index].dis = d.second;
      _index ++;
      //printf("v %d dis %lf\n" , d.first, d.second);
      //if (fabs(dis - d.second) > 1e-6) {
      //  fprintf(stderr,"dis %lf origin %lf\n", dis, d.second);
      //}
      dests.push_back(make_pair(d.first,d.second));
    }
    delete algorithm;
    algorithm = NULL;

    std::sort(covered_points.begin(), covered_points.end());
    std::map<int, int> mp;
    for(int i = 0; i < covered_points.size(); ++i){
      mp[covered_points[i].id] = i;
    }
    BodyHeadOfSVG body_header(source_index , dests.size());
    output_file.write((char*)&body_header , sizeof(body_header));

    vector<BodyPartOfSVGWithK> body_parts(dests.size());
    for(int i = 0; i < dests.size(); ++i) {
      BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , mp[dests[i].first]);
      //BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , i);
      body_parts[i] = body_part;
    }
    sort(body_parts.begin() , body_parts.end());
    for (int i = 0; i < body_parts.size(); ++i) {
      BodyPartOfSVG b = body_parts[i];
      output_file.write( (char*)&b , sizeof(b));
    }
  }
  time_once.printTime("time past ");
  output_file.close();
}

void svg_precompute_mmp(const string& input_obj_name, const int fixed_k, string& svg_file_name)
{
	std::vector<double> points;
	std::vector<unsigned> faces;
	std::vector<int> realIndex;
	int originalVertNum = 0;

	clock_t start = clock();
	bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(), points, faces, realIndex, originalVertNum);
	if (!success)
	{
		fprintf(stderr, "something is wrong with the input file");
		return;
	}
	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);		//create internal
	clock_t end = clock();
	fprintf(stderr, "loading model took %lf seconds", (double)(end - start) / (double)CLOCKS_PER_SEC);

	svg_file_name = input_obj_name.substr(0, input_obj_name.length() - 4)
		+ "_SVGMMP_k" + to_string(fixed_k) + ".binary";

	int begin_vertex_index = 0;

	int end_vertex_index = points.size() / 3 - 1;

	ofstream output_file(svg_file_name.c_str(), ios::out | ios::binary);
	int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
	HeadOfSVG head_of_svg(begin_vertex_index, end_vertex_index, num_of_vertex);
	output_file.write((char*)&head_of_svg, sizeof(head_of_svg));
	ElapasedTime time_once;

	double dis_time_total(0);
	double past_time(0);
	//end_vertex_index = 10;

	for (int tmp_source = begin_vertex_index; tmp_source <= end_vertex_index; ++tmp_source) {

		if (time_once.getTime() - past_time > 5) {
			past_time = time_once.getTime();
			char buf[128];
			sprintf(buf, "Computed %.0lf percent", (double)tmp_source  * 100. / (end_vertex_index - begin_vertex_index));
			time_once.printTime(buf);
		}
		ElapasedTime dis_time;
		int source_index = tmp_source;
		vector<int> srcs;
		srcs.push_back(source_index);
		vector<geodesic::SurfacePoint> sources;
		for (unsigned i = 0; i < srcs.size(); ++i)
		{
			srcs[i] %= originalVertNum;
			srcs[i] = realIndex[srcs[i]];
			sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
		}

		geodesic::GeodesicAlgorithmBase *algorithm;

		algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh);

		double step = 1.0;
		const double eta = 100;

		algorithm->step = step;
		algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
		map<int, double> fixedDests;
		algorithm->propagate_local(sources, fixed_k, fixedDests);

		vector<pair<int, double>> dests;

		struct node {
			int id;
			double dis;
			int operator<(const node & other) const{
				return dis < other.dis;
			}
		};

		//printf("tmp_source %d dest size %d\n" , tmp_source, dest_verts.size());
		vector<node> covered_points;
		covered_points.resize(fixedDests.size());
		int _index = 0;
		for (auto& d : fixedDests) {
			geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
			double dis;
			//printf("line 274\n");
			algorithm->best_source(dest_p, dis);
			//printf("line 276\n");
			covered_points[_index].id = d.first;
			covered_points[_index].dis = d.second;
			_index++;
			if (fabs(dis - d.second) > 1e-6) {
				fprintf(stderr, "dis %lf origin %lf\n", dis, d.second);
			}
			dests.push_back(make_pair(d.first, d.second));
		}
		delete algorithm;
		algorithm = NULL;

		std::sort(covered_points.begin(), covered_points.end());
		std::map<int, int> mp;
		for (int i = 0; i < covered_points.size(); ++i){
			mp[covered_points[i].id] = i;
		}
		BodyHeadOfSVG body_header(source_index, dests.size());
		output_file.write((char*)&body_header, sizeof(body_header));

		vector<BodyPartOfSVGWithK> body_parts(dests.size());
		for (int i = 0; i < dests.size(); ++i) {
			BodyPartOfSVGWithK body_part(dests[i].first, dests[i].second, mp[dests[i].first]);
			//BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , i);
			body_parts[i] = body_part;
		}
		sort(body_parts.begin(), body_parts.end());
		for (int i = 0; i < body_parts.size(); ++i) {
			BodyPartOfSVG b = body_parts[i];
			output_file.write((char*)&b, sizeof(b));
		}
	}

	time_once.printTime("time past ");

	output_file.close();


}

void svg_precompute(const string& input_obj_name, const int fixed_k, string& svg_file_name) {
  //Step 1. Initialize models
  CRichModel model(input_obj_name);
  model.Preprocess();

  int begin_vertex_index = 0;

  int end_vertex_index = model.GetNumOfVerts() - 1;

  svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
    + "_SVG_k" + to_string(fixed_k) +  ".binary";

  ofstream output_file (svg_file_name.c_str() , ios::out | ios::binary);
  int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
  HeadOfSVG head_of_svg(begin_vertex_index , end_vertex_index , num_of_vertex );     
  output_file.write((char*)&head_of_svg , sizeof(head_of_svg));

  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //#pragma omp parallel for
  for (int tmp_source = begin_vertex_index;tmp_source <= end_vertex_index;++tmp_source) {

    if (time_once.getTime() -  past_time > 5 ) {
      past_time = time_once.getTime();
      char buf[128];

	  double percent = (double)tmp_source  * 100. / double(end_vertex_index - begin_vertex_index + 1);
	  double current_time = time_once.getTime();
	  double remain_time = current_time / percent * (100 - percent);
	  printf("Computed %.0lf percent, time %lf, estimate_remain_time %lf\n",
		  percent, current_time, remain_time);
    }
    ElapasedTime dis_time;
    int source_index = tmp_source;
    double farthest = 0.0;
    //Step 2: Construct a group of source points;
    vector<int> sources;
    sources.push_back(source_index);
    //Step 3: Construct a new algorithm object
    CICHWithFurtherPriorityQueue alg(model, sources);
    //Step 4: Locally propagate wavefronts and stop at the prescribed geodesic distance.
    set<int> fixedDests;
    //The first parameter is the distance threshold, 
    //and the second is to return those vertices where the geodesic distance makes sense.
    double max_radius = 1e10;
    alg.ExecuteLocally_SVG(max_radius, fixedDests,fixed_k);
    dis_time_total += dis_time.getTime();
    //printf("Totally collected: %d\n", fixedDests.size());
    set<int>::iterator itr;
    int cnt = 0;
    vector<pair<int,double>> dests;

    struct node {
      int id;
      double dis;
      int operator<(const node & other) const{
        return dis < other.dis;
      }
    };
    vector<node> covered_points;
    covered_points.resize(fixedDests.size());
    int _index = 0;
    for (itr = fixedDests.begin(); itr != fixedDests.end(); ++itr) {
      int v = *itr;
      map<int, CICHWithFurtherPriorityQueue::InfoAtVertex>::const_iterator it = alg.m_InfoAtVertices.find(v);
      int indexOfParent = it->second.indexOfParent;
      int parent = 0;
      covered_points[_index].id = v;
      covered_points[_index].dis = it->second.disUptodate;
      _index ++;
      if (farthest < it->second.disUptodate){
        farthest = it->second.disUptodate;
      }
      if (it->second.fParentIsPseudoSource) {
        parent = indexOfParent;
      } else {
        parent = it->second.indexOfRootVertOfParent;
      }
      if (parent == sources[0]) {
        dests.push_back(pair<int,double>(v , alg.m_InfoAtVertices[v].disUptodate));
        ++cnt;
      }
    }
    std::sort(covered_points.begin(), covered_points.end());
    std::map<int, int> mp;
    for(int i = 0; i < covered_points.size(); ++i){
      mp[covered_points[i].id] = i;
    }

    BodyHeadOfSVG body_header(source_index , dests.size());
    output_file.write((char*)&body_header , sizeof(body_header));

    vector<BodyPartOfSVGWithK> body_parts(dests.size());
    for(int i = 0; i < dests.size(); ++i) {
      BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , mp[dests[i].first]);
      body_parts[i] = body_part;
    }
    sort(body_parts.begin() , body_parts.end());
	//printf("source %d, neigh_size %d\n", tmp_source, body_parts.size());
    for (int i = 0; i < body_parts.size(); ++i) {
      BodyPartOfSVG b = body_parts[i];
      output_file.write( (char*)&b , sizeof(b));
    }
  }

  time_once.printTime("time past ");
  output_file.close();

}

void svg_precompute_debug(const string& input_obj_name, const int fixed_k) {
	//Step 1. Initialize models
	CRichModel model(input_obj_name);
	model.Preprocess();

	int begin_vertex_index = 0;

	int end_vertex_index = 10;

	int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
	HeadOfSVG head_of_svg(begin_vertex_index, end_vertex_index, num_of_vertex);

	ElapasedTime time_once;
	double dis_time_total(0);
	double past_time(0);
	//#pragma omp parallel for
	for (int tmp_source = begin_vertex_index; tmp_source <= end_vertex_index; ++tmp_source) {

		if (time_once.getTime() - past_time > 5) {
			past_time = time_once.getTime();
			char buf[128];
			sprintf(buf, "Computed %.0lf percent", (double)tmp_source  * 100. / (end_vertex_index - begin_vertex_index));
			time_once.printTime(buf);
		}
		ElapasedTime dis_time;
		int source_index = tmp_source;
		double farthest = 0.0;
		//Step 2: Construct a group of source points;
		vector<int> sources;
		sources.push_back(source_index);
		//Step 3: Construct a new algorithm object
		CICHWithFurtherPriorityQueue alg(model, sources);
		//Step 4: Locally propagate wavefronts and stop at the prescribed geodesic distance.
		set<int> fixedDests;
		//The first parameter is the distance threshold, 
		//and the second is to return those vertices where the geodesic distance makes sense.
		double max_radius = 1e10;
		alg.ExecuteLocally_SVG(max_radius, fixedDests, fixed_k);
		dis_time_total += dis_time.getTime();
		//printf("Totally collected: %d\n", fixedDests.size());
		set<int>::iterator itr;
		int cnt = 0;
		vector<pair<int, double>> dests;

		struct node {
			int id;
			double dis;
			int operator<(const node & other) const{
				return dis < other.dis;
			}
		};
		vector<node> covered_points;
		covered_points.resize(fixedDests.size());
		int _index = 0;
		for (itr = fixedDests.begin(); itr != fixedDests.end(); ++itr) {
			int v = *itr;
			map<int, CICHWithFurtherPriorityQueue::InfoAtVertex>::const_iterator it = alg.m_InfoAtVertices.find(v);
			int indexOfParent = it->second.indexOfParent;
			int parent = 0;
			covered_points[_index].id = v;
			covered_points[_index].dis = it->second.disUptodate;
			_index++;
			if (farthest < it->second.disUptodate){
				farthest = it->second.disUptodate;
			}
			if (it->second.fParentIsPseudoSource) {
				parent = indexOfParent;
			}
			else {
				parent = it->second.indexOfRootVertOfParent;
			}
			if (parent == sources[0]) {
				dests.push_back(pair<int, double>(v, alg.m_InfoAtVertices[v].disUptodate));
				++cnt;
			}
		}
		std::sort(covered_points.begin(), covered_points.end());
		std::map<int, int> mp;
		for (int i = 0; i < covered_points.size(); ++i){
			mp[covered_points[i].id] = i;
		}

		vector<BodyPartOfSVGWithK> body_parts(dests.size());
		for (int i = 0; i < dests.size(); ++i) {
			BodyPartOfSVGWithK body_part(dests[i].first, dests[i].second, mp[dests[i].first]);
			body_parts[i] = body_part;
		}
		sort(body_parts.begin(), body_parts.end());
		printf("source %d body_parts %d\n", tmp_source, body_parts.size());

		{
			vector<int> dests;
			for (auto& b : body_parts) {
				dests.push_back(b.dest_index);
			}
			CylinderPath cylinder_path(0.01);
			cylinder_path.cntGeodesicPaths(model, tmp_source, dests);
		}

	}
	time_once.printTime("time past ");

}


void svg_precompute_hy(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta)
{

  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  fprintf(stderr,"******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;
  CRichModel model(input_obj_name);
  model.Preprocess();

  clock_t start = clock();
  bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(),points,faces, realIndex, originalVertNum);
  if(!success)
  {
     fprintf(stderr, "something is wrong with the input file" );
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  fprintf(stderr,  "loading model took %lf seconds" , (double)(end - start) / (double)CLOCKS_PER_SEC );

  //svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
  //  + "_HY" + to_string(eps_vg) +  ".binary";
  svg_file_name = get_DGG_filename(input_obj_name, "HY", eps_vg, const_for_theta);

  int begin_vertex_index = 0;

  int end_vertex_index = points.size() / 3 - 1;

  ofstream output_file (svg_file_name.c_str() , ios::out | ios::binary);
  int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
  HeadOfSVG head_of_svg(begin_vertex_index , end_vertex_index , num_of_vertex );     
  output_file.write((char*)&head_of_svg , sizeof(head_of_svg));
  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //end_vertex_index = 0;

  for (int tmp_source = begin_vertex_index;tmp_source <= end_vertex_index;++tmp_source) {

    if (time_once.getTime() -  past_time > 5 ) {
      past_time = time_once.getTime();
      char buf[128];
      sprintf(buf, "Computed %.0lf percent", (double) tmp_source  * 100. / (end_vertex_index - begin_vertex_index));
      time_once.printTime(buf );
    }
    ElapasedTime dis_time;
    int source_index = tmp_source;
    vector<int> srcs;
    srcs.push_back(source_index);
    vector<geodesic::SurfacePoint> sources;
    for (unsigned i = 0; i < srcs.size(); ++i)
    {
      srcs[i] %= originalVertNum;
      srcs[i] = realIndex[srcs[i]];
      sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
    }

    geodesic::GeodesicAlgorithmBase *algorithm;

    algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh); 

    double step = 1.0;
    const double eta = 100;

    algorithm->step = step;
    algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
    map<int,double> fixedDests;
    algorithm->propagate_vg(sources, eps_vg, fixedDests);



    vector<pair<int,double>> dests;

    struct node {
      int id;
      double dis;
      int operator<(const node & other) const{
        return dis < other.dis;
      }
    };

    //printf("tmp_source %d dest size %d\n" , tmp_source, fixedDests.size());
    vector<node> covered_points;
    covered_points.resize(fixedDests.size());
    int _index = 0;
    for (auto& d:fixedDests) {
      geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
      double dis;
      //printf("line 274\n");
      algorithm->best_source(dest_p, dis);
      //printf("line 276\n");
      covered_points[_index].id = d.first;
      covered_points[_index].dis = d.second;
      _index ++;
      //printf("v %d dis %lf\n" , d.first, d.second);
      //if (fabs(dis - d.second) > 1e-6) {
      //  fprintf(stderr,"dis %lf origin %lf\n", dis, d.second);
      //}
      dests.push_back(make_pair(d.first,d.second));
    }

    std::sort(covered_points.begin(), covered_points.end());
    std::map<int, int> mp;
    for(int i = 0; i < covered_points.size(); ++i){
      mp[covered_points[i].id] = i;
    }
    BodyHeadOfSVG body_header(source_index , dests.size());
    output_file.write((char*)&body_header , sizeof(body_header));

    vector<BodyPartOfSVGWithK> body_parts(dests.size());
    for(int i = 0; i < dests.size(); ++i) {
      BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , mp[dests[i].first]);
      //BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , i);
      body_parts[i] = body_part;
    }
    sort(body_parts.begin() , body_parts.end());

    vector<double> angles(body_parts.size());
    for (int i = 0; i < body_parts.size(); ++i) {
      //double distance;
      //auto p(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
      //unsigned best_source = algorithm->best_source(p, distance);	
      auto dest_vert(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
      vector<geodesic::SurfacePoint> path; 
      algorithm->trace_back(dest_vert, path);
      geodesic::SurfacePoint p;
      if (path.size() > 1) {
        p = *(path.rbegin()+1);
      } else {
        p = (path[0]);
      }
 //         VERTEX,
 //   EDGE,
 //   FACE,
	//UNDEFINED_POINT
      auto& neighs = model.Neigh(source_index);
      vector<double> sum_angle(neighs.size()+1);
      sum_angle[0] = 0;
      for (int j = 1; j <= neighs.size(); ++j) {
        sum_angle[j] = sum_angle[j-1] + neighs[j-1].second;
      }

      double angle = 0;
      if ( p.type() == geodesic::VERTEX) {
        //printf("vertex\n");
        bool flag_found = false;
        for (int j = 0; j < neighs.size(); ++j) {
          auto& neigh = neighs[j];
          if ( p.base_element()->id() == model.Edge(neigh.first).indexOfRightVert) {
            //printf("yes\n");
            flag_found = true;
            angle = sum_angle[j];
            break;
          }
        }
        if (!flag_found) {
          angle = 0;
          //printf("vertex %d source %d\n" , p.base_element()->id(), source_index);
        }

      } else if(p.type() == geodesic::EDGE) {
        
        //printf("edge\n");
        int v0 = p.base_element()->adjacent_vertices()[0]->id();
        int v1 = p.base_element()->adjacent_vertices()[1]->id();
        bool flag = false;
        for (int j = 0; j < neighs.size(); ++j) {
          auto& neigh = neighs[j];
          if (v0 == model.Edge(neigh.first).indexOfRightVert) {
            int jminus1 = (j-1+neighs.size())%neighs.size();
            int vjminus1 = model.Edge(neighs[jminus1].first).indexOfRightVert;
            int jplus1 = (j+1)%neighs.size();
            int vjplus1 = model.Edge(neighs[jplus1].first).indexOfRightVert;
            //printf("v1 %d j -1 %d j + 1 %d\n" , v1, vjminus1, vjplus1); 
            CPoint3D p_cpoint3d(p.x(),p.y(),p.z());
            
            if (v1 == vjminus1) {//v1 first
              double l = model.Edge(neighs[jminus1].first).length;
              double r = (model.Vert(source_index) - p_cpoint3d).Len();
              double b = (model.Vert(vjminus1) - p_cpoint3d).Len();
              angle = sum_angle[jminus1] + acos((l * l + r * r - b * b) / (2 * l * r));
            } else if (v1 == vjplus1) {//v0 first
              double l = model.Edge(neighs[j].first).length;
              double r = (model.Vert(source_index) - p_cpoint3d).Len();
              double b = (model.Vert(v0) - p_cpoint3d).Len();
              angle = sum_angle[j] + acos((l * l + r * r - b * b) / (2 * l * r));
            }else{
              fprintf(stderr,"fuck error line 680\n");
            }
            flag = true;
            break;
          }
        }
        if (!flag) {
          fprintf(stderr,"flag %d\n" , flag);
        }

      } else{
        fprintf(stderr,"fuck error face\n");
      }
    
      angles[i] = angle;
    }

    vector<BodyPartOfSVGWithAngle> body_parts_with_angle(body_parts.size());
    for (int i = 0; i < body_parts.size(); ++i) {
      BodyPartOfSVG b = body_parts[i];
      BodyPartOfSVGWithAngle b_with_angle(b.dest_index, b.dest_dis, angles[i],0,0);
      body_parts_with_angle[i] = b_with_angle;
      //output_file.write( (char*)&b_with_angle , sizeof(b_with_angle));
    }
    sort(body_parts_with_angle.begin(), body_parts_with_angle.end());
    //for (auto& b:body_parts_with_angle) {
    //  fprintf(stderr,"b %lf " , b.angle);
    //}
    //printf("\n");
    float angle_sum = model.AngleSum(source_index);
    vector<double> tmp_angles(body_parts_with_angle.size()*2);
    for (int i = 0; i < body_parts_with_angle.size(); ++i) {
      tmp_angles[i] = body_parts_with_angle[i].angle;
    }
    for (int i = body_parts_with_angle.size(); i < tmp_angles.size(); ++i) {
      tmp_angles[i] = body_parts_with_angle[i-body_parts_with_angle.size()].angle + angle_sum;
    }
    //if (source_index == 0) {
    //  for (auto& p:tmp_angles) {
    //    fprintf(stderr,"p %lf " , p);
    //  }
    //}
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

      int start_pos = lower_bound(tmp_angles.begin(),tmp_angles.end(),start_angle) - tmp_angles.begin();
      if (start_pos > 0) start_pos--;
      int end_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), end_angle) - tmp_angles.begin();
      //printf("start_angle %lf start_pos_angle %lf end_angle  %lf end_pos_angle %lf\n" , start_angle, tmp_angles[start_pos], end_angle, tmp_angles[end_pos]);
      if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
      if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
      body_parts_with_angle[i].begin_pos = start_pos;
      body_parts_with_angle[i].end_pos = end_pos;
    }
    for (auto& b:body_parts_with_angle) {
      output_file.write((char*)&b , sizeof(b));    
    }

    delete algorithm;
    algorithm = NULL;
  }
  time_once.printTime("time past ");

  output_file.close();


}

void getDisAndAngles(const CRichModel& model, int source, double eps_vg, vector<int>& dests, vector<double>& angles, vector<double>& dis)
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
		if (source == 0 && d == 20) {
			printf("source %d dest %d parent %d dis %.10lf\n", source, d, parent, it->second.disUptodate);
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

void getFanOutput(const vector<int>& dests, const vector<double>& angles,
				  const vector<double>& dis, const CRichModel& model,  
				  double theta, int source ,
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


void svg_precompute_ich_vert(int source , int total_vert, CRichModel& model,
						   	 double eps_vg, double theta, WxnBuffer& wxn_buffer,
							 bool is_debug_mode, ElapasedTime& time_once,
							 double& past_time, double& average_degree)
{
	if (time_once.getTime() - past_time > 5) {
		past_time = time_once.getTime();
		double percent = (double)source  * 100. / double(total_vert);
		double current_time = time_once.getTime();
		double remain_time = current_time / percent * (100 - percent);
		printf("Computed %.0lf percent, time %lf, estimate_remain_time %lf\n",
					percent, current_time, remain_time);
		//time_once.printTime(buf);
	}

	vector<int> dests;
	vector<double> angles;
	vector<double> dis;
	getDisAndAngles(model, source, eps_vg, dests, angles, dis);



	vector<BodyPartOfSVGWithAngle> body_parts_with_angle;
	getFanOutput(dests, angles, dis, model, theta, source, body_parts_with_angle);
	if (is_debug_mode) {
		CylinderPath cylinder_path(0.01);
		cylinder_path.cntGeodesicPaths(model, source, dests);
		printf("neigh sz %d\n", body_parts_with_angle.size());
	}
	BodyHeadOfSVG body_header(source, body_parts_with_angle.size());
	//output_file.write((char*)&body_header, sizeof(body_header));
	wxn_buffer.addStruct((void*)&body_header, sizeof(body_header));
	for (auto& b : body_parts_with_angle) {
		//output_file.write((char*)&b, sizeof(b));
		wxn_buffer.addStruct((void*)&b, sizeof(b));
	}
	average_degree += body_parts_with_angle.size();
}



void svg_precompute_ich_debug(const string& input_obj_name, const string& debug_svg_filename, const string& neigh_filename)
{
	//void readInputFile(const string& svg_file_name,
	int node_number;
	std::ifstream input_file(debug_svg_filename, std::ios::in | std::ios::binary);
	HeadOfSVG head_of_svg;
	input_file.read((char*)&head_of_svg, sizeof(head_of_svg));
	head_of_svg.print();
	node_number = head_of_svg.num_of_vertex;
	//std::ofstream output_file(neigh_filename, std::ios::out);
	FILE* output_file = fopen(neigh_filename.c_str(), "w");
	ElapasedTime time_once;
	CRichModel model(input_obj_name);
	model.Preprocess();

	double past_time = 0;
	double average_degree = 0;
	for (int i = 0; i < node_number; ++i) {
		if (time_once.getTime() - past_time > 5) {
			past_time = time_once.getTime();
			char buf[128];
			sprintf(buf, "Computed %.0lf percent", (double)i  * 100. / double(node_number));
			time_once.printTime(buf);
		}

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
		//printf("nei %d\n", body_head.neighbor_num);
		average_degree += body_head.neighbor_num;
		fprintf(output_file, "%d %d\n", i, number_of_neighbor);
		vector<int> dests;
		for (auto& b : body_parts) {
			dests.push_back(b.dest_index);
		}
		//if (i < 10) {
		if (false) {
			if (i == 0) {
				printf("dests : ");
				for (auto b : body_parts) {
					printf("%d %lf ", b.dest_index, b.dest_dis);
				}
				printf("\n");
			}
			CylinderPath cylinder_path(0.01);
			cylinder_path.cntGeodesicPaths(model, i, dests);
		}
	}
	input_file.close();
	fclose(output_file);
	printf("Average_degree_before %lf\n", average_degree / node_number);
}


void svg_precompute_ich(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta, bool is_debug_mode)
{
	ElapasedTime total_t;
	double theta = asin(sqrt(eps_vg));
	theta *= const_for_theta;
	fprintf(stderr, "******** eps %.10lf const %lf theta %lf du\n", eps_vg, const_for_theta, theta / M_PI * 180.0);
	CRichModel model(input_obj_name);
	model.Preprocess();

	svg_file_name = get_DGG_filename(input_obj_name, "DGGICH", eps_vg, const_for_theta);

	printf("binary filename generated\n");
	int begin_vertex_index = 0;
	int end_vertex_index = model.GetNumOfVerts();
	int num_of_vertex = end_vertex_index - begin_vertex_index + 1;

	WxnBuffer wxn_buffer(svg_file_name);
	HeadOfSVG head_of_svg(begin_vertex_index, end_vertex_index, num_of_vertex);

	//output_file.write((char*)&head_of_svg, sizeof(head_of_svg));
	wxn_buffer.addStruct(&head_of_svg, sizeof(head_of_svg));
	ElapasedTime time_once;

	double past_time;
	double average_degree = 0;
	//if (!is_debug_mode) {
	if (true) {
		int start_vert = 0;
		int end_vert = model.GetNumOfVerts() - 1;
		int total_vert = end_vert - start_vert + 1;
		for (int source = start_vert; source <= end_vert; ++source) {
			svg_precompute_ich_vert(source, total_vert, model,
				eps_vg, theta, wxn_buffer,
				is_debug_mode, time_once,
				past_time, average_degree);
		}
	}
	else {
		std::random_device rd;
		std::mt19937 gen(0);
		std::uniform_int_distribution<> rnd_vert(0, model.GetNumOfVerts() - 1);
		int total_vert = 10;
		for (int i = 0; i < total_vert; ++i) {
			int source = rnd_vert(gen);
			printf("i %d\n", i);
			svg_precompute_ich_vert(i, total_vert, model,
				eps_vg, theta, wxn_buffer,
				is_debug_mode, time_once,
				past_time, average_degree);
		}
	}
	wxn_buffer.close();
	printf("average_degree %lf\n", average_degree / model.GetNumOfVerts());
	time_once.printTime("ich_time");
	double ich_time = time_once.getTime();
	string output_filename;
	double prune_time;
	//if (!is_debug_mode) {
	if (true) {
		JIAJUN_DGG_PRUNING::dgg_pruning(svg_file_name, eps_vg, output_filename, prune_time);
	}
	svg_file_name = output_filename;
	fprintf(stderr, "prunning time %lf\n", prune_time);
	fprintf(stderr, "total_time_and_pruning %lf\n", ich_time + prune_time);
}

void svg_precompute_hy_pruning(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta)
{
  ElapasedTime total_t;
  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  fprintf(stderr,"******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;
  CRichModel model(input_obj_name);      
  model.Preprocess();

  clock_t start = clock();
  bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(),points,faces, realIndex, originalVertNum);
  if(!success)
  {
    fprintf(stderr, "something is wrong with the input file\n" );
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  fprintf(stderr, "loading model took %.2lf secodns\n" , (double)(end - start) / (double)CLOCKS_PER_SEC);;

  svg_file_name = get_DGG_filename(input_obj_name, "DGG", eps_vg, const_for_theta);

  int begin_vertex_index = 0;
  int end_vertex_index = points.size() / 3 - 1;

  ofstream output_file (svg_file_name.c_str() , ios::out | ios::binary);
  int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
  HeadOfSVG head_of_svg(begin_vertex_index , end_vertex_index , num_of_vertex );     
  output_file.write((char*)&head_of_svg , sizeof(head_of_svg));
  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //end_vertex_index = 0;

  double t_propagate=0;
  double t_backtrace=0;
  double t_sort=0;
  double t_write = 0;

  for (int tmp_source = begin_vertex_index;tmp_source <= end_vertex_index;++tmp_source) {

    if (time_once.getTime() -  past_time > 5 ) {
      past_time = time_once.getTime();
      char buf[128];
      sprintf(buf, "Computed %.0lf percent", (double) tmp_source  * 100. / (end_vertex_index - begin_vertex_index));
      time_once.printTime(buf );
    }
    ElapasedTime dis_time;
    int source_index = tmp_source;
    vector<int> srcs;
    srcs.push_back(source_index);
    vector<geodesic::SurfacePoint> sources;
    for (unsigned i = 0; i < srcs.size(); ++i)
    {
      srcs[i] %= originalVertNum;
      srcs[i] = realIndex[srcs[i]];
      sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
    }

    ElapasedTime tm_propagate;

    geodesic::GeodesicAlgorithmBase *algorithm;

    algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh); 

    double step = 1.0;
    const double eta = 100;

    algorithm->step = step;
    algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
    map<int,double> fixedDests;
    algorithm->propagate_vg(sources, eps_vg, fixedDests);
	if (tmp_source == 0) {
		printf("dests: %d \n", fixedDests.size());
		for (auto f : fixedDests) {
			printf("%d ", f.first);
		}
		printf("\n");
	}
    t_propagate += tm_propagate.getTime();

    ElapasedTime tm_backtrace;

    vector<pair<int,double>> dests;

    struct node {
      int id;
      double dis;
      int operator<(const node & other) const{
        return dis < other.dis;
      }
    };

    //printf("tmp_source %d dest size %d\n" , tmp_source, fixedDests.size());
    vector<node> covered_points;
    covered_points.resize(fixedDests.size());
    int _index = 0;
    for (auto& d:fixedDests) {
      geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
      double dis;
      //printf("line 274\n");
      algorithm->best_source(dest_p, dis);
      //printf("line 276\n");
      covered_points[_index].id = d.first;
      covered_points[_index].dis = d.second;
      _index ++;
      //printf("v %d dis %lf\n" , d.first, d.second);
      //if (fabs(dis - d.second) > 1e-6) {
      //  fprintf(stderr,"dis %lf origin %lf\n", dis, d.second);
      //}
      dests.push_back(make_pair(d.first,d.second));
    }

    t_backtrace += tm_backtrace.getTime();
    ElapasedTime tm_sort;

    std::sort(covered_points.begin(), covered_points.end());
    std::map<int, int> mp;
    for(int i = 0; i < covered_points.size(); ++i){
      mp[covered_points[i].id] = i;
    }
    BodyHeadOfSVG body_header(source_index , dests.size());  
    output_file.write((char*)&body_header , sizeof(body_header));

    vector<BodyPartOfSVGWithK> body_parts(dests.size());
    for(int i = 0; i < dests.size(); ++i) {
      BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , mp[dests[i].first]);
      //BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , i);
      body_parts[i] = body_part;
    }
    sort(body_parts.begin() , body_parts.end());

    vector<double> angles(body_parts.size());
    for (int i = 0; i < body_parts.size(); ++i) {
      //double distance;
      //auto p(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
      //unsigned best_source = algorithm->best_source(p, distance);	
      auto dest_vert(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
      vector<geodesic::SurfacePoint> path; 
      algorithm->trace_back(dest_vert, path);
      geodesic::SurfacePoint p;
      if (path.size() > 1) {
        p = *(path.rbegin()+1);
      } else {
        p = (path[0]);
      }
 //         VERTEX,
 //   EDGE,
 //   FACE,
	//UNDEFINED_POINT
      auto& neighs = model.Neigh(source_index);
      vector<double> sum_angle(neighs.size()+1);
      sum_angle[0] = 0;
      for (int j = 1; j <= neighs.size(); ++j) {
        sum_angle[j] = sum_angle[j-1] + neighs[j-1].second;
      }

      double angle = 0;
      if ( p.type() == geodesic::VERTEX) {
        //printf("vertex\n");
        bool flag_found = false;
        for (int j = 0; j < neighs.size(); ++j) {
          auto& neigh = neighs[j];
          if ( p.base_element()->id() == model.Edge(neigh.first).indexOfRightVert) {
            //printf("yes\n");
            flag_found = true;
            angle = sum_angle[j];
            break;
          }
        }
        if (!flag_found) {
          angle = 0;
        }
      } else if(p.type() == geodesic::EDGE) {
        
        //printf("edge\n");
        int v0 = p.base_element()->adjacent_vertices()[0]->id();
        int v1 = p.base_element()->adjacent_vertices()[1]->id();
        bool flag = false;
        for (int j = 0; j < neighs.size(); ++j) {
          auto& neigh = neighs[j];
          if (v0 == model.Edge(neigh.first).indexOfRightVert) {
            int jminus1 = (j-1+neighs.size())%neighs.size();
            int vjminus1 = model.Edge(neighs[jminus1].first).indexOfRightVert;
            int jplus1 = (j+1)%neighs.size();
            int vjplus1 = model.Edge(neighs[jplus1].first).indexOfRightVert;
            //printf("v1 %d j -1 %d j + 1 %d\n" , v1, vjminus1, vjplus1); 
            CPoint3D p_cpoint3d(p.x(),p.y(),p.z());
            
            if (v1 == vjminus1) {//v1 first
              double l = model.Edge(neighs[jminus1].first).length;
              double r = (model.Vert(source_index) - p_cpoint3d).Len();
              double b = (model.Vert(vjminus1) - p_cpoint3d).Len();
              angle = sum_angle[jminus1] + acos((l * l + r * r - b * b) / (2 * l * r));
            } else if (v1 == vjplus1) {//v0 first
              double l = model.Edge(neighs[j].first).length;
              double r = (model.Vert(source_index) - p_cpoint3d).Len();
              double b = (model.Vert(v0) - p_cpoint3d).Len();
              angle = sum_angle[j] + acos((l * l + r * r - b * b) / (2 * l * r));
            }else{
              fprintf(stderr,"fuck error line 680\n");
            }
            flag = true;
            break;
          }
        }
        if (!flag) {
          fprintf(stderr,"flag %d\n" , flag);
        }

      } else{
        fprintf(stderr,"fuck error face\n");
      }
    
      angles[i] = angle;
    }

    vector<BodyPartOfSVGWithAngle> body_parts_with_angle(body_parts.size());
    for (int i = 0; i < body_parts.size(); ++i) {
      BodyPartOfSVG b = body_parts[i];
      BodyPartOfSVGWithAngle b_with_angle(b.dest_index, b.dest_dis, angles[i],0,0);
      body_parts_with_angle[i] = b_with_angle;
      //output_file.write( (char*)&b_with_angle , sizeof(b_with_angle));
    }
    sort(body_parts_with_angle.begin(), body_parts_with_angle.end());
    //for (auto& b:body_parts_with_angle) {
    //  fprintf(stderr,"b %lf " , b.angle);
    //}
    //printf("\n");
    float angle_sum = model.AngleSum(source_index);
    vector<double> tmp_angles(body_parts_with_angle.size()*2);
    for (int i = 0; i < body_parts_with_angle.size(); ++i) {
      tmp_angles[i] = body_parts_with_angle[i].angle;
    }
    for (int i = body_parts_with_angle.size(); i < tmp_angles.size(); ++i) {
      tmp_angles[i] = body_parts_with_angle[i-body_parts_with_angle.size()].angle + angle_sum;
    }
    //if (source_index == 0) {
    //  for (auto& p:tmp_angles) {
    //    fprintf(stderr,"p %lf " , p);
    //  }
    //}
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

      int start_pos = lower_bound(tmp_angles.begin(),tmp_angles.end(),start_angle) - tmp_angles.begin();
      if (start_pos > 0) start_pos--;
      int end_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), end_angle) - tmp_angles.begin();
      //if (source_index == 0 ) {
      //  fprintf(stderr,"dest %d dis %lf start_angle %lf start_pos_angle %lf end_angle  %lf end_pos_angle %lf\n" , body_parts_with_angle[i].dest_index, body_parts_with_angle[i].dest_dis,start_angle, tmp_angles[start_pos], end_angle, tmp_angles[end_pos]);
      //}
      if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
      if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
      body_parts_with_angle[i].begin_pos = start_pos;
      body_parts_with_angle[i].end_pos = end_pos;
    }

    t_sort += tm_sort.getTime();
    for (auto& b:body_parts_with_angle) {
      output_file.write((char*)&b , sizeof(b));    
    }

    delete algorithm;
    algorithm = NULL;
  }

  double t = time_once.getTime();

  output_file.close();


  string output_filename;
  double prune_time;
  JIAJUN_DGG_PRUNING::dgg_pruning(svg_file_name, eps_vg, output_filename, prune_time);

  //time_once.printTime("time_past ");
  //total_t.printTime("total_time ");
  fprintf(stderr,"propagate time %lf precent %lf%%\n" , t_propagate , t_propagate / t * 100.0);
  fprintf(stderr,"backtrace time %lf percent %lf%%\n" , t_backtrace, t_backtrace / t * 100.0);
  fprintf(stderr,"sort time %lf percent %lf%%\n" , t_sort, t_sort / t * 100.0);

  fprintf(stderr,"prunning time %lf\n" , prune_time);
  fprintf(stderr,"total_time_and_pruning %lf\n" , t + prune_time);

}


void dggPropagateCompute(const int source_index, geodesic::Mesh& mesh, double eps_vg,
	double theta, const CRichModel& model, BodyHeadOfSVG& body_header, vector<BodyPartOfSVGWithAngle>& body_parts_with_angle)
{
	vector<int> srcs;
	srcs.push_back(source_index);
	vector<geodesic::SurfacePoint> sources;
	for (unsigned i = 0; i < srcs.size(); ++i)
	{
		sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
	}

	ElapasedTime tm_propagate;

	geodesic::GeodesicAlgorithmBase *algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh);

	double step = 1.0;
	const double eta = 100;

	algorithm->step = step;
	algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
	map<int, double> fixedDests;
	algorithm->propagate_vg(sources, eps_vg, fixedDests);

	ElapasedTime tm_backtrace;

	vector<pair<int, double>> dests;

	struct node {
		int id;
		double dis;
		int operator<(const node & other) const{
			return dis < other.dis;
		}
	};

	vector<node> covered_points;
	covered_points.resize(fixedDests.size());
	int _index = 0;
	for (auto& d : fixedDests) {
		geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
		double dis;
		//printf("line 274\n");
		algorithm->best_source(dest_p, dis);
		//printf("line 276\n");
		covered_points[_index].id = d.first;
		covered_points[_index].dis = d.second;
		_index++;
		dests.push_back(make_pair(d.first, d.second));
	}

	ElapasedTime tm_sort;

	std::sort(covered_points.begin(), covered_points.end());
	std::map<int, int> mp;
	for (int i = 0; i < covered_points.size(); ++i){
		mp[covered_points[i].id] = i;
	}
	//BodyHeadOfSVG body_header(source_index, dests.size());
	body_header = BodyHeadOfSVG(source_index, dests.size());
	

	vector<BodyPartOfSVGWithK> body_parts(dests.size());
	for (int i = 0; i < dests.size(); ++i) {
		BodyPartOfSVGWithK body_part(dests[i].first, dests[i].second, mp[dests[i].first]);
		//BodyPartOfSVGWithK body_part(dests[i].first , dests[i].second , i);
		body_parts[i] = body_part;
	}
	sort(body_parts.begin(), body_parts.end());

	vector<double> angles(body_parts.size());
	for (int i = 0; i < body_parts.size(); ++i) {
		//double distance;
		//auto p(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
		//unsigned best_source = algorithm->best_source(p, distance);	
		auto dest_vert(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
		vector<geodesic::SurfacePoint> path;
		algorithm->trace_back(dest_vert, path);
		geodesic::SurfacePoint p;
		if (path.size() > 1) {
			p = *(path.rbegin() + 1);
		}
		else {
			p = (path[0]);
		}
		//         VERTEX,
		//   EDGE,
		//   FACE,
		//UNDEFINED_POINT
		auto& neighs = model.Neigh(source_index);
		vector<double> sum_angle(neighs.size() + 1);
		sum_angle[0] = 0;
		for (int j = 1; j <= neighs.size(); ++j) {
			sum_angle[j] = sum_angle[j - 1] + neighs[j - 1].second;
		}

		double angle = 0;
		if (p.type() == geodesic::VERTEX) {
			//printf("vertex\n");
			bool flag_found = false;
			for (int j = 0; j < neighs.size(); ++j) {
				auto& neigh = neighs[j];
				if (p.base_element()->id() == model.Edge(neigh.first).indexOfRightVert) {
					//printf("yes\n");
					flag_found = true;
					angle = sum_angle[j];
					break;
				}
			}
			if (!flag_found) {
				angle = 0;
				//printf("vertex %d source %d\n" , p.base_element()->id(), source_index);
			}

		}
		else if (p.type() == geodesic::EDGE) {

			//printf("edge\n");
			int v0 = p.base_element()->adjacent_vertices()[0]->id();
			int v1 = p.base_element()->adjacent_vertices()[1]->id();
			bool flag = false;
			for (int j = 0; j < neighs.size(); ++j) {
				auto& neigh = neighs[j];
				if (v0 == model.Edge(neigh.first).indexOfRightVert) {
					int jminus1 = (j - 1 + neighs.size()) % neighs.size();
					int vjminus1 = model.Edge(neighs[jminus1].first).indexOfRightVert;
					int jplus1 = (j + 1) % neighs.size();
					int vjplus1 = model.Edge(neighs[jplus1].first).indexOfRightVert;
					//printf("v1 %d j -1 %d j + 1 %d\n" , v1, vjminus1, vjplus1); 
					CPoint3D p_cpoint3d(p.x(), p.y(), p.z());

					if (v1 == vjminus1) {//v1 first
						double l = model.Edge(neighs[jminus1].first).length;
						double r = (model.Vert(source_index) - p_cpoint3d).Len();
						double b = (model.Vert(vjminus1) - p_cpoint3d).Len();
						angle = sum_angle[jminus1] + acos((l * l + r * r - b * b) / (2 * l * r));
					}
					else if (v1 == vjplus1) {//v0 first
						double l = model.Edge(neighs[j].first).length;
						double r = (model.Vert(source_index) - p_cpoint3d).Len();
						double b = (model.Vert(v0) - p_cpoint3d).Len();
						angle = sum_angle[j] + acos((l * l + r * r - b * b) / (2 * l * r));
					}
					else{
						fprintf(stderr, "fuck error line 680\n");
					}
					flag = true;
					break;
				}
			}
			if (!flag) {
				fprintf(stderr, "flag %d\n", flag);
			}

		}
		else{
			fprintf(stderr, "fuck error face\n");
		}

		angles[i] = angle;
	}

	//vector<BodyPartOfSVGWithAngle> body_parts_with_angle(body_parts.size());
	body_parts_with_angle.resize(body_parts.size());
	for (int i = 0; i < body_parts.size(); ++i) {
		BodyPartOfSVG b = body_parts[i];
		BodyPartOfSVGWithAngle b_with_angle(b.dest_index, b.dest_dis, angles[i], 0, 0);
		body_parts_with_angle[i] = b_with_angle;
	}
	sort(body_parts_with_angle.begin(), body_parts_with_angle.end());
	//for (auto& b:body_parts_with_angle) {
	//  fprintf(stderr,"b %lf " , b.angle);
	//}
	//printf("\n");
	float angle_sum = model.AngleSum(source_index);
	vector<double> tmp_angles(body_parts_with_angle.size() * 2);
	for (int i = 0; i < body_parts_with_angle.size(); ++i) {
		tmp_angles[i] = body_parts_with_angle[i].angle;
	}
	for (int i = body_parts_with_angle.size(); i < tmp_angles.size(); ++i) {
		tmp_angles[i] = body_parts_with_angle[i - body_parts_with_angle.size()].angle + angle_sum;
	}
	//if (source_index == 0) {
	//  for (auto& p:tmp_angles) {
	//    fprintf(stderr,"p %lf " , p);
	//  }
	//}
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
		//if (source_index == 0 ) {
		//  fprintf(stderr,"dest %d dis %lf start_angle %lf start_pos_angle %lf end_angle  %lf end_pos_angle %lf\n" , body_parts_with_angle[i].dest_index, body_parts_with_angle[i].dest_dis,start_angle, tmp_angles[start_pos], end_angle, tmp_angles[end_pos]);
		//}
		if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
		if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
		body_parts_with_angle[i].begin_pos = start_pos;
		body_parts_with_angle[i].end_pos = end_pos;
	}

	delete algorithm;
	algorithm = NULL;




}

void dggPropagate(const int source_index, geodesic::Mesh& mesh, double eps_vg,
                  double theta, FILE* output_file, const CRichModel& model)
{
	BodyHeadOfSVG body_header;
	vector<BodyPartOfSVGWithAngle> body_parts_with_angle;
	dggPropagateCompute(source_index, mesh, eps_vg, theta, model, body_header, body_parts_with_angle);
	fwrite(&body_header, sizeof(body_header), 1, output_file);
	for (auto& b : body_parts_with_angle) {
		fwrite(&b, sizeof(b), 1, output_file);
	}
}


void dggPropagateHead(const HeadOfSVG& head, const string& part_svg_filename, 
											geodesic::Mesh& mesh, double eps_vg,
											double theta, const CRichModel& model)
{
	//ofstream output_file (part_svg_filename.c_str() , ios::out | ios::binary);
	FILE* output_file = fopen(part_svg_filename.c_str(), "wb");
	//WxnBuffer buffer(output_file);
	//output_file.write((char*)&head , sizeof(head));
	fwrite(&head, sizeof(head), 1, output_file);
	//buffer.addStruct((void*)&head, sizeof(head));
	ElapasedTime time_total;
	double last_t;
	for (int i = head.begin_vertex_index; i <= head.end_vertex_index; ++i) {
		if (time_total.getTime() - last_t > 5) {
			//time_total.printTime("time");
			double percent = (i - head.begin_vertex_index) / (double)(head.end_vertex_index - head.begin_vertex_index) * 100.0;
			last_t = time_total.getTime();
			double remain_hours = last_t / percent * (100.0-percent) / 3600.0;
			printf("current %.2lf percent, time %lf s, remain %.2lf hours\n", percent, last_t, remain_hours);
		}
		dggPropagate(i, mesh, eps_vg, theta, output_file, model);
	}
	time_total.printTime("part time");
	//output_file.close();
	fclose(output_file);
}

void combinePartPrecomputeFiles(const vector<string>& part_filenames,const string& svg_filename, int num_of_vertex, int thread_num)
{
	HeadOfSVG head;
	head.begin_vertex_index = 0; head.end_vertex_index = num_of_vertex - 1;
	head.num_of_vertex = num_of_vertex;
	ofstream output_file (svg_filename.c_str() , ios::out | ios::binary);
	output_file.write((char*)&head , sizeof(head));

	for (int thread_id = 0; thread_id < thread_num; ++thread_id) {
		const string& part_svg_filename = part_filenames[thread_id];
		std::ifstream input_file (part_svg_filename, std::ios::in | std::ios::binary);
		HeadOfSVG head_of_svg;
		input_file.read( (char*)&head_of_svg , sizeof(head_of_svg));
		//		printf("head %d %d vert_sz %d\n" , head_of_svg.begin_vertex_index , head_of_svg.end_vertex_index, head_of_svg.num_of_vertex);

		for (int i = head_of_svg.begin_vertex_index; i <= head_of_svg.end_vertex_index;++i) {
			//printf("i%d\n" , i);
			BodyHeadOfSVG body_head;
			input_file.read((char*)&body_head , sizeof(body_head));
			//printf("readed head\n");
			vector<BodyPartOfSVGWithAngle> body_parts;
			//printf("neigh %d\n" , body_head.neighbor_num);
			body_parts.reserve(body_head.neighbor_num);
			//printf("neigh %d\n" , body_head.neighbor_num);
			for(int j = 0; j < body_head.neighbor_num;++j) {
				//printf("j%d\n" , j);
				BodyPartOfSVGWithAngle body_part;
				input_file.read((char*)&body_part , sizeof(body_part));
				body_parts.push_back(body_part);
			}
			output_file.write((char*)&body_head , sizeof(body_head));
			for (auto& b:body_parts) {
				output_file.write((char*)&b , sizeof(b));
			}
		}
		input_file.close();
	}
	output_file.close();
}


void ichPropogateHead(const HeadOfSVG& head, const string& part_svg_filename, double eps_vg, double theta, const CRichModel& model) 
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
		getDisAndAngles(model, source, eps_vg, dests, angles, dis);

		vector<BodyPartOfSVGWithAngle> body_parts_with_angle;
		getFanOutput(dests, angles, dis, model, theta, source, body_parts_with_angle);

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

const double JiajunMaxDist(1e5);

template <class T>
void dijkstra_pruning(const vector<vector<int>>&  graph_neighbor,
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
				if (src == 0 ) {
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
		//if (i > 0 && i % (head_of_svg.num_of_vertex / 10) == 0){
		//	std::cerr << "read " << i * 100 / head_of_svg.num_of_vertex << " percent \n";
		//}
	}
	input_file.close();
}


void test(vector<int>& a, int thread_id)
{
	a[thread_id] = thread_id;
	
}

template<class T>
void dijkstraPruningThread(int thread_id, int thread_num, int node_number, 
						   const vector<vector<int>>& graph_neighbor,
						   const vector<vector<T>>& graph_neighbor_dis,
						   vector<vector<bool>>& graph_neighbor_deleted, double eps_vg)
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

	for (int i = begin; i <= end; ++i) {
		dijkstra_pruning<T>(graph_neighbor, 
		graph_neighbor_dis, graph_neighbor_deleted[i], i,
		eps_vg, dis, mark);
	}
}

template<class T> 
void wxn_pruning(const string& svg_file_name, double eps_vg, string& test_output_filename, int thread_num, double& prune_time)
{
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
			tt[thread_id] = std::thread(&dijkstraPruningThread<T>, thread_id,
										thread_num, node_number, ref(graph_neighbor),
										ref(graph_neighbor_dis), ref(graph_neighbor_deleted),
										eps_vg);

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
				} else {
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


void svg_precompute_ich_multithread_before_pruning(const string& input_obj_name, double eps_vg,
												   string& svg_file_name, double const_for_theta, 
												   int thread_num, double& ich_multi_time)
{

	char buf[1024];
	sprintf(buf, "%s_DGGICH%.10lf_c%.0lf.binary", input_obj_name.substr(0, input_obj_name.length() - 4).c_str(), eps_vg, const_for_theta);
	svg_file_name = string(buf);
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
		tt[i] = std::thread(&ichPropogateHead, heads[i], svg_part_file_names[i], eps_vg, theta, std::ref(model));
	}
	for (int i = 0; i < thread_num; ++i) {
		tt[i].join();
	}
	ich_multi_time = time_multi.getTime();
	time_multi.printTime("ich_multi_time");

	ElapasedTime combine_time;
	combinePartPrecomputeFiles(svg_part_file_names, svg_file_name, model.GetNumOfVerts(), thread_num);
	combine_time.printTime("combine time");

}


void svg_precompute_ich_multithread(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta, int thread_num)
{
	double ich_multi_time = 0;
	svg_precompute_ich_multithread_before_pruning(input_obj_name, eps_vg, svg_file_name, const_for_theta, thread_num, ich_multi_time);
	string wxn_output_filename(svg_file_name.substr(0, svg_file_name.length() - 7) + "_pruning.binary");
	double prune_time;
	ElapasedTime wxn_pruning_time;
	wxn_pruning<double>(svg_file_name, eps_vg, wxn_output_filename, thread_num, prune_time);
	wxn_pruning_time.printTime("wxn pruning time");
	fprintf(stderr, "prunning time %lf\n", prune_time);
	fprintf(stderr, "total_time_and_pruning %lf\n", ich_multi_time + prune_time);
}


void  svg_precompute_hy_multithread(const string& input_obj_name, double eps_vg, string& output_filename, double const_for_theta, int thread_num)
{
  ElapasedTime total_t;
  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  fprintf(stderr,"******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;
  CRichModel model(input_obj_name);
  model.Preprocess();

  clock_t start = clock();
  bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(),points,faces, realIndex, originalVertNum);
  if(!success)
  {
    fprintf(stderr, "something is wrong with the input file\n" );
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  fprintf(stderr, "loading model took %.2lf secodns\n" , (double)(end - start) / (double)CLOCKS_PER_SEC);;


	//seperate vertices to parts
	vector<HeadOfSVG> heads;
	vector<string> svg_part_file_names;
	int part_size = model.GetNumOfVerts() / thread_num;
	for (int i = 0; i < thread_num; ++i) {
		HeadOfSVG head;
		head.num_of_vertex = model.GetNumOfVerts();
		head.begin_vertex_index =  i * part_size;
		if (i != thread_num - 1) {
			head.end_vertex_index = (i+1)*part_size - 1;
		}else{
			head.end_vertex_index = model.GetNumOfVerts() - 1;
		}
		heads.push_back(head);
		printf("head %d %d vert_sz %d\n" , head.begin_vertex_index , head.end_vertex_index, head.num_of_vertex);
		char buf[1024];
		sprintf(buf, "%s_DGG%lf_c%.0lf_part%d.binary", input_obj_name.substr(0,input_obj_name.length() - 4 ).c_str(), eps_vg, const_for_theta,i);
		svg_part_file_names.push_back((string)buf);
	}
	
	std::thread *tt = new std::thread[thread_num];

	for (int thread_id = 0; thread_id < thread_num; ++thread_id) {
		tt[thread_id] = std::thread(&dggPropagateHead, std::ref(heads[thread_id]), std::ref(svg_part_file_names[thread_id]), 
														std::ref(mesh), eps_vg, theta, std::ref(model));

		//void dggPropagateHead(const HeadOfSVG& head, const string& part_svg_filename, 
		//									geodesic::Mesh& mesh, double eps_vg,
		//									double theta, const CRichModel& model)

		//dggPropagateHead(heads[thread_id], svg_part_file_names[thread_id], 
		//												mesh, eps_vg, theta, model);
	}
	for (int i = 0; i < thread_num; ++i){
		tt[i].join();
	}
		delete[] tt;
	//combine svg_files
	char buf[1024];
	sprintf(buf,"%s_DGG%lf_c%.0lf.binary", input_obj_name.substr(0,input_obj_name.length() - 4 ).c_str(), eps_vg, const_for_theta);
	string svg_file_name = string(buf);
	ElapasedTime combine_time;
	combinePartPrecomputeFiles(svg_part_file_names, svg_file_name, model.GetNumOfVerts(), thread_num);
	combine_time.printTime("combine time");

  double prune_time;
  JIAJUN_DGG_PRUNING::dgg_pruning(svg_file_name, eps_vg, output_filename, prune_time);

  fprintf(stderr,"prunning time %lf\n" , prune_time);
  fprintf(stderr,"total time %lf\n" , total_t.getTime());

}

  void svg_precompute_hy_new(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta)
{
  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  fprintf(stderr,"******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;
  CRichModel model(input_obj_name);
  model.Preprocess();

  clock_t start = clock();
  bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(),points,faces, realIndex, originalVertNum);
  if(!success)
  {
    fprintf(stderr,"something is wrong with the input file\n");;
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  fprintf(stderr, "loading model took %lf seconds" , (double)(end - start) / (double)CLOCKS_PER_SEC );

  char buf[1024];
  sprintf(buf,"%s_DGG%lf_c%.0lf.binary", input_obj_name.substr(0,input_obj_name.length() - 4 ).c_str(), eps_vg, const_for_theta);
  svg_file_name = string(buf);
  //svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
  //  + "_HY" + to_string(eps_vg) +  ".binary";


  int num_of_vertex = points.size() / 3;
  ofstream output_file (svg_file_name.c_str() , ios::out | ios::binary);
  HeadOfSVG head_of_svg(0 , num_of_vertex - 1 , num_of_vertex );     
  output_file.write((char*)&head_of_svg , sizeof(head_of_svg));
  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //end_vertex_index = 0;

  vector<vector<SVGEdgeWithAngle>> graph_with_angle(num_of_vertex);
  for (int tmp_source = 0;tmp_source < num_of_vertex;++tmp_source) {
    if (time_once.getTime() -  past_time > 5 ) {
      past_time = time_once.getTime();
      char buf[128];
      sprintf(buf, "Computed %.0lf percent", (double) tmp_source  * 100. / num_of_vertex);
      time_once.printTime(buf );
    }
    ElapasedTime dis_time;
    int source_index = tmp_source;
    vector<int> srcs;
    srcs.push_back(source_index);
    vector<geodesic::SurfacePoint> sources;
    for (unsigned i = 0; i < srcs.size(); ++i)
    {
      srcs[i] %= originalVertNum;
      srcs[i] = realIndex[srcs[i]];
      sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
    }

    geodesic::GeodesicAlgorithmBase *algorithm;

    algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh); 

    double step = 1.0;
    const double eta = 100;

    algorithm->model_ptr_ = &model;
    algorithm->step = step;
    algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
    map<int,double> fixedDests;
    algorithm->propagate_vg(sources, eps_vg, fixedDests);

    for (auto& d:fixedDests) {//d.first is dest vertex
      geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
      double dis;
      double angle;
      algorithm->best_source_with_angle(dest_p, dis,angle);
      if (d.first != source_index) {
        graph_with_angle[d.first].push_back(SVGEdgeWithAngle(source_index,dis,angle));      
      }
      //if (source_index == 10) {
      //  fprintf(stderr,"source %d dest %d dis %lf angle %lf\n" , source_index, d.first, d.second, angle);
      //}
      //if (d.first == 10) {
      //  fprintf(stderr,"source %d dest %d dis %lf angle %lf\n" , source_index, d.first, d.second, angle);
      //}
    }
    delete algorithm;
    algorithm = NULL;
  }

  for (int source_index = 0; source_index < num_of_vertex; ++source_index) {

    auto& graph_body_part = graph_with_angle[source_index];
    sort(graph_body_part.begin(), graph_body_part.end());
    int neigh_size = graph_body_part.size();
    vector<BodyPartOfSVGWithRange> body_parts_with_range(neigh_size);
    for (int i = 0; i < neigh_size; ++i) {
      body_parts_with_range[i] = BodyPartOfSVGWithRange(graph_body_part[i].dest_index,graph_body_part[i].dest_dis,-1,-1);
    }

    float angle_sum = model.AngleSum(source_index);
    vector<double> tmp_angles(neigh_size*2);
    for (int i = 0; i < neigh_size; ++i) {
      tmp_angles[i] = graph_body_part[i].angle;
    }
    for (int i = neigh_size; i < tmp_angles.size(); ++i) {
      tmp_angles[i] = graph_body_part[i-neigh_size].angle + angle_sum;
    }
    for (int i = 0; i < graph_body_part.size(); ++i) {//assume i is father
      float father_angle = graph_body_part[i].angle;
      //based on father_angle as 0
      float start_angle = M_PI - theta + father_angle;
      float end_angle = angle_sum - (M_PI - theta) + father_angle;
      if (start_angle > end_angle) {
        body_parts_with_range[i].begin_pos = -1;
        body_parts_with_range[i].end_pos = -1;
        continue;
      }

      int start_pos = lower_bound(tmp_angles.begin(),tmp_angles.end(),start_angle) - tmp_angles.begin();
      if (start_pos > 0) start_pos--;
      int end_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), end_angle) - tmp_angles.begin();
      //if (source_index == 0) {
      //  fprintf(stderr,"dest %d dis %lf start_angle %lf start_pos_angle %lf end_angle  %lf end_pos_angle %lf\n" ,body_parts_with_range[i].dest_index , body_parts_with_range[i].dest_dis,  start_angle, tmp_angles[start_pos], end_angle, tmp_angles[end_pos]);
      //}
      if (start_pos >= neigh_size) start_pos -= neigh_size;
      if (end_pos >= neigh_size) end_pos -= neigh_size;
      body_parts_with_range[i].begin_pos = start_pos;
      body_parts_with_range[i].end_pos = end_pos;
    }
    BodyHeadOfSVG body_header(source_index , neigh_size);  
    output_file.write((char*)&body_header , sizeof(body_header));

    for (auto& b:body_parts_with_range) {
      output_file.write((char*)&b , sizeof(b));    
    }

  }
  double t = time_once.getTime();
  output_file.close();
}

#if 0
void svg_precompute_hy_fast(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta)
{
  ElapasedTime t;
  svg_precompute_hy_new(input_obj_name,eps_vg,svg_file_name, const_for_theta);
  string output_filename;
  double prune_time;
  JIAJUN_DGG_PRUNING::dgg_pruning_new(svg_file_name, eps_vg, output_filename, prune_time);

  t.printTime("time_past ");

 // fprintf(stderr,"total_time_and_pruning %lf\n" , t + prune_time);
}
#endif

void findPointUsingPathTracer(const CRichModel& model, const int source_vert, 
							  const int direction_edge,	const double dis,
							  const double angle, YXPathTracer& path_tracer,
							  CPoint3D& p, int& face_idx)
{
	const auto& e = model.Edge(direction_edge);
	int v1 = e.indexOfLeftVert;
	int v2 = e.indexOfRightVert;
	int e_id_yx = path_tracer.metric.getEdgeFrom2Verts(v1, v2);
	YXPath tmp_path;
	path_tracer.computeSinglePathArbitraryStart(e_id_yx, angle, dis, tmp_path, source_vert);
	//center_3d = CenterPoint(path_tracer.convertToYXPoint3D(tmp_path.back()).toCPoint3D());
	p = path_tracer.convertToYXPoint3D(tmp_path.back()).toCPoint3D();
	face_idx = path_tracer.getFaceIndex(tmp_path.back());
}


void getPointOnEdge(const vector<double>& sum_angle, double angle, const CRichModel& model,
					const vector<pair<int, double>>& neighs , int source_vert,
					CPoint3D& p, int& v0_id, int& v1_id)
{
	auto it = lower_bound(sum_angle.begin(), sum_angle.end(), angle);
	it = prev(it);
	int pos_in_neigh = it - sum_angle.begin();
	//vert pos_in_neigh and pos_in_neigh + 1
	v0_id = model.Edge(neighs[pos_in_neigh].first).indexOfRightVert;
	v1_id = model.Edge(neighs[(pos_in_neigh + 1) % neighs.size()].first).indexOfRightVert;
	auto& v0 = model.Vert(v0_id);
	auto& v1 = model.Vert(v1_id);
	auto& o = model.Vert(source_vert);
	double l = (v0 - v1).Len();
	double r = (v0 - o).Len();
	double b = (v1 - o).Len();
	double alpha = acos((l * l + r * r - b * b) / (2 * l * r));
	double theta = angle - *it;
	double l2 = r / sin(M_PI - theta - alpha) * sin(theta);
	p = v0 * (1 - l2 / l) + v1 * l2 / l;
}

//auto it = lower_bound(sum_angle.begin(), sum_angle.end(), divide_angles[k]);
//it = prev(it);
//int pos_in_neigh = it - sum_angle.begin();
////vert pos_in_neigh and pos_in_neigh + 1
//int v0_id = model.Edge(neighs[pos_in_neigh].first).indexOfRightVert;
//int v1_id = model.Edge(neighs[(pos_in_neigh + 1) % neighs.size()].first).indexOfRightVert;
//auto& v0 = model.Vert(v0_id);
//auto& v1 = model.Vert(v1_id);
//auto& o = model.Vert(i);
//double l = (v0 - v1).Len();
//double r = (v0 - o).Len();
//double b = (v1 - o).Len();
//double alpha = acos((l * l + r * r - b * b) / (2 * l * r));
//double theta = divide_angles[k] - *it;
//double l2 = r / sin(M_PI - theta - alpha) * sin(theta);
//CPoint3D p = v0 * (1 - l2 / l) + v1 * l2 / l;


void cnt_percent(const CRichModel& model, SparseGraph<float>* s_graph)
{
	int cnt = 0;
	int percent[20] = { 0 };
	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		auto& angles = s_graph->graph_neighbor_angle[i];
		vector<double> angles_diff(angles.size());
		bool flag = false;
		for (int j = 0; j < angles.size() - 1; ++j) {
			angles_diff[j] = angles[j + 1] - angles[j];
			//printf("a%lf ", angles[j]);
			if (angles_diff[j] >= 2 * 2 * M_PI / angles.size()) {
				//printf("%lf ", angles_diff[j] / (2 * M_PI / angles.size()));
				int pos_in_percent = int(angles_diff[j] / (2 * M_PI / angles.size()));
				percent[pos_in_percent]++;
				flag = true;
				cnt++;
			}
		}
		//printf("\n");
		if (flag) {
			//cnt++;
			//printf("\n");
		}
	}
	printf("cnt %d total %d \n", cnt, model.GetNumOfVerts());
	for (int i = 0; i <= 10; ++i) {
		printf("%9d ", i);
	}
	printf("\n");
	for (int i = 0; i <= 10; ++i) {
		printf("%.7lf ", percent[i] / (double)cnt);
	}
	printf("\n");
}

void svg_precompute_LiuYongjin_fixing(const string& input_file_name, double eps_vg, const string& svg_file_name)
{
	SparseGraph<float>* s_graph = NULL;
	s_graph = new LC_HY<float>();
	s_graph->read_svg_file_with_angle((string)svg_file_name);
	CRichModel model(input_file_name);
	model.Preprocess();

	YXPathTracer path_tracer;
	printf("model_file_name %s\n", input_file_name.c_str());
	path_tracer.init(input_file_name.c_str());

	dynamic_cast<LC_HY<float>*>(s_graph)->setModel(model);


	vector<int> origin_of_new_faces;
	ModelForSubdivide sub_model(model, origin_of_new_faces);
	vector<CPoint3D> added_points_list;
	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		auto& geo_dises = s_graph->graphNeighborDis(i);
		auto& angles = s_graph->graph_neighbor_angle[i];
		auto& neighs = model.Neigh(i);
		vector<double> sum_angle(neighs.size() + 1);
		sum_angle[0] = 0;
		for (int j = 1; j < sum_angle.size(); ++j) {
			sum_angle[j] = sum_angle[j - 1] + neighs[j - 1].second;
		}
		
		for (int j = 0; j < angles.size(); ++j) {
			double angles_diff{ 0 };
			double dis{ 0 };
			if (j != angles.size() - 1) {
				angles_diff = angles[j + 1] - angles[j];
				dis = (geo_dises[j] + geo_dises[j + 1]) / 2;
			} else {
				angles_diff = sum_angle.back() - angles[j];
				dis = (geo_dises[j] + geo_dises[0]) / 2;
			}
			int pos_in_percent = int(angles_diff / (2 * M_PI / angles.size()));
			if (pos_in_percent >= 2) {
				vector<double> divide_angles(pos_in_percent - 1);
				for (int k = 0; k < divide_angles.size(); ++k) {
					divide_angles[k] = angles[j] + (double)(k + 1) / (double)pos_in_percent * angles_diff;
					CPoint3D p;
					int face_index;
					
					findPointUsingPathTracer(model, i, neighs[0].first, dis,
											 divide_angles[k], path_tracer, p, face_index);
					

					//printf("face_index %d\n", face_index);

					sub_model.addVertexOnFace(p, face_index);
					added_points_list.push_back(p);
				}
			}
		}
	}

	printBallToObj(added_points_list, "added_points_list.obj", 0.001);


	sub_model.subdivide();
	sub_model.WriteToFile("fertility_nf10k_ani2_lyj.obj");


}