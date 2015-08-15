#include "stdafx.h"
#include "wxnTime.h"
#include "svg_precompute.h"
#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"

#include "MMP\geodesic_algorithm_vg_mmp.h"
#include "JIAJUN\dgg_pruning.h"
#include "svg_definition.h"
#include <thread>

bool flag_first_output = true;
struct WxnBuffer{
	char* buf;
	int len;
	int capacity;
	FILE* file;
	WxnBuffer(FILE* _file) :file(_file) {
		capacity = 16 * 1024 * 1024;
		buf = new char[capacity];
		len = 0;
	}
	void addStruct(const void* ptr, int struct_size)
	{
		//fwrite(ptr, struct_size, 1, file);
		//return;
		if (len + struct_size > capacity) {
			fwrite(buf, sizeof(char), len, file);
			len = 0;
		}

		if (struct_size > capacity) {
			fprintf(stderr, "str too large!");
			exit(1);
		}
		memcpy((void*)(buf + len), ptr, struct_size);
		len += struct_size;
	}
	void addText(const char* str, int str_len) {
		if (len + str_len > capacity) {
			fwrite(buf, sizeof(char), len, file);
			len = 0;
		}
		if (str_len > capacity) {
			fprintf(stderr, "str too large!");
			exit(1);
		}
		for (int i = 0; i < str_len; ++i) {
			buf[len++] = str[i];
		}
	}

	~WxnBuffer(){
		delete[] buf;
	}
};


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


//#pragma pack(4)
struct SVGEdge{
	int v;
//	bool deleted;
	double dis;
};



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
    + "_SVGMMP_k" + to_string(fixed_k) +  ".binary";

  int begin_vertex_index = 0;

  int end_vertex_index = points.size() / 3 - 1;

  ofstream output_file (svg_file_name.c_str() , ios::out | ios::binary);
  int num_of_vertex = end_vertex_index - begin_vertex_index + 1;
  HeadOfSVG head_of_svg(begin_vertex_index , end_vertex_index , num_of_vertex );     
  output_file.write((char*)&head_of_svg , sizeof(head_of_svg));
  ElapasedTime time_once;

  double dis_time_total(0);
  double past_time(0);
  //end_vertex_index = 10;

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
    algorithm->propagate_local(sources, fixed_k, fixedDests);



    vector<pair<int,double>> dests;

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
    for (auto& d:fixedDests) {
      geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
      double dis;
      //printf("line 274\n");
      algorithm->best_source(dest_p, dis);
      //printf("line 276\n");
      covered_points[_index].id = d.first;
      covered_points[_index].dis = d.second;
      _index ++;
      if (fabs(dis - d.second) > 1e-6) {
        fprintf(stderr,"dis %lf origin %lf\n", dis, d.second);
      }
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
    
    

    for (int i = 0; i < body_parts.size(); ++i) {
      BodyPartOfSVG b = body_parts[i];
      output_file.write( (char*)&b , sizeof(b));
    }
  }

  time_once.printTime("time past ");
  output_file.close();

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

  char buf[1024];
  sprintf(buf,"%s_HY%lf_c%.0lf.binary", input_obj_name.substr(0,input_obj_name.length() - 4 ).c_str(), eps_vg, const_for_theta);
  svg_file_name = string(buf);
  //svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
  //  + "_HY" + to_string(eps_vg) +  ".binary";

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

  char buf[1024];
  sprintf(buf,"%s_DGG%lf_c%.0lf.binary", input_obj_name.substr(0,input_obj_name.length() - 4 ).c_str(), eps_vg, const_for_theta);
  svg_file_name = string(buf);
  //svg_file_name = input_obj_name.substr(0,input_obj_name.length() - 4 ) 
  //  + "_HY" + to_string(eps_vg) +  ".binary";

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


void dggPropagate(const int source_index, geodesic::Mesh& mesh, double eps_vg,
                  double theta, WxnBuffer& buffer, const CRichModel& model)
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
    map<int,double> fixedDests;
    algorithm->propagate_vg(sources, eps_vg, fixedDests);

    ElapasedTime tm_backtrace;

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
    for (auto& d:fixedDests) {
      geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
      double dis;
      //printf("line 274\n");
      algorithm->best_source(dest_p, dis);
      //printf("line 276\n");
      covered_points[_index].id = d.first;
      covered_points[_index].dis = d.second;
      _index ++;
      dests.push_back(make_pair(d.first,d.second));
    }

    ElapasedTime tm_sort;

    std::sort(covered_points.begin(), covered_points.end());
    std::map<int, int> mp;
    for(int i = 0; i < covered_points.size(); ++i){
      mp[covered_points[i].id] = i;
    }
    BodyHeadOfSVG body_header(source_index , dests.size());  
    //output_file.write((char*)&body_header , sizeof(body_header));
	//fwrite(&body_header, sizeof(body_header), 1, output_file);
	buffer.addStruct(&body_header, sizeof(body_header));

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
      //if (source_index == 0 ) {
      //  fprintf(stderr,"dest %d dis %lf start_angle %lf start_pos_angle %lf end_angle  %lf end_pos_angle %lf\n" , body_parts_with_angle[i].dest_index, body_parts_with_angle[i].dest_dis,start_angle, tmp_angles[start_pos], end_angle, tmp_angles[end_pos]);
      //}
      if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
      if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
      body_parts_with_angle[i].begin_pos = start_pos;
      body_parts_with_angle[i].end_pos = end_pos;
    }

    for (auto& b:body_parts_with_angle) {
      //output_file.write((char*)&b , sizeof(b));    
		//fwrite(&b, sizeof(b), 1, output_file);
		buffer.addStruct(&b, sizeof(b));
	}
    delete algorithm;
    algorithm = NULL;
}


void dggPropagateHead(const HeadOfSVG& head, const string& part_svg_filename, 
											geodesic::Mesh& mesh, double eps_vg,
											double theta, const CRichModel& model)
{
	//ofstream output_file (part_svg_filename.c_str() , ios::out | ios::binary);
	FILE* output_file = fopen(part_svg_filename.c_str(), "wb");
	WxnBuffer buffer(output_file);
	//output_file.write((char*)&head , sizeof(head));
	//fwrite(&head, sizeof(head), 1, output_file);
	buffer.addStruct((void*)&head, sizeof(head));
	ElapasedTime time_total;
	double last_t;
	for (int i = head.begin_vertex_index; i <= head.end_vertex_index; ++i) {
		if (time_total.getTime() - last_t > 5) {
			time_total.printTime("time");
			last_t = time_total.getTime();
		}
		dggPropagate(i, mesh, eps_vg, theta, buffer, model);
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
														mesh, eps_vg, theta, std::ref(model));

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

    //vector<BodyPartOfSVGWithAngle> body_parts_with_angle(body_parts.size());
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