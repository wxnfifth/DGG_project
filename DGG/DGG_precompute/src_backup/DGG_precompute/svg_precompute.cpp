#include "stdafx.h"
#include "wxnTime.h"
#include "svg_precompute.h"
#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"

#include "MMP\geodesic_algorithm_vg_mmp.h"
#include "JIAJUN\dgg_pruning.h"
#include "svg_definition.h"

bool flag_first_output = true;
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
    cout << "something is wrong with the input file" << endl;
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  cout << "loading model took " << (double)(end - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;

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
      //  printf("dis %lf origin %lf\n", dis, d.second);
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
    cout << "something is wrong with the input file" << endl;
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  cout << "loading model took " << (double)(end - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;

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
      //  printf("dis %lf origin %lf\n", dis, d.second);
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
		cout << "something is wrong with the input file" << endl;
		return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
	cout << "loading model took " << (double)(end - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;

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
        printf("dis %lf origin %lf\n", dis, d.second);
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
  printf("******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

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
    cout << "something is wrong with the input file" << endl;
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  cout << "loading model took " << (double)(end - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;

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
      //  printf("dis %lf origin %lf\n", dis, d.second);
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
              printf("fuck error line 680\n");
            }
            flag = true;
            break;
          }
        }
        if (!flag) {
          printf("flag %d\n" , flag);
        }

      } else{
        printf("fuck error face\n");
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
    //  printf("b %lf " , b.angle);
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
    //    printf("p %lf " , p);
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

  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  printf("******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

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
    cout << "something is wrong with the input file" << endl;
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  cout << "loading model took " << (double)(end - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;

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
      //  printf("dis %lf origin %lf\n", dis, d.second);
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
              printf("fuck error line 680\n");
            }
            flag = true;
            break;
          }
        }
        if (!flag) {
          printf("flag %d\n" , flag);
        }

      } else{
        printf("fuck error face\n");
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
    //  printf("b %lf " , b.angle);
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
    //    printf("p %lf " , p);
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

  double t = time_once.getTime();

  output_file.close();


  string output_filename;
  double prune_time;
  JIAJUN_DGG_PRUNING::dgg_pruning(svg_file_name, eps_vg, output_filename, prune_time);

  time_once.printTime("time_past ");

  printf("total_time_and_pruning %lf\n" , t + prune_time);

}


void svg_precompute_hy_fast(const string& input_obj_name, double eps_vg, string& svg_file_name, double const_for_theta)
{

  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  printf("******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

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
    cout << "something is wrong with the input file" << endl;
    return;
  }
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);		//create internal
  clock_t end = clock();
  cout << "loading model took " << (double)(end - start) / (double)CLOCKS_PER_SEC << " seconds." << endl;

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
      //  printf("dis %lf origin %lf\n", dis, d.second);
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
              printf("fuck error line 680\n");
            }
            flag = true;
            break;
          }
        }
        if (!flag) {
          printf("flag %d\n" , flag);
        }

      } else{
        printf("fuck error face\n");
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
    //  printf("b %lf " , b.angle);
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
    //    printf("p %lf " , p);
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

  double t = time_once.getTime();

  output_file.close();


  string output_filename;
  double prune_time;
  JIAJUN_DGG_PRUNING::dgg_pruning(svg_file_name, eps_vg, output_filename, prune_time);

  time_once.printTime("time_past ");

  printf("total_time_and_pruning %lf\n" , t + prune_time);





}