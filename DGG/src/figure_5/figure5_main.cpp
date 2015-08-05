#include "stdafx.h"
#include "WXN\wxn_path_helper.h"
#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"
#include "MMP\geodesic_algorithm_vg_mmp.h"

void findDest(CRichModel& model, int source, int dest)
{
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

  geodesic::GeodesicAlgorithmBase *algorithm;

  algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh); 

  double step = 1.0;
  const double eta = 100;

  algorithm->step = step;
  algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
  map<int,double> fixedDests;
  algorithm->propagate_vg(sources, eps_vg, fixedDests);





}

void figure_6()
{
  string input_obj_name = "fertility_nf30k_anisotropic.obj";
  CRichModel model();
  model.Preprocess();
  int source = 3465;//10848;//8056;
  //int dest = 5761;//14229;//5924;
  int dest;
  findDest(model,source,dest);


  double cylinder_radius = 0.002;
  CylinderPath cylinder_path(cylinder_radius);
  cylinder_path.addGeodesicPath(model,source,dest);
  cylinder_path.write_to_file("fertility_path_staight.obj");


  vector<int> sources;
  sources.push_back(source);
  CICHWithFurtherPriorityQueue algo_source(model,sources);
  algo_source.Execute();

  vector<int> dests;
  dests.push_back(dest);
  CICHWithFurtherPriorityQueue algo_dest(model,dests);
  algo_dest.Execute();

  double dis_exact = algo_source.m_InfoAtVertices[dest].disUptodate;
  set<int> path_verts;
  vector<IntersectionWithPath> result_path;
  algo_source.FindSourceVertex(dest,result_path);
  path_verts.insert(source);
  path_verts.insert(dest);
  for (auto& p:result_path) {
    if (p.isVertex) {
      path_verts.insert(p.index);
    }
  }

  vector<pair<int,double>> vert_dis;
  vector<pair<int,double>> vert_error;
  for (int i = 0; i < model.GetNumOfVerts(); ++i) {
    if (path_verts.find(i)!=path_verts.end()) continue;
    double d0 = algo_source.m_InfoAtVertices[i].disUptodate;
    double d1 = algo_dest.m_InfoAtVertices[i].disUptodate;
    vert_dis.push_back(make_pair(i,d0+d1));
    vert_error.push_back(make_pair(i,(d0+d1-dis_exact)/dis_exact));
  }
  sort(vert_error.begin(), vert_error.end(),
     [](const pair<int, double>& lhs, const pair<int, double>& rhs) {
             return lhs.second < rhs.second; } );
  for (int i = 0; i < 20; ++i) {
    printf("vert(i %d) error %lf \n" , vert_error[i].first, vert_error[i].second);
  }

  int mid_vert = vert_error[10].first;
  CylinderPath cylinder_path_eps(cylinder_radius);
  cylinder_path_eps.addGeodesicPath(model,source,mid_vert);
  cylinder_path_eps.addGeodesicPath(model,mid_vert,dest);
  cylinder_path_eps.write_to_file("fertility_path_eps.obj");
  
  vector<CPoint3D> points;
  points.push_back(model.Vert(source));
  points.push_back(model.Vert(dest));
  points.push_back(model.Vert(mid_vert));
  printBallToObj(points, "points_2.obj", cylinder_radius * 5);

}


int main()
{
  
  figure_6();






  return 0;

}