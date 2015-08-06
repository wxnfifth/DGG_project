#include "wxn\wxn_path_helper.h"
#include "mmp\geodesic_algorithm_vg_mmp.h"
#include "svg_definition.h"


struct BodyPartOfSVGWithPoint: public BodyPartOfSVGWithAngle{
  CPoint3D p;
    BodyPartOfSVGWithPoint(){}
    BodyPartOfSVGWithPoint(int _dest_index , float _dest_dis, float _angle, int _begin_pos, int _end_pos, const CPoint3D& _p):
      BodyPartOfSVGWithAngle(_dest_index,_dest_dis,_angle,_begin_pos,_end_pos),p(_p){}
};


void computeDests(CRichModel& model, const string& input_obj_name,
                    const double eps_vg, const int source, 
                    set<int>& dest_verts) 
{
  double const_for_theta = 5;
  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  fprintf(stderr,"******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;

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

  int source_index = source;
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

  vector<node> covered_points;
  covered_points.resize(fixedDests.size());
  int _index = 0;
  for (auto& d:fixedDests) {
    geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
    double dis;
    algorithm->best_source(dest_p, dis);
    covered_points[_index].id = d.first;
    covered_points[_index].dis = d.second;
    _index ++;
    dests.push_back(make_pair(d.first,d.second));
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

  vector<double> angles(body_parts.size());
  vector<CPoint3D> fan_pts(body_parts.size());
  for (int i = 0; i < body_parts.size(); ++i) {
    auto dest_vert(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
    vector<geodesic::SurfacePoint> path; 
    algorithm->trace_back(dest_vert, path);
    geodesic::SurfacePoint p;
    if (path.size() > 1) {
      p = *(path.rbegin()+1);
    } else {
      p = (path[0]);
    }
    //   VERTEX,
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
    fan_pts[i] = CPoint3D(p.x(),p.y(),p.z());

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

  vector<BodyPartOfSVGWithPoint> body_parts_with_angle(body_parts.size());
  for (int i = 0; i < body_parts.size(); ++i) {
    BodyPartOfSVG b = body_parts[i];
    BodyPartOfSVGWithPoint b_with_angle(b.dest_index, b.dest_dis, angles[i],0,0,fan_pts[i]);
    body_parts_with_angle[i] = b_with_angle;
    //output_file.write( (char*)&b_with_angle , sizeof(b_with_angle));
  }
  sort(body_parts_with_angle.begin(), body_parts_with_angle.end());
  float angle_sum = model.AngleSum(source_index);
  vector<double> tmp_angles(body_parts_with_angle.size()*2);
  for (int i = 0; i < body_parts_with_angle.size(); ++i) {
    tmp_angles[i] = body_parts_with_angle[i].angle;
  }
  for (int i = body_parts_with_angle.size(); i < tmp_angles.size(); ++i) {
    tmp_angles[i] = body_parts_with_angle[i-body_parts_with_angle.size()].angle + angle_sum;
  }

  auto compareFunc = [](pair<BodyPartOfSVGWithPoint,int>& a, pair<BodyPartOfSVGWithPoint,int>& b) {return a.first.dest_dis > b.first.dest_dis;};
  typedef priority_queue<pair<BodyPartOfSVGWithPoint,int>, vector<pair<BodyPartOfSVGWithPoint,int>>, decltype(compareFunc)> myque_type;
  myque_type que(compareFunc);

  //for (auto& a:tmp_angles) {
  //  printf("%lf " , a);
  //}
  //printf("\n");

  for (int i = 0; i < body_parts_with_angle.size(); ++i) {//assume i is father
    que.push(make_pair(body_parts_with_angle[i],i));
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
    if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
    if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
    body_parts_with_angle[i].begin_pos = start_pos;
    body_parts_with_angle[i].end_pos = end_pos;
  }






}
void computeDestFan(CRichModel& model, const string& input_obj_name,
                    const double eps_vg, const int source, 
                    int& dest1, CPoint3D&fan1_p1, CPoint3D& fan1_p2, 
                    int& dest2, CPoint3D&fan2_p1, CPoint3D& fan2_p2)
{
  double const_for_theta = 5;
  double theta = asin(sqrt(eps_vg));
  theta *= const_for_theta;
  fprintf(stderr,"******** eps %lf const %lf theta %lf du\n" , eps_vg, const_for_theta, theta  / M_PI * 180.0);

  std::vector<double> points;	
  std::vector<unsigned> faces;
  std::vector<int> realIndex;
  int originalVertNum = 0;

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

  int source_index = source;
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

  vector<node> covered_points;
  covered_points.resize(fixedDests.size());
  int _index = 0;
  for (auto& d:fixedDests) {
    geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
    double dis;
    algorithm->best_source(dest_p, dis);
    covered_points[_index].id = d.first;
    covered_points[_index].dis = d.second;
    _index ++;
    dests.push_back(make_pair(d.first,d.second));
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

  vector<double> angles(body_parts.size());
  vector<CPoint3D> fan_pts(body_parts.size());
  for (int i = 0; i < body_parts.size(); ++i) {
    auto dest_vert(geodesic::SurfacePoint(&mesh.vertices()[body_parts[i].dest_index]));
    vector<geodesic::SurfacePoint> path; 
    algorithm->trace_back(dest_vert, path);
    geodesic::SurfacePoint p;
    if (path.size() > 1) {
      p = *(path.rbegin()+1);
    } else {
      p = (path[0]);
    }
    //   VERTEX,
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
    fan_pts[i] = CPoint3D(p.x(),p.y(),p.z());

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

  vector<BodyPartOfSVGWithPoint> body_parts_with_angle(body_parts.size());
  for (int i = 0; i < body_parts.size(); ++i) {
    BodyPartOfSVG b = body_parts[i];
    BodyPartOfSVGWithPoint b_with_angle(b.dest_index, b.dest_dis, angles[i],0,0,fan_pts[i]);
    body_parts_with_angle[i] = b_with_angle;
    //output_file.write( (char*)&b_with_angle , sizeof(b_with_angle));
  }
  sort(body_parts_with_angle.begin(), body_parts_with_angle.end());
  float angle_sum = model.AngleSum(source_index);
  vector<double> tmp_angles(body_parts_with_angle.size()*2);
  for (int i = 0; i < body_parts_with_angle.size(); ++i) {
    tmp_angles[i] = body_parts_with_angle[i].angle;
  }
  for (int i = body_parts_with_angle.size(); i < tmp_angles.size(); ++i) {
    tmp_angles[i] = body_parts_with_angle[i-body_parts_with_angle.size()].angle + angle_sum;
  }

  auto compareFunc = [](pair<BodyPartOfSVGWithPoint,int>& a, pair<BodyPartOfSVGWithPoint,int>& b) {return a.first.dest_dis > b.first.dest_dis;};
  typedef priority_queue<pair<BodyPartOfSVGWithPoint,int>, vector<pair<BodyPartOfSVGWithPoint,int>>, decltype(compareFunc)> myque_type;
  myque_type que(compareFunc);

  //for (auto& a:tmp_angles) {
  //  printf("%lf " , a);
  //}
  //printf("\n");

  for (int i = 0; i < body_parts_with_angle.size(); ++i) {//assume i is father
    que.push(make_pair(body_parts_with_angle[i],i));
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
    if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
    if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
    body_parts_with_angle[i].begin_pos = start_pos;
    body_parts_with_angle[i].end_pos = end_pos;
  }
   for (int i = 0; i < 5; ++i) {
 // for (int i = 0; i < 5; ++i) {
    //printf("d %lf \n" , que.top().dest_dis);
    que.pop();
  }


  int result_pos1 = que.top().second;
  auto p1 = body_parts_with_angle[result_pos1];
  dest1 = p1.dest_index;
  printf("origin_pos %d begin_pos %d end_pos %d\n" , result_pos1, p1.begin_pos, p1.end_pos);
  printf("dis %lf\n" , p1.dest_dis);
  fan1_p1 = body_parts_with_angle[p1.begin_pos].p;
  fan1_p2 = body_parts_with_angle[p1.end_pos].p;

  int result_pos2 = result_pos1 + body_parts_with_angle.size() / 4;
  while(!que.empty() ){
    if ( abs(que.top().second - result_pos2) < 10) {
      break;
    }
    que.pop();
  }
  result_pos2 = que.top().second;

  auto p = body_parts_with_angle[result_pos2];
  dest2 = p.dest_index;
  printf("origin_pos %d begin_pos %d end_pos %d\n" , result_pos2, p.begin_pos, p.end_pos);
  printf("dis %lf\n" , p.dest_dis);
  fan2_p1 = body_parts_with_angle[p.begin_pos].p;
  fan2_p2 = body_parts_with_angle[p.end_pos].p;


}

void figure8_a()
{
  string model_filename = "bunny_nf36k.obj";
  CRichModel model(model_filename);
  model.Preprocess();
  int source = 14706;
  int dest1;
  CPoint3D fan1_p1,fan1_p2;
  int dest2;
  CPoint3D fan2_p1,fan2_p2;
  double eps_vg = 0.001;
  double cylinder_radius = 0.0002;
  computeDests(model, model_filename,eps_vg,source,dest1,fan1_p1,fan1_p2,dest2,fan2_p1,fan2_p2);
  CylinderPath cylinder1(cylinder_radius);
  cylinder1.addGeodesicPath(model,source,dest1);
  cylinder1.addLine(model.Vert(source),fan1_p1);
  cylinder1.addLine(model.Vert(source),fan1_p2);
  cylinder1.write_to_file("bunny_fan1.obj");
  CylinderPath cylinder2(cylinder_radius);
  cylinder2.addGeodesicPath(model,source,dest2);
  cylinder2.addLine(model.Vert(source),fan2_p1);
  cylinder2.addLine(model.Vert(source),fan2_p2);
  cylinder2.write_to_file("bunny_fan2.obj");
  
  vector<CPoint3D> source_verts;
  source_verts.push_back(model.Vert(source));
  printBallToObj(source_verts, "bunny_source.obj", cylinder_radius * 6);
  vector<CPoint3D> dest_verts;
  dest_verts.push_back(model.Vert(dest1));
  dest_verts.push_back(model.Vert(dest2));
  printBallToObj(dest_verts, "bunny_dest.obj", cylinder_radius * 4);
  


}


void figure8_b()
{
  string model_filename = "bunny_nf36k.obj";
  CRichModel model(model_filename);
  model.Preprocess();
  int source = 14706;
  int dest1;
  CPoint3D fan1_p1,fan1_p2;
  int dest2;
  CPoint3D fan2_p1,fan2_p2;
  double eps_vg = 0.001;
  double cylinder_radius = 0.0002;
  computeDestFan(model, model_filename,eps_vg,source,dest1,fan1_p1,fan1_p2,dest2,fan2_p1,fan2_p2);
  CylinderPath cylinder1(cylinder_radius);
  cylinder1.addGeodesicPath(model,source,dest1);
  cylinder1.addLine(model.Vert(source),fan1_p1);
  cylinder1.addLine(model.Vert(source),fan1_p2);
  cylinder1.write_to_file("bunny_fan1.obj");
  CylinderPath cylinder2(cylinder_radius);
  cylinder2.addGeodesicPath(model,source,dest2);
  cylinder2.addLine(model.Vert(source),fan2_p1);
  cylinder2.addLine(model.Vert(source),fan2_p2);
  cylinder2.write_to_file("bunny_fan2.obj");
  
  vector<CPoint3D> source_verts;
  source_verts.push_back(model.Vert(source));
  printBallToObj(source_verts, "bunny_source.obj", cylinder_radius * 6);
  vector<CPoint3D> dest_verts;
  dest_verts.push_back(model.Vert(dest1));
  dest_verts.push_back(model.Vert(dest2));
  printBallToObj(dest_verts, "bunny_dest.obj", cylinder_radius * 4);
  


}


int main()
{
  figure8_b();

}