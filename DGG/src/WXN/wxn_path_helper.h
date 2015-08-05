#ifndef _WXN_PATH_HELPER_H_
#define _WXN_PATH_HELPER_H_
#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"

void cylinder(CPoint3D start_p, CPoint3D end_p);

void generateCylinder(CPoint3D start_p, CPoint3D end_p, vector<CPoint3D>& verts, vector<CBaseModel::CFace>& faces, double radius);

void generateArrow(CPoint3D start_p, CPoint3D end_p, vector<CPoint3D>& verts, vector<CBaseModel::CFace>& faces, double radius);

void printBallToObj(const vector<CPoint3D>& vertex_list, const string& file_name, double scale);

void output_faces(const CRichModel* model_ptr, const set<int>& faces, const string& output_filename);

void output_cylinder(const string& filename, const vector<CPoint3D>& verts, const vector<CBaseModel::CFace>& faces);

class CylinderPath{
  vector<CPoint3D> verts_;
  vector<CBaseModel::CFace> faces_;
  double radius_;
public:
  CylinderPath(double radius):radius_(radius){}
  void addGeodesicPath(CRichModel& model, int v0, int v1)
  {
    vector<int> sources;
    sources.push_back(v0);
    CICHWithFurtherPriorityQueue ich_algoritm(model,sources);
    ich_algoritm.Execute();

    vector<CPoint3D> path_points;
    vector<IntersectionWithPath> paths;
    ich_algoritm.FindSourceVertex(v1, paths);
    for (auto& v:paths) {
      path_points.push_back(v.GetPosition(model));
    }

    //CylinderPath cylinder_path(0.002);
    for (int i = 0; i < path_points.size() - 1; ++i) {
      addLine(path_points[i],path_points[i+1]);
    }
  }

  void addLine(const CPoint3D&p0, const CPoint3D& p1)
  {
    generateCylinder(p0, p1, verts_,faces_, radius_);
  }
  void write_to_file(const string& filename) {
    output_cylinder(filename, verts_, faces_);
  }


};

#endif
