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

	void addGeodesicPaths(CRichModel& model, int v0, const vector<int>& vts)
  {
    vector<int> sources;
    sources.push_back(v0);
    CICHWithFurtherPriorityQueue ich_algoritm(model,sources);
    ich_algoritm.Execute();

		for (auto& v:vts) {
			vector<CPoint3D> path_points;
			vector<IntersectionWithPath> paths;
			ich_algoritm.FindSourceVertex(v, paths);
			for (auto& p:paths) {
				path_points.push_back(p.GetPosition(model));
			}

			//CylinderPath cylinder_path(0.002);
			for (int i = 0; i < path_points.size() - 1; ++i) {
				addLine(path_points[i],path_points[i+1]);
			}
		}
  }

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
	void addGeodesicPaths(CRichModel& model, vector<int>& vts)
  {
    for (int i = 0; i < vts.size() - 1; ++i) {
      addGeodesicPath(model,vts[i],vts[i+1]);
    }
  }

  void addLine(const CPoint3D&p0, const CPoint3D& p1)
  {
    generateCylinder(p0, p1, verts_,faces_, radius_);
  }
  void addLine(const CPoint3D&p0, const CPoint3D& p1,const double len)
  {
    auto& p_end = p0 + (p1-p0).Normalize() * len;
    generateCylinder(p0, p_end, verts_,faces_, radius_);
  }
  void addLines(const vector<CPoint3D>& pts)
  {
    for (int i = 0; i < pts.size() - 1; ++i) {
      addLine(pts[i], pts[i+1]);
    }
  }
  void write_to_file(const string& filename) {
    output_cylinder(filename, verts_, faces_);
  }


};

#endif
