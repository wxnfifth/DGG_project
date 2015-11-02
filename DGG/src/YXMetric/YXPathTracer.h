#ifndef _YXPATHTRACER_H_
#define _YXPATHTRACER_H_
#include "YXMetric.h"
using std::vector;


struct YXPathPoint{
  //edge: v1-->v2, on triangle v1v2v3, then the point is : v1*a+v2*b+v3*(1-a-b)
  int edgeId;
  double a, b;
  YXPathPoint() : edgeId(0), a(0.0), b(0.0) {}
  YXPathPoint(int e, double a, double b) : edgeId(e), a(a), b(b) {}
};

typedef vector<YXPathPoint> YXPath;


class YXPathTracer{
public:

  YXMesh3D mesh;
  YXMetric metric;

  void init(const char * filename) {
    int rt_value = mesh.load(filename);
    assert( rt_value == 0);
    metric.buildFromMesh(mesh);
    metric.assignNormalsForMesh(mesh);
  }

  void unfoldStraightLine(int base_edge, double position, double dis1, double dis2, double remain_length, YXPath & path) const;


  void computeSinglePath(int start_edge, double angle, double length, YXPath & path) const;
  void computeSinglePathArbitraryStart(int start_edge, double angle, double length, YXPath & path , int source_index = -1) const;


  void computePatch(int start_vertex, int numOfPaths, double length, vector<YXPath> &paths) const; 

  YXPoint3D convertToYXPoint3D(const YXPathPoint & p) const;
  int getFaceIndex(const YXPathPoint& p) const;

  void generateObj(const vector<YXPath> &paths, const char * filename) const;

  double approxCurvature(int v, const YXPoint3D & B, const YXPoint3D & C);
};
#endif