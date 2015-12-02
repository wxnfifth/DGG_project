#include "stdafx.h"
#include "YXPathTracer.h"

int g_processindex = 0;

void YXPathTracer::unfoldStraightLine(int base_edge, double position, double dis1, double dis2, double remain_length, YXPath & path) const
{
  if(remain_length < EPS) return;
  base_edge = metric.Edge(base_edge).reverseEdgeIndex;
  double base_len = metric.Edge(base_edge).length;
  position = base_len - position;
  YXPoint2D intersectPoint(position, 0.0);
  const YXPoint2D v1(0.0, 0.0);
  YXPoint2D v2(base_len, 0.0);
  YXPoint2D & v3 = metric.ptToEdge2D[base_edge];
  YXPoint2D srcPoint = YXGeometry::compute3Dto2D(base_len, dis2, dis1);
  srcPoint.y = - srcPoint.y;
  if((v3 - srcPoint).crossProduct(intersectPoint - srcPoint) >= 0.0 ){
    // on left'
    int res;
    double b = YXGeometry::computeLineIntersection(v3, v1, srcPoint, intersectPoint, res);
    double a = 1.0 - b;
    YXPoint2D newIntersectionPoint = v3.scale(a);
    double current_seg_length = (intersectPoint - newIntersectionPoint).length();
    if(current_seg_length > remain_length){
      YXPoint2D endpoint = intersectPoint + (newIntersectionPoint - intersectPoint).scale(remain_length / current_seg_length);
      double area = fabs((v3-v2).crossProduct(v1-v2));
      a = fabs((v3-v2).crossProduct(endpoint - v2)) / area;
      b = fabs((v3-v1).crossProduct(endpoint-v1)) / area;
      path.push_back(YXPathPoint(base_edge, a, b));
    }else{
      base_edge = metric.getPrevEdgeAroundFace(base_edge);
      position = metric.Edge(base_edge).length * b;
      dis1 = (srcPoint - v3).length();
      path.push_back(YXPathPoint(base_edge, a, b));
      unfoldStraightLine(base_edge, position, dis1, dis2, remain_length - current_seg_length, path);
    }
  }else{
    // on right
    int res;
    double b = YXGeometry::computeLineIntersection(v2, v3, srcPoint, intersectPoint, res);
    double a = 1.0 - b;
    YXPoint2D newIntersectionPoint = v2.scale(a) + v3.scale(b);
    double current_seg_length = (intersectPoint - newIntersectionPoint).length();
    if(current_seg_length > remain_length){
      YXPoint2D endpoint = intersectPoint + (newIntersectionPoint - intersectPoint).scale(remain_length / current_seg_length);
      double area = fabs((v3-v2).crossProduct(v1-v2));
      a = fabs((v3-v2).crossProduct(endpoint - v2)) / area;
      b = fabs((v3-v1).crossProduct(endpoint-v1)) / area;
      path.push_back(YXPathPoint(base_edge, a, b));
    }else{
      base_edge = metric.getNextEdgeAroundFace(base_edge);
      position = metric.Edge(base_edge).length * b;
      dis2 = (srcPoint - v3).length();
      path.push_back(YXPathPoint(base_edge, a, b));
      unfoldStraightLine(base_edge, position, dis1, dis2, remain_length - current_seg_length, path);
    }
  }
}


void YXPathTracer::computeSinglePath(int start_edge, double angle, double length, YXPath & path) const {
  path.push_back(YXPathPoint(start_edge, 1.0, 0.0));
  int cross_edge = metric.getNextEdgeAroundFace(start_edge);
  double beta = metric.edgeAngle[cross_edge];
  double gamma = YX_PI - beta - angle;
  double base_len = metric.Edge(start_edge).length;
  double sin_gamma = sin(gamma);
  double position = sin(angle) * base_len / sin_gamma;
  double current_seg_length = sin(beta) * base_len / sin_gamma;
  double b = position / metric.Edge(cross_edge).length;
  double a = 1.0 - b;
  if(current_seg_length > length){
    double s = length / current_seg_length;
    path.push_back(YXPathPoint(cross_edge, a*s, b*s));
  }else{
    path.push_back(YXPathPoint(cross_edge, a, b));
    unfoldStraightLine(cross_edge, position, /*dis1*/base_len, /*dis2*/metric.Edge(metric.getPrevEdgeAroundFace(start_edge)).length, 
      length - current_seg_length, path);
  }
}


void YXPathTracer::computeSinglePathArbitraryStart(int start_edge,
                                                   double angle, double length,
                                                   YXPath & path,
                                                   int source_index
                                                   ) const
{
  //if (source_index == 27) {
    //int temp_e = start_edge;
    //while (true) {
    //  printf("e %d %lf v1 %d v2 %d source %d\n" , temp_e, metric.edgeAngle[temp_e],
    //                                              metric.Edge(temp_e).v1, metric.Edge(temp_e).v2,
    //                                              source_index + 5002);
    //  temp_e = metric.getNextEdgeAroundVertex(temp_e);
      //if (temp_e == start_edge) {
      //  break;
      //}
    //}
    //printf("\n");
  //}
  int temp_e = start_edge;
  vector<int> edge_sequece;
  while (true) {
    edge_sequece.push_back(temp_e);
    temp_e = metric.getNextEdgeAroundVertex(temp_e);
    if (temp_e == start_edge) {
      break;
    }
  }
  std::reverse(edge_sequece.begin() + 1 , edge_sequece.end());



  double current_angle = angle;
  int e;
  for (auto itr = edge_sequece.begin(); itr != edge_sequece.end(); ++itr) {
    e = *itr;
    if (current_angle > metric.edgeAngle[e]) {
      current_angle -= metric.edgeAngle[e];
    } else{
      break;
    }
  }
  //if (source_index == 27) {
  //  printf(" angle %lf current angle %lf\n" , angle , current_angle);
  //}
  //double current_angle = angle;
  //int e = start_edge;
  //while (current_angle > metric.edgeAngle[e]) {
  //  current_angle -= metric.edgeAngle[e];
  //  e = metric.getNextEdgeAroundVertex(e);
  //}
  //printf("current angle %lf\n", angle);
  //for (auto e : edge_sequece) {
	 // printf("tmp_e %d v0 %d v1 %d\n", e, metric.Edge(e).v1, metric.Edge(e).v2);
  //}
  computeSinglePath(e, current_angle, length, path);
}

void YXPathTracer::computePatch(int start_vertex, int numOfPaths, double length, vector<YXPath> &paths) const {
  double total_angle = metric.angle[start_vertex];
  double angle_inc = total_angle / numOfPaths;
  int e = metric.getAnyEdgeStartFromVert(start_vertex);
  double current_angle = 0.1;
  for (int i = 0; i < numOfPaths; ++i) {
    //printf("computing path %d\n", i);
    g_processindex = i;
    while (current_angle > metric.edgeAngle[e]) {
      current_angle -= metric.edgeAngle[e];
      e = metric.getNextEdgeAroundVertex(e);
    }
    YXPath tmp;
    computeSinglePath(e, current_angle, length, tmp);
    paths.push_back(tmp);
    current_angle += angle_inc;
  }
}

void YXPathTracer::generateObj(const vector<YXPath> &paths, const char * filename) const{
  FILE * fp = fopen(filename, "w");
  fprintf(fp, "g Curve\n");
  for(int i = 0; i < paths.size(); ++i){
    g_processindex = i;
    YXPoint3D prevPt;
    double len = 0.0, remain = 1.0;
    for(int j = 0; j < paths[i].size(); ++j){
      YXPoint3D pt = convertToYXPoint3D(paths[i][j]);
      if(j != 0){
        len += (pt - prevPt).length();
        remain -= (pt-prevPt).length();
      }
      prevPt = pt;
      fprintf(fp, "v %.9lf %.9lf %.9lf\n", pt.x, pt.y, pt.z);
    }
    //printf("total length = %lf\n", len);
  }
  int t = 0;
  for(int i = 0; i < paths.size(); ++i){
    for(int j = 1; j < paths[i].size(); ++j){
      fprintf(fp, "l %d %d\n", t+j, t+j+1);
    }
    t += paths[i].size();
  }
  fclose(fp);
}

double YXPathTracer::approxCurvature(int v, const YXPoint3D & B, const YXPoint3D & C) {
  YXPoint3D & center = mesh.vert[v];
  double a2 = (B-C).length2(), b2 = (center-C).length2(), c2 = (center - B).length2();
  double cosA = (b2 + c2 - a2)/2.0/sqrt(b2*c2);
  double curvature = 2.0 * sqrt(1.0 - cosA*cosA) / sqrt(a2);
  YXPoint3D t = YXGeometry::projectPointToLine(B, C, center);
  if( YXGeometry::projectScalar(center - t, mesh.vertNormal[v]) < 0.0){
    return -curvature;
  }
  return curvature;
}

YXPoint3D YXPathTracer::convertToYXPoint3D(const YXPathPoint & p) const {
  int v1 = metric.Edge(p.edgeId).v1;
  int v2 = metric.Edge(p.edgeId).v2;
  int v3 = metric.Edge(metric.getNextEdgeAroundFace(p.edgeId)).v2;
  return mesh.vert[v1].scale(p.a) + mesh.vert[v2].scale(p.b) + mesh.vert[v3].scale(1.0-p.a-p.b);
}

int YXPathTracer::getFaceIndex(const YXPathPoint& p) const {
	return metric.faceIndex(p.edgeId);
}


