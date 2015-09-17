#ifndef _WXN_TRIANGLE_H_
#define _WXN_TRIANGLE_H_
#include "stdafx.h"
#include "wxn_geometry.h"

void triangulateSegLoop(const vector<CPoint3D>& ptInList,vector<CPoint3D>& ptOutList,vector<CBaseModel::CFace>& triList,vector<vector<int>>* segListPt=NULL,char* configCmd="pzQ");
void triangulateSegLoop(const vector<CPoint3D>& ptInList,const vector<vector<int>> segInList,vector<CPoint3D>& ptOutList,vector<CBaseModel::CFace>& triList,vector<vector<int>>& segOutList,char* configCmd="pzQ");

void triangulateSegLoop(const vector<Point2D>& ptInList,const vector<vector<int>> segInList,vector<Point2D>& ptOutList,vector<CBaseModel::CFace>& triList,vector<vector<int>>& segOutList,char* configCmd="pzQ");

//pzQ  pDS10zQ;
#endif
