#include "stdafx.h"
#include "ich\Point3D.h"
#include "ich\BaseModel.h"
#include "wxn_triangle.h"
#include "triangle.h"

void triangulateSegLoop(const vector<CPoint3D>& ptInList,vector<CPoint3D>& ptOutList,vector<CBaseModel::CFace>& triList,vector<vector<int>>* segListPt,char* configCmd)//pzQ  
	//if want better result can use,pDS10zQ;
{
	triangulateio ptIn,ptOut,vorout;
	ptIn.numberofpoints = ptInList.size()+0;
	ptIn.numberofpointattributes = 0;
	ptIn.pointlist = new REAL[ptIn.numberofpoints * 2];
	ptIn.pointattributelist = NULL;
	ptIn.pointmarkerlist = new int[ptIn.numberofpoints];
	ptIn.numberoftriangles = 0;

	ptIn.numberofsegments = ptInList.size();
	ptIn.segmentlist = new int [ptIn.numberofsegments * 2];
	ptIn.segmentmarkerlist = NULL;
	ptIn.numberofregions = 0;
	ptIn.regionlist = NULL;
	ptIn.numberofholes = 0;

	for(int i = 0; i < (int)ptInList.size();++i){
		ptIn.pointlist[i*2] = ptInList[i].x;
		ptIn.pointlist[i*2+1] = ptInList[i].y;
	}
	for(int i = 0; i < (int)ptIn.numberofpoints;++i){
		ptIn.pointmarkerlist[i] = 1;
	}
	for(int i = 0; i < (int)ptInList.size();++i){
		ptIn.segmentlist[i*2] = i;
		ptIn.segmentlist[i*2+1] = (i+1)%ptInList.size();
	}
	ptOut.pointlist = NULL;
	ptOut.pointmarkerlist = NULL;
	ptOut.segmentlist = NULL;
	ptOut.segmentmarkerlist = NULL;
	ptOut.trianglelist = NULL;

	//if( printReport ){
	//report(&ptIn,0,1,0,1,0,0);
	//}
	//"pDS10zQ"
	triangulate(configCmd,&ptIn,&ptOut,&vorout);
	//if( printReport ){
	//	report(&ptOut,0,1,0,1,0,1);
	//}

	ptOutList.resize(ptOut.numberofpoints);
	for(int i = 0; i < ptOut.numberofpoints;++i){
		ptOutList[i]=CPoint3D(ptOut.pointlist[i*2],ptOut.pointlist[i*2+1],0);
	}

	if( ptOut.numberofcorners != 3 ){
		printf( "error ptOut.numberofcorners != 3\n" );
		return;
	}

	triList.resize(ptOut.numberoftriangles);
	for(int i = 0; i < ptOut.numberoftriangles; ++i){
		for(int j = 0; j < ptOut.numberofcorners;++j){
			triList[i][j] = ptOut.trianglelist[i*ptOut.numberofcorners+j];
		}
	}
	if( segListPt != NULL ){
		(*segListPt).resize(ptOut.numberofsegments);
		for(int i = 0; i < ptOut.numberofsegments;++i){
			for(int j = 0; j < 2;++j){
				(*segListPt)[i].push_back(ptOut.segmentlist[i*2+j]);
			}
		}
	}
}
void triangulateSegLoop(const vector<CPoint3D>& ptInList,const vector<vector<int>> segInList,vector<CPoint3D>& ptOutList,vector<CBaseModel::CFace>& triList,vector<vector<int>>& segOutList,char* configCmd)
{
	triangulateio ptIn,ptOut,vorout;
	ptIn.numberofpoints = ptInList.size()+0;
	ptIn.numberofpointattributes = 0;
	ptIn.pointlist = new REAL[ptIn.numberofpoints * 2];
	ptIn.pointattributelist = NULL;
	ptIn.pointmarkerlist = new int[ptIn.numberofpoints];
	ptIn.numberoftriangles = 0;
	
	ptIn.numberofsegments = segInList.size();
	ptIn.segmentlist = new int [ptIn.numberofsegments * 2];
	ptIn.segmentmarkerlist = NULL;
	ptIn.numberofregions = 0;
	ptIn.regionlist = NULL;
	ptIn.numberofholes = 0;

	for(int i = 0; i < ptInList.size();++i){
		ptIn.pointlist[i*2] = ptInList[i].x;
		ptIn.pointlist[i*2+1] = ptInList[i].y;
	}
	for(int i = 0; i < ptIn.numberofpoints;++i){
		ptIn.pointmarkerlist[i] = 1;
	}

	/*for(int i = 0; i < num;++i){
		ptIn.segmentlist[i*2] = i;
		ptIn.segmentlist[i*2+1] = (i+1)%(num);
	}*/
	for(int i = 0; i < segInList.size();++i){
		ptIn.segmentlist[i*2] = segInList[i][0];
		ptIn.segmentlist[i*2+1] = segInList[i][1];
	}

	ptOut.pointlist = NULL;
	ptOut.pointmarkerlist = NULL;
	ptOut.segmentlist = NULL;
	ptOut.segmentmarkerlist = NULL;
	ptOut.trianglelist = NULL;

	//"pDS10zQ"
	triangulate(configCmd,&ptIn,&ptOut,&vorout);


	ptOutList.resize(ptOut.numberofpoints);
	for(int i = 0; i < ptOut.numberofpoints;++i){
		ptOutList[i]=CPoint3D(ptOut.pointlist[i*2],ptOut.pointlist[i*2+1],0);
	}

	if( ptOut.numberofcorners != 3 ){
		printf( "error ptOut.numberofcorners != 3\n" );
		return;
	}

	triList.resize(ptOut.numberoftriangles);
	for(int i = 0; i < ptOut.numberoftriangles; ++i){
		for(int j = 0; j < ptOut.numberofcorners;++j){
			triList[i][j] = ptOut.trianglelist[i*ptOut.numberofcorners+j];
		}
	}
	segOutList.resize(ptOut.numberofsegments);
	for(int i = 0; i < ptOut.numberofsegments;++i){
		for(int j = 0; j < 2;++j){
			segOutList[i].push_back(ptOut.segmentlist[i*2+j]);
		}
	}
}
void triangulateSegLoop(const vector<Point2D>& ptInList,const vector<vector<int>> segInList,vector<Point2D>& ptOutList,vector<CBaseModel::CFace>& triList,vector<vector<int>>& segOutList,char* configCmd)
{
	triangulateio ptIn,ptOut,vorout;
	ptIn.numberofpoints = ptInList.size()+0;
	ptIn.numberofpointattributes = 0;
	ptIn.pointlist = new REAL[ptIn.numberofpoints * 2];
	ptIn.pointattributelist = NULL;
	ptIn.pointmarkerlist = new int[ptIn.numberofpoints];
	ptIn.numberoftriangles = 0;
	
	ptIn.numberofsegments = segInList.size();
	ptIn.segmentlist = new int [ptIn.numberofsegments * 2];
	ptIn.segmentmarkerlist = NULL;
	ptIn.numberofregions = 0;
	ptIn.regionlist = NULL;
	ptIn.numberofholes = 0;

	for(int i = 0; i < ptInList.size();++i){
		ptIn.pointlist[i*2] = ptInList[i].x;
		ptIn.pointlist[i*2+1] = ptInList[i].y;
	}
	for(int i = 0; i < ptIn.numberofpoints;++i){
		ptIn.pointmarkerlist[i] = 1;
	}

	/*for(int i = 0; i < num;++i){
		ptIn.segmentlist[i*2] = i;
		ptIn.segmentlist[i*2+1] = (i+1)%(num);
	}*/
	for(int i = 0; i < segInList.size();++i){
		ptIn.segmentlist[i*2] = segInList[i][0];
		ptIn.segmentlist[i*2+1] = segInList[i][1];
	}

	ptOut.pointlist = NULL;
	ptOut.pointmarkerlist = NULL;
	ptOut.segmentlist = NULL;
	ptOut.segmentmarkerlist = NULL;
	ptOut.trianglelist = NULL;

	//"pDS10zQ"
	triangulate(configCmd,&ptIn,&ptOut,&vorout);


	ptOutList.resize(ptOut.numberofpoints);
	for(int i = 0; i < ptOut.numberofpoints;++i){
		ptOutList[i]=Point2D(ptOut.pointlist[i*2],ptOut.pointlist[i*2+1]);
	}

	if( ptOut.numberofcorners != 3 ){
		printf( "error ptOut.numberofcorners != 3\n" );
		return;
	}

	triList.resize(ptOut.numberoftriangles);
	for(int i = 0; i < ptOut.numberoftriangles; ++i){
		for(int j = 0; j < ptOut.numberofcorners;++j){
			triList[i][j] = ptOut.trianglelist[i*ptOut.numberofcorners+j];
		}
	}
	segOutList.resize(ptOut.numberofsegments);
	for(int i = 0; i < ptOut.numberofsegments;++i){
		for(int j = 0; j < 2;++j){
			segOutList[i].push_back(ptOut.segmentlist[i*2+j]);
		}
	}
}