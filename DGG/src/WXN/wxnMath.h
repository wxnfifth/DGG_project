#ifndef _WXNMATH_H_
#define _WXNMATH_H_

#include <cmath>
#include <GL/glu.h>
#include "ich\Point3D.h"
#include "ich\BaseModel.h"
#include "ich\RichModel.h"
#include "wxn_geometry.h"
#include "wxn_triangle.h"

//color reference http://www.tayloredmktg.com/rgb/
const CPoint3D color_black(0,0,0);
const CPoint3D color_white(1,1,1);
const CPoint3D color_666699(0x66/255.0,0x66/255.0,0x99/255.0);
const CPoint3D color_ff9999(0xff/255.0,0x99/255.0,0x99/255.0);
const CPoint3D color_hot_pink(255.0/255.0,105.0/255.0,180.0/255.0);
const CPoint3D color_pink(255.0/255.0,192.0/255.0,203.0/255.0);//	ffc0cb
const CPoint3D color_deep_pink(255.0/255.0,20.0/255.0,147.0/255.0);//	ff1493
const CPoint3D color_dodger_blue(30.0/255.0,144.0/255.0,255.0/255.0);//	1e90ff
const CPoint3D color_dark_sea_green(143/255.0,188/255.0,143/255.0);
const CPoint3D color_dark_violet(148.0/255.0,0.0/255.0,211.0/255.0);
const CPoint3D color_blue(0,0,1);
const CPoint3D color_green(0,1,0);
const CPoint3D color_red(1,0,0);
const CPoint3D color_yellow(1,1,0);
const CPoint3D color_indian_red(205.0/255.0,92.0/255.0,92.0/255.0);
const CPoint3D color_forest_green(34.0/255.0,139.0/255.0,34.0/255.0);
const CPoint3D color_sky_blue(135/255.0,206.0/255.0,250/255.0);
const CPoint3D color_saddle_brown(139.0/255.0,69.0/255.0,19.0/255.0);


CPoint3D GetOGLPos(int x, int y);


#define IS_NAN(x) ((x)!=(x))

//void calculateCoveredFaces(const CRichModel& model , const map<int, Point2D>& vertex_set , set<int>& face_set );
void renderTorus(const CPoint3D& p , const CPoint3D& color , const double radius = 0.003 , const int& sub_times = 10);
void renderSphere(const CPoint3D& p , const CPoint3D& color = color_red , const double radius = 0.003 , const int& sub_times = 10);

void renderCylinder(float x1, float y1, float z1, float x2,float y2, float z2, float radius,int subdivisions,GLUquadricObj *quadric);


void renderCylinder_convenient(float x1, float y1, float z1, float x2,float y2, float z2, float radius,int subdivisions);

void renderCylinder_convenient(CPoint3D p1, CPoint3D p2 , float radius , int subdivisions=5);
double cos_from_edges(double const a,			//compute the cosine of the angle given the lengths of the edges
							 double const b,
							 double const c);


bool calculateBarycentric(const CPoint3D v[],const CPoint3D& p,CPoint3D& barycentric_cordinate);
bool calculateBarycentric(const CPoint3D v[],const CPoint3D& p,double barycentric_cordinate[]);


void solveQuadraticEquation(const double& a , const double& b ,const double& c,vector<double>& result);

double calcMidPoint(const double& dS1V1,const double& dS1V2,const double& dS2V1,const double &dS2V2,const double&dV1V2,double& distance_s1_middle );


//http://en.wikipedia.org/wiki/Circumscribed_circle
//refer to the url above for the equations
CPoint3D calcCircumcenter(CPoint3D& p1 , CPoint3D& p2,CPoint3D& p3);



//https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas
//see the site for reference

CPoint3D rotatePoint(const CPoint3D& ptv1,const CPoint3D& norm,const CPoint3D& ptv2,const double d1, const double& cosTheta);

void map3DTriangleTo2D(const CPoint3D pt[],CPoint3D pt2D[]);
void map3DTriangleTo2D(const CPoint3D pt[],Point2D pt2D[]);


void calculate3DTriangleCoordinateIn2D(const CPoint3D p[] , Point2D vertex_2d[]);
void calculate2DTriangleCoordinate(const double length[],Point2D vertex_2d[]);
double calculateOppositeVertexInquadrangle(const double length[]);
double calculateCosFromTriangleLength(const double& a,const double& b , const double& c,int* cnt=NULL);

CPoint3D ptIn3DEdgeTo2D(const CPoint3D& ptIn,const CPoint3D& EP3D1,const CPoint3D& EP3D2,const CPoint3D& EP2D1,const CPoint3D& EP2D2);
CPoint3D ptIn2DTriangleToBarycentric(const CPoint3D v2D[],const CPoint3D& p2D);

//http://www.blackpawn.com/texts/pointinpoly/default.html
CPoint3D pointIn3DTriangleToBarycentric(const CPoint3D v3D[],const CPoint3D& p3D);

CPoint3D pointIn2DTriangleToBarycentricRMS(const Point2D v2D[],const Point2D& p2D);
CPoint3D pointIn2DTriangleToBarycentric(const Point2D v2D[],const Point2D& p2D);
CPoint3D pointIn2DTriangleTo3D(const CPoint3D v[],const Point2D v2D[],const Point2D& p2D);
Point2D ptIn3DTriangleTo2D(const CPoint3D v3D[],const Point2D v2D[],const CPoint3D& p3D);

double ptInEdgeToProp(const CPoint3D& ptIn,const CPoint3D& EP1,const CPoint3D& EP2);

CPoint3D ptIn2DTriangleTo3D(const CPoint3D v[],const CPoint3D v2D[],const CPoint3D& p2D);
double cos_from_edges(double const a,			//compute the cosine of the angle given the lengths of the edges
							 double const b,
							 double const c);
double calculateTriangleQuality(const double& a ,const double& b,const double& c);

double calcAngleFromTriangle(const double& a , const double& b , const double& c);

void drawLine(const CPoint3D& p1,const CPoint3D& p2,vector<CPoint3D>&  vList,vector<vector<int>>& triList);


CPoint3D findU1(const CPoint3D vx2D,const CPoint3D& p1,double lP1Midr1,double theta);

CPoint3D findLineIntersectionPt(const CPoint3D& p12D,const CPoint3D& u1,const CPoint3D& p22D,const CPoint3D& u2);
double calculateOppositeVertexDistanceInquadrangle(const double length[]);


template<typename  T> bool vectorSortEqual(const vector<T>& va,const vector<T>& vb)
{
	if( va.size() != vb.size() ) {
		//printf( "vector size not equal\n");
		return false;
	}
	vector<T>a(va);
	vector<T>b(vb);
	sort(a.begin(),a.end());
	sort(b.begin(),b.end());
	for(int i = 0; i < a.size();++i){
		if( a[i] != b[i] ){
			return false;
		}
	}
	return true;
}
template<typename  T> bool vectorSortEqual(const T va[],const T vb[],int n)
{
	vector<T>a(va,va+n);
	vector<T>b(vb,vb+n);
	sort(a.begin(),a.end());
	sort(b.begin(),b.end());
	for(int i = 0; i < a.size();++i){
		if( a[i] != b[i] ){
			return false;
		}
	}
	return true;
}
template<typename T> bool vectorEqual(const vector<T>&va,const vector<T>&vb)
{
	if( va.size() != vb.size() ) return false;
	for(int i = 0; i < (int)va.size();++i){
		if( va[i] != vb[i] ){
			return false;
		}
	}
	return true;
}


template<typename T> void printTwoDimentionArray(const vector<vector<T>>& v,const string& print_infomation = "")
{
	cout << print_infomation << "\n";
	for(int i = 0; i < v.size();++i){
		for(int j = 0; j < v[i].size();++j){
			cout << v[i][j] << " ";
		}
		cout << "\n";
	}
}

template<template <typename> class P = std::less >
struct compare_pair_second {
    template<class T1, class T2> bool operator()(const std::pair<T1, T2>& left, const std::pair<T1, T2>& right) {
        return P<T2>()(left.second, right.second);
    }
};

class BasicModel{
public:
	vector<CPoint3D> vertex_list_;
	vector<CPoint3D> normal_list_;
	vector<CBaseModel::CFace> face_list_;
	vector<pair<int,int>> edge_to_be_removed_;
	BasicModel(){}
	BasicModel(const vector<CPoint3D>& _vertex_list_,const vector<CBaseModel::CFace>& _face_list_):vertex_list_(_vertex_list_),face_list_(_face_list_){

	}
	BasicModel(const vector<CPoint3D>& _vertex_list_):vertex_list_(_vertex_list_){
		face_list_.clear();
	}

	int addVertex(const CPoint3D& p){
		vertex_list_.push_back(p);
		return vertex_list_.size()-1;
	}

	void addFace(const CBaseModel::CFace& face){
		if(face[0] != face[1] && face[1] != face[2] && face[2] != face[0] ){

			double length[3] = {0}; 
			for(int i = 0; i <3; ++i){
				length[i] = (vertex_list_[face[i]] - vertex_list_[face[(i+1)%3]]).Len();
			}
			double angle[3]={0};
			for(int i = 0; i < 3; ++i){
				angle[i] = acos(cos_from_edges(length[i],length[(i+1)%3],length[(i+2)%3]));
				if(angle[i] < 1e-5 ){
					printf("angle %lf\n" , angle[i]);
					pair<int,int> edge = pair<int,int>(face[i],face[(i+1)%3]);
					if(edge.first > edge.second ){
						swap(edge.first,edge.second);
					}
					bool flag_found = false;
					for(int j = 0; j < edge_to_be_removed_.size();++j){
						if(edge_to_be_removed_[j] == edge ){
							flag_found = true;
							break;
						}
					}
					if( !flag_found  ){
						edge_to_be_removed_.push_back(edge);
					}
					break;
				}
			}

			face_list_.push_back(face);
		}
	}

	void printStatic(){

		double min_angel = 1e50;
		double max_angel = -1.0;
		double min_triangle_quality=1e50;
		double average_triangle_quality=0.0;
		double min_average_angel = 0.0;
		for(int i = 0; i < face_list_.size();++i){
			CPoint3D* vertex[3];
			for(int j = 0; j < 3; ++j){
				vertex[j] = &vertex_list_[face_list_[i][j]];
			}
			double length[3] = {0}; 
			for(int j = 0; j <3; ++j){
				length[j] = (*vertex[j] - *vertex[(j+1)%3]).Len();
			}
			double angle[3]={0};
			double min_angle_in_triangle = 1e50;
			for(int j = 0; j < 3; ++j){
				angle[j] = acos(cos_from_edges(length[j],length[(j+1)%3],length[(j+2)%3]));
				min_angle_in_triangle = min(min_angle_in_triangle,angle[j]);
				min_angel = min(min_angel,angle[j]);
				max_angel = max(max_angel,angle[j]);
			}
			min_average_angel += min_angle_in_triangle;
			double triangle_quality = calculateTriangleQuality(length[0],length[1],length[2]);
			min_triangle_quality = min(min_triangle_quality,triangle_quality);
			average_triangle_quality += triangle_quality;
		}
		min_average_angel /= face_list_.size();
		average_triangle_quality /= face_list_.size();
		cout << "______________________statistic" << endl;
		cout << "averaget_triangle_quality " << average_triangle_quality << " min triangle quality " << min_triangle_quality << "\n";
		printf("min average angle %lf min angle %lf\n" , min_average_angel / PI * 180.0 , min_angel / PI * 180.0 );
		cout << "______________________statistics end" << endl;
	}

	void WriteToFile(const string& file_name ){
		
		vector<CBaseModel::CFace> temp_face_list(face_list_);
		map<int,int> vertex_to_be_removed;
		for(int i = 0; i < edge_to_be_removed_.size();++i){
			vertex_list_[edge_to_be_removed_[i].second] = CPoint3D::MAX_VALUE;
			if(vertex_to_be_removed.find(edge_to_be_removed_[i].second) != vertex_to_be_removed.end()){

			}else{
				vertex_to_be_removed[edge_to_be_removed_[i].second] = edge_to_be_removed_[i].first;
			}
			printf("%d %d\n" , edge_to_be_removed_[i].first,edge_to_be_removed_[i].second);
			for(int j = 0; j < face_list_.size();++j){
				vector<int> vertexs;
				for(int k = 0; k < 3; ++k){
					if( face_list_[j][k] == edge_to_be_removed_[i].first ||
						face_list_[j][k] == edge_to_be_removed_[i].second ){
							vertexs.push_back(k);
					}
				}

				if( vertexs.size() == 2 ){
					printf("face id %d vertexs size 2 \n" , j );
					temp_face_list[j] = CBaseModel::CFace(-1,-1,-1);
				}else if(vertexs.size() == 1){
					printf("face %d %d %d\n",face_list_[j][0],face_list_[j][1],face_list_[j][2]);
					printf("face id %d vertexs size 1 edge_to_be_removed_[i].first %d edge_to_be_removed_[i].second %d\n" , j ,edge_to_be_removed_[i].first,edge_to_be_removed_[i].second  );
					if(face_list_[j][0] >= 0 ){
						temp_face_list[j][vertexs[0]] = edge_to_be_removed_[i].first;
					}
				}
			}
		}
		face_list_.swap(temp_face_list);
		vector<int> vertex_index(vertex_list_.size());
		int pos = 0;
		for(int i = 0; i < vertex_index.size();++i){
			if( vertex_list_[i] != CPoint3D::MAX_VALUE ){
				vertex_index[i] = pos++;
			}else{
				vertex_index[i] = -1;
			}
		}
		for(int i = 0; i < face_list_.size();++i){
			if( face_list_[i][0] >= 0 ){
				int temp[3];
				copy ( face_list_[i].verts , face_list_[i].verts + 3 ,temp  );
				face_list_[i][0] = vertex_index[face_list_[i][0]];
				face_list_[i][1] = vertex_index[face_list_[i][1]];
				face_list_[i][2] = vertex_index[face_list_[i][2]];
				if(face_list_[i][0] < 0 || face_list_[i][1] < 0 || face_list_[i][2] < 0 ){
					printf("what the fuck face[%d] , %d %d %d origin %d %d %d\n" , i , face_list_[i][0],face_list_[i][1],face_list_[i][2] , temp[0],temp[1],temp[2]);
				}
			}
			/*if(face_list_[i][0] == 0 || face_list_[i][1] == 0 || face_list_[i][2] == 0 ){
				printf("what the fuck 0 face[%d] , %d %d %d\n" , i , face_list_[i][0],face_list_[i][1],face_list_[i][2]);
			}*/

		}

		cout << "write to file " << file_name << "\n";
		FILE* file_obj = fopen(file_name.c_str(),"w");
		if( file_obj == NULL ){
			printf("error in write to file %s\n", file_name.c_str());
			return ;
		}
		for(int j = 0; j < vertex_list_.size();++j){
			if( vertex_list_[j] != CPoint3D::MAX_VALUE ){
				fprintf(file_obj,"v %lf %lf %lf\n" , vertex_list_[j].x,vertex_list_[j].y,vertex_list_[j].z);
			}
		}
		for(int j = 0; j < face_list_.size();++j){
			if( face_list_[j][0] >= 0 ){
				fprintf(file_obj,"f %d %d %d\n" , face_list_[j][0]+1,face_list_[j][1]+1,face_list_[j][2]+1);
			}
			
			//if(face_list_[j][0] == 0 || face_list_[j][1] == 0 || face_list_[j][2] == 0 ){
			//	printf("what the fuck 0 face[%d] , %d %d %d\n" , j , face_list_[j][0],face_list_[j][1],face_list_[j][2]);
			//}

		}
		fclose(file_obj);
	}
};

class WXNModel{
public:
	vector<CPoint3D> vertex_list_;
	vector<CBaseModel::CFace> face_list_;
	const CRichModel& model_;
	//vector<vector<CPoint3D>> 
public:
	WXNModel(const CRichModel& model):model_(model){
		vertex_list_.assign(model.m_Verts.begin(),model.m_Verts.end());
		face_list_.assign(model.m_Faces.begin(),model.m_Faces.end());
	}
	void updateFace(int edge_pos,int face_index,vector<pair<int,int>>& new_faces,CPoint3D& vertex_for_insert){
		int vertex_index = vertex_list_.size();
		int edge_index = model_.GetEdgeIndexFromFace(face_index,edge_pos);
		int neighbor_face_index = model_.GetNeighborFaceIndexFromFace(face_index,edge_pos);
		int neighbor_vertex_index = model_.Edge(model_.Edge(edge_index).indexOfReverseEdge).indexOfOppositeVert;
		int v0 = face_list_[face_index][edge_pos];
		int v1 = face_list_[face_index][(edge_pos+1)%3];
		int v2 = face_list_[face_index][(edge_pos+2)%3];
		int v3 = neighbor_vertex_index;
		face_list_[face_index] = CBaseModel::CFace(v0,vertex_index,v2);
		new_faces.push_back( pair<int,int>( face_index,face_list_.size()) );
		face_list_.push_back(CBaseModel::CFace(vertex_index,v1,v2));
		face_list_[neighbor_face_index] = CBaseModel::CFace(v0,v3,vertex_index);
		new_faces.push_back( pair<int,int>( neighbor_face_index,face_list_.size()) );
		face_list_.push_back(CBaseModel::CFace(v3,v1,vertex_index));

	}
	void Subdivide(const CPoint3D& vertex,const int& face_index,vector<pair<int,int>>& new_faces,int& source_index,vector<int>& temp_vector=vector<int>()){//new faces's first is origin face
		
		if( face_index != -1 ){
			int vertex_index = vertex_list_.size();
			CPoint3D vertex_for_insert = vertex;
			if( face_index > face_list_.size() ){
				printf("error in Subdivide in WXNModel\n");
				return;
			}
			int v_index[3];
			CPoint3D v_3D[3];
			for(int j = 0; j < 3; ++j){
				v_index[j] = face_list_[face_index][j];
				v_3D[j] = model_.Vert(v_index[j]);
			}
			CPoint3D barycentric_coordinate;
			calculateBarycentric(v_3D,vertex,barycentric_coordinate);
			double area = ((v_3D[0]-v_3D[1])*(v_3D[0]-v_3D[2])).Len();
			double area_01 = ((v_3D[0]-v_3D[1])*(v_3D[0]-vertex)).Len();
			double area_12 = ((v_3D[1]-v_3D[2])*(v_3D[1]-vertex)).Len();
			double area_20 = ((v_3D[2]-v_3D[0])*(v_3D[2]-vertex)).Len();
			
			if(face_index == 1518){
				v_3D[0].Print("v_3D[0]");
				v_3D[1].Print("v_3D[1]");
				v_3D[2].Print("v_3D[2]");
				printf("area_01 %g area_12 %g area_20 %g\n" , area_01,area_12,area_20);
				barycentric_coordinate.Print("barycentric_...");
				printf("\n");
			}
			face_list_[face_index] = CBaseModel::CFace(v_index[0],v_index[1],vertex_index);
			new_faces.push_back(pair<int,int>(face_index,face_list_.size()));
			face_list_.push_back(CBaseModel::CFace(v_index[1],v_index[2],vertex_index));

			new_faces.push_back(pair<int,int>(face_index,face_list_.size()));
			face_list_.push_back(CBaseModel::CFace(v_index[2],v_index[0],vertex_index));

			vertex_list_.push_back(vertex);
			source_index = vertex_list_.size() - 1;
			//return vertex_list_.size()-1;
		}else{
			new_faces.clear();
			

		}
	}
	void WriteToFile(const string& file_name ){
				cout << "write to file " << file_name << "\n";
		FILE* file_obj = fopen(file_name.c_str(),"w");
		if( file_obj == NULL ){
			printf("error in write to file %s\n", file_name.c_str());
			return ;
		}
		for(int j = 0; j < vertex_list_.size();++j){
			fprintf(file_obj,"v %lf %lf %lf\n" , vertex_list_[j].x,vertex_list_[j].y,vertex_list_[j].z);
		}
		for(int j = 0; j < face_list_.size();++j){
			fprintf(file_obj,"f %d %d %d\n" , face_list_[j][0]+1,face_list_[j][1]+1,face_list_[j][2]+1);
		}
		fclose(file_obj);
	}
};

class ModelForSubdivide{
public:
	vector<CPoint3D> vertex_list_;
	vector<CBaseModel::CFace> face_list_;
	vector< vector<int> > vertex_on_simple_edge;
	vector< vector<int> > vertex_in_face;
	vector< int >& origin_of_new_face;
	const CRichModel& model_;
	//vector<vector<CPoint3D>> 
public:
	ModelForSubdivide(const CRichModel& model, vector<int>&_origin_of_new_face) :model_(model), origin_of_new_face(_origin_of_new_face){
		vertex_list_.assign(model.m_Verts.begin(), model.m_Verts.end());
		face_list_.assign(model.m_Faces.begin(), model.m_Faces.end());
		vertex_on_simple_edge.resize(model_.GetNumOfTotalUndirectedEdges());
		vertex_in_face.resize(model_.GetNumOfFaces());
	}

	void addVertexOnFace(const CPoint3D& vertex, const int face_index) {
		if (face_index < 0) return;
			int v_index[3];
			int simple_edge_index[3];
			CPoint3D v_3D[3];
			for (int j = 0; j < 3; ++j){
				v_index[j] = face_list_[face_index][j];
				v_3D[j] = model_.Vert(v_index[j]);
				simple_edge_index[j] = model_.GetEdgeFromFace(face_index, j).indexOfSimpleEdge;
			}
			int vertex_index = vertex_list_.size();
			vertex_list_.push_back(vertex);
			vertex_in_face[face_index].push_back(vertex_index);
	}

	void addVertexOnEdge(const CPoint3D& vertex, int simple_edge_id) {
		int vertex_index = vertex_list_.size();
		vertex_list_.push_back(vertex);
		vector<int>& vertex_on_edge_id = vertex_on_simple_edge[simple_edge_id];
		vector<pair<int, double>> vertex_on_a_simple_edge;
		vertex_on_a_simple_edge.push_back(pair<int, double>(vertex_index,
			(model_.Vert(model_.SimpleEdge(simple_edge_id).v1) -
			vertex).Len()));
		bool need_to_insert = true;
		for (int i = 0; i < vertex_on_edge_id.size(); ++i){
			pair<int, double> temp_vertex;//first is index ,second is distance to v1 on simple edge
			temp_vertex.first = vertex_on_edge_id[i];
			temp_vertex.second = (model_.Vert(model_.SimpleEdge(simple_edge_id).v1) -
				vertex_list_[temp_vertex.first]).Len();
			//if ((temp_vertex.second - vertex_on_a_simple_edge[0].second) / model_.GetEdgeFromFace(face_index, on_edge_id).length < 1e-3) {
			//	need_to_insert = false;
			//	break;
			//}
			vertex_on_a_simple_edge.push_back(temp_vertex);
		}
		if (need_to_insert){
			sort(vertex_on_a_simple_edge.begin(), vertex_on_a_simple_edge.end(), compare_pair_second<std::less>());
			vertex_on_edge_id.clear();
			for (int i = 0; i < vertex_on_a_simple_edge.size(); ++i){
				vertex_on_edge_id.push_back(vertex_on_a_simple_edge[i].first);
			}
		}
		else{
			vertex_list_.pop_back();

		}
	}


	void addVertex(const CPoint3D& vertex, const int& face_index, int& source_index, int on_edge_id = -1){
		if (face_index >= 0){
			int v_index[3];
			int simple_edge_index[3];
			CPoint3D v_3D[3];
			for (int j = 0; j < 3; ++j){
				v_index[j] = face_list_[face_index][j];
				v_3D[j] = model_.Vert(v_index[j]);
				simple_edge_index[j] = model_.GetEdgeFromFace(face_index, j).indexOfSimpleEdge;
			}
			if (on_edge_id == -1){
				double barycentric_coordinate[3];
				calculateBarycentric(v_3D, vertex, barycentric_coordinate);
				for (int i = 0; i < 3; ++i){
					if (barycentric_coordinate[i] < 1e-6){
						on_edge_id = (i + 1) % 3;
						break;
					}
				}
			}
			int vertex_index = vertex_list_.size();
			source_index = vertex_index;
			vertex_list_.push_back(vertex);
			if (on_edge_id == -1){
				vertex_in_face[face_index].push_back(vertex_index);
			}
			else {
				vector<int>& vertex_on_edge_id = vertex_on_simple_edge[simple_edge_index[on_edge_id]];
				vector<pair<int, double>> vertex_on_a_simple_edge;
				vertex_on_a_simple_edge.push_back(pair<int, double>(vertex_index,
					(model_.Vert(model_.SimpleEdge(simple_edge_index[on_edge_id]).v1) -
					vertex).Len()));
				bool need_to_insert = true;
				for (int i = 0; i < vertex_on_edge_id.size(); ++i){
					pair<int, double> temp_vertex;//first is index ,second is distance to v1 on simple edge
					temp_vertex.first = vertex_on_edge_id[i];
					temp_vertex.second = (model_.Vert(model_.SimpleEdge(simple_edge_index[on_edge_id]).v1) -
						vertex_list_[temp_vertex.first]).Len();
					if ((temp_vertex.second - vertex_on_a_simple_edge[0].second) / model_.GetEdgeFromFace(face_index, on_edge_id).length < 1e-3){
						need_to_insert = false;
						break;
					}
					vertex_on_a_simple_edge.push_back(temp_vertex);
				}
				if (need_to_insert){
					sort(vertex_on_a_simple_edge.begin(), vertex_on_a_simple_edge.end(), compare_pair_second<std::less>());
					vertex_on_edge_id.clear();
					for (int i = 0; i < vertex_on_a_simple_edge.size(); ++i){
						vertex_on_edge_id.push_back(vertex_on_a_simple_edge[i].first);
					}
				}
				else{
					vertex_list_.pop_back();

				}
			}
		}
	}

	void subdivide(){
		double cnt = 0;
		printf("_______________cnt %d proption %lf \n", cnt, (double)cnt / (double)model_.GetNumOfFaces());

		for (int face_index = 0; face_index < model_.GetNumOfFaces(); ++face_index){
			int v_index[3];
			int simple_edge_index[3];
			CPoint3D v_3D[3];
			for (int j = 0; j < 3; ++j){
				v_index[j] = face_list_[face_index][j];
				v_3D[j] = model_.Vert(v_index[j]);
				simple_edge_index[j] = model_.GetEdgeFromFace(face_index, j).indexOfSimpleEdge;
			}
			if (vertex_in_face[face_index].empty() &&
				vertex_on_simple_edge[simple_edge_index[0]].empty() &&
				vertex_on_simple_edge[simple_edge_index[1]].empty() &&
				vertex_on_simple_edge[simple_edge_index[2]].empty()){
				continue;
			}
			if (vertex_on_simple_edge[simple_edge_index[0]].empty() &&
				vertex_on_simple_edge[simple_edge_index[1]].empty() &&
				vertex_on_simple_edge[simple_edge_index[2]].empty() &&
				vertex_in_face[face_index].size() == 1
				) {
				vector<int> input_point_index;
				vector<CBaseModel::CFace> triangle_list;
				for (int j = 0; j < 3; ++j){
					input_point_index.push_back(v_index[j]);
				}
				input_point_index.push_back(vertex_in_face[face_index][0]);
				triangle_list.push_back(CBaseModel::CFace(0, 1, 3));
				triangle_list.push_back(CBaseModel::CFace(1, 2, 3));
				triangle_list.push_back(CBaseModel::CFace(2, 0, 3));
				for (int j = 0; j < triangle_list.size(); ++j){
					triangle_list[j][0] = input_point_index[triangle_list[j][0]];
					triangle_list[j][1] = input_point_index[triangle_list[j][1]];
					triangle_list[j][2] = input_point_index[triangle_list[j][2]];
				}
				face_list_[face_index] = triangle_list[0];
				for (int j = 1; j < triangle_list.size(); ++j){
					origin_of_new_face.push_back(face_index);
					face_list_.push_back(triangle_list[j]);
				}
				continue;
			}

			Point2D v_2D[3];
			vector<Point2D> input_point_list_2D;
			vector<int> input_point_index;
			map3DTriangleTo2D(v_3D, v_2D);
			for (int j = 0; j < 3; ++j){
				input_point_list_2D.push_back(v_2D[j]);
				input_point_index.push_back(v_index[j]);
				if (model_.SimpleEdge(simple_edge_index[j]).v1 == v_index[j]){
					for (int k = 0; k < vertex_on_simple_edge[simple_edge_index[j]].size(); ++k){
						CPoint3D& temp_vertex_3d = vertex_list_[vertex_on_simple_edge[simple_edge_index[j]][k]];
						Point2D temp_vertex_2d = ptIn3DTriangleTo2D(v_3D, v_2D, temp_vertex_3d);
						input_point_list_2D.push_back(temp_vertex_2d);
						input_point_index.push_back(vertex_on_simple_edge[simple_edge_index[j]][k]);
					}
				}
				else {
					for (int k = vertex_on_simple_edge[simple_edge_index[j]].size() - 1; k >= 0; --k){
						CPoint3D& temp_vertex_3d = vertex_list_[vertex_on_simple_edge[simple_edge_index[j]][k]];
						Point2D temp_vertex_2d = ptIn3DTriangleTo2D(v_3D, v_2D, temp_vertex_3d);
						input_point_list_2D.push_back(temp_vertex_2d);
						input_point_index.push_back(vertex_on_simple_edge[simple_edge_index[j]][k]);
					}
				}
			}
			vector<vector<int>> input_segment_list;
			for (int j = 0; j < input_point_list_2D.size(); ++j){
				vector<int> segment(2);
				segment[0] = j;
				segment[1] = (j + 1) % input_point_list_2D.size();
				input_segment_list.push_back(segment);
			}
			for (int j = 0; j < vertex_in_face[face_index].size(); ++j){
				CPoint3D& temp_vertex_3d = vertex_list_[vertex_in_face[face_index][j]];
				Point2D temp_vertex_2d = ptIn3DTriangleTo2D(v_3D, v_2D, temp_vertex_3d);
				input_point_list_2D.push_back(temp_vertex_2d);
				input_point_index.push_back(vertex_in_face[face_index][j]);
			}
			vector<Point2D> output_point_list;
			vector<CBaseModel::CFace> triangle_list;
			vector<vector<int>> output_segment_list;
			triangulateSegLoop(input_point_list_2D, input_segment_list, output_point_list, triangle_list, output_segment_list, "pzQS");
			//assert(output_point_list.size() == input_point_list_2D.size());
			if (output_point_list.size() > input_point_list_2D.size()) {
				for (int i = input_point_list_2D.size(); i < output_point_list.size(); ++i) {
					input_point_index.push_back(vertex_list_.size());
					vertex_list_.push_back(pointIn2DTriangleTo3D(v_3D, v_2D, output_point_list[i]));
				}
			}
			for (int j = 0; j < triangle_list.size(); ++j) {
				int a = triangle_list[j][0];
				int b = triangle_list[j][1];
				int c = triangle_list[j][2];
				triangle_list[j][0] = input_point_index[triangle_list[j][0]];
				triangle_list[j][1] = input_point_index[triangle_list[j][1]];
				triangle_list[j][2] = input_point_index[triangle_list[j][2]];
				if (triangle_list[j][0] < 0 || triangle_list[j][1] < 0 || triangle_list[j][2] < 0 || 
					triangle_list[j][0] > 100000 || triangle_list[j][1] > 100000 || triangle_list[j][2] > 100000) {
					//printf("sz%d szoutput_point_list %d input_point_list_2D%d ", input_point_index.size(), output_point_list.size(), input_point_list_2D.size());
					//for (auto p : input_point_index) {
					//	printf(" %d", p);
					//}
					//printf("\n");
					//printf("inx %d %d %d\n", a, b, c);
					//printf("%d %d %d\n", triangle_list[j][0], triangle_list[j][1], triangle_list[j][2]);
				}
			}
			face_list_[face_index] = triangle_list[0];
			for (int j = 1; j < triangle_list.size(); ++j){
				origin_of_new_face.push_back(face_index);
				face_list_.push_back(triangle_list[j]);
			}
		}
	}
	void WriteToFile(const string& file_name){
		cout << "write to file " << file_name << "\n";
		FILE* file_obj = fopen(file_name.c_str(), "w");
		if (file_obj == NULL){
			printf("error in write to file %s\n", file_name.c_str());
			return;
		}
		for (int j = 0; j < vertex_list_.size(); ++j){
			fprintf(file_obj, "v %lf %lf %lf\n", vertex_list_[j].x, vertex_list_[j].y, vertex_list_[j].z);
		}
		for (int j = 0; j < face_list_.size(); ++j){
			fprintf(file_obj, "f %d %d %d\n", face_list_[j][0] + 1, face_list_[j][1] + 1, face_list_[j][2] + 1);
		}
		fclose(file_obj);
	}

};
#endif
