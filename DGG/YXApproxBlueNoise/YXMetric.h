#ifndef __YX_METRIC_H__
#define __YX_METRIC_H__

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "windows.h"

#pragma warning(disable:4996)

#define IS_NAN(x) ((x)!=(x))
#define YX_MAX(x,y) ((x)>(y)?(x):(y))
#define YX_MIN(x,y) ((x)<(y)?(x):(y))

const float YX_PI =  3.14159265358979323846;
const float YX_E =   2.71828182845904523536;
const float EPS = 1e-8;
const float YX_2PI = 6.283185307179586476925;

class YXTimer{
public:
	void start(){
        QueryPerformanceCounter(&t1);
    }
	float gettime(){
		QueryPerformanceCounter(&t2);
		QueryPerformanceFrequency(&f); 
		return (float)(t2.QuadPart - t1.QuadPart)/f.QuadPart;
	}
private:
	LARGE_INTEGER t1,t2, f;
};

// 3D
struct YXPoint3D{
	float x, y, z;
	YXPoint3D(){};
	YXPoint3D(float a, float b) : x(a), y(b), z(0.0) { };
	YXPoint3D(float a, float b, float c) : x(a), y(b), z(c){ };
	YXPoint3D operator+(const YXPoint3D & o) const {return YXPoint3D(x+o.x, y+o.y, z+o.z);};
	YXPoint3D operator-(const YXPoint3D & o) const {return YXPoint3D(x-o.x, y-o.y, z-o.z);};
	YXPoint3D scale(float s)const {return YXPoint3D(x*s, y*s, z*s); };
	float length()const { return sqrt(x*x+y*y+z*z); };
	float length2()const { return x*x+y*y+z*z; };
	YXPoint3D crossProduct(const YXPoint3D & o)const { return YXPoint3D(y*o.z-z*o.y, z*o.x-x*o.z, x*o.y-y*o.x); };
	float dotProduct(const YXPoint3D & o)const { return x*o.x + y*o.y + z*o.z; };
};

// 2D
struct YXPoint2D{
	float x, y;
	YXPoint2D(){};
	YXPoint2D(float a, float b) : x(a), y(b){ };
	YXPoint2D operator+(const YXPoint2D & o) const {return YXPoint2D(x+o.x, y+o.y);};
	YXPoint2D operator-(const YXPoint2D & o) const {return YXPoint2D(x-o.x, y-o.y);};
	YXPoint2D scale(float s)const {return YXPoint2D(x*s, y*s); };
	float length()const { return sqrt(x*x+y*y); };
	float length2()const { return x*x+y*y; };
	//2D product returns a scalor.
	float crossProduct(const YXPoint2D & o) const { return x*o.y - y*o.x; };
	float dotProduct(const YXPoint2D & o) const {return x*o.x + y*o.y; };
	YXPoint3D toPoint3D() const { return YXPoint3D(x, y, 0.0); };
};

// 4D
struct YXPoint4D{
	float x, y, z, w;
	YXPoint4D(){};
	YXPoint4D(float a, float b, float c, float d) : x(a), y(b), z(c), w(d) { };
	YXPoint4D operator+(const YXPoint4D & o) const {return YXPoint4D(x+o.x, y+o.y, z+o.z, w+o.w);};
	YXPoint4D operator-(const YXPoint4D & o) const {return YXPoint4D(x-o.x, y-o.y, z-o.z, w-o.w);};
	YXPoint4D scale(float s)const {return YXPoint4D(x*s, y*s, z*s, w*s); };
	float length()const { return sqrt(x*x+y*y+z*z+w*w); };
	float length2()const { return x*x+y*y+z*z+w*w; };
};

struct YXFace{
	int v[3];
	YXFace(){};
	YXFace(int v1, int v2, int v3) { v[0] = v1; v[1] = v2; v[2] = v3; };
	void reset(int a, int b, int c){
		v[0] = a; v[1] = b; v[2] = c;
	}
};

struct YXEdge{
	int v1;
	int v2; //actually, only v1 is enough.
	//int myEdgeIndex;
	int reverseEdgeIndex;
	float length;
};

struct YXTetrahedron{
	int v[4];
	void reset(int a, int b, int c, int d){
		v[0] = a; v[1] = b; v[2] = c; v[3] = d;
	}
};



class YXDuaMap{
	typedef unsigned long long LL;
	std::map<LL, int> mp;
public:
	void set(int key1, int key2, int value){
		if(key1 > key2) std::swap(key1, key2);
		LL l = ((LL)key1 << 32) + key2;
		mp[l] = value;
	}
	bool has(int key1, int key2){
		if(key1 > key2) std::swap(key1, key2);
		LL l = ((LL)key1 << 32) + key2;
		return mp.find(l) != mp.end();
	}
	int get(int key1, int key2){
		if(key1 > key2) std::swap(key1, key2);
		LL l = ((LL)key1 << 32) + key2;
		return mp[l];
	}
};

struct YXGeometry{
	static float randomfloat(unsigned int &seed){
		const unsigned int P = 1103515245;
		const unsigned int T = 12345;
		float r = 1.0;
		seed = (P * seed + T) / 65535 % 32768;
		r = (float)seed / 32768.0 + r / 32768.0;
		seed = (P * seed + T) / 65535 % 32768;
		r = (float)seed / 32768.0 + r / 32768.0;
		return r;
	}
	//compute triangle area.
	static float computeTriangleArea(float len1, float len2, float len3){
		float p = (len1 + len2 + len3) * 0.5;
		return sqrt( p * (p-len1) * (p-len2) * (p-len3) );
	}
	//convert 3D point to 2D 
	static YXPoint2D compute3Dto2D(float baseLength, float len1, float len2){
		YXPoint2D res;
		res.y = 2.0 * computeTriangleArea(baseLength, len1, len2) / baseLength;
		if(IS_NAN(res.y)) res.y = 0.0;
		res.x = sqrt(len1*len1 - res.y*res.y);
		if(len2*len2 > baseLength*baseLength + len1*len1){
			res.x = -res.x;
		}
		return res;
	}
	//project scalar: project v onto w, ==> scalar is s, thus the projection is s*w
	static float projectScalar(const YXPoint3D & v, const YXPoint3D & w){
		float num = v.dotProduct(w);
		if(fabs(num) < 1e-6)return 0.0;
		return num / w.dotProduct(w);
	}
	//project point to line
	static YXPoint3D projectPointToLine(const YXPoint3D & A, const YXPoint3D & B, const YXPoint3D & P){
		return A + (B-A).scale( projectScalar(P-A, B-A) );
	}
	//unfold point P to plane ABC
	static YXPoint3D unfoldToPlane(const YXPoint3D & A, const YXPoint3D & B, const YXPoint3D & C, const YXPoint3D & P){
		YXPoint3D dir = projectPointToLine(A, B, C) - C;
		YXPoint3D pos = projectPointToLine(A, B, P);
		return dir.scale( sqrt( (pos-P).length2() / dir.length2() ) ) + pos;
	}
	//compute line intersection point. all points are 2D.
	//两条边：v1-->v2, v3-->v4，返回在v1-->v2上的交点t，坐标为：v1+(v2-v1)*t
	static float computeLineIntersection(const YXPoint2D & v1, const YXPoint2D & v2, const YXPoint2D &v3, const YXPoint2D &v4, int & res){
		if ( fabs((v2-v1).x * (v4-v3).y - (v2-v1).y * (v4-v3).x) < EPS ) {
			res = -1; //两直线平行
			return 0.0; 
		}
		if( fabs(v3.x - v4.x) < EPS ){ //斜率不存在
			res = 0;
			return (v3.x - v1.x) / (v2.x - v1.x);
		}
		float k1 = (v2.x - v1.x) / (v4.x - v3.x);
		float k2 = (v1.x - v3.x) / (v4.x - v3.x);
		res = 0;
		return ( (v4.y - v3.y) * k2 + v3.y - v1.y ) / (v2.y - v1.y - (v4.y - v3.y) * k1);
	}
	//WangRui's Integral Geodesic Approximation [Bowers 2010], Siggraph Asia
	static float approximateGeodesicDistance(const YXPoint3D & a, const YXPoint3D & b, const YXPoint3D & anormal, const YXPoint3D & bnormal){
		float euc_dis = (a - b).length();
		YXPoint3D v = (b - a).scale(1.0/euc_dis);
		float c1 = anormal.dotProduct(v);
		float c2 = bnormal.dotProduct(v);
		if(fabs(c1-c2) < EPS) return euc_dis / sqrt(1.0 - c1 * c2);
		return (asin(c1) - asin(c2)) / (c1 - c2) * euc_dis;
	}
};

class YXMesh3D{
public:
	std::vector <YXPoint3D> vert;	// vertices
	std::vector <YXFace> face;		// faces
	std::vector <YXPoint3D> faceNormal;
	std::vector <YXPoint3D> vertNormal;
	int nface;	// number of faces
	int nvert;	// number of vertices
	float scale;
	YXMesh3D() : nvert(0), nface(0){};
	float pointDistance( int a, int b ) const{
		return (vert[a]-vert[b]).length();
	}
	void computeBoundingBox(){
		float xmax, ymax, zmax;
		float xmin, ymin, zmin;
		xmax = ymax = zmax = 1e-99;
		xmin = ymin = zmin = 1e99;
		for(int i = 0; i < vert.size(); ++i){
			xmax = YX_MAX(xmax, vert[i].x);
			ymax = YX_MAX(ymax, vert[i].y);
			zmax = YX_MAX(zmax, vert[i].z);
			xmin = YX_MIN(xmin, vert[i].x);
			ymin = YX_MIN(ymin, vert[i].y);
			zmin = YX_MIN(zmin, vert[i].z);
		}
		YXPoint3D tmp(xmax-xmin, ymax-ymin, zmax-zmin);
		YXPoint3D zero = tmp.scale(0.5);
		scale = tmp.length() / sqrt(3.0);
		for(int i = 0; i < vert.size(); ++i){
			vert[i] = (vert[i] - zero).scale(1.0/scale);
		}
	}
	void pushVertex(float x, float y, float z){
		vert.push_back( YXPoint3D(x, y, z) );
		nvert = (int) vert.size();
	}
	void pushFace(int v1, int v2, int v3){
		face.push_back(YXFace(v1, v2, v3));
		nface = (int) face.size();
	}
	int load(const char * filename){
		vert.clear(); face.clear(); nface = nvert = 0;
		FILE * fp = fopen(filename, "r"); if(fp == NULL) return -1;
		char buf[4096]; YXPoint3D vx; YXFace fc; int vid;
		while(fgets(buf, 4096, fp) != NULL){
			if(strncmp(buf, "Vertex", 6) == 0){
				sscanf(buf, "Vertex %d %f %f %f", &vid, &vx.x, &vx.y, &vx.z);
				++nvert; vert.push_back(vx);
			}
			else if(strncmp(buf, "Face", 4) == 0){
				sscanf(buf, "Face %d %d %d %d", &vid, &fc.v[0], &fc.v[1], &fc.v[2]);
				--fc.v[0]; --fc.v[1]; --fc.v[2]; ++nface;
				face.push_back(fc);
			}
		}
		fclose(fp); return 0;
	}
	//void print(const char * filename){
	//	FILE * fp = fopen(filename, "w");
	//	for(int i = 0; i < nvert; ++i)
	//		fprintf(fp, "Vertex %d %.9lf %.9lf %.9lf\n", i+1, vert[i].x, vert[i].y, vert[i].z);
	//	for(int i = 0; i < nface; ++i)
	//		fprintf(fp, "Face %d %d %d %d\n", i+1, face[i].v[0]+1, face[i].v[1]+1, face[i].v[2]+1);
	//	fclose(fp);
	//}
	void print_to_obj(const char * filename, YXPoint2D * texture = NULL){
		FILE * fp = fopen(filename, "w");
		static int signature = 0;
		fprintf(fp, "g 3D_object_%04d_%p\n", ++signature, this);
		for(int i = 0; i < nvert; ++i){
			fprintf(fp, "v %f %f %f\n", vert[i].x, vert[i].y, vert[i].z);
		}
		
			for(int i = 0; i < nvert; ++i){
				if(texture != NULL)fprintf(fp, "vt %f %f\n", texture[i].x, texture[i].y);
				else fprintf(fp, "vt 0.5 0.5\n");
			}
		
		for(int i = 0; i < nface; ++i){
			fprintf(fp, "f %d/%d %d/%d %d/%d\n", face[i].v[0]+1, face[i].v[0]+1, 
				face[i].v[1]+1, face[i].v[1]+1, face[i].v[2]+1, face[i].v[2]+1);
		}
		fclose(fp);
	}
	void print_to_obj(const char * filename, float * texture){
		YXPoint2D * p = new YXPoint2D[nvert];
		for(int i = 0; i < nvert; ++i) p[i].x = p[i].y = texture[i];
		print_to_obj(filename, p);
	}
	int get_edge_midpoint(YXDuaMap &mp, int v1, int v2){
		if(mp.has(v1, v2)){
			return mp.get(v1, v2);
		}
		int res = vert.size();
		mp.set(v1, v2, res);
		YXPoint3D v = (vert[v1] + vert[v2]).scale(0.5);
		//const float shift[] = {1.000001, 0.999999, 1.0};
		vert.push_back(v);
		return res;
	}
	void subdivide(){
		YXDuaMap mp;
		for(int i = 0; i < nface; ++i){
			int ov3 = get_edge_midpoint(mp, face[i].v[0], face[i].v[1]);
			int ov2 = get_edge_midpoint(mp, face[i].v[0], face[i].v[2]);
			int ov1 = get_edge_midpoint(mp, face[i].v[2], face[i].v[1]);
			YXFace f;
			f.reset(face[i].v[0], ov3, ov2);
			face.push_back(f);
			f.reset(face[i].v[1], ov1, ov3);
			face.push_back(f);
			f.reset(face[i].v[2], ov2, ov1);
			face.push_back(f);
			face[i].reset(ov1, ov2, ov3);
		}
		nface = (int)face.size();
		nvert = (int)vert.size();
	}
};

class YXMetric{
	typedef unsigned long long uint64_t;
public:
	YXEdge * edge;
	float * angle;
	float * area;
	float * edgeAngle;
	YXPoint2D * ptToEdge2D;
	float totalArea;
	int nvert, nface, nedge;

	int * anyEdgeStartFromVert;
	
	int getNextEdgeAroundFace(int edgeIndex) const{
		return edgeIndex - edgeIndex % 3 + (edgeIndex+1) % 3;
	}
	int getPrevEdgeAroundFace(int edgeIndex) const{
		return edgeIndex - edgeIndex % 3 + (edgeIndex+2) % 3;
	}
	int getNextEdgeAroundVertex(int edgeIndex) const {
		if(edge[edgeIndex].reverseEdgeIndex == -1) return -1;
		return getNextEdgeAroundFace( edge[edgeIndex].reverseEdgeIndex );
	}
	int getPrevEdgeAroundVertex(int edgeIndex) const {
		return edge[ getPrevEdgeAroundFace(edgeIndex) ].reverseEdgeIndex;
	}
	int getAnyEdgeStartFromVert(int vertIndex)const{
		return anyEdgeStartFromVert[vertIndex];
	}
	const YXEdge & Edge(int edgeIndex) const{
		return edge[edgeIndex];
	}
	int faceIndex(int edgeIndex) const{
		return edgeIndex / 3;
	}
	bool isConcaveVertex(int vertIndex) const{
		return angle[vertIndex] > YX_2PI;
	}
	void buildFromMesh(const YXMesh3D & mesh){
		nvert = mesh.nvert;
		nface = mesh.nface;
		nedge = nface * 3;
		edge = new YXEdge[nface*3];
		anyEdgeStartFromVert = new int [nvert];
		for(int i = 0; i < nvert; ++i){
			anyEdgeStartFromVert[i] = -1;
		}
		for(int i = 0; i < nface; ++i){
			for(int j = 0; j < 3; ++j){
				int idx = i * 3 + j;
				int v1 = mesh.face[i].v[j];
				int v2 = mesh.face[i].v[(j+1)%3];
				if(anyEdgeStartFromVert[v1] < 0) anyEdgeStartFromVert[v1] = idx;
				edge[idx].v1 = v1;
				edge[idx].v2 = v2; //optional.
				edge[idx].length = mesh.pointDistance( v1, v2 );
				//edge[idx].myEdgeIndex = idx; //optional.
			}
		}
		buildReverseEdge();
		computeAnglesForVertex();
		computePointToEdge2D();
	}
	void assignNormalsForMesh(YXMesh3D & mesh) const{
		mesh.faceNormal.resize(nface);
		mesh.vertNormal.resize(nvert);
		for(int i = 0; i < nface; ++i){
			const int * v = mesh.face[i].v;
			mesh.faceNormal[i] = (mesh.vert[v[0]] - mesh.vert[v[1]]).crossProduct(mesh.vert[v[1]] - mesh.vert[v[2]]);
		}
		for(int i = 0; i < nvert; ++i){
			mesh.vertNormal[i] = YXPoint3D(0.0, 0.0, 0.0);
			int e = getAnyEdgeStartFromVert(i);
			int e_idx = e;
			while(1){
				mesh.vertNormal[i] = mesh.vertNormal[i] + mesh.faceNormal[ faceIndex(e_idx) ];
				e_idx = getNextEdgeAroundVertex(e_idx);
				if(e_idx == e || e_idx == -1) break;
			}
			mesh.vertNormal[i] = mesh.vertNormal[i].scale( 1.0 / mesh.vertNormal[i].length() );
			
		}
		for(int i = 0; i < nface; ++i){
			mesh.faceNormal[i] = mesh.faceNormal[i].scale( 1.0 / mesh.faceNormal[i].length() );
		}
	}
protected:
	void buildReverseEdge(){
		//build Reverse Edge
		std::map< uint64_t, int > mp;
		std::map< uint64_t, int > :: iterator itr;
		uint64_t hash;
		int a, b;
		for(int i = 0; i < nedge; ++i){
			a = edge[i].v1;
			b = edge[i].v2;
			if(edge[i].v1 > edge[i].v2) std::swap(a, b);
			hash = ( (uint64_t)a << (uint64_t)32 ) + (uint64_t) b;
			itr = mp.find(hash);
			if(itr != mp.end()){
				int reverseIndex = itr->second;
				edge[i].reverseEdgeIndex = reverseIndex;
				edge[reverseIndex].reverseEdgeIndex = i;
				mp.erase(itr);
			}else{
				mp[hash] = i;
				edge[i].reverseEdgeIndex = -1;
			}
		}
	}
	void computeAnglesForVertex(){
		angle = new float [nvert];
		area = new float [nface];
		edgeAngle = new float[nedge];
		memset(angle, 0, sizeof(float) * nvert);
		totalArea = 0.0;
		for(int i = 0; i < nface; ++i){
			int a = i * 3;
			int b = i * 3 + 1;
			int c = i * 3 + 2;
			area[i] = YXGeometry::computeTriangleArea(edge[a].length, edge[b].length, edge[c].length);
			totalArea += area[i];
			edgeAngle[b] = asin( 2.0 * area[i] / edge[a].length / edge[b].length );
			angle[edge[b].v1] += edgeAngle[b];
			edgeAngle[c] = asin( 2.0 * area[i] / edge[b].length / edge[c].length );
			angle[edge[c].v1] += edgeAngle[c];
			edgeAngle[a] = asin( 2.0 * area[i] / edge[c].length / edge[a].length );
			angle[edge[a].v1] += edgeAngle[a];
		}
	}
	void computePointToEdge2D(){
		ptToEdge2D = new YXPoint2D [nedge];
		for(int i = 0; i < nface; ++i){
			int a = i * 3;
			int b = i * 3 + 1;
			int c = i * 3 + 2;
			ptToEdge2D[a] = YXGeometry::compute3Dto2D(edge[a].length, edge[c].length, edge[b].length);
			ptToEdge2D[b] = YXGeometry::compute3Dto2D(edge[b].length, edge[a].length, edge[c].length);
			ptToEdge2D[c] = YXGeometry::compute3Dto2D(edge[c].length, edge[b].length, edge[a].length);
		}
	}
};




#endif