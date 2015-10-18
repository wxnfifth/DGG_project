#ifndef __YX_METRIC_H__
#define __YX_METRIC_H__

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#ifdef _MSC_VER
#include "windows.h"
#include "ICH\Point3D.h"
#endif

//#define IS_GPU_CODE

#ifndef CPU_AND_GPU
#ifdef IS_GPU_CODE
#define CPU_AND_GPU __host__ __device__ 
#else
#define CPU_AND_GPU
#endif
#endif

//#pragma warning(disable:4996)

#define IS_NAN(x) ((x)!=(x))
#define YX_MAX(x,y) ((x)>(y)?(x):(y))
#define YX_MIN(x,y) ((x)<(y)?(x):(y))

#ifdef IS_GPU_CODE
__device__ const double YX_PI =  3.14159265358979323846;
__device__ const double YX_E =   2.71828182845904523536;
__device__ const double EPS = 1e-12;
__device__ const double YX_2PI = 6.283185307179586476925;
#else
CPU_AND_GPU const double YX_PI =  3.14159265358979323846;
CPU_AND_GPU const double YX_E =   2.71828182845904523536;
CPU_AND_GPU const double EPS = 1e-12;
CPU_AND_GPU const double YX_2PI = 6.283185307179586476925;
#endif

#ifdef _MSC_VER
class YXTimer{
public:
	void start(){
        QueryPerformanceCounter(&t1);
    }
	double gettime(){
		QueryPerformanceCounter(&t2);
		QueryPerformanceFrequency(&f); 
		return (double)(t2.QuadPart - t1.QuadPart)/f.QuadPart;
	}
private:
	LARGE_INTEGER t1,t2, f;
};
#else
class YXTimer{
	clock_t t;
public:
	void start(){
		t = clock();
	}
	double gettime(){
		return double(clock() - t) / CLOCKS_PER_SEC;
	}
};
#endif

// 3D

struct YXPoint3D{
	double x, y, z;
	CPU_AND_GPU YXPoint3D(){};
	CPU_AND_GPU YXPoint3D(double a, double b) : x(a), y(b), z(0.0) { };
	CPU_AND_GPU YXPoint3D(double a, double b, double c) : x(a), y(b), z(c){ };
	CPU_AND_GPU YXPoint3D(const CPoint3D& p) :x(p.x), y(p.y), z(p.z) {};
	CPU_AND_GPU YXPoint3D operator+(const YXPoint3D & o) const { return YXPoint3D(x + o.x, y + o.y, z + o.z); };
	CPU_AND_GPU YXPoint3D operator-(const YXPoint3D & o) const { return YXPoint3D(x - o.x, y - o.y, z - o.z); };
	CPU_AND_GPU YXPoint3D operator*(double o) const { return YXPoint3D(x * o, y * o, z * o); };
	CPU_AND_GPU YXPoint3D scale(double s)const { return YXPoint3D(x*s, y*s, z*s); };
	CPU_AND_GPU double length()const { return sqrt(x*x + y*y + z*z); };
	CPU_AND_GPU double length2()const { return x*x + y*y + z*z; };
	CPU_AND_GPU YXPoint3D crossProduct(const YXPoint3D & o)const { return YXPoint3D(y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x); };
	CPU_AND_GPU double dotProduct(const YXPoint3D & o)const { return x*o.x + y*o.y + z*o.z; };
	CPU_AND_GPU CPoint3D toCPoint3D() { return CPoint3D(x, y, z); }
	CPU_AND_GPU YXPoint3D toUnitLength() {
		double len = length();
		x /= len; y /= len; z /= len;
		return *this;
	}
  

};

// 2D
struct YXPoint2D{
	double x, y;
	CPU_AND_GPU YXPoint2D(){};
	CPU_AND_GPU YXPoint2D(double a, double b) : x(a), y(b){ };
	CPU_AND_GPU YXPoint2D operator+(const YXPoint2D & o) const {return YXPoint2D(x+o.x, y+o.y);};
	CPU_AND_GPU YXPoint2D operator-(const YXPoint2D & o) const {return YXPoint2D(x-o.x, y-o.y);};
	CPU_AND_GPU YXPoint2D scale(double s)const {return YXPoint2D(x*s, y*s); };
	CPU_AND_GPU double length()const { return sqrt(x*x+y*y); };
	CPU_AND_GPU double length2()const { return x*x+y*y; };
	//2D product returns a scalor.
	CPU_AND_GPU double crossProduct(const YXPoint2D & o) const { return x*o.y - y*o.x; };
	CPU_AND_GPU double dotProduct(const YXPoint2D & o) const {return x*o.x + y*o.y; };
	CPU_AND_GPU YXPoint3D toPoint3D() const { return YXPoint3D(x, y, 0.0); };
};

// 4D
struct YXPoint4D{
	double x, y, z, w;
	CPU_AND_GPU YXPoint4D(){};
	CPU_AND_GPU YXPoint4D(double a, double b, double c, double d) : x(a), y(b), z(c), w(d) { };
	CPU_AND_GPU YXPoint4D operator+(const YXPoint4D & o) const {return YXPoint4D(x+o.x, y+o.y, z+o.z, w+o.w);};
	CPU_AND_GPU YXPoint4D operator-(const YXPoint4D & o) const {return YXPoint4D(x-o.x, y-o.y, z-o.z, w-o.w);};
	CPU_AND_GPU YXPoint4D scale(double s)const {return YXPoint4D(x*s, y*s, z*s, w*s); };
	CPU_AND_GPU double length()const { return sqrt(x*x+y*y+z*z+w*w); };
	CPU_AND_GPU double length2()const { return x*x+y*y+z*z+w*w; };
};

struct YXFace{
	int v[3];
	CPU_AND_GPU YXFace(){};
	CPU_AND_GPU YXFace(int v1, int v2, int v3) { v[0] = v1; v[1] = v2; v[2] = v3; };
	CPU_AND_GPU void reset(int a, int b, int c){
		v[0] = a; v[1] = b; v[2] = c;
	}
};

struct YXEdge{
	int v1;
	int v2; //actually, only v1 is enough.
	//int myEdgeIndex;
	int reverseEdgeIndex;
	double length;
};

struct YXTetrahedron{
	int v[4];
	CPU_AND_GPU void reset(int a, int b, int c, int d){
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
	CPU_AND_GPU static double randomDouble(unsigned int &seed){
		const unsigned int P = 1103515245;
		const unsigned int T = 12345;
		double r = 1.0;
		seed = (P * seed + T) / 65535 % 32768;
		r = (double)seed / 32768.0 + r / 32768.0;
		seed = (P * seed + T) / 65535 % 32768;
		r = (double)seed / 32768.0 + r / 32768.0;
		return r;
	}
	//compute triangle area.
	CPU_AND_GPU static double computeTriangleArea(double len1, double len2, double len3){
		double p = (len1 + len2 + len3) * 0.5;
		return sqrt( p * (p-len1) * (p-len2) * (p-len3) );
	}
	//convert 3D point to 2D 
	CPU_AND_GPU static YXPoint2D compute3Dto2D(double baseLength, double len1, double len2){
		//YXPoint2D res;
		//res.x = ((len1 * len1 - len2 * len2) / baseLength + baseLength) / 2.0;
		//res.y = -sqrt(len1 * len1 - res.x * res.x);
		//return res;

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
	CPU_AND_GPU static double projectScalar(const YXPoint3D & v, const YXPoint3D & w){
		double num = v.dotProduct(w);
		if(fabs(num) < 1e-6)return 0.0;
		return num / w.dotProduct(w);
	}
	CPU_AND_GPU static double projectScalar(const YXPoint2D & v, const YXPoint2D & w){
		double num = v.dotProduct(w);
		if(fabs(num) < 1e-9)return 0.0;
		return num / w.dotProduct(w);
	}
	//project point to line
	CPU_AND_GPU static YXPoint3D projectPointToLine(const YXPoint3D & A, const YXPoint3D & B, const YXPoint3D & P){
		return A + (B-A).scale( projectScalar(P-A, B-A) );
	}
	//unfold point P to plane ABC
	CPU_AND_GPU static YXPoint3D unfoldToPlane(const YXPoint3D & A, const YXPoint3D & B, const YXPoint3D & C, const YXPoint3D & P){
		YXPoint3D dir = projectPointToLine(A, B, C) - C;
		YXPoint3D pos = projectPointToLine(A, B, P);
		return dir.scale( sqrt( (pos-P).length2() / dir.length2() ) ) + pos;
	}
	//compute line intersection point. all points are 2D.
	//two segments: v1-->v2, v3-->v4, we will return the scalar t on v1-->v2, the result ineterseciont point is v1+(v2-v1)*t
	CPU_AND_GPU static double computeLineIntersection(const YXPoint2D & v1, const YXPoint2D & v2, const YXPoint2D &v3, const YXPoint2D &v4, int & res){
		if ( fabs((v2-v1).x * (v4-v3).y - (v2-v1).y * (v4-v3).x) < EPS ) {
			res = -1; //two segments are parallel
			return 0.0; 
		}
		if( fabs(v3.x - v4.x) < EPS ){ //the segment is parallel to y-axes
			res = 0;
			return (v3.x - v1.x) / (v2.x - v1.x);
		}
		double k1 = (v2.x - v1.x) / (v4.x - v3.x);
		double k2 = (v1.x - v3.x) / (v4.x - v3.x);
		res = 0;
		return ( (v4.y - v3.y) * k2 + v3.y - v1.y ) / (v2.y - v1.y - (v4.y - v3.y) * k1);
	}
	//WangRui's Integral Geodesic Approximation [Bowers 2010], Siggraph Asia
	CPU_AND_GPU static double approximateGeodesicDistance(const YXPoint3D & a, const YXPoint3D & b, const YXPoint3D & anormal, const YXPoint3D & bnormal){
		double euc_dis = (a - b).length();
		YXPoint3D v = (b - a).scale(1.0/euc_dis);
		double c1 = anormal.dotProduct(v);
		double c2 = bnormal.dotProduct(v);
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
	double scale;
	YXMesh3D() : nface(0), nvert(0){};
	double pointDistance( int a, int b ) const{
		return (vert[a]-vert[b]).length();
	}
	void computeBoundingBox(){
		double xmax, ymax, zmax;
		double xmin, ymin, zmin;
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
	void pushVertex(double x, double y, double z){
		vert.push_back( YXPoint3D(x, y, z) );
		nvert = (int) vert.size();
	}
	void pushFace(int v1, int v2, int v3){
		face.push_back(YXFace(v1, v2, v3));
		nface = (int) face.size();
	}

  int load_obj(const string& file_name) {
    vert.clear(); face.clear();
    nface = 0; nvert = 0;
    FILE* fp = fopen(file_name.c_str(), "r");
    if (fp == NULL) return -1;
    char buf[4096];
    YXPoint3D vx; 
    YXFace fc; 
    while(fgets(buf, 4096,fp) != NULL) {
      if (strncmp(buf, "v ", 2) == 0) {
        sscanf(buf, "v %lf %lf %lf", &vx.x, &vx.y, &vx.z);
        ++nvert;
        vert.push_back(vx);
      } else if(strncmp(buf, "f ", 2) == 0) {
        sscanf(buf, "f %d %d %d", &fc.v[0], &fc.v[1], &fc.v[2]);
        --fc.v[0]; --fc.v[1]; --fc.v[2]; ++nface;
        face.push_back(fc);
      }
    }
    fclose(fp);
    return 0;
  }

  int load_m(const string& file_name) {
    vert.clear(); face.clear(); nface = nvert = 0;
    FILE * fp = fopen(file_name.c_str(), "r"); if(fp == NULL) return -1;
    char buf[4096]; YXPoint3D vx; YXFace fc; int vid;
    while (fgets(buf, 4096, fp) != NULL) {
      if(strncmp(buf, "Vertex", 6) == 0){
        sscanf(buf, "Vertex %d %lf %lf %lf", &vid, &vx.x, &vx.y, &vx.z);
        ++nvert; vert.push_back(vx);
      }
      else if (strncmp(buf, "Face", 4) == 0) {
        sscanf(buf, "Face %d %d %d %d", &vid, &fc.v[0], &fc.v[1], &fc.v[2]);
        --fc.v[0]; --fc.v[1]; --fc.v[2]; ++nface;
        face.push_back(fc);
      }
    }
    fclose(fp); return 0;
  }

	int load(const char * filename) {
    string file_name_str = string(filename);
    int find_pos = file_name_str.rfind('.');
    if (find_pos == string::npos) {
      printf("Unkown type for %s\n", filename);
      return -1;
    }
    string post_fix = file_name_str.substr(find_pos+1);
    //PRINT_VAR(post_fix);
    if (post_fix == string("obj")) {
      return load_obj(filename);
    } else {
      return load_m(filename);
    }
	}
	void print(const char * filename){
		FILE * fp = fopen(filename, "w");
		for(int i = 0; i < nvert; ++i)
			fprintf(fp, "Vertex %d %.9lf %.9lf %.9lf\n", i+1, vert[i].x, vert[i].y, vert[i].z);
		for(int i = 0; i < nface; ++i)
			fprintf(fp, "Face %d %d %d %d\n", i+1, face[i].v[0]+1, face[i].v[1]+1, face[i].v[2]+1);
		fclose(fp);
	}
	void print_to_obj(const char * filename, YXPoint2D * texture = NULL){
		FILE * fp = fopen(filename, "w");
		static int signature = 0;
		fprintf(fp, "g 3D_object_%04d_%p\n", ++signature, this);
		for(int i = 0; i < nvert; ++i){
			fprintf(fp, "v %lf %lf %lf\n", vert[i].x, vert[i].y, vert[i].z);
		}
		
			for(int i = 0; i < nvert; ++i){
				if(texture != NULL)fprintf(fp, "vt %lf %lf\n", texture[i].x, texture[i].y);
				else fprintf(fp, "vt 0.5 0.5\n");
			}
		
		for(int i = 0; i < nface; ++i){
			fprintf(fp, "f %d/%d %d/%d %d/%d\n", face[i].v[0]+1, face[i].v[0]+1, 
				face[i].v[1]+1, face[i].v[1]+1, face[i].v[2]+1, face[i].v[2]+1);
		}
		fclose(fp);
	}
	void print_to_obj(const char * filename, double * texture){
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
		//const double shift[] = {1.000001, 0.999999, 1.0};
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
	double * angle;
	double * area;
	double * edgeAngle;
	YXPoint2D * ptToEdge2D;
	double totalArea;
	int nvert, nface, nedge;

	int * anyEdgeStartFromVert;

	YXMetric() : edge(nullptr), angle(nullptr), area(nullptr), edgeAngle(nullptr), 
		ptToEdge2D(nullptr), totalArea(0.0), nvert(0), nface(0), nedge(0), anyEdgeStartFromVert(nullptr){}

	~YXMetric() {
		if(edge != nullptr) delete [] edge;
		if(angle != nullptr) delete [] angle;
		if(area != nullptr) delete [] area;
		if(edgeAngle != nullptr) delete [] edgeAngle;
		if(ptToEdge2D != nullptr) delete [] ptToEdge2D;
		if(anyEdgeStartFromVert != nullptr) delete [] anyEdgeStartFromVert;
	}

	#ifdef IS_GPU_CODE
	void releaseGPU(){
		if(edge != NULL) cudaFree(edge);
		if(angle != NULL) cudaFree(angle);
		if(area != NULL) cudaFree(area);
		if(edgeAngle != NULL) cudaFree(edgeAngle);
		if(ptToEdge2D != NULL) cudaFree(ptToEdge2D);
		if(anyEdgeStartFromVert != NULL) cudaFree(anyEdgeStartFromVert);
	}
	#endif
	
	CPU_AND_GPU int getNextEdgeAroundFace(int edgeIndex) const{
		return edgeIndex - edgeIndex % 3 + (edgeIndex+1) % 3;
	}
	CPU_AND_GPU int getPrevEdgeAroundFace(int edgeIndex) const{
		return edgeIndex - edgeIndex % 3 + (edgeIndex+2) % 3;
	}

  CPU_AND_GPU int getEdgeFrom2Verts(int v1, int v2) const{
    int start_e = getAnyEdgeStartFromVert(v1);
    int e = start_e;
    int result = -1;
    //printf("edge size %d\n", nedge);

    while (e != -1) {
      //printf("e %d line 419\n", e);
      if (Edge(e).v2 == v2) {
        result = e;
        break;
      }
      e = getNextEdgeAroundVertex(e);
      if (e == start_e) break;
    }
    if (result >= 0) {
      assert(v1 == Edge(result).v1 && v2 == Edge(result).v2);
    }
    return result;
  }
	CPU_AND_GPU int getNextEdgeAroundVertex(int edgeIndex) const {
		if(edge[edgeIndex].reverseEdgeIndex == -1) return -1;
		return getNextEdgeAroundFace( edge[edgeIndex].reverseEdgeIndex );
	}
	CPU_AND_GPU int getPrevEdgeAroundVertex(int edgeIndex) const {
		return edge[ getPrevEdgeAroundFace(edgeIndex) ].reverseEdgeIndex;
	}
	CPU_AND_GPU int getAnyEdgeStartFromVert(int vertIndex)const{
		return anyEdgeStartFromVert[vertIndex];
	}
	CPU_AND_GPU const YXEdge & Edge(int edgeIndex) const{
		return edge[edgeIndex];
	}
	CPU_AND_GPU int faceIndex(int edgeIndex) const{
		return edgeIndex / 3;
	}
	CPU_AND_GPU bool isConcaveVertex(int vertIndex) const{
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

	#ifdef IS_GPU_CODE
	void buildOnGPU(const YXMesh3D & mesh){
		nvert = mesh.nvert;
		nface = mesh.nface;
		nedge = nface * 3;
		cudaMalloc(&edge, nface * 3 * sizeof(YXEdge));

	}
	#endif

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
		angle = new double [nvert];
		area = new double [nface];
		edgeAngle = new double[nedge];
		memset(angle, 0, sizeof(double) * nvert);
		totalArea = 0.0;
		for (int i = 0; i < nface; ++i) {
			int a = i * 3;
			int b = i * 3 + 1;
			int c = i * 3 + 2;
			area[i] = YXGeometry::computeTriangleArea(edge[a].length, edge[b].length, edge[c].length);
			totalArea += area[i];
      double a_len2 = edge[a].length * edge[a].length;
      double b_len2 = edge[b].length * edge[b].length;
      double c_len2 = edge[c].length * edge[c].length;
			//edgeAngle[b] = asin( 2.0 * area[i] / edge[a].length / edge[b].length );
      edgeAngle[b] = acos((a_len2 + b_len2 - c_len2) / (2.0 * edge[a].length * edge[b].length));
			angle[edge[b].v1] += edgeAngle[b]; 
			//edgeAngle[c] = asin( 2.0 * area[i] / edge[b].length / edge[c].length );
      edgeAngle[c] = acos((b_len2 + c_len2 - a_len2) / (2.0 * edge[b].length * edge[c].length));
			angle[edge[c].v1] += edgeAngle[c];
			//edgeAngle[a] = asin( 2.0 * area[i] / edge[c].length / edge[a].length );
      edgeAngle[a] = acos((c_len2 + a_len2 - b_len2) / (2.0 * edge[c].length * edge[a].length));
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
