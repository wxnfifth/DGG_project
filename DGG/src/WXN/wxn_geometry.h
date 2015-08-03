#ifndef _WXN_GEOMETRY_
#define _WXN_GEOMETRY_
#include <cmath>
#include "ich\Point3D.h"
#include "ich\RichModel.h"
struct PointOnFace{
    int face_id_;
    CPoint3D v;
    //CPoint3D barycentric_coordinate;
    PointOnFace(){}
    PointOnFace(int _face_id , CPoint3D _v):face_id_(_face_id),v(_v){}
    void print(const string& s){
        fprintf(stderr,"PointOnFace %s: %d , %lf %lf %lf\n" , s.c_str() , face_id_ , v.x , v.y , v.z);
    }
};

struct Point2D{
	double x, y;
	Point2D(){x = 0.0; y = 0.0;};
	Point2D(double a, double b) : x(a), y(b){};
	Point2D operator+(const Point2D & o) const {return Point2D(x+o.x, y+o.y);};
	Point2D& operator +=(const Point2D& pt)
	{
		x += pt.x;
		y += pt.y;
		return *this;
	}
	Point2D& operator -=(const Point2D& pt)
	{
		x -= pt.x;
		y -= pt.y;
		return *this;
	}
    void print(const string& information="")const{
		printf( "%s: point2d %g %g\n" , information.c_str() , x , y);
	}
	Point2D& operator /=(double times)
	{
		x /= times;
		y /= times;
		return *this;
	}
	bool operator !=(const Point2D o) const{
		return fabs(x - o.x) > DOUBLE_EPSILON || fabs(y - o.y) > DOUBLE_EPSILON;
	}

	Point2D operator-(const Point2D & o) const {return Point2D(x-o.x, y-o.y);};
	Point2D operator*(const double d)const{
		return Point2D(x*d,y*d);
	}
	Point2D scale(double s)const {return Point2D(x*s, y*s); };
	double length()const { return sqrt(x*x+y*y); };
	double length2()const { return x*x+y*y; };
	//2D product returns a scalor.
	double crossProduct(const Point2D & o) const { return x*o.y - y*o.x; }
	double dotProduct(const Point2D & o) const {return x*o.x + y*o.y; }
    void normalize() {  
        double len = length();
        if( fabs(len) < DOUBLE_EPSILON ) return;
        x /= len;
        y /= len; 
    }
    double find_angle(const Point2D&  p) const {
        double cos_theta = this->dotProduct(p);
        double sin_theta = this->crossProduct(p);
        double angle = atan2(sin_theta , cos_theta);
        return angle;
    }
	
};

// Class template for computing the barycentric coordinates of a point
// relative to a given triangle. Typically, this class will be
// instantiated for either <float> or <double> type, depending
// on the desired precision.
//
// This code precomputes and stores data that speed up computing
// the barycentric coordinates. This is particularly useful when
// repeatedly evaluating different points on the same triangle.
//
// (C) 2002 Ivan Neulander

template <class T> class Barycentric {
  int _proj; // which dimension we discard

  T _v0[3];  // original position of point that we
             // moved to the origin

  T _v1[2];  // the other two points with one dimension discarded
  T _v2[2]; 

  T _bv1[2]; // precomputed precursor to barycentric
  T _bv2[2];

public:
  Barycentric() { Init(); }

  Barycentric(const T *v0, const T *v1, const T *v2) 
  { Init(v0, v1, v2); }

  static void Diff(const T *a, const T *b, T *res) {
    res[0] = a[0] - b[0];
    res[1] = a[1] - b[1];
    res[2] = a[2] - b[2];
  }

  void Init() {
    _proj = -1;
    _v0[0] = _v0[1] = _v0[2] = 0;
    _v1[0] =_v1[1] = 0;
    _v2[0] =_v2[1] = 0;
    _bv1[0] =_bv1[1] = 0;
    _bv2[0] =_bv2[1] = 0;
  }
  
  void Init(const T *iv0, const T *iv1, const T *iv2);
  // Call this to specify the 3 verts of the triangle for which
  // barycentric coordinates are to be computed.


  bool IsValid() const { return _proj >= 0; }
  // True iff a non-degenerate triangle has been specified
  // using constructor or Init().
  
  bool Eval(const T *pos, T *bary) {
    // Attempts to compute barycentric coords of point 'pos' relative
    // to the triangle specified using constructor or Init().  This
    // method works even if the point lies outside triangle; however
    // the point must always be coplanar with the triangle. Returns
    // true on success, false on failure (i.e. triangle is
    // degenerate or uninitialized).

    T p0, p1;
    switch(_proj) {
    case 0:
      p0 = pos[1] - _v0[1];
      p1 = pos[2] - _v0[2];
      break;
    case 1:
      p0 = pos[2] - _v0[2];
      p1 = pos[0] - _v0[0];
      break;
    case 2:
      p0 = pos[0] - _v0[0];
      p1 = pos[1] - _v0[1];
      break;
    default:
      return false;
    }

    bary[1] = _bv2[1]*p0 - _bv2[0]*p1;
    bary[2] = _bv1[0]*p1 - _bv1[1]*p0;
    bary[0] = 1 - bary[1] - bary[2];

    return true;
  }

  bool Test(const T *pos, T *bary, T min = 0, T max = 1) {
    // This is an optimized version of Eval(). It only computes
    // barycentric coordinates if point 'pos' lies within the
    // triangle. If so, it returns true. Otherwise it returns
    // false. As with Eval(), it assumes that 'pos' is coplanar with
    // triangle.

    T p0, p1;

    switch(_proj) {
    case 0:
      p0 = pos[1] - _v0[1];
      p1 = pos[2] - _v0[2];
      break;
    case 1:
      p0 = pos[2] - _v0[2];
      p1 = pos[0] - _v0[0];
      break;
    case 2:
      p0 = pos[0] - _v0[0];
      p1 = pos[1] - _v0[1];
      break;
    default:
      return false;
    }

    bary[1] = _bv2[1]*p0 - _bv2[0]*p1;
    if (bary[1] < min || bary[1] > max) return false;
 
    bary[2] = _bv1[0]*p1 - _bv1[1]*p0;
    if (bary[2] < min || bary[2] > max) return false;

    bary[0] = 1 - bary[1] - bary[2];
    return bary[0] >= min;
  }
};

template <class T> void Barycentric<T>::
 Init(const T *iv0, const T *iv1, const T *iv2) { 
  // Invoke this method before evaluating or testing any barycentric
  // coords. It reads in the verts (iv0, iv1, iv2) of the triangle to
  // be used as a basis for barycentric coordinates, and precomputes
  // needed information. Returns 'true' if triangle is valid, false
  // otherwise.
  
  T v1[3], v2[3];
  
  Init();
  
  // copy iv0
  _v0[0] = iv0[0];
  _v0[1] = iv0[1];
  _v0[2] = iv0[2];
  
  // translate everything to place iv0 at origin
  Diff(iv1, iv0, v1);
  Diff(iv2, iv0, v2);
  
  // The system we're solving is overdetermined; We solve using one of
  // the following projections: (y,z), (z,x), (x,y): namely, the one
  // that yields the biggest determinant. So we just need to decide
  // which dimension to discard.
  
  T det0 = v1[1]*v2[2] - v2[1]*v1[2];
  T det1 = v1[2]*v2[0] - v2[2]*v1[0];
  T det2 = v1[0]*v2[1] - v2[0]*v1[1];
  
  T adet0 = fabs(det0);
  T adet1 = fabs(det1);
  T adet2 = fabs(det2);
  
  // Pick proj according to the biggest among adet0, adet1, adet2.
  // Note that if they're all zero, then proj will be assigned 2, so
  // that's the only case where we need to check before dividing by
  // adet2.
  int proj = adet0 > adet1 ? 
    (adet0 > adet2 ? 0 : 2) : (adet1 > adet2 ? 1 : 2);
  
  T detInv;
  switch(proj) {
  case 0:
    detInv = 1/det0;
    _proj = 0;
    
    _bv1[0] = detInv*v1[1];
    _bv1[1] = detInv*v1[2];
    _bv2[0] = detInv*v2[1];
    _bv2[1] = detInv*v2[2];
    
    _v1[0] = v1[1];
    _v1[1] = v1[2];
     _v2[0] = v2[1];
    _v2[1] = v2[2];
    
    break;
    
  case 1:
    detInv = 1/det1;
    _proj = 1;
    
    _bv1[0] = detInv*v1[2];
    _bv1[1] = detInv*v1[0];
    _bv2[0] = detInv*v2[2];
    _bv2[1] = detInv*v2[0];

    _v1[0] = v1[2];
    _v1[1] = v1[0];
    _v2[0] = v2[2];
    _v2[1] = v2[0];
    
    break;
    
  case 2:
    if (adet2 == 0) return;
    
    detInv = 1/det2;
    _proj = 2;
    
    _bv1[0] = detInv*v1[0];
    _bv1[1] = detInv*v1[1];
    _bv2[0] = detInv*v2[0];
    _bv2[1] = detInv*v2[1];
    
    _v1[0] = v1[0];
    _v1[1] = v1[1];
     _v2[0] = v2[0];
    _v2[1] = v2[1];
    
    break;
  }
}

Point2D operator*(double times, const Point2D& pt);
Point2D operator/(const Point2D& pt,const double& times);

class FindCircleInGraph{
public:
	const vector<vector<int>>& graph_;
	vector<bool>visited_;
	vector<int> father_index_;
	vector<vector<int>> circles_;
	FindCircleInGraph(const vector<vector<int>>& _graph_):graph_(_graph_){
		visited_.resize(graph_.size());
		fill(visited_.begin(),visited_.end(),false);
		father_index_.resize(graph_.size());
		fill(father_index_.begin(),father_index_.end(),-1);
	}
	void executeAlgorithm(){
		for(int u = 0; u < graph_.size();++u){
			if( visited_[u] ) continue;
			dfs(u);
		}
	}
	void printCircle(){
		for(int i = 0; i < circles_.size();++i){
			printf("circle ");
			for(int j = 0; j < circles_[i].size();++j){
				printf("%d " , circles_[i][j]+1);
			}
 			printf("\n");
		}
	}
	void printGraph(){
		printf("graph size %d\n" , graph_.size());
		set<int> edges;
		for(int i = 0; i < graph_.size();++i){
			for(int j = 0; j < graph_[i].size();++j){
				int u = i;
				int v = graph_[i][j];
				if( u > v ){
					swap(u,v);
				}
				int id = u + v * graph_.size();
				if( edges.find(id) != edges.end() ) continue;
				edges.insert(id);
				printf("edge %d %d\n" , u +1, v +1);
			}
		}

	}
	void maxCircle(vector<int>& max_circle){
		int max_circle_size = 0;
		int max_circle_pos = 0;
		if( circles_.size() == 0 ) return;
		for(int i = 0; i < circles_.size();++i){
			if( circles_[i].size() >= max_circle_size ){
				max_circle_pos = i;
				max_circle_size = circles_[i].size();
			}
		}
		//printf("max %d\n" , max_circle_pos);
		max_circle.assign(circles_[max_circle_pos].begin(),circles_[max_circle_pos].end());
	}
    void dfs(int start_node) {
        visited_[start_node] = true;
        for(int i = 0; i < graph_[start_node].size();++i){
            int next_node = graph_[start_node][i];
            if( visited_[next_node] ){
                vector<int>circle;
                //printf("start ");
                if( father_index_[start_node] != next_node ){
                    //printf("find one back end %d %d\n",start_node+1,next_node+1);
                    int temp_node = start_node;
                    while( temp_node != -1 && temp_node != next_node ){
                        circle.push_back(temp_node);
                        //printf("%d " , temp_node+1);
                        temp_node = father_index_[temp_node];
                    }
                    //printf("%d " , temp_node+1); 
                    //printf("\n");
                    if( temp_node == -1 ){
                    }else{
                        circle.push_back(temp_node);
                        circles_.push_back(circle);
                    }
                }
                continue;
            }else{
                visited_[next_node] = true;
                father_index_[next_node] = start_node;
                dfs(next_node);
            }
        }
    }

};

struct TexturedFaces{
    struct Box{
        Point2D min_p , max_p;
        Box(){
            min_p = Point2D(1e10 , 1e10);
            max_p = Point2D(-1e10,-1e10);
        }
        void update(const Point2D& p) {
            min_p.x = std::min(min_p.x , p.x);
            min_p.y = std::min(min_p.y , p.y);
            max_p.x = std::max(max_p.x , p.x);
            max_p.y = std::max(max_p.y , p.y);
        }
        bool inBox(const Point2D& p) const {
            return min_p.x <= p.x &&
                   min_p.y <= p.y &&
                   max_p.x >= p.x &&
                   max_p.y >= p.y;
        }
    };
private:
    int source_index;
    double radius;
    double rotate_angle;
public:
    Point2D center_2d;
    PointOnFace center_on_surface;
    map<int, Point2D> planar_coordinates;
    vector<int> textured_faces;
    vector<Box> box_of_faces;
    const CRichModel* model_pt;
    TexturedFaces();
    TexturedFaces(const CRichModel* _model_pt , const int _source_index , const double _radius , const double _rotate_angle) 
        :model_pt(_model_pt) ,
        source_index(_source_index),
        radius(_radius),
        rotate_angle(_rotate_angle){}
    TexturedFaces(const CRichModel* _model_pt , const int _source_index , const double _rotate_angle) 
        :model_pt(_model_pt) ,
        source_index(_source_index),
        radius(-1.0),
        rotate_angle(_rotate_angle){}
    void setSourceIndex(int source){ source_index = source;}
    int getSourceIndex() const{ return source_index;}
    double getRadius() const {return radius;}
    void setRadius(const double _radius) { radius = _radius;}
    void setRotateAngle (const double _rotate_angle) { rotate_angle = _rotate_angle;}
    double getRotateAngle () const {return rotate_angle;}
    
    void calculateCoveredFaces();
    
    void mapPointToModel(const Point2D& p , PointOnFace& point_on_face) const ;
    void mapModelPoint2Plane(const PointOnFace& point_3D , Point2D& point_2D) const;
    void RenderFaces(const GLenum& mode) const;
    void dumpTexturedFaces() const;
};


struct TexturedPoint{
public:
    Point2D p_2d;
    PointOnFace p_3d;
public:
    TexturedPoint(){}
    TexturedPoint(const Point2D& _p_2d , const PointOnFace& _p_3d) :
        p_2d(_p_2d) , p_3d(_p_3d){}
    const Point2D& Point_2D()const{
        return p_2d;
    }
    const CPoint3D& Point_3D()const{
        return p_3d.v;
    }
};


struct TexturedFacesWithControlPoints :TexturedFaces
{
    struct ControlEdge{
        int e1 , e2;
        ControlEdge(){}
        ControlEdge(int _e1 , int _e2):
            e1(_e1) , e2(_e2) {}
    };
    vector<Point2D> origin_control_points;
    vector<ControlEdge> origin_control_edges;
    vector<TexturedPoint> control_points;
    vector<ControlEdge> control_edges;
    vector<TexturedPoint> temp_control_points;
    vector<ControlEdge> temp_control_edges;

    vector<vector<TexturedPoint>> control_edge_lines;
    //vector<vector<Point2D>> control_polygons;
    //vector<vector<PointOnFace>> control_polygons_3d;
    map<int, Point2D> deformed_coordinates;
    double control_point_radius;
    int select_control_id;

    TexturedFacesWithControlPoints(){}
    TexturedFacesWithControlPoints(const CRichModel* _model_pt , const int _source_index , const double _radius,const double _rotate_angle)
        :TexturedFaces(_model_pt , _source_index ,_radius,_rotate_angle) {
            select_control_id = -1;
            control_point_radius = 0.008;
    }

    void SetControlPoints(const vector<Point2D>& _control_points ,
                          const vector<int>& _edges);
    void addControlPoint(const PointOnFace& hit_point , const bool& flag_mouse_first_point);
    void endAddControlPoint();

    void renderControlPoints()const;
    void UpdateControlPoints();
    void updateDeformedCoordinates();
    void RenderFaces(const GLenum& mode) const;
    void selectControlPoint(const CPoint3D& hit_point);
    int  getClosestControlPoint(const CPoint3D&p);
    void  moveSelection(const PointOnFace& last_p ,
                        const PointOnFace& p);
};


#endif