#ifndef _EXACTGEODESIC_H_
#define _EXACTGEODESIC_H_
#include "Mesh.h"
#include <list>
#include <vector>
#include <map>
#include <time.h>
#include <cassert>

#include "PriorityQueue.h"

using namespace std;

namespace MMP_HalfEdge
{
	struct SurfacePoint
	{
		SurfacePoint() {}
		SurfacePoint(UINT f, Vector3D p) {faceIdx = f; pos = p;}
		UINT faceIdx;
		Vector3D pos;
	};

	//线源结构
	struct LineSource
	{
		LineSource() {}
		LineSource(UINT f, Vector3D pos, Vector3D color, UINT s)
		{
			points_on_plat.clear(); colors.clear();
			Face_idx = f; points_on_plat.push_back(pos); colors.push_back(color); srcId = s;
		}
		LineSource(const SurfacePoint& sp, UINT sId)
		{
			Face_idx = sp.faceIdx; points_on_plat.push_back(sp.pos); colors.push_back(Vector3D(1.0, 0.0, 0.0)); srcId = sId;
		}
		UINT Face_idx;												//面的编号
		vector<Vector3D> points_on_plat;
		vector<Vector3D> colors;
		UINT srcId;
	};

	//平面上直线Ax+By+Cz=0
	struct Line
	{
		double A;
		double B;
		double C;
	};

	struct GeodesicKeyNode
	{
		bool isVert;
		UINT vert_idx;
		UINT edge_idx;
		double x;
	};

	struct GeodesicPath
	{
		Vector3D src; UINT srcFace; UINT srcVert;
		Vector3D dst; UINT dstFace; UINT dstVert;
		list<GeodesicKeyNode> keypoints;
		list<unsigned> passedFaces;
	};

	struct VertOccupyInfo {
		bool isVert;
		unsigned elementIndex;
		double pos;
	};

	typedef list<Window> WindowList;

	class ExactGeodesic
	{
	public:
		ExactGeodesic();
		ExactGeodesic(ExactGeodesic* exactGeodesic);
		~ExactGeodesic();
	public:
		Priority_Queue windowPriorityQueue;
		CMesh* mesh;
		vector<UINT> vertex_src_idx;
		int winNum;
		vector<double> vertexDist;
		double maxDist;
		vector<LineSource> polySource;
		vector<LineSource> curve;

	public:
		GeodesicPath geodesicpath;
		bool pathbuilt;

		vector<WindowList> windowlist;
		vector <VertOccupyInfo> vertOccupyied;

	public:
		void AssignMesh(CMesh* _mesh);
		void AssignSrcVertex(UINT src_idx);
		void RemoveSrcVertex(UINT src_idx);

		void AddPolySourceVertex(LineSource s);
		//建立已选定源点的距离场
		void BuildGeodesicDistField(std::vector<SurfacePoint> *stopPoints = NULL);
		//建立从源点到点point的测地线
		double BuildGeodesicPath(const Vector3D& point, const UINT& face_idx, UINT vertex_id);

		double CalcGeodesicLenOnly(const Vector3D& point, const UINT& face_idx, UINT vertex_id, UINT &srcId);

		double CalcGeodesicDistToVertex(UINT vertex_id, UINT &srcId);

		void Clear();

		void ClearGeodesicPath();

		void FillInVertDist(map<UINT, vector<UINT> > &region2Faces);

		void SplitWindows(unsigned edgeIndex);

		void UpdateAroundVert(unsigned centerVert);

		//对点窗口win进行传播
		void PropogatePointWindow(const Window& win);
		//对线窗口win进行传播
		void PropogateLineWindow(const Window& win);

	private:

		//检查顶点到所在窗口伪源点是否有路径与所在三角形其他两边相交
		bool CheckVertexOnEdge(const Window& w, const double& pos);
		//为point在面face_idx上寻找窗口
		void SearchWindowOnFace(const UINT& face_idx, const Vector3D& point, double& pathlen, UINT& edge_idx, double& minpos);
		//寻找通过窗口回到源点最短的位置
		void FindShortestPathThroughWindow(WindowList::iterator win, double x2, double y2, double& pos);
		//对窗口win进行传播
		void Propogate(const Window& win);
		
		//添加一个伪源点pseudo_src_idx，其到源点vertex_src_idx的距离为sigma
		void AddPseudoSource(UINT pseudo_src_idx, double sigma);

		void AddPseudoSource(UINT face_idx, Vector3D position, double sigma, UINT srcId, UINT pseuSrcId);

		void AddPseudoSource(UINT pseudo_src_idx, double sigma, UINT* edges, UINT srcId, UINT pseuSrcId);
		//以线段为源，添加窗口
		void AddPseudoSource(LineSource ls);
		//为窗口win设置值
		bool BuildUpWindow(Window& win, UINT Edge_idx, double b1, double b2, double d1, double d2, double sigma, bool tao, Window_Type winType, UINT srcId, UINT pseuSrcId);
		//为两个窗口相交解一元二次方程；w1和w2是两个窗口的相交部分。解放在p1和p2中。返回true表示有解，返回false表示无解
		void SolveEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//为两个点窗口相交解一元二次方程
		void SolvePPEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//为一个点窗口和一个线窗口相交解一元二次方程
		void SolvePLEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//为两个线窗口相交解一元一次方程
		void SolveLLEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//解平面直线相交坐标
		bool SolvePlanarEquation(const Line& l1, const Line& l2, double& x, double& y);
		//判断边上的点p更“靠近”s0和s1哪个伪源点，true为s0，false为s1
		bool CloserToWindow(const Window& w0, const Window& w1, double p, double& pathlen0, double& pathlen1);
		bool CloserToWindow(const Window& w0, const Window& w1, double p);
		//将新窗口加入到其边上，处理相交的情况
		void AddWindow(Window& win);

		//将窗口win从leftb到rightb切分为一个窗口，其中leftb和rightb为所在边上的坐标(而非相对于窗口端点的坐标)
		void CutWindow(const Window& win, double leftb, double rightb, Window& subwin);
		//计算一个窗口的Frontier到源点的距离
		double CalcMinDist(const Window& win);

		bool StopConditionSatisfied(std::vector<SurfacePoint> *stopPoints);

		bool toLeft(Vector2D& v0, Vector2D& v1, Vector2D& v2);
		bool InTriangle(Vector2D& v0, Vector2D& v1, Vector2D& v2, Vector2D& v3);

		void ReverseWinToTwinEdge(Window &win);
	};

}
#endif