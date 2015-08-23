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

	//��Դ�ṹ
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
		UINT Face_idx;												//��ı��
		vector<Vector3D> points_on_plat;
		vector<Vector3D> colors;
		UINT srcId;
	};

	//ƽ����ֱ��Ax+By+Cz=0
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
		//������ѡ��Դ��ľ��볡
		void BuildGeodesicDistField(std::vector<SurfacePoint> *stopPoints = NULL);
		//������Դ�㵽��point�Ĳ����
		double BuildGeodesicPath(const Vector3D& point, const UINT& face_idx, UINT vertex_id);

		double CalcGeodesicLenOnly(const Vector3D& point, const UINT& face_idx, UINT vertex_id, UINT &srcId);

		double CalcGeodesicDistToVertex(UINT vertex_id, UINT &srcId);

		void Clear();

		void ClearGeodesicPath();

		void FillInVertDist(map<UINT, vector<UINT> > &region2Faces);

		void SplitWindows(unsigned edgeIndex);

		void UpdateAroundVert(unsigned centerVert);

		//�Ե㴰��win���д���
		void PropogatePointWindow(const Window& win);
		//���ߴ���win���д���
		void PropogateLineWindow(const Window& win);

	private:

		//��鶥�㵽���ڴ���αԴ���Ƿ���·�����������������������ཻ
		bool CheckVertexOnEdge(const Window& w, const double& pos);
		//Ϊpoint����face_idx��Ѱ�Ҵ���
		void SearchWindowOnFace(const UINT& face_idx, const Vector3D& point, double& pathlen, UINT& edge_idx, double& minpos);
		//Ѱ��ͨ�����ڻص�Դ����̵�λ��
		void FindShortestPathThroughWindow(WindowList::iterator win, double x2, double y2, double& pos);
		//�Դ���win���д���
		void Propogate(const Window& win);
		
		//���һ��αԴ��pseudo_src_idx���䵽Դ��vertex_src_idx�ľ���Ϊsigma
		void AddPseudoSource(UINT pseudo_src_idx, double sigma);

		void AddPseudoSource(UINT face_idx, Vector3D position, double sigma, UINT srcId, UINT pseuSrcId);

		void AddPseudoSource(UINT pseudo_src_idx, double sigma, UINT* edges, UINT srcId, UINT pseuSrcId);
		//���߶�ΪԴ����Ӵ���
		void AddPseudoSource(LineSource ls);
		//Ϊ����win����ֵ
		bool BuildUpWindow(Window& win, UINT Edge_idx, double b1, double b2, double d1, double d2, double sigma, bool tao, Window_Type winType, UINT srcId, UINT pseuSrcId);
		//Ϊ���������ཻ��һԪ���η��̣�w1��w2���������ڵ��ཻ���֡������p1��p2�С�����true��ʾ�н⣬����false��ʾ�޽�
		void SolveEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//Ϊ�����㴰���ཻ��һԪ���η���
		void SolvePPEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//Ϊһ���㴰�ں�һ���ߴ����ཻ��һԪ���η���
		void SolvePLEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//Ϊ�����ߴ����ཻ��һԪһ�η���
		void SolveLLEquation(const Window& w0, const Window& w1, double& p1, double& p2);
		//��ƽ��ֱ���ཻ����
		bool SolvePlanarEquation(const Line& l1, const Line& l2, double& x, double& y);
		//�жϱ��ϵĵ�p����������s0��s1�ĸ�αԴ�㣬trueΪs0��falseΪs1
		bool CloserToWindow(const Window& w0, const Window& w1, double p, double& pathlen0, double& pathlen1);
		bool CloserToWindow(const Window& w0, const Window& w1, double p);
		//���´��ڼ��뵽����ϣ������ཻ�����
		void AddWindow(Window& win);

		//������win��leftb��rightb�з�Ϊһ�����ڣ�����leftb��rightbΪ���ڱ��ϵ�����(��������ڴ��ڶ˵������)
		void CutWindow(const Window& win, double leftb, double rightb, Window& subwin);
		//����һ�����ڵ�Frontier��Դ��ľ���
		double CalcMinDist(const Window& win);

		bool StopConditionSatisfied(std::vector<SurfacePoint> *stopPoints);

		bool toLeft(Vector2D& v0, Vector2D& v1, Vector2D& v2);
		bool InTriangle(Vector2D& v0, Vector2D& v1, Vector2D& v2, Vector2D& v3);

		void ReverseWinToTwinEdge(Window &win);
	};

}
#endif