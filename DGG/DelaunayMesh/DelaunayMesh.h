#include "ExactGeodesic.h"
#include <vector>
#include <set>

class DelaunayMesh
{
public:

	struct EdgeInPQ
	{
		unsigned edgeIndex;
		double priority;

		EdgeInPQ() {edgeIndex = -1; priority = 1e30;}
		EdgeInPQ(unsigned e, double p) {edgeIndex = e; priority = p;}

		bool operator< (const EdgeInPQ& right) const
		{
			return priority < right.priority;
		}
	};

	typedef std::multiset<EdgeInPQ> PriorityQueue;
	typedef PriorityQueue::iterator PQHandle;

	DelaunayMesh() 
	{ 
		mesh = NULL; redunMesh = NULL; originalMesh = NULL; delta = 1.0; 
		fourNonDelaunay = 0;
		threeNonDelaunay = 0;
		twoNonDelaunay = 0;
		oneNonDelaunay = 0;
		zeroNonDelaunay = 0;
		twoEdgesDelaunay = 0;
		oneEdgesDelaunay = 0;
		zeroEdgesDelaunay = 0;
		perpFootFallInFourInterval = 0; 
		FourIntervalMid = 0;
		perpFootFallInThreeInterval = 0;
		ThreeIntervalMid = 0;
		justPerpFoot = 0;
		surroundExistZeroDelaunayEdge = 0;
		surroundExistOneDelaunayEdge = 0;
		surroundExistTwoDelaunayEdge = 0;
		surroundExistThreeDelaunayEdge = 0;
		surroundExistFourDelaunayEdge = 0;
		flipMakeFourEdgeNonDelaunay = 0;
		flipMakeThreeEdgeNonDelaunay = 0;
		flipMakeTwoEdgeNonDelaunay = 0;
		flipMakeOneEdgeNonDelaunay = 0;
		flipMakeZeroEdgeNonDelaunay = 0;
		guanranteeDelaunayNum = 0;
	}
	~DelaunayMesh() 
	{ 
		mesh = NULL; 
		if (redunMesh) { delete redunMesh; redunMesh = NULL; }
		if (originalMesh) { delete originalMesh; originalMesh = NULL; }
	}

	void AssignMesh(CMesh *_mesh);
	void ConstructDelaunayMesh();
	void MergeTwoNeighPointOnPhysicalEdge();
	void MergeThreeNeighPointOnPhysicalEdge();
	void RemoveTooCloseVerts();

	unsigned NonDelaunayEdgeNum();

private:

	void Split(unsigned edgeIndex, double pos);

	void Flip(unsigned edgeIndex);

	unsigned ViolatingVertEdgePair();
	bool isEdgeLegal(CMesh *mesh, unsigned curEdge);
	double calcOpAngleSum(CMesh *mesh, unsigned curEdge);

	void FindSplitPos(unsigned violatedEdge, unsigned &edgeIndex, double&pos);

	void FindPerpendicularFootFarthest(unsigned curEdge, unsigned twinEdge, unsigned &edgeIndex, double &pos);

	void FindTheBestPosInBounder(double lowerBound, double upperBound, 
		unsigned curEdge, unsigned twinEdge, double len[5], double cosAngle1, 
		unsigned &edgeIndex, double &pos);
	void FindPerpendicularFoot(unsigned edgeIndex, double &pos);

	void FindPosNearestToMidPoint(unsigned edgeIndex, double &pos);
	void FindPosNearestToPerpendicular(unsigned edgeIndex, double &pos);
	void FindPosMakeSurroundEdgeLegal(unsigned baseEdge, unsigned surroundEdge, double &pos);

	void BounderIntersection(vector< pair<double, double> > &bounds, vector<unsigned> &exceptList, double &lowerBound, double &upperBound);
	bool BoundNotEmpty(double lowerBound, double upperBound, unsigned curEdge);

	void MergeTwoNeighPointOnPhysicalEdge(unsigned iEdge, double leftBound, double rightBound);
	void CalculateCircumCircle(vector<Vector2D> &points, Vector2D &center, double &R);
	void IntersectWithFeasibleIntervals(list< pair<double, double> > &feasibleIntervals, list< pair<double, double> > &intervals);
	void MergeThreeNeighPointOnPhysicalEdge(unsigned iEdge0, unsigned iEdge1, double leftBound, double rightBound);

	bool InInterval(double pos, pair<double, double> interval);

public:
	CMesh *redunMesh;
	string baseFileName;

	unsigned type;
	double delta;

	unsigned fourNonDelaunay;
	unsigned threeNonDelaunay;
	unsigned twoNonDelaunay;
	unsigned oneNonDelaunay;
	unsigned zeroNonDelaunay;
	unsigned twoEdgesDelaunay, oneEdgesDelaunay, zeroEdgesDelaunay;

	unsigned perpFootFallInFourInterval;
	unsigned FourIntervalMid;
	unsigned perpFootFallInThreeInterval;
	unsigned ThreeIntervalMid;
	unsigned justPerpFoot;

	unsigned mergeVertCnt;
	unsigned surroundExistZeroDelaunayEdge;
	unsigned surroundExistOneDelaunayEdge;
	unsigned surroundExistTwoDelaunayEdge;
	unsigned surroundExistThreeDelaunayEdge;
	unsigned surroundExistFourDelaunayEdge;

	unsigned flipMakeFourEdgeNonDelaunay;
	unsigned flipMakeThreeEdgeNonDelaunay;
	unsigned flipMakeTwoEdgeNonDelaunay;
	unsigned flipMakeOneEdgeNonDelaunay;
	unsigned flipMakeZeroEdgeNonDelaunay;

	unsigned guanranteeDelaunayNum;

private:
	CMesh *mesh;
	CMesh *originalMesh;
	double rhoV, rhoE;

	vector<unsigned> fatherEdge;
	vector<bool> originalEdgeHasBeenSplitted;
	vector< pair<unsigned, unsigned> > originalEdgeEndVerts;
	unsigned originalMeshEdgeNum;
	unsigned originalMeshVertNum;
};