///////////////////////class CMesh////////////////////////////////////
//		可以处理任意拓扑，任意多边形的二流型网格（2-manifode mesh)		//
//		可以读取一般的obj,smf,wrl等文件类型							//
//		根据点表和面表生成完整的连接关系								//
//																	//
//////////////////////////////////////////////////////////////////////
#ifndef MESH_H
#define MESH_H

#include "model.h"
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

typedef unsigned int UINT;
typedef std::list<Vector3D>	_VECTORLIST;
typedef std::list<UINT>		_UINTLIST;

#define MAX_VERTEX_PER_FACE 20

class CVertex
{
public:
	Vector3D	m_vPosition;		//点的坐标
	UINT*		m_piEdge;			//从该点发出的半边,要根据点的度数动态创建
	_UINTLIST	m_lEdgeList;		//用来构造m_piEdge的临时链表
	short		m_nValence;			//点的度数
	Vector3D	m_vNormal;			//顶点法向，由附近面法向平均得到
	bool		m_bIsBoundary;		//是否在边界上
	int			m_nCutValence;
	UINT		m_color;			//用于标记可行面或颜色等信息;
	bool m_bValid;

public:
	//constructions
	CVertex() { m_piEdge=NULL; m_nValence=0; m_nCutValence = 0; m_bIsBoundary = false; m_color = 0; m_bValid = true;}
	CVertex(double x,double y,double z) {m_vPosition=Vector3D(x,y,z);m_piEdge=NULL; m_nValence=0; m_bIsBoundary = false; m_nCutValence = 0; m_bValid = true;}
	CVertex(Vector3D v) {m_vPosition=v;m_piEdge=NULL; m_nValence=0; m_bIsBoundary = false; m_nCutValence = 0; m_bValid = true;}
	virtual ~CVertex();

	//operations
	CVertex& operator = (CVertex& v);

};

class CTexture
{
public:
	Vector2D m_vPosition;
public:
	CTexture() {m_vPosition = Vector2D(0, 0);}
	CTexture(double x, double y) {m_vPosition = Vector2D(x, y);}
	CTexture(Vector2D v) {m_vPosition = v;}
};

class CEdge
{
public:
	UINT	m_iVertex[2];		//边的两端点，Vertex0－>Vertex1

	UINT	m_iTwinEdge;		//与该边方向相反的另一条边，如果为-1则该边为边界
	UINT	m_iNextEdge;		//沿逆时针方向的下一条边
	UINT	m_iFace;			//该边所属的面，应该在它的左边
	UINT	m_color;			//用于标记可行面或颜色等信息;
	double  m_length;			//边长度;
	bool m_bValid;

public:		
	bool	m_bCut;
	int		m_nCutTag;

public:
	//constructions
	CEdge() {
		m_iVertex[0]=m_iVertex[1]=m_iTwinEdge=m_iNextEdge=m_iFace=-1; 
		m_bCut = false; m_nCutTag = 0;m_color = 0;m_length = 0; m_bValid = true;
	}
	CEdge(UINT iV0, UINT iV1) { m_iVertex[0]=iV0; m_iVertex[1]=iV1;m_iTwinEdge=m_iNextEdge=m_iFace=-1; m_bCut = false; m_nCutTag = 0; m_bValid = true;}
	virtual ~CEdge();

	//operations
	CEdge& operator = (const CEdge& e);
};

class CFace
{
public:
	short	m_nType;		//几边形
	UINT*	m_piVertex;		//所有点
	UINT*	m_piEdge;		//所有边
	double* m_pdAngle;
	Vector3D m_vNormal;		//法向
	Vector3D m_vMassPoint;	//法向
	double	m_dArea;		//面积
	UINT	m_color;		//用于标记可行面或颜色等信息;
	bool m_bValid;

public:
	//constructions
	CFace() {m_nType=0;m_piVertex=m_piEdge=NULL;m_vNormal=Vector3D(0.0,0.0,1.0);m_dArea=0.0;m_color = 0; m_bValid = true;}
	CFace(short s);
	virtual ~CFace();

	//operations
	void Create(short s);
	CFace& operator = (const CFace& f);
};

class CMesh :public CModel
{
public:
	UINT		m_nVertex;				//点数
	CVertex*	m_pVertex;				//点表

	std::vector<CTexture> m_pTexture;
	CTexture maxTex;

	UINT		m_nEdge;				//边数
	CEdge*		m_pEdge; 				//边表
	UINT		m_nFace;	 			//面数
	CFace*		m_pFace;				//面表
	UINT		m_nVertexCapacity;		//当前顶点列表容量
	UINT		m_nEdgeCapacity;		//当前边表容量
	UINT		m_nFaceCapacity;		//当前面表容量

	unsigned m_nValidVertNum;
	unsigned m_nValidFaceNum;
	unsigned m_nValidEdgeNum;

	std::vector<Vector3D> isolatedPoints;

	double *m_pAngles;

	std::string		filename;				//三角网格文件名
	double scaleD;
	Vector3D origin;

	//temp
	_UINTLIST m_lFocusEdge;
	_UINTLIST m_lFocusVertex;
	_UINTLIST m_lFocusFace;
	UINT	m_iPickedFace;
	UINT	m_iPickedEdge;
	UINT	m_iPickedVertex;

	bool	m_bClosed;

public:
	CMesh() {
		m_nVertex=m_nEdge=m_nFace=0;
		m_pVertex=NULL;m_pEdge=NULL;m_pFace=NULL;m_pAngles=NULL;
		m_iPickedFace=m_iPickedEdge=m_iPickedVertex=-1;
	}
	CMesh(CMesh* pMesh);
	virtual ~CMesh();

public:
	bool	Load(const char* sFileName);	// load from file
	bool	Save(const char* sFileName);	// save to file

	bool	construct();// construct connectivity

	CMesh*	clone();

	//将iEdge相邻的两个面（若为边界边则只有一个面）以iEdge进行细分。
	UINT	split(UINT iEdge, double posPercent = -1.0);

	//flip edge whose two adjacent facws are in the same plane
	void flip(unsigned iEdge);

	void collapse(unsigned iEdge, Vector3D newPos);

	void removeVert(unsigned iVert, unsigned startVert);

	void moveVertTo(unsigned iVert, Vector3D newPos);

	//计算所有边的边长，保存在CEdge->m_length中
	void	calcAllEdgeLength();
private:
	void	clear();
	bool	reConstruct();// construct connectivity from current mesh
	bool	loadFromSMF(const char* sFileName);
	bool	saveToSMF(const char* sFileName);

	//计算面i的法向量
	void	calFaceNormal(UINT i);
	//计算顶点i的法向量，由其相邻面的法向平均得到
	void	calVertexNormal(UINT i);
	//当现有点边面列表可能面临长度不够时，将空间加倍。
	void	expandCapacity();
	//计算每个顶点周围的角度
	void calcVertexAngle();
};

#endif //MESH_H