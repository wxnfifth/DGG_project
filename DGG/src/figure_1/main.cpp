#include "DelaunayMesh.h"

#include <iostream>
#include <string>
#include <iomanip>
#include <ctime>

//#define TESTCODEENABLE

using namespace std;
using namespace MMP_HalfEdge;

CMesh *mesh = NULL;
ExactGeodesic *exactGeodesic = NULL;
DelaunayMesh *delaunayMesh = NULL;

int main(int argc, char **argv)
{
#ifndef TESTCODEENABLE
	if (argc < 3) 
	{
		cout << "USAGE: DelaunayMesh.exe [in.obj] [type]" << endl;
		cout << "[type]: " << endl; 
		cout << "\t0. Perpendicular" << endl;
		cout << "\t1. MidPoint" << endl;
		cout << "\t2: Dense sample" << endl;
		cout << "\t3: power(2, x) sample" << endl;
		cout << "\t4: Make the no new violating if possible, or skip the one with minimal opposite angle" << endl;
		cout << "\t5: Make the new violating minimal" << endl;
		return -1;
	}
	string baseFileName = argv[1];
	baseFileName = baseFileName.substr(baseFileName.rfind("\\")+1, baseFileName.rfind(".") - baseFileName.rfind("\\")-1);
	char meshOutFile[255];
	sprintf(meshOutFile, "%s.delaunay%d.obj", baseFileName.c_str(), atoi(argv[2]));

	mesh = new CMesh();
	if (!mesh->Load(argv[1]))
	{
		cout << "Cannot load mesh " << baseFileName << endl;
		return -2;
	}

	delaunayMesh = new DelaunayMesh();
	delaunayMesh->AssignMesh(mesh);
	delaunayMesh->type = atoi(argv[2]);
	delaunayMesh->baseFileName = baseFileName;

	cout << "Start to building delaunay mesh..." << endl;
	cout << "NonDelaunayEdgeNum: " << delaunayMesh->NonDelaunayEdgeNum() << endl;
	clock_t start = clock();
	delaunayMesh->ConstructDelaunayMesh();
	clock_t end = clock();
	cout << "Delaunay mesh built." << endl;
	cout << "Remove too close points..." << endl;
	/*delaunayMesh->RemoveTooCloseVerts();*/
	cout << "Done." << endl;

	if (delaunayMesh->type == 5)
	{
		cout << "Merge two neighbour verts on physical edges..." << endl;
		delaunayMesh->mergeVertCnt = 0;
		unsigned prevMergeVertCnt = -1;
		while (delaunayMesh->mergeVertCnt != prevMergeVertCnt)
		{
			prevMergeVertCnt = delaunayMesh->mergeVertCnt;
			delaunayMesh->MergeTwoNeighPointOnPhysicalEdge();
		}

		// 	cout << "Merge three neighbour verts on physical edges..." << endl;
		// 	delaunayMesh->mergeVertCnt = 0;
		// 	unsigned prevMergeVertCnt = -1;
		// 	while (delaunayMesh->mergeVertCnt != prevMergeVertCnt)
		// 	{
		// 		prevMergeVertCnt = delaunayMesh->mergeVertCnt;
		// 		delaunayMesh->MergeThreeNeighPointOnPhysicalEdge();
		// 	}

		cout << "Done." << endl;
	}
	cout << "Time: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;

	// output resulted mesh
// 	ofstream output(meshOutFile);
// 	CMesh *redunMesh = mesh;
// 	for (unsigned i = 0; i < redunMesh->m_nVertex; ++i)
// 		output << "v " << setprecision(20) << redunMesh->m_pVertex[i].m_vPosition.x << " " 
// 		<< redunMesh->m_pVertex[i].m_vPosition.y << " " 
// 		<< redunMesh->m_pVertex[i].m_vPosition.z << endl;
// 	for (unsigned i = 0; i < redunMesh->m_nFace; ++i)
// 		output << "f " << redunMesh->m_pFace[i].m_piVertex[0]+1 << " " << redunMesh->m_pFace[i].m_piVertex[1]+1 << " " << redunMesh->m_pFace[i].m_piVertex[2]+1 << endl;
	mesh->Save(meshOutFile);

	delete mesh;
	return 0;
#else
	mesh = new CMesh();
	mesh->Load(argv[1]);
	unsigned flipEdge = mesh->m_pVertex[1].m_piEdge[0];
	flipEdge = mesh->m_pEdge[flipEdge].m_iNextEdge;

	mesh->flip(flipEdge);

	ofstream output("test_flip.obj");
	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
		output << "v " << mesh->m_pVertex[i].m_vPosition.x * mesh->scaleD + mesh->origin.x << " " 
		<< mesh->m_pVertex[i].m_vPosition.y * mesh->scaleD + mesh->origin.y << " " 
		<< mesh->m_pVertex[i].m_vPosition.z * mesh->scaleD + mesh->origin.z << endl;
	for (unsigned i = 0; i < mesh->m_nFace; ++i)
		output << "f " << mesh->m_pFace[i].m_piVertex[0]+1 << " " << mesh->m_pFace[i].m_piVertex[1]+1 << " " << mesh->m_pFace[i].m_piVertex[2]+1 << endl;
	
	delete mesh;
	return 0;
#endif
}