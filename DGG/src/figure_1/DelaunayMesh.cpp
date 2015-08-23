#include "DelaunayMesh.h"
#include <stack>
#include <set>
#include <iostream>
#include <assert.h>
#include <algorithm>

using namespace MMP_HalfEdge;
using namespace std;

//#define OUTPUTINNERMESH
//#define DIRECTFOOT

void DelaunayMesh::AssignMesh(CMesh *_mesh)
{
	mesh = _mesh;
	originalMesh = new CMesh(mesh);

	//min edge length and min angle
	double Lmin = 1e30, Thetamin = 1e30, Lmax = -1e30;
	unsigned minEdge = -1, maxEdge = -1;
	for (unsigned i = 0; i < mesh->m_nEdge; ++i)
	{
		minEdge = mesh->m_pEdge[i].m_length < Lmin ? i : minEdge;
		maxEdge = mesh->m_pEdge[i].m_length > Lmax ? i : maxEdge;
		Lmin = mesh->m_pEdge[minEdge].m_length;
		Lmax = mesh->m_pEdge[maxEdge].m_length;
	}
	for (unsigned i = 0; i < mesh->m_nFace; ++i)
		for (unsigned j = 0; j < mesh->m_pFace[i].m_nType; ++j)
			Thetamin = mesh->m_pFace[i].m_pdAngle[j] < Thetamin ? mesh->m_pFace[i].m_pdAngle[j] : Thetamin;

	rhoV = min((Lmin*sin(Thetamin)) / (0.5 + sin(Thetamin)), Lmin / 2.0);
	rhoE = 2.0 * rhoV * sin(Thetamin);

	unsigned samplePointNumPreCalc = 0;
	for (unsigned i = 0; i < mesh->m_nEdge; ++i)
	{
		unsigned twinEdge = mesh->m_pEdge[i].m_iTwinEdge;
		if (i > twinEdge) continue;
		if (mesh->m_pEdge[i].m_length > 2 * rhoV)
		{
			samplePointNumPreCalc += 2;
			unsigned curSampleCandidateNum = (mesh->m_pEdge[i].m_length - 2 * rhoV) / rhoE;
			samplePointNumPreCalc += curSampleCandidateNum;
		}
		else if (mesh->m_pEdge[i].m_length > rhoV)
			++samplePointNumPreCalc;
	}

	cout << "Min Edge Length: " << Lmin << 
		" (" << mesh->m_pEdge[minEdge].m_iVertex[0] << ", " << mesh->m_pEdge[minEdge].m_iVertex[1] << ")" << endl;
	cout << "Max Edge Length: " << Lmax << 
		" (" << mesh->m_pEdge[maxEdge].m_iVertex[0] << ", " << mesh->m_pEdge[maxEdge].m_iVertex[1] << ")" << endl;

	cout << "Min sinTheta: " << sin(Thetamin) << endl;
	cout << "RhoV: " << rhoV << endl;
	cout << "RhoE: " << rhoE << endl;
	cout << "Predicted sample point number: " << samplePointNumPreCalc << endl;

	originalMeshEdgeNum = mesh->m_nEdge;
	originalMeshVertNum = mesh->m_nVertex;

	fatherEdge.resize(mesh->m_nEdge);
	for (unsigned i = 0; i < fatherEdge.size(); ++i) fatherEdge[i] = i;

	originalEdgeHasBeenSplitted.resize(mesh->m_nEdge, false);

	originalEdgeEndVerts.resize(mesh->m_nEdge);
	for (unsigned i = 0; i < originalEdgeEndVerts.size(); ++i)
		originalEdgeEndVerts[i] = make_pair(mesh->m_pEdge[i].m_iVertex[0], mesh->m_pEdge[i].m_iVertex[1]);
}

void DelaunayMesh::ConstructDelaunayMesh()
{
	unsigned iteration = 0, samplePointNum = 0;
	PriorityQueue edgeStack;
	vector<PQHandle> edgeInStack;

	edgeInStack.resize(mesh->m_nEdge, edgeStack.end());
	for (unsigned i = 0; i < mesh->m_nEdge; ++i)
	{
		if (i < mesh->m_pEdge[i].m_iTwinEdge)				// boundary edges are push into the stack
		{
			double opSumAngle = calcOpAngleSum(mesh, i);
			if (opSumAngle < 0) edgeInStack[i] = edgeStack.insert(EdgeInPQ(i, opSumAngle));
		}
	}

	while (!edgeStack.empty())
	{
// 		if (iteration % 1000 == 0)
// 		{
// 			unsigned violatingPair = ViolatingVertEdgePair();
// 			cout << "\rIteration " << iteration << "\tViolating Count: " <<  violatingPair << "\t";
// 			if (violatingPair == 0) break;
// 			/*system("pause");*/
// #ifdef OUTPUTINNERMESH
// 			cout << "Current iteration done. Output mesh? (y/n)" << endl;
// 			char ans = 'y';
// 			/*cin >> ans;*/
// 			if (ans == 'y')
// 			{
// 				char fileName[255];
// 				sprintf(fileName, "internalMesh%d.obj", 0);
// 				ofstream output(fileName);
// 				for (unsigned i = 0; i < mesh->m_nVertex; ++i)
// 				{
// 					output << "v " << mesh->m_pVertex[i].m_vPosition.x << " "
// 						<< mesh->m_pVertex[i].m_vPosition.y << " "
// 						<< mesh->m_pVertex[i].m_vPosition.z << endl;
// 				}
// 				for (unsigned i = 0; i < mesh->m_nFace; ++i)
// 				{
// 					output << "f " << mesh->m_pFace[i].m_piVertex[0] + 1 << " "
// 						<< mesh->m_pFace[i].m_piVertex[1] + 1 << " "
// 						<< mesh->m_pFace[i].m_piVertex[2] + 1 << endl;
// 				}
// 				output.close();
// 				system("pause");
// 			}
// #endif
// 		}
		++iteration;

		PQHandle iter = edgeStack.begin();
		unsigned curEdge = edgeStack.begin()->edgeIndex;
		edgeStack.erase(iter);
		edgeInStack[curEdge] = edgeStack.end();
		unsigned twinEdge = mesh->m_pEdge[curEdge].m_iTwinEdge;
		unsigned face0 = mesh->m_pEdge[curEdge].m_iFace;
		unsigned face1 = -1;
		if (twinEdge != -1) face1 = mesh->m_pEdge[twinEdge].m_iFace;
		if (isEdgeLegal(mesh, curEdge)) continue;
		unsigned surroundingEdge[4];
		surroundingEdge[0] = mesh->m_pEdge[curEdge].m_iNextEdge;
		surroundingEdge[1] = mesh->m_pEdge[surroundingEdge[0]].m_iNextEdge;
		if (twinEdge == -1)
		{
			surroundingEdge[2] = -1;
			surroundingEdge[3] = -1;
		}
		else
		{
			surroundingEdge[2] = mesh->m_pEdge[twinEdge].m_iNextEdge;
			surroundingEdge[3] = mesh->m_pEdge[surroundingEdge[2]].m_iNextEdge;
		}

		Vector3D n0 = mesh->m_pFace[face0].m_vNormal;
		n0.normalize();
		Vector3D n1(0, 0, 0);
		if (face1 != -1) { n1 = mesh->m_pFace[face1].m_vNormal; n1.normalize(); }

		unsigned e[4] = {-1};
		e[0] = mesh->m_pEdge[curEdge].m_iNextEdge;
		e[1] = mesh->m_pEdge[e[0]].m_iNextEdge;
		if (twinEdge != -1)
		{
			e[2] = mesh->m_pEdge[twinEdge].m_iNextEdge;
			e[3] = mesh->m_pEdge[e[2]].m_iNextEdge;
		}

		if (fatherEdge[curEdge] >= originalMeshEdgeNum && 
			twinEdge != -1 && fatherEdge[twinEdge] >= originalMeshEdgeNum
			/*|| EQUALZERO(n0 * n1 - 1)*/)
		{
			Flip(curEdge);

// 			unsigned nonDelaunayEdgesAround = 0;
// 			for (unsigned i = 0; i < 4; ++i)
// 			{
// 				if(e[i] == -1) continue;
// 				if (!isEdgeLegal(mesh, e[i])) ++nonDelaunayEdgesAround;
// 			}
// 
// 			switch (nonDelaunayEdgesAround)
// 			{
// 			case 0: ++flipMakeZeroEdgeNonDelaunay; break;
// 			case 1: ++flipMakeOneEdgeNonDelaunay; break;
// 			case 2: ++flipMakeTwoEdgeNonDelaunay; break;
// 			case 3: ++flipMakeThreeEdgeNonDelaunay; break;
// 			case 4: ++flipMakeFourEdgeNonDelaunay; break;
// 			}
		}
		else
		{
			++samplePointNum;
			unsigned splitEdge; double splitPos;

			FindSplitPos(curEdge, splitEdge, splitPos);
			/*if (splitPos < 0.0) continue;*/
			Split(splitEdge, splitPos);

			edgeInStack.resize(mesh->m_nEdge, edgeStack.end());
			double opAngleSum = calcOpAngleSum(mesh, curEdge);
			if (opAngleSum < 0) edgeInStack[curEdge] = edgeStack.insert(EdgeInPQ(curEdge, opAngleSum));
			if (twinEdge != -1)
			{
				if (edgeInStack[twinEdge] != edgeStack.end())
				{ edgeStack.erase(edgeInStack[twinEdge]); edgeInStack[twinEdge] = edgeStack.end(); }
				opAngleSum = calcOpAngleSum(mesh, twinEdge);
				if (opAngleSum < 0) edgeInStack[twinEdge] = edgeStack.insert(EdgeInPQ(twinEdge, opAngleSum));
			}

			// statistic
// 			unsigned nonDelaunayEdgesAround = 0;
// 			for (unsigned i = 0; i < 4; ++i)
// 			{
// 				if(e[i] == -1) continue;
// 				if (!isEdgeLegal(mesh, e[i])) ++nonDelaunayEdgesAround;
// 			}
// 
// 			switch (nonDelaunayEdgesAround)
// 			{
// 			case 0: ++zeroNonDelaunay; break;
// 			case 1: ++oneNonDelaunay; break;
// 			case 2: ++twoNonDelaunay; break;
// 			case 3: ++threeNonDelaunay; break;
// 			case 4: ++fourNonDelaunay; break;
// 			default: break;
// 			}
// 
// 			bool curEdgeDelaunay = isEdgeLegal(mesh, curEdge);
// 			bool twinEdgeDelaunay;
// 			if (twinEdge != -1) twinEdgeDelaunay = isEdgeLegal(mesh, twinEdge);
// 			else twinEdgeDelaunay = isEdgeLegal(mesh, mesh->m_nEdge-3);
// 
// 			if (curEdgeDelaunay && twinEdgeDelaunay)
// 				++twoEdgesDelaunay;
// 			else if (curEdgeDelaunay && !twinEdgeDelaunay || 
// 				!curEdgeDelaunay && twinEdgeDelaunay)
// 				++oneEdgesDelaunay;
// 			else
// 				++zeroEdgesDelaunay;
		}

		for (unsigned i = 0; i < 4; ++i) 
		{
			if (surroundingEdge[i] == -1) continue;
			if (surroundingEdge[i] > mesh->m_pEdge[surroundingEdge[i]].m_iTwinEdge) 
				surroundingEdge[i] = mesh->m_pEdge[surroundingEdge[i]].m_iTwinEdge;
			if (edgeInStack[surroundingEdge[i]] != edgeStack.end()) 
			{ edgeStack.erase(edgeInStack[surroundingEdge[i]]); edgeInStack[surroundingEdge[i]] = edgeStack.end(); }
			double opAngleSum = calcOpAngleSum(mesh, surroundingEdge[i]);
			if (opAngleSum < 0) edgeInStack[surroundingEdge[i]] = edgeStack.insert(EdgeInPQ(surroundingEdge[i], opAngleSum));
		}
	}
	cout << endl;
	cout << "Total iteration: " << iteration << endl;
	cout << "Total Sample Number: " << samplePointNum << endl << endl;
// 	cout << "Of all the sampling points: " << endl;
// 	cout << "\t " << fourNonDelaunay << " points make four edges around non-Delaunay." << endl;
// 	cout << "\t " << threeNonDelaunay << " points make three edges around non-Delaunay." << endl;
// 	cout << "\t " << twoNonDelaunay << " points make two edges around non-Delaunay." << endl;
// 	cout << "\t " << oneNonDelaunay << " points make one edges around non-Delaunay." << endl;
// 	cout << "\t " << zeroNonDelaunay << " points make zero edges around non-Delaunay." << endl;
// 	cout << "Of all the sampling points: " << endl;
// 	cout << "\t " << twoEdgesDelaunay << " points destroy the current non-Delaunay edge ." << endl;
// 	cout << "\t " << oneEdgesDelaunay << " points left on current non-Delaunay edge." << endl;
// 	cout << "\t " << zeroEdgesDelaunay << " points make an extra non-Delaunay edge." << endl;
// 
// 	cout << endl;
// 	cout << "Four interval is non-empty and perpendicular foot fall in it: " << perpFootFallInFourInterval << endl;
// 	cout << "Four interval is non-empty but perpendicular foot doesn't fall in it: " << FourIntervalMid << endl;
// 	cout << "Three interval is non-empty and perpendicular foot fall in it: " << perpFootFallInThreeInterval << endl;
// 	cout << "Three interval is non-empty but perpendicular foot doesn't fall in it: " << ThreeIntervalMid << endl;
// 	cout << "Three interval is even empty so just take perpendicular foot: " << justPerpFoot << endl;
// 
// 	cout << endl;
// 	cout << "Surrounding edges existing zero non-Delaunay: " << surroundExistZeroDelaunayEdge << endl;
// 	cout << "Surrounding edges existing one non-Delaunay: " << surroundExistOneDelaunayEdge << endl;
// 	cout << "Surrounding edges existing two non-Delaunay: " << surroundExistTwoDelaunayEdge << endl;
// 	cout << "Surrounding edges existing three non-Delaunay: " << surroundExistThreeDelaunayEdge << endl;
// 	cout << "Surrounding edges existing four non-Delaunay: " << surroundExistFourDelaunayEdge << endl;
// 
// 	cout << endl;
// 	cout << "Of all flipping: " << endl;
// 	cout << "Make 4 surrounding edges non-Delaunay: " << flipMakeFourEdgeNonDelaunay << endl;
// 	cout << "Make 3 surrounding edges non-Delaunay: " << flipMakeThreeEdgeNonDelaunay << endl;
// 	cout << "Make 2 surrounding edges non-Delaunay: " << flipMakeTwoEdgeNonDelaunay << endl;
// 	cout << "Make 1 surrounding edges non-Delaunay: " << flipMakeOneEdgeNonDelaunay << endl;
// 	cout << "Make 0 surrounding edges non-Delaunay: " << flipMakeZeroEdgeNonDelaunay << endl;
// 
// 	cout << endl;
// 	cout << "Minimized the Delaunay edges opAngleSum: " << guanranteeDelaunayNum << endl;
}

void DelaunayMesh::MergeTwoNeighPointOnPhysicalEdge()
{
	for (unsigned i = 0; i < originalMesh->m_nEdge; ++i)
	{
		/*if (mergeVertCnt >= 1) break;*/
		unsigned v0 = originalMesh->m_pEdge[i].m_iVertex[0];
		unsigned v1 = originalMesh->m_pEdge[i].m_iVertex[1];

		unsigned startEdge = -1;
		unsigned vi = -1, vj = -1, prevVi = -1, nextVj = -1;
		for (unsigned j = 0; j < mesh->m_pVertex[v0].m_nValence; ++j)
		{
			unsigned curEdge = mesh->m_pVertex[v0].m_piEdge[j];
			if (fatherEdge[curEdge] == i) 
			{
				startEdge = curEdge;
				break;
			}
		}
		if (startEdge == -1) continue;
		
		while (mesh->m_pEdge[startEdge].m_iVertex[1] != v1)
		{
			prevVi = vi;
			vi = mesh->m_pEdge[startEdge].m_iVertex[0];
			vj = mesh->m_pEdge[startEdge].m_iVertex[1];

			unsigned nextEdge = -1;
			for (unsigned j = 0; j < mesh->m_pVertex[vj].m_nValence; ++j)
			{
				nextEdge = mesh->m_pVertex[vj].m_piEdge[j];
				if (fatherEdge[nextEdge] == i) break;
			}

			nextVj = mesh->m_pEdge[nextEdge].m_iVertex[1];

			if (vi != v0 && vj != v1)
			{
				double leftBound = (mesh->m_pVertex[v0].m_vPosition - mesh->m_pVertex[prevVi].m_vPosition).length();
				double rightBound = (mesh->m_pVertex[v0].m_vPosition - mesh->m_pVertex[nextVj].m_vPosition).length();
				MergeTwoNeighPointOnPhysicalEdge(startEdge, leftBound, rightBound);
			}

			startEdge = nextEdge;
		}
	}

	cout << "Merge vertices number: " << mergeVertCnt << endl;
}

void DelaunayMesh::MergeThreeNeighPointOnPhysicalEdge()
{
	for (unsigned i = 0; i < originalMesh->m_nEdge; ++i)
	{
		/*if (mergeVertCnt >= 1) break;*/
		if (i > originalMesh->m_pEdge[i].m_iTwinEdge) continue;
		unsigned v0 = originalMesh->m_pEdge[i].m_iVertex[0];
		unsigned v1 = originalMesh->m_pEdge[i].m_iVertex[1];

		unsigned startEdge = -1;
		unsigned vi = -1, vj = -1, vk = -1, prevVi = -1, nextVj = -1;
		for (unsigned j = 0; j < mesh->m_pVertex[v0].m_nValence; ++j)
		{
			unsigned curEdge = mesh->m_pVertex[v0].m_piEdge[j];
			if (fatherEdge[curEdge] == i) 
			{
				startEdge = curEdge;
				break;
			}
		}
		if (startEdge == -1) continue;

		while (mesh->m_pEdge[startEdge].m_iVertex[1] != v1)
		{
			prevVi = vi;
			vi = mesh->m_pEdge[startEdge].m_iVertex[0];
			vj = mesh->m_pEdge[startEdge].m_iVertex[1];

			unsigned nextEdge = -1;
			for (unsigned j = 0; j < mesh->m_pVertex[vj].m_nValence; ++j)
			{
				nextEdge = mesh->m_pVertex[vj].m_piEdge[j];
				if (fatherEdge[nextEdge] == i) break;
			}
			vk = mesh->m_pEdge[nextEdge].m_iVertex[1];

			unsigned nextNextEdge = -1;
			for (unsigned j = 0; j < mesh->m_pVertex[vk].m_nValence; ++j)
			{
				unsigned curEdge = mesh->m_pVertex[vk].m_piEdge[j];
				if (fatherEdge[curEdge] == i) 
				{
					nextNextEdge = curEdge;
					break;
				}
			}
			if (nextNextEdge == -1) break;

			nextVj = mesh->m_pEdge[nextNextEdge].m_iVertex[1];

			if (vi != v0 && vk != v1)
			{
				double leftBound = (mesh->m_pVertex[v0].m_vPosition - mesh->m_pVertex[prevVi].m_vPosition).length();
				double rightBound = (mesh->m_pVertex[v0].m_vPosition - mesh->m_pVertex[nextVj].m_vPosition).length();
				MergeThreeNeighPointOnPhysicalEdge(startEdge, nextEdge, leftBound, rightBound);
			}

			startEdge = nextNextEdge;
		}
	}

	cout << "Merge vertices number: " << mergeVertCnt << endl;
}

void DelaunayMesh::RemoveTooCloseVerts()
{
	std::map<Vector3D, UINT> redunVertex;
	std::vector<UINT> realIndex;
	UINT vertexIndex = 0;

	vector<Vector3D> verts;
	vector<unsigned> faces;

	redunMesh = new CMesh();
	redunMesh->scaleD = mesh->scaleD;
	redunMesh->origin = mesh->origin;

	string buf;
	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
	{
		map<Vector3D, UINT>::iterator iter = redunVertex.find(mesh->m_pVertex[i].m_vPosition);
		if (iter == redunVertex.end())
		{
			redunVertex[mesh->m_pVertex[i].m_vPosition] = vertexIndex;
			realIndex.push_back(vertexIndex++);
			verts.push_back(mesh->m_pVertex[i].m_vPosition);
		}
		else
			realIndex.push_back(iter->second);
	}

	for (unsigned i = 0; i < mesh->m_nFace; ++i)
	{
		unsigned faceVerts[3];
		for (unsigned j = 0; j < 3; ++j) faceVerts[j] = mesh->m_pFace[i].m_piVertex[j];

		if (realIndex[faceVerts[0]] == realIndex[faceVerts[1]] || realIndex[faceVerts[1]] == realIndex[faceVerts[2]] || realIndex[faceVerts[2]] == realIndex[faceVerts[0]]) continue;
		faces.push_back(realIndex[faceVerts[0]]); faces.push_back(realIndex[faceVerts[1]]); faces.push_back(realIndex[faceVerts[2]]); 
	}

	redunMesh->m_nVertex = redunVertex.size();
	redunMesh->m_nFace = faces.size() / 3;

	redunMesh->m_pVertex = new CVertex[redunMesh->m_nVertex];
	redunMesh->m_pFace = new CFace[redunMesh->m_nFace];

	for (unsigned i = 0; i < redunMesh->m_nVertex; ++i)
		redunMesh->m_pVertex[i].m_vPosition = verts[i];
	for (unsigned i = 0; i < redunMesh->m_nFace; ++i)
	{
		redunMesh->m_pFace[i].Create(3);
		redunMesh->m_pFace[i].m_piVertex[0] = faces[i*3];
		redunMesh->m_pFace[i].m_piVertex[1] = faces[i*3+1];
		redunMesh->m_pFace[i].m_piVertex[2] = faces[i*3+2];
	}
}

unsigned DelaunayMesh::NonDelaunayEdgeNum()
{
	unsigned nonDelaunayEdgeNum = 0;
	for (unsigned i = 0; i < mesh->m_nEdge; ++i)
	{
		if (i > mesh->m_pEdge[i].m_iTwinEdge) continue;
		if (isEdgeLegal(mesh, i)) continue;
		++nonDelaunayEdgeNum;
	}
	return nonDelaunayEdgeNum;
}

void DelaunayMesh::Split(unsigned edgeIndex, double pos)
{
	unsigned twinEdgeIndex = mesh->m_pEdge[edgeIndex].m_iTwinEdge;

	unsigned prevEdgeNum = mesh->m_nEdge;

	mesh->split(edgeIndex, pos);

	unsigned addEdgeNum = mesh->m_nEdge - prevEdgeNum;
	fatherEdge.resize(mesh->m_nEdge);
	if (addEdgeNum == 6)
	{
		fatherEdge[prevEdgeNum] = fatherEdge[twinEdgeIndex];
		fatherEdge[prevEdgeNum+1] = fatherEdge[edgeIndex];
		for (unsigned i = 2; i < 6; ++i) fatherEdge[prevEdgeNum+i] = prevEdgeNum+i;
		originalEdgeHasBeenSplitted[fatherEdge[edgeIndex]] = true;
		originalEdgeHasBeenSplitted[fatherEdge[twinEdgeIndex]] = true;
	}
	else if (addEdgeNum == 3)
	{
		fatherEdge[prevEdgeNum] = fatherEdge[edgeIndex];
		for (unsigned i = 1; i < 3; ++i) fatherEdge[prevEdgeNum+i] = prevEdgeNum+i;
		originalEdgeHasBeenSplitted[fatherEdge[edgeIndex]] = true;
	}
}

void DelaunayMesh::Flip(unsigned edgeIndex)
{
	mesh->flip(edgeIndex);
}

unsigned DelaunayMesh::ViolatingVertEdgePair()
{
	unsigned violatingCnt = 0;
	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
	{
		unsigned violatingVert = i;
		for (unsigned j = 0; j < mesh->m_pVertex[i].m_nValence; ++j)
		{
			unsigned opEdge = mesh->m_pEdge[mesh->m_pVertex[i].m_piEdge[j]].m_iNextEdge;
			unsigned twinEdge = mesh->m_pEdge[opEdge].m_iTwinEdge;

			double edgeLen = mesh->m_pEdge[opEdge].m_length;
			double a, b, c, cot1, cot2;
			double cos1, sin1, cos2, sin2;
			a = edgeLen;
			b = mesh->m_pEdge[mesh->m_pEdge[opEdge].m_iNextEdge].m_length;
			c = mesh->m_pEdge[mesh->m_pVertex[i].m_piEdge[j]].m_length;
			/*cos1 = (b*b+c*c-a*a) / (2*b*c); sin1 = sqrt(fabs(1-cos1*cos1)); cot1 = cos1 / sin1;*/
			cot1 = (b*b+c*c-a*a) / sqrt((b+c+a)*(b+c-a)*(a+b-c)*(a-b+c));

			if (twinEdge == -1) 
			{
				if (cot1 < 0) ++violatingVert;
				continue;
			}

			b = mesh->m_pEdge[mesh->m_pEdge[twinEdge].m_iNextEdge].m_length;
			c = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[twinEdge].m_iNextEdge].m_iNextEdge].m_length;
			/*cos2 = (b*b+c*c-a*a) / (2*b*c); sin2 = sqrt(fabs(1-cos2*cos2)); cot2 = cos2 / sin2;*/
			cot2 = (b*b+c*c-a*a) / sqrt((b+c+a)*(b+c-a)*(a+b-c)*(a-b+c));

			if (cot1 + cot2 < 0)  ++ violatingCnt;
		}
	}

	return violatingCnt;
}

bool DelaunayMesh::isEdgeLegal(CMesh *mesh, unsigned curEdge)
{
	double a = mesh->m_pEdge[curEdge].m_length;
	double b = mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iNextEdge].m_length;
	double c = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iNextEdge].m_iNextEdge].m_length;
	double cotAlpha1 = (b*b+c*c-a*a) / sqrt((b+c+a)*(b+c-a)*(a+b-c)*(a-b+c));

	unsigned twinEdge = mesh->m_pEdge[curEdge].m_iTwinEdge;

	if (twinEdge == -1) 
		return cotAlpha1 >= 0;
	
	b = mesh->m_pEdge[mesh->m_pEdge[twinEdge].m_iNextEdge].m_length;
	c = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[twinEdge].m_iNextEdge].m_iNextEdge].m_length;
	double cotAlpha2 = (b*b+c*c-a*a) / sqrt((b+c+a)*(b+c-a)*(a+b-c)*(a-b+c));

	if (_fpclass(cotAlpha1) == _FPCLASS_PINF || _fpclass(cotAlpha2) == _FPCLASS_PINF) return true;
	if (cotAlpha1 != cotAlpha1 || cotAlpha2 != cotAlpha2) return true;
	return cotAlpha1 + cotAlpha2 >= 0;
}

double DelaunayMesh::calcOpAngleSum(CMesh *mesh, unsigned curEdge)
{
	double a = mesh->m_pEdge[curEdge].m_length;
	double b = mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iNextEdge].m_length;
	double c = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iNextEdge].m_iNextEdge].m_length;
	double cotAlpha1 = (b*b+c*c-a*a) / sqrt((b+c+a)*(b+c-a)*(a+b-c)*(a-b+c));

	unsigned twinEdge = mesh->m_pEdge[curEdge].m_iTwinEdge;

	if (twinEdge == -1) 
		return cotAlpha1;

	b = mesh->m_pEdge[mesh->m_pEdge[twinEdge].m_iNextEdge].m_length;
	c = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[twinEdge].m_iNextEdge].m_iNextEdge].m_length;
	double cotAlpha2 = (b*b+c*c-a*a) / sqrt((b+c+a)*(b+c-a)*(a+b-c)*(a-b+c));

	return cotAlpha1 + cotAlpha2;
}

void DelaunayMesh::FindSplitPos(unsigned violatedEdge, unsigned &edgeIndex, double&pos)
{
	unsigned curEdge = violatedEdge;
	unsigned twinEdge = mesh->m_pEdge[curEdge].m_iTwinEdge;

	if (type == 0)
	{
		double pos1 = 1e30, pos2 = 1e30;
		double edgeLen = mesh->m_pEdge[curEdge].m_length;
		FindPerpendicularFoot(curEdge, pos1);
		FindPerpendicularFoot(twinEdge, pos2);

		double toBound1 = min(pos1/edgeLen, 1 - pos1/edgeLen);
		double toBound2 = min(pos2/edgeLen, 1 - pos2/edgeLen);
		if (toBound1 < 0 && toBound2 < 0) 
		{
			pos = -1.0;
			return;
		}
		/*assert(!(toBound1 < 0 && toBound2 < 0));*/
		if (toBound1 > toBound2) 
		{
			edgeIndex = curEdge;
			pos = pos1;
		}
		else
		{
			edgeIndex = twinEdge;
			pos = pos2;
		}
	}
	else if (type == 1)
	{
		edgeIndex = curEdge;
		pos = -1.0;
	}
	else if (type == 2)
	{
		FindPerpendicularFootFarthest(curEdge, twinEdge, edgeIndex, pos);
	}
	else if (type == 3)
	{
		if (mesh->m_pEdge[curEdge].m_iVertex[0] < originalMeshVertNum || twinEdge == -1)
		{
			edgeIndex = curEdge;
			FindPosNearestToMidPoint(edgeIndex, pos);
		}
		else
		{
			edgeIndex = twinEdge;
			FindPosNearestToMidPoint(twinEdge, pos);
		}
		if (!(pos > 0 && pos < mesh->m_pEdge[edgeIndex].m_length))
		{
			cout << pos << endl;
		}
		assert(pos > 0 && pos < mesh->m_pEdge[edgeIndex].m_length);
	}
	else if (type == 4)
	{
		edgeIndex = curEdge;

		vector< pair<double, double> > bounds;
		vector<unsigned> exceptList;

		double lowerBound, upperBound;

		unsigned e[5] = {0};
		e[0] = curEdge;
		e[1] = mesh->m_pEdge[e[0]].m_iNextEdge;
		e[2] = mesh->m_pEdge[e[1]].m_iNextEdge;
		if (twinEdge != -1)
		{
			e[3] = mesh->m_pEdge[twinEdge].m_iNextEdge;
			e[4] = mesh->m_pEdge[e[3]].m_iNextEdge;
		}
		double len[5];
		for (unsigned i = 0; i < 5; ++i)
			len[i] = mesh->m_pEdge[e[i]].m_length;
		double cosOutAngle[4];
		bool flippableEdges[4];
		for (unsigned i = 1; i < 5; ++i)
		{
			unsigned e0 = mesh->m_pEdge[e[i]].m_iTwinEdge;
			if (e0 == -1)
			{
				// if it's a boundary edge, a non-obtuse angle is required for non-Delaunay condition
				cosOutAngle[i-1] = 0.0;
				flippableEdges[i-1] = false;
				continue;
			}
			unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
			unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;
			double l0 = mesh->m_pEdge[e0].m_length;
			double l1 = mesh->m_pEdge[e1].m_length;
			double l2 = mesh->m_pEdge[e2].m_length;
			cosOutAngle[i-1] = (l1*l1 + l2*l2 - l0*l0) / (2*l1*l2);

			flippableEdges[i-1] =(fatherEdge[e[i]] >= originalMeshEdgeNum && 
				e0 != -1 && fatherEdge[e0] >= originalMeshEdgeNum);
		}

		double cosAngle1 = (len[0]*len[0] + len[2]*len[2] - len[1]*len[1]) / (2*len[0]*len[2]);
		double cosAngle2 = (len[0]*len[0] + len[3]*len[3] - len[4]*len[4]) / (2*len[0]*len[3]);
		double cosAngle3 = (len[0]*len[0] + len[1]*len[1] - len[2]*len[2]) / (2*len[0]*len[1]);
		double cosAngle4 = (len[0]*len[0] + len[4]*len[4] - len[3]*len[3]) / (2*len[0]*len[4]);

		double h1 = mesh->m_pFace[mesh->m_pEdge[curEdge].m_iFace].m_dArea *2.0 / len[1];
		double h2 = mesh->m_pFace[mesh->m_pEdge[curEdge].m_iFace].m_dArea *2.0 / len[2];
		
		double h3 = 1.0, h4 = 1.0;
		if (twinEdge != -1)
		{
			h3 = mesh->m_pFace[mesh->m_pEdge[twinEdge].m_iFace].m_dArea *2.0 / len[3];
			h4 = mesh->m_pFace[mesh->m_pEdge[twinEdge].m_iFace].m_dArea *2.0 / len[4];
		}

		if (twinEdge != -1)
		{
			upperBound = (len[2]*h3 + len[3]*h2) / (h3*cosAngle1 + h2*cosAngle2);
			lowerBound = (len[1]*h4 + len[4]*h1) / (h4*cosAngle3 + h1*cosAngle4);
			lowerBound = len[0] - lowerBound;
		}
		else
		{
			upperBound = len[2] / cosAngle1;
			lowerBound = len[1] / cosAngle3;
			lowerBound = len[0] - lowerBound;
		}
		// Make the current edge Delaunay
		bounds.push_back(make_pair(lowerBound, upperBound));
		
		// Make the upper two edges Delaunay
		if (flippableEdges[0]) upperBound = len[0];
		else FindPosMakeSurroundEdgeLegal(curEdge, e[1], upperBound);
		if (flippableEdges[1]) lowerBound = 0.0;
		else FindPosMakeSurroundEdgeLegal(curEdge, e[2], lowerBound);
		bounds.push_back(make_pair(0.0, upperBound));
		bounds.push_back(make_pair(lowerBound, len[0]));

		if (twinEdge != -1)
		{
			if (flippableEdges[2]) lowerBound = len[0];
			else FindPosMakeSurroundEdgeLegal(twinEdge, e[3], lowerBound);
			if (flippableEdges[3]) upperBound = 0.0;
			else FindPosMakeSurroundEdgeLegal(twinEdge, e[4], upperBound);
			lowerBound = len[0] - lowerBound;
			upperBound = len[0] - upperBound;
		}
		else
		{
			lowerBound = 0; upperBound = len[0];
		}

		// Make the downer two edges (if there are) Delaunay
		bounds.push_back(make_pair(0.0, upperBound));
		bounds.push_back(make_pair(lowerBound, len[0]));

		for (unsigned i = 0; i < bounds.size(); ++i)
		{
			if (bounds[i].first < 0) bounds[i].first = 0;
			if (bounds[i].first > len[0]) bounds[i].first = len[0];
			if (bounds[i].second < 0) bounds[i].second = 0;
			if (bounds[i].second > len[0]) bounds[i].second = len[0];
		}

		BounderIntersection(bounds, exceptList, lowerBound, upperBound);

		if (BoundNotEmpty(lowerBound, upperBound, curEdge))
		{
			FindTheBestPosInBounder(lowerBound, upperBound, curEdge, twinEdge, len, cosAngle1, edgeIndex, pos);
			double posTmp; unsigned edgeTmp;
#ifdef DIRECTFOOT
			FindPerpendicularFoot(curEdge, posTmp); edgeTmp = curEdge;
#else
			FindPerpendicularFootFarthest(curEdge, twinEdge, edgeTmp, posTmp);
#endif
			if (edgeTmp == curEdge && pos == posTmp || edgeTmp == twinEdge && pos == len[0] - posTmp)
				++perpFootFallInFourInterval;
			else
				++FourIntervalMid;
		}
		else
		{
			unsigned smallestAngleIdx = -1;
			double smallestAngle = -1e30;
			for (unsigned i = 0; i < (twinEdge == -1 ? 2 : 4); ++i)
			{
				if (cosOutAngle[i] < smallestAngle) continue;
				if (flippableEdges[i]) continue;
				smallestAngle = cosOutAngle[i];
				smallestAngleIdx = i;
			}

			exceptList.push_back(smallestAngleIdx+1);
			for (unsigned i = 0; i < (twinEdge == -1 ? 2 : 4); ++i)
				if (flippableEdges[i]) exceptList.push_back(i+1);
			BounderIntersection(bounds, exceptList, lowerBound, upperBound);

			if (BoundNotEmpty(lowerBound, upperBound, curEdge))
			{
				FindTheBestPosInBounder(lowerBound, upperBound, curEdge, twinEdge, len, cosAngle1, edgeIndex, pos);
				double posTmp; unsigned edgeTmp;
#ifdef DIRECTFOOT
				FindPerpendicularFoot(curEdge, posTmp); edgeTmp = curEdge;
#else
				FindPerpendicularFootFarthest(curEdge, twinEdge, edgeTmp, posTmp);
#endif
				if (edgeTmp == curEdge && pos == posTmp || edgeTmp == twinEdge && pos == len[0] - posTmp)
					++perpFootFallInThreeInterval;
				else
					++ThreeIntervalMid;
			}
			else
			{
				FindPerpendicularFootFarthest(curEdge, twinEdge, edgeIndex, pos);
				++justPerpFoot;
			}
		}
	}
	else if (type == 5)
	{
		vector< pair<double, double> > bounds;
		double lowerBound, upperBound;

		unsigned e[5] = {0};
		e[0] = curEdge;
		e[1] = mesh->m_pEdge[e[0]].m_iNextEdge;
		e[2] = mesh->m_pEdge[e[1]].m_iNextEdge;
		if (twinEdge != -1)
		{
			e[3] = mesh->m_pEdge[twinEdge].m_iNextEdge;
			e[4] = mesh->m_pEdge[e[3]].m_iNextEdge;
		}

		unsigned nonDelaunayCnt = 0;
		for (unsigned i = 1; i < 5; ++i)
			if (!isEdgeLegal(mesh, e[i])) ++nonDelaunayCnt;
		switch(nonDelaunayCnt)
		{
		case 0: ++surroundExistZeroDelaunayEdge; break;
		case 1: ++surroundExistOneDelaunayEdge; break;
		case 2: ++surroundExistTwoDelaunayEdge; break;
		case 3: ++surroundExistThreeDelaunayEdge; break;
		case 4: ++surroundExistFourDelaunayEdge; break;
		}

		double len[5];
		for (unsigned i = 0; i < 5; ++i)
			len[i] = mesh->m_pEdge[e[i]].m_length;
		
		bool flippableEdges[4];
		for (unsigned i = 1; i < 5; ++i)
		{
			unsigned e0 = mesh->m_pEdge[e[i]].m_iTwinEdge;
			if (e0 == -1)
			{
				// if it's a boundary edge, a non-obtuse angle is required for non-Delaunay condition
				flippableEdges[i-1] = false;
				continue;
			}
			unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
			unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;
			double l0 = mesh->m_pEdge[e0].m_length;
			double l1 = mesh->m_pEdge[e1].m_length;
			double l2 = mesh->m_pEdge[e2].m_length;

			flippableEdges[i-1] =(fatherEdge[e[i]] >= originalMeshEdgeNum && 
				e0 != -1 && fatherEdge[e0] >= originalMeshEdgeNum);
		}

		double cosAngle1 = (len[0]*len[0] + len[2]*len[2] - len[1]*len[1]) / (2*len[0]*len[2]);
		double cosAngle2 = (len[0]*len[0] + len[3]*len[3] - len[4]*len[4]) / (2*len[0]*len[3]);
		double cosAngle3 = (len[0]*len[0] + len[1]*len[1] - len[2]*len[2]) / (2*len[0]*len[1]);
		double cosAngle4 = (len[0]*len[0] + len[4]*len[4] - len[3]*len[3]) / (2*len[0]*len[4]);

		double h1 = mesh->m_pFace[mesh->m_pEdge[curEdge].m_iFace].m_dArea *2.0 / len[1];
		double h2 = mesh->m_pFace[mesh->m_pEdge[curEdge].m_iFace].m_dArea *2.0 / len[2];

		double h3 = 1.0, h4 = 1.0;
		if (twinEdge != -1)
		{
			h3 = mesh->m_pFace[mesh->m_pEdge[twinEdge].m_iFace].m_dArea *2.0 / len[3];
			h4 = mesh->m_pFace[mesh->m_pEdge[twinEdge].m_iFace].m_dArea *2.0 / len[4];
		}

		if (twinEdge != -1)
		{
			upperBound = (len[2]*h3 + len[3]*h2) / (h3*cosAngle1 + h2*cosAngle2);
			lowerBound = (len[1]*h4 + len[4]*h1) / (h4*cosAngle3 + h1*cosAngle4);
			lowerBound = len[0] - lowerBound;
		}
		else
		{
			upperBound = len[2] / cosAngle1;
			lowerBound = len[1] / cosAngle3;
			lowerBound = len[0] - lowerBound;
		}
		// Make the current edge Delaunay
		bounds.push_back(make_pair(lowerBound, upperBound));

		// Make the upper two edges Delaunay
		/*if (flippableEdges[0]) upperBound = len[0];*/
		/*else */FindPosMakeSurroundEdgeLegal(curEdge, e[1], upperBound);
		/*if (flippableEdges[1]) lowerBound = 0.0;*/
		/*else */FindPosMakeSurroundEdgeLegal(curEdge, e[2], lowerBound);
		bounds.push_back(make_pair(0.0, upperBound));
		bounds.push_back(make_pair(lowerBound, len[0]));

		if (twinEdge != -1)
		{
			/*if (flippableEdges[2]) lowerBound = len[0];*/
			/*else */FindPosMakeSurroundEdgeLegal(twinEdge, e[3], lowerBound);
			/*if (flippableEdges[3]) upperBound = 0.0;*/
			/*else */FindPosMakeSurroundEdgeLegal(twinEdge, e[4], upperBound);
			lowerBound = len[0] - lowerBound;
			upperBound = len[0] - upperBound;
		}
		else
		{
			lowerBound = 0; upperBound = len[0];
		}

		// Make the downer two edges (if there are) Delaunay
		bounds.push_back(make_pair(0.0, upperBound));
		bounds.push_back(make_pair(lowerBound, len[0]));

		for (unsigned i = 0; i < bounds.size(); ++i)
		{
			if (bounds[i].first < 0) bounds[i].first = 0;
			if (bounds[i].first > len[0]) bounds[i].first = len[0];
			if (bounds[i].second < 0) bounds[i].second = 0;
			if (bounds[i].second > len[0]) bounds[i].second = len[0];
		}

		set<double> endPoints;
		for (unsigned i = 1; i < bounds.size(); ++i)
		{
			endPoints.insert(bounds[i].first); endPoints.insert(bounds[i].second);
		}

		lowerBound = bounds[0].first;
		upperBound = bounds[0].second;
		unsigned maxDelaunayNum = 0;
		double optLB = lowerBound, optUB = upperBound;
		for (set<double>::iterator iter = endPoints.begin(); iter != endPoints.end(); ++iter)
		{
			if (*iter < bounds[0].first) continue;
			if (*iter > bounds[0].second) break;

			double mid = (lowerBound + *iter) / 2.0;

			unsigned curDelaunayNum = 0;
			for (unsigned i = 1; i < bounds.size(); ++i)
				if (InInterval(mid, bounds[i])) ++curDelaunayNum;

			if (curDelaunayNum > maxDelaunayNum)
			{
				maxDelaunayNum = curDelaunayNum;
				optLB = lowerBound; optUB = *iter;
			}
			lowerBound = *iter;
		}

		double mid = (lowerBound + upperBound) / 2.0;
		unsigned curDelaunayNum = 0;
		for (unsigned i = 1; i < bounds.size(); ++i)
			if (InInterval(mid, bounds[i])) ++curDelaunayNum;

		if (curDelaunayNum > maxDelaunayNum)
		{
			maxDelaunayNum = curDelaunayNum;
			optLB = lowerBound; optUB = upperBound;
		}
		unsigned perpenEdge = -1; double perpenPos = -1;
		FindPerpendicularFootFarthest(curEdge, twinEdge, perpenEdge, perpenPos);
		if (perpenEdge == twinEdge) 
		{
			perpenEdge = curEdge;
			perpenPos = len[0] - perpenPos;
		}
		if (InInterval(perpenPos, make_pair(optLB, optUB)))
		{
			edgeIndex = perpenEdge;
			pos = perpenPos;
		}
		else
		{
			edgeIndex = curEdge;
			pos = (optLB + optUB) / 2.0;
		}
	}
}

void DelaunayMesh::FindPerpendicularFootFarthest(unsigned curEdge, unsigned twinEdge, unsigned &edgeIndex, double &pos)
{
	double pos1 = 1e30, pos2 = 1e30;
	double edgeLen = mesh->m_pEdge[curEdge].m_length;
	FindPosNearestToPerpendicular(curEdge, pos1);
	FindPosNearestToPerpendicular(twinEdge, pos2);

	double toBound1 = min(pos1/edgeLen, 1 - pos1/edgeLen);
	double toBound2 = min(pos2/edgeLen, 1 - pos2/edgeLen);
	if (toBound1 < 0 && toBound2 < 0) 
	{
		pos = -1.0;
		return;
	}
	/*assert(!(toBound1 < 0 && toBound2 < 0));*/
	if (toBound1 > toBound2) 
	{
		edgeIndex = curEdge;
		pos = pos1;
	}
	else
	{
		edgeIndex = twinEdge;
		pos = pos2;
	}
}

void DelaunayMesh::FindTheBestPosInBounder(double lowerBound, double upperBound, 
	unsigned curEdge, unsigned twinEdge, double len[5], double cosAngle, 
	unsigned &edgeIndex, double &pos)
{
#ifdef DIRECTFOOT
 	double perpendicularFoot;
	FindPerpendicularFoot(curEdge, perpendicularFoot);
	
	if (perpendicularFoot >= lowerBound && perpendicularFoot <= upperBound)
	{
		/*FindPerpendicularFootFarthest(curEdge, twinEdge, edgeIndex, pos);*/
		edgeIndex = curEdge; pos = perpendicularFoot;
	}
	else 
	{
		pos = (lowerBound+upperBound) / 2.0;
	}
#else
	// the optimal interval is not empty
	FindPerpendicularFootFarthest(curEdge, twinEdge, edgeIndex, pos);
	if (edgeIndex == twinEdge)
	{
		edgeIndex = curEdge;
		pos = len[0] - pos;
	}

	if (pos < lowerBound || pos > upperBound)
	{
		pos = (lowerBound+upperBound) / 2.0;
	}
#endif
}

void DelaunayMesh::FindPerpendicularFoot(unsigned edgeIndex, double &pos)
{
	if (edgeIndex == -1) {pos = 1e30; return;}
	unsigned e0 = edgeIndex;
	unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
	unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;
	double len0 = mesh->m_pEdge[e0].m_length;
	double len1 = mesh->m_pEdge[e1].m_length;
	double len2 = mesh->m_pEdge[e2].m_length;

	pos = (len0*len0+len2*len2 - len1*len1) / (2*len0);
}

void DelaunayMesh::FindPosNearestToMidPoint(unsigned edgeIndex, double &pos)
{
	pos = log10(mesh->m_pEdge[edgeIndex].m_length / 2.0 / delta)/log10(2.0);
	int k1 = ceil(pos), k2 = floor(pos);
// 	if (k1 == INT_MAX || k2 == INT_MAX)
// 	{
// 		cout << "k1: " << k1 << " k2: " << k2 << endl;
// 		system("pause");
// 	}
	double pos1 = pow(2.0, (double)k1)*delta, pos2 = pow(2.0, (double)k2)*delta;
	if ((mesh->m_pEdge[edgeIndex].m_length - pos1) / mesh->m_pEdge[edgeIndex].m_length < LHDOUBLE_EPS)
		pos = pos2;
	else 
	{
// 		Vector3D p0 = mesh->m_pVertex[mesh->m_pEdge[edgeIndex].m_iVertex[0]].m_vPosition;
// 		Vector3D d0 = mesh->m_pVertex[mesh->m_pEdge[edgeIndex].m_iVertex[1]].m_vPosition - 
// 			mesh->m_pVertex[mesh->m_pEdge[edgeIndex].m_iVertex[0]].m_vPosition;
// 		d0.normalize();

		if (abs(mesh->m_pEdge[edgeIndex].m_length/2.0 - pos1) < abs(mesh->m_pEdge[edgeIndex].m_length/2.0 - pos2))
			pos = pos1;
		else pos = pos2;
// 		Vector3D p1 = p0 + d0 * pos1, p2 = p0 + d0 * pos2;
// 		Vector3D q = mesh->m_pVertex[mesh->m_pEdge[mesh->m_pEdge[edgeIndex].m_iNextEdge].m_iVertex[1]].m_vPosition;
// 
// 		if ((p1-q).length() < (p2-q).length()) pos = pos1;
// 		else pos = pos2;
	}
}

void DelaunayMesh::FindPosNearestToPerpendicular(unsigned edgeIndex, double &pos)
{
	FindPerpendicularFoot(edgeIndex, pos);
	if (edgeIndex == -1 || pos < 0 || pos > mesh->m_pEdge[edgeIndex].m_length) return;
	/*assert(mesh->m_pEdge[edgeIndex].m_length > rhoV);*/
	double leftLen = (mesh->m_pVertex[mesh->m_pEdge[edgeIndex].m_iVertex[0]].m_vPosition - 
		mesh->m_pVertex[originalEdgeEndVerts[fatherEdge[edgeIndex]].first].m_vPosition).length();
	double kL = (leftLen+pos - rhoV) / rhoE;
	kL = kL < 0 ? 0 : kL;
	
	if (kL < INT_MAX) kL = (int)kL;
	
	double kU = kL + 1;
	double pos1 = kL*rhoE + rhoV - leftLen, pos2 = kU*rhoE + rhoV - leftLen;
	/*assert(pos1 > 0 && pos1 < mesh->m_pEdge[edgeIndex].m_length);*/
	if (pos1 < rhoV) { pos1 = rhoV; pos2 = rhoV + rhoE; }

	if ((mesh->m_pEdge[edgeIndex].m_length-pos2) / mesh->m_pEdge[edgeIndex].m_length < LHDOUBLE_EPS)
		pos = pos1;
	else
	{
		Vector3D p0 = mesh->m_pVertex[mesh->m_pEdge[edgeIndex].m_iVertex[0]].m_vPosition;
		Vector3D d0 = mesh->m_pVertex[mesh->m_pEdge[edgeIndex].m_iVertex[1]].m_vPosition - 
			mesh->m_pVertex[mesh->m_pEdge[edgeIndex].m_iVertex[0]].m_vPosition;
		d0.normalize();
		Vector3D p1 = p0 + d0 * pos1, p2 = p0 + d0 * pos2;
		Vector3D q = mesh->m_pVertex[mesh->m_pEdge[mesh->m_pEdge[edgeIndex].m_iNextEdge].m_iVertex[1]].m_vPosition;

		if ((p1-q).length() < (p2-q).length()) pos = pos1;
		else pos = pos2;
	}
}

void DelaunayMesh::FindPosMakeSurroundEdgeLegal(unsigned baseEdge, unsigned surroundEdge, double &pos)
{
	unsigned surroundTwinEdge = mesh->m_pEdge[surroundEdge].m_iTwinEdge;
	unsigned e[3]; double len[3];
	if (surroundTwinEdge == -1)
	{
		FindPerpendicularFoot(baseEdge, pos);
		return;
	}
	e[0] = surroundTwinEdge;
	e[1] = mesh->m_pEdge[e[0]].m_iNextEdge;
	e[2] = mesh->m_pEdge[e[1]].m_iNextEdge;

	len[0] = mesh->m_pEdge[e[0]].m_length;
	len[1] = mesh->m_pEdge[e[1]].m_length;
	len[2] = mesh->m_pEdge[e[2]].m_length;

	double cosAngle1 = - (len[1]*len[1] + len[2]*len[2] - len[0]*len[0]) / (2*len[1]*len[2]);	
	double angle1 = acos(cosAngle1);

	e[0] = baseEdge;
	e[1] = mesh->m_pEdge[e[0]].m_iNextEdge;
	e[2] = mesh->m_pEdge[e[1]].m_iNextEdge;

	len[0] = mesh->m_pEdge[e[0]].m_length;
	len[1] = mesh->m_pEdge[e[1]].m_length;
	len[2] = mesh->m_pEdge[e[2]].m_length;

	double cosAngle2;
	assert(surroundEdge == e[1] || surroundEdge == e[2]);
	if (surroundEdge == e[1]) cosAngle2 = (len[0]*len[0]+len[1]*len[1]-len[2]*len[2]) / (2*len[0]*len[1]);
	else cosAngle2 = (len[0]*len[0]+len[2]*len[2]-len[1]*len[1]) / (2*len[0]*len[2]);
	double angle2 = acos(cosAngle2);

	double angle4 = acos((len[1]*len[1]+len[2]*len[2]-len[0]*len[0]) / (2*len[1]*len[2]));

	double angle3 = PI - angle1 - angle2;
	// if angle3 < 0 or angle3 > angle4, the pos will exceeds[0, len[0]]
	pos = sin(angle3) * mesh->m_pEdge[surroundTwinEdge].m_length / sin(angle1);
	if (surroundEdge == e[1])
		pos = len[0] - pos;
}

void DelaunayMesh::BounderIntersection(vector< pair<double, double> > &bounds, vector<unsigned> &exceptList, double &lowerBound, double &upperBound)
{
	lowerBound = 0.0; upperBound = 1e30;
	for (unsigned i = 0; i < bounds.size(); ++i)
	{
		lowerBound = max(lowerBound, bounds[i].first);
		upperBound = min(upperBound, bounds[i].second);
	}
	
	for (unsigned i = 0; i < bounds.size(); ++i)
	{
		bool inExceptList = false;
		for (unsigned j = 0; j < exceptList.size(); ++j)
		{
			if (exceptList[j] == i)
			{
				inExceptList = true;
				break;
			}
		}
		if (inExceptList) continue;

		lowerBound = max(lowerBound, bounds[i].first);
		upperBound = min(upperBound, bounds[i].second);
	}
}

bool DelaunayMesh::BoundNotEmpty(double lowerBound, double upperBound, unsigned curEdge)
{
	if (lowerBound < 0) lowerBound = 0; if (lowerBound > mesh->m_pEdge[curEdge].m_length) lowerBound = mesh->m_pEdge[curEdge].m_length;
	if (upperBound < 0) upperBound = 0; if (upperBound > mesh->m_pEdge[curEdge].m_length) upperBound = mesh->m_pEdge[curEdge].m_length;
	return (lowerBound - upperBound) / mesh->m_pEdge[curEdge].m_length < -1e-2;
}

void DelaunayMesh::MergeTwoNeighPointOnPhysicalEdge(unsigned iEdge, double leftBound, double rightBound)
{
	unsigned iv0 = mesh->m_pEdge[iEdge].m_iVertex[0];
	unsigned iv1 = mesh->m_pEdge[iEdge].m_iVertex[1];
	
	Vector3D originFinal = mesh->m_pVertex[originalEdgeEndVerts[fatherEdge[iEdge]].first].m_vPosition;
	Vector3D endFinal = mesh->m_pVertex[originalEdgeEndVerts[fatherEdge[iEdge]].second].m_vPosition;
	Vector3D XaxisFinal = endFinal - originFinal; XaxisFinal.normalize();

	list< pair<double, double> > feasibleIntervals;
	feasibleIntervals.push_back(make_pair(leftBound, rightBound));

	unsigned iTwinEdge = mesh->m_pEdge[iEdge].m_iTwinEdge;
	unsigned oE0 = originalMesh->m_pEdge[fatherEdge[iEdge]].m_iNextEdge;
	unsigned oE1 = originalMesh->m_pEdge[oE0].m_iNextEdge;
	unsigned oE2 = -1, oE3 = -1;
	if (iTwinEdge != -1)
	{
		oE2 = originalMesh->m_pEdge[fatherEdge[iTwinEdge]].m_iNextEdge;
		oE3 = originalMesh->m_pEdge[oE2].m_iNextEdge;
	}
	double lFinal = originalMesh->m_pEdge[fatherEdge[iEdge]].m_length;
	double l0 = originalMesh->m_pEdge[oE0].m_length;
	double l1 = originalMesh->m_pEdge[oE1].m_length;
	double l2 = 0, l3 = 0;
	if (iTwinEdge != -1)
	{
		l2 = originalMesh->m_pEdge[oE2].m_length;
		l3 = originalMesh->m_pEdge[oE3].m_length;
	}

	unsigned commonFace0 = mesh->m_pEdge[iEdge].m_iFace;
	unsigned commonFace1 = -1;
	if (iTwinEdge != -1) commonFace1 = mesh->m_pEdge[iTwinEdge].m_iFace;

	unsigned vSpecial0 = -1, vSpecial1 = -1;

	Vector3D original0 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE0].m_iVertex[1]].m_vPosition;
	Vector3D Xaxis0 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE0].m_iVertex[0]].m_vPosition - 
		originalMesh->m_pVertex[originalMesh->m_pEdge[oE0].m_iVertex[1]].m_vPosition);
	Xaxis0.normalize();

	Vector3D original1 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE1].m_iVertex[1]].m_vPosition;
	Vector3D Xaxis1 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE1].m_iVertex[0]].m_vPosition - 
		originalMesh->m_pVertex[originalMesh->m_pEdge[oE1].m_iVertex[1]].m_vPosition);
	Xaxis1.normalize();

	Vector3D original2, Xaxis2, original3, Xaxis3;
	if (iTwinEdge != -1)
	{
		original2 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE2].m_iVertex[0]].m_vPosition;
		Xaxis2 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE2].m_iVertex[1]].m_vPosition - 
			originalMesh->m_pVertex[originalMesh->m_pEdge[oE2].m_iVertex[0]].m_vPosition);
		Xaxis2.normalize();

		original3 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE3].m_iVertex[0]].m_vPosition;
		Xaxis3 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE3].m_iVertex[1]].m_vPosition - 
			originalMesh->m_pVertex[originalMesh->m_pEdge[oE3].m_iVertex[0]].m_vPosition);
		Xaxis3.normalize();
	}

	Vector2D Xaxis1_2D, Xaxis2_2D;
	Xaxis1_2D.x = (l1*l1 + lFinal*lFinal - l0*l0) / (2*lFinal);
	Xaxis1_2D.y = sqrt(l1*l1 - Xaxis1_2D.x*Xaxis1_2D.x);
	Xaxis2_2D.x = (l2*l2 + lFinal*lFinal - l3*l3) / (2*lFinal);
	Xaxis2_2D.y = -sqrt(l2*l2 - Xaxis2_2D.x*Xaxis2_2D.x);
	Xaxis1_2D.normalize(); Xaxis2_2D.normalize();

	vector< unsigned > edgesOnPhysicalEdges;
	for (unsigned i = 0; i < mesh->m_pVertex[iv0].m_nValence; ++i)
	{
		unsigned adjEdge = mesh->m_pVertex[iv0].m_piEdge[i];
		unsigned opEdge = mesh->m_pEdge[adjEdge].m_iNextEdge;
		unsigned prevEdge = mesh->m_pEdge[opEdge].m_iNextEdge;
		unsigned prevTwinEdge = mesh->m_pEdge[prevEdge].m_iTwinEdge;

		if (mesh->m_pEdge[adjEdge].m_iFace == commonFace0 || 
			mesh->m_pEdge[adjEdge].m_iFace == commonFace1)
			continue;

		edgesOnPhysicalEdges.push_back(adjEdge);
		edgesOnPhysicalEdges.push_back(opEdge);
		if (prevTwinEdge == -1) edgesOnPhysicalEdges.push_back(prevEdge);
		else if (mesh->m_pEdge[prevTwinEdge].m_iFace == commonFace1)
			vSpecial0 = mesh->m_pEdge[adjEdge].m_iVertex[1];
	}
	for (unsigned i = 0; i < mesh->m_pVertex[iv1].m_nValence; ++i)
	{
		unsigned adjEdge = mesh->m_pVertex[iv1].m_piEdge[i];
		unsigned opEdge = mesh->m_pEdge[adjEdge].m_iNextEdge;
		unsigned prevEdge = mesh->m_pEdge[opEdge].m_iNextEdge;
		unsigned prevTwinEdge = mesh->m_pEdge[prevEdge].m_iTwinEdge;

		if (mesh->m_pEdge[adjEdge].m_iFace == commonFace0 || 
			mesh->m_pEdge[adjEdge].m_iFace == commonFace1)
			continue;

		edgesOnPhysicalEdges.push_back(adjEdge);
		edgesOnPhysicalEdges.push_back(opEdge);
		if (prevTwinEdge == -1) edgesOnPhysicalEdges.push_back(prevEdge);
		else if (mesh->m_pEdge[prevTwinEdge].m_iFace == commonFace0)
			vSpecial1 = mesh->m_pEdge[adjEdge].m_iVertex[1];
	}

	// find the interval that makes the edgesOnPhysicalEdges all Delaunay
	for (unsigned i = 0; i < edgesOnPhysicalEdges.size(); ++i)
	{
		unsigned curEdge = edgesOnPhysicalEdges[i];
		unsigned twinEdge = mesh->m_pEdge[curEdge].m_iTwinEdge;

		unsigned e0 = mesh->m_pEdge[curEdge].m_iNextEdge;
		unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
		unsigned e2 = -1, e3 = -1;
		if (twinEdge != -1)
		{
			e2 = mesh->m_pEdge[twinEdge].m_iNextEdge;
			e3 = mesh->m_pEdge[e2].m_iNextEdge;
		}
		unsigned v[4] = {-1};
		v[0] = mesh->m_pEdge[curEdge].m_iVertex[0];
		v[1] = mesh->m_pEdge[curEdge].m_iVertex[1];
		v[2] = mesh->m_pEdge[e0].m_iVertex[1];
		if (e2 != -1) v[3] = mesh->m_pEdge[e2].m_iVertex[1];

		if (twinEdge != -1 && mesh->m_pEdge[twinEdge].m_iFace == commonFace0)
			v[3] = vSpecial1;
		else if (twinEdge == -1 && commonFace1 == -1 || twinEdge != -1 && mesh->m_pEdge[twinEdge].m_iFace == commonFace1)
			v[3] = vSpecial0;

		vector<Vector2D> pointsOnCircle;
		Vector2D center; double R;

		if (fatherEdge[curEdge] == oE0)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				if (v[j] == -1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original0;
				double x = vec0 * Xaxis0, y = (vec0 ^ Xaxis0).normalize();
				double cosAngle = Xaxis0*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// clock wise
				double x0 = cosAngle*x + sinAngle*y;
				double y0 = -sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0) + Xaxis1_2D*l1);
			}
		}
		else if (fatherEdge[curEdge] == oE1)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				if (v[j] == -1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original1;
				double x = vec0 * Xaxis1, y = (vec0 ^ Xaxis1).normalize();
				double cosAngle = Xaxis1*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// counter clock wise
				double x0 = cosAngle*x - sinAngle*y;
				double y0 = sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0));
			}
		}
		else if (fatherEdge[curEdge] == oE2)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				if (v[j] == -1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original2;
				double x = vec0 * Xaxis2, y = -(vec0 ^ Xaxis2).normalize();
				double cosAngle = Xaxis2*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// clock wise
				double x0 = cosAngle*x + sinAngle*y;
				double y0 = -sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0));
			}
		}
		else if (fatherEdge[curEdge] == oE3)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				if (v[j] == -1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original3;
				double x = vec0 * Xaxis3, y = -(vec0 ^ Xaxis3).normalize();
				double cosAngle = Xaxis3*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// counter clock wise
				double x0 = cosAngle*x - sinAngle*y;
				double y0 = sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0) + Xaxis2_2D*l2);
			}
		}
		else
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				if (v[j] == -1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - originFinal;
				double x = vec0 * XaxisFinal, y = (vec0 ^ XaxisFinal).normalize();
				if ((vec0 ^ XaxisFinal) * (Xaxis1 ^ XaxisFinal) < 0) y = -y;
				pointsOnCircle.push_back(Vector2D(x, y));
			}
		}

		CalculateCircumCircle(pointsOnCircle, center, R);
		if (fabs(center.y) < R)
		{
			list< pair<double, double> > curIntervals;
			double width = sqrt(R*R - center.y*center.y);
			if (mesh->m_pEdge[curEdge].m_iVertex[0] == iv0 || 
				mesh->m_pEdge[curEdge].m_iVertex[1] == iv0 || 
				mesh->m_pEdge[curEdge].m_iVertex[0] == iv1 || 
				mesh->m_pEdge[curEdge].m_iVertex[1] == iv1)
			{
				curIntervals.push_back(make_pair(center.x - width, center.x + width));
			}
			else
			{
				curIntervals.push_back(make_pair(leftBound, center.x - width));
				curIntervals.push_back(make_pair(center.x + width, rightBound));
			}
			IntersectWithFeasibleIntervals(feasibleIntervals, curIntervals);
		}
	}

	// if the result interval is not empty
	if (!feasibleIntervals.empty())
	{
		// readjust the postion of iv0
		++mergeVertCnt;
		cout << "Merge vert " << iv0 << " " << iv1 << endl;
		double pos = (feasibleIntervals.begin()->first + feasibleIntervals.begin()->second) / 2.0;
		Vector3D newPos = originFinal + XaxisFinal * pos;
		mesh->collapse(iEdge, newPos);
		/*mesh->moveVertTo(iv0, newPos);*/
	}
}

void DelaunayMesh::CalculateCircumCircle(vector<Vector2D> &points, Vector2D &center, double &R)
{
	vector<Vector3D> points3D;
	points3D.resize(points.size());
	for (unsigned i = 0; i < points.size(); ++i)
		points3D[i] = Vector3D(points[i].x, points[i].y, 0);

	double S = ((points3D[1] - points3D[0]) ^ (points3D[2] - points3D[0])).length();
	double alpha = (points3D[1] - points3D[2]).length2() * (points3D[0] - points3D[1]) * (points3D[0] - points3D[2]);
	double beta = (points3D[0] - points3D[2]).length2() * (points3D[1] - points3D[0]) * (points3D[1] - points3D[2]);
	double gamma = (points3D[0] - points3D[1]).length2() * (points3D[2] - points3D[0]) * (points3D[2] - points3D[1]);

	alpha /= (2 * S * S); beta /= (2 * S * S); gamma /= (2 * S * S); 

	center = alpha * points[0] + beta * points[1] + gamma * points[2];
	R = (points[0] - points[1]).length() * (points[1] - points[2]).length() * (points[2] - points[0]).length() / (2*S);
}

void DelaunayMesh::IntersectWithFeasibleIntervals(list< pair<double, double> > &feasibleIntervals, list< pair<double, double> > &intervals)
{
	unsigned feasibleSize = feasibleIntervals.size();
	unsigned cnt = 0;

	for (list< pair<double, double> >::iterator iter = feasibleIntervals.begin(); cnt < feasibleSize; )
	{
		for (list< pair<double, double> >::iterator iter2 = intervals.begin(); iter2 != intervals.end(); ++iter2)
		{
			double L = max(iter->first, iter2->first);
			double R = min(iter->second, iter2->second);
			if (L < R) feasibleIntervals.push_back(make_pair(L, R));
		}
		list< pair<double, double> >::iterator iterNext = iter; ++iterNext;
		feasibleIntervals.erase(iter);
		iter = iterNext;
		++cnt;
	}
}

void DelaunayMesh::MergeThreeNeighPointOnPhysicalEdge(unsigned iEdge0, unsigned iEdge1, double leftBound, double rightBound)
{
	assert(mesh->m_pEdge[iEdge0].m_iVertex[1] == mesh->m_pEdge[iEdge1].m_iVertex[0]);
	assert(fatherEdge[iEdge0] == fatherEdge[iEdge1]);

	unsigned iv0 = mesh->m_pEdge[iEdge0].m_iVertex[0];
	unsigned iv1 = mesh->m_pEdge[iEdge1].m_iVertex[1];

	Vector3D originFinal = mesh->m_pVertex[originalEdgeEndVerts[fatherEdge[iEdge0]].first].m_vPosition;
	Vector3D endFinal = mesh->m_pVertex[originalEdgeEndVerts[fatherEdge[iEdge0]].second].m_vPosition;
	Vector3D XaxisFinal = endFinal - originFinal; XaxisFinal.normalize();

	list< pair<double, double> > feasibleIntervals;
	feasibleIntervals.push_back(make_pair(leftBound, rightBound));

	unsigned iTwinEdge0 = mesh->m_pEdge[iEdge0].m_iTwinEdge;
	unsigned iTwinEdge1 = mesh->m_pEdge[iEdge1].m_iTwinEdge;
	unsigned oE0 = originalMesh->m_pEdge[fatherEdge[iEdge0]].m_iNextEdge;
	unsigned oE1 = originalMesh->m_pEdge[oE0].m_iNextEdge;
	unsigned oE2 = -1, oE3 = -1;
	if (iTwinEdge0 != -1)
	{
		oE2 = originalMesh->m_pEdge[fatherEdge[iTwinEdge0]].m_iNextEdge;
		oE3 = originalMesh->m_pEdge[oE2].m_iNextEdge;
	}
	double lFinal = originalMesh->m_pEdge[fatherEdge[iEdge0]].m_length;
	double l0 = originalMesh->m_pEdge[oE0].m_length;
	double l1 = originalMesh->m_pEdge[oE1].m_length;
	double l2 = 0, l3 = 0;
	if (iTwinEdge0 != -1)
	{
		l2 = originalMesh->m_pEdge[oE2].m_length;
		l3 = originalMesh->m_pEdge[oE3].m_length;
	}

	unsigned commonFace00 = mesh->m_pEdge[iEdge0].m_iFace;
	unsigned commonFace10 = mesh->m_pEdge[iEdge1].m_iFace;
	unsigned commonFace01 = -1, commonFace11 = -1;
	if (iTwinEdge0 != -1) commonFace01 = mesh->m_pEdge[iTwinEdge0].m_iFace;
	if (iTwinEdge1 != -1) commonFace11 = mesh->m_pEdge[iTwinEdge1].m_iFace;

	unsigned vSpecial0 = -1, vSpecial1 = -1;

	Vector3D original0 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE0].m_iVertex[1]].m_vPosition;
	Vector3D Xaxis0 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE0].m_iVertex[0]].m_vPosition - 
		originalMesh->m_pVertex[originalMesh->m_pEdge[oE0].m_iVertex[1]].m_vPosition);
	Xaxis0.normalize();

	Vector3D original1 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE1].m_iVertex[1]].m_vPosition;
	Vector3D Xaxis1 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE1].m_iVertex[0]].m_vPosition - 
		originalMesh->m_pVertex[originalMesh->m_pEdge[oE1].m_iVertex[1]].m_vPosition);
	Xaxis1.normalize();

	Vector3D original2, Xaxis2, original3, Xaxis3;
	if (iTwinEdge0 != -1)
	{
		original2 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE2].m_iVertex[0]].m_vPosition;
		Xaxis2 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE2].m_iVertex[1]].m_vPosition - 
			originalMesh->m_pVertex[originalMesh->m_pEdge[oE2].m_iVertex[0]].m_vPosition);
		Xaxis2.normalize();

		original3 = originalMesh->m_pVertex[originalMesh->m_pEdge[oE3].m_iVertex[0]].m_vPosition;
		Xaxis3 = (originalMesh->m_pVertex[originalMesh->m_pEdge[oE3].m_iVertex[1]].m_vPosition - 
			originalMesh->m_pVertex[originalMesh->m_pEdge[oE3].m_iVertex[0]].m_vPosition);
		Xaxis3.normalize();
	}

	Vector2D Xaxis1_2D, Xaxis2_2D;
	Xaxis1_2D.x = (l1*l1 + lFinal*lFinal - l0*l0) / (2*lFinal);
	Xaxis1_2D.y = sqrt(l1*l1 - Xaxis1_2D.x*Xaxis1_2D.x);
	Xaxis2_2D.x = (l2*l2 + lFinal*lFinal - l3*l3) / (2*lFinal);
	Xaxis2_2D.y = -sqrt(l2*l2 - Xaxis2_2D.x*Xaxis2_2D.x);
	Xaxis1_2D.normalize(); Xaxis2_2D.normalize();

	vector< unsigned > edgesOnPhysicalEdges;
	for (unsigned i = 0; i < mesh->m_pVertex[iv0].m_nValence; ++i)
	{
		unsigned adjEdge = mesh->m_pVertex[iv0].m_piEdge[i];
		unsigned opEdge = mesh->m_pEdge[adjEdge].m_iNextEdge;
		unsigned prevEdge = mesh->m_pEdge[opEdge].m_iNextEdge;
		unsigned prevTwinEdge = mesh->m_pEdge[prevEdge].m_iTwinEdge;

		if (mesh->m_pEdge[adjEdge].m_iFace == commonFace00 || 
			mesh->m_pEdge[adjEdge].m_iFace == commonFace01)
			continue;

		edgesOnPhysicalEdges.push_back(adjEdge);
		edgesOnPhysicalEdges.push_back(opEdge);
		if (prevTwinEdge == -1) edgesOnPhysicalEdges.push_back(prevEdge);
		else if (mesh->m_pEdge[prevTwinEdge].m_iFace == commonFace01)
			vSpecial0 = mesh->m_pEdge[adjEdge].m_iVertex[1];
	}
	for (unsigned i = 0; i < mesh->m_pVertex[iv1].m_nValence; ++i)
	{
		unsigned adjEdge = mesh->m_pVertex[iv1].m_piEdge[i];
		unsigned opEdge = mesh->m_pEdge[adjEdge].m_iNextEdge;
		unsigned prevEdge = mesh->m_pEdge[opEdge].m_iNextEdge;
		unsigned prevTwinEdge = mesh->m_pEdge[prevEdge].m_iTwinEdge;

		if (mesh->m_pEdge[adjEdge].m_iFace == commonFace10 || 
			mesh->m_pEdge[adjEdge].m_iFace == commonFace11)
			continue;

		edgesOnPhysicalEdges.push_back(adjEdge);
		edgesOnPhysicalEdges.push_back(opEdge);
		if (prevTwinEdge == -1) edgesOnPhysicalEdges.push_back(prevEdge);
		else if (mesh->m_pEdge[prevTwinEdge].m_iFace == commonFace10)
			vSpecial1 = mesh->m_pEdge[adjEdge].m_iVertex[1];
	}

	// find the interval that makes the edgesOnPhysicalEdges all Delaunay
	for (unsigned i = 0; i < edgesOnPhysicalEdges.size(); ++i)
	{
		unsigned curEdge = edgesOnPhysicalEdges[i];
		unsigned twinEdge = mesh->m_pEdge[curEdge].m_iTwinEdge;

		unsigned e0 = mesh->m_pEdge[curEdge].m_iNextEdge;
		unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
		unsigned e2 = -1, e3 = -1;
		if (twinEdge != -1)
		{
			e2 = mesh->m_pEdge[twinEdge].m_iNextEdge;
			e3 = mesh->m_pEdge[e2].m_iNextEdge;
		}
		unsigned v[4] = {-1};
		v[0] = mesh->m_pEdge[curEdge].m_iVertex[0];
		v[1] = mesh->m_pEdge[curEdge].m_iVertex[1];
		v[2] = mesh->m_pEdge[e0].m_iVertex[1];
		if (e2 != -1) v[3] = mesh->m_pEdge[e2].m_iVertex[1];

		if (mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iTwinEdge].m_iFace == commonFace00)
			v[3] = vSpecial1;
		else if (mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iTwinEdge].m_iFace == commonFace11)
			v[3] = vSpecial0;

		vector<Vector2D> pointsOnCircle;
		Vector2D center; double R;

		if (fatherEdge[curEdge] == oE0)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original0;
				double x = vec0 * Xaxis0, y = (vec0 ^ Xaxis0).normalize();
				double cosAngle = Xaxis0*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// clock wise
				double x0 = cosAngle*x + sinAngle*y;
				double y0 = -sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0) + Xaxis1_2D*l1);
			}
		}
		else if (fatherEdge[curEdge] == oE1)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original1;
				double x = vec0 * Xaxis1, y = (vec0 ^ Xaxis1).normalize();
				double cosAngle = Xaxis1*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// counter clock wise
				double x0 = cosAngle*x - sinAngle*y;
				double y0 = sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0));
			}
		}
		else if (fatherEdge[curEdge] == oE2)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original2;
				double x = vec0 * Xaxis2, y = -(vec0 ^ Xaxis2).normalize();
				double cosAngle = Xaxis2*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// clock wise
				double x0 = cosAngle*x + sinAngle*y;
				double y0 = -sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0));
			}
		}
		else if (fatherEdge[curEdge] == oE3)
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - original3;
				double x = vec0 * Xaxis3, y = -(vec0 ^ Xaxis3).normalize();
				double cosAngle = Xaxis3*XaxisFinal;
				double sinAngle = sqrt(fabs(1-cosAngle*cosAngle));
				// counter clock wise
				double x0 = cosAngle*x - sinAngle*y;
				double y0 = sinAngle*x + cosAngle*y;
				pointsOnCircle.push_back(Vector2D(x0, y0) + Xaxis2_2D*l2);
			}
		}
		else
		{
			for (unsigned j = 0; j < 4; ++j)
			{
				if (v[j] == iv0 || v[j] == iv1) continue;
				Vector3D curVertPos = mesh->m_pVertex[v[j]].m_vPosition;
				Vector3D vec0 = curVertPos - originFinal;
				double x = vec0 * XaxisFinal, y = (vec0 ^ XaxisFinal).normalize();
				if ((vec0 ^ XaxisFinal) * (Xaxis1 ^ XaxisFinal) < 0) y = -y;
				pointsOnCircle.push_back(Vector2D(x, y));
			}
		}

		CalculateCircumCircle(pointsOnCircle, center, R);
		if (fabs(center.y) < R)
		{
			list< pair<double, double> > curIntervals;
			double width = sqrt(R*R - center.y*center.y);
			if (mesh->m_pEdge[curEdge].m_iVertex[0] == iv0 || 
				mesh->m_pEdge[curEdge].m_iVertex[1] == iv0 || 
				mesh->m_pEdge[curEdge].m_iVertex[0] == iv1 || 
				mesh->m_pEdge[curEdge].m_iVertex[1] == iv1)
			{
				curIntervals.push_back(make_pair(center.x - width, center.x + width));
			}
			else
			{
				curIntervals.push_back(make_pair(leftBound, center.x - width));
				curIntervals.push_back(make_pair(center.x + width, rightBound));
			}
			IntersectWithFeasibleIntervals(feasibleIntervals, curIntervals);
		}
	}

	// if the result interval is not empty
	if (!feasibleIntervals.empty())
	{
		// readjust the postion of iv0
		++mergeVertCnt;
		unsigned ivM = mesh->m_pEdge[iEdge0].m_iVertex[1];
		cout << "Merge vert " << iv0 << " "  << ivM << " " << iv1 << endl;
		double pos = (feasibleIntervals.begin()->first + feasibleIntervals.begin()->second) / 2.0;
		Vector3D newPos = originFinal + XaxisFinal * pos;

		mesh->collapse(iEdge0, (mesh->m_pVertex[iv0].m_vPosition + mesh->m_pVertex[ivM].m_vPosition) / 2);
		mesh->collapse(iEdge1, newPos);
		/*mesh->moveVertTo(iv0, newPos);*/
	}
}

bool DelaunayMesh::InInterval(double pos, pair<double, double> interval)
{
	return pos > interval.first && pos < interval.second;
}