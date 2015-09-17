#include "ExactGeodesic.h"

using namespace MMP_HalfEdge;

ExactGeodesic::ExactGeodesic()
{
	mesh = NULL;
	vertex_src_idx.clear();

	pathbuilt = false;
	winNum = 0;
}

ExactGeodesic::ExactGeodesic(ExactGeodesic* exactGeodesic)
{
	mesh = exactGeodesic->mesh;
	vertex_src_idx.resize(exactGeodesic->vertex_src_idx.size());
	polySource.resize(exactGeodesic->polySource.size());

	for(UINT i = 0;i < vertex_src_idx.size();++i)
		vertex_src_idx[i] = exactGeodesic->vertex_src_idx[i];
	for(UINT i = 0;i < polySource.size();++i)
		polySource[i] = exactGeodesic->polySource[i];

	windowlist.resize(mesh->m_nEdge);
	for(int i = 0;i < mesh->m_nEdge;++i)
		windowlist[i] = exactGeodesic->windowlist[i];
}

ExactGeodesic::~ExactGeodesic()
{
	mesh = NULL;
	windowlist.clear();
}

void ExactGeodesic::AssignMesh(CMesh* _mesh)
{
	mesh = _mesh;

	windowlist.clear();
	windowlist.resize(mesh->m_nEdge);

	vertOccupyied.clear();
	vertOccupyied.resize(mesh->m_nEdge);

	vertexDist.clear();
	vertexDist.resize(mesh->m_nVertex, 1e30);
}

void ExactGeodesic::AssignSrcVertex(UINT src_idx)
{
	vertex_src_idx.push_back(src_idx);
	vertexDist[src_idx] = 0;
}

void ExactGeodesic::RemoveSrcVertex(UINT src_idx)
{
	for (vector<UINT>::iterator iter = vertex_src_idx.begin();
		iter != vertex_src_idx.end(); ++iter)
	{
		if (*iter == src_idx)
		{
			vertex_src_idx.erase(iter);
			break;
		}
	}
}

void ExactGeodesic::AddPolySourceVertex(LineSource s)
{
	polySource.push_back(s);
}

void ExactGeodesic::BuildGeodesicDistField(std::vector<SurfacePoint> *stopPoints)
{
	for (unsigned i = 0; i < mesh->m_nVertex; ++i) vertexDist[i] = 1e30;

	//添加指定的源点(顶点)
	for(UINT i = 0;i < vertex_src_idx.size();++i)
		AddPseudoSource(vertex_src_idx[i], 0);

	for (UINT i = 0; i < polySource.size(); ++i)
	{
		for (UINT j = 0; j < polySource[i].points_on_plat.size(); ++j)
			AddPseudoSource(polySource[i].Face_idx, polySource[i].points_on_plat[j], 0, polySource[i].srcId, polySource[i].srcId);
		/*AddPseudoSource(polySource[i]);*/
	}

	for (UINT i = 0; i < curve.size(); ++i)
	{
		for (UINT j = 0; j < curve[i].points_on_plat.size(); ++j)
			AddPseudoSource(curve[i].Face_idx, curve[i].points_on_plat[j], 0, curve[i].srcId, curve[i].srcId);
/*		AddPseudoSource(curve[i]);*/
	}

	unsigned checkPeriod = 10;
	unsigned iterations = 0;
	while(!windowPriorityQueue.empty())
	{
		Window nextWindow = windowPriorityQueue.pop();

		Propogate(nextWindow);

		if (++iterations % checkPeriod == 0)
			if (StopConditionSatisfied(stopPoints))
				break;
	}

	//Fill out gaps
	for (UINT i = 0; i < mesh->m_nEdge; ++i)
	{
		for (WindowList::iterator iter = windowlist[i].begin();
			iter != windowlist[i].end(); ++iter)
		{
			WindowList::iterator iterNext = iter; ++iterNext;
			double nextStart = mesh->m_pEdge[i].m_length;
			if (iterNext != windowlist[i].end()) nextStart = iterNext->b1;
			if (iter->b2 != nextStart)
			{
				if (iterNext != windowlist[i].end())
				{
					iter->b2 = (iter->b2 + nextStart) / 2;
					iterNext->b1 = iter->b2;
				}
				else
				{
					iter->b2 = nextStart;
				}
			}
		}
	}
	
	//输出windowlist的内容
	
// 	ofstream output("windowlist_mine.txt");
// 	for(UINT i = 0;i < mesh->m_nEdge;++i)
// 	{
// 		if (! (mesh->m_pEdge[i].m_iFace == 34855 && mesh->m_pEdge[mesh->m_pEdge[i].m_iTwinEdge].m_iFace == 34860 || 
// 			mesh->m_pEdge[i].m_iFace == 34860 && mesh->m_pEdge[mesh->m_pEdge[i].m_iTwinEdge].m_iFace == 34855)) continue;
// 		WindowList::iterator it = windowlist[i].begin();
// 		for(;it != windowlist[i].end();++it)
// 		{
// 			output<<"srcId: " << it->srcId <<" idx: "<<it->Edge_idx<<" b1: "<<it->b1<<" b2: "<<it->b2<<" d1: "<<it->d1<<" d2: "<<it->d2<<" sigma: "<<it->sigma<<" windowtype: "<<it->winType<<" ;  "<<endl;
// 		}
// 		output<<endl;
// 	}
// 	ofstream output("windowCalculate.txt", ios::app);
// 	double totalWindowNum = 0;
// 	for(UINT i = 0;i < mesh->m_nEdge;++i)
// 	{
// 		for(WindowList::iterator it = windowlist[i].begin();
// 			it != windowlist[i].end();++it)
// 			++totalWindowNum;
// 	}
// 	output<<"Total window num: "<<totalWindowNum<<endl;
// 	output<<"Average window num: "<<totalWindowNum / mesh->m_nEdge<<endl;
}

double ExactGeodesic::BuildGeodesicPath(const Vector3D& point, const UINT& face_idx, UINT vertex_id)
{
	double pathLen = 0;
	UINT edge_idx;
	double minpos;

	WindowList::iterator cur_window;

	geodesicpath.dst = point;
	geodesicpath.keypoints.clear();
	geodesicpath.passedFaces.clear();
	geodesicpath.srcFace = -1; geodesicpath.srcVert = -1;
	geodesicpath.dstFace = -1; geodesicpath.dstVert = -1;

	bool cur_is_vertex = false;
	UINT vertex_idx;

	Vector3D cur_point;																														//从目的点回溯的当前点
	vector<Vector3D> src_point;
	vector<unsigned> src_vertIds; vector<unsigned> src_faceIds;
	for(unsigned int i = 0;i < vertex_src_idx.size();++i)
	{
		src_point.push_back(mesh->m_pVertex[vertex_src_idx[i]].m_vPosition);											//源点
		src_vertIds.push_back(vertex_src_idx[i]);
	}
	for(unsigned int i = 0;i < curve.size();++i)
	{
		for(unsigned int j = 0;j < curve[i].points_on_plat.size();++j)
		{
			src_point.push_back(curve[i].points_on_plat[j]);
			src_faceIds.push_back(curve[i].Face_idx);
		}
	}

	for(unsigned int i = 0;i < polySource.size();++i)
	{
		for(unsigned int j = 0;j < polySource[i].points_on_plat.size();++j)
		{
			src_point.push_back(polySource[i].points_on_plat[j]);
			src_faceIds.push_back(polySource[i].Face_idx);
		}
	}

	if(vertex_id != -1)																											//如果point是vertex
	{
		cur_is_vertex = true;
		cur_point = point;
		vertex_idx = vertex_id;
		geodesicpath.dst = mesh->m_pVertex[vertex_id].m_vPosition;
		geodesicpath.dstVert = vertex_id;
	}
	else																																				//否则，可能在边上，也可能在面内
	{
		double curPathLen = 1e30;
		SearchWindowOnFace(face_idx, point, curPathLen, edge_idx, minpos);											//无论哪种情况，均可找到测地线与面face_idx的某一边上的交点(尽管point在边上时就是point)
		Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[0]].m_vPosition;
		Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[1]].m_vPosition;
		Vector3D deltav = v1 - v0; deltav.normalize();
		cur_point = v0 + deltav * minpos;
		
		GeodesicKeyNode tmpkeypoint;
		if (EQUALZERO(minpos) || EQUALZERO(mesh->m_pEdge[edge_idx].m_length-minpos)) {
			tmpkeypoint.isVert = true;
			tmpkeypoint.vert_idx = EQUALZERO(minpos) ? mesh->m_pEdge[edge_idx].m_iVertex[0] : 
				mesh->m_pEdge[edge_idx].m_iVertex[1];
			cur_is_vertex = true;
			vertex_idx = tmpkeypoint.vert_idx;
		}
		else  {
			tmpkeypoint.isVert = false;
			cur_is_vertex = false;
		}
		
		tmpkeypoint.edge_idx = edge_idx; tmpkeypoint.x = minpos; geodesicpath.dstFace = face_idx;
		geodesicpath.keypoints.push_back(tmpkeypoint);
		geodesicpath.passedFaces.push_back(mesh->m_pEdge[edge_idx].m_iFace);
	}

	bool isCurPointCloseEnough = false;
	while(!isCurPointCloseEnough){ 
		WindowList::iterator it;

		while (cur_is_vertex) {
			edge_idx = 0; minpos = 0;

			UINT new_vertex_idx = -1;
			cur_is_vertex = vertOccupyied[vertex_idx].isVert;

			if (vertOccupyied[vertex_idx].elementIndex == -1) {
				double dist = 1e30;
				for (unsigned i = 0; i < src_point.size(); ++i) {
					double curDist = (src_point[i] - cur_point).length();
					if (curDist < dist) {
						dist = curDist;
						geodesicpath.src = src_point[i];
						break;
					}
				}
				isCurPointCloseEnough = true;
			}
			if (isCurPointCloseEnough) break;
			if (vertOccupyied[vertex_idx].isVert) {
				new_vertex_idx = vertOccupyied[vertex_idx].elementIndex;
			}
			else {
				edge_idx = vertOccupyied[vertex_idx].elementIndex;
				minpos = vertOccupyied[vertex_idx].pos;
			}

			GeodesicKeyNode gkn;
			gkn.isVert = cur_is_vertex; gkn.vert_idx = new_vertex_idx; gkn.edge_idx = edge_idx; gkn.x = minpos;
			geodesicpath.keypoints.push_back(gkn);
			if (cur_is_vertex) {
				geodesicpath.passedFaces.push_back(mesh->m_pEdge[mesh->m_pVertex[new_vertex_idx].m_piEdge[0]].m_iFace);
			}
			else {
				geodesicpath.passedFaces.push_back(mesh->m_pEdge[edge_idx].m_iFace);
			}

			if (cur_is_vertex) vertex_idx = new_vertex_idx;
		}
		if (isCurPointCloseEnough) break;

		for(it = windowlist[edge_idx].begin(); it != windowlist[edge_idx].end(); ++it)
		{
			if(it->b2 + DOUBLE_EPS < minpos) continue;
			else break;
		}
		if(it == windowlist[edge_idx].end()) {
			Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[0]].m_vPosition;
			Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[1]].m_vPosition;
			Vector3D deltav = v1 - v0; deltav.normalize();
			Vector3D tp = v0 + deltav * minpos;
			double dist = 1e30;
			for (unsigned i = 0; i <src_point.size(); ++i) {
				double curdist = (src_point[i] - tp).length();
				if (curdist < dist) {
					dist = curdist;
					geodesicpath.src = src_point[i];
				}
			}
			break;
		}

		UINT next_edge = mesh->m_pEdge[edge_idx].m_iNextEdge;
		UINT prev_edge = mesh->m_pEdge[next_edge].m_iNextEdge;
		double length1 = mesh->m_pEdge[prev_edge].m_length;
		double length2 = mesh->m_pEdge[edge_idx].m_length;
		double length3 = mesh->m_pEdge[next_edge].m_length;
		double vx = (length1*length1+length2*length2-length3*length3)/(2*length2);
		double vy = sqrt(fabs(length1*length1 - vx*vx));

		Line l2; l2.A = vy; l2.B = length2 - vx; l2.C = -length2*vy;
		Line l3; l3.A = vy; l3.B = - vx; l3.C = 0;

		double wlen = it->b2 - it->b1;
		double kx = minpos;
		double ky = 0;

		//源点和边上点的连线，或者线窗口里的平行线
		Line l1;

		if(it->winType == POINTWINDOW) {
			double sx = (it->d1*it->d1+wlen*wlen-it->d2*it->d2)/(2*wlen) + it->b1;
			double sy = sqrt(fabs(it->d1*it->d1 - (sx-it->b1)*(sx-it->b1)));
			l1.A = sy - ky; l1.B = kx - sx; l1.C = sx*ky - kx*sy;
			unsigned curFace = mesh->m_pEdge[edge_idx].m_iFace;
			unsigned curVerts[3];
			for (unsigned i = 0; i < 3; ++i) curVerts[i] = mesh->m_pFace[curFace].m_piVertex[i];

			bool srcReached = false;
			for (unsigned i = 0; i < src_vertIds.size(); ++i) {
				if (src_vertIds[i] == curVerts[0] || src_vertIds[i] == curVerts[1] || src_vertIds[i] == curVerts[2]) {
					srcReached = true;
					break;
				}
			}
			for (unsigned i = 0; i < src_faceIds.size(); ++i) {
				if (src_faceIds[i] == curFace) {
					srcReached = true;
					break;
				}
			}

			if (it->sigma == 0 && srcReached) {
				Vector3D vecx = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[1]].m_vPosition 
					- mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[0]].m_vPosition;
				vecx.normalize();

				Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[prev_edge].m_iVertex[0]].m_vPosition 
					- mesh->m_pVertex[mesh->m_pEdge[prev_edge].m_iVertex[1]].m_vPosition;
				Vector3D vecy = v0 - (v0*vecx)*vecx;
				vecy.normalize();

				Vector3D pseudoSource = sx*vecx + sy*vecy + mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[0]].m_vPosition;
				UINT vert1 = mesh->m_pEdge[edge_idx].m_iVertex[0];
				UINT vert2 = mesh->m_pEdge[edge_idx].m_iVertex[1];
				UINT opVert = mesh->m_pEdge[next_edge].m_iVertex[1];

				geodesicpath.src = pseudoSource;
				break;
 			}
		}
// 		else if(it->winType == LINEWINDOW) {
// 			double sina = fabs(it->d1-it->d2)/(it->b2-it->b1);
// 			double cosa = sqrt(1-sina*sina);
// 			Vector2D end1, end2;
// 			
// 			if(it->d1 > it->d2) {
// 				l1.A = cosa; l1.B = -sina; l1.C = -cosa*minpos;
// 				end1.x = it->b1 + it->d1/(sqrt(1+l1.A*l1.A/(l1.B*l1.B))); end1.y = -l1.A/l1.B*end1.x-l1.C/l1.B;
// 				end2.x = it->b2 + it->d2/(sqrt(1+l1.A*l1.A/(l1.B*l1.B))); end2.y = -l1.A/l1.B*end2.x-l1.C/l1.B;
// 			}
// 			else {
// 				l1.A = -cosa; l1.B = -sina; l1.C = cosa*minpos;
// 				end1.x = it->b1 - it->d1/(sqrt(1+l1.A*l1.A/(l1.B*l1.B))); end1.y = -l1.A/l1.B*end1.x-l1.C/l1.B;
// 				end2.x = it->b2 - it->d2/(sqrt(1+l1.A*l1.A/(l1.B*l1.B))); end2.y = -l1.A/l1.B*end2.x-l1.C/l1.B;
// 			}
// 			if (InTriangle(end1, Vector2D(0, 0), Vector2D(length2, 0), Vector2D(vx, vy)) && 
// 				InTriangle(end2, Vector2D(0, 0), Vector2D(length2, 0), Vector2D(vx, vy))) {
// 					Line ls; ls.A = end2.y - end1.y; ls.B = end1.x - end2.x; ls.C = end2.x*end1.y-end1.x*end2.y;
// 					double x, y; SolvePlanarEquation(l1, ls, x, y); 
// 
// 					Vector3D vecx = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[1]].m_vPosition 
// 						- mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[0]].m_vPosition;
// 					vecx.normalize();
// 
// 					Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[prev_edge].m_iVertex[0]].m_vPosition 
// 						- mesh->m_pVertex[mesh->m_pEdge[prev_edge].m_iVertex[1]].m_vPosition;
// 					Vector3D vecy = v0 - (v0*vecx)*vecx;
// 					vecy.normalize();
// 
// 					geodesicpath.src = x*vecx + y*vecy + mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[0]].m_vPosition;
// 					break;
// 			}
// 		}

		Vector3D next_point = cur_point;
		double x1, y1, x2, y2 = 0;
		double tmp_minpos;
		UINT tmp_edge_idx = -1;
		double angle1 = 0, angle2 = 0;

		if(SolvePlanarEquation(l1, l2, x1, y1) && y1 > -DOUBLE_EPS && y1 < vy + DOUBLE_EPS) {
			tmp_edge_idx = mesh->m_pEdge[next_edge].m_iTwinEdge;
			tmp_minpos = (vy-y1)/vy*length3;
			Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[tmp_edge_idx].m_iVertex[0]].m_vPosition;
			Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[tmp_edge_idx].m_iVertex[1]].m_vPosition;
			Vector3D deltav = v1 - v0; deltav.normalize();
			next_point = v0 + deltav * tmp_minpos;
		}
		else if(SolvePlanarEquation(l1, l3, x2, y2) && y2 > -DOUBLE_EPS && y2 < vy + DOUBLE_EPS) {
			tmp_edge_idx = mesh->m_pEdge[prev_edge].m_iTwinEdge;
			tmp_minpos = y2/vy*length1;
			Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[tmp_edge_idx].m_iVertex[0]].m_vPosition;
			Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[tmp_edge_idx].m_iVertex[1]].m_vPosition;
			Vector3D deltav = v1 - v0; deltav.normalize();
			next_point = v0 + deltav * tmp_minpos;
		}

		Vector2D vec0(minpos-vx, -vy); Vector2D vec1(-vx, -vy); Vector2D vec2(length2-vx, -vy);
		vec0.normalize(); vec1.normalize(); vec2.normalize();
		angle1 = acos(vec0*vec1); angle2 = acos(vec0*vec2);

		assert(tmp_edge_idx != -1);
		cur_point = next_point;
		minpos = tmp_minpos;
		edge_idx = tmp_edge_idx;
		vertex_idx = -1;

		if(EQUALZERO(minpos) || EQUALZERO(minpos-mesh->m_pEdge[edge_idx].m_length)) {
			cur_is_vertex = true;
			if(EQUALZERO(minpos)) vertex_idx = mesh->m_pEdge[edge_idx].m_iVertex[0];
			else vertex_idx = mesh->m_pEdge[edge_idx].m_iVertex[1];
		}
		else cur_is_vertex = false;

		GeodesicKeyNode tmpkeypoint;
		tmpkeypoint.isVert = cur_is_vertex;
		tmpkeypoint.vert_idx = vertex_idx;
		tmpkeypoint.edge_idx = edge_idx;
		tmpkeypoint.x = minpos;
		geodesicpath.keypoints.push_back(tmpkeypoint);
		geodesicpath.passedFaces.push_back(mesh->m_pEdge[edge_idx].m_iFace);

		if (cur_is_vertex) {
			for (unsigned i = 0; i < src_vertIds.size(); ++i) {
				if (vertex_idx == src_vertIds[i]) {
					isCurPointCloseEnough = true;
					geodesicpath.src = mesh->m_pVertex[vertex_idx].m_vPosition;
					break;
				}
			}
		}
	}

	pathbuilt = true;

	pathLen = 0;
	for (list<GeodesicKeyNode>::iterator iter = geodesicpath.keypoints.begin();
		iter != geodesicpath.keypoints.end(); ++iter) {
			list<GeodesicKeyNode>::iterator iterNext = iter; ++iterNext;

			Vector3D r0;
			if (iter->isVert) {
				r0 = mesh->m_pVertex[iter->vert_idx].m_vPosition;
			}
			else {
				Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[iter->edge_idx].m_iVertex[0]].m_vPosition;
				Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[iter->edge_idx].m_iVertex[1]].m_vPosition;
				Vector3D unv = v1 - v0; unv.normalize();
				r0 = v0 + iter->x * unv;
			}

			if (iter == geodesicpath.keypoints.begin()) {
				pathLen += (geodesicpath.dst-r0).length();
			}

			if (iterNext != geodesicpath.keypoints.end()) {
				Vector3D r1;
				if (iterNext->isVert) {
					r1 = mesh->m_pVertex[iterNext->vert_idx].m_vPosition;
				}
				else {
					Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[iterNext->edge_idx].m_iVertex[0]].m_vPosition;
					Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[iterNext->edge_idx].m_iVertex[1]].m_vPosition;
					Vector3D unv = v1 - v0; unv.normalize();
					r1 = v0 + iterNext->x * unv;
				}
				pathLen += (r0-r1).length();
			}
			else {
				pathLen += (geodesicpath.src-r0).length();
			}
	}

	return pathLen;
}

double ExactGeodesic::CalcGeodesicLenOnly(const Vector3D& point, const UINT& face_idx, UINT vertex_id, UINT &srcId)
{
	double pathlen = 1e30;
	UINT edge_idx;
	double minpos;

	WindowList::iterator cur_window;

	if(vertex_id != -1)																											//如果point是vertex
	{
		for(int i = 0;i < mesh->m_pVertex[vertex_id].m_nValence;++i)				//就要从所有以顶点为端点的窗口中选一个最小的出来
		{
			edge_idx = mesh->m_pVertex[vertex_id].m_piEdge[i];
			UINT twin_edge_idx = mesh->m_pEdge[edge_idx].m_iTwinEdge;

			cur_window = windowlist[edge_idx].begin();
			if(!windowlist[edge_idx].empty())
			{
				double tmplen = cur_window->d1 + cur_window->sigma;
				if(tmplen + DOUBLE_EPS < pathlen)	
				{
					pathlen = tmplen;
					srcId = cur_window->srcId;
				}
			}

			if(twin_edge_idx == -1 || windowlist[twin_edge_idx].empty()) continue;
			cur_window = windowlist[twin_edge_idx].end();
			--cur_window;
			double tmplen = cur_window->d2 + cur_window->sigma;
			if(tmplen + DOUBLE_EPS < pathlen)	
			{
				pathlen = tmplen;
				srcId = cur_window->srcId;
			}
		}
	}
	else																																				//否则，可能在边上，也可能在面内
	{
		SearchWindowOnFace(face_idx, point, pathlen, edge_idx, minpos);											//无论哪种情况，均可找到测地线与面face_idx的某一边上的交点(尽管point在边上时就是point)
	}
	return pathlen;
}

double ExactGeodesic::CalcGeodesicDistToVertex(UINT vertex_id, UINT &srcId)
{
	if (vertex_id >= mesh->m_nVertex) return 1e30;

	return CalcGeodesicLenOnly(mesh->m_pVertex[vertex_id].m_vPosition, 
		0, 
		vertex_id, srcId);
}

void ExactGeodesic::Clear()
{
	windowPriorityQueue.clear();

	windowlist.clear();
	windowlist.resize(mesh->m_nEdge);
	vertOccupyied.clear();
	vertOccupyied.resize(mesh->m_nVertex);
	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
		vertexDist[i] = 1e30;

// 	vertex_src_idx.clear();
// 	polySource.clear();
// 	curve.clear();

	ClearGeodesicPath();
}

void ExactGeodesic::ClearGeodesicPath()
{
	geodesicpath.keypoints.clear();
	geodesicpath.passedFaces.clear();
	pathbuilt = false;
}

void ExactGeodesic::FillInVertDist(map<UINT, vector<UINT> > &region2Faces)
{
// 	maxDist = 0;
// 	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
// 	{
// 		if (vertexDist[i] > 1e20)
// 		{
// 			for (unsigned j = 0; j < mesh->m_pVertex[i].m_nValence; ++j)
// 			{
// 				unsigned adjEdge = mesh->m_pVertex[i].m_piEdge[j];
// 				unsigned opVert = mesh->m_pEdge[adjEdge].m_iVertex[0] == i ? 
// 					mesh->m_pEdge[adjEdge].m_iVertex[1] : mesh->m_pEdge[adjEdge].m_iVertex[0];
// 				vertexDist[i] = vertexDist[i] < vertexDist[opVert] + mesh->m_pEdge[adjEdge].m_length ? 
// 					vertexDist[i] : vertexDist[opVert] + mesh->m_pEdge[adjEdge].m_length;
// 			}
// 		}
// 		maxDist = maxDist < vertexDist[i] ? vertexDist[i] : maxDist;
// 	}
	vertexDist.resize(mesh->m_nVertex);
	maxDist = 0;
	UINT srcId;
	vector<UINT> vert2Region; vert2Region.resize(mesh->m_nVertex);
	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
	{
		vertexDist[i] = CalcGeodesicDistToVertex(i, srcId);
		vert2Region[i] = srcId;
		if (vertexDist[i] < 1e20 && vertexDist[i] > maxDist) maxDist = vertexDist[i];
	}
	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
	{
		if (vertexDist[i] > 1e20)
		{
			unsigned minNeigh = -1;
			for (unsigned j = 0; j < mesh->m_pVertex[i].m_nValence; ++j)
			{
				unsigned adjEdge = mesh->m_pVertex[i].m_piEdge[j];
				unsigned opVert = mesh->m_pEdge[adjEdge].m_iVertex[0] == i ? 
					mesh->m_pEdge[adjEdge].m_iVertex[1] : mesh->m_pEdge[adjEdge].m_iVertex[0];
				if (vertexDist[i] > vertexDist[opVert] + mesh->m_pEdge[adjEdge].m_length)
				{
					vertexDist[i] = vertexDist[opVert] + mesh->m_pEdge[adjEdge].m_length;
					vert2Region[i] = vert2Region[opVert];
				}
			}
		}
	}

	for (unsigned i = 0; i < mesh->m_nFace; ++i) {
		vector<unsigned> addRegion;
		for (unsigned j = 0; j < 3; ++j) {
			bool flag = false;
			for (unsigned k = 0; k < addRegion.size(); ++k) {
				if (addRegion[k] == vert2Region[mesh->m_pFace[i].m_piVertex[j]]) {
					flag = true;
					break;
				}
			}
			if (flag) continue;
			addRegion.push_back(vert2Region[mesh->m_pFace[i].m_piVertex[j]]);
			region2Faces[vert2Region[mesh->m_pFace[i].m_piVertex[j]]].push_back(i);
		}
	}
}

void ExactGeodesic::SplitWindows(unsigned edgeIndex)
{
	windowlist.resize(mesh->m_nEdge);
	unsigned subEdge = mesh->m_pEdge[mesh->m_pEdge[edgeIndex].m_iNextEdge].m_iTwinEdge;
	subEdge = mesh->m_pEdge[subEdge].m_iNextEdge;
	unsigned subTwinEdge = mesh->m_pEdge[subEdge].m_iTwinEdge;
	unsigned twinEdge = mesh->m_pEdge[edgeIndex].m_iTwinEdge;
	double pos = mesh->m_pEdge[edgeIndex].m_length;

	for (WindowList::iterator iter = windowlist[edgeIndex].begin(); iter != windowlist[edgeIndex].end(); )
	{
		if (iter->b2 < pos)
		{
			++iter; continue;
		}
		WindowList::iterator iterNext = iter; ++iterNext;

		if (iter->b1 < pos)
		{
			Window headWin, tailWin;
			CutWindow(*iter, iter->b1, pos, headWin);
			CutWindow(*iter, pos, iter->b2, tailWin);
			tailWin.Edge_idx = subEdge; tailWin.b1 -= pos; tailWin.b2 -= pos;

			windowlist[edgeIndex].insert(iterNext, headWin);
			windowlist[subEdge].push_back(tailWin);
			windowlist[edgeIndex].erase(iter);
			iter = iterNext;
			continue;
		}

		iter->Edge_idx = subEdge; iter->b1 -= pos; iter->b2 -= pos;
		windowlist[subEdge].push_back(*iter);
		windowlist[edgeIndex].erase(iter);
		iter = iterNext;
	}

	// reverse windows on edge subEdge to its twin edge
	if (subTwinEdge == -1) return;

	pos = mesh->m_pEdge[subTwinEdge].m_length;
	for (WindowList::iterator iter = windowlist[subTwinEdge].begin(); iter != windowlist[subTwinEdge].end(); )
	{
		if (iter->b2 < pos)
		{
			++iter; continue;
		}
		WindowList::iterator iterNext = iter; ++iterNext;

		if (iter->b1 < pos)
		{
			Window headWin, tailWin;
			CutWindow(*iter, iter->b1, pos, headWin);
			CutWindow(*iter, pos, iter->b2, tailWin);
			tailWin.Edge_idx = twinEdge; tailWin.b1 -= pos; tailWin.b2 -= pos;

			windowlist[subTwinEdge].insert(iterNext, headWin);
			windowlist[twinEdge].push_back(tailWin);
			windowlist[subTwinEdge].erase(iter);
			iter = iterNext;
			continue;
		}

		iter->Edge_idx = twinEdge; iter->b1 -= pos; iter->b2 -= pos;
		windowlist[twinEdge].push_back(*iter);
		windowlist[subTwinEdge].erase(iter);
		iter = iterNext;
	}
}

void ExactGeodesic::UpdateAroundVert(unsigned centerVert)
{
	vertexDist.resize(mesh->m_nVertex, 1e30);
	AddPseudoSource(centerVert, 0.0);
	while(!windowPriorityQueue.empty())
	{
		Window nextWindow = windowPriorityQueue.pop();
		Propogate(nextWindow);
	}
}

bool ExactGeodesic::CheckVertexOnEdge(const Window& w, const double& pos)
{
	UINT edge_idx = w.Edge_idx;
	UINT next_idx = mesh->m_pEdge[edge_idx].m_iNextEdge;
	UINT prev_idx = mesh->m_pEdge[next_idx].m_iNextEdge;

	for (UINT i = 0; i < curve.size(); ++i)
		if (mesh->m_pEdge[edge_idx].m_iFace == curve[i].Face_idx) return true;
	for (unsigned i = 0; i < polySource.size(); ++i)
		if (mesh->m_pEdge[edge_idx].m_iFace == polySource[i].Face_idx) return true;

	double wlen = w.b2 - w.b1;
	double sx = (w.d1*w.d1 + wlen*wlen - w.d2*w.d2)/(2*wlen) + w.b1;
	double sy = sqrt(fabs(w.d1*w.d1 - (sx-w.b1)*(sx-w.b1)));
	double length1 = mesh->m_pEdge[prev_idx].m_length;
	double length2 = mesh->m_pEdge[edge_idx].m_length;
	double length3 = mesh->m_pEdge[next_idx].m_length;
	double vx = (length1*length1 + length2*length2 - length3*length3)/(2*length2);
	double vy = sqrt(fabs(length1*length1-vx*vx));
	Line l1, l2;
	l1.A = vy, l1.B = -vx, l1.C = 0;
	l2.A = vy, l2.B = length2 - vx, l2.C = -length2*vy;
	Line ls;
	ls.A = sy;ls.B = pos-sx;ls.C = -pos*sy;
	/*
	if(EQUALZERO(pos))
	{
		ls.A = sy;
		ls.B = -sx;
		ls.C = 0;
	}
	else if(EQUALZERO(length2-pos))
	{
		ls.A = sy;
		ls.B = length2-sx;
		ls.C = -length2*sy;
	}
	*/
	double xx, yy;
	if(SolvePlanarEquation(ls, l1, xx, yy))
	{
// 		if(yy < vy + DOUBLE_EPS && yy > -DOUBLE_EPS)
// 			return true;
		if(EQUALZERO(pos))
		{
			if(yy < vy + DOUBLE_EPS && yy > DOUBLE_EPS)
				return true;
		}
		else
		{
			if(yy < vy + DOUBLE_EPS && yy > -DOUBLE_EPS)
				return true;
		}
	}
	if(SolvePlanarEquation(ls, l2, xx, yy))
	{
// 		if(yy < vy + DOUBLE_EPS && yy > -DOUBLE_EPS)
// 			return true;
		if(EQUALZERO(pos-length2))
		{
			if(yy < vy + DOUBLE_EPS && yy > DOUBLE_EPS)
				return true;
		}
		else
		{
			if(yy < vy + DOUBLE_EPS && yy > -DOUBLE_EPS)
				return true;
		}
	}
	return false;
}

void ExactGeodesic::SearchWindowOnFace(const UINT& face_idx, const Vector3D& point, double& pathlen, UINT& edge_idx, double& minpos)
{
	for(int i = 0;i <  mesh->m_pFace[face_idx].m_nType;++i)
	{
		int cur_edge = mesh->m_pFace[face_idx].m_piEdge[i];
		//edgev为x轴坐标
		Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[cur_edge].m_iVertex[0]].m_vPosition;
		Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[cur_edge].m_iVertex[1]].m_vPosition;
		Vector3D edgev = v1 - v0;
		//point在所在面上的向量坐标
		Vector3D tmpv = point - v0;
		//ex, ey为所在面上point的坐标
		double ex = (tmpv * edgev) / edgev.length();
		double ey = (tmpv ^ edgev).length() / edgev.length();

		WindowList::iterator it;
		for(it = windowlist[cur_edge].begin();it != windowlist[cur_edge].end();++it)
		{
			//计算point和当前窗口对应伪源点之间的连线与窗口所在边的交点
			double x_pro = (it->d1*it->d1 + (it->b2-it->b1)*(it->b2-it->b1) - it->d2*it->d2)/(2*(it->b2-it->b1));
			double sx = it->b1 + x_pro;
			double sy = sqrt(fabs(it->d1*it->d1 - x_pro*x_pro));																						//???

			double curlen, pos;
			FindShortestPathThroughWindow(it, ex, ey, pos);
			if(it->winType == POINTWINDOW)
				curlen = sqrt((pos-sx)*(pos-sx)+sy*sy) + sqrt((pos-ex)*(pos-ex)+ey*ey) + it->sigma;
			else
				curlen = sqrt((pos-ex)*(pos-ex)+ey*ey)+(it->d2-it->d1)/(it->b2-it->b1)*pos+(it->b2*it->d1-it->b1*it->d2)/(it->b2-it->b1);
			if(curlen + DOUBLE_EPS < pathlen)
			{
				pathlen = curlen;
				minpos = pos;
				edge_idx = cur_edge;
			}
		}

		int cur_twin_edge = mesh->m_pEdge[cur_edge].m_iTwinEdge;
		ex = mesh->m_pEdge[cur_edge].m_length - ex;																																//转换到twin edge的坐标
		ey = -ey;

		for(it = windowlist[cur_twin_edge].begin();it != windowlist[cur_twin_edge].end();++it)
		{
			//计算point和当前窗口对应伪源点之间的连线与窗口所在边的交点
			double x_pro = (it->d1*it->d1 + (it->b2-it->b1)*(it->b2-it->b1) - it->d2*it->d2)/(2*(it->b2-it->b1));
			double sx = it->b1 + x_pro;
			double sy = sqrt(fabs(it->d1*it->d1 - x_pro*x_pro));

			double curlen, pos;
			FindShortestPathThroughWindow(it, ex, ey, pos);
			if(it->winType == POINTWINDOW)
				curlen = sqrt((pos-sx)*(pos-sx)+sy*sy) + sqrt((pos-ex)*(pos-ex)+ey*ey) + it->sigma;
			else
				curlen = sqrt((pos-ex)*(pos-ex)+ey*ey)+(it->d2-it->d1)/(it->b2-it->b1)*pos+(it->b2*it->d1-it->b1*it->d2)/(it->b2-it->b1);
			if(curlen + DOUBLE_EPS< pathlen)
			{
				pathlen = curlen;
				minpos = pos;
				edge_idx = cur_twin_edge;
			}
		}
	}
}

void ExactGeodesic::FindShortestPathThroughWindow(WindowList::iterator win, double x2, double y2, double& pos)
{
	if(win->winType == POINTWINDOW)																																				//圆窗口
	{
		double x_pro = (win->d1*win->d1 + (win->b2-win->b1)*(win->b2-win->b1) - win->d2*win->d2)/(2*(win->b2-win->b1));
		double x1 = win->b1 + x_pro;
		double y1 = sqrt(fabs(win->d1*win->d1 - x_pro*x_pro));

		double A = y2*y2 - y1*y1;
		double B = -2*(x1*y2*y2-x2*y1*y1);
		double C = x1*x1*y2*y2 - x2*x2*y1*y1;

		if(EQUALZERO(A))
		{
			if(!EQUALZERO(B))
			{
				double r = -C/B;
				if(r > win->b1 && r < win->b2)
				{
					double sax = (r+win->b1)/2;
					double fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					sax = (r+win->b2)/2;
					double fd2 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd1 > 0 && fd2 > 0)
						pos = win->b1;
					else if(fd1 > 0 && fd2 < 0)
					{
						double v1 = sqrt((win->b1-x1)*(win->b1-x1)+y1*y1) + sqrt((win->b1-x2)*(win->b1-x2)+y2*y2);
						double v2 = sqrt((win->b2-x1)*(win->b2-x1)+y1*y1) + sqrt((win->b2-x2)*(win->b2-x2)+y2*y2);
						if(v1 < v2) pos = win->b1;
						else pos = win->b2;
					}
					else if(fd1 < 0 && fd2 > 0)
						pos = r;
					else if(fd1 < 0 && fd2 < 0)
						pos = win->b2;
				}
				else
				{
					double sax = (win->b1+win->b2)/2;
					double fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd1 >= 0)
						pos = win->b1;
					else
						pos = win->b2;
				}
			}
			else if(EQUALZERO(C))																			//处处相等
			{
				pos = (win->b1+win->b2)/2;
			}
			else
			{
				double sax = (win->b1+win->b2)/2;
				double fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
				if(fd1 >= 0)
					pos = win->b1;
				else
					pos = win->b2;
			}
		}
		else
		{
			double r1 = (-B - sqrt(fabs(B*B - 4*A*C)))/(2*A);
			double r2 = (-B + sqrt(fabs(B*B - 4*A*C)))/(2*A);

			if (EQUALZERO(r1-r2)) 
			{
				if(r1 > win->b1 && r1 < win->b2)
				{
					double sax = (r1+win->b1)/2;
					double fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					sax = (r1+win->b2)/2;
					double fd2 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd1 > 0 && fd2 > 0)
						pos = win->b1;
					else if(fd1 > 0 && fd2 < 0)
					{
						double v1 = sqrt((win->b1-x1)*(win->b1-x1)+y1*y1) + sqrt((win->b1-x2)*(win->b1-x2)+y2*y2);
						double v2 = sqrt((win->b2-x1)*(win->b2-x1)+y1*y1) + sqrt((win->b2-x2)*(win->b2-x2)+y2*y2);
						if(v1 < v2) pos = win->b1;
						else pos = win->b2;
					}
					else if(fd1 < 0 && fd2 > 0)
						pos = r1;
					else if(fd1 < 0 && fd2 < 0)
						pos = win->b2;
				}
				else
				{
					double sax = (win->b1+win->b2)/2;
					double fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd1 >= 0)
						pos = win->b1;
					else
						pos = win->b2;
				}
				return;
			}

			if(r1 > r2)
			{
				double tmp = r1;
				r1 = r2;
				r2 = tmp;
			}

			double sax, fd, fd1, fd2, fd3 = 0;

			if(r1 <= win->b1)
			{
				if(r2 <= win->b1 || r2 >= win->b2)
				{
					sax = (win->b1+win->b2)/2;
					fd = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd >= 0)
						pos = win->b1;
					else
						pos = win->b2;
				}
				else if(r2 > win->b1 && r2 < win->b2)
				{
					sax = (win->b1+r2)/2;
					fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					sax = (r2+win->b2)/2;
					fd2 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd1 > 0 && fd2 > 0)
						pos = win->b1;
					else if(fd1 > 0 && fd2 < 0)
					{
						double v1 = sqrt((win->b1-x1)*(win->b1-x1)+y1*y1) + sqrt((win->b1-x2)*(win->b1-x2)+y2*y2);
						double v2 = sqrt((win->b2-x1)*(win->b2-x1)+y1*y1) + sqrt((win->b2-x2)*(win->b2-x2)+y2*y2);
						if(v1 < v2) pos = win->b1;
						else pos = win->b2;
					}
					else if(fd1 < 0 && fd2 > 0)
						pos = r2;
					else if(fd1 < 0 && fd2 < 0)
						pos = win->b2;
				}
			}
			else if(r1 > win->b1 && r1 < win->b2)
			{
				if(r2 > win->b1 && r2 < win->b2)
				{
					sax = (win->b1+r1)/2;
					fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					sax = (r1+r2)/2;
					fd2 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					sax = (r2+win->b2)/2;
					fd3 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd1 > 0 && fd2 > 0 && fd3 > 0)
					{
						pos = win->b1;
					}
					else if(fd1 > 0 && fd2 > 0 && fd3 < 0 || fd1 > 0 && fd2 < 0 && fd3 < 0)
					{
						double v1 = sqrt((win->b1-x1)*(win->b1-x1)+y1*y1) + sqrt((win->b1-x2)*(win->b1-x2)+y2*y2);
						double v2 = sqrt((win->b2-x1)*(win->b2-x1)+y1*y1) + sqrt((win->b2-x2)*(win->b2-x2)+y2*y2);
						if(v1 < v2) pos = win->b1;
						else pos = win->b2;
					}
					else if(fd1 > 0 && fd2 < 0 && fd3 > 0)
					{
						double v1 = sqrt((win->b1-x1)*(win->b1-x1)+y1*y1) + sqrt((win->b1-x2)*(win->b1-x2)+y2*y2);
						double v2 = sqrt((r2-x1)*(r2-x1)+y1*y1) + sqrt((r2-x2)*(r2-x2)+y2*y2);
						if(v1 < v2) pos = win->b1;
						else pos = r2;
					}
					else if(fd1 < 0 && fd2 > 0 && fd3 > 0)
					{
						pos = r1;
					}
					else if(fd1 < 0 && fd2 > 0 && fd3 < 0)
					{
						double v1 = sqrt((r1-x1)*(r1-x1)+y1*y1) + sqrt((r1-x2)*(r1-x2)+y2*y2);
						double v2 = sqrt((win->b2-x1)*(win->b2-x1)+y1*y1) + sqrt((win->b2-x2)*(win->b2-x2)+y2*y2);
						if(v1 < v2) pos = r1;
						else pos = win->b2;
					}
					else if(fd1 < 0 && fd2 < 0 && fd3 > 0)
					{
						pos = r2;
					}
					else if(fd1 < 0 && fd2 < 0 && fd3 < 0)
					{
						pos = win->b2;
					}
				}
				else
				{
					sax = (win->b1+r1)/2;
					fd1 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					sax = (r1+win->b2)/2;
					fd2 = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
					if(fd1 > 0 && fd2 > 0)
						pos = win->b1;
					else if(fd1 > 0 && fd2 < 0)
					{
						double v1 = sqrt((win->b1-x1)*(win->b1-x1)+y1*y1) + sqrt((win->b1-x2)*(win->b1-x2)+y2*y2);
						double v2 = sqrt((win->b2-x1)*(win->b2-x1)+y1*y1) + sqrt((win->b2-x2)*(win->b2-x2)+y2*y2);
						if(v1 < v2) pos = win->b1;
						else pos = win->b2;
					}
					else if(fd1 < 0 && fd2 > 0)
						pos = r1;
					else if(fd1 < 0 && fd2 < 0)
						pos = win->b2;
				}
			}
			else
			{
				sax = (win->b1+win->b2)/2;
				fd = (sax-x1)/sqrt((sax-x1)*(sax-x1)+y1*y1) + (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2);
				if(fd >= 0)
					pos = win->b1;
				else
					pos = win->b2;
			}
		}
	}
	else																																																//线窗口
	{
		double m = (win->d2-win->d1)/(win->b2-win->b1);
		double n = (win->b2*win->d1-win->b1*win->d2)/(win->b2-win->b1);
		double r1 = x2+sqrt(m*m/(1-m*m))*y2;
		double r2 = x2-sqrt(m*m/(1-m*m))*y2;

		if(r1 > r2)
		{
			double tmp = r2;
			r2 = r1;
			r1 = tmp;
		}

		if(r1 < win->b1)
		{
			if(r2 < win->b1 || r2 > win->b2)
			{
				double sax = (win->b1+win->b2)/2;
				double fd = (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2)+m;
				if(fd < 0) pos = win->b2;
				else pos = win->b1;
			}
			else if(r2 > win->b1 && r2 < win->b2)
			{
				double sax1 = (win->b1+r2)/2;
				double sax2 = (r2+win->b2)/2;
				double fd1 = (sax1-x2)/sqrt((sax1-x2)*(sax1-x2)+y2*y2)+m;
				double fd2 = (sax2-x2)/sqrt((sax2-x2)*(sax2-x2)+y2*y2)+m;
				if(fd1 < 0 && fd2 < 0) pos = win->b2;
				else if(fd1 < 0 && fd2 > 0) pos = r2;
				else if(fd1 > 0 && fd2 < 0)
					pos = sqrt((win->b1-x2)*(win->b1-x2)+y2*y2)+m*win->b1+n > sqrt((win->b2-x2)*(win->b2-x2)+y2*y2)+m*win->b2+n ? win->b2 : win->b1;
				else if(fd1 > 0 && fd2 > 0)
					pos = win->b1;
			}
		}
		else if(r1 > win->b1 && r1 < win->b2)
		{
			if(r2 < win->b2)
			{
				double sax1 = (win->b1+r1)/2;
				double sax2 = (r1+r2)/2;
				double sax3 = (r2+win->b2)/2;

				double fd1 = (sax1-x2)/sqrt((sax1-x2)*(sax1-x2)+y2*y2)+m;
				double fd2 = (sax2-x2)/sqrt((sax2-x2)*(sax2-x2)+y2*y2)+m;
				double fd3 = (sax3-x2)/sqrt((sax3-x2)*(sax3-x2)+y2*y2)+m;

				if(fd1 < 0 && fd2 < 0 && fd3 < 0)
					pos = win->b2;
				else if(fd1 < 0 && fd2 < 0 && fd3 > 0)
					pos = r2;
				else if(fd1 < 0 && fd2 > 0 && fd3 < 0)
					pos = sqrt((r1-x2)*(r1-x2)+y2*y2)+m*r1+n > sqrt((win->b2-x2)*(win->b2-x2)+y2*y2)+m*win->b2+n ? win->b2 : r1;
				else if(fd1 < 0 && fd2 > 0 && fd3 > 0)
					pos = r1;
				else if(fd1 > 0 && fd2 < 0 && fd3 < 0 || fd1 > 0 && fd2 > 0 && fd3 < 0)
					pos = sqrt((win->b1-x2)*(win->b1-x2)+y2*y2)+m*win->b1+n > sqrt((win->b2-x2)*(win->b2-x2)+y2*y2)+m*win->b2+n ? win->b2 : win->b1;
				else if(fd1 > 0 && fd2 < 0 && fd3 > 0)
					pos = sqrt((win->b1-x2)*(win->b1-x2)+y2*y2)+m*win->b1+n > sqrt((r2-x2)*(r2-x2)+y2*y2)+m*r2+n ? r2 : win->b1;
				else if(fd1 > 0 && fd2 > 0 && fd3 > 0)
					pos = win->b1;
			}
			else
			{
				double sax1 = (win->b1+r1)/2;
				double sax2 = (r1+win->b2)/2;
				double fd1 = (sax1-x2)/sqrt((sax1-x2)*(sax1-x2)+y2*y2)+m;
				double fd2 = (sax2-x2)/sqrt((sax2-x2)*(sax2-x2)+y2*y2)+m;
				if(fd1 < 0 && fd2 < 0) pos = win->b2;
				else if(fd1 < 0 && fd2 > 0) pos = r1;
				else if(fd1 > 0 && fd2 < 0)
					pos = sqrt((win->b1-x2)*(win->b1-x2)+y2*y2)+m*win->b1+n > sqrt((win->b2-x2)*(win->b2-x2)+y2*y2)+m*win->b2+n ? win->b2 : win->b1;
				else if(fd1 > 0 && fd2 > 0)
					pos = win->b1;
			}
		}
		else
		{
			double sax = (win->b1+win->b2)/2;
			double fd = (sax-x2)/sqrt((sax-x2)*(sax-x2)+y2*y2)+m;
			if(fd < 0) pos = win->b2;
			else pos = win->b1;
		}
	}
}

void ExactGeodesic::Propogate(const Window& win)
{
	switch(win.winType)
	{
	case POINTWINDOW:
		PropogatePointWindow(win);
		break;
	case LINEWINDOW:
		PropogateLineWindow(win);
		break;
	}
}

void ExactGeodesic::PropogatePointWindow(const Window& win)
{
	double s = win.b2 - win.b1;

	double sx = win.b1 + (win.d1*win.d1 + s*s - win.d2*win.d2) / (2*s);
	double sy = sqrt(fabs(win.d1*win.d1 - (sx-win.b1)*(sx-win.b1)));

	//连接win.b1和伪源点的直线Ax+By+C=0
	Line r1, r2;
	r1.A = sy;
	r1.B = win.b1 - sx;
	r1.C = -win.b1*sy;

	r2.A = sy;
	r2.B = win.b2 - sx;
	r2.C = -win.b2*sy;

	int edge_idx = mesh->m_pEdge[win.Edge_idx].m_iTwinEdge;

	if(edge_idx == -1) return;

	UINT next_edge_idx = mesh->m_pEdge[edge_idx].m_iNextEdge;
	UINT next_next_edge_idx = mesh->m_pEdge[next_edge_idx].m_iNextEdge;

	UINT vert[3];
	vert[0] = mesh->m_pEdge[edge_idx].m_iVertex[1];
	vert[1] = mesh->m_pEdge[next_edge_idx].m_iVertex[1];
	vert[2] = mesh->m_pEdge[edge_idx].m_iVertex[0];

	double len1 = mesh->m_pEdge[edge_idx].m_length;
	double len2 = mesh->m_pEdge[next_edge_idx].m_length;
	double len3 = mesh->m_pEdge[next_next_edge_idx].m_length;

	double x_project = (len1*len1 + len2*len2 - len3*len3) / (2*len1);
	double y = -sqrt(fabs(len2*len2 - x_project*x_project));


	Line e1, e2;
	e1.A = y; e1.B = -x_project; e1.C = 0;
	e2.A = y; e2.B = len1 - x_project; e2.C = -len1*y;

	Line l;
	l.A = sy - y; l.B = x_project - sx; l.C = sx*y - x_project*sy;
	double x = -l.C/l.A;				//源点与三角形另一顶点连线与x轴的交点(l.A可能为0吗?)

	double xx, yy;
	double a, b, d1, d2;
	Window newwindow;


	UINT edges[2];
	edges[0] = next_edge_idx; edges[1] = next_next_edge_idx;
	/*
	if(isSaddle[mesh->m_pEdge[next_edge_idx].m_iVertex[0]] &&
	EQUALZERO(win.b1) && 
	x < win.b1 + DOUBLE_EPS)
	*/
	if(mesh->m_pAngles[vert[0]] > 2*PI - DOUBLE_EPS && EQUALZERO(win.b1) && 
		//!EQUALZERO(win.d1) && 
		(win.d1*win.d1+(win.b2-win.b1)*(win.b2-win.b1)-win.d2*win.d2)/(win.d1*(win.b2-win.b1))+(len1*len1+len2*len2-len3*len3)/(len1*len2)<DOUBLE_EPS)
	{
		AddPseudoSource(mesh->m_pEdge[next_edge_idx].m_iVertex[0], win.sigma + win.d1, edges, win.srcId, mesh->m_pEdge[next_edge_idx].m_iVertex[0]);
	}
	/*
	else if(isSaddle[mesh->m_pEdge[next_next_edge_idx].m_iVertex[1]] &&
	EQUALZERO(win.b2-len1) && 
	x + DOUBLE_EPS > win.b2)
	*/
	else if(mesh->m_pAngles[vert[2]] > 2*PI - DOUBLE_EPS && EQUALZERO(win.b2-len1) && 
		//!EQUALZERO(win.d2) && 
		((win.b2-win.b1)*(win.b2-win.b1)+win.d2*win.d2-win.d1*win.d1)/(win.d2*(win.b2-win.b1))+(len1*len1+len3*len3-len2*len2)/(len1*len3)<DOUBLE_EPS)
	{
		AddPseudoSource(mesh->m_pEdge[next_next_edge_idx].m_iVertex[1], win.sigma + win.d2, edges, win.srcId, mesh->m_pEdge[next_next_edge_idx].m_iVertex[1]);
	}

	if(x < win.b1 + DOUBLE_EPS && !EQUALZERO(sy))								//r1和r2只与e2相交
	{
		if (SolvePlanarEquation(r1, e2, xx, yy))
		{
			if (!EQUALZERO(y)) a = (y-yy)/y*len3;
			else a = (xx-x_project)/(len1-x_project)*len3;
			d1 = sqrt((sx-xx)*(sx-xx) + (sy-yy)*(sy-yy));

			if (SolvePlanarEquation(r2, e2, xx, yy))
			{
				if (!EQUALZERO(y)) b = (y-yy)/y*len3;
				else b = (xx-x_project)/(len1-x_project)*len3;
				d2 = sqrt((sx-xx)*(sx-xx) + (sy-yy)*(sy-yy));
				if(BuildUpWindow(newwindow, next_next_edge_idx, a, b, d1, d2, win.sigma, true, POINTWINDOW, win.srcId, win.pseuSrcId))
					AddWindow(newwindow);
			}
		}
	}
	//else if((x > b2 + DOUBLE_EPS || EQUALZERO(x-b2)))						//r1和r2只与e1相交
	else if(x > win.b2 - DOUBLE_EPS && !EQUALZERO(sy))						//r1和r2只与e1相交
	{
		if (SolvePlanarEquation(r1, e1, xx, yy))
		{
			if (!EQUALZERO(y)) a = yy/y*len2;
			else a = xx/x_project*len2;
			d1 = sqrt((sx-xx)*(sx-xx) + (sy-yy)*(sy-yy));

			if (SolvePlanarEquation(r2, e1, xx, yy))
			{
				if (!EQUALZERO(y)) b = yy/y*len2;
				else b = xx/x_project*len2;
				d2 = sqrt((sx-xx)*(sx-xx) + (sy-yy)*(sy-yy));
				if(BuildUpWindow(newwindow, next_edge_idx, a, b, d1, d2, win.sigma, true, POINTWINDOW, win.srcId, win.pseuSrcId))
					AddWindow(newwindow);
			}
		}
	}
	else if(x  > win.b1 + DOUBLE_EPS && x + DOUBLE_EPS < win.b2)							//r1和r2与e1和e2均相交;此时不需要考虑sy是否为0
	{
		if (SolvePlanarEquation(r1, e1, xx, yy))
		{
			if (!EQUALZERO(y)) a = yy/y*len2;
			else a = xx/x_project*len2;
			d1 = sqrt((sx-xx)*(sx-xx) + (sy-yy)*(sy-yy));
			b = len2;
			d2 = sqrt((sx-x_project)*(sx-x_project) + (sy-y)*(sy-y));
			if(BuildUpWindow(newwindow, next_edge_idx, a, b, d1, d2, win.sigma, true, POINTWINDOW, win.srcId, win.pseuSrcId))
				AddWindow(newwindow);
		}

		if (SolvePlanarEquation(r2, e2, xx, yy))
		{
			if (!EQUALZERO(y)) b = (y-yy)/y*len3;
			else b = (xx-x_project)/(len1-x_project)*len3;
			d2 = sqrt((sx-xx)*(sx-xx) + (sy-yy)*(sy-yy));
			a = 0;
			d1 = sqrt((sx-x_project)*(sx-x_project) + (sy-y)*(sy-y));
			if(BuildUpWindow(newwindow, next_next_edge_idx, a, b, d1, d2, win.sigma, true, POINTWINDOW, win.srcId, win.pseuSrcId))
				AddWindow(newwindow);
		}
	}

	if(x < win.b1 + DOUBLE_EPS) {
		double curUpdateDist = win.sigma + win.d1 + sqrt((win.b1-x_project)*(win.b1-x_project)+y*y)/*sqrt(sx*sx+sy*sy) + len2*/;
		if (curUpdateDist < vertexDist[vert[1]]) {
			vertexDist[vert[1]] = curUpdateDist;
			VertOccupyInfo voi; 
			voi.isVert = true; voi.elementIndex = vert[0]; voi.pos = 0;
			vertOccupyied[vert[1]] = voi;
		}
	}
	else if(x > win.b2 - DOUBLE_EPS) {
		double curUpdateDist = win.sigma + win.d2 + sqrt((win.b2-x_project)*(win.b2-x_project)+y*y)/*sqrt((sx-len1)*(sx-len1)+sy*sy) + len3*/;
		if (curUpdateDist < vertexDist[vert[1]]) {
			vertexDist[vert[1]] = curUpdateDist;
			VertOccupyInfo voi; 
			voi.isVert = true; voi.elementIndex = vert[2]; voi.pos = len1;
			vertOccupyied[vert[1]] = voi;
		}
	}
	else {
		double curUpdateDist = win.sigma + sqrt((sx-x_project)*(sx-x_project)+(sy-y)*(sy-y));
		if (curUpdateDist < vertexDist[vert[1]]) {
			vertexDist[vert[1]] = curUpdateDist;
			VertOccupyInfo voi; 
			voi.isVert = false; voi.elementIndex = win.Edge_idx; voi.pos = x;
			vertOccupyied[vert[1]] = voi;
		}
	}
}

void ExactGeodesic::PropogateLineWindow(const Window& win)
{
	UINT edge_idx = mesh->m_pEdge[win.Edge_idx].m_iTwinEdge;
	UINT next_edge_idx = mesh->m_pEdge[edge_idx].m_iNextEdge;
	UINT next_next_edge_idx = mesh->m_pEdge[next_edge_idx].m_iNextEdge;

	double len1 = mesh->m_pEdge[edge_idx].m_length;
	double len2 = mesh->m_pEdge[next_edge_idx].m_length;
	double len3 = mesh->m_pEdge[next_next_edge_idx].m_length;

	double x_project = (len1*len1 + len2*len2 - len3*len3) / (2*len1);
	double y = -sqrt(fabs(len2*len2 - x_project*x_project));

	Line e1, e2;
	e1.A = y;
	e1.B = -x_project;
	e1.C = 0;

	e2.A = y;
	e2.B = len1 - x_project;
	e2.C = -len1*y;

	Line l1, l2;

	double sina = fabs(win.d1-win.d2)/(win.b2-win.b1);
	double cosa = sqrt(1-sina*sina);

	if(win.d1 > win.d2)
	{
		l1.A = cosa; l1.B = -sina; l1.C = -cosa*win.b1;
	}
	else
	{
		l1.A = -cosa; l1.B = -sina; l1.C = cosa*win.b1;
	}
	l2.A = l1.A; l2.B = l1.B; l2.C = -l2.A*win.b2;

	double x1, y1, x2, y2 = 1e30;
	double b1, b2, d1, d2;

	SolvePlanarEquation(l1, e1, x1, y1);
	SolvePlanarEquation(l2, e1, x2, y2);

	if(x1 < x_project)
	{
		d1 = (win.d2-win.d1)/(win.b2-win.b1)*win.b1+(win.b2*win.d1-win.b1*win.d2)/(win.b2-win.b1)+sqrt((win.b1-x1)*(win.b1-x1)+y1*y1);
		double px = 0;
		if(x2 > x_project)
		{
			x2 = x_project;
			y2 = y;
			Line tmpl;
			tmpl.A = l1.A; tmpl.B = l1.B; tmpl.C = -l1.A*x_project-l1.B*y;
			px = -tmpl.C/tmpl.A;
		}
		else
		{
			px = win.b2;
		}
		d2 = (win.d2-win.d1)/(win.b2-win.b1)*px+(win.b2*win.d1-win.b1*win.d2)/(win.b2-win.b1)+sqrt((px-x2)*(px-x2)+y2*y2);

		b1 = y1/y*len2;
		b2 = y2/y*len2;
		if(!EQUALZERO(b1-b2))
		{
			Window newwindow;
			if(BuildUpWindow(newwindow, next_edge_idx, b1, b2, d1, d2, win.sigma, true, LINEWINDOW, win.srcId, win.pseuSrcId))
				AddWindow(newwindow);
		}
	}

	SolvePlanarEquation(l1, e2, x1, y1);
	SolvePlanarEquation(l2, e2, x2, y2);
	if(x2 > x_project)
	{
		d2 = (win.d2-win.d1)/(win.b2-win.b1)*win.b2+(win.b2*win.d1-win.b1*win.d2)/(win.b2-win.b1)+sqrt((win.b2-x2)*(win.b2-x2)+y2*y2);
		double px = 0;
		if(x1 < x_project)
		{
			x1 = x_project;
			y1 = y;
			Line tmpl;
			tmpl.A = l1.A; tmpl.B = l1.B; tmpl.C = -l1.A*x_project-l1.B*y;
			px = -tmpl.C/tmpl.A;
		}
		else
		{
			px = win.b1;
		}
		d1 = (win.d2-win.d1)/(win.b2-win.b1)*px+(win.b2*win.d1-win.b1*win.d2)/(win.b2-win.b1)+sqrt((px-x1)*(px-x1)+y1*y1);

		b1 = (y-y1)/y*len3;
		b2 = (y-y2)/y*len3;
		if(!EQUALZERO(b1-b2))
		{
			Window newwindow;
			if(BuildUpWindow(newwindow, next_next_edge_idx, b1, b2, d1, d2, win.sigma, true, LINEWINDOW, win.srcId, win.pseuSrcId))
				AddWindow(newwindow);
		}
	}
}

void ExactGeodesic::AddPseudoSource(UINT pseudo_src_idx, double sigma)
{
	//从伪源点扩散到以其为顶点的所有的面上
	Window newWindow;
	int faceNum = mesh->m_pVertex[pseudo_src_idx].m_nValence;

	UINT Edge_idx, next_Edge_idx, nextnext_Edge_idx;
	double Edgelength, nextEdgelength, nextnextEdgelength;

	if (sigma == 0) 
	{
		vertOccupyied[pseudo_src_idx].elementIndex = -1;
		vertexDist[pseudo_src_idx] = sigma;
	}

	for(int i = 0;i < faceNum;++i)
	{
		Edge_idx = mesh->m_pVertex[pseudo_src_idx].m_piEdge[i];
		next_Edge_idx = mesh->m_pEdge[Edge_idx].m_iNextEdge;
		nextnext_Edge_idx = mesh->m_pEdge[next_Edge_idx].m_iNextEdge;

		Edgelength = mesh->m_pEdge[Edge_idx].m_length;
		nextEdgelength = mesh->m_pEdge[next_Edge_idx].m_length;
		nextnextEdgelength = mesh->m_pEdge[nextnext_Edge_idx].m_length;

		//构造在边Edge_idx上的新窗口
		if(BuildUpWindow(newWindow, Edge_idx, 0, Edgelength, 0, Edgelength, sigma, true, POINTWINDOW, pseudo_src_idx, pseudo_src_idx))
			AddWindow(newWindow);

		//构造在边next_Edge_idx上的新窗口
		if(BuildUpWindow(newWindow, next_Edge_idx, 0, nextEdgelength, Edgelength, nextnextEdgelength, sigma, true, POINTWINDOW, pseudo_src_idx, pseudo_src_idx))
			AddWindow(newWindow);

		//构造在边nextnext_Edge_idx上的新窗口
		if(BuildUpWindow(newWindow, nextnext_Edge_idx,0, nextnextEdgelength, nextnextEdgelength, 0, sigma, true, POINTWINDOW, pseudo_src_idx, pseudo_src_idx))
			AddWindow(newWindow);
	}
}

void ExactGeodesic::AddPseudoSource(UINT face_idx, Vector3D position, double sigma, UINT srcId, UINT pseuSrcId)
{
	Window newWindow;
	UINT cur_edge, v0, v1;
	double d0, d1;

	for(int i = 0;i < mesh->m_pFace[face_idx].m_nType;++i)
	{
		cur_edge = mesh->m_pFace[face_idx].m_piEdge[i];
		v0 = mesh->m_pEdge[cur_edge].m_iVertex[0];
		v1 = mesh->m_pEdge[cur_edge].m_iVertex[1];

		d0 = (position-mesh->m_pVertex[v0].m_vPosition).length();
		d1 = (position-mesh->m_pVertex[v1].m_vPosition).length();

		if(BuildUpWindow(newWindow, cur_edge, 0, mesh->m_pEdge[cur_edge].m_length, d0, d1, sigma, true, POINTWINDOW, srcId, pseuSrcId))
			AddWindow(newWindow);
		vertOccupyied[v0].elementIndex = -1;
	}
}

void ExactGeodesic::AddPseudoSource(UINT pseudo_src_idx, double sigma, UINT* edges, UINT srcId, UINT pseuSrcId)
{
	Window newWindow;
	Vector3D position = mesh->m_pVertex[pseudo_src_idx].m_vPosition;

	UINT cur_edge, v0, v1;
	double d0, d1;

	for(int i = 0;i < 2;++i)
	{
		cur_edge = edges[i];
		v0 = mesh->m_pEdge[cur_edge].m_iVertex[0];
		v1 = mesh->m_pEdge[cur_edge].m_iVertex[1];

		d0 = pseudo_src_idx == v0 ? 0 : (position-mesh->m_pVertex[v0].m_vPosition).length();
		d1 = pseudo_src_idx == v1 ? 0 : (position-mesh->m_pVertex[v1].m_vPosition).length();

		if(BuildUpWindow(newWindow, cur_edge, 0, mesh->m_pEdge[cur_edge].m_length, d0, d1, sigma, true, POINTWINDOW, srcId, pseuSrcId))
			AddWindow(newWindow);
	}
}

void ExactGeodesic::AddPseudoSource(LineSource ls)
{
	for(int i = 0;i < mesh->m_pFace[ls.Face_idx].m_nType;++i)
	{
		UINT edge_idx = mesh->m_pFace[ls.Face_idx].m_piEdge[i];
		Vector3D v0 = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[0]].m_vPosition;
		Vector3D v1 = mesh->m_pVertex[mesh->m_pEdge[edge_idx].m_iVertex[1]].m_vPosition;
		Vector3D vx = v1 - v0;
		vx.normalize();

		for(int j = 0;j < (int)ls.points_on_plat.size();++j)
		{
			Vector3D tmppoint = ls.points_on_plat[j] - v0;
			double ex1 = tmppoint * vx;										//(ex1,ey1)和(ex2,ey2)为第j和第j+1个顶点在平面上的坐标
			double ey1 = (tmppoint ^ vx).length();

			tmppoint = ls.points_on_plat[(j+1)%ls.points_on_plat.size()] - v0;
			double ex2 = tmppoint * vx;
			double ey2 = (tmppoint ^ vx).length();

			if(ex1 > ex2)
			{
				double tmp = ex1;
				ex1 = ex2;
				ex2 = tmp;

				tmp = ey1;
				ey1 = ey2;
				ey2 = tmp;
			}

			if(!EQUALZERO(ex1-ex2))												//只有这样该段线源在该边上才会有可能有窗口
			{
				Line l;																		//线源所在直线
				l.A = ey2-ey1; l.B = ex1-ex2; l.C = ex2*ey1-ex1*ey2;

				Line l1;																		//经过x1y1垂直于l的直线
				l1.A = l.B; l1.B = -l.A; l1.C = l.A*ey1-l.B*ex1;
				double b1 = -l1.C/l1.A;

				Line l2;																		//经过x2y2垂直琎l的直线
				l2.A = l.B; l2.B = -l.A; l2.C = l.A*ey2-l.B*ex2;
				double b2 = -l2.C/l2.A;

				double d1, d2 = 0;

				if(b1 < 0 && b2 < 0 || b1 > mesh->m_pEdge[edge_idx].m_length 
					&& b2 > mesh->m_pEdge[edge_idx].m_length)			//没有在该边上投射出线窗口
					continue;
				else if(b1 < 0 && b2 > 0 && b2 < mesh->m_pEdge[edge_idx].m_length)
				{
					b1 = 0;
					Line tmpl;
					tmpl.A = l.B; tmpl.B = -l.A;tmpl.C = 0;
					double tmpx, tmpy;
					SolvePlanarEquation(tmpl, l, tmpx, tmpy);
					d1 = sqrt(tmpx*tmpx+tmpy*tmpy);
					d2 = sqrt((b2-ex2)*(b2-ex2)+ey2*ey2);
				}
				else if(b1 > 0 && b1 < mesh->m_pEdge[edge_idx].m_length && b2 > mesh->m_pEdge[edge_idx].m_length)
				{
					b2 = mesh->m_pEdge[edge_idx].m_length;
					Line tmpl;
					tmpl.A = l.B; tmpl.B = -l.A; tmpl.C = -l.B*b2;
					double tmpx, tmpy;
					SolvePlanarEquation(tmpl, l, tmpx, tmpy);
					d1 = sqrt((b1-ex1)*(b1-ex1)+ey1*ey1);
					d2 = sqrt((tmpx-b2)*(tmpx-b2)+tmpy*tmpy);
				}
				else if(b1 > 0 && b1 < mesh->m_pEdge[edge_idx].m_length && b2 > 0 && b2 < mesh->m_pEdge[edge_idx].m_length)
				{
					d1 = sqrt((b1-ex1)*(b1-ex1)+ey1*ey1);
					d2 = sqrt((b2-ex2)*(b2-ex2)+ey2*ey2);
				}

				Window window;
				UINT srcId = 0;
				if(BuildUpWindow(window, edge_idx, b1, b2, d1, d2, 0, true, LINEWINDOW, ls.srcId, ls.srcId))
					AddWindow(window);
			}
		}
	}
}

bool ExactGeodesic::BuildUpWindow(Window& win, UINT Edge_idx, double b1, double b2, double d1, double d2, double sigma, bool tao, Window_Type winType, UINT srcId, UINT pseuSrcId)
{
	if(b1 > b2 || EQUALZERO(b1-b2))						//左边界不小于右边界的窗口非法
		return false;
	if(b1 < -DOUBLE_EPS || b2 < -DOUBLE_EPS)
		return false;
	win.Edge_idx = Edge_idx;
	win.b1 = b1;
	win.b2 = b2;
	win.d1 = d1;
	win.d2 = d2;
	win.sigma = sigma;
	win.tao = tao;
	win.winType = winType;
	win.idx_in_pq = -1;
	win.srcId = srcId;
	win.pseuSrcId = pseuSrcId;

	return true;
}

void ExactGeodesic::SolveEquation(const Window& w0, const Window& w1, double& p1, double& p2)
{
	if(w0.winType == POINTWINDOW && w1.winType == POINTWINDOW)
		SolvePPEquation(w0, w1, p1, p2);
	else if(w0.winType == POINTWINDOW && w1.winType == LINEWINDOW ||
		w0.winType == LINEWINDOW && w1.winType == POINTWINDOW)
		SolvePLEquation(w0, w1, p1, p2);
	else if(w0.winType == LINEWINDOW && w1.winType == LINEWINDOW)
		SolveLLEquation(w0, w1, p1, p2);
	else
		p1 = p2 = 1e30;
}

void ExactGeodesic::SolvePPEquation(const Window& w0, const Window& w1, double& p1, double& p2)
{
	double s0 = w0.b2 - w0.b1;
	double s1 = w1.b2 - w1.b1;
	double s0x = (w0.d1*w0.d1 + s0*s0 - w0.d2*w0.d2) / (2*s0) + w0.b1;
	double s1x = (w1.d1*w1.d1 + s1*s1 - w1.d2*w1.d2) / (2*s1) + w1.b1;

	double s0y = sqrt(fabs(w0.d1*w0.d1 - (s0x-w0.b1)*(s0x-w0.b1)));
	double s1y = sqrt(fabs(w1.d1*w1.d1 - (s1x-w1.b1)*(s1x-w1.b1)));

	double alpha = s1x - s0x;
	double beta = w1.sigma - w0.sigma;
	double gama = (s0x*s0x+s0y*s0y) - (s1x*s1x+s1y*s1y) - beta*beta;

	double A = alpha*alpha - beta*beta;
	double B = gama*alpha + 2*s1x*beta*beta;
	double C = 0.25*gama*gama - (s1x*s1x+s1y*s1y)*beta*beta;

	if(EQUALZERO(A))
	{
		if(!EQUALZERO(B))
		{
			p1 = -C/B;
			p2 = 1e30;
		}
		else if(EQUALZERO(C))
		{
			p1 = 1e30;
			p2 = 1e30;
		}
		else
			return;
	}
	else
	{
		double delta = B*B - 4*A*C;
		if(delta < -DOUBLE_EPS)							//delta < 0
			return;
		else if(delta < DOUBLE_EPS)						//delta == 0
		{
			p1 = -B/(2*A);
			p2 = 1e30;
		}
		else															//delta > 0
		{
			p1 = (-B + sqrt(delta)) / (2*A);
			p2 = (-B - sqrt(delta)) / (2*A);
		}
	}

	if(p1 > p2)
	{
		double tmp = p1;
		p1 = p2;
		p2 = tmp;
	}
	return;
}

void ExactGeodesic::SolvePLEquation(const Window& w0, const Window& w1, double& p1, double& p2)
{
	double m, n, sx, sy;
	if(w0.winType == LINEWINDOW)
	{
		m = (w0.d2 - w0.d1)/(w0.b2-w0.b1);
		n = (w0.b2*w0.d1-w0.b1*w0.d2)/(w0.b2-w0.b1);
		sx = (w1.d1*w1.d1+(w1.b2-w1.b1)*(w1.b2-w1.b1)-w1.d2*w1.d2)/(2*(w1.b2-w1.b1))+w1.b1;
		sy = sqrt(w1.d1*w1.d1 - (sx-w1.b1)*(sx-w1.b1));
	}
	else
	{
		m = (w1.d2 - w1.d1)/(w1.b2-w1.b1);
		n = (w1.b2*w1.d1-w1.b1*w1.d2)/(w1.b2-w1.b1);
		sx = (w0.d1*w0.d1+(w0.b2-w0.b1)*(w0.b2-w0.b1)-w0.d2*w0.d2)/(2*(w0.b2-w0.b1))+w0.b1;
		sy = sqrt(w0.d1*w0.d1 - (sx-w0.b1)*(sx-w0.b1));
	}

	double A = 1-m*m;
	double B = -2*(sx+m*n);
	double C = sx*sx+sy*sy-n*n;

	if(EQUALZERO(A))
	{
		if(!EQUALZERO(B))
		{
			p1 = -C/B;
			p2 = 1e30;
		}
		else if(EQUALZERO(C))
		{
			p1 = -1e30;
			p2 = -1e30;
		}
		else
		{
			p1 = 1e30;
			p2 = 1e30;
		}
	}
	else
	{
		double delta = B*B - 4*A*C;
		if(delta < -DOUBLE_EPS)							//delta < 0
		{
			p1 = 1e30;
			p2 = 1e30;
		}
		else if(delta < DOUBLE_EPS)						//delta == 0
		{
			p1 = -B/(2*A);
			p2 = 1e30;
		}
		else															//delta > 0
		{
			p1 = (-B + sqrt(delta)) / (2*A);
			p2 = (-B - sqrt(delta)) / (2*A);
		}
	}

	if(p1 > p2)
	{
		double tmp = p1;
		p1 = p2;
		p2 = tmp;
	}
}

void ExactGeodesic::SolveLLEquation(const Window& w0, const Window& w1, double& p1, double& p2)
{
	double m0 = (w0.d2 - w0.d1)/(w0.b2 - w0.b1);
	double n0 = (w0.b2*w0.d1-w0.b1*w0.d2)/(w0.b2-w0.b1);
	double m1 = (w1.d2 - w1.d1)/(w1.b2 - w1.b1);
	double n1 = (w1.b2*w1.d1-w1.b1*w1.d2)/(w1.b2-w1.b1);

	if(EQUALZERO(m1-m0))
	{
		p1 = p2 = 1e30;
	}
	else
	{
		p1 = p2 = (n1-n0)/(m0-m1);
	}
}

bool ExactGeodesic::SolvePlanarEquation(const Line& l1, const Line& l2, 	double& x, double& y)
{
	//两条直线平行，返回false
	if(EQUALZERO(l1.A*l2.B - l2.A*l1.B))
	{
		x = 1e30;
		y = 1e30;
		return false;
	}

	x = (l1.B*l2.C-l2.B*l1.C)/(l1.A*l2.B-l2.A*l1.B);
	y = (l2.A*l1.C-l1.A*l2.C)/(l1.A*l2.B-l2.A*l1.B);

	return true;
}

bool ExactGeodesic::CloserToWindow(const Window& w0, const Window& w1, double p, double& pathlen0, double& pathlen1)
{
	if(w0.winType == POINTWINDOW)
	{
		double s0 = w0.b2 - w0.b1;
		double s0x = (w0.d1*w0.d1 + s0*s0 - w0.d2*w0.d2) / (2*s0) + w0.b1;
		double s0y = sqrt(fabs(w0.d1*w0.d1 - (s0x-w0.b1)*(s0x-w0.b1)));

		pathlen0 = sqrt((s0x-p)*(s0x-p) + s0y*s0y) + w0.sigma;
	}
	else
	{
		pathlen0 = p*(w0.d2-w0.d1)/(w0.b2-w0.b1)+(w0.b2*w0.d1-w0.b1*w0.d2)/(w0.b2-w0.b1);
	}
	if(w1.winType == POINTWINDOW)
	{
		double s1 = w1.b2 - w1.b1;
		double s1x = (w1.d1*w1.d1 + s1*s1 - w1.d2*w1.d2) / (2*s1) + w1.b1;
		double s1y = sqrt(fabs(w1.d1*w1.d1 - (s1x-w1.b1)*(s1x-w1.b1)));

		pathlen1 = sqrt((s1x-p)*(s1x-p) + s1y*s1y) + w1.sigma;
	}
	else
	{
		pathlen1 = p*(w1.d2-w1.d1)/(w1.b2-w1.b1)+(w1.b2*w1.d1-w1.b1*w1.d2)/(w1.b2-w1.b1);
	}
	return pathlen0 + DOUBLE_EPS < pathlen1;
}

bool ExactGeodesic::CloserToWindow(const Window& w0, const Window& w1, double p)
{
	double pathlen0, pathlen1;
	if(w0.winType == POINTWINDOW)
	{
		double s0 = w0.b2 - w0.b1;
		double s0x = (w0.d1*w0.d1 + s0*s0 - w0.d2*w0.d2) / (2*s0) + w0.b1;
		double s0y = sqrt(fabs(w0.d1*w0.d1 - (s0x-w0.b1)*(s0x-w0.b1)));

		pathlen0 = sqrt((s0x-p)*(s0x-p) + s0y*s0y) + w0.sigma;
	}
	else
	{
		pathlen0 = p*(w0.d2-w0.d1)/(w0.b2-w0.b1)+(w0.b2*w0.d1-w0.b1*w0.d2)/(w0.b2-w0.b1);
	}
	if(w1.winType == POINTWINDOW)
	{
		double s1 = w1.b2 - w1.b1;
		double s1x = (w1.d1*w1.d1 + s1*s1 - w1.d2*w1.d2) / (2*s1) + w1.b1;
		double s1y = sqrt(fabs(w1.d1*w1.d1 - (s1x-w1.b1)*(s1x-w1.b1)));

		pathlen1 = sqrt((s1x-p)*(s1x-p) + s1y*s1y) + w1.sigma;
	}
	else
	{
		pathlen1 = p*(w1.d2-w1.d1)/(w1.b2-w1.b1)+(w1.b2*w1.d1-w1.b1*w1.d2)/(w1.b2-w1.b1);
	}
	return pathlen0 + DOUBLE_EPS < pathlen1;
}

void ExactGeodesic::AddWindow(Window& win)
{
	UINT edgeIdx = win.Edge_idx;
	if (edgeIdx == -1) return;
	if(EQUALZERO(win.b1 - win.b2)) return;

	Window* newWindow = &win;
	WindowList::iterator it = windowlist[win.Edge_idx].begin();
	bool shouldIterDecr = false;

	Window intersectWindow;
	Window updateWindow;
	WindowList::iterator existWindow, ptr, pos;

	double left, right, p1, p2, p;
	bool chooselabel;

	double bounds[6];
	bool belongExist[5];
	int subIdx;

	double path0, path1;

	int winCounter = 0;
	int label = 0;

	for(; it != windowlist[edgeIdx].end();)
	{
		existWindow = it++;
		if(existWindow->b1 > newWindow->b2 + DOUBLE_EPS || EQUALZERO(existWindow->b1-newWindow->b2))					//如果当前最左已存在窗口的b1都已大于新窗口newWindow的b2，则后面也不会有相交了
		{
			shouldIterDecr = true;
			break;
		}

		if(newWindow->b1 > existWindow->b2 + DOUBLE_EPS || EQUALZERO(newWindow->b1-existWindow->b2))					//当前existWindow还未与newWindow相交，进行下一个existWindow的判断
			continue;

		if(EQUALZERO(newWindow->b1 - newWindow->b2))					//如果本轮的newWindow已经是点窗口了，则退出循环
			break;

		bool existWindowIsInPQ = (existWindow->idx_in_pq > 0 && existWindow->idx_in_pq < (int)windowPriorityQueue.size());				//当前windowlist上的窗口是否存在于PQ中

		updateWindow.b1 = 0;													//将其初始化为点窗口
		updateWindow.b2 = 0;													//而这个条件则可保证如果没有可更新的newWindow，下一轮循环就退出了

		/*
		以下处理不相交的部分作为新窗口插入
		*/
		subIdx = 0;
		//existWindow左侧超出的部分被截为新窗口
		if(existWindow->b1 + DOUBLE_EPS < newWindow->b1)
		{
			bounds[subIdx] = existWindow->b1;
			bounds[subIdx+1] = newWindow->b1;
			belongExist[subIdx] = true;
			++subIdx;
		}
		//或者newWindow左侧超出的部分被截为新窗口
		else if(newWindow->b1 + DOUBLE_EPS < existWindow->b1)
		{
			bounds[subIdx] = newWindow->b1;
			bounds[subIdx+1] = existWindow->b1;
			belongExist[subIdx] = false;
			++subIdx;
		}
		else bounds[subIdx] = existWindow->b1;

		/*
		以下处理相交部分
		相交的部分要处理不同窗口的情况
		*/
		left = existWindow->b1 > newWindow->b1 ? existWindow->b1 : newWindow->b1;
		right = existWindow->b2 < newWindow->b2 ? existWindow->b2 : newWindow->b2;

		p1 = 1e30, p2 = 1e30;
		double p[4]; int psize = 0;
		p[psize] = left;
		++psize;

		SolveEquation(*existWindow, *newWindow, p1, p2);

		if (p1 > left + DOUBLE_EPS && p1 + DOUBLE_EPS < right)
		{
			p[psize] = p1;
			++psize;
		}
		if (p2 > left + DOUBLE_EPS && p2 + DOUBLE_EPS < right)
		{
			p[psize] = p2;
			++psize;
		}
		p[psize] = right;
		++psize;

		for (int i = 0;i < psize-1;++i)
		{
			bounds[subIdx+1] = p[i+1];
			CloserToWindow(*existWindow, *newWindow, (p[i]+p[i+1])/2, path0, path1);
			if (path0 < path1 + DOUBLE_EPS) belongExist[subIdx] = true;
			else if (path1 < path0 + DOUBLE_EPS) belongExist[subIdx] = false;
			++subIdx;
		}

		//existWindow右侧超出的部分被截为新窗口
		if (existWindow->b2 > newWindow->b2 + DOUBLE_EPS)
		{
			bounds[subIdx+1] = existWindow->b2;
			belongExist[subIdx] = true;
			++subIdx;
		}
		//或者newWindow右侧超出的部分被截为新窗口，则切出来留作下次的newWindow
		else if (newWindow->b2 > existWindow->b2 + DOUBLE_EPS)
		{
			bounds[subIdx+1] = newWindow->b2;
			belongExist[subIdx] = false;
			++subIdx;
		}

		//Merge all adjacent windows
		for (int i = 0;i < subIdx-1;)
		{
			if (belongExist[i] == belongExist[i+1])
			{
				bounds[i+1] = bounds[i+2];
				
				for (int j = i+1; j < subIdx-1;++j)
				{
					bounds[j+1] = bounds[j+2];
					belongExist[j] = belongExist[j+1];
				}
				--subIdx;
			}
			else ++i;
		}

		winCounter += subIdx;

		bool isExistWindowNotCut = false;
		for (int i = 0;i < subIdx;++i)
		{
			if (EQUALZERO(bounds[i] - existWindow->b1) && EQUALZERO(bounds[i+1]-existWindow->b2)
				&& belongExist[i]){
					isExistWindowNotCut = true;
					continue;
			}
			//if ((bounds[i+1]-bounds[i]) / mesh->m_pEdge[edgeIdx].m_length < 0.01) continue;
			if (belongExist[i]) CutWindow(*existWindow, bounds[i], bounds[i+1], intersectWindow);
			else CutWindow(*newWindow, bounds[i], bounds[i+1], intersectWindow);
			/*
			CloserToWindow(*existWindow, *newWindow, (bounds[i]+bounds[i+1])/2, path0, path1);
			if (!(path0 < path1 + DOUBLE_EPS) && !(path1 < path0 + DOUBLE_EPS))
				intersectWindow.srcId = existWindow->srcId | newWindow->srcId;
				*/
			if (i != subIdx-1 || belongExist[i])
			{
				ptr = existWindow;
				if (isExistWindowNotCut)
					windowlist[edgeIdx].insert(++ptr, intersectWindow);
				else
					windowlist[edgeIdx].insert(ptr, intersectWindow);

				if (!belongExist[i] || existWindowIsInPQ)
				{
					windowPriorityQueue.push(intersectWindow, CalcMinDist(intersectWindow), --ptr);
				}
			}
			else
				updateWindow = intersectWindow;
		}

		if (!isExistWindowNotCut)
		{
			if(existWindowIsInPQ)
				windowPriorityQueue.remove(existWindow->idx_in_pq);
			windowlist[existWindow->Edge_idx].erase(existWindow);
		}

		*newWindow = updateWindow;
	}

	if(newWindow->b1 + DOUBLE_EPS	 < newWindow->b2)													//如果是点窗口，则不将其插入
	{
		if(shouldIterDecr) --it;
		windowlist[newWindow->Edge_idx].insert(it, *newWindow);
		windowPriorityQueue.push(*newWindow, CalcMinDist(*newWindow), --it);
		++winCounter;
	}

	//CombineSameWindow(edgeIdx, winCounter);

	winNum += winCounter;

	//delete newWindow;
}

void ExactGeodesic::CutWindow(const Window& win, double leftb, double rightb, Window& subwin)
{
	subwin.Edge_idx = win.Edge_idx;
	subwin.b1 = leftb;
	subwin.b2 = rightb;
	subwin.sigma = win.sigma;
	subwin.tao = win.tao;
	subwin.winType = win.winType;
	subwin.idx_in_pq = -1;
	subwin.srcId = win.srcId;
	subwin.pseuSrcId = win.pseuSrcId;

	if(win.winType == POINTWINDOW)
	{
		subwin.d1 = sqrt(fabs((win.d1*win.d1*(win.b2-leftb)+win.d2*win.d2*(leftb-win.b1)-(leftb-win.b1)*(win.b2-leftb)*(win.b2-win.b1))/(win.b2-win.b1)));
		subwin.d2 = sqrt(fabs((win.d1*win.d1*(win.b2-rightb)+win.d2*win.d2*(rightb-win.b1)-(rightb-win.b1)*(win.b2-rightb)*(win.b2-win.b1))/(win.b2-win.b1)));
	}
	else
	{
		subwin.d1 = (win.d2-win.d1)/(win.b2-win.b1)*leftb+(win.b2*win.d1-win.b1*win.d2)/(win.b2-win.b1);
		subwin.d2 = (win.d2-win.d1)/(win.b2-win.b1)*rightb+(win.b2*win.d1-win.b1*win.d2)/(win.b2-win.b1);
	}
}

double ExactGeodesic::CalcMinDist(const Window& win)
{
	if(win.winType == POINTWINDOW)
	{
		double s = win.b2 - win.b1;
		if(win.d1*win.d1 + s*s - win.d2*win.d2 <= 0)
			return win.d1 + win.sigma;
		else if(win.d2*win.d2 + s*s - win.d1*win.d1 <= 0)
			return win.d2 + win.sigma;
		else
		{
			double x = (win.d1*win.d1 + s*s - win.d2*win.d2)/(2*s);
			return sqrt(fabs(win.d1*win.d1 - x*x)) + win.sigma;
		}
	}
	else if(win.winType == LINEWINDOW)												//线窗口的最小距离必然是b1或b2之一
	{
		if(win.d1 < win.d2) return win.d1;
		else return win.d2;
	}
	return 0;
}

bool ExactGeodesic::StopConditionSatisfied(std::vector<SurfacePoint> *stopPoints)
{
	if (!stopPoints || stopPoints->empty()) return false;
	double curQueueTop = windowPriorityQueue.topKey();
	UINT srcId;

	for (UINT i = 0; i < stopPoints->size(); ++i)
	{
		unsigned faceIdx = (*stopPoints)[i].faceIdx;
		Vector3D pos = (*stopPoints)[i].pos;

		double curMinDist = 1e30;
		for (UINT j = 0; j < mesh->m_pFace[faceIdx].m_nType; ++j)
		{
			UINT curVert = mesh->m_pFace[faceIdx].m_piVertex[j];
			double dist = CalcGeodesicDistToVertex(curVert, srcId) + (pos - mesh->m_pVertex[curVert].m_vPosition).length();
			curMinDist = curMinDist < dist ? curMinDist : dist;
		}

		if (curMinDist > curQueueTop) return false;
	}

	return true;
}

bool ExactGeodesic::toLeft(Vector2D& v0, Vector2D& v1, Vector2D& v2)
{
	return ((v2 - v0) ^ (v1 - v0)) < DOUBLE_EPS;
}

bool ExactGeodesic::InTriangle(Vector2D& v0, Vector2D& v1, Vector2D& v2, Vector2D& v3)
{
	bool res1 = toLeft(v0, v1, v2);
	bool res2 = toLeft(v0, v2, v3);
	bool res3 = toLeft(v0, v3, v1);
	return (res1 && res2 && res3) || (!res1 && !res2 && !res3);
}

void ExactGeodesic::ReverseWinToTwinEdge(Window &win)
{
	unsigned twinEdge = mesh->m_pEdge[win.Edge_idx].m_iTwinEdge;
	double b1 = mesh->m_pEdge[twinEdge].m_length - win.b2;
	double b2 = mesh->m_pEdge[twinEdge].m_length - win.b1;
	win.b1 = b1; win.b2 = b2;
	swap(win.d1, win.d2);
	win.Edge_idx = twinEdge;
	win.tao = !win.tao;
}