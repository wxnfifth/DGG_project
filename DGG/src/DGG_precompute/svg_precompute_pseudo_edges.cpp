#include "svg_precompute_pseudo_edges.h"
#include "wxn\wxn_dijstra.h"
#include "wxn\wxn_geometry.h"
#include "MMP\geodesic_algorithm_vg_mmp.h"
#include "YXMetric\YXPathTracer.h"
#include "wxn\wxn_path_helper.h"
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include <windows.h>

#define WXN_PRINT_DEBUG_LINE(msg) fprintf(stderr, ""msg" at %s, line %d.\n", __FILE__, __LINE__)


void initial_heads_and_parts(vector<BodyHeadOfSVG>& total_heads,
	vector<vector<BodyPartOfSVGWithAngle>>& total_parts,
	CRichModel& model, SparseGraph<float>* s_graph)
{
	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		auto& geo_dises = s_graph->graphNeighborDis(i); //i点的到各个邻居的geodesic distance 
		auto& angles = s_graph->graphNeighborAngle(i); //i点到各个邻居的angle 
		total_parts[i].reserve(geo_dises.size());
		for (int j = 0; j < geo_dises.size(); ++j) {
			total_parts[i].push_back(BodyPartOfSVGWithAngle(s_graph->graphNeighbor(i)[j], geo_dises[j], angles[j], 0, 0));
		}
		total_heads[i].source_index = i;
		total_heads[i].neighbor_num = geo_dises.size();


		if (i == 1700) {
			CylinderPath path(0.0005);
			path.addGeodesicPaths(model, i, s_graph->graphNeighbor(i));
			path.write_to_file("paths_without_pseudo_edges.obj");
			printBallToObj(vector < CPoint3D > {model.Vert(i)}, "paths_source.obj", 0.003);
		}

	}
}


void find_pseudo_point_using_path_tracer(const CRichModel& model, const int source_vert,
	const int direction_edge, const double dis,
	const double angle, YXPathTracer& path_tracer,
	CPoint3D& p, int& face_idx)
{
	//double yx_angle = model.AngleSum(source_vert) - angle;
	double yx_angle = angle;
	const auto& e = model.Edge(direction_edge);
	int v1 = e.indexOfLeftVert;
	int v2 = e.indexOfRightVert;
	int e_id_yx = path_tracer.metric.getEdgeFrom2Verts(v1, v2);
	YXPath tmp_path;
	path_tracer.computeSinglePathArbitraryStart(e_id_yx, yx_angle, dis, tmp_path, source_vert);
	p = path_tracer.convertToYXPoint3D(tmp_path.back()).toCPoint3D();

	{
		YXPathPoint tmp_path_back = tmp_path.back();
		int v0 = path_tracer.metric.Edge(tmp_path_back.edgeId).v1;
		int v1 = path_tracer.metric.Edge(tmp_path_back.edgeId).v2;
		int v2 = path_tracer.metric.Edge(path_tracer.metric.getNextEdgeAroundFace(tmp_path_back.edgeId)).v2;
		//printf("v0 %d v1 %d v2 %d\n", v0, v1, v2);
		face_idx = model.GetFaceIndexFromTreeVertices(v0, v1, v2);
		//printf("f %d v0 %d v1 %d v2 %d\n", face_idx, model.Face(face_idx)[0], model.Face(face_idx)[1], model.Face(face_idx)[2]);
		//face_idx = path_tracer.getFaceIndex(tmp_path.back());
	}
	if (false) {
		static int a = 0;
		CylinderPath path(0.001);
		vector<CPoint3D> pts;
		for (int i = 0; i < tmp_path.size(); ++i) {
			auto pts_c3d = path_tracer.convertToYXPoint3D(tmp_path[i]).toCPoint3D();
			pts.push_back(pts_c3d);
			printf("pts a %lf b %lf\n", tmp_path[i].a, tmp_path[i].b);
			printf("pts xyz %lf %lf %lf\n", pts_c3d.x, pts_c3d.y, pts_c3d.z);
		}
		path.addLines(pts);
		path.write_to_file("path_tracer_" + to_string(a++) + ".obj");
		printf("tmp_path.size() %d \n", tmp_path.size() );
		printf("\n");
	}
	//const CRichModel& model = *refined_model_;
	//double geodesic_dis = center_2d.length();
	//double angle_on_2d = 2.0 * M_PI - center_2d.angle();
	//double angle_sum_of_vert = model.AngleSum(source_vert);
	//double angle_on_mesh = angle_on_2d / (2.0 * M_PI) * angle_sum_of_vert;
	//const auto& e = model.Edge(direction_edge);
	//int v1 = e.indexOfLeftVert;
	//int v2 = e.indexOfRightVert;
	//int e_id_yx = path_tracer.metric.getEdgeFrom2Verts(v1, v2);
	//YXPath tmp_path;
	//path_tracer.computeSinglePathArbitraryStart(e_id_yx, angle_on_mesh, geodesic_dis, tmp_path);
	//center_3d = CenterPoint(path_tracer.convertToYXPoint3D(tmp_path.back()).toCPoint3D());
	//

}


double compute_Angle_Dest_From_MMP(geodesic::Mesh& mesh, int dest_index, geodesic::GeodesicAlgorithmBase* algorithm,
	const CRichModel& model, const CPoint3D& p_source)
{
	auto dest_vert(geodesic::SurfacePoint(&mesh.vertices()[dest_index]));
	vector<geodesic::SurfacePoint> path;
	algorithm->trace_back(dest_vert, path);
	geodesic::SurfacePoint p_next_to_dest;
	if (path.size() > 1) {
		p_next_to_dest = *(path.begin() + 1);
	}
	else {
		p_next_to_dest = (path[0]);
	}
	double angle = 0;
	auto& neighs = model.Neigh(dest_index);
	vector<int> neigh_verts;
	vector<double> neigh_angles;
	neigh_verts.reserve(neighs.size());
	neigh_angles.reserve(neigh_angles.size());
	for (auto& neigh : neighs) {
		neigh_verts.push_back(model.Edge(neigh.first).indexOfRightVert);
		neigh_angles.push_back(neigh.second);
	}
	auto& sum_angle = model.NeighAngleSum(dest_index);

	//printf("_______type\n");
	if (p_next_to_dest.type() == geodesic::VERTEX) {
		//printf(" vertex \n");
	}
	else if (p_next_to_dest.type() == geodesic::EDGE) {
		//printf(" edge \n");
	}
	else if (p_next_to_dest.type() == geodesic::FACE) {
		//printf(" face \n");
	}
	else if (p_next_to_dest.type() == geodesic::UNDEFINED_POINT) {
		//printf(" undefined_point \n");
	}

	if (p_next_to_dest.type() == geodesic::VERTEX) {
		bool flag_found = false;
		for (int j = 0; j < neigh_verts.size(); ++j) {
			if (p_next_to_dest.base_element()->id() == neigh_verts[j]) {
				flag_found = true;
				angle = sum_angle[j];
				break;
			}
		}
		if (!flag_found) {
			angle = 0;
		}
	}
	else if (p_next_to_dest.type() == geodesic::EDGE) {
		int v0 = p_next_to_dest.base_element()->adjacent_vertices()[0]->id();
		int v1 = p_next_to_dest.base_element()->adjacent_vertices()[1]->id();
		bool flag = false;
		for (int j = 0; j < neigh_verts.size(); ++j) {
			if (v0 == neigh_verts[j]) {
				int jminus1 = (j - 1 + neighs.size()) % neigh_verts.size();
				int vjminus1 = neigh_verts[jminus1];
				int jplus1 = (j + 1) % neigh_verts.size();
				int vjplus1 = neigh_verts[jplus1];
				CPoint3D p_cpoint3d(p_next_to_dest.x(), p_next_to_dest.y(), p_next_to_dest.z());
				if (v1 == vjminus1) {//v1 first
					double l = model.Edge(neighs[jminus1].first).length;
					double r = (model.Vert(dest_index) - p_cpoint3d).Len();
					double b = (model.Vert(vjminus1) - p_cpoint3d).Len();
					angle = sum_angle[jminus1] + acos((l * l + r * r - b * b) / (2 * l * r));
				}
				else if (v1 == vjplus1) {//v0 first
					double l = model.Edge(neighs[j].first).length;
					double r = (model.Vert(dest_index) - p_cpoint3d).Len();
					double b = (model.Vert(v0) - p_cpoint3d).Len();
					angle = sum_angle[j] + acos((l * l + r * r - b * b) / (2 * l * r));
				}
				else {
					fprintf(stderr, "fuck error line 680\n");
				}
				flag = true;
				break;
			}
		}
		if (!flag) {
			fprintf(stderr, "flag %d\n", flag);
		}
	}
	else if (p_next_to_dest.type() == geodesic::FACE) {
		angle = 0;
	}
	else {
		fprintf(stderr, "fuck error!\n");
	}
	//printf("angle %lf\n", angle);
	return angle;
}


double compute_Angle_From_MMP(geodesic::Mesh& mesh, int dest_index, geodesic::GeodesicAlgorithmBase* algorithm,
	const vector<int>& neigh_verts, const vector<double>& neigh_angles,
	const vector<double>& sum_angle, const CRichModel& model, const CPoint3D& p_source)
{
	auto dest_vert(geodesic::SurfacePoint(&mesh.vertices()[dest_index]));
	vector<geodesic::SurfacePoint> path;
	algorithm->trace_back(dest_vert, path);
	geodesic::SurfacePoint p_last_point;
	if (path.size() > 1) {
		p_last_point = *(path.rbegin() + 1);
	}
	else {
		p_last_point = (path[0]);
	}

	double angle = 0;
	if (p_last_point.type() == geodesic::VERTEX) {
		//printf("vertex\n");
		bool flag_found = false;
		for (int j = 0; j < neigh_verts.size(); ++j) {
			int neigh_vert = neigh_verts[j];
			if (p_last_point.base_element()->id() == neigh_vert) {
				//printf("yes\n");
				flag_found = true;
				angle = sum_angle[j];
				break;
			}
		}
		if (!flag_found) {
			angle = 0;
		}
	}
	else if (p_last_point.type() == geodesic::EDGE) {
		//printf("edge\n");
		int v0 = p_last_point.base_element()->adjacent_vertices()[0]->id();
		int v1 = p_last_point.base_element()->adjacent_vertices()[1]->id();
		bool flag = false;
		for (int j = 0; j < neigh_verts.size(); ++j) {
			double neigh_vert = neigh_verts[j];
			if (v0 == neigh_vert) {
				int jminus1 = (j - 1 + neigh_verts.size()) % neigh_verts.size();
				int vjminus1 = neigh_verts[jminus1];
				int jplus1 = (j + 1) % neigh_verts.size();
				int vjplus1 = neigh_verts[jplus1];
				CPoint3D p_cpoint3d(p_last_point.x(), p_last_point.y(), p_last_point.z());

				if (v1 == vjminus1) {//v1 first
					double l = (p_source - model.Vert(vjminus1)).Len();//  model.Edge(neighs[jminus1].first).length;
					double r = (p_source - p_cpoint3d).Len();
					double b = (model.Vert(vjminus1) - p_cpoint3d).Len();
					angle = sum_angle[jminus1] + acos((l * l + r * r - b * b) / (2 * l * r));
				}
				else if (v1 == vjplus1) {//v0      
					double l = (p_source - model.Vert(neigh_vert)).Len();//    model.Edge(neighs[j].first).length;
					double r = (p_source - p_cpoint3d).Len();
					double b = (model.Vert(v0) - p_cpoint3d).Len();
					angle = sum_angle[j] + acos((l * l + r * r - b * b) / (2 * l * r));
				}
				else{
					fprintf(stderr, "fuck error line 680\n");
				}
				flag = true;
				break;
			}
		}
		if (!flag) {
			fprintf(stderr, "flag %d\n", flag);
		}

	}
	else {
		fprintf(stderr, "fuck error face\n");
	}
	return angle;
}



void compute_point_on_face_neighbors(int current_source_index, const CPoint3D& p_source, int face_index_source, // input 
	geodesic::Mesh& mesh, double eps_vg, double theta, const CRichModel& model,// input
	BodyHeadOfSVG& body_header, vector<BodyPartOfSVGWithAngle>& body_parts_with_angle, vector<double>& dest_angles) // output
{
	vector<geodesic::SurfacePoint> sources;
	sources.push_back(geodesic::SurfacePoint(&mesh.faces()[face_index_source], p_source));
	//printf("line 3023\n");
	geodesic::GeodesicAlgorithmVGMMP algorithm(&mesh);
	double step = 1.0;
	const double eta = 100;
	algorithm.step = step;
	algorithm.binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
	map<int, double> fixedDests;
	//printf("line 3031\n");
	algorithm.propagate_vg(sources, eps_vg, fixedDests);

	vector<pair<int, double>> dests;

	for (auto& d : fixedDests) {
		geodesic::SurfacePoint dest_p = geodesic::SurfacePoint(&mesh.vertices()[d.first]);
		double dis;
		//printf("line 274\n");
		algorithm.best_source(dest_p, dis);
		dests.push_back(make_pair(d.first, d.second));
	}
	//printf("line 3043\n");
	body_header = BodyHeadOfSVG(current_source_index, dests.size());

	vector<int> neigh_verts;
	vector<double> neigh_angles;
	model.PointOnFaceNeigh(face_index_source, p_source, neigh_verts, neigh_angles);
	vector<double> sum_angle(neigh_angles.size() + 1);
	sum_angle[0] = 0;
	for (int j = 1; j <= neigh_angles.size(); ++j) {
		sum_angle[j] = sum_angle[j - 1] + neigh_angles[j - 1];
	}

	vector<double> angles(dests.size());
	for (int i = 0; i < dests.size(); ++i) {
		angles[i] = compute_Angle_From_MMP(mesh, dests[i].first, &algorithm,
			neigh_verts, neigh_angles,
			sum_angle, model, p_source);
	}

	dest_angles.clear();
	dest_angles.resize(dests.size());
	for (int i = 0; i < dests.size(); ++i) {
		dest_angles[i] = compute_Angle_Dest_From_MMP(mesh, dests[i].first, &algorithm, model, p_source);
	}
	//printf("line 3012 computed dest_angles \n");


	body_parts_with_angle.resize(dests.size());
	for (int i = 0; i < dests.size(); ++i) {
		auto b = dests[i];
		BodyPartOfSVGWithAngle b_with_angle(b.first, b.second, angles[i], 0, 0);
		body_parts_with_angle[i] = b_with_angle;
	}
	sort(body_parts_with_angle.begin(), body_parts_with_angle.end());

	//printf("line 3025\n");
	double angle_sum = sum_angle.back();
	vector<double> tmp_angles(body_parts_with_angle.size() * 2);
	for (int i = 0; i < body_parts_with_angle.size(); ++i) {
		tmp_angles[i] = body_parts_with_angle[i].angle;
	}
	for (int i = body_parts_with_angle.size(); i < tmp_angles.size(); ++i) {
		tmp_angles[i] = body_parts_with_angle[i - body_parts_with_angle.size()].angle + angle_sum;
	}
	for (int i = 0; i < body_parts_with_angle.size(); ++i) {//assume i is father
		double father_angle = body_parts_with_angle[i].angle;
		//based on father_angle as 0
		double start_angle = M_PI - theta + father_angle;
		double end_angle = angle_sum - (M_PI - theta) + father_angle;
		if (start_angle > end_angle) {
			body_parts_with_angle[i].begin_pos = -1;
			body_parts_with_angle[i].end_pos = -1;
			continue;
		}

		int start_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), start_angle) - tmp_angles.begin();
		if (start_pos > 0) start_pos--;
		int end_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), end_angle) - tmp_angles.begin();

		if (start_pos >= body_parts_with_angle.size()) start_pos -= body_parts_with_angle.size();
		if (end_pos >= body_parts_with_angle.size()) end_pos -= body_parts_with_angle.size();
		body_parts_with_angle[i].begin_pos = start_pos;
		body_parts_with_angle[i].end_pos = end_pos;
	}

	//delete algorithm;
	//algorithm = NULL;
}


void write_to_output_file(CRichModel& model, const vector<BodyHeadOfSVG>& total_heads, const vector<vector<BodyPartOfSVGWithAngle>>& total_parts, const string& output_file_name)
{
		int num_of_vertex = total_heads.size();
		ofstream output_file(output_file_name.c_str(), ios::out | ios::binary);
		HeadOfSVG head_of_svg(0, num_of_vertex - 1, num_of_vertex);
		output_file.write((char*)&head_of_svg, sizeof(head_of_svg));
		for (int i = 0; i < total_heads.size(); ++i) {
			output_file.write((char*)&total_heads[i], sizeof(total_heads[i]));
			for (auto& b : total_parts[i]) {
				output_file.write((char*)&b, sizeof(b));
			}
		}

		output_file.close();
}


void add_points_list(CRichModel& model, SparseGraph<float>* s_graph, double eps_vg, YXPathTracer& path_tracer, geodesic::Mesh& mesh, double theta, vector<BodyHeadOfSVG>& total_heads, vector<vector<BodyPartOfSVGWithAngle>>& total_parts, const string& svg_file_name)
{
	double min_angle = asin(sqrt(eps_vg)) * 2;
	bool flag_first_time_dump = false;
	int flag_dump_cnt = 0;
	vector<PointOnFace> added_points_list;

	vector<geodesic::SurfacePoint> surface_pts;
	bool flag_first = true;
	
	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		if (false) {
			if (i == -1 && flag_first) {
				i = 1700;
			}
		}
		printf("i%d ", i);
		const vector<float>& geo_dises = s_graph->graphNeighborDis(i); //i点的到各个邻居的//geodesic distance 
		const std::vector<float>&  angles = s_graph->graphNeighborAngle(i); //i点到各个邻居
		int neighor_size = angles.size();
		//的angle 
		//第一条边
		if (model.Neigh(i).size() == 0) continue;
		int first_edge_id = model.Neigh(i)[0].first;
		double angle_sum = model.AngleSum(i);
		double max_dis = *max_element(geo_dises.begin(), geo_dises.end());
		if (false) {
			CPoint3D p;
			int face_index;
			for (int j = 0; j < model.Neigh(i).size(); ++j) {
				int e_id = model.Neigh(i)[j].first;
				printf("e_id %d v0 %d v1 %d\n", e_id, model.Edge(e_id).indexOfLeftVert, model.Edge(e_id).indexOfRightVert);
			}
			for (int j = 0; j < 10; ++j) {
				printf("j xxxx %d\n", j);
				
				find_pseudo_point_using_path_tracer(model, i, first_edge_id, max_dis * (double)(j+1) / 10,
					(double)j/10 * 2 * PI, path_tracer, p, face_index);
				CylinderPath path_tracer(0.001);
				auto source = geodesic::SurfacePoint(&mesh.vertices()[i]);
				auto dest = geodesic::SurfacePoint(&mesh.faces()[face_index], p);
				printf("v0 %d v1 %d v2 %d\n", model.Face(face_index)[0], model.Face(face_index)[1], model.Face(face_index)[2]);
				printf("xyz %lf %lf %lf\n", dest.x(), dest.y(), dest.z());
				path_tracer.addGeodesicPath(mesh, source, dest);
				path_tracer.write_to_file("test_traced_edges" + to_string(j) + ".obj");
			}
			break;
		}
		if (false) {
			CPoint3D p;
			int face_index;
			for (int j = 0; j < neighor_size; ++j) {
				printf("j xxxx %d\n", j);
				find_pseudo_point_using_path_tracer(model, i, first_edge_id, max_dis,
					angles[j], path_tracer, p, face_index);

				CylinderPath path(0.001);
				int d0 = s_graph->graphNeighbor(i)[j];
				path.addGeodesicPath(model, i, d0);
				path.write_to_file("neighor_edges" + to_string(j) + ".obj");
				CylinderPath path_tracer(0.001);
				auto source = geodesic::SurfacePoint(&mesh.vertices()[i]);
				auto dest = geodesic::SurfacePoint(&mesh.faces()[face_index], p);
				path_tracer.addGeodesicPath(mesh, source, dest);
				path_tracer.write_to_file("traced_edges" + to_string(j) + ".obj");

			}
			break;
		}

		for (int j = 0; j < neighor_size; ++j) {
			double angles_diff{ 0 };
			if (j != angles.size() - 1) {
				angles_diff = angles[j + 1] - angles[j];
			} else {
				angles_diff = angle_sum + angles[0] - angles[j];
			}
			int pos_in_percent = int(angles_diff / min_angle);
			if (pos_in_percent >= 2) {
				vector<double> divided_angles(pos_in_percent - 1);
				for (int k = 0; k < divided_angles.size(); ++k) {
					divided_angles[k] = angles[j] + (double)(k + 1) / (double)pos_in_percent * angles_diff;
					if (divided_angles[k] > angle_sum) {
						divided_angles[k] -= angle_sum;
					}
					if (i == 1700) {
						printf("angles[%d]=%lf angles[%d]=%lf, k %d angles %lf\n", j, angles[j], j + 1, angles[j + 1], k, divided_angles[k]);
					}

					CPoint3D p;
					int face_index;
					find_pseudo_point_using_path_tracer(model, i, first_edge_id, max_dis,
						divided_angles[k], path_tracer, p, face_index);
					int current_source_index = model.GetNumOfVerts() + added_points_list.size();
					added_points_list.push_back(PointOnFace(face_index, p));
					total_heads[i].neighbor_num++;
					total_parts[i].push_back(BodyPartOfSVGWithAngle(current_source_index, max_dis, divided_angles[k], 0, 0));

					if (i == 1700) {
						surface_pts.push_back(geodesic::SurfacePoint(&mesh.faces()[face_index], p));
					}
					BodyHeadOfSVG body_head;
					std::vector<BodyPartOfSVGWithAngle> body_parts_with_angle;
					vector<double> dest_angles;
					compute_point_on_face_neighbors(current_source_index, p, face_index,	mesh, eps_vg, theta, model,// input
						body_head, body_parts_with_angle, dest_angles); // output
					total_heads.push_back(body_head);
					total_parts.push_back(body_parts_with_angle);
					for (int t = 0; t < body_parts_with_angle.size(); ++t) {
						int dest = body_parts_with_angle[t].dest_index;
						double dis = body_parts_with_angle[t].dest_dis;
						total_heads[dest].neighbor_num++;
						total_parts[dest].push_back(BodyPartOfSVGWithAngle(current_source_index, dis, dest_angles[t], 0, 0));
					}
				}
			}
		}

		if (false) {
			if (i == 1700) {
				//surface_pts.push_back(geodesic::SurfacePoint(&mesh.faces()[face_index], p));
				CylinderPath path(0.0005);
				path.addGeodesicPath(mesh, geodesic::SurfacePoint(&mesh.vertices()[i]), surface_pts);
				path.write_to_file("path_pseudo_edges.obj");
			}

			if (i == 1700 && flag_first) {
				i = -1;
				flag_first = false;
			}
		}
	}



	printf("\narrange angeles\n");
	for (int i = 0; i < model.GetNumOfVerts(); ++i)
	{
		printf("i %d ", i);
		auto& current_head = total_heads[i];
		auto& current_parts = total_parts[i];
		auto& one_ring_sum_angle = model.NeighAngleSum(i);
		sort(current_parts.begin(), current_parts.end());
		double angle_sum = one_ring_sum_angle.back();
		vector<double> tmp_angles(current_parts.size() * 2);
		for (int i = 0; i < current_parts.size(); ++i) {
			tmp_angles[i] = current_parts[i].angle;
		}
		for (int i = current_parts.size(); i < tmp_angles.size(); ++i) {
			tmp_angles[i] = current_parts[i - current_parts.size()].angle + angle_sum;
		}
		for (int i = 0; i < current_parts.size(); ++i) {//assume i is father
			double father_angle = current_parts[i].angle;
			//based on father_angle as 0
			double start_angle = M_PI - theta + father_angle;
			double end_angle = angle_sum - (M_PI - theta) + father_angle;
			if (start_angle > end_angle) {
				current_parts[i].begin_pos = -1;
				current_parts[i].end_pos = -1;
				continue;
			}

			int start_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), start_angle) - tmp_angles.begin();
			if (start_pos > 0) start_pos--;
			int end_pos = lower_bound(tmp_angles.begin(), tmp_angles.end(), end_angle) - tmp_angles.begin();
			if (start_pos >= current_parts.size()) start_pos -= current_parts.size();
			if (end_pos >= current_parts.size()) end_pos -= current_parts.size();
			current_parts[i].begin_pos = start_pos;
			current_parts[i].end_pos = end_pos;
		}
	}


	string output_file_name = svg_file_name.substr(0, svg_file_name.length() - 7) + "_fixed.binary";

	write_to_output_file(model, total_heads, total_parts, output_file_name);

	printf("added_points / verts_num %lf\n", (double)added_points_list.size() / model.GetNumOfVerts());
	
	//printf("added_edges / origin_edges %lf\n" ,    )
	//printBallToObj(added_points_list, "added_points_list.obj", 0.001);
	//delete s_graph;

}


void svg_precompute_add_pseudo_edges_main(const string& input_file_name, double eps_vg, double const_for_theta, const string& svg_file_name)
{
	double theta = asin(sqrt(eps_vg)) * const_for_theta;    //compute theta
	//initial model
	CRichModel model(input_file_name);
	model.Preprocess();
	//initial s_graph
	SparseGraph<float>* s_graph = new LC_HY<float>();
	s_graph->read_svg_file_with_angle((string)svg_file_name);
	WXN_PRINT_DEBUG_LINE("reading s_graph done");
	//initial mesh
	geodesic::Mesh mesh;
	std::vector<double> points;
	std::vector<unsigned> faces;
	std::vector<int> realIndex;
	int originalVertNum = 0;
	clock_t start = clock();
	bool success = geodesic::read_mesh_from_file(input_file_name.c_str(), points, faces, realIndex, originalVertNum);
	if (!success)
	{
		fprintf(stderr, "something is wrong with the input file");
		return;
	}
	mesh.initialize_mesh_data(points, faces);		//create internal
	WXN_PRINT_DEBUG_LINE("initial mesh done");

	YXPathTracer path_tracer;
	path_tracer.init(input_file_name.c_str());
	WXN_PRINT_DEBUG_LINE("path tracer done");

	vector<BodyHeadOfSVG> total_heads(model.GetNumOfVerts());
	vector<vector<BodyPartOfSVGWithAngle>> total_parts(model.GetNumOfVerts());

	initial_heads_and_parts(total_heads, total_parts, model, s_graph);
	WXN_PRINT_DEBUG_LINE("initial heads and parts done");

	add_points_list(model, s_graph, eps_vg, path_tracer, mesh, theta, total_heads, total_parts, svg_file_name);

	delete s_graph;
}

void svg_precompute_add_pseudo_edges(const string& input_file_name, double eps_vg, double const_for_theta, const string& svg_file_name)
{
	svg_precompute_add_pseudo_edges_main(input_file_name, eps_vg, const_for_theta, svg_file_name);
	_CrtDumpMemoryLeaks();
	system("pause");
}




