#include "stdafx.h"
#include "DGG_precompute\svg_precompute.h"
#include "MMP\geodesic_algorithm_vg_mmp.h"
#include "ICH\RichModel.h"
#include "WXN\wxn_path_helper.h"

void getDests(string& input_obj_name, const CRichModel& model, int source_index, double eps_vg, map<int, double>& fixedDests)
{
	std::vector<double> points;
	std::vector<unsigned> faces;
	std::vector<int> realIndex;
	int originalVertNum = 0;

	clock_t start = clock();
	bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(), points, faces, realIndex, originalVertNum);
	if (!success)
	{
		fprintf(stderr, "something is wrong with the input file\n");
		return;
	}
	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);		//create internal
	vector<int> srcs;
	srcs.push_back(source_index);
	vector<geodesic::SurfacePoint> sources;
	for (unsigned i = 0; i < srcs.size(); ++i)
	{
		srcs[i] %= originalVertNum;
		srcs[i] = realIndex[srcs[i]];
		sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
	}
	geodesic::GeodesicAlgorithmBase *algorithm;

	algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh);

	double step = 1.0;
	const double eta = 100;

	algorithm->step = step;
	algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
	algorithm->propagate_vg(sources, eps_vg, fixedDests);

}

void figure5()
{
	string input_obj_name = "bunny_nf10k.obj";
	int source_index = 4183;
	double eps_vg = 1e-3;
	int dest_index = 4169;
	//目的地 4169
	int middle_index = 1134;
	//中间点1134
	map<int, double> fixedDests;
	CRichModel model(input_obj_name);
	model.Preprocess();
	getDests(input_obj_name, model, source_index, eps_vg, fixedDests);
	double cylinder_radius = 0.0002;
	CylinderPath cylinder1(cylinder_radius);
	vector<int> dests;
	for (auto& d : fixedDests) {
		dests.push_back(d.first);
	}
	cylinder1.addGeodesicPaths(model, source_index, dests);
	cylinder1.write_to_file("bunny_dgg_edges.obj");
	
	vector<CPoint3D> saddle_pts;
	for (int d = 0; d < model.GetNumOfVerts(); ++d) {
		if (model.AngleSum(d) > 2 * PI + 1e-6) {
			saddle_pts.push_back(model.Vert(d));
		}
	}
	printBallToObj(saddle_pts, "bunny_saddles.obj", 0.0012);

	{
		vector<int> sources;
		sources.push_back(source_index);
		CICHWithFurtherPriorityQueue algo_source(model, sources);
		algo_source.Execute();

		vector<int> dests;
		dests.push_back(dest_index);
		CICHWithFurtherPriorityQueue algo_dest(model, dests);
		algo_dest.Execute();

		double dis_exact = algo_source.m_InfoAtVertices[dest_index].disUptodate;
		
		double d0 = algo_source.m_InfoAtVertices[middle_index].disUptodate;
		double d1 = algo_dest.m_InfoAtVertices[middle_index].disUptodate;
		//vert_dis.push_back(make_pair(i, d0 + d1));
		//vert_error.push_back(make_pair(i, (d0 + d1 - dis_exact) / dis_exact));
		
		double error = fabs(dis_exact - d0 + d1) / dis_exact;
		printf("error %lf\n", error);

		CylinderPath cylinder_path_eps(cylinder_radius);
		cylinder_path_eps.addGeodesicPath(model, source_index, middle_index);
		cylinder_path_eps.addGeodesicPath(model, middle_index, dest_index);
		cylinder_path_eps.write_to_file("bunny_path_eps.obj");

		CylinderPath cylinder_path(cylinder_radius);
		cylinder_path.addGeodesicPath(model, source_index, dest_index);
		cylinder_path.write_to_file("bunny_path_staight.obj");

	}




}

int findDest(CRichModel& model, string& input_obj_name, int source,
	double eps_vg, int edge_v0, int edge_v1, double cylinder_radius)
{
	std::vector<double> points;
	std::vector<unsigned> faces;
	std::vector<int> realIndex;
	int originalVertNum = 0;
	//CRichModel model(input_obj_name);
	//model.Preprocess();

	clock_t start = clock();
	bool success = geodesic::read_mesh_from_file(input_obj_name.c_str(), points, faces, realIndex, originalVertNum);
	if (!success)
	{
		fprintf(stderr, "something is wrong with the input file\n");
		exit(1);
	}
	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);		//create internal


	int source_index = source;
	vector<int> srcs;
	srcs.push_back(source_index);
	vector<geodesic::SurfacePoint> sources;
	for (unsigned i = 0; i < srcs.size(); ++i)
	{
		srcs[i] %= originalVertNum;
		srcs[i] = realIndex[srcs[i]];
		sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[srcs[i]]));
	}


	geodesic::GeodesicAlgorithmBase *algorithm;

	algorithm = new geodesic::GeodesicAlgorithmVGMMP(&mesh);

	double step = 1.0;
	const double eta = 100;

	algorithm->step = step;
	algorithm->binWidth = mesh.avg_edge() / sqrt((double)mesh.vertices().size()) * eta;
	map<int, double> fixedDests;
	algorithm->propagate_vg(sources, eps_vg, fixedDests);

	//int dest = fixedDests.rbegin()->first;
	//int dest = fixedDests.rbegin()->first;
	auto& itr = fixedDests.rbegin();
	while (model.AngleSum(itr->first) + 1e-5 <= 2 * M_PI) {
		itr++;
		printf("__angle = %lf * 2PI\n", model.AngleSum(itr->first) / 2.0 / M_PI );
	}
	int dest = itr->first;

	std::vector<geodesic::SurfacePoint> trace_path;
	geodesic::interval_pointer result_interval = NULL;
	algorithm->trace_back_find_edge(geodesic::SurfacePoint(&mesh.vertices()[dest]),
									trace_path, edge_v0, edge_v1, result_interval);
	{
		std::vector<geodesic::SurfacePoint> path;
		geodesic::interval_pointer last_interval;
		algorithm->trace_back_interval(geodesic::SurfacePoint(&mesh.vertices()[dest]),
			trace_path, last_interval);
		auto& e_ptr = last_interval->edge();
		double x0 = last_interval->start() / e_ptr->length();
		double x1 = last_interval->stop() / e_ptr->length();
		geodesic::SurfacePoint e_p0(e_ptr, x0);
		geodesic::SurfacePoint e_p1(e_ptr, x1);

		algorithm->trace_back(e_p0, path);
		CylinderPath cylinder_path(cylinder_radius * 2);
		//4803,1489
		geodesic::edge_pointer edge_d0 = mesh.get_edge(1489, 4803);
		geodesic::SurfacePoint edge_d0_midp(edge_d0, 0.8);
		cylinder_path.addGeodesicPath(mesh, sources[0], e_p0);
		geodesic::edge_pointer edge_d1 = mesh.get_edge(1487, 4803);
		geodesic::SurfacePoint edge_d1_midp(edge_d1, 0.8);
		cylinder_path.addGeodesicPath(mesh, sources[0], e_p1);
		cylinder_path.write_to_file("bunny_a_window_fan.obj");

		vector<CPoint3D> points{ e_p0.xyz(), e_p1.xyz() };
		printf("___________x0 %lf x1 %lf edge %d\n", x0, x1, e_ptr->id());
		printBallToObj(points, "bunny_a_interval_points.obj", cylinder_radius * 10);
	}

	//result_interval->start();

	if (result_interval != NULL) {

		auto& e_ptr = result_interval->edge();
		CPoint3D v0(e_ptr->v0()->xyz());
		CPoint3D v1(e_ptr->v1()->xyz());
		double x0 = result_interval->start() / e_ptr->length();
		double x1 = result_interval->stop() / e_ptr->length();
		geodesic::SurfacePoint e_p0(e_ptr, x0);
		geodesic::SurfacePoint e_p1(e_ptr, x1);
		//x0 = 0.1;
		//x1 = 0.6;
		printf("_______interval pos x0 %lf x1 %lf\n", x0, x1);
		//edge_p0 = (1 - x0)*v0 + x0 * v1;
		//(1 - a)*v0 + a*v1;
		//edge 1489,4803
		geodesic::edge_pointer edge_d0 = mesh.get_edge(1489, 4803);
		double len = (model.Vert(1489) - model.Vert(4803)).Len();
		double dis_to_4803 = 0.026;
		double ratio = 0.2388;
		geodesic::SurfacePoint edge_d0_midp(edge_d0, ratio);
		CylinderPath cylinder_edge_path(cylinder_radius * 2);
		cylinder_edge_path.addGeodesicPath(mesh, sources[0], edge_d0_midp);
		//edge 



		//4803,1487
		geodesic::edge_pointer edge_d1 = mesh.get_edge(1391, 1487);
		geodesic::SurfacePoint edge_d1_midp(edge_d1, 0.575195878072);
		cylinder_edge_path.addGeodesicPath(mesh, sources[0], edge_d1_midp);
		cylinder_edge_path.write_to_file("bunny_b_window_fan.obj");
		//vector<geodesic::SurfacePoint> path;
		//algorithm->trace_back(edge_d0_midp, path);




		geodesic::SurfacePoint edge_sp0(e_ptr, x0);
		geodesic::SurfacePoint edge_sp1(e_ptr, x1);
		//edge_p1 = (1 - x1)*v0 + x1 * v1;
		//interval_points.push_back(make_pair(p0, p1));
		CylinderPath cylinder_path(cylinder_radius * 2);
		cylinder_path.addGeodesicPath(mesh, sources[0], edge_sp0);
		cylinder_path.addGeodesicPath(mesh, sources[0], edge_sp1);
		cylinder_path.write_to_file("bunny_windows_fan.obj");
		//printf("____angle %lf\n", model.AngleSum(dest));
		printBallToObj(vector < CPoint3D > {edge_sp0.xyz(), edge_sp1.xyz()},
			"bunny_interval_points.obj", cylinder_radius * 2);

	}
	return dest;

}

void writeEdges(CRichModel& model, int source, int dest, double eps_vg, double cylinder_radius, const string& edge_filename)
{
	map<int, double> fixedDests;
	getDests(model.GetFullPathAndFileName(), model, source, eps_vg, fixedDests);
	CylinderPath cylinder1(cylinder_radius);
	vector<int> dests;
	for (auto& d : fixedDests) {
		if (d.first != dest) {
			dests.push_back(d.first);
		}
	}
	cylinder1.addGeodesicPaths(model, source, dests);
	cylinder1.write_to_file(edge_filename);
}

int getMidVert(CICHWithFurtherPriorityQueue& algo_source, CICHWithFurtherPriorityQueue& algo_dest,
			   CRichModel& model, int source, int dest)
{
	double dis_exact = algo_source.m_InfoAtVertices[dest].disUptodate;
	set<int> path_verts;
	vector<IntersectionWithPath> result_path;
	algo_source.FindSourceVertex(dest, result_path);
	path_verts.insert(source);
	path_verts.insert(dest);
	for (auto& p : result_path) {
		if (p.isVertex) {
			path_verts.insert(p.index);
		}
	}

	vector<pair<int, double>> vert_dis;
	vector<pair<int, double>> vert_error;
	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		if (path_verts.find(i) != path_verts.end()) continue;
		double d0 = algo_source.m_InfoAtVertices[i].disUptodate;
		double d1 = algo_dest.m_InfoAtVertices[i].disUptodate;
		vert_dis.push_back(make_pair(i, d0 + d1));
		vert_error.push_back(make_pair(i, (d0 + d1 - dis_exact) / dis_exact));
	}
	sort(vert_error.begin(), vert_error.end(),
		[](const pair<int, double>& lhs, const pair<int, double>& rhs) {
		return lhs.second < rhs.second; });
	for (int i = 0; i < 20; ++i) {
		printf("vert(i %d) error %lf \n", vert_error[i].first, vert_error[i].second);
	}
	int mid_vert_pos = 5;
	auto& mid_error = vert_error[mid_vert_pos];
	int mid_vert = mid_error.first;
	printf("choose vert(i %d) error %lf \n", mid_error.first, mid_error.second);
	return mid_vert;
}

void figure_5_b_bunny()
{
	string input_obj_name = "bunny_nf10k.obj";
	CRichModel model(input_obj_name);
	model.Preprocess();
	int source = 4183;
	double eps_vg = 0.001;
	double cylinder_radius = 0.0004;
	int dest;
	//get window pass edge (1393,1394)
	//int edge_v0 = 1393;
	//int edge_v1 = 1394;
	int edge_v0 = 4180;
	int edge_v1 = 4586;
	//CPoint3D edge_p0;
	//CPoint3D edge_p1;

	dest = findDest(model, input_obj_name, source, eps_vg, edge_v0, edge_v1, cylinder_radius);

	string edge_filename("bunny_dgg_edges.obj");
	writeEdges(model, source, dest, eps_vg, cylinder_radius, edge_filename);

	vector<int> sources;
	sources.push_back(source);
	CICHWithFurtherPriorityQueue algo_source(model, sources);
	algo_source.Execute();

	vector<int> dests{ dest };
	CICHWithFurtherPriorityQueue algo_dest(model, dests);
	algo_dest.Execute();

	int mid_vert = getMidVert(algo_source, algo_dest, model, source, dest);

	CylinderPath cylinder_path(cylinder_radius * 2);
	cylinder_path.addGeodesicPath(model, source, dest);
	cylinder_path.write_to_file("bunny_path_staight.obj");
	CylinderPath cylinder_path_eps(cylinder_radius * 2);
	cylinder_path_eps.addGeodesicPath(model, source, mid_vert);
	cylinder_path_eps.addGeodesicPath(model, mid_vert, dest);
	cylinder_path_eps.write_to_file("bunny_path_eps.obj");

	vector<CPoint3D> points{ model.Vert(source), model.Vert(dest) };
	printBallToObj(points, "bunny_source_dest_points.obj", cylinder_radius * 10);
	printBallToObj(vector < CPoint3D > {model.Vert(mid_vert)}, "bunny_mid_points.obj", cylinder_radius * 5);

	set<int> faces{ 1996 };
	output_faces(&model, faces, "bunny_a_face.obj");

}


void figure_5_b_bunny_2d()
{
	string input_obj_name = "bunny_nf10k.obj";
	CRichModel model(input_obj_name);
	model.Preprocess();
	int source = 4183;
	//1 based
	vector<int> face_list{ 7616, 7602, 7619, 7618, 7620, 8726, 8727, 8729, 8730, 8731,
		8728, 1973, 1972, 1971, 8384, 8385, 9229, 9228, 9362,8391 };//9362
	//for (int f_id : face_list) {
	for (int i = 0; i < face_list.size(); ++i) {
		int f_id = face_list[i];
		auto& f = model.Face(f_id - 1);
		double len[3];
		for (int j = 0; j < 3; ++j) {
			len[j] = (model.Vert(f[j]) - model.Vert(f[(j + 1) % 3])).Len();
		}
		if (i > 0){
			auto& pre_f = model.Face(face_list[i - 1] - 1);
			int new_v = -1;
			for (int j = 0; j < 3; ++j) {
				if (f[j] != pre_f[0] && f[j] != pre_f[1] && f[j] != pre_f[2]) {
					new_v = j;
					break;
				}
			}
			int old_v0, old_v1;
			double len_0, len_1;
			if (new_v == 0) {
				old_v0 = 1; old_v1 = 2;
				len_0 = len[0]; len_1 = len[2];
			}
			else if (new_v == 1) {
				old_v0 = 0; old_v1 = 2;
				len_0 = len[0]; len_1 = len[1];
			}
			else if (new_v == 2) {
				old_v0 = 0; old_v1 = 1;
				len_0 = len[2]; len_1 = len[1];
			}
			printf("face %d: v %d new_v %d %lf v %d new_v %d %lf\n",
				f_id, f[old_v0], f[new_v], len_0, f[old_v1], f[new_v], len_1);
		}
		else{
			printf("face %d: e %d %d %lf e %d %d %lf e %d %d %lf\n", f_id - 1,
				f[0], f[1], len[0],
				f[1], f[2], len[1],
				f[2], f[0], len[2]);
		}
	}

}


int main()
{
	figure_5_b_bunny();
	//figure_5_b_bunny_2d();





	return 0;
}