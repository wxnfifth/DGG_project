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

int findDest(CRichModel& model, string& input_obj_name, int source, double eps_vg)
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

	return fixedDests.rbegin()->first;

}

void figure_pruning_bunny()
{
	string input_obj_name = "bunny_nf10k.obj";
	CRichModel model(input_obj_name);
	model.Preprocess();
	int source = 4183;
	double eps_vg = 0.001;
	double cylinder_radius = 0.0004;
	int dest;
	dest = findDest(model, input_obj_name, source, eps_vg);

	{
		map<int, double> fixedDests;
		CRichModel model(input_obj_name);
		model.Preprocess();
		getDests(input_obj_name, model, source, eps_vg, fixedDests);
		CylinderPath cylinder1(cylinder_radius);
		vector<int> dests;
		for (auto& d : fixedDests) {
			if (d.first != dest) {
				dests.push_back(d.first);
			}
		}
		cylinder1.addGeodesicPaths(model, source, dests);
		cylinder1.write_to_file("bunny_dgg_edges.obj");
	}

	vector<int> sources;
	sources.push_back(source);
	CICHWithFurtherPriorityQueue algo_source(model, sources);
	algo_source.Execute();

	vector<int> dests;
	dests.push_back(dest);
	CICHWithFurtherPriorityQueue algo_dest(model, dests);
	algo_dest.Execute();

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
	int mid_vert_pos = 3;
	auto mid_error = vert_error[mid_vert_pos];
	int mid_vert = mid_error.first;
	printf("choose vert(i %d) error %lf \n", mid_error.first, mid_error.second);

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

}


int main()
{
	figure_pruning_bunny();


	return 0;
}