#include "ICH\BaseModel.h"
#include "ICH\RichModel.h"
#include "wxn\wxnmath.h"
#include "wxn\wxnTime.h"
#include "wxn\wxn_path_helper.h"
#include "YXMetric\YXApproxBlueNoise.h"
#include "wxn\wxn_dijstra.h"

void testApproxBlueNoise_tmp(){
	YXApproxBlueNoise bn;
	YXTimer t;
	t.start();
	bn.load("C:\\Users\\YingXiang\\GamelabSVN\\YingXiang\\3DModels\\bunny\\bunny_nf144k.m");
	printf("Time used for loading: %f seconds\n", t.gettime());
	int N, P;
	//t.start();
	//N = 10000;
	//P = 10;
	//bn.sampling(N, P);
	//printf("Time used for sampling: %lf seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
	//	t.gettime(), N, P, bn.result.size());

	t.start();
	N = 50000;
	P = 10;
	bn.sampling(N, P);
	printf("Time used for sampling: %f seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
		t.gettime(), N, P, bn.result.size());

	t.start();
	N = 200000;
	P = 10;
	bn.sampling(N, P);
	printf("Time used for sampling: %f seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
		t.gettime(), N, P, bn.result.size());

	t.start();
	N = 1000000;
	P = 5;
	bn.sampling(N, P);
	printf("Time used for sampling: %f seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
		t.gettime(), N, P, bn.result.size());
}

CPoint3D YXPoint2Point3D(const YXMesh3D& mesh_3d,const YXMeshPoint& mesh_point)
{
	auto& f = mesh_3d.face[mesh_point.faceid];
	double a = mesh_point.a;
	double b = mesh_point.b;
	double c = 1 - mesh_point.a - mesh_point.b;
	YXPoint3D p = mesh_3d.vert[f.v[0]] * a + mesh_3d.vert[f.v[1]] * b + mesh_3d.vert[f.v[2]] * c;
	return CPoint3D(p.x, p.y, p.z);
}

void add_sampling_noise(string& input_name)
{
	YXApproxBlueNoise bn;
	vector<PointOnFace> sampled_pts;
	bn.load(input_name.c_str());
	int N = 10000;
	int P = 10;
	bn.sampling(N, P);
	//bn.result.size();
	sampled_pts.reserve(bn.result.size());
	for (auto& p : bn.result) {
		sampled_pts.push_back(PointOnFace(p.faceid, YXPoint2Point3D(bn.mesh, p)));
	}
	vector<CPoint3D> pts_for_output;
	for (auto& s : sampled_pts) {
		pts_for_output.push_back(s.v);
	}

	printBallToObj(pts_for_output, "fertility_sampled_pts.obj", 0.001);
	CRichModel model(input_name);
	model.Preprocess();
	vector<int> origin_faces;
	ModelForSubdivide sub_model(model, origin_faces);
	for (auto& p : sampled_pts) {
		int s;
		sub_model.addVertex(p.v, p.face_id_, s);
	}
	sub_model.subdivide();
	sub_model.WriteToFile("fertility_nf10k_ani2_sub.obj");

}

void addNoiseByAddingPointOnEdge(string& input_name)
{
	double ratio = 1.01;
	CRichModel model(input_name);
	model.Preprocess();
	double mean_edge_length = model.GetMeanEdgeLength();
	vector<pair<CPoint3D,int>> points_on_edge;
	for (int i = 0; i < model.GetNumOfSimpleEdges(); ++i) {
		auto& e = model.SimpleEdge(i);
		double len = (model.Vert(e.v1) - model.Vert(e.v2)).Len();
		if (len > mean_edge_length * ratio) {
			int num_of_points = int(len / mean_edge_length + 0.5) ;
			
			double seg_len = 1.0 / (num_of_points + 1);
			for (int j = 0; j < num_of_points; ++j) {
				double tmp_ratio = (j + 1) * seg_len;
				CPoint3D p = model.Vert(e.v1) * tmp_ratio + model.Vert(e.v2) * (1 - tmp_ratio);
				points_on_edge.push_back(make_pair(p, i));
			}
		}
	}


	vector<int> origin_faces;
	ModelForSubdivide sub_model(model, origin_faces);
	for (auto& p_on_edge : points_on_edge) {
		int s;
		sub_model.addVertexOnEdge(p_on_edge.first, p_on_edge.second);
	}
	sub_model.subdivide();
	sub_model.WriteToFile("fertility_nf10k_ani2_edgepoint.obj");


}

void printEdgeVariance()
{
	CRichModel model_an0("fertility_nf10k_ani0.obj");
	model_an0.Preprocess();
	printf("model ani0 %lf\n", model_an0.GetEdgeVariance());
	CRichModel model_an1("fertility_nf10k_ani1.obj");
	model_an1.Preprocess();
	printf("model ani1 %lf\n", model_an1.GetEdgeVariance());
	CRichModel model_an2("fertility_nf10k_ani2.obj");
	model_an2.Preprocess();
	printf("model ani2 %lf\n", model_an2.GetEdgeVariance());
	CRichModel model_an2_edgepoint("fertility_nf10k_ani2_edgepoint.obj");
	model_an2_edgepoint.Preprocess();
	printf("model ani2 %lf\n", model_an2_edgepoint.GetEdgeVariance());


}

void getDGGDistance(const CRichModel& model, int source_vert, const string& dgg_filename, vector<double>& dgg_dis)
{
	SparseGraph<float>* s_graph = NULL;
	s_graph = new LC_HY<float>();
	s_graph->read_svg_file_with_angle((string)dgg_filename);
	dynamic_cast<LC_HY<float>*>(s_graph)->setModel(model);
	ElapasedTime t;
	s_graph->findShortestDistance(source_vert);
	t.printTime("dgg time");
	dgg_dis.reserve(model.GetNumOfVerts());

	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		double dis = s_graph->distanceToSource(i);
		//printf("%lf ", dis);
		if (_finite(dis) && dis < 1e10) {
			dgg_dis.push_back(dis);
		}
		else {
			dgg_dis.push_back(0);
		}
	}
	delete s_graph;
}

void outputDistanceField(CRichModel& model, const vector<double>& dis, const string& output_filename)
{
	assert(model.GetNumOfVerts() == dis.size());

	vector<pair<double, double>> ich_texture;
	double max_dis = *std::max_element(dis.begin(), dis.end());
	printf("max_dis %lf\n", max_dis);
	ich_texture.reserve(model.GetNumOfVerts());
	for (double d : dis) {
		ich_texture.push_back(make_pair(d / max_dis, d / max_dis));
		//ich_texture.push_back(make_pair(d, d));
	}
	model.FastSaveObjFile(output_filename, ich_texture);

}


void getICHDistance(const CRichModel& rich_model, int source_vert, vector<double>& correct_dis)
{

	ElapasedTime t;
	vector<int> sources;
	sources.push_back(source_vert);
	CICHWithFurtherPriorityQueue ich_algorithm(rich_model, sources);
	ich_algorithm.Execute();
	t.printTime("ich");
	correct_dis.resize(rich_model.GetNumOfVerts());
	for (int i = 0; i < rich_model.GetNumOfVerts(); ++i) {
		double d = ich_algorithm.m_InfoAtVertices[i].disUptodate;
		if (d == d) {
			correct_dis[i] = d;
		}
		else{
			correct_dis[i] = 0;
		}
	}

}

void computeError(const vector<double>& correct_dis, const vector<double>& dis, vector<double>& error)
{
	assert(correct_dis.size() == dis.size());
	error.resize(dis.size());
	for (int i = 0; i < correct_dis.size(); ++i) {
		if (fabs(correct_dis[i]) < 1e-9 || !isfinite(correct_dis[i]) || !isfinite((dis[i]))) {
			error[i] = 0;
		}
		else{
			error[i] = fabs(correct_dis[i] - dis[i]) / correct_dis[i];
		}
	}
}

void computeStatics(const vector<double>& error)
{
	double min_error = 1e10;
	double max_error = -1e10;
	double average_error = 0;
	for (auto& e : error) {
		min_error = min(min_error, e);
		max_error = max(max_error, e);
		average_error += e;
	}
	average_error /= error.size();
	printf("min error is %lf\n", min_error);
	printf("max error is %lf\n", max_error);
	printf("average error is %lf\n", average_error);
}


void computeAndOutputDistanceField(const string& model_name, const string& dgg_name, int source, const string& distance_field_name)
{
	CRichModel model(model_name);
	model.Preprocess();
	vector<double> dgg_dis;
	getDGGDistance(model, source, dgg_name, dgg_dis);
	outputDistanceField(model, dgg_dis, distance_field_name);

	vector<double> correct_dis;
	getICHDistance(model, source, correct_dis);
	vector<double> dgg_error;
	computeError(correct_dis, dgg_dis, dgg_error);
	computeStatics(dgg_error);

}

void findNearestSource(const string& model_name_1, int source_1, string& model_name_2, int& source_2)
{
	CRichModel model_1(model_name_1);
	model_1.Preprocess();
	auto& p1 = model_1.Vert(source_1);
	CRichModel model_2(model_name_2);
	model_2.Preprocess();
	source_2 = 0;
	double min_dis = 1e10;
	for (int i = 0; i < model_2.GetNumOfVerts(); ++i) {
		double dis = (p1 - model_2.Vert(i)).Len();
		if (dis < min_dis) {
			min_dis = dis;
			source_2 = i;
		}
	}
	printf("min_dis %lf\n", min_dis);
}

void printBuddhaResult()
{
	string buddha_small_name = "buddha_nf40k\\buddha_nf40k.obj";
	string buddha_small_dgg_name = "buddha_nf40k\\buddha_nf40k_DGG0.010000_c20_pruning.binary";
	string buddha_small_distance_name = "buddha_nf40k\\buddha_nf40k_distance.obj";
	string buddha_middle_name = "buddha_nf300k\\buddha_nf300k.obj";
	string buddha_middle_dgg_name = "buddha_nf300k\\buddha_nf300k_DGG0.010000_c20_pruning.binary";
	string buddha_middle_distance_name = "buddha_nf300k\\buddha_nf300k_distance.obj";
	string buddha_large_name = "buddha_nf600k\\buddha_nf600k.obj";
	string buddha_large_dgg_name = "buddha_nf600k\\buddha_nf600k_DGG0.010000_c20_pruning.binary";
	string buddha_large_distance_name = "buddha_nf600k\\buddha_nf600k_distance.obj";
	int source_small = 5104;
	int source_middle;
	int source_large;
	findNearestSource(buddha_small_name, source_small, buddha_middle_name, source_middle);
	findNearestSource(buddha_small_name, source_small, buddha_large_name, source_large);

	computeAndOutputDistanceField(buddha_small_name, buddha_small_dgg_name, source_small, buddha_small_distance_name);
	computeAndOutputDistanceField(buddha_middle_name, buddha_middle_dgg_name, source_middle, buddha_middle_distance_name);
	computeAndOutputDistanceField(buddha_large_name, buddha_large_dgg_name, source_large, buddha_large_distance_name);
}


int main()
{
	//string input_name = "fertility_nf10k_ani2.obj";
	//add_sampling_noise(input_name);
	
	//addNoiseByAddingPointOnEdge(input_name);

	//printEdgeVariance();

	//5104

	printBuddhaResult();



	return 0;
}