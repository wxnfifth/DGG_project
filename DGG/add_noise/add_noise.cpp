#include "ICH\BaseModel.h"
#include "ICH\RichModel.h"
#include "wxn\wxnmath.h"
#include "wxn\wxn_path_helper.h"
#include "YXMetric\YXApproxBlueNoise.h"

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
	int N = 5000;
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
	double ratio = 1.6;
	CRichModel model(input_name);
	model.Preprocess();
	double mean_edge_length = model.GetMeanEdgeLength();
	vector<pair<CPoint3D,int>> points_on_edge;
	for (int i = 0; i < model.GetNumOfSimpleEdges(); ++i) {
		auto& e = model.SimpleEdge(i);
		double len = (model.Vert(e.v1) - model.Vert(e.v2)).Len();
		if (len > mean_edge_length * ratio) {
			int num_of_points = int(len / mean_edge_length + 0.5) - 1;
			double seg_len = 1.0 / (num_of_points + 1);
			for (int j = 0; j < num_of_points; ++j) {
				double ratio = (j + 1) * seg_len;
				CPoint3D p = model.Vert(e.v1) * ratio + model.Vert(e.v2) * (1 - ratio);
				points_on_edge.push_back(make_pair(p, i));
			}
		}
	}


	vector<int> origin_faces;
	ModelForSubdivide sub_model(model, origin_faces);
	for (auto& p_on_edge : points_on_edge) {
		int s;
		sub_model.addVertex(p_on_edge.first, 0, s, p_on_edge.second);
	}
	sub_model.subdivide();
	sub_model.WriteToFile("fertility_nf10k_ani2_edgepoint.obj");


}


int main()
{
	string input_name = "fertility_nf10k_ani2.obj";
	//add_sampling_noise(input_name);
	
	addNoiseByAddingPointOnEdge(input_name);



	return 0;
}