
#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"
#include "wxn\wxn_dijstra.h"
#include "wxn\wxnTime.h"
#include "MMP\geodesic_algorithm_exact.h"
#include "YXMetric\YXMetric.h"
#include "wxn\wxn_path_helper.h"
#include <random>
//#include "figure_1\DelaunayMesh.h"

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
		if (d==d) {
			correct_dis[i] = d;
		}
		else{
			correct_dis[i] = 0;
		}
	}

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
		if (_finite(dis)) {
			dgg_dis.push_back(dis);
		} else {
			dgg_dis.push_back(0);
		}
	}
	delete s_graph;
}

void getSVGDistance(const CRichModel& model, int source_vert, const string& svg_filename, vector<double>& svg_dis)
{
	SparseGraph<float>* s_graph = NULL;
	s_graph = new Dijstra_vector<float>();
	s_graph->read_svg_file_binary((string)svg_filename);
	ElapasedTime t;
	s_graph->findShortestDistance(source_vert);
	t.printTime("svg_time");
	svg_dis.reserve(model.GetNumOfVerts());
	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		double dis = s_graph->distanceToSource(i);
		if (_finite(dis)) {
			svg_dis.push_back(dis);
		}
		else {
			svg_dis.push_back(0);
		}
	}
	delete s_graph;
}


void outputDistanceField(CRichModel& model, const vector<double>& dis, const string& output_filename)
{
	assert(model.GetNumOfVerts() == dis.size());

	vector<pair<double, double>> ich_texture;
	double max_dis = std::max(max_dis, *std::max_element(dis.begin(), dis.end()));
	ich_texture.reserve(model.GetNumOfVerts());
	for (double d : dis) {
		//ich_texture.push_back(make_pair(d / max_dis, d / max_dis));
		ich_texture.push_back(make_pair(d, d));
	}
	model.FastSaveObjFile(output_filename, ich_texture);

}

void computeError(const vector<double>& correct_dis, const vector<double>& dis, vector<double>& error)
{
	assert(correct_dis.size() == dis.size());
	error.resize(dis.size());
	for (int i = 0; i < correct_dis.size(); ++i) {
		if (fabs(correct_dis[i]) < 1e-9 || !isfinite(correct_dis[i]) || !isfinite((dis[i]))) {
			error[i] = 0;
		} else{
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
	printf("min error is %lf\n" , min_error);
	printf("max error is %lf\n", max_error);
	printf("average error is %lf\n", average_error);
}

void writeErrorFIle(vector<double>& error, const string& filename)
{
	FILE* file = fopen(filename.c_str(), "w");
	assert(file != NULL);
	for (auto& e : error) {
		fprintf(file, "%lf\n", e);
	}
	fclose(file);
}

void normalizeError(vector<double>& svg_error, double min_error, double max_error)
{
//	for (auto& e : svg_error) {
	for (int i = 0; i < svg_error.size(); ++i) {
		svg_error[i] = (svg_error[i] - min_error) / (max_error - min_error);
		//printf("'%lf %lf %lf'", svg_error[i], min_error, max_error);
		if (svg_error[i] >= 1){
			svg_error[i] = 1;
		}
	}
	//printf("\n");
}

void logNormalizeErrorAndOutput(vector<double>& svg_error, double min_error, double max_error, CRichModel& model, const string& output_filename)
{
	//	for (auto& e : svg_error) {
	for (int i = 0; i < svg_error.size(); ++i) {
		if (svg_error[i] < min_error) {
			svg_error[i] = min_error;
			//printf("min_error ");
		}
		if (svg_error[i] >= max_error){
			svg_error[i] = max_error;
			//printf("max_error ");
		}
		svg_error[i] = (log(svg_error[i]) - log(min_error)) / (log(max_error) - log(min_error));
		//printf("'%lf %lf %lf'", svg_error[i], min_error, max_error);
	}
	//printf("\n");
	outputDistanceField(model, svg_error, output_filename);
}



void figure_1(int source_vert)
{
	string obj_file_name = "gargoyle_nf700k.obj";
	string dgg_filename = "gargoyle_nf700k_DGG0.000100_c20_pruning.binary";
	string svg_filename = "gargoyle_nf700k_SVG_k483.binary";

	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();

	vector<double> correct_dis;

	getICHDistance(model, source_vert, correct_dis);
	outputDistanceField(model, correct_dis, obj_prefix + "_ich.obj");

	vector<double> dgg_dis;
	getDGGDistance(model, source_vert, dgg_filename, dgg_dis);
	outputDistanceField(model, dgg_dis, obj_prefix + "_dgg_dis.obj");

	vector<double> svg_dis;
	getSVGDistance(model, source_vert, svg_filename, svg_dis);
	outputDistanceField(model, svg_dis, obj_prefix + "_svg_dis.obj");

	printf("dgg_statis:\n");
	vector<double> dgg_error;
	computeError(correct_dis, dgg_dis, dgg_error);
	writeErrorFIle(dgg_error, "dgg_error.txt");
	computeStatics(dgg_error);

	printf("svg_statis:\n");
	vector<double> svg_error;
	computeError(correct_dis, svg_dis, svg_error);
	writeErrorFIle(svg_error, "svg_error.txt");
	computeStatics(svg_error);

	double min_error = 0;
	double max_error = 0.0001;
	normalizeError(dgg_error, min_error, max_error);
	outputDistanceField(model,dgg_error, obj_prefix + "dgg_error.obj");

	normalizeError(svg_error, min_error, max_error);
	outputDistanceField(model,svg_error, obj_prefix + "svg_error.obj");
}


void getObjDistance(CRichModel& model, string& fmm_filename , vector<double>& fmm_dis)
{
	FILE* file = fopen(fmm_filename.c_str(), "r");
	char buf[1024];
	
	while (fgets(buf, 1024, file) != NULL) {
		if (buf[0] == 'v' && buf[1] == 't') {
			double d;
			sscanf(buf + 3, "%lf", &d);
			fmm_dis.push_back(d);
		}
	}
	//assert(fmm_dis.size() == model.GetNumOfVerts());
	fclose(file);
}

void test_fast_marching(int source_vert)
{
	string obj_file_name = "gargoyle_nf700k.obj";
	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();
	vector<double> correct_dis;

	getICHDistance(model, source_vert, correct_dis);

	string fmm_filename = "gargoyle_nf700k_fmm_linear_dis.obj";
	vector<double> fmm_dis;
	getObjDistance(model, fmm_filename, fmm_dis);

	vector<double> fmm_error;
	computeError(correct_dis, fmm_dis, fmm_error);
	writeErrorFIle(fmm_error, "fmm_error.txt");
	computeStatics(fmm_error);

}

void normalizeDis(vector<double>& dis)
{
	double max_dis = *max_element(dis.begin(), dis.end());
	for (auto& d : dis){
		d /= max_dis;
	}
}

void normalizeDis(const vector<double>& dis, vector<double>& normalized_dis)
{
	normalized_dis.reserve(dis.size());
	double max_dis = *max_element(dis.begin(), dis.end());
	for (auto d : dis){
		d /= max_dis;
		normalized_dis.push_back(d);
	}
}


void test_heat(int source_vert)
{
	string obj_file_name = "gargoyle_nf700k.obj";
	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();
	vector<double> correct_dis;

	getICHDistance(model, source_vert, correct_dis);
	normalizeDis(correct_dis);

	string heat_filename = "gargoyle_nf700k_heat_s32042_m5.00.obj";
	vector<double> heat_dis;
	getObjDistance(model, heat_filename, heat_dis);
	normalizeDis(heat_dis);

	vector<double> heat_error;
	computeError(correct_dis, heat_dis, heat_error);
	writeErrorFIle(heat_error, "heat_error.txt");
	computeStatics(heat_error);

}

void test_delaunay_mesh()
{
	string obj_file_name = "gargoyle_nf700k.obj";
	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();
	int source_vert = 32042;

	string delaunay_obj_filename = "gargoyle_nf700k.delaunay5.obj";
	CRichModel delaunay_model(delaunay_obj_filename);
	delaunay_model.Preprocess();

	for (int i = 0; i < model.GetNumOfVerts(); ++i) {
		auto& p0 = model.Vert(i);
		auto& p1 = delaunay_model.Vert(i);
		if ((p0 - p1).Len() > 1e-6) {
			printf("i %d diff %lf\n", i , (p0 - p1).Len());
		}
	}

}

void test_heat_delaunay(int source_vert)
{
	string obj_file_name = "gargoyle_nf700k.obj";
	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();

	vector<double> correct_dis;
	getICHDistance(model, source_vert, correct_dis);
	normalizeDis(correct_dis);

	vector<string> heat_filenames = { "gargoyle_nf700k.delaunay2_heat_s32042_m3.00.obj",
									  "gargoyle_nf700k.delaunay3_heat_s32042_m3.00.obj",
									  "gargoyle_nf700k.delaunay4_heat_s32042_m3.00.obj" ,
									  "gargoyle_nf700k_heat_s32042_m5.00.obj"};
	vector<string> error_filenames = { "delaunay2_heat_error.txt",
									  "delaunay3_heat_error.txt",
							    	  "delaunay4_heat_error.txt",
									  "heat_erro.txt"};
	for (int type = 0; type < 4; ++type) 
	{
		//string heat_filename = "gargoyle_nf700k.delaunay4_heat_s32042_m3.00.obj";
		string heat_filename = heat_filenames[type];
		cout << heat_filename << "\n";
		vector<double> heat_dis;
		getObjDistance(model, heat_filename, heat_dis);
		//for (int i = 0; i < 100; ++i) {
		//	printf("%lf ", heat_dis[i]);
		//}
		//printf("\n");
		normalizeDis(heat_dis);
		heat_dis = vector<double>(heat_dis.begin(), heat_dis.begin() + model.GetNumOfVerts());

		vector<double> heat_error;
		computeError(correct_dis, heat_dis, heat_error);
		writeErrorFIle(heat_error, error_filenames[type]);
		computeStatics(heat_error);
	}
}

int main_old()
{
	int source_vert = 93598;
	string obj_filename = "gargoyle_nf700k.obj";
	string anisotropic_filename = "gargoyle_nf700k_anisotropic.obj";
	CRichModel model(obj_filename);
	model.Preprocess();
	CRichModel model_anisotropic(anisotropic_filename);
	model_anisotropic.Preprocess();
	double min_dis = 1e100;
	int min_pos = 0;
	for (int i = 0; i < model_anisotropic.GetNumOfVerts(); ++i) {
		double dis = (model_anisotropic.Vert(i) - model.Vert(source_vert)).Len();
		if (dis < min_dis) {
			min_dis = dis;
			min_pos = i;
		}
	}
	printf("min_pos %d min_dis %lf\n", min_pos, min_dis);//91011
	return 0;

}


int main(int argc, char** argv) {
	string input_file_name = "fertility_nf1k_anisotropic.obj";

	CRichModel model(input_file_name);
	model.Preprocess();
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
		return 0;
	}
	mesh.initialize_mesh_data(points, faces);		//create internal

	int sample_num = 200; 

	vector<CPoint3D> result_list;
	int source_index = atoi(argv[1]);
	result_list.push_back(model.Vert(source_index));
	geodesic::SurfacePoint source(&mesh.vertices()[source_index]);
	vector<geodesic::SurfacePoint> sources{ source };
	unsigned int seed = time(0);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> uniform_d(0, 1);//random noise len

	for (int i = 1; i < sample_num; ++i) {
		printf("i %d\n", i);
		geodesic::GeodesicAlgorithmBase *algorithm;

		algorithm = new geodesic::GeodesicAlgorithmExact(&mesh);

		algorithm->propagate(sources);

		int farthest_vert;
		double farthest_dis = 0;
		for (int j = 0; j < model.GetNumOfVerts(); ++j) {
			double dis;
			algorithm->best_source(geodesic::SurfacePoint(&mesh.vertices()[j]),dis);
			if (dis > farthest_dis) {
				farthest_dis = dis;
				farthest_vert = j;
			}
		}

		//neighbor faces
		int farthest_neighbor_v;
		int farthest_neighbor_e;
		double farthest_neighbor_dis = 0;
		for (auto n : model.Neigh(farthest_vert)) {
			int v = (model.Edge(n.first).indexOfRightVert);
			double dis;
			algorithm->best_source(geodesic::SurfacePoint(&mesh.vertices()[v]), dis);
			if (dis > farthest_neighbor_dis) {
				farthest_neighbor_dis = dis;
				farthest_neighbor_v = v;
				farthest_neighbor_e = n.first;
			}
		}

		int face_id;
		{
			int third_v0 = model.Edge(farthest_neighbor_e).indexOfOppositeVert;
			int third_v1 = model.Edge(model.Edge(farthest_neighbor_e).indexOfReverseEdge).indexOfOppositeVert;
			double third_dis0;
			double third_dis1;
			algorithm->best_source(geodesic::SurfacePoint(&mesh.vertices()[third_v0]), third_dis0);
			algorithm->best_source(geodesic::SurfacePoint(&mesh.vertices()[third_v1]), third_dis1);
			if (third_dis0 > third_dis1) {
				face_id = model.Edge(farthest_neighbor_e).indexOfFrontFace;
			}
			else{
				face_id = model.Edge(model.Edge(farthest_neighbor_e).indexOfReverseEdge).indexOfFrontFace;
			}
		}

		//random times = 100
		int random_times = 500;
		double max_dis = 0;
		geodesic::SurfacePoint max_surface_point;
		CPoint3D p;
		double max_a, max_b;
		for (int j = 0; j < random_times; ++j) {

			float alfa = (float)sqrt(uniform_d(gen));
			float beta = (float)uniform_d(gen);
			double a = 1.0 - alfa;
			double b = alfa * (1.0 - beta);
			if (a < 1e-6 || b < 1e-6) continue;
			geodesic::SurfacePoint dest(&mesh.faces()[face_id], a, b);
			double dis;
			algorithm->best_source(dest,dis);
			if (dis > max_dis) {
				max_a = a;
				max_b = b;
				max_dis = dis;
				max_surface_point = dest;
				p = model.Vert(model.Face(face_id)[0]) * a + model.Vert(model.Face(face_id)[1]) * b + model.Vert(model.Face(face_id)[2]) * (1 - a - b);
			}
		}
		printf("a %lf b %lf\n", max_a, max_b);
		//source = max_surface_point;
		sources.push_back(max_surface_point);
		result_list.push_back(p);
		//vector<geodesic::SurfacePoint> path;
		//algorithm->trace_back(dest, path);



		delete algorithm;
	}

	printBallToObj(result_list, "result_balls.obj", 0.01);
	
	return 0;

}


//sampling fast marching
int main_sampling_fmm(int argc, char** argv) {
	string input_file_name = "fertility_nf1k_anisotropic.obj";

	CRichModel model(input_file_name);
	model.Preprocess();
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
		return 0;
	}
	mesh.initialize_mesh_data(points, faces);		//create internal

	int sample_num = 200;

	vector<CPoint3D> result_list;
	int source_index = atoi(argv[1]);
	result_list.push_back(model.Vert(source_index));
	geodesic::SurfacePoint source(&mesh.vertices()[source_index]);
	vector<geodesic::SurfacePoint> sources{ source };
	unsigned int seed = time(0);

	for (int i = 1; i < sample_num; ++i) {
		printf("i %d\n", i);
		geodesic::GeodesicAlgorithmBase *algorithm;

		algorithm = new geodesic::GeodesicAlgorithmExact(&mesh);

		algorithm->propagate(sources);

		int farthest_vert;
		double farthest_dis = 0;
		for (int j = 0; j < model.GetNumOfVerts(); ++j) {
			double dis;
			algorithm->best_source(geodesic::SurfacePoint(&mesh.vertices()[j]), dis);
			if (dis > farthest_dis) {
				farthest_dis = dis;
				farthest_vert = j;
			}
		}

		geodesic::SurfacePoint max_surface_point(&mesh.vertices()[farthest_vert]);

		//source = max_surface_point;
		sources.push_back(max_surface_point);
		result_list.push_back(model.Vert(farthest_vert));
		//vector<geodesic::SurfacePoint> path;
		//algorithm->trace_back(dest, path);



		delete algorithm;
	}

	printBallToObj(result_list, "result_balls_s200_fmm.obj", 0.01);

	return 0;

}


int main_armardillo() {
	int source_vert = 89834;
	string obj_file_name = "armadillo_nf346k.obj";

	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();

	vector<double> correct_dis;

	getICHDistance(model, source_vert, correct_dis);
	outputDistanceField(model, correct_dis, obj_prefix + "_ich.obj");


}

int main_figure_backup()
{
	int source_vert = 93598;
	//int source_vert = 32042;
	string obj_file_name = "gargoyle_nf700k.obj";
	string dgg_filename = "gargoyle_nf700k_DGG0.000100_c20_pruning.binary";
	string svg_filename = "gargoyle_nf700k_SVG_k483.binary";

	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();
	double min_error = 0.00001;
	double max_error = 0.05;

	vector<double> correct_dis;
	getICHDistance(model, source_vert, correct_dis);
	outputDistanceField(model, correct_dis, obj_prefix + "_ich.obj");

	if (true) {
		vector<double> dgg_dis;
		getDGGDistance(model, source_vert, dgg_filename, dgg_dis);
		outputDistanceField(model, dgg_dis, obj_prefix + "_dgg_dis.obj");
		printf("dgg_statis:\n");
		vector<double> dgg_error;
		computeError(correct_dis, dgg_dis, dgg_error);
		writeErrorFIle(dgg_error, "dgg_error.txt");
		computeStatics(dgg_error);
		logNormalizeErrorAndOutput(dgg_error, min_error, max_error, model, obj_prefix + "dgg_log_error.obj");
	}

	if (true) {
		vector<double> svg_dis;
		getSVGDistance(model, source_vert, svg_filename, svg_dis);
		outputDistanceField(model, svg_dis, obj_prefix + "_svg_dis.obj");
		printf("svg_statis:\n");
		vector<double> svg_error;
		computeError(correct_dis, svg_dis, svg_error);
		writeErrorFIle(svg_error, "svg_error.txt");
		computeStatics(svg_error);
		//normalizeError(svg_error, min_error, max_error);
		logNormalizeErrorAndOutput(svg_error, min_error, max_error, model, obj_prefix + "svg_log_error.obj");
	}

	string fmm_filename = "gargoyle_nf700k_fmm_linear_dis.obj";
	vector<double> fmm_dis;
	getObjDistance(model, fmm_filename, fmm_dis);
	vector<double> fmm_error;
	computeError(correct_dis, fmm_dis, fmm_error);
	writeErrorFIle(fmm_error, "fmm_error.txt");
	computeStatics(fmm_error);
	logNormalizeErrorAndOutput(fmm_error, min_error, max_error, model, obj_prefix + "fmm_log_error.obj");



	string heat_filename = "gargoyle_heat.0.obj";
	vector<double> heat_dis;
	getObjDistance(model, heat_filename, heat_dis);
	vector<double> heat_error;
	computeError(correct_dis, heat_dis, heat_error);
	writeErrorFIle(heat_error, "heat_error.txt");
	computeStatics(heat_error);
	logNormalizeErrorAndOutput(heat_error, min_error, max_error, model, obj_prefix + "heat_log_error.obj");

	//vector<double> normalized_correct_dis;
	//normalizeDis(correct_dis, normalized_correct_dis);

	//vector<string> heat_filenames = { 
	//	"gargoyle_nf700k.delaunay4_heat_s" + to_string(source_vert) + "_m3.00.obj",
	//	"gargoyle_nf700k_heat_s" + to_string(source_vert) + "_m3.00.obj" };
	//vector<string> error_filename_prefix = { 
	//	"delaunay4_heat_error",
	//	"heat_erro" };
	//for (int type = 0; type < 2; ++type)
	//{
	//	//string heat_filename = "gargoyle_nf700k.delaunay4_heat_s32042_m3.00.obj";
	//	string heat_filename = heat_filenames[type];
	//	cout << heat_filename << "\n";
	//	vector<double> heat_dis;
	//	getObjDistance(model, heat_filename, heat_dis);
	//	normalizeDis(heat_dis);
	//	heat_dis = vector<double>(heat_dis.begin(), heat_dis.begin() + model.GetNumOfVerts());

	//	vector<double> heat_error;
	//	computeError(normalized_correct_dis, heat_dis, heat_error);
	//	//writeErrorFIle(heat_error, error_filename_prefix[type]+".txt");
	//	computeStatics(heat_error);
	//	logNormalizeErrorAndOutput(heat_error, min_error, max_error, model, error_filename_prefix[type] + ".obj");
	//}
	return 0;
}