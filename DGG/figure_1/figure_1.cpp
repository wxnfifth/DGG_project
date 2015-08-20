#include "ICH\RichModel.h"
#include "ICH\ICHWithFurtherPriorityQueue.h"
#include "wxn\wxn_dijstra.h"
#include "wxn\wxnTime.h"

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
	s_graph->findShortestDistance(source_vert);
	
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
	s_graph->findShortestDistance(source_vert);

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


void outputDistanceField(CRichModel& model, vector<double>& dis, const string& output_filename)
{
	assert(model.GetNumOfVerts() == dis.size());

	vector<pair<double, double>> ich_texture;
	double max_dis = std::max(max_dis, *std::max_element(dis.begin(), dis.end()));
	ich_texture.reserve(model.GetNumOfVerts());
	for (double d : dis) {
		ich_texture.push_back(make_pair(d / max_dis, d / max_dis));
	}
	model.FastSaveObjFile(output_filename, ich_texture);

}

void computeError(const vector<double>& correct_dis, const vector<double>& dis, vector<double>& error)
{
	assert(correct_dis.size() == dis.size());
	error.resize(dis.size());
	for (int i = 0; i < correct_dis.size(); ++i) {
		if (fabs(correct_dis[i]) < 1e-7) {
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
	for (auto& e : svg_error) {
		e = (e - min_error) / (max_error - min_error);
	}
}


void figure_1()
{
	string obj_file_name = "gargoyle_nf700k.obj";
	string dgg_filename = "gargoyle_nf700k_DGG0.000100_c20_pruning.binary";
	string svg_filename = "gargoyle_nf700k_SVG_k483.binary";

	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();

	int source_vert = 319073;

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

int main()
{

  figure_1();



}