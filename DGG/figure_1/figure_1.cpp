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
		if (isfinite(d)) {
			correct_dis[i] = d;
		}
		else{
			correct_dis[i] = 0;
		}
	}

}


void figure_1()
{
	string obj_file_name = "gargoyle_nf700k.obj";
	string dgg_filename = "gargoyle_nf700k_DGG0.000100_c20_pruning.binary";
	
	string obj_prefix = obj_file_name.substr(0, obj_file_name.length() - 4);
	CRichModel model(obj_file_name);
	model.Preprocess();
	SparseGraph<float>* s_graph = NULL;
	s_graph = new LC_HY<float>();
	s_graph->read_svg_file_with_angle((string)dgg_filename);
	dynamic_cast<LC_HY<float>*>(s_graph)->setModel(model);

	int source_vert = 319073;

	vector<double> correct_dis;

	getICHDistance(model, source_vert, correct_dis);

	double max_dis = std::max(max_dis, *std::max_element(correct_dis.begin(), correct_dis.end()));

	vector<pair<double, double>> ich_texture;
	ich_texture.reserve(model.GetNumOfVerts());
	for (double d : correct_dis) {
		ich_texture.push_back(make_pair(d / max_dis, d / max_dis));
	}
	model.FastSaveObjFile(obj_prefix + "_ich.obj", ich_texture);

}

int main()
{





}