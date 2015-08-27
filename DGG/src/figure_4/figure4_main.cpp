#include "wxn\wxn_dijstra.h"
#include "svg_definition.h"
#include "wxn\wxn_path_helper.h"

void figure_4()
{
  string model_name = "bunny_nf10k.obj";
  string svg_filename = "bunny_nf10k_DGG0.001000_c12_pruning.binary";
  //string svg_filename = "bunny_nf10k_DGG0.003162_c7_pruning.binary";
  int source_vertex = 1397;
  CRichModel model(model_name);
  model.Preprocess();
  SparseGraph<float>* s_graph = NULL;
  s_graph = new LC_HY<float>();
  s_graph->read_svg_file_with_angle((string)svg_filename);
  dynamic_cast<LC_HY<float>*>(s_graph)->setModel(model);



  s_graph->findShortestDistance(source_vertex);  //find the geodesic distance from a 
  
  vector<int> sources;
  sources.push_back(source_vertex);
  CICHWithFurtherPriorityQueue ich_algoritm(model,sources);
  ich_algoritm.Execute();
  vector<int> result_path_nodes;
  int dest = 0;
  for (int i = 0; i < model.GetNumOfVerts(); ++i) {
    vector<int> path_nodes;
    s_graph->getPath(i,path_nodes);
    if (path_nodes.size() == 5) {
      //printf("path
      double exact_dis = ich_algoritm.m_InfoAtVertices[i].disUptodate;
      double dgg_dis = s_graph->distanceToSource(i);
      double error = fabs(dgg_dis - exact_dis) / exact_dis;
      printf("dis_exact %lf dis_dgg %lf error %lf\n" , exact_dis, dgg_dis, error);
      dest = i;
      result_path_nodes = path_nodes;
      if (exact_dis > 0.5  ) break;
     // break;
    }
  }
  CylinderPath path_approximate(0.0005);
  path_approximate.addGeodesicPaths(model,result_path_nodes);
  path_approximate.write_to_file("bunny_DGG_path_0.001.obj");

  vector<CPoint3D> pts;
  for (int i = 1; i < result_path_nodes.size() -1; ++i) {
    pts.push_back(model.Vert(result_path_nodes[i]));
  }
  printBallToObj(pts,"bunny_dgg_path_pts.obj",0.0015);
  vector<CPoint3D> source_dest;
  source_dest.push_back(model.Vert(source_vertex));
  source_dest.push_back(model.Vert(dest));
  printBallToObj(source_dest,"bunny_dgg_source_pts.obj",0.0015);
  
  CylinderPath path_exact(0.0005);
  path_exact.addGeodesicPath(model,source_vertex,dest);
  path_exact.write_to_file("bunny_exact_path_0.001.obj");
   
  //
  double exact_dis = ich_algoritm.m_InfoAtVertices[dest].disUptodate;
  double approximate_dis = s_graph->distanceToSource(dest);
  double error = fabs(approximate_dis - exact_dis) / exact_dis;
  printf("eps 0.001, source %d dest %d dis %lf error %lf\n" , 
    source_vertex, dest, exact_dis,  error);


}


void dump_edge_lens(const string& model_name, const string& svg_filename, const string& dgg_filename)
{
	//string model_name = "bunny_nf10k.obj";
	//string svg_filename = "bunny_nf10k_DGG0.001000_c12_pruning.binary";
	CRichModel model(model_name);
	model.Preprocess();
	string prefix = model_name.substr(0, model_name.length() - 4);
	{
		SparseGraph<float>* s_graph = NULL;
		s_graph = new Dijstra_vector<float>();
		s_graph->read_svg_file_binary((string)svg_filename);

		string svg_edge_filename = prefix + "_svg_edge_lens.txt";
		FILE* svg_file = fopen(svg_edge_filename.c_str(), "w");
		for (int i = 0; i < s_graph->NodeNum(); ++i) {
			auto& neigh_dis = s_graph->graphNeighborDis(i);
			for (auto& d : neigh_dis) {
				fprintf(svg_file, "%lf\n", d);
			}
		}
		fclose(svg_file);
		delete s_graph;
	}
	{
		SparseGraph<float>* s_graph = NULL;
		s_graph = new LC_HY<float>();
		s_graph->read_svg_file_with_angle((string)dgg_filename);
		dynamic_cast<LC_HY<float>*>(s_graph)->setModel(model);
		string dgg_edge_filename = prefix + "_dgg_edge_lens.txt";
		FILE* dgg_file = fopen(dgg_edge_filename.c_str(), "w");
		for (int i = 0; i < s_graph->NodeNum(); ++i) {
			auto& neigh_dis = s_graph->graphNeighborDis(i);
			for (auto& d : neigh_dis) {
				fprintf(dgg_file, "%lf\n", d);
			}
		}
		fclose(dgg_file);
		delete s_graph;
	}





}

int main(int argc, char** argv)
{
	//figure_4();

	if (argc < 3) {
		printf("error input, example xx.exe xx.obj svg.binary dgg.binary");
		exit(1);
	}

	dump_edge_lens(argv[1],argv[2],argv[3]);


  return 0;
}