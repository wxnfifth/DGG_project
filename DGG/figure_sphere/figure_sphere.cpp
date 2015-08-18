#include "ICH\RichModel.h"
#include <random>


void figure_sphere(string sphere_filename, double noise_percent)
{

	CRichModel model(sphere_filename);
	model.Preprocess();
	string prefix = sphere_filename.substr(0, sphere_filename.length() - 4);
	double noise_len = model.GetMeanEdgeLength() * 0.1;
	vector<int> vert_index(model.GetNumOfVerts());
	for (int i = 0; i < vert_index.size(); ++i) {
		vert_index[i] = i;
	}


	random_shuffle(vert_index.begin(), vert_index.end());
	int noise_vert_num = noise_percent * model.GetNumOfVerts();

	printf("noise_vert_num %d\n", noise_vert_num);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<double> d(1, 1);
	for (int i = 0; i < noise_vert_num; ++i) {
		int vert = vert_index[i];
		CPoint3D direct(d(gen), d(gen), d(gen));
		direct.Normalize();	 direct *= noise_len;
		model.m_Verts[vert] += direct;
	}
	char output_filename[256];
	sprintf(output_filename, "%s_noise%.2lf.obj", prefix.c_str(), noise_percent);
	model.FastSaveObjFile(output_filename);
}

int main(int argc, char** argv)
{
	if (argc < 3) {
		fprintf(stderr, "insufficient arguments!\n");
		fprintf(stderr, "usage xx.exe sphere.obj 0.01\n");
		exit(1);
	}
	string obj_name = argv[1];
	double noise_percent = atof(argv[2]);
	figure_sphere(obj_name, noise_percent);

}