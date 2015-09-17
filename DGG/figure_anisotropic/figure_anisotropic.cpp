#include "ich\RichModel.h"
#include "wxn\wxnMath.h"

void print_anisotropic(string& model_name)
{
	CRichModel model(model_name);
	model.Preprocess();
	double ave_anisotropic = 0;
	for (int i = 0; i < model.GetNumOfFaces(); ++i) {
		auto& f = model.Face(i);
		double len[3];
		for (int j = 0; j < 3; ++j) {
			len[j] = (model.Vert(f[j]) - model.Vert(f[(j + 1) % 3])).Len();
		}
		double quality = calculateTriangleQuality(len[0], len[1], len[2]);
		ave_anisotropic += 1 / quality;
	}
	ave_anisotropic /= model.GetNumOfFaces();
	printf("average anisotropic = %lf\n", ave_anisotropic);
}

int main(int argc, char** argv)
{
	if (argc < 2) {
		printf( "error! example xx.exe xx.obj\n" );
		exit(1);
	}
	string model_name = argv[1];
	print_anisotropic(model_name);





	return 0;
}