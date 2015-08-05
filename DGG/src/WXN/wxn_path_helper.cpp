
#include "wxn_path_helper.h"
#include "ICH\RichModel.h"
#include <vector>
using std::vector;


void cylinder(CPoint3D start_p, CPoint3D end_p) 
{
	double radius = 0.015;
	vector<CPoint3D> up_points;
	vector<CPoint3D> down_points;

	CPoint3D axis = (end_p - start_p).Normalize();
	CPoint3D p0_vec = (axis * CPoint3D(0,0,1)).Normalize();
	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	//p0_vec = radius * p0_vec;
	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	CPoint3D p4_vec = -1.0 * p0_vec;
	CPoint3D p6_vec = (axis * p0_vec).Normalize();
	CPoint3D p2_vec = -1.0 * p6_vec;
	CPoint3D p7_vec = (p0_vec + p6_vec).Normalize();
	CPoint3D p1_vec = (p0_vec + p2_vec).Normalize();
	CPoint3D p3_vec = (p2_vec + p4_vec).Normalize();
	CPoint3D p5_vec = (p4_vec + p6_vec).Normalize();

	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	//printf("p7 len %lf p0 len %lf\n" , p7_vec.Len(), p0_vec.Len());

	up_points.push_back(p0_vec);
	up_points.push_back(p1_vec);
	up_points.push_back(p2_vec );
	up_points.push_back(p3_vec );
	up_points.push_back(p4_vec );
	up_points.push_back(p5_vec );
	up_points.push_back(p6_vec );
	up_points.push_back(p7_vec );

	vector<CPoint3D> up_points_2;
	for (int i = 0; i < up_points.size(); ++i) {
		up_points_2.push_back(up_points[i]);
		CPoint3D p = (up_points[i] + up_points[(i+1)%up_points.size()]).Normalize();
		up_points_2.push_back(p);
	}
	swap(up_points, up_points_2);

	for (auto& p:up_points) {
		p *= radius;
		p += start_p;
	}
  FILE* file_tmp;
	fopen_s(&file_tmp, "up_disk.obj","w");
	fprintf(file_tmp , "v %lf %lf %lf\n" , start_p.x, start_p.y, start_p.z);
	for (int i = 0; i < up_points.size(); ++i) {
		fprintf(file_tmp, "v %lf %lf %lf\n" , up_points[i].x, up_points[i].y, up_points[i].z);
	}
	for (int i = 0; i < up_points.size(); ++i) {
		fprintf(file_tmp, "f 1 %d %d\n" , i + 2, (i+1) % up_points.size() + 2);
	}

	fclose(file_tmp);


	for (auto& p:up_points) {
		down_points.push_back(p + end_p - start_p);
	}

	FILE* file_cylinder(NULL);
  fopen_s(&file_cylinder,"cylinder.obj", "w");
	for (auto& p:up_points) {
		fprintf(file_cylinder, "v %lf %lf %lf\n" , p.x, p.y, p.z);
	}
	for (auto& p:down_points) {
		fprintf(file_cylinder, "v %lf %lf %lf\n" , p.x, p.y, p.z);
	}
	for (int i = 0; i < up_points.size(); ++i) {
		int up_v0 = i + 1;
		int up_v1 = (i+1)% up_points.size() + 1;
		int down_v0 = i + up_points.size() + 1;
		int down_v1 = (i+1)%up_points.size() + up_points.size() + 1;
		fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v1, up_v1);
		fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v0, down_v1);
	}
	fclose(file_cylinder);
}

void generateCylinder(CPoint3D start_p, CPoint3D end_p, vector<CPoint3D>& verts, vector<CBaseModel::CFace>& faces, double radius) 
{
	vector<CPoint3D> up_points;
	vector<CPoint3D> down_points;

	CPoint3D axis = (end_p - start_p).Normalize();
	CPoint3D p0_vec;
	if (axis.equal(CPoint3D(1,0,0))) {
		p0_vec = (axis * CPoint3D(0,0,1)).Normalize();
	}else{
		p0_vec = (axis * CPoint3D(1,0,0)).Normalize();
	}
	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	//p0_vec = radius * p0_vec;
	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	CPoint3D p4_vec = -1.0 * p0_vec;
	CPoint3D p6_vec = (axis * p0_vec).Normalize();
	CPoint3D p2_vec = -1.0 * p6_vec;
	CPoint3D p7_vec = (p0_vec + p6_vec).Normalize();
	CPoint3D p1_vec = (p0_vec + p2_vec).Normalize();
	CPoint3D p3_vec = (p2_vec + p4_vec).Normalize();
	CPoint3D p5_vec = (p4_vec + p6_vec).Normalize();

	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	//printf("p7 len %lf p0 len %lf\n" , p7_vec.Len(), p0_vec.Len());

	up_points.push_back(p0_vec);
	up_points.push_back(p1_vec);
	up_points.push_back(p2_vec );
	up_points.push_back(p3_vec );
	up_points.push_back(p4_vec );
	up_points.push_back(p5_vec );
	up_points.push_back(p6_vec );
	up_points.push_back(p7_vec );

	vector<CPoint3D> up_points_2;
	for (int i = 0; i < up_points.size(); ++i) {
		up_points_2.push_back(up_points[i]);
		CPoint3D p = (up_points[i] + up_points[(i+1)%up_points.size()]).Normalize();
		up_points_2.push_back(p);
	}
	swap(up_points, up_points_2);

	for (auto& p:up_points) {
		p *= radius;
		p += start_p;
	}

	for (auto& p:up_points) {
		down_points.push_back(p + end_p - start_p);
	}

	//FILE* file_cylinder = fopen("cylinder.obj", "w");
	//for (auto& p:up_points) {
	//	fprintf(file_cylinder, "v %lf %lf %lf\n" , p.x, p.y, p.z);
	//}
	//for (auto& p:down_points) {
	//	fprintf(file_cylinder, "v %lf %lf %lf\n" , p.x, p.y, p.z);
	//}
	//for (int i = 0; i < up_points.size(); ++i) {
	//	int up_v0 = i + 1;
	//	int up_v1 = (i+1)% up_points.size() + 1;
	//	int down_v0 = i + up_points.size() + 1;
	//	int down_v1 = (i+1)%up_points.size() + up_points.size() + 1;
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v1, up_v1);
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v0, down_v1);
	//}
	//fclose(file_cylinder);
	int verts_old_sz = verts.size();
	for (auto& p:up_points) {
		verts.push_back(p);
	}
	for (auto& p:down_points) {
		verts.push_back(p);
	}
	for (int i = 0; i < up_points.size(); ++i) {
		int up_v0 = verts_old_sz + i + 1;
		int up_v1 = verts_old_sz + (i+1) % up_points.size() + 1;
		int down_v0 = verts_old_sz + i + up_points.size() + 1;
		int down_v1 = verts_old_sz + (i+1) %up_points.size() + up_points.size() + 1;
		faces.push_back(CBaseModel::CFace(up_v0, down_v1, up_v1));
		faces.push_back(CBaseModel::CFace(up_v0, down_v0, down_v1));
	}

	//for (int i = 0; i < up_points.size(); ++i) {
	//	int up_v0 = i + 1;
	//	int up_v1 = (i+1)% up_points.size() + 1;
	//	int down_v0 = i + up_points.size() + 1;
	//	int down_v1 = (i+1)%up_points.size() + up_points.size() + 1;
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v1, up_v1);
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v0, down_v1);
	//}



}

void generateArrow(CPoint3D start_p, CPoint3D end_p, vector<CPoint3D>& verts, vector<CBaseModel::CFace>& faces, double radius) 
{
	vector<CPoint3D> up_points;
	vector<CPoint3D> down_points;

	CPoint3D axis = (end_p - start_p).Normalize();
	CPoint3D p0_vec;
	if (axis.equal(CPoint3D(1,0,0))) {
		p0_vec = (axis * CPoint3D(0,0,1)).Normalize();
	}else{
		p0_vec = (axis * CPoint3D(1,0,0)).Normalize();
	}
	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	//p0_vec = radius * p0_vec;
	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	CPoint3D p4_vec = -1.0 * p0_vec;
	CPoint3D p6_vec = (axis * p0_vec).Normalize();
	CPoint3D p2_vec = -1.0 * p6_vec;
	CPoint3D p7_vec = (p0_vec + p6_vec).Normalize();
	CPoint3D p1_vec = (p0_vec + p2_vec).Normalize();
	CPoint3D p3_vec = (p2_vec + p4_vec).Normalize();
	CPoint3D p5_vec = (p4_vec + p6_vec).Normalize();

	//printf("p0 %lf %lf %lf\n" , p0_vec.x, p0_vec.y, p0_vec.z);
	//printf("p7 len %lf p0 len %lf\n" , p7_vec.Len(), p0_vec.Len());

	up_points.push_back(p0_vec);
	up_points.push_back(p1_vec);
	up_points.push_back(p2_vec );
	up_points.push_back(p3_vec );
	up_points.push_back(p4_vec );
	up_points.push_back(p5_vec );
	up_points.push_back(p6_vec );
	up_points.push_back(p7_vec );

	vector<CPoint3D> up_points_2;
	for (int i = 0; i < up_points.size(); ++i) {
		up_points_2.push_back(up_points[i]);
		CPoint3D p = (up_points[i] + up_points[(i+1)%up_points.size()]).Normalize();
		up_points_2.push_back(p);
	}
	swap(up_points, up_points_2);

	for (auto& p:up_points) {
		p *= radius;
		p += start_p;
	}

	for (auto& p:up_points) {
		down_points.push_back(end_p);
	}

	//FILE* file_cylinder = fopen("cylinder.obj", "w");
	//for (auto& p:up_points) {
	//	fprintf(file_cylinder, "v %lf %lf %lf\n" , p.x, p.y, p.z);
	//}
	//for (auto& p:down_points) {
	//	fprintf(file_cylinder, "v %lf %lf %lf\n" , p.x, p.y, p.z);
	//}
	//for (int i = 0; i < up_points.size(); ++i) {
	//	int up_v0 = i + 1;
	//	int up_v1 = (i+1)% up_points.size() + 1;
	//	int down_v0 = i + up_points.size() + 1;
	//	int down_v1 = (i+1)%up_points.size() + up_points.size() + 1;
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v1, up_v1);
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v0, down_v1);
	//}
	//fclose(file_cylinder);
	int verts_old_sz = verts.size();
	for (auto& p:up_points) {
		verts.push_back(p);
	}
	for (auto& p:down_points) {
		verts.push_back(p);
	}
	for (int i = 0; i < up_points.size(); ++i) {
		int up_v0 = verts_old_sz + i + 1;
		int up_v1 = verts_old_sz + (i+1) % up_points.size() + 1;
		int down_v0 = verts_old_sz + i + up_points.size() + 1;
		int down_v1 = verts_old_sz + (i+1) %up_points.size() + up_points.size() + 1;
		faces.push_back(CBaseModel::CFace(up_v0, down_v1, up_v1));
		faces.push_back(CBaseModel::CFace(up_v0, down_v0, down_v1));
	}

	//for (int i = 0; i < up_points.size(); ++i) {
	//	int up_v0 = i + 1;
	//	int up_v1 = (i+1)% up_points.size() + 1;
	//	int down_v0 = i + up_points.size() + 1;
	//	int down_v1 = (i+1)%up_points.size() + up_points.size() + 1;
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v1, up_v1);
	//	fprintf(file_cylinder, "f %d %d %d\n" , up_v0, down_v0, down_v1);
	//}



}

void printBallToObj(const vector<CPoint3D>& vertex_list, const string& file_name, double scale)
{
	CBaseModel sphere_model("c:\\util\\sphere.obj");

	sphere_model.LoadModel();
	FILE* ball_file(NULL);
  fopen_s(&ball_file,file_name.c_str(), "w");
	int cnt = 1;
	for(int i = 0; i < vertex_list.size();++i) {
		CBaseModel temp_model(sphere_model);
		temp_model.Scale(scale);
		temp_model.Translate(vertex_list[i]);
		fprintf(ball_file, "g ball_%d\n", i % 2);
		for (int vert_id = 0; vert_id < temp_model.GetNumOfVerts(); ++vert_id) {
			fprintf(ball_file, "v %lf %lf %lf\n", temp_model.Vert(vert_id).x, temp_model.Vert(vert_id).y, temp_model.Vert(vert_id).z);
		}
		for(int face_id = 0; face_id < temp_model.GetNumOfFaces(); ++face_id) {
			fprintf(ball_file, "f %d %d %d\n", temp_model.Face(face_id)[0] + cnt, temp_model.Face(face_id)[1] + cnt, temp_model.Face(face_id)[2] + cnt);
		}
		cnt += temp_model.GetNumOfVerts();
	}
	fclose(ball_file);

}


void output_faces(const CRichModel* model_ptr, const set<int>& faces, const string& output_filename)
{
	FILE* output_file;
  fopen_s(&output_file,output_filename.c_str(), "w");
	if (output_file == nullptr) {
		printf("%s not found!\n");
		exit(1);
	}
	int cnt = 1;
	for (auto f:faces) {
		for (int i = 0; i < 3; ++i) {
			auto& p = model_ptr->Vert(model_ptr->Face(f)[i]);
			fprintf(output_file, "v %lf %lf %lf\n" ,  p.x, p.y, p.z);
		}
		fprintf(output_file, "f %d %d %d\n" , cnt,cnt+1,cnt+2);
		cnt += 3;
	}
	fclose(output_file);
}

void output_cylinder(const string& filename, const vector<CPoint3D>& verts, const vector<CBaseModel::CFace>& faces)
{
	FILE* cylinder_file;
  fopen_s(&cylinder_file,filename.c_str(), "w");
	for (auto& p:verts) {
		fprintf(cylinder_file, "v %lf %lf %lf\n" , p.x, p.y, p.z);
	}
	for (auto& f:faces) {
		fprintf(cylinder_file, "f %d %d %d\n" , f[0] , f[1], f[2]);
	}
	fclose(cylinder_file);
}

