#ifndef __YX_WHITE_NOISE_H__
#define __YX_WHITE_NOISE_H__

#include "YXMetric.h"

struct YXMeshPoint{
	float a;
	float b;
	int faceid;
};

class YXWhiteNoise{
public:
	YXMeshPoint * pt;
	int * ptCountOnFace;
	YXWhiteNoise() : pt(NULL), ptCountOnFace(NULL){};
	void sampling(const YXMetric * metric, int number){
		pt = new YXMeshPoint[number];
		ptCountOnFace = new int [metric->nface + 1];
		float partArea = 0.0;
		ptCountOnFace[0] = 0;
		unsigned int seed = (unsigned int)time(0);
		//seed = 100;
		for(int i = 0; i < metric->nface; ++i){
			partArea += metric->area[i];
			ptCountOnFace[i+1] = int( (partArea / metric->totalArea * number) + 0.5 );
			for(int j = ptCountOnFace[i]; j < ptCountOnFace[i+1]; ++j){
				float alfa = (float)sqrt(YXGeometry::randomfloat(seed));
				float beta = (float)YXGeometry::randomfloat(seed) ;
				pt[j].a = 1.0 - alfa;
				pt[j].b = alfa * (1.0 - beta);
				pt[j].faceid = i;
			}
		}
	}
	YXMeshPoint * getPointOnFace(int faceIndex) const{
		return pt + ptCountOnFace[faceIndex];
	}
	int size(int faceIndex) const {
		return ptCountOnFace[faceIndex+1] - ptCountOnFace[faceIndex];
	}
	void clear(){
		if(pt) delete [] pt;
		if(ptCountOnFace) delete [] ptCountOnFace;
		pt = NULL;
		ptCountOnFace = NULL;
	}
};

/*********
//test code : you can run this code to generate a white noise distribution
void testWhiteNoise(){
	YXTimer t;
	t.start();
	YXMesh3D mesh;
	mesh.load( "D:\\models\\bunny_nf144k.m");
	// or :
	//mesh.pushVertex( x, y, z );
	//mesh.pushFace(v1, v2, v3); // start from 0
	YXMetric metric;
	metric.buildFromMesh( mesh );
	printf("Preprocess time: %lf seconds\n", t.gettime());

	t.start();
	YXWhiteNoise whiteNoiseSampling;
	whiteNoiseSampling.sampling( &metric, 1000000 ); // metric and number of samples
	printf("Sampling time : %lf seconds\n", t.gettime());

	// get samples on face face_id:
	int face_id = 99;
	YXMeshPoint * pt = whiteNoiseSampling.getPointOnFace( face_id );
	for(int i = 0; i < whiteNoiseSampling.size( face_id ); ++i){
		printf("[Face %d]Sample #%d : {face_id=%d, (a,b)=(%f,%f)}\n", face_id, i+1, pt[i].faceid, pt[i].a, pt[i].b);
	}

	whiteNoiseSampling.clear();
}
*********/

#endif //__YX_WHITE_NOISE_H__