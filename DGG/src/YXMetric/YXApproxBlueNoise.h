#ifndef __YX_APPROX_BLUE_NOISE_H__
#define __YX_APPROX_BLUE_NOISE_H__

#include "YXMetric.h"
#include "YXWhiteNoise.h"
#include "YXSimpleHashMap.h"

#ifndef UINT64_DEFINED
#define UINT64_DEFINED
typedef unsigned long long uint64;
#endif


struct YXBlueNoisePoint{
	uint64 cellID;
	int status;
	int index;
	int faceid;
	YXPoint3D pt;
	bool operator<(YXBlueNoisePoint & other){
		return this->cellID < other.cellID;
	}
};

class YXApproxBlueNoise{
	enum{IDLE,REJECTED,ACCEPTED};
public:
	/*
	 * approxSamplingNumber : approximate number of blue noise samples.
	 * preSampleRate : rate of white noise.  the presampled number of white noise = preSampleRate * approxSamplingNumber
	 * the results of sampling are stored in vector<> result.
	 * return : the minimum distance between two samples. (approximate)
	*/
	float radius;
	std::vector <YXMeshPoint> result;
	YXMesh3D mesh;

	void load(const char * filename) {
		string filename_str(filename);
		int found = filename_str.rfind(".");
		string post_fix = filename_str.substr(found+1);
		cout << post_fix << "\n";
		if (post_fix == "m") {
			mesh.load_m(filename);
		}
		else if (post_fix == "obj"){
			mesh.load_obj(filename);
		}
		else{
			printf("filetype %s can not recognize!\n", filename_str.c_str());
			exit(1);
		}

		metric.buildFromMesh(mesh);
		metric.assignNormalsForMesh(mesh);
		xmax = ymax = zmax = 1e-99;
		xmin = ymin = zmin = 1e99;
		for(int i = 0; i < (int)mesh.vert.size(); ++i){
			xmax = YX_MAX(xmax, mesh.vert[i].x);
			ymax = YX_MAX(ymax, mesh.vert[i].y);
			zmax = YX_MAX(zmax, mesh.vert[i].z);
			xmin = YX_MIN(xmin, mesh.vert[i].x);
			ymin = YX_MIN(ymin, mesh.vert[i].y);
			zmin = YX_MIN(zmin, mesh.vert[i].z);
		}
	}
	float sampling(int approxSamplingNumber, int preSampleRate = 10){
		YXWhiteNoise whiteNoise;
		int WN = preSampleRate * approxSamplingNumber;
		whiteNoise.sampling(&metric, WN);
		result.clear();
		// 0.74 is the radius statistics 
		radius = 0.74 * sqrt(metric.totalArea / 2.0 / sqrt(3.0) / approxSamplingNumber) * 2;
		computeGridSize(); 
		assignCells(&whiteNoise, WN);
		result.clear();
		randomOrderSampling(&whiteNoise, WN);
		whiteNoise.clear();
		return radius;
	}
protected:
	YXMetric metric;
	std::vector <YXBlueNoisePoint> ptInCells;
	int gridSizeX, gridSizeY, gridSizeZ;
	float gridXmin, gridYmin, gridZmin;
	float xmax, ymax, zmax;
	float xmin, ymin, zmin;
	YXSimpleHashMap <uint64, int> beginIndex, endIndex;
	uint64 position_to_cellID(YXPoint3D & pt){
		int x = int((pt.x - gridXmin) / radius);
		int y = int((pt.y - gridYmin) / radius);
		int z = int((pt.z - gridZmin) / radius);
		return (uint64)((uint64)x * gridSizeY + y) * gridSizeZ + z;
	};
	void cellID_to_xyz(uint64 id, int &x, int &y, int &z){
		z = int(id % gridSizeZ);
		id /= uint64(gridSizeZ);
		y = int(id % gridSizeY);
		x = int(id / gridSizeY);
	}
	void computeGridSize(){
		gridSizeX = int ((xmax - xmin) / radius) + 1;
		gridSizeY = int ((ymax - ymin) / radius) + 1;
		gridSizeZ = int ((zmax - zmin) / radius) + 1;
		gridXmin = xmin;
		gridYmin = ymin;
		gridZmin = zmin;
	}
	void assignCells(YXWhiteNoise * wn, int N){
		YXMeshPoint * pt = wn->pt;
		ptInCells.resize(N);
		for(int i = 0; i < N; ++i){
			const YXFace & f = mesh.face[pt[i].faceid];
			ptInCells[i].pt = mesh.vert[f.v[0]].scale( pt[i].a ) + 
				mesh.vert[f.v[1]].scale( pt[i].b ) + mesh.vert[f.v[2]].scale( 1.0 - pt[i].a - pt[i].b );
			ptInCells[i].cellID = position_to_cellID(ptInCells[i].pt);
			ptInCells[i].status = IDLE;
			ptInCells[i].index = i;
			ptInCells[i].faceid = pt[i].faceid;
		}
		std::sort(ptInCells.begin(), ptInCells.end());
		int hashsize = N/4;
		if(hashsize > 100000) hashsize = 100001;
		beginIndex.init( hashsize );
		endIndex.init( hashsize );
		beginIndex.set(ptInCells[0].cellID, 0);
		for(int i = 1; i < N; ++i){
			if(ptInCells[i].cellID != ptInCells[i-1].cellID) beginIndex.set(ptInCells[i].cellID, i);
		}
		for(int i = 0; i < N-1; ++i){
			if(ptInCells[i].cellID != ptInCells[i+1].cellID) endIndex.set(ptInCells[i].cellID, i);
		}
		endIndex.set(ptInCells[N-1].cellID, N-1);
	}

	void randomOrderSampling(YXWhiteNoise * wn, int N){
		std::vector <int> randomIndex;
		randomIndex.resize(N);
		for(int i = 0; i < N; ++i) randomIndex[i] = i;
		std::random_shuffle(randomIndex.begin(), randomIndex.end());
		for(int i = 0; i < N; ++i){
			int idx = randomIndex[i];
			if(ptInCells[idx].status != IDLE) continue;
			rejectNeighbors(idx);
			result.push_back( wn->pt[ ptInCells[idx].index ] );
		}
	}
	void rejectNeighbors(int idx){
		int x, y, z;
		cellID_to_xyz(ptInCells[idx].cellID, x, y, z);
		for(int dx = x-1; dx <= x+1; ++dx)
			for(int dy = y-1; dy <= y+1; ++dy)
				for(int dz = z-1; dz <= z+1; ++dz){
					if(dx < 0 || dy < 0 || dz < 0 || dx >= gridSizeX || dy >= gridSizeY || dz >= gridSizeZ)continue;
					uint64 neighborCell = (uint64)((uint64)dx * gridSizeY + dy) * gridSizeZ + dz;
					int * begin = beginIndex.get(neighborCell);
					if(begin == NULL)continue;
					int * end = endIndex.get(neighborCell);
					for(int i = *begin; i <= *end; ++i){
						if(ptInCells[i].status != IDLE) continue;
						if(YXGeometry::approximateGeodesicDistance(ptInCells[i].pt, ptInCells[idx].pt, 
							mesh.faceNormal[ ptInCells[i].faceid ], mesh.faceNormal[ ptInCells[idx].faceid ]) < radius){
								ptInCells[i].status = REJECTED;
						}
					}
				}
	}
};


#endif // __YX_APPROX_BLUE_NOISE_H__
