 
#include "stdafx.h"
#include "YXPathTracer.h"
#include "YXMetric.h"
#include <stdio.h>
using std::vector;

template <int N>
struct NVector{
	double d[N];
	double & operator[](int idx){
		return d[idx];
	}
	NVector operator+(const NVector & o) const {
		NVector tmp;
		for(int i = 0; i < N; ++i) tmp[i] = d[i] + o.d[i];
		return tmp;
	};
	NVector operator-(const NVector & o) const {
		NVector tmp;
		for(int i = 0; i < N; ++i) tmp[i] = d[i] - o.d[i];
		return tmp;
	};
	double length() const {
		double res = 0.0;
		for(int i = 0; i < N; ++i) res += d[i] * d[i];
		return sqrt(res);
	}
	double length2() const {
		double res = 0.0;
		for(int i = 0; i < N; ++i) res += d[i] * d[i];
		return res;
	}
	double shiftDistance(const NVector & o) const{
		double dis = 1e99;
		for(int s = 0; s < N; ++s){
			double tmpdis = 0.0;
			for(int i = 0; i < N; ++i){
				tmpdis += (d[i] - o.d[i])*(d[i] - o.d[i]);
			}
			if(tmpdis < dis) dis = tmpdis;
		}
		return dis;
	}
	double dotProduct(const NVector & o) const {
		double res = 0.0;
		for(int i = 0; i < N; ++i) res += d[i] * o.d[i];
		return res;
	}
	NVector shiftToLowest() const {
		const double * p = std::min_element(d, d+N);
		int shift = (int)(p - &d[0]);
		NVector tmp;
		for(int i = 0; i < N; ++i){
			tmp[i] = d[(i+shift)%N];
		}
		return tmp;
	}
	void print() const {
		printf("(");
		for(int i = 0; i < N; ++i) printf(" %g ", d[i]);
		printf(")\n");
	}
};


void read_samples(const char * filename, vector<int> & sample){
	FILE * fp = fopen(filename, "r");
	char buf[1024];
	int p;
	while(fgets(buf, 1024, fp) != nullptr){
		sscanf(buf, "%d", &p);
		sample.push_back(p);
	}
	fclose(fp);
}



typedef std::pair<int, double> dispair;
	struct CMP{
		int operator()(dispair a, dispair b){
			return a.second < b.second;
		}
	};

int main_backup(int argc, char * argv[]){

	const int K=36;
	printf("K=%d\n", K);
	const int K2 = K/2;
	typedef NVector<K2> CurvatureDescriptor;
	double radius;

	YXPathTracer pathtr;
	pathtr.init(argv[1]);	
	vector<int> samples;
	read_samples(argv[2], samples);
	sscanf(argv[3], "%lf", &radius);
	vector<YXPath> paths;
	vector<CurvatureDescriptor> descriptors;
	for(int i = 0; i < samples.size(); ++i){
		paths.clear();
		pathtr.computePatch(samples[i], K, radius, paths);
		CurvatureDescriptor c;
		for(int j = 0; j < K2; ++j){
			c[j] = pathtr.approxCurvature(samples[i], 
				pathtr.convertToYXPoint3D(*paths[j].rbegin()), 
				pathtr.convertToYXPoint3D(*paths[j+K2].rbegin()));
		}
		c = c.shiftToLowest();
		descriptors.push_back(c);
		char buf[1024];
		sprintf(buf, "curve/%d.obj", i);
		pathtr.generateObj(paths, buf);
		//c.print();
		//printf("max element at: %ld\n", std::max_element(c.d, c.d+K2)-c.d);
	}

	for(int i = 0; i < samples.size(); ++i){
		vector<dispair> dis;
		for(int j = 0; j < samples.size(); ++j){
			dispair a;
			a.first = j;
			//a.second = (descriptors[i] - descriptors[j]).length2();
			a.second = descriptors[i].shiftDistance(descriptors[j]) ;
			dis.push_back(a);
		}
		std::sort(dis.begin(), dis.end(), CMP());
		printf("%d : ", i);
		for(int j = 1; j <= 10; ++j){
			printf("%3d\t", dis[j].first);
		}
		printf("\n");
	}

	return 0;
}

