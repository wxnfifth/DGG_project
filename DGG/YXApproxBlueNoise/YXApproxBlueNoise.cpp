// YXApproxBlueNoise.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "YXApproxBlueNoise.h"

void testApproxBlueNoise(){
	YXApproxBlueNoise bn;
	YXTimer t;
	t.start();
	bn.load("C:\\Users\\YingXiang\\GamelabSVN\\YingXiang\\3DModels\\bunny\\bunny_nf144k.m");
	printf("Time used for loading: %f seconds\n", t.gettime());
	int N, P;
	//t.start();
	//N = 10000;
	//P = 10;
	//bn.sampling(N, P);
	//printf("Time used for sampling: %lf seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
	//	t.gettime(), N, P, bn.result.size());

	t.start();
	N = 50000;
	P = 10;
	bn.sampling(N, P);
	printf("Time used for sampling: %f seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
		t.gettime(), N, P, bn.result.size());

	t.start();
	N = 200000;
	P = 10;
	bn.sampling(N, P);
	printf("Time used for sampling: %f seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
		t.gettime(), N, P, bn.result.size());

	t.start();
	N = 1000000;
	P = 5;
	bn.sampling(N, P);
	printf("Time used for sampling: %f seconds. \n Want %d samples with %d xRate, generated %d sampels\n\n",
		t.gettime(), N, P, bn.result.size());
}

int _tmain(int argc, _TCHAR* argv[])
{
	testApproxBlueNoise();
	return 0;
}

