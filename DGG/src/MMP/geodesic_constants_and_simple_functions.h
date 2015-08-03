#ifndef GEODESIC_CONSTANTS
#define GEODESIC_CONSTANTS

// some constants and simple math functions

#include <assert.h>
#include <math.h>
#include <limits>
#include <fstream>
#include <cstdio>
#include <string>
#include <map>

namespace geodesic{

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//double const GEODESIC_INF = std::numeric_limits<double>::max();
double const GEODESIC_INF = 1e100;

//in order to avoid numerical problems with "infinitely small" intervals,
//we drop all the intervals smaller than SMALLEST_INTERVAL_RATIO*edge_length
double const SMALLEST_INTERVAL_RATIO = 1e-6;		
//double const SMALL_EPSILON = 1e-10;


inline double cos_from_edges(double const a,			//compute the cosine of the angle given the lengths of the edges
							 double const b,
							 double const c)
{
	assert(a>1e-50);
	assert(b>1e-50);
	assert(c>1e-50);

	double result = (b*b + c*c - a*a)/(2.0*b*c);
	result = std::max(result, -1.0);
	return std::min(result, 1.0);
}

inline double angle_from_edges(double const a,			//compute the cosine of the angle given the lengths of the edges
							   double const b,
							   double const c)
{
	return acos(cos_from_edges(a,b,c));
}

template<class Points, class Faces>
inline bool read_mesh_from_file(const char* filename,
								Points& points,
								Faces& faces, 
								std::vector<int> &realIndex, 
								int &originalVertNum)
{
	std::ifstream file(filename);
	assert(file.is_open());
	if(!file.is_open()) return false;
	/*
	unsigned num_points;
	file >> num_points;
	assert(num_points>=3);

	unsigned num_faces;
	file >> num_faces;

	points.resize(num_points*3);
	for(typename Points::iterator i=points.begin(); i!=points.end(); ++i)
	{
		file >> *i;
	}

	faces.resize(num_faces*3);
	for(typename Faces::iterator i=faces.begin(); i!=faces.end(); ++i)
	{
		file >> *i;
	}
	file.close();
	*/
	char type;
	std::string curLine;
	double coord[3];
	unsigned int vtxIdx[3];
	std::map<string, int> mapForDuplicate;
	originalVertNum = 0;

	while(getline(file, curLine))
	{
		if (curLine.size() < 2) continue;
		if (curLine[0] == 'v' && curLine[1] != 't')
		{
			std::map<string, int>::iterator pos = mapForDuplicate.find(curLine);
			if (pos == mapForDuplicate.end())
			{
				int oldSize = mapForDuplicate.size();
				realIndex.push_back(oldSize);
				mapForDuplicate[curLine] = oldSize;
				sscanf(curLine.c_str(), "v %lf %lf %lf", &coord[0], &coord[1], &coord[2]);
				for (int i = 0;i < 3;++i) points.push_back(coord[i]);
			}
			else
			{
				realIndex.push_back(pos->second);
			}
			++originalVertNum;
		}
		else if (curLine[0] == 'f')
		{
			unsigned tex;
			if (curLine.find('/') != std::string::npos)
				sscanf(curLine.c_str(), "f %d/%d %d/%d %d/%d", &vtxIdx[0], &tex, &vtxIdx[1], &tex, &vtxIdx[2], &tex);
			else
				sscanf(curLine.c_str(), "f %d %d %d", &vtxIdx[0], &vtxIdx[1], &vtxIdx[2]);
			
			vtxIdx[0] = realIndex[vtxIdx[0]-1];
			vtxIdx[1] = realIndex[vtxIdx[1]-1];
			vtxIdx[2] = realIndex[vtxIdx[2]-1];
			if (vtxIdx[0] == vtxIdx[1] || vtxIdx[0] == vtxIdx[2] || vtxIdx[1] == vtxIdx[2]) continue;

			for (int i = 0;i < 3;++i) faces.push_back(vtxIdx[i]);
		}
	}
	file.close();
#ifndef DOTEST
	printf("There are %d non-coincidence vertices.\n", points.size() / 3);
#endif
	return true;
}

} //geodesic

#endif	//GEODESIC_CONSTANTS
