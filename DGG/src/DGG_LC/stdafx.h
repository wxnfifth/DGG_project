// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once
#include <Windows.h>
//#define _AFXDLL
//#include <afx.h>
#include "targetver.h"
#include <stdio.h>
#include <tchar.h>


#define _USE_MATH_DEFINES
#include <math.h>
#include <ctime>
#include <iomanip>
#include <vector>
#include <map>
#include <stack>
#include <set>
#include <list>
#include <vector>
#include <queue>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include <limits>
#include <string>
#include <cassert>
#include <gl\gl.h>
#include <gl\glu.h>
#include <GL\glut.h>
#include <omp.h>
#define MAKE_STRING(a) (# a)
#define PRINT_VAR(a) (std::cerr << MAKE_STRING(a) << " " << (a) << "\n")

//#include <ompassem.h>
#include <iomanip>
//#pragma comment(lib, "OPENGL32.LIB")
//#pragma comment(lib, "GLU32.LIB")
#pragma warning(disable:4018)

//#pragma comment(lib, "glut32.lib")
//using namespace std;

const double DOUBLE_EPSILON = 1e-10;
const double LENGTH_EPSILON_CONTROL = 1e-6;
const double WXN_LENGTH_EPSILON = 1e-8;
const double PI = 3.1415926535897932384626;
const double RateOfNormalShift = 1.5e-3;
const double ToleranceOfConvexAngle = 5e-3;

// TODO: reference additional headers your program requires here
