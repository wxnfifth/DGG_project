
#include "stdafx.h"
#include "wxnMath.h"
#include "ich\Point3D.h"
#include "gl\glut.h"

double cos_from_edges(double const a,			//compute the cosine of the angle given the lengths of the edges
							 double const b,
							 double const c)
{
	if( a < 1e-50 || b < 1e-50 || c < 1e-50 ){
		return 1.0;
	}

	double result = (b*b + c*c - a*a)/(2.0*b*c);
	result = max(result, -1.0);
	return min(result, 1.0);
}


//void renderTorus(const CPoint3D& p , const CPoint3D& color , const double radius , const int& sub_times)
//{
//    glDisable(GL_LIGHTING);
//    glPushMatrix();
//    glTranslatef( p.x,p.y,p.z);
//    color.SetColor();
//    glutSolidTorus(radius * 0.2 , radius * 1.1 , sub_times , sub_times); 
//    glPopMatrix();
//    glEnable(GL_LIGHTING);
//    //void glutSolidTorus(GLdouble innerRadius,
//    //                GLdouble outerRadius,
//    //                GLint nsides, GLint rings);
//}

void renderSphere(const CPoint3D& p , const CPoint3D& color , const double radius , const int& sub_times)
{
    glDisable(GL_LIGHTING);
    glPushMatrix();
    glTranslatef( p.x,p.y,p.z);
    color.SetColor();
    glutSolidSphere(radius , sub_times , sub_times); 
    glPopMatrix();
    glEnable(GL_LIGHTING);
}


Point2D operator*(double times, const Point2D& pt)
{
	return Point2D(times*pt.x,times*pt.y);
}
Point2D operator/(const Point2D& pt,const double& times)
{
	return Point2D(pt.x/times,pt.y/times);
}

CPoint3D GetOGLPos(int x, int y)
{
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLfloat winX, winY, winZ;
    GLdouble posX, posY, posZ;
 
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );
 
    winX = (float)x;
    winY = (float)viewport[3] - (float)y;
    glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
    //printf("winZ %lf\n" , winZ);
 
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
 
    return CPoint3D(posX, posY, posZ);
}

void renderCylinder_convenient(CPoint3D p1, CPoint3D p2 , float radius , int subdivisions)
{
    renderCylinder_convenient(p1.x , p1.y , p1.z , p2.x , p2.y , p2.z , radius , subdivisions);
}

void renderCylinder_convenient(float x1 , float y1 , float z1, float x2 , float y2, float z2, float radius , int subdivisions)
{
    //the same quadric can be re-used for drawing many cylinders
    GLUquadricObj *quadric=gluNewQuadric();
    gluQuadricNormals(quadric, GLU_SMOOTH);
    renderCylinder(x1,y1,z1,x2,y2,z2,radius,subdivisions,quadric);
    gluDeleteQuadric(quadric);
}


void renderCylinder(float x1, float y1, float z1, float x2,float y2, float z2, float radius,int subdivisions,GLUquadricObj *quadric)
{
    float vx = x2-x1;
    float vy = y2-y1;
    float vz = z2-z1;

    //handle the degenerate case of z1 == z2 with an approximation
    if(vz == 0)
        vz = .0001;

    float v = sqrt( vx*vx + vy*vy + vz*vz );
    float ax = 57.2957795*acos( vz/v );
    if ( vz < 0.0 )
        ax = -ax;
    float rx = -vy*vz;
    float ry = vx*vz;
    glPushMatrix();

    //draw the cylinder body
    glTranslatef( x1,y1,z1 );
    glRotatef(ax, rx, ry, 0.0);
    gluQuadricOrientation(quadric,GLU_OUTSIDE);
    gluCylinder(quadric, radius, radius, v, subdivisions, 1);

    ////draw the first cap
    gluQuadricOrientation(quadric,GLU_INSIDE);
    gluDisk( quadric, 0.0, radius, subdivisions, 1);
    glTranslatef( 0,0,v );

    ////draw the second cap
    gluQuadricOrientation(quadric,GLU_OUTSIDE);
    gluDisk( quadric, 0.0, radius, subdivisions, 1);

    glPopMatrix();
}

bool calculateBarycentric(const CPoint3D v[],const CPoint3D& p,CPoint3D& barycentric_cordinate)
{
	double tri[3][3];
	double pos[3];
	double bary_coordinate[3];
	for(int i = 0; i < 3; ++i){
		tri[i][0] = v[i].x;
		tri[i][1] = v[i].y;
		tri[i][2] = v[i].z;
	}
	pos[0] = p.x;
	pos[1] = p.y;
	pos[2] = p.z;
	Barycentric<double> bary(tri[0], tri[1], tri[2]);
	if (!bary.IsValid()) {
		printf("degenerated triangle!\n");
		barycentric_cordinate = CPoint3D(1.0/3.0,1.0/3.0,1.0/3.0);
		return false;
	}
	bary.Eval(pos, bary_coordinate);
	barycentric_cordinate = CPoint3D(bary_coordinate[0],bary_coordinate[1],bary_coordinate[2]);
	//bary_point.Print("bary");
	return true;
}



void solveQuadraticEquation(const double& a , const double& b ,const double& c,vector<double>& result)
{
	if( b * b - 4 * a * c < 0 ) return;
	double delta = sqrt( b * b - 4 * a * c );
	result.push_back( ( - b + delta ) / 2.0 / a );
	result.push_back( ( - b - delta ) / 2.0 / a );
}

double calcMidPoint(const double& dS1V1,const double& dS1V2,const double& dS2V1,const double &dS2V2,const double&dV1V2,double& distance_s1_middle){
	double dV1V2Sqr = dV1V2 * dV1V2;
	double dS1V1Sqr = dS1V1 * dS1V1;
	double dS1V2Sqr = dS1V2 * dS1V2;
	double dS2V2Sqr = dS2V2 * dS2V2;
	double dS2V1Sqr = dS2V1 * dS2V1;
	double cosTheta1;
	const double eps = 1e-16;
	if( fabs(dS1V1) < eps ){
		cosTheta1 = 0;
	}else{
		cosTheta1 = ( dV1V2Sqr + dS1V1Sqr - dS1V2Sqr ) / 2.0 / dV1V2 / dS1V1;
	}
	double cosTheta2;
	if( fabs(dS2V2) < eps ){
		cosTheta2 = 0;
	}else{
		cosTheta2 = ( dV1V2Sqr + dS2V2Sqr - dS2V1Sqr ) / 2.0 / dV1V2 / dS2V2;
	}

	//double d = ( dS1V1Sqr - dS2V2Sqr - dV1V2Sqr ) / 2.0 /
	//		   ( dS1V1 * cosTheta1 - dS2V2 * cosTheta2 - dV1V2 );
	double distance_v1_middle = ( dV1V2Sqr + dS2V2Sqr - dS1V1Sqr - 2 * dS2V2 * cosTheta2 * dV1V2 ) / 2.0 /
		( dV1V2 - dS1V1 * cosTheta1 - dS2V2 * cosTheta2 );

	distance_s1_middle = sqrt( dS1V1Sqr + distance_v1_middle * distance_v1_middle - 2 * dS1V1 * distance_v1_middle * cosTheta1 );
	double distance_s2_middle = sqrt( dS2V2Sqr + (dV1V2 - distance_v1_middle) * ( dV1V2 - distance_v1_middle ) - 2 * dS2V2 * ( dV1V2 - distance_v1_middle ) * cosTheta2 );
	if( fabs(distance_s1_middle - distance_s2_middle) > 1e-12 ){
		cout << "error in calcMidPoint" << distance_s1_middle << " " << distance_s2_middle << "\n";
	}

	return distance_v1_middle;
}
//http://en.wikipedia.org/wiki/Circumscribed_circle
//refer to the url above for the equations
CPoint3D calcCircumcenter(CPoint3D& p1 , CPoint3D& p2,CPoint3D& p3){
	double alpha, beta,gamma;
	alpha = ((p2-p3).LenSqr())*((p1-p2)^(p1-p3)) /
			2.0 / ((p1-p2)*(p2-p3)).LenSqr();
	beta  = ((p1-p3).LenSqr())*((p2-p1)^(p2-p3)) / 
			2.0 / ((p1-p2)*(p2-p3)).LenSqr();
	gamma = ((p1-p2).LenSqr())*((p3-p1)^(p3-p2)) /
			2.0 / ((p1-p2)*(p2-p3)).LenSqr();
	return alpha * p1 + beta * p2 + gamma * p3;
}

//https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas
//see the site for reference


CPoint3D rotatePoint(const CPoint3D& ptv1,const CPoint3D& normin,const CPoint3D& ptv2,const double d1, const double& cosTheta)
{
	//normin.Print("normin");
	CPoint3D norm = normin / normin.Len();
	//norm.Print("norm");
	double a = ptv1.x;
	double b = ptv1.y;
	double c = ptv1.z;
	double u = norm.x;
	double v = norm.y;
	double w = norm.z;
	double u2 = u * u;
	double v2 = v * v;
	double w2 = w * w;
	double x = ptv2.x;
	double y = ptv2.y;
	double z = ptv2.z;
	double cosT = cosTheta;
	double oneMinusCosT = 1 - cosT;
	double sinT = sqrt( 1.0 - cosT * cosT );
	double dv1v2 = (ptv1-ptv2).Len();
	CPoint3D p;
	//printf("a %lf b %lf c %lf u %lf v %lf w %f u2 %lf v2 %lf w2 %lf x %lf y %lf z %lf cosT %lf\n" , a,b,c,u,v,w,u2,v2,w2,x,y,z,cosT);
    p.x = (a*(v2 + w2) - u*(b*v + c*w - u*x - v*y - w*z)) * oneMinusCosT + x*cosT + (-c*v + b*w - w*y + v*z)*sinT;
	p.y = (b*(u2 + w2) - v*(a*u + c*w - u*x - v*y - w*z)) * oneMinusCosT + y*cosT + (c*u - a*w + w*x - u*z)*sinT;
    p.z = (c*(u2 + v2) - w*(a*u + b*v - u*x - v*y - w*z)) * oneMinusCosT + z*cosT + (-b*u + a*v - v*x + u*y)*sinT;
	//printf("d1 %lf %lf\n",d1,dv1v2);
	//p.Print("p");
	//ptv1.Print("ptv1");
	return p;

}

void calculate3DTriangleCoordinateIn2D(const CPoint3D p[] , Point2D vertex_2d[])
{
    double length[3];
    length[0] = (p[0] - p[1]).Len();
    length[1] = (p[1] - p[2]).Len();
    length[2] = (p[2] - p[0]).Len();
    calculate2DTriangleCoordinate(length , vertex_2d);
}


void calculate2DTriangleCoordinate(const double length[] , Point2D vertex_2d[])
{
	double costheta = (length[0] *length[0] + length[2] * length[2] - length[1] * length[1] ) /2.0 / length[0] / length[2];
	double sintheta = sqrt( 1 - costheta * costheta );
	vertex_2d[0] = Point2D(0,0);
	vertex_2d[1] = Point2D(length[0],0);
	vertex_2d[2] = Point2D(length[2] * costheta , length[2] * sintheta);
}

double calculateCosFromTriangleLength(const double& a,const double& b , const double& c,int* cnt)
{
	double temp_eps = 1e-14;
	if( fabs(a) < temp_eps || fabs(b) < temp_eps || fabs(c) < temp_eps ){
		return 0.0;
	}
	double result = ( a * a + b * b - c * c ) / 2.0 / a / b;
	if( result > 1.0 ){
		//		std::cout.precision(10);
		//std::cout << std::scientific;
		//	cout << "result " << result << "\n";
		if( cnt != NULL ) (*cnt)++;
		result = 1.0;
	}else if( result < -1.0 ){
		if( cnt != NULL ) (*cnt)++;
		result = -1.0;
	}
	return result;
}

double calculateOppositeVertexInquadrangle(const double length[])
{
	double costheta1 = calculateCosFromTriangleLength(length[0],length[1],length[2]);
	double sintheta1 = sqrt(1.0 - costheta1 * costheta1 );
	double costheta2 = calculateCosFromTriangleLength(length[1],length[3],length[4]);
	double sintheta2 = sqrt(1.0 - costheta2 * costheta2 );
	double costheta = costheta1 * costheta2 - sintheta1 * sintheta2;
	return length[0] * length[0] + length[3] * length[3] - 2 * length[0] * length[3] * costheta;
}


void map3DTriangleTo2D(const CPoint3D pt[],CPoint3D pt2D[])
{
	double len[3];
	double costheta;
	double sintheta;
	for(int i = 0; i < 3; ++i){
		len[i] = (pt[i] - pt[(i+1)%3]).Len();
		
	}
	costheta = (len[0] *len[0] + len[2] * len[2] - len[1] * len[1] ) /2.0 / len[0] / len[2];
	sintheta = sqrt( 1 - costheta * costheta );
	pt2D[0] = CPoint3D(0,0,0);
	pt2D[1] = CPoint3D(len[0],0,0);
	pt2D[2] = CPoint3D(len[2] * costheta , len[2] * sintheta , 0 );
}



double ptInEdgeToProp(const CPoint3D& ptIn,const CPoint3D& EP1,const CPoint3D& EP2)
{
	double prop1;
	double prop2;
	double diffInX,diffInY,diffInZ;
	diffInX = fabs(EP1.x-EP2.x);
	diffInY = fabs(EP1.y-EP2.y);
	diffInZ = fabs(EP1.z-EP2.z);
	if( diffInX > diffInY && diffInX > diffInZ ){
		prop1 = (ptIn.x - EP1.x) / (EP2.x - EP1.x);
	}else if( diffInY > diffInX && diffInY > diffInZ ){
		prop1 = (ptIn.y - EP1.y) / (EP2.y - EP1.y);
	}else{
		prop1 = (ptIn.z - EP1.z) / (EP2.z - EP1.z);
	}
	prop2 = 1 - prop1;
	if( !(EP1 * prop2 + EP2 * prop1 == ptIn)){
		double len = (EP2 - EP1).Len(); 
		ptIn.Print("ptIn2DEdgeToProp");
		printf( "prop1 %lf len %lf len1 + len2 %lf\n" , prop1 , len , (ptIn-EP1).Len() + (EP2 - ptIn).Len() );
		printf( "diff %g\n" , (EP1 * prop2 + EP2 * prop1 - ptIn).Len());
		printf("\n");
	}
	return prop1;
}

CPoint3D ptIn3DEdgeTo2D(const CPoint3D& ptIn,const CPoint3D& EP3D1,const CPoint3D& EP3D2,const CPoint3D& EP2D1,const CPoint3D& EP2D2)
{
	//no use for now
//http://hewiki.heroengine.com/wiki/Technical:Cartesian_To_Barycentric_Coordinates
	//CPoint3D v0 = pt3D[0] - pt3D[2];
	//CPoint3D v1 = pt3D[1] - pt3D[2];
	//CPoint3D v2 = ptIn - pt3D[2];
	//double d00 = v0 ^ v0;
	//double d01 = v0 ^ v1;
	//double d11 = v1 ^ v1;
	//double d20 = v2 ^ v0;
	//double d21 = v2 ^ v1;
	//double denom = d00 * d11 - d01 * d01;
	//double a , b , c;
	//if( denom != 0 ){
	//	b = ( d11 * d20 - d01 * d21 ) / denom;
	//	c = ( d00 * d21 - d01 * d20 ) / denom;
	//}
	//if( (ptIn - (b * pt3D[0] + (1-b-c)*pt3D[1] +c*pt3D[2]) ).Len() > 1e-10){
	//	ptIn.Print("pt");
	//	(b * pt3D[0] + (1-b-c)*pt3D[1] +c*pt3D[2]).Print("calc pt");
	//	printf( "%lf %lf %lf\n" , b , 1-b-c,c);
	//	printf( "error!\n");

	//}

	//return CPoint3D::MAX_VALUE;
	////if( b >= 0 && c >= 0 && ( b + c ) <= 1 ){
	////	return CPoint3D(b,1-b-c,c);
	////}

	double prop1 , prop2;
	prop1 = ptInEdgeToProp(ptIn,EP3D1,EP3D2);
	prop2 = 1 - prop1;
	//prop2 = (EP3D2 - ptIn).Len() / len;
	//if( (prop1 + prop2) != 1 ){
	//	printf( "prop1 + prop2 = %lf \n" , prop1 + prop2 -1 );
	//}
	return EP2D1 * prop2 + EP2D2 * prop1;

}

double calcAngleFromTriangle(const double& a , const double& b , const double& c)
{
	return acos((a*a+b*b-c*c)/2.0/a/b);
}
//template<typename  T> bool vectorEqual(const vector<T>& va,const vector<T>& vb);
//template<typename  T> bool vectorEqual(const T va[],const T vb[],int n);


CPoint3D findU1(const CPoint3D vx2D,const CPoint3D& p1,double lP1Midr1,double theta)
{
	CPoint3D u1 = (vx2D - p1).Normalize() * lP1Midr1;
	double costheta = cos(theta);
	double sintheta = sin(theta);
	CPoint3D u2 = CPoint3D(u1.x * costheta - u1.y * sintheta , u1.x * sintheta + u1.y * costheta ,0);
	return u2;

}

CPoint3D ptIn2DTriangleTo3D(const CPoint3D v[],const CPoint3D v2D[],const CPoint3D& p2D)
{
	CPoint3D v01 = v2D[1] - v2D[0];
	CPoint3D v02 = v2D[2] - v2D[0];
	CPoint3D v0p = p2D - v2D[0];
	double lambda1,lambda2;
	lambda1 = (v0p.x * v02.y - v0p.y * v02.x ) / ( v01.x * v02.y - v02.x * v01.y );
	lambda2 = (v0p.x * v01.y - v0p.y * v01.x ) / ( v02.x * v01.y - v01.x * v02.y );
	if( p2D != ( 1 - lambda1 - lambda2 ) * v2D[0] + lambda1 * v2D[1] + lambda2 * v2D[2] ){
		printf("convert 2d pt to 3d error!\n");
	}
	if( !(lambda1 == lambda1) ){
		p2D.Print("p2D");
		printf( "tmp %lf \n"  ,v01.x * v02.y - v02.x * v01.y);
		printf( "tmp up %lf\n" , (v0p.x * v02.y - v0p.y * v02.x));
		printf("lambda %lf %lf \n" , lambda1 , lambda2 );
	}
	return 	( 1 - lambda1 - lambda2 ) * v[0] + lambda1 * v[1] + lambda2 * v[2];
}
CPoint3D ptIn2DTriangleToBarycentric(const CPoint3D v2D[],const CPoint3D& p2D)
{
	CPoint3D v01 = v2D[1] - v2D[0];
	CPoint3D v02 = v2D[2] - v2D[0];
	CPoint3D v0p = p2D - v2D[0];
	double lambda1,lambda2;
	lambda1 = (v0p.x * v02.y - v0p.y * v02.x ) / ( v01.x * v02.y - v02.x * v01.y );
	lambda2 = (v0p.x * v01.y - v0p.y * v01.x ) / ( v02.x * v01.y - v01.x * v02.y );
	if( p2D != ( 1 - lambda1 - lambda2 ) * v2D[0] + lambda1 * v2D[1] + lambda2 * v2D[2] ){
		printf("ptIn2DTriangleToBarycentric error!\n");
	}
	if( !(lambda1 == lambda1) ){
		p2D.Print("p2D");
		printf( "tmp %lf \n"  ,v01.x * v02.y - v02.x * v01.y);
		printf( "tmp up %lf\n" , (v0p.x * v02.y - v0p.y * v02.x));
		printf("lambda %lf %lf \n" , lambda1 , lambda2 );
	}
	return 	CPoint3D( 1 - lambda1 - lambda2 , lambda1 , lambda2);
}
//CPoint3D pointIn2DTriangleToBarycentric(const Point2D v2D[],const Point2D& p2D)
//{
//	Point2D A = v2D[0];
//	Point2D B = v2D[1];
//	Point2D C = v2D[2];
//	Point2D v0 = C - A;
//	Point2D v1 = B - A;
//	Point2D v2 = p2D - A;
//	double dot00 = v0.dotProduct(v0);
//	double dot01 = v0.dotProduct(v1);
//	double dot02 = v0.dotProduct(v2);
//	double dot11 = v1.dotProduct(v1);
//	double dot12 = v1.dotProduct(v2);
//	double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
//	if( invDenom < LENGTH_EPSILON_CONTROL ){
//		Point2D A = v2D[2];
//		Point2D B = v2D[0];
//		Point2D C = v2D[1];
//		Point2D v0 = C - A;
//		Point2D v1 = B - A;
//		Point2D v2 = p2D - A;
//		double dot00 = v0.dotProduct(v0);
//		double dot01 = v0.dotProduct(v1);
//		double dot02 = v0.dotProduct(v2);
//		double dot11 = v1.dotProduct(v1);
//		double dot12 = v1.dotProduct(v2);
//		double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
//		double lambda_u = (dot11 * dot02 - dot01 * dot12) * invDenom;
//		double lambda_v = (dot00 * dot12 - dot01 * dot02) * invDenom;
//		if( A *(1-lambda_u-lambda_v) + B*lambda_v+C*lambda_u != p2D){
//			printf("changed\n");
//			v0.print("v0");
//			v1.print("v1");
//			v2.print("v2");
//			v2D[0].print("v2D[0]");
//			v2D[1].print("v2D[1]");
//			v2D[2].print("v2D[2]");
//			p2D.print("p2D");
//			printf("lambda_u %g lambda_v %g invDenom %g\n" , lambda_u,lambda_v,invDenom);
//			(A *(1-lambda_u-lambda_v) + B*lambda_v+C*lambda_u).print("barycentric");
//			p2D.print("p2D");
//			printf("error in pointIn2DTriangleToBarycentric\n");
//		}
//		return CPoint3D(lambda_v,lambda_u,1-lambda_u-lambda_v);
//	}
//	double lambda_u = (dot11 * dot02 - dot01 * dot12) * invDenom;
//	double lambda_v = (dot00 * dot12 - dot01 * dot02) * invDenom;
//	if( A *(1-lambda_u-lambda_v) + B*lambda_v+C*lambda_u != p2D){
//		v2D[0].print("v2D[0]");
//		v2D[1].print("v2D[1]");
//		v2D[2].print("v2D[2]");
//		p2D.print("p2D");
//		printf("lambda_u %g lambda_v %g invDenom %g\n" , lambda_u,lambda_v,invDenom);
//		(A *(1-lambda_u-lambda_v) + B*lambda_v+C*lambda_u).print("barycentric");
//		p2D.print("p2D");
//		printf("error in pointIn2DTriangleToBarycentric\n");
//	}
//	return CPoint3D(1-lambda_u-lambda_v,lambda_v,lambda_u);
//}


//template <class Real>
//void rms::BarycentricCoords( const Vector3<Real> & vTriVtx1, 
//							 const Vector3<Real> & vTriVtx2,
//							 const Vector3<Real> & vTriVtx3,
//							 const Vector3<Real> & vVertex,
//							 Real & fBary1, Real & fBary2, Real & fBary3 )
//{
//
//	Wml::Vector3<Real> kV02 = vTriVtx1 - vTriVtx3;
//    Wml::Vector3<Real> kV12 = vTriVtx2 - vTriVtx3;
//    Wml::Vector3<Real> kPV2 = vVertex - vTriVtx3;
//
//    Real fM00 = kV02.Dot(kV02);
//    Real fM01 = kV02.Dot(kV12);
//    Real fM11 = kV12.Dot(kV12);
//    Real fR0 = kV02.Dot(kPV2);
//    Real fR1 = kV12.Dot(kPV2);
//    Real fDet = fM00*fM11 - fM01*fM01;
////    assert( Wml::Math<Real>::FAbs(fDet) > (Real)0.0 );
//    Real fInvDet = ((Real)1.0)/fDet;
//
//    fBary1 = (fM11*fR0 - fM01*fR1)*fInvDet;
//    fBary2 = (fM00*fR1 - fM01*fR0)*fInvDet;
//    fBary3 = (Real)1.0 - fBary1 - fBary2;
//}

//http://www.gamedev.net/topic/451357-computing-barycentric-coordinates/
CPoint3D pointIn2DTriangleToBarycentricRMS(const Point2D v2D[],const Point2D& p2D)
{
	const Point2D& v1 = v2D[0];
	const Point2D& v2 = v2D[1];
	const Point2D& v3 = v2D[2];
    Point2D kV02 = v1 - v3;
    Point2D kV12 = v2 - v3;
    Point2D kPV2 = p2D - v3;
    double fM00 = kV02.dotProduct(kV02);
    double fM01 = kV02.dotProduct(kV12);
    double fM11 = kV12.dotProduct(kV12);
    double fR0 = kV02.dotProduct(kPV2);
    double fR1 = kV12.dotProduct(kPV2);
    double fDet = fM00*fM11 - fM01*fM01;
    double fInvDet = (1.0)/fDet;
    CPoint3D rt;
    rt.x = (fM11*fR0 - fM01*fR1)*fInvDet;
    rt.y = (fM00*fR1 - fM01*fR0)*fInvDet;
    rt.z = 1.0 - rt.x - rt.y;

	return rt;
}
CPoint3D pointIn2DTriangleToBarycentric(const Point2D v2D[],const Point2D& p2D)
{
	const Point2D& v1 = v2D[0];
	const Point2D& v2 = v2D[1];
	const Point2D& v3 = v2D[2];
	double F  = (v3-p2D).crossProduct(v1-p2D) / (v2-v1).crossProduct(v3-v1);
	double G  = (v1-p2D).crossProduct(v2-p2D) / (v2-v1).crossProduct(v3-v1);
	if( (v2-v1).crossProduct(v3-v1) > DOUBLE_EPSILON ){
		if( (1-F-G)*v1+F*v2+G*v3 != p2D ){
			v1.print("v1");
			v2.print("v2");
			v3.print("v3");
			//(v2-v1).print("v2-v1");
			//(v3-v1).print("v3-v1"); 
			printf("error %g %g %g!\n",F,G,1-F-G);
            return CPoint3D(1/3.0,1/3.0,1/3.0);
		}
	}
	return CPoint3D(1-F-G,F,G);
}


//http://www.blackpawn.com/texts/pointinpoly/default.html
CPoint3D pointIn3DTriangleToBarycentric(const CPoint3D v3D[],const CPoint3D& p3D)
{
	CPoint3D A = v3D[0];
	CPoint3D B = v3D[1];
	CPoint3D C = v3D[2];
	CPoint3D v0 = C - A;
	CPoint3D v1 = B - A;
	CPoint3D v2 = p3D - A;
	double dot00 = v0 ^ v0;
	double dot01 = v0 ^ v1;
	double dot02 = v0 ^ v2;
	double dot11 = v1 ^ v1;
	double dot12 = v1 ^ v2;
	double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	double lambda_u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	double lambda_v = (dot00 * dot12 - dot01 * dot02) * invDenom;
	
	if( invDenom > 1e8 ){
		
	}
	if( ( (A *(1-lambda_u-lambda_v) + B*lambda_v+C*lambda_u - p3D).Len() )  > 0.01 ){
		printf("divide result %lf %lf %lf\n" , v0.x / v1.x , v0.y / v1.y , v0.z / v1.z );
		A.Print("A");
		B.Print("B");
		C.Print("C");
		v0.Print("v0");
		v1.Print("v1");
		v2.Print("v2");
		p3D.Print("p3D");
		printf("coef %e %e %e\n" , 1-lambda_u-lambda_v,lambda_v,lambda_u);
		printf("invDenom %e\n" ,invDenom);
		printf("dot11 * dot02 - dot01 * dot12 %lf\n" , dot11 * dot02 - dot01 * dot12);
		printf("dot00 * dot12 - dot01 * dot02 %lf\n" , dot00 * dot12 - dot01 * dot02); 
		printf("diff %e  " , ( (A *(1-lambda_u-lambda_v) + B*lambda_v+C*lambda_u - p3D).Len() ));
		
		//(A *(1-lambda_u-lambda_v) + B*lambda_v+C*lambda_u).Print("strange");
		//p3D.Print("p3D");
		//printf("invDenom %lf\n" , invDenom);
		printf("error in pointIn3DTriangleToBarycentric\n");
	}

	return CPoint3D(1-lambda_u-lambda_v,lambda_v,lambda_u);
}
//CPoint3D pointIn3DTriangleToBarycentric(const CPoint3D v[],const CPoint3D& p)
//{
//	
//	double triArea = ((v[1] - v[0]) * (v[2] - v[0])).Len();
//	double ud = ((v[1]-p) * (v[2] - p)).Len() / triArea;
//	double vd = ((v[0]-p) * (v[2] - p)).Len() / triArea;
//	double wd = ((v[0]-p) * (v[1] - p)).Len() / triArea;
//
//	if( (ud * v[0] + vd * v[1] + wd * v[2] - p ).Len() >1e-1 ){
//		printf("triArea %lf\n" , triArea);
//	}
//	return CPoint3D(ud,vd,wd);
//}

//Vector3d Tri::barycentric(Vector3d p) {
//    double triArea = (p1 - p0).cross(p2 - p0).norm() * 0.5;
//    double u = ((p1 - p).cross(p2 - p).norm() * 0.5) / triArea;
//    double v = ((p0 - p).cross(p2 - p).norm() * 0.5) / triArea;
//    double w = ((p0 - p).cross(p1 - p).norm() * 0.5) / triArea;
//    return Vector3d(u,v,w);
//}
Point2D ptIn3DTriangleTo2D(const CPoint3D v3D[],const Point2D v2D[],const CPoint3D& p3D)
{
	CPoint3D barycentric_coordnate;// = pointIn3DTriangleToBarycentric(v3D,p3D);

	if( calculateBarycentric(v3D,p3D,barycentric_coordnate) ){

	}else{
		printf("ptIn3DTriangleTo2D error !\n");
	}
	return barycentric_coordnate.x * v2D[0] + barycentric_coordnate.y * v2D[1] +
		barycentric_coordnate.z * v2D[2];
}




CPoint3D findLineIntersectionPt(const CPoint3D& p12D,const CPoint3D& u1,const CPoint3D& p22D,const CPoint3D& u2)
{
	double lambda1 = ( (p12D.x - p22D.x) * u2.y - (p12D.y - p22D.y) * u2.x ) / ( u1.y * u2.x - u1.x * u2.y );
	return p12D + lambda1 * u1;
}

void drawLine(const CPoint3D& p1,const CPoint3D& p2,vector<CPoint3D>&  vList,vector<vector<int>>& triList)
{
	CPoint3D tmp(1,0,0);
	CPoint3D direct = ((p1-p2)*tmp).Normalize() * 0.5;
	vList.push_back(p1);
	vList.push_back(p2);
	CPoint3D p0 = p1 + direct;
	int num = 3;
	CPoint3D p = p0;
	for(int i = 0; i < num; ++i){
		vList.push_back(p);
		CPoint3D tmpp;
		tmpp = rotatePoint(p1,p2-p1,p,1.0,cos( PI * 2.0 / num) );
		p = tmpp;
	}
	for(int i = 0; i < num; ++i){
		vList.push_back(vList[2+i] + p2 - p1 );
	}
	for(int i = 0; i < num; ++i){
		vector<int> tri;
		tri.push_back(0);
		tri.push_back(i+2);
		tri.push_back((i+1)%num+2);
		triList.push_back(tri);
	}
	for(int i = 0; i < num; ++i){
		vector<int> tri;
		tri.push_back(1);
		tri.push_back(i+2+num);
		tri.push_back((i+1)%10+2+num);
		triList.push_back(tri);
	}
	for(int i = 0; i < num; ++i){
		vector<int> tri(3);
		tri[0] = i + 2;
		tri[1] = ( i + 1 ) % num + 2;
		tri[2] = i + 2 + num;
		triList.push_back(tri);
		tri[0] = i + 2 + num;
		tri[1] = ( i + 1 ) % num + 2 + num;
		tri[2] = ( i + 1 ) % num + 2;
		triList.push_back(tri);
	}
}


double calculateTriangleQuality(const double& a ,const double& b,const double& c){
	double  s = (a+b+c)/2.0;
	double in_circle_radius = sqrt( (s-a)*(s-b)*(s-c)/s );
	return 6.0 / sqrt(3.0) * in_circle_radius / max(max(a,b),c);
}

double calculateOppositeVertexDistanceInquadrangle(const double length[])
{
	double costheta1 = calculateCosFromTriangleLength(length[0],length[1],length[2]);
	double sintheta1 = sqrt(1.0 - costheta1 * costheta1 );
	double costheta2 = calculateCosFromTriangleLength(length[1],length[3],length[4]);
	double sintheta2 = sqrt(1.0 - costheta2 * costheta2 );
	double costheta = costheta1 * costheta2 - sintheta1 * sintheta2;
	return max(0.,length[0] * length[0] + length[3] * length[3] - 2 * length[0] * length[3] * costheta);
}

void map3DTriangleTo2D(const CPoint3D pt[],Point2D pt2D[])
{
	double len[3];
	double costheta;
	double sintheta;
	for(int i = 0; i < 3; ++i){
		len[i] = (pt[i] - pt[(i+1)%3]).Len();
		
	}
	costheta = (len[0] *len[0] + len[2] * len[2] - len[1] * len[1] ) /2.0 / len[0] / len[2];
	sintheta = sqrt( 1 - costheta * costheta );
	pt2D[0] = Point2D(0,0);
	pt2D[1] = Point2D(len[0],0);
	pt2D[2] = Point2D(len[2] * costheta , len[2] * sintheta );
}
bool calculateBarycentric(const CPoint3D v[],const CPoint3D& p,double barycentric_cordinate[])
{
	double tri[3][3];
	double pos[3];
	for(int i = 0; i < 3; ++i){
		tri[i][0] = v[i].x;
		tri[i][1] = v[i].y;
		tri[i][2] = v[i].z;
	}
	pos[0] = p.x;
	pos[1] = p.y;
	pos[2] = p.z;
	Barycentric<double> bary(tri[0], tri[1], tri[2]);
	if (!bary.IsValid()) {
		printf("degenerated triangle!\n");
		//barycentric_cordinate = CPoint3D(1.0/3.0,1.0/3.0,1.0/3.0);
		barycentric_cordinate[0] = 1.0 / 3.0;
		barycentric_cordinate[1] = 1.0 / 3.0;
		barycentric_cordinate[2] = 1.0 / 3.0;
		return false;
	}
	bary.Eval(pos, barycentric_cordinate);
	//barycentric_cordinate = CPoint3D(bary_coordinate[0],bary_coordinate[1],bary_coordinate[2]);
	//bary_point.Print("bary");
	return true;
}
CPoint3D pointIn2DTriangleTo3D(const CPoint3D v[],const Point2D v2D[],const Point2D& p2D)
{
	CPoint3D barycentric = pointIn2DTriangleToBarycentric(v2D,p2D);
	return barycentric.x * v[0] + barycentric.y * v[1] + barycentric.z * v[2];
}


//CPoint3D tmp2dPt = ptIn3DTriangleTo2D(vertexList[vertexOnEdge[edgeId[j]]],pt3D,pt2D);