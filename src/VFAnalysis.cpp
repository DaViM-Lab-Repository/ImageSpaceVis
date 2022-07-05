

//#include "VFAnalysis.h"

#include "Skeleton.h"

#include <math.h>

double r1, r2, i1, i2;  //Make these Jaccobian Matrix components global !!!!!
icMatrix3x3 FieldMatrix;
double x_cp, y_cp;      //use to get the singularity coordinates 3/25/06

std::vector<Singularity> singularity_list;
std::vector<int> singular_tris;

std::vector<Trajectory> trajs;
std::vector<int> incoming_seps;
std::vector<int> outgoing_seps;

extern Polyhedron *poly;

unsigned char Integrator_opt = 0;  // the user specified integrator type (0 - Euler, 1 - RK2, 2 - RK4 ...)


void capture_singular_tris(void)
{
	unsigned int i, j;
	Triangle *f;

	icVector2 vec[3];      //vectors at three vertices
    double  theta[3];      //for storing the angles between two vector for Gauss circle judgement

	double  vec_ang[3];  //the angle for each vector under the polar frame of current triangle
	double  ang_sum;

	////Initialize
	singularity_list.clear();
	singular_tris.clear();


	////Calculate the Gaussian angle
	for (i = 0; i <poly->ntris; i++) {
		f = poly->tlist[i];
		ang_sum = 0;

		//For each triangle, calculate the vector for all vertices
		for (j=0; j < f->nverts; j++) {
			////using the vectors stored in the vertices
			vec[j] =  f->verts[j]->vec;

			vec_ang[j] = atan2(vec[j].entry[1], vec[j].entry[0]);

			if(vec_ang[j] < 0) vec_ang[j] += 2 * PI;
		}

		for(j = 0; j < f->nverts; j++)
		{
			theta[j] = vec_ang[(j+1)%3] - vec_ang[j];

			if( theta[j] < -PI)
				theta[j] += 2 * PI;
			
			if( theta[j] > PI)
				theta[j] -= 2 * PI;


			ang_sum += theta[j];
		}

		if(fabs(ang_sum) >= (2 * PI - 0.5))
		{
			//The triangle must have singularities inside, mark it as yellow color
			//Still need to judge whether it is one of current singularities or not
			singular_tris.push_back(i);
		}
	}

	////Calculate the coordinates of being found singularities
	if(!singular_tris.empty())
	{
		compute_fixedpt_info();
	}
}


void get_loc_Jacobian(int tri, double Jacobian[2][2], double &c, double &f)
{
	double x[3], y[3], vx[3], vy[3];
	int i, j;
	Triangle *face;

	///These information has been obtained before
	///Can we reuse them and enhance the performance 
    face = poly->tlist[tri];

	//For each triangle, calculate the vector for all vertices
	for (j=0; j<face->nverts; j++) {
		x[j] = face->verts[j]->x;
		y[j] = face->verts[j]->y;

		vx[j] = face->verts[j]->vx;
		vy[j] = face->verts[j]->vy;
	}

	/////Second way to calculate the a b c d e f values

	double coord[3][3],  *inver_coord ;  //inver_coord[3][3];

	for(i = 0; i < 3; i++)
	{
		coord[0][i] = x[i];
		coord[1][i] = y[i];
		coord[2][i] = 1.;
	}

	inver_coord = MatrixOpp((double*)coord, 3, 3);

	icMatrix3x3 result, rightM;
	result.set(vx[0],vx[1],vx[2],  vy[0],vy[1],vy[2],  1,1,1);
	rightM.set(inver_coord[0], inver_coord[1], inver_coord[2],
				inver_coord[3], inver_coord[4], inver_coord[5],
				inver_coord[6], inver_coord[7], inver_coord[8]);


	result.rightMultiply(rightM);

	FieldMatrix.set(result);

	double a, b, d, e;

	a = result.entry[0][0];
	b = result.entry[0][1];
	c = result.entry[0][2];
	d = result.entry[1][0];
	e = result.entry[1][1];
	f = result.entry[1][2];

	//need to store it as a part of the information of the elements or unknown singularities
	Jacobian[0][0] = a;
	Jacobian[0][1] = b;
	Jacobian[1][0] = d;
	Jacobian[1][1] = e;

	//use to calculate the coordinates of the singularity 3/25/06
	x_cp = (f*b - c*e)/(a*e - d*b);
	y_cp = (c*d - a*f)/(a*e - b*d);
}


void compute_fixedpt_info(void)
{
	int i, j;
	Triangle *tri;
	double x_loc, y_loc; 
	icVector2 glob;
	double Jacobian[2][2];
	double alpha[3];

	////Initialize the unknown singularities link

	for(i = 0; i < singular_tris.size(); i++)
	{
		//For each being captured triangle, compute the coordinates of singularities inside it
		tri = poly->tlist[singular_tris[i]];

		double c, f;

		get_loc_Jacobian(singular_tris[i], Jacobian, c, f);

		Singularity onesing;

		onesing.Jacobian.set(Jacobian);
		onesing.gpos.set (x_cp, y_cp);
		onesing.TriangleID = singular_tris[i];
		onesing.type = get_singularity_type(Jacobian);

		if (onesing.type == SADDLE)
		{
			onesing.outgoing.entry[0] = -onesing.Jacobian.entry[0][1]; //outgoing direction
			onesing.outgoing.entry[1] = onesing.Jacobian.entry[0][0]- r1;

			onesing.incoming.entry[0] = -onesing.Jacobian.entry[0][1]; //incoming direction
			onesing.incoming.entry[1] = onesing.Jacobian.entry[0][0]- r2;

			normalize(onesing.outgoing);
			normalize(onesing.incoming);
		}

		singularity_list.push_back(onesing);
	}
}


int get_singularity_type(double jacobian[2][2])
{
	double a,b,d,e;

	a = jacobian[0][0];
	b = jacobian[0][1];
	d = jacobian[1][0];
	e = jacobian[1][1];

	//need to store it as a part of the information of the elements or unknow singularities

	//Getting the eigenvalue here

	double A, B, C, delta;

	A = 1;
	B = -(a + e);
	C = (a * e - b * d);

	delta =  B*B - 4*A*C;

	if( delta >= 0)
	{
		i1 = i2 = 0;
		r1 = (-B + sqrt(delta))/2;
		r2 = (-B - sqrt(delta))/2;
	}
	else
	{
		r1 = r2 = -B/2.;
		i1 = sqrt(-delta)/2;
		i2 = -sqrt(-delta)/2;
	}

	////////////////////////////////////////////////////////////////
	if( r1 > 0 && r2 > 0 && i1 == 0 && i2 == 0)
		return 1; //it is a source

	else if( r1 < 0 && r2 < 0 && i1 == 0 && i2 == 0)
		return 2; //it is a sink

	else if( r1 * r2 < 0 && i1 == 0 && i2 == 0)
		return 3; //it is a saddle

	else if( r1 == 0 && r2 == 0 && i1 != 0 && i2 != 0)
	//else if( fabs(r1) <= 1e-4 && fabs(r2) <= 1e-4 && i1 != 0 && i2 != 0)
		return 4; //it is a center

	else if( r1 < 0 && r2 < 0 && i1 != 0 && i2 != 0)
		return 6; //it is an attracting focus

	else if( r1 > 0 && r2 > 0 && i1 != 0 && i2 != 0)
		return 7; //it is a repelling focus

	else
		return 0; //Unknow, there must be some error here!
}


/*----------------------------------------------------------------------*/
////Routines for matrix operations

/****************************************************
Maybe better for 2X2 and 3X3 matrix inversion
****************************************************/

double *MatrixOpp(double A[],int m,int n) //inverse 
{ 
    int i,j,x,y,k; 
    double *SP=NULL,*AB=NULL,*B=NULL,XX,*C; 
    SP=(double *)malloc(m*n*sizeof(double)); 
    AB=(double *)malloc(m*n*sizeof(double)); 
    B=(double *)malloc(m*n*sizeof(double)); 
    
    XX=Surplus(A,m,n); 
    XX=1/XX; 
    
    for(i=0;i<m;i++) 
	for(j=0;j<n;j++) 
	{
		for(k=0;k<m*n;k++) 
			B[k]=A[k]; 
			{
				for(x=0;x<n;x++) 
				B[i*n+x]=0; 
				for(y=0;y<m;y++) 
				B[m*y+j]=0; 
				B[i*n+j]=1; 
				SP[i*n+j]=Surplus(B,m,n); 
				AB[i*n+j]=XX*SP[i*n+j]; 
			} 
	} 

    C=MatrixInver(AB,m,n); 

	free(SP);
	free(AB);
	free(B);
    
    return C; 
} 
    
double * MatrixInver(double A[],int m,int n) //zhuanzhi
{ 
    int i,j; 
    double *B=NULL; 
    B=(double *)malloc(m*n*sizeof(double)); 
    
    for(i=0;i<n;i++) 
	for(j=0;j<m;j++) 
		B[i*m+j]=A[j*n+i]; 
    
    return B; 
} 
    
double Surplus(double A[],int m,int n) //hanglieshi
{ 
    
    int i,j,k,p,r; 
    double XX,temp=1,temp1=1,s=0,s1=0; 
    
    if(n==2) 
    {for(i=0;i<m;i++) 
    for(j=0;j<n;j++) 
    if((i+j)%2) temp1*=A[i*n+j]; 
    else temp*=A[i*n+j]; 
    XX=temp-temp1;} 
    else{ 
    for(k=0;k<n;k++) 
    {for(i=0,j=k;i<m,j<n;i++,j++) 
    temp*=A[i*n+j]; 
    if(m-i) 
    {for(p=m-i,r=m-1;p>0;p--,r--) 
    temp*=A[r*n+p-1];} 
    s+=temp; 
    temp=1; 
    } 
    
    for(k=n-1;k>=0;k--) 
    {for(i=0,j=k;i<m,j>=0;i++,j--) 
    temp1*=A[i*n+j]; 
    if(m-i) 
    {for(p=m-1,r=i;r<m;p--,r++) 
    temp1*=A[r*n+p];} 
    s1+=temp1; 
    temp1=1; 
    } 
    
    XX=s-s1;} 
    return XX; 
} 

#include "glut.h"

void DrawSolidCircle(double cx, double cy)
{
	int i;
	//double R = 0.0001;
	double R = 0.008/*/(2*zoom_factor)*/;
	double theta, deta ;
	deta = 2 * PI/49.;
	double x, y;
	theta = 0.;
	glBegin(GL_POLYGON);
	for(i = 0; i < 50; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}

void draw_hollow_circle(double cx, double cy)
{
	int i;
	//double R = 0.0001;
	double R = 0.0085/*/(2*zoom_factor)*/;
	double theta, deta ;
	deta = 2 * PI/99.;
	double x, y;
	theta = 0.;
	glLineWidth(2.);
	glBegin(GL_LINE_LOOP);
	for(i = 0; i < 100; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}


void SetColorByType(int type)
{
	if(type == SOURCE){
		glColor4f(0.,1., 0.,1);
	}

	else if(type == SINK){
		glColor4f(1.,0., 0.,1);
	}

	else if(type == SADDLE){
		glColor4f(0.,0.2,1., 1);
	}

	else if(type == CWCENTER){
		glColor4f(1.,0., 1.,1);
	}

	else if(type == CCWCENTER){
		glColor4f(0.3, 1., 1.,1);
	}

	else if(type == AFOCUS){
		//glColor4f(1., 1., 0., 1);
		//glColor4f(1.,0.5, 0, 1);
		glColor4f(1.,0, 0, 1);
	}

	else if(type == RFOCUS){
		//glColor4f(0, 1, 0.7, 1);
		glColor4f(0.,1., 0, 1);
	}
}

void 
display_singularities()
{
	int i;

	//glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
	glEnable(GL_BLEND);
	glEnable(GL_POLYGON_SMOOTH);

	for (i=0; i<singularity_list.size(); i++)
	{
		SetColorByType (singularity_list[i].type);
		DrawSolidCircle(singularity_list[i].gpos.entry[0], singularity_list[i].gpos.entry[1]);
		glColor3f(0, 0, 0);
		draw_hollow_circle(singularity_list[i].gpos.entry[0], singularity_list[i].gpos.entry[1]);
	}
}

/////////////////////////////////////////////////////////////////////////
////// trace out streamline
//
///***************************************************************
//New method to calculate the intersection of two line segments
//****************************************************************/
//
///* meaning of return value
// 0----Intersection dosn't exists                                                   
// 1----Intersection exists.                                                        
// 2----two line segments are parallel.                                         
// 3----two line segments are collinear, but not overlap.                      
// 4----two line segments are collinear, and share one same end point.       
// 5----two line segments are collinear, and overlap.                           
//*/    
//
//int Trajectory::get_intersection(double PointA[2], double PointB[2], 
//								 double PointC[2], double PointD[2], double t[2])
//{
//
//    double delta;
//    double t1,t2;
//    double a,b,c,d;
//    double xba,yba,xdc,ydc,xca,yca;
//
//    xba=PointB[0]-PointA[0];    yba=PointB[1]-PointA[1];
//    xdc=PointD[0]-PointC[0];    ydc=PointD[1]-PointC[1];
//    xca=PointC[0]-PointA[0];    yca=PointC[1]-PointA[1];
//
//    delta=xba*ydc-yba*xdc;
//    t1=xca*ydc-yca*xdc;
//    t2=xca*yba-yca*xba;
//
//    if(delta!=0)
//    {
//        t[0]=t1/delta;   t[1]=t2/delta;
//        /*two segments intersect (including intersect at end points)*/
//        //if ( t[0]<=1 && t[0]>=0 && t[1]<=1 && t[1]>=0 ) return 1;
//         if ( t[0]<=1 && (t[0]>=0 || fabs(t[0]) < 1e-8)
//			 && t[1]<=1 && (t[1]>=0 || fabs(t[1]) < 1e-8 )) 
//			 return 1;
//       else return 0; 
//    }
//
//    else
//    {       
//        /* AB & CD are parallel. */
//        if ( (t1!=0) && (t2!=0) ) return 2;
//
//        /* when AB & CD are collinear */
//
//        /*if AB isn't a vertical line segment, project to x-axis */
//        if(PointA[0]!=PointB[0])   
//        {
//            a=min(PointA[0],PointB[0]); b=max(PointA[0],PointB[0]);
//            c=min(PointC[0],PointD[0]); d=max(PointC[0],PointD[0]);
//
//            if ( (d<a) || (c>b) ) return  3;
//            else if( (d==a) || (c==b) ) return 4;  
//            else return 5;
//        }
//
//        else         /* if AB is a vertical line segment, project to y-axis */  
//        {
//
//            a=min(PointA[1],PointB[1]); b=max(PointA[1],PointB[1]);
//            c=min(PointC[1],PointD[1]); d=max(PointC[1],PointD[1]); 
//
//            if( (d<a) || (c>b) ) return  3;
//            else if( (d==a) || (c==b) ) return 4;
//            else return 5;
//        }
//    }
//}
//
//
//void  Trajectory::cross_boundary(double pre[2], double cur[2], 
//								 int face_id, double alpha[3], 
//								 int &which_edge, double t[2])
//{
//	Triangle *face = poly->tlist[face_id];
//
//	if(alpha[0] < 0 && alpha[1] < 0)
//	{
//		double p0[2] = {0, 0};
//		double p1[2] = {face->x2, face->y2};
//		if(get_intersection(pre, cur, p0, p1, t)==1)
//		{
//			which_edge = 1;
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//			return;
//		}
//		
//		//Calculate the intersection with edge v1v2
//
//		p0[0] = face->x1; p0[1] = 0;
//
//		if(get_intersection(pre, cur, p0, p1, t)==1)
//		{
//			which_edge = 0;
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//			return;
//		}
//	}
//
//	else if(alpha[0] < 0 && alpha[2] < 0)
//	{
//		//Calculate the intersection with edge v0v1
//		double p0[2] = {0, 0};
//		double p1[2] = {face->x1, 0};
//
//		if(get_intersection(pre, cur, p0, p1, t)==1)
//		{
//			which_edge = 2;
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//			return;
//		}
//		
//		//Calculate the intersection with edge v1v2
//		p0[0] = p1[0]; p0[1] = p1[1];
//		p1[0] = face->x2; p1[1] = face->y2;
//		
//		if(get_intersection(pre, cur, p0, p1, t)==1)
//		{
//			which_edge = 0;
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//			return;
//		}
//	}
//
//	else if(alpha[1] < 0 && alpha[2] < 0)
//	{
//		//Calculate the intersection with edge v0v1
//		double p0[2] = {0, 0};
//		double p1[2] = {face->x1, 0};
//
//		if(get_intersection(pre, cur, p0, p1, t)==1)
//		{
//			which_edge = 2;
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//			return;
//		}
//		
//		//Calculate the intersection with edge v0v2
//		p1[0] = face->x2; p1[1] = face->y2;
//	
//		if(get_intersection(pre, cur, p0, p1, t)==1)
//		{
//			which_edge = 1;
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//			return;
//		}
//	}
//
//	else if(alpha[0] < 0)
//	{
//		double p0[2] = {face->x1, 0};
//		double p1[2] = {face->x2, face->y2};
//
//		which_edge = 0;
//		get_intersection(pre, cur, p0, p1, t);
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//		return;
//	}
//
//	else if(alpha[1] < 0)
//	{
//		double p0[2] = {face->x2, face->y2};
//		double p1[2] = {0, 0};
//		which_edge = 1;
//		get_intersection(pre, cur, p0, p1, t);
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//		return;
//	}
//
//	else if(alpha[2] < 0)
//	{
//		double p0[2] = {0, 0};
//		double p1[2] = {face->x1, 0};
//		which_edge = 2;
//		get_intersection(pre, cur, p0, p1, t);
//			cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
//			cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
//		return;
//	}
//}
//
//
//void 
//Trajectory::pass_vertex(int vert_id, int &theone, int type)
//{
//	Vertex *vert = poly->vlist[vert_id];
//
//	int i;
//
//	double vang;
//	//Face *face;
//	Corner *c;
//	int NewTID = -1;
//
//    int orient;
//
//	////for 3D surfaces, we should not use the following angle directly!
//	vang = atan2(vert->vec.entry[1], vert->vec.entry[0]);
//
//	if(type == 1)
//		vang += PI;
//
//	if(vang < 0)
//		vang += 2*PI;
//	
//	else if(vang > 2*PI)
//		vang -= 2*PI;
//
//	if(vert->corners[0]->orient)
//		orient=1;
//	else
//		orient=0;
//
//	for( i = 0; i < vert->ncorners; i++)
//	{
//
//		c = vert->corners[i];
//
//		////first, we check which angle area the vector on the vertex will fall in
//		if(orient > 0)
//		{
//			if(c->BeginAng > c->EndAng)
//			{
//				if((vang >= c->BeginAng && vang < 2*PI)|| (vang <= c->EndAng && vang >= 0))
//				{
//					NewTID = i;
//					break;
//				}
//			}
//			else{
//				if(vang >= c->BeginAng && vang <= c->EndAng)
//				{
//					NewTID = i;
//					break;
//				}
//			}
//		}
//		else{
//			if(c->BeginAng < c->EndAng)
//			{
//				if((vang <= c->BeginAng && vang >= 0)|| (vang >= c->EndAng && vang < 2*PI))
//				{
//					NewTID = i;
//					break;
//				}
//			}
//			else{
//				if(vang <= c->BeginAng && vang >= c->EndAng)
//				{
//					NewTID = i;
//					break;
//				}
//			}
//		}
//	}
//
//	if(NewTID < 0) //avoid crash!
//	{
//		theone = -1;
//		return;
//	}
//	
//	theone = vert->corners[NewTID]->t;
//}
//
//
//void Trajectory::cross_a_vertex(int &face_id, double cur_p[2], double pre_p[2],int type, int &passornot)
//{
//	double vert[2];
//	int alpha_index = 0;
//	double max_alpha ;
//    int newtriangleid = 0;
//	Vertex* crossVert;
//
//	Triangle *face = poly->tlist[face_id];
//
//
//	////New way to get the possible crossing vertex
//	icVector2 test_dir;
//	test_dir.entry[0] = cur_p[0] ;
//	test_dir.entry[1] = cur_p[1] ;
//	max_alpha = length(test_dir);
//	alpha_index = 0;
//
//	test_dir.entry[0] = cur_p[0] - face->x1;
//	test_dir.entry[1] = cur_p[1] ;
//	if(length(test_dir) < max_alpha)
//	{
//		max_alpha = length(test_dir);
//	    alpha_index = 1;
//	}
//
//	test_dir.entry[0] = cur_p[0] - face->x2;
//	test_dir.entry[1] = cur_p[1] - face->y2;
//	if(length(test_dir) < max_alpha)
//	{
//		max_alpha = length(test_dir);
//	    alpha_index = 2;
//	}
//
//	crossVert = face->verts[alpha_index];
//
//	if(alpha_index == 0)
//	{
//		vert[0] = vert[1] = 0;
//	}
//
//	else if(alpha_index == 1)
//	{
//		vert[0] = face->x1;
//		vert[1] = 0;
//	}
//
//	else
//	{
//		vert[0] = face->x2;
//		vert[1] = face->y2;
//	}
//
//
//	double A, B, C;
//	A = pre_p[1] - cur_p[1];
//	B = cur_p[0] - pre_p[0];
//	C = (pre_p[0]*cur_p[1] - cur_p[0]*pre_p[1]);
//
//	double pending = A*vert[0] + B*vert[1] + C;
//
//	if(fabs(pending) <= 1e-8) ////passing the vertex
//	{
//		pass_vertex(crossVert->index, newtriangleid, type);
//		face_id = newtriangleid;	
//		passornot = alpha_index+1;
//		return;
//	}
//
//	passornot = 0;
//}
//
//
//void Trajectory::get_next_triangle(int &face_id, double pre[2], double cur[2], double param_t[2], int type, 
//					 int &PassVertornot, double alpha[3])
//{
//	int which_edge = -1;
//
//	int prev_face_id = face_id;
//
//	Triangle *prev_face = poly->tlist[face_id];
//
//	Vertex *vert = NULL;
//
//	PassVertornot = 0;
//	
//	////We should put pass vertex testing here before testing crossing edge
//	cross_a_vertex(face_id, cur, pre, type, PassVertornot);
//	if(PassVertornot > 0)
//	{
//		return ;
//	}
//
//	face_id = prev_face_id;  //////added on 06/08/05
//
//	cross_boundary(pre, cur, face_id, alpha, which_edge, param_t);
//
//
//	if(param_t[0] == -1 && param_t[1] == -1)
//	{
//		face_id = prev_face_id;   ////something wrong here
//		return;
//	}
//
//	////if not passing a vertex, judge which triangle it will enter later
//	pass_edge(face_id, which_edge);
//}
//
//
///********************************************************************
//If we have already judge which edge the curve will cross
//We can use the edge information to get next triangle
//********************************************************************/
//void Trajectory::pass_edge(int &face_id, int which_edge)
//{
//	////Using edge information to get next triangle
//
//	Triangle *face = poly->tlist[face_id];
//
//	Edge *theone ;
//	Vertex *v;
//
//	for(int i = 0; i < 3; i++)
//	{
//		v = face->verts[which_edge];
//
//		theone = face->edges[i];
//
//		if(theone->verts[0] != v && theone->verts[1] != v)
//			break;
//	}
//
//	Triangle *next_tri;
//	if(theone->tris[0] != face)
//		next_tri = theone->tris[0];
//	else
//		next_tri = theone->tris[1];
//
//	if (theone->ntris<2)
//		face_id = -1;
//	else
//		face_id = next_tri->index;
//
//}
//
//
//int Trajectory::trace_in_triangle(int &face_id, double globalp[3], int type, int &flag)
//{
//	int i;
//	double alpha[3];
//	double cur_point[2], pre_point[2];
//	double vert0[3];
//	icVector3 VP, globalv;
//
//	Triangle *face = poly->tlist[face_id];
//
//	Triangle *pre_f = face;
//	
//	////Temporary curve point array
//
//	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 100);
//	int NumPoints = 0;
//
//	////get the local coordinates of the starting point
//
//	VP.entry[0] = globalp[0] - face->verts[0]->x;
//	VP.entry[1] = globalp[1] - face->verts[0]->y;
//	VP.entry[2] = globalp[2] - face->verts[0]->z;
//
//	pre_point[0] = cur_point[0] = dot(VP, face->LX);
//	pre_point[1] = cur_point[1] = dot(VP, face->LY);
//
//	vert0[0] = face->verts[0]->x;   ////for update the global point
//	vert0[1] = face->verts[0]->y;
//	vert0[2] = face->verts[0]->z;
//
//	//globalface = face_id;
//
//	////////////////////////////////////////////////////
//    for(i = 0; i < 100; i++)
//	{
//		////1. calculate the barycentric coordinates for current point
//
//		face->get_barycentric_coordinates_loc(cur_point, alpha);
//
//		////2. if current point is inside current triangle
//		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
//			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
//			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
//		{
//			////store the point into the temp curve points list
//
//			temp_point_list[NumPoints].gpx = globalp[0];
//			temp_point_list[NumPoints].gpy = globalp[1];
//			temp_point_list[NumPoints].gpz = globalp[2];
//			temp_point_list[NumPoints].lpx = cur_point[0];
//			temp_point_list[NumPoints].lpy = cur_point[1];
//			temp_point_list[NumPoints].triangleid = face->index;  
//			NumPoints++;
//
//			pre_point[0] = cur_point[0];
//			pre_point[1] = cur_point[1];
//
//			//if(cal_next_point_euler1(pre_point, cur_point, face_id, alpha, type))
//			//if(cal_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
//			//if(cal_nextpt_RK4(pre_point, cur_point, face_id, alpha, type))
//			
//			if (get_next_pt(pre_point, cur_point, face_id, alpha, type, unsigned char(Integrator_opt)))
//			{
//				////update the global point
//
//				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;
//
//				//Get global coordinates of the point
//				globalp[0] = vert0[0] + globalv.entry[0];
//				globalp[1] = vert0[1] + globalv.entry[1];
//				globalp[2] = vert0[2] + globalv.entry[2];
//
//			}
//
//			else{  ////the curve reach a singularity
//				flag = 1;
//
//				////Store the record into global line segment array
//				//store_to_global_line_segs(temp_point_list, NumPoints);
//				//free(temp_point_list);
//				//return face_id;
//				break;
//			}
//		}
//
//		////3. if the point is out of current triangle
//		else{
//			double t[2] = {0.};
//
//			int PassVertornot = 0;
//            
//			get_next_triangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);
//
//			////update the globalpoint here (Modified on 01/30/07)
//			if(PassVertornot > 0)
//			{
//				if (face_id < 0) 
//				{
//					//free(temp_point_list);
//					//return -1;
//					face_id = -1;
//					break;
//				}
//
//				//we first need to know which vertex it is in the new triangle 01/30/07
//				Vertex* vertid = pre_f->verts[PassVertornot-1];
//				Triangle *cur_f = poly->tlist[face_id];
//				int vert_new = 0;
//				for(int k = 0; k < 3; k++)
//				{
//					if(cur_f->verts[k] == vertid)
//					{
//						vert_new = k;
//						break;
//					}
//				}
//
//				alpha[vert_new]=1-0.001;	
//				alpha[(vert_new+1)%3]=0.0005;
//				alpha[(vert_new+2)%3]=0.0005;
//
//
//				/* Get the new cur_point */
//				cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
//				cur_point[1] = alpha[2]*cur_f->y2;
//
//				globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;
//
//				globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
//				globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
//				globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
//				face=cur_f;
//			}
//
//			else{
//				//// transfer it to the global coordinates
//				globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;
//
//				globalp[0] = vert0[0] + globalv.entry[0];
//				globalp[1] = vert0[1] + globalv.entry[1];
//				globalp[2] = vert0[2] + globalv.entry[2];
//			}
//
//			////Add the intersection point to the temporary points' list
//			temp_point_list[NumPoints].gpx = globalp[0];
//			temp_point_list[NumPoints].gpy = globalp[1];
//			temp_point_list[NumPoints].gpz = globalp[2];
//			temp_point_list[NumPoints].lpx = cur_point[0];
//			temp_point_list[NumPoints].lpy = cur_point[1];
//			temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
//			NumPoints++;
//
//			//if(NumPoints > 1){
// 		//		////Store the record into global line segment array
//   //            store_to_global_line_segs(temp_point_list, NumPoints);
//			//}
//
//			//free(temp_point_list);
//			//return face_id;
//
//			break;
//		}
//
//	}
//
//	if(NumPoints > 1)
//		store_to_global_line_segs(temp_point_list, NumPoints);
//
//	free(temp_point_list);
//	return face_id;
//}
//
//
//void Trajectory::cal_one_traj(int face_id, double x, double y, double z, int type)
//{
//	int i;
//	int flag = 0;
//	double globalp[3];
//	int pre_face, cur_face;
//
//	pre_face = cur_face = face_id;
//
//	globalp[0] = x;   globalp[1] = y;    globalp[2] = z;
//	int NUMTRACETRIS = (int)sqrt((double)poly->ntris);
//
//	for(i = 0; i < 10*NUMTRACETRIS; i++)
//	{
//
//		if(cur_face == -1)
//		{
//			return;
//		}
//
//		pre_face = cur_face;
//		cur_face = trace_in_triangle(cur_face, globalp, type, flag); ////0 means always forward tracing here
//
//		if(flag > 0 || pre_face == cur_face ) 
//		{
//			return;
//		}
//
//	}
//}
//
//
//
//void Trajectory::store_to_global_line_segs(CurvePoints *temp, int num)
//{
//	int i;
//	icVector3 dis_vec;
//
//	LineSeg oneseg;
//
//	for (i=0; i<num; i++)
//	{
//		oneseg.gstart.entry[0] = temp[i].gpx;
//		oneseg.gstart.entry[1] = temp[i].gpy;
//		oneseg.gstart.entry[2] = temp[i].gpz;
//		oneseg.start.entry[0] = temp[i].lpx;
//		oneseg.start.entry[1] = temp[i].lpy;
//
//		oneseg.gend.entry[0] = temp[i+1].gpx;
//		oneseg.gend.entry[1] = temp[i+1].gpy;
//		oneseg.gend.entry[2] = temp[i+1].gpz;
//		oneseg.end.entry[0] = temp[i+1].lpx;
//		oneseg.end.entry[1] = temp[i+1].lpy;
//
//		////Use local coordinates to calculate the length
//		dis_vec.entry[0] = temp[i+1].lpx - temp[i].lpx;
//		dis_vec.entry[1] = temp[i+1].lpy - temp[i].lpy;
//		dis_vec.entry[2] = 0;
//
//		oneseg.length = length(dis_vec);
//
//		oneseg.Triangle_ID = temp[i].triangleid;
//
//		/*save to the global list*/
//		linesegs.push_back(oneseg);
//	}
//}
//
//
//bool 
//Trajectory::get_next_pt(double first[2], 
//						double second[2], 
//						int &face_id, 
//						double alpha[3], 
//						int type,
//						unsigned char opt)
//{
//	switch (opt)
//	{
//	case 0:
//		return cal_next_point_euler1(first, second, face_id, alpha, type);
//	case 1:
//		return cal_nextpt_2ndeuler(first, second, face_id, alpha, type);
//	case 2:
//		return cal_nextpt_RK4(first, second, face_id, alpha, type);
//
//	default:
//		return false;
//	}
//}
//
//
//bool Trajectory::cal_next_point_euler1(double first[2], double second[2], int &face_id, double alpha[3], int type)
//{
//
//	Triangle *t = poly->tlist[face_id];
//
//	icVector2 VecAtPoint = t->get_vector_at_pt(first);
//
//	if(length(VecAtPoint) < 1e-15 ) return false;
//	
//	//scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
//	normalize(VecAtPoint);
//
//	//VecAtPoint = 0.4*length(t->verts[0]->vec)*VecAtPoint;  // for the normalized vector field 11/19/2010
//	//VecAtPoint = 0.004*length(t->verts[0]->t_vec)*VecAtPoint;  // for the non-normalized vector field 11/19/2010
//
//	////Using first order Euler method to test 12/11/05
//	if(type == 0)
//	{
//		second[0] = first[0] + eulerstep_scalar * VecAtPoint.entry[0];
//		second[1] = first[1] + eulerstep_scalar * VecAtPoint.entry[1];
//	}
//	else
//	{
//		second[0] = first[0] - eulerstep_scalar * VecAtPoint.entry[0];
//		second[1] = first[1] - eulerstep_scalar * VecAtPoint.entry[1];
//	}
//
//	return true;
//}
//
//
//
//bool Trajectory::cal_nextpt_2ndeuler(double first[2], double second[2], int &face_id, double alpha[3], int type)
//{
//	////Using first order Euler method to get the next point
//
//	Triangle *t = poly->tlist[face_id];
//
//	icVector2 VecAtPoint = t->get_vector_at_pt(first);
//
//	if(length(VecAtPoint) < 1e-15) return false;
//
//	//scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
//	normalize(VecAtPoint);
//
//	//VecAtPoint = 0.4*length(t->verts[0]->vec)*VecAtPoint; /*evaluate in very small step now*/
//	//VecAtPoint = 0.004*length(t->verts[0]->t_vec)*VecAtPoint; /*evaluate in very small step now*/
//
//	double temp[2] = {0.};
//	double euler2nd_step = eulerstep_scalar; 
//	
//	if(type == 0)
//	{
//		temp[0] = first[0] + euler2nd_step/2.*VecAtPoint.entry[0];
//		temp[1] = first[1] + euler2nd_step/2*VecAtPoint.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - euler2nd_step/2*VecAtPoint.entry[0];
//		temp[1] = first[1] - euler2nd_step/2*VecAtPoint.entry[1];
//	}
//
//	////get the vector at next point
//	double alpha1[3] = {0.};
//	
//	t->get_barycentric_coordinates_loc(temp, alpha1);
//
//	if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1 
//			&& (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1 
//			&& (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)
//			
//	{
//		icVector2 VecAtPoint2 = t->get_vector_at_pt(temp);
//		normalize(VecAtPoint2);
//
//		icVector2 total_v;
//		total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
//		total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
//		normalize(total_v);
//		total_v = 0.4*length(t->verts[0]->vec)*total_v;          // for normalized vector field
//
//		if(type == 0)
//		{
//			second[0] = first[0] + euler2nd_step*total_v.entry[0];
//			second[1] = first[1] + euler2nd_step*total_v.entry[1];
//		}
//		else
//		{
//			second[0] = first[0] - euler2nd_step*total_v.entry[0];
//			second[1] = first[1] - euler2nd_step*total_v.entry[1];
//		}
//	}
//
//	else
//	{
//		if(type == 0)
//		{
//			second[0] = first[0] + euler2nd_step*VecAtPoint.entry[0];
//			second[1] = first[1] + euler2nd_step*VecAtPoint.entry[1];
//		}
//		else
//		{
//			second[0] = first[0] - euler2nd_step*VecAtPoint.entry[0];
//			second[1] = first[1] - euler2nd_step*VecAtPoint.entry[1];
//		}
//	}
//
//	return true;
//}
//

//bool 
//Trajectory::cal_nextpt_RK4(double first[2], double second[2], int &face_id, double alpha[3], int type)
//{
//	////Using 4th Runge-Kutta method to get the next point
//	icVector2 t_v;
//	double temp[2] = {0.};
//
//	
//	Triangle *t = poly->tlist[face_id];
//
//	/*compute K1*/
//	icVector2 V1 = t->get_vector_at_pt(first);
//
//	if(length(V1) < 1e-14) return false;
//
//
//	double RK4_step = eulerstep_scalar/2.;  // for normalized
//	//RK4_step = eulerstep_scalar/16.; // for diesel slice
//
//	t_v = V1;
//	normalize(t_v);
//	//t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/
//
//	if(type == 0)
//	{
//		temp[0] = first[0] + RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] + RK4_step/2*t_v.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] - RK4_step/2*t_v.entry[1];
//	}
//
//	/*compute K2*/
//	////Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
//	////icVector2 V2 = GetVectorAtPoints(face_id, alpha);
//	//local_To_global(face_id, temp, t_gp);
//	//object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
//	//icVector2 V2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
//	icVector2 V2 = t->get_vector_at_pt(temp);
//	t_v = V2;
//	normalize(t_v);
//	//t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/
//	
//	if(type == 0)
//	{
//		temp[0] = first[0] + RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] + RK4_step/2*t_v.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] - RK4_step/2*t_v.entry[1];
//	}
//
//	/*compute K3*/
//	//local_To_global(face_id, temp, t_gp);
//	//object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
//	//icVector2 V3 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
//	icVector2 V3 = t->get_vector_at_pt(temp);
//	t_v = V3;
//	normalize(t_v);
//	//t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/
//	
//	if(type == 0)
//	{
//		temp[0] = first[0] + RK4_step*t_v.entry[0];
//		temp[1] = first[1] + RK4_step*t_v.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - RK4_step*t_v.entry[0];
//		temp[1] = first[1] - RK4_step*t_v.entry[1];
//	}
//	
//	/*compute K4*/
//	//local_To_global(face_id, temp, t_gp);
//	//object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
//	//icVector2 V4 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
//	icVector2 V4 = t->get_vector_at_pt(temp);
//
//	icVector2 total_v = 1./6.*(V1+2*V2+2*V3+V4);
//	t_v = total_v;
//	normalize(t_v);
//	//t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/
//
//	if(type == 0)
//	{
//		second[0] = first[0] + RK4_step*t_v.entry[0];
//		second[1] = first[1] + RK4_step*t_v.entry[1];
//	}
//
//	else
//	{
//		second[0] = first[0] - RK4_step*t_v.entry[0];
//		second[1] = first[1] - RK4_step*t_v.entry[1];
//	}
//
//	return true;
//}
//
//const double SEPARATRIXSTEP = 0.09;
//void Trajectory::cal_startpt_sep(int triangleID, icVector3 sep_vector, double saddle_ce[3], double newpos[3])
//{
//	int count = 2;
//	Triangle *t = poly->tlist[triangleID];
//
//	////Get the first position of the beginning point
//	newpos[0] = saddle_ce[0] + SEPARATRIXSTEP * sep_vector.entry[0];
//	newpos[1] = saddle_ce[1] + SEPARATRIXSTEP * sep_vector.entry[1];
//	newpos[2] = saddle_ce[2] + SEPARATRIXSTEP * sep_vector.entry[2];
//
//	while(count <= 500)
//	{
//		if(t->fall_in_the_tri(newpos))
//			return;
//	    
//		newpos[0] = saddle_ce[0] + SEPARATRIXSTEP/(/*2**/count) * sep_vector.entry[0];
//		newpos[1] = saddle_ce[1] + SEPARATRIXSTEP/(/*2**/count) * sep_vector.entry[1];
//	    newpos[2] = saddle_ce[2] + SEPARATRIXSTEP/(/*2**/count) * sep_vector.entry[2];
//
//		count+=1;
//	}
//
//
//}

////////////////////////////////////////////////////////////////////////
////  Calculate separatrices starting from and ending at saddles
///

//void  cal_seps()
//{
//	int i;
//
//	icVector3 g_out, g_in, sing_pos;
//	Triangle *t;
//	double newpos[3]={0.};
//	int start_tri;
//	int flag = -1;
//
//	trajs.clear();
//	incoming_seps.clear();
//	outgoing_seps.clear();
//
//	for (i=0; i<singularity_list.size(); i++)
//	{
//		if (singularity_list[i].type != SADDLE)
//			continue;
//
//		t = poly->tlist[singularity_list[i].TriangleID];
//		
//		g_out = singularity_list[i].outgoing.entry[0] * t->LX
//			+ singularity_list[i].outgoing.entry[1] * t->LY;
//
//		g_in = singularity_list[i].incoming.entry[0] * t->LX
//			+ singularity_list[i].incoming.entry[1] * t->LY;
//
//		/*1. Positive outgoing separatrix */
//		Trajectory traj;
//		traj.eulerstep_scalar = t->edges[0]->length/10;
//
//		sing_pos.set(singularity_list[i].gpos.entry[0], singularity_list[i].gpos.entry[1], 0);
//		start_tri = singularity_list[i].TriangleID;
//		traj.cal_startpt_sep(start_tri, g_out, sing_pos.entry, newpos);
//		traj.cal_one_traj(start_tri, newpos[0], newpos[1], newpos[2], 0);
//		traj.index = trajs.size()-1;
//		trajs.push_back(traj);
//		outgoing_seps.push_back(trajs.size()-1);
//
//		/*2. Negative outgoing separatrix */
//		traj.linesegs.clear();
//		sing_pos.set(singularity_list[i].gpos.entry[0], singularity_list[i].gpos.entry[1], 0);
//		start_tri = singularity_list[i].TriangleID;
//		traj.cal_startpt_sep(start_tri, -g_out, sing_pos.entry, newpos);
//		traj.cal_one_traj(start_tri, newpos[0], newpos[1], newpos[2], 0);
//
//		traj.index = trajs.size()-1;
//		trajs.push_back(traj);
//		outgoing_seps.push_back(trajs.size()-1);
//
//
//		/*3. Positive incoming separatrix */
//		traj.linesegs.clear();
//		sing_pos.set(singularity_list[i].gpos.entry[0], singularity_list[i].gpos.entry[1], 0);
//		start_tri = singularity_list[i].TriangleID;
//		traj.cal_startpt_sep(start_tri, g_in, sing_pos.entry, newpos);
//		traj.cal_one_traj(start_tri, newpos[0], newpos[1], newpos[2], 1);
//
//		traj.index = trajs.size()-1;
//		trajs.push_back(traj);
//		incoming_seps.push_back(trajs.size()-1);
//
//		/*4. Negative outgoing separatrix */
//		traj.linesegs.clear();
//		sing_pos.set(singularity_list[i].gpos.entry[0], singularity_list[i].gpos.entry[1], 0);
//		start_tri = singularity_list[i].TriangleID;
//		traj.cal_startpt_sep(start_tri, -g_in, sing_pos.entry, newpos);
//		traj.cal_one_traj(start_tri, newpos[0], newpos[1], newpos[2], 1);
//
//		traj.index = trajs.size()-1;
//		trajs.push_back(traj);
//		incoming_seps.push_back(trajs.size()-1);
//	}
//}


//void 
//display_one_streamline(Trajectory traj)
//{
//	int i;
//
//	glBegin(GL_LINE_STRIP);
//	for (i=0; i<traj.linesegs.size(); i++)
//	{
//		glVertex2d (traj.linesegs[i].gstart.entry[0], traj.linesegs[i].gstart.entry[1]);
//	}
//	glEnd();
//}


//void 
//display_separatrices()
//{
//	int i;
//
//	glLineWidth(2.);
//
//	glColor3f (1, 0, 0);  // connect to sinks
//	for (i=0; i<outgoing_seps.size(); i++)
//	{
//		display_one_streamline(trajs[outgoing_seps[i]]);
//	}
//
//	glColor3f (0, 1, 0);  // connect to sources
//	for (i=0; i<incoming_seps.size(); i++)
//		display_one_streamline(trajs[incoming_seps[i]]);
//}