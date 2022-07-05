
////VFAnalysis.h

/*This file provide routines fulfilling the analysis of current vector field (generally based on the vector field
stored in 2D local frame)

The inputs of this module are underneath mesh, vector field stored at each vertex
The output is the results, such as singularities, limit cycles, separatrices

Most of the calculations are performed here, especially local tracing calculation for trajectory
*/
#pragma once

#include "icVector.H"
#include "icMatrix.H"

#include <vector>
using namespace std;


enum SINGULARElEM{
	NOTHING,
	SOURCE,
	SINK,
	SADDLE,
	CWCENTER,
	CCWCENTER,
	AFOCUS,
	RFOCUS,
	REGULAR
};

class Singularity{
public:
	bool connected;
	unsigned char nnodes;
	unsigned char  num_connect_cycles;
	unsigned char type;                    //the type of the singularity
	int  index;                            //the index of the singularity in the list
	int  TriangleID;                       //the triangle containing the singularity
	//int *connect_cycles;
	int separatrices;                      //the index of the separatrix group belongs to the saddle 10/13/05
	int  node_index;                      //the index of the node in ECG graph
	
	icVector2  gpos;                       //global coordinates of singularity
	//icVector2  lpos;                       //local coordinates of singularity
	icVector2  incoming, outgoing;       //the local e-vectors of the singularity (esp. saddle)
	icMatrix2x2 Jacobian;                  //the jacobian associated with the fixed point

	/*-------for ECG-graph --------if I can find a way to build the graph during analysis, we can safely remove these variables */
	//Singularity **connect_singulars;         //singularities connecting with the singularity
	//unsigned char num_connect_singulars;
	//LimitCycle **connect_cycles;             //limit cycles connecting with the singularity
	//ECG_Node **nodes;                     //the list of nodes connecting with the singularity in ECG

	/*member functions for updating the connection information involved this fixed point*/
	void update_list_to_PO(int limitcycle);

}; // end of Singularity class


/* Line segment data structure */
typedef struct LineSeg{
	int Triangle_ID;            //which triangle this line segment locates in
	icVector2 start, end;    //local coordinates for start and end points
	icVector3 gstart, gend;  //global coordinates for start and end points
	double length;              //we may need to store the length of the line segment
} LineSeg;


/*For temporary trajetory*/ 
typedef struct CurvePoints{
	int triangleid;
	double lpx, lpy;
	double gpx, gpy, gpz;
	double length;
}CurvePoints;


class Trajectory{
public:
	int index;
	//int  nlinesegs;
	//int  curMaxNumLinesegs;
	//LineSeg *linesegs;
	std::vector<LineSeg> linesegs;

	int saddleID;          /*which saddle this trajectory belongs to*/

	double eulerstep_scalar;

	/*Construct the trajectory*/
	Trajectory(int index, int curMaxNum = 200);

	Trajectory(){}

	~Trajectory()
	{
		linesegs.clear();
	}

	//get the length of the trajectory
	double get_length();

	//get any line segment according to the input index
	LineSeg *get_line_seg(int index);

	/*------------ routines for calculating streamlines -------------------------------*/
	//all the tracing calculation should be fulfilled here
	bool cal_next_point_euler1(double first[2], double second[2], int &face_id, double alpha[3], int type);
	bool cal_nextpt_2ndeuler(double first[2], double second[2], int &face_id, double alpha[3], int type);
	bool cal_nextpt_RK4(double first[2], double second[2], int &face_id, double alpha[3], int type);

	void store_to_global_line_segs(CurvePoints *temp, int num);
	void local_To_global(int faceid, double locpos[2], icVector3 &glpos);
	int trace_in_triangle(int &face_id, double globalp[3], int type, int &flag);
	void get_next_triangle(int &face_id, double pre[2], double cur[2], double param_t[2], int type, 
					 int &PassVertornot, double alpha[3]);
	void cross_a_vertex(int &face_id, double cur_p[2], double pre_p[2],int type, int &passornot);
	void pass_edge(int &face_id, int which_edge);
	void pass_vertex(int vert_id, int &theone, int type);
    void cross_boundary(double pre[2], double cur[2], 
								 int face_id, double alpha[3], 
								 int &which_edge, double t[2]);
	int get_intersection(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2]);

	void cal_one_traj(int face_id, double x, double y, double z, int type);
	
	bool get_next_pt(double first[2], double second[2], int &face_id, double alpha[3], int type, unsigned char opt);

	void cal_startpt_sep(int triangleID, icVector3 sep_vector, double saddle_ce[3], double newpos[3]);

}; //end of Trajectory class



////
//The followings are some possible useful routines from previous program

////singularities extraction
void capture_singular_tris(void);
void compute_fixedpt_info(void);
int GetSingType(int MarkTriangleID);

int get_singularity_type(double jacobian[2][2]);
void get_loc_Jacobian(int tri, double Jacobian[2][2], double &c, double &f);
int get_singularity_type_loc(int tri, double Jacobian[2][2], double x_loc, double y_loc, double alpha[3]);
void compute_fixedpt_info(void);

////functions for calculate trajectory
void CalSingleTrajectory(double x, double y);
void CalSingleSeparatrix(double x, double y, int inout);
void CalSeparatrices();

//// display singularity

void  display_singularities();

////limit cycle detection

//// compute separatrices

void  cal_seps();

void  display_separatrices();







//////functions for matrix inversion
////These routines should be moved to another library
double *MatrixOpp(double A[],int m,int n); //inverse 
double *MatrixInver(double A[],int m,int n); //transpose
double Surplus(double A[],int m,int n); //


double cal_determinant2x2(double a[][2]);
bool cal_inverse2x2(double a[][2], double inverse_a[][2]);
//void cal_eigenvectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2]);
