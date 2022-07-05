

// An image space based Morse decomposition

#ifndef __MORSEDECOMP_H__
#define __MORSEDECOMP_H__

#include "directed_graph.h"
#include <stack>
#include <vector>
using namespace std;

#define RANDOMINTEGRATION

typedef unsigned char vec3ch[3];
union MyColor {
   vec3ch val;
   struct {
      unsigned char r;
      unsigned char g;
      unsigned char b;
   };
};

class Img_MorseDecomp
{
public:
	DirGraph *dg;
	vector<SCComponent> scclist;
	int integral_steps;
	//float x_interval, y_interval;
	int XDIM, YDIM;  // XDIM - how many columns! YDIM - how many rows!
	int nsamples;
	MyColor **vec_img;
	bool **fixedPt_pixels;

	vector<vector<int>> end_pixels; //store all the end pixel positions for each starting pixel

	float min_vx, max_vx, min_vy, max_vy;
	double predict_stepsize;

	int valid_MS;

	int rand_integrator, rand_stepsize;

	Img_MorseDecomp(int X=0, int Y=0, unsigned char ***in_vec_img = NULL)
		: XDIM(X), YDIM(Y)
	{
		dg = NULL;
		scclist.clear();
		integral_steps = 1;
		predict_stepsize = 1.;

		//x_interval = 1./(XDIM-1);
		//y_interval = 1./(YDIM-1);

		if (in_vec_img == NULL)
		{
			vec_img = NULL;
		}

		else
		{
			int i, j;

			vec_img = new MyColor *[Y];

			for (i=0; i<Y; i++)
			{
				vec_img[i] = new MyColor[X];
			}

			for (i=0; i<Y; i++)
			{
				for (j=0; j<X; j++)
				{
					vec_img[i][j].val[0] = in_vec_img[i][j][0];
					vec_img[i][j].val[1] = in_vec_img[i][j][1];
					vec_img[i][j].val[2] = in_vec_img[i][j][2];
				}
			}
		}

		nsamples = 25; //16; //9; //4; //10  // we can use FTLE results to guide the sampling distribution!!!
		valid_MS = 0;
		rand_integrator = 1;
		rand_stepsize = 2;

		if (X>0 && Y>0)
		{
			fixedPt_pixels = new bool *[Y];

			int i;
			for (i=0; i<Y; i++)
				fixedPt_pixels[i] = new bool[X];
			//initialize the end_pixels vector
			end_pixels.resize(XDIM*YDIM);

			//for (i=0; i<X*Y; i++)
			//	end_pixels[i].resize(nsamples);
		}

		else
			fixedPt_pixels = NULL;


	}

	~Img_MorseDecomp()
	{
		if (dg != NULL)
			delete dg;
		dg = NULL;

		if (vec_img != NULL)
			delete [] vec_img;
		vec_img = NULL;

		if (fixedPt_pixels != NULL)
			delete [] fixedPt_pixels;
		fixedPt_pixels = NULL;

		scclist.clear();

		end_pixels.clear();
	}

	void finalize()
	{
		if (dg != NULL)
			delete dg;
		dg = NULL;

		if (vec_img != NULL)
			delete [] vec_img;
		vec_img = NULL;

		if (fixedPt_pixels != NULL)
			delete [] fixedPt_pixels;
		fixedPt_pixels = NULL;

		scclist.clear();

		end_pixels.clear();
	}
	
	void init_graph();

	void build_DG_geometry (); // construct the flow combinatorialization using a geometric-based method (need to think of it for the image-space method)

	void build_DG_tau(int steps);  // we use the number of the integration steps to replace time

	void trace_all_pixels (int nsamples, int steps, bool forward_backward);

	void trace_one_pixel (int row, int column, int nsamples, int steps, bool forward_backward);

	void generate_samples_in_one_pixel (int row, int column, int nsamples, vector<pair<float, float>> &samples); // compute n sub-samples within a pixel
	
	void generate_samples_in_one_pixel_regular (int row, int column, int nsamples, vector<pair<float, float>> &samples); // regularly (or evenly) sample the pixel

	void generate_samples_in_one_pixel_random (int row, int column, int nsamples, vector<pair<float, float>> &samples); // randomly sample the pixel

	void get_vec_at_regular_grid_image_based(float x, float y, /*int xdim, int ydim,*/ float &vx, float &vy);
	bool get_vec_at_regular_grid_image_based2(double xy[2], double vxy[2], bool backward);

	double bilinear_interpolate(double a, double b, 
						  double f00, double f01, double f10, double f11);

	void image_space_trace_forward (float sx, float sy, int nsteps, pair<int, int> &end_pixel);

	void image_space_trace_backward (float sx, float sy, int nsteps, pair<int, int> &end_pixel);

	void image_space_trace_forward_Euler (float sx, float sy, int nsteps, pair<int, int> &end_pixel);
	void image_space_trace_backward_Euler (float sx, float sy, int nsteps, pair<int, int> &end_pixel);
	
	void image_space_trace_forward_RK4 (float sx, float sy, int nsteps, pair<int, int> &end_pixel);
	void image_space_trace_backward_RK4 (float sx, float sy, int nsteps, pair<int, int> &end_pixel);

	bool get_nextpt_RK23_quad(double first[2], double second[2], double offset[2], bool type);
	bool get_nextpt_RK45_quad(double first[2], double second[2], double offset[2], bool type);
	
	bool RK23_2d(double pre_p[2], double next_p[2], double offset[2], double &hstep_loc, double &hnext,
			 double eps, double &eps_did, bool backward /*bool deriv(double xy[2], double vxy[2], bool backward) void deriv(double cur_p[2], double vec[2])*/
			 );
	bool RK45_2d(double pre_p[2], double next_p[2], double offset[2], double &hstep_loc, double &hnext,
			 double eps, double &eps_did, bool backward //, bool deriv(double xy[2], int xdim, int ydim, double vxy[2], bool backward) /*void deriv(double cur_p[2], double vec[2])*/
			 );

	void build_edges();

	void build_one_edge(int node1, int node2);

	//void set_vec_img(unsigned char ***in_vec_img, float xmin, float xmax, float ymin, float ymax);
	void set_vec_img(unsigned char *in_vec_img, float xmin, float xmax, float ymin, float ymax);

	//void set_fixedPt_mask(bool **in_fixedPt);
	void set_fixedPt_mask(bool *in_fixedPt);

	//void init_all_edges();

	void cal_Img_MorseDecomp(int nsteps);

	void build_SCCElemList();
	void mark_all_valid_SCCS();

	//classify Morse sets
	int get_Type_MorseSet(int scc_id);


	// visualization
	void mark_valid_SCC_pixels (unsigned char *tex, unsigned char *out); // this is a test code for visualizing the obtained Morse sets
};

#endif //__MORSEDECOMP_H__