
// The implementation of the image-space Morse decomposition

#include "Img_MorseDecomp.h"

#include <math.h>
#include "icVector.H"
#include <omp.h>


void init_graph();

void Img_MorseDecomp::build_DG_geometry () // construct the flow combinatorialization using a geometric-based method (need to think of it for the image-space method)
{
}

void Img_MorseDecomp::build_DG_tau(int steps)  // we use the number of the integration steps to replace time
{
	//forward
	trace_all_pixels (nsamples, steps, false);

//	int i, j;
//#pragma omp parallel for 
//	for (i=0; i<XDIM; i++)
//	{
//		for (j=0; j<YDIM; j++)
//		{
//			trace_one_pixel (i, j, nsamples, steps, false);
//		}
//	}

	//backward
	trace_all_pixels (nsamples, steps, true);

	//build edges
	//build_edges();
}


void Img_MorseDecomp::trace_all_pixels (int nsamples, int steps, bool forward_backward)
{
	int i, j;

	//omp_set_num_threads (100);

//#pragma omp parallel for 
	for (i=0; i<XDIM; i++)
	{
		for (j=0; j<YDIM; j++)
		{
			trace_one_pixel (i, j, nsamples, steps, forward_backward);
		}
	}
	//for (i=0; i<XDIM*YDIM; i++)
	//{
	//	int row = i/XDIM;
	//	int col = i%XDIM;
	//	trace_one_pixel (row, col, nsamples, steps, forward_backward);
	//}
}

bool repeated_int(std::vector<int> int_array, int elem)
{
	for (int i=0; i<int_array.size(); i++)
	{
		if (int_array[i] == elem)
			return true;
	}

	return false;
}

void Img_MorseDecomp::trace_one_pixel (int row, int column, int nsamples, int steps, bool forward_backward)
{
	// we trace from the sampled positions within a pixel

	// get samples
	vector<pair<float, float>> samples;

	generate_samples_in_one_pixel (row, column, nsamples, samples);

	int i;
	
	int start_pixel_id = row*XDIM+column;

	if (!forward_backward)
	{
		for (i=0; i<samples.size(); i++)
		{
			pair<int, int> end_pixel;

#ifdef  RANDOMINTEGRATION
			if (rand_integrator == 0)
				image_space_trace_forward_Euler (samples[i].first, samples[i].second,
				   steps, end_pixel);
			else if (rand_integrator == 1)
				image_space_trace_forward (samples[i].first, samples[i].second,
				   steps, end_pixel);
			else
				image_space_trace_forward_RK4 (samples[i].first, samples[i].second,
				   steps, end_pixel);
#else
			image_space_trace_forward (samples[i].first, samples[i].second,
				steps, end_pixel);
#endif
			int end_pixel_id = end_pixel.first*XDIM+end_pixel.second;

			//#pragma omp critical
			//{
			//if (!repeated_int(end_pixels[start_pixel_id], end_pixel_id))
			//    end_pixels[start_pixel_id].push_back(end_pixel_id);
			//}
		
			#pragma omp critical
			{
			if (!dg->is_repeated_edge_2(start_pixel_id, end_pixel_id))
				build_one_edge(start_pixel_id, end_pixel_id);
			}
		}
	}

	else
	{
		for (i=0; i<samples.size(); i++)
		{
			pair<int, int> end_pixel;

#ifdef RANDOMINTEGRATION
			if (rand_integrator == 0)
				image_space_trace_backward_Euler (samples[i].first, samples[i].second,
				   steps, end_pixel);
			else if (rand_integrator == 1)
				image_space_trace_backward (samples[i].first, samples[i].second,
				   steps, end_pixel);
			else
				image_space_trace_backward_RK4 (samples[i].first, samples[i].second,
				   steps, end_pixel);
#else
			image_space_trace_backward (samples[i].first, samples[i].second,
				steps, end_pixel);
#endif
			int end_pixel_id = end_pixel.first*XDIM+end_pixel.second;
			
			//#pragma omp critical	
			//{
			//if (!repeated_int(end_pixels[end_pixel_id], start_pixel_id))
			//	end_pixels[end_pixel_id].push_back(start_pixel_id);
			//}

			#pragma omp critical
			{
			if (!dg->is_repeated_edge_2(end_pixel_id, start_pixel_id))
				build_one_edge(end_pixel_id, start_pixel_id);
			}
			
		}
	}

	//samples.resize(0);
	samples.clear();
	//samples.shrink_to_fit();
	//vector<pair<float, float>> tmp;
	vector<pair<float, float>>().swap(samples);
}

void Img_MorseDecomp::generate_samples_in_one_pixel_regular (int row, int column, int nsamples, vector<pair<float, float>> &samples)
{
	int nseeds = sqrt((float) nsamples);

	double x_sub_interval = /*x_interval*/1./(nseeds);
	double y_sub_interval = /*y_interval*/1./(nseeds);
	
	int i, j;

	float start_x = column /** x_interval*/ + x_sub_interval/2.;
	float start_y = row /** y_interval*/ + y_sub_interval/2.;

//#pragma omp parallel for nowait  (shared variable "samples")
	for (i=0; i<nseeds; i++)
	{
		for (j=0; j<nseeds; j++)
		{
			pair <float, float> sample;
			sample.first = start_x + i * x_sub_interval;
			sample.second = start_y + j * y_sub_interval;

			samples.push_back(sample);
		}
	}
}

void Img_MorseDecomp::generate_samples_in_one_pixel_random (int row, int column, int nsamples, vector<pair<float, float>> &samples)
{

	// we try two random sample strategies here

	int i;

	// since we are in image-space, the coordinates are integer

	float start_x = column  /** x_interval*/;
	float start_y = row /** y_interval*/ ;

//#pragma omp parallel for nowait  (shared variable "samples")
	for (i=0; i<nsamples; i++)
	{
		float x_off = (float)rand()/(float)INT_MAX;
		float y_off = (float)rand()/(float)INT_MAX;

		pair <float, float> sample;
		sample.first = start_x + /*x_interval **/ x_off;
		sample.second = start_y + /*y_interval * */y_off;

		samples.push_back(sample);
	}

	// strategy 2
	//int nseeds = sqrt((float) nsamples);
	//double x_sub_interval = x_interval/(nseeds-1);
	//double y_sub_interval = y_interval/(nseeds-1);
}


void Img_MorseDecomp::generate_samples_in_one_pixel (int row, int column, int nsamples, vector<pair<float, float>> &samples)
{
	float check = sqrt((float)nsamples);

	if ((int) check == check)
	{
		generate_samples_in_one_pixel_regular(row, column, nsamples, samples);
	}

	else
	{
		generate_samples_in_one_pixel_random (row, column, nsamples, samples);
	}
}


double Img_MorseDecomp::bilinear_interpolate(double a, double b, 
						  double f00, double f01, double f10, double f11)
{
	return (f00*(1-a)*(1-b)+f10*a*(1-b)+f01*(1-a)*b+f11*a*b);
}


void Img_MorseDecomp::get_vec_at_regular_grid_image_based(float x, float y,/* int xdim, int ydim,*/ float &vx, float &vy)
{
	int lrow = (int) y;
	int lcol = (int) x;

	if (lrow<0 || lrow>=YDIM-1 || lcol<0 || lcol>=XDIM-1) // reaching the boundary
	{
		vx = vy = 0.;
		return;
	}

	// Do we need to consider the center of the pixel instead, which is (lcol+.5, lrow+.5)
	double a = (x - lcol);
	double b = (y - lrow);

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
	double f00 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol].val[0]/255.;
	double f01 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol].val[0]/255.;
	double f10 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol+1].val[0]/255.;
	double f11 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol+1].val[0]/255.;
	vx = bilinear_interpolate(a, b, f00, f01, f10, f11);

	f00 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol].val[1]/255.;
	f01 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol].val[1]/255.;
	f10 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol+1].val[1]/255.;
	f11 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol+1].val[1]/255.;
	vy = bilinear_interpolate(a, b, f00, f01, f10, f11);
}


bool Img_MorseDecomp::get_vec_at_regular_grid_image_based2(double xy[2], double vxy[2], bool backward)
{
	double x = xy[0];
	double y = xy[1];
	//int lrow = (int) y;
	//int lcol = (int) x;

	int lrow = (int) y;
	int lcol = (int) x;

	if (lrow<0 || lrow>=YDIM-1 || lcol<0 || lcol>=XDIM-1) // reaching the boundary
	{
		vxy[0] = vxy[1] = 0.;
		return false;
	}

	// Do we need to consider the center of the pixel instead, which is (lcol+.5, lrow+.5)
	double a = (x - lcol);
	double b = (y - lrow);

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
	double f00 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol].val[0]/255.;
	double f01 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol].val[0]/255.;
	double f10 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol+1].val[0]/255.;
	double f11 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol+1].val[0]/255.;
	vxy[0] = bilinear_interpolate(a, b, f00, f01, f10, f11);


	f00 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol].val[1]/255.;
	f01 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol].val[1]/255.;
	f10 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol+1].val[1]/255.;
	f11 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol+1].val[1]/255.;
	vxy[1] = bilinear_interpolate(a, b, f00, f01, f10, f11);

	if (backward) {vxy[0] =-vxy[0]; vxy[1] = -vxy[1];}
	return true;
}

extern bool 
RK2_integrate_multi_steps(double trace_gp[2], int N, bool forward_backward);

extern bool 
Euler_integrate_multi_steps(double trace_gp[2], int N, bool forward_backward);

extern bool 
RK4_integrate_multi_steps(double trace_gp[2], int N, bool forward_backward);


void Img_MorseDecomp::image_space_trace_forward (float sx, float sy, int nsteps, pair<int, int> &end_pixel)
{
	float vec[2], vec2[2];

	float x = sx, y = sy;

	//get_vec_at_regular_grid_image_based (x, y, vec[0], vec[1]);

	int i = (int) sy;  //which_row;
	int j = (int) sx;  //which_column;

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	//double pre_ang, cur_ang, ang_diff;
	//pre_ang = atan2(vec[1], vec[0]);

	//closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (step = 0; step < nsteps; step ++)
	{
		//if (fixedPt_pixels[i][j]) break;
		//get_nextpt_RK23_quad(trace_gp, trace_gp, offset, false);

		if (!RK2_integrate_multi_steps(trace_gp, rand_stepsize, false)) break;

//#ifdef RANDOMINTEGRATION
//		if (rand_integrator==0)
//		{
//			if (!Euler_integrate_multi_steps(trace_gp, rand_stepsize, false)) break;
//		}
//		else if (rand_integrator == 1)
//		{
//			if (!RK2_integrate_multi_steps(trace_gp, rand_stepsize, false)) break;
//		}
//		else
//		{
//			if (!RK4_integrate_multi_steps(trace_gp, rand_stepsize, false)) break;
//		}
//#else
//		if (!RK2_integrate_multi_steps(trace_gp, 2, false)) break;
//#endif
		int jj = (int) trace_gp[0];
		int ii = (int) trace_gp[1];

		if (ii<=0 || ii>=YDIM-1 || jj<=0 || jj>=XDIM-1)
			break;

		j = jj;
		i = ii;
	}

	end_pixel.first = i;
	end_pixel.second = j;
}

void Img_MorseDecomp::image_space_trace_backward (float sx, float sy, int nsteps, pair<int, int> &end_pixel)
{

	float vec[2], vec2[2];

	float x = sx, y = sy;

	//get_vec_at_regular_grid_image_based (x, y, vec[0], vec[1]);

	int i = (int) sy;  //which_row;
	int j = (int) sx;  //which_column;

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	//double pre_ang, cur_ang, ang_diff;
	//pre_ang = atan2(vec[1], vec[0]);

	//closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (step = 0; step < nsteps; step ++)
	{
		//if (fixedPt_pixels[i][j]) break;
		//get_nextpt_RK23_quad(trace_gp, trace_gp, offset, true);
		
		if (!RK2_integrate_multi_steps(trace_gp, rand_stepsize, true)) break;


//#ifdef RANDOMINTEGRATION
//		if (rand_integrator==0)
//		{
//			if (!Euler_integrate_multi_steps(trace_gp, rand_stepsize, true)) break;
//		}
//		else if (rand_integrator == 1)
//		{
//			if (!RK2_integrate_multi_steps(trace_gp, rand_stepsize, true)) break;
//		}
//		else
//		{
//			if (!RK4_integrate_multi_steps(trace_gp, rand_stepsize, true)) break;
//		}
//#else
//		if (!RK2_integrate_multi_steps(trace_gp, 2, true)) break;
//#endif
		
		//j = (int) trace_gp[0];
		//i = (int) trace_gp[1];

		//if (i<=0 || i>=YDIM-1 || j<=0 || j>=XDIM-1)
		//	break;

		int jj = (int) trace_gp[0];
		int ii = (int) trace_gp[1];

		if (ii<=0 || ii>=YDIM-1 || jj<=0 || jj>=XDIM-1)
			break;

		j = jj;
		i = ii;
	}

	end_pixel.first = i;
	end_pixel.second = j;
}



void Img_MorseDecomp::image_space_trace_forward_Euler (float sx, float sy, int nsteps, pair<int, int> &end_pixel)
{
	float vec[2], vec2[2];

	float x = sx, y = sy;

	//get_vec_at_regular_grid_image_based (x, y, vec[0], vec[1]);

	int i = (int) sy;  //which_row;
	int j = (int) sx;  //which_column;

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	//double pre_ang, cur_ang, ang_diff;
	//pre_ang = atan2(vec[1], vec[0]);

	//closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (step = 0; step < nsteps; step ++)
	{

		if (!Euler_integrate_multi_steps(trace_gp, rand_stepsize, false)) break;

		//j = (int) trace_gp[0];
		//i = (int) trace_gp[1];
		//
		//if (i<=0 || i>=YDIM-1 || j<=0 || j>=XDIM-1)
		//	break;
		int jj = (int) trace_gp[0];
		int ii = (int) trace_gp[1];

		if (ii<=0 || ii>=YDIM-1 || jj<=0 || jj>=XDIM-1)
			break;

		j = jj;
		i = ii;
	}

	end_pixel.first = i;
	end_pixel.second = j;
}

void Img_MorseDecomp::image_space_trace_backward_Euler (float sx, float sy, int nsteps, pair<int, int> &end_pixel)
{
	float vec[2], vec2[2];

	float x = sx, y = sy;

	//get_vec_at_regular_grid_image_based (x, y, vec[0], vec[1]);

	int i = (int) sy;  //which_row;
	int j = (int) sx;  //which_column;

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	//double pre_ang, cur_ang, ang_diff;
	//pre_ang = atan2(vec[1], vec[0]);

	//closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (step = 0; step < nsteps; step ++)
	{

		if (!Euler_integrate_multi_steps(trace_gp, rand_stepsize, true)) break;
		
		//j = (int) trace_gp[0];
		//i = (int) trace_gp[1];

		//if (i<=0 || i>=YDIM-1 || j<=0 || j>=XDIM-1)
		//	break;

		int jj = (int) trace_gp[0];
		int ii = (int) trace_gp[1];

		if (ii<=0 || ii>=YDIM-1 || jj<=0 || jj>=XDIM-1)
			break;

		j = jj;
		i = ii;
	}

	end_pixel.first = i;
	end_pixel.second = j;
}
	
void Img_MorseDecomp::image_space_trace_forward_RK4 (float sx, float sy, int nsteps, pair<int, int> &end_pixel)
{
	float vec[2], vec2[2];

	float x = sx, y = sy;

	//get_vec_at_regular_grid_image_based (x, y, vec[0], vec[1]);

	int i = (int) sy;  //which_row;
	int j = (int) sx;  //which_column;

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	//double pre_ang, cur_ang, ang_diff;
	//pre_ang = atan2(vec[1], vec[0]);

	//closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (step = 0; step < nsteps; step ++)
	{

		if (!RK4_integrate_multi_steps(trace_gp, rand_stepsize, false)) break;

		//j = (int) trace_gp[0];
		//i = (int) trace_gp[1];
		//
		//if (i<=0 || i>=YDIM-1 || j<=0 || j>=XDIM-1)
		//	break;

		int jj = (int) trace_gp[0];
		int ii = (int) trace_gp[1];

		if (ii<=0 || ii>=YDIM-1 || jj<=0 || jj>=XDIM-1)
			break;

		j = jj;
		i = ii;
	}

	end_pixel.first = i;
	end_pixel.second = j;
}

void Img_MorseDecomp::image_space_trace_backward_RK4 (float sx, float sy, int nsteps, pair<int, int> &end_pixel)
{
	float vec[2], vec2[2];

	float x = sx, y = sy;

	//get_vec_at_regular_grid_image_based (x, y, vec[0], vec[1]);

	int i = (int) sy;  //which_row;
	int j = (int) sx;  //which_column;

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	//double pre_ang, cur_ang, ang_diff;
	//pre_ang = atan2(vec[1], vec[0]);

	//closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (step = 0; step < nsteps; step ++)
	{

		if (!RK4_integrate_multi_steps(trace_gp, rand_stepsize, true)) break;
		
		//j = (int) trace_gp[0];
		//i = (int) trace_gp[1];

		//if (i<=0 || i>=YDIM-1 || j<=0 || j>=XDIM-1)
		//	break;

		int jj = (int) trace_gp[0];
		int ii = (int) trace_gp[1];

		if (ii<=0 || ii>=YDIM-1 || jj<=0 || jj>=XDIM-1)
			break;

		j = jj;
		i = ii;
	}

	end_pixel.first = i;
	end_pixel.second = j;
}



bool Img_MorseDecomp::get_nextpt_RK23_quad(double first[2], double second[2], double offset[2], bool type)
{

	double t_vec[2] = {0.};
	if (!get_vec_at_regular_grid_image_based2(first, t_vec, type)) return false;
	icVector2 vec (t_vec);;
	//vec.entry[0] = t_vec[0];
	//vec.entry[1] = t_vec[1];
	if(length(vec) < 5.e-6) return false;

	double eps_did, eps = 1.e-9;
	int i;
	double hstep_loc, hnext;

	/*
	   normalize the predict_stepsize 08/04/2010
	*/
	//if (predict_stepsize < 0.1)
	//	predict_stepsize = 0.1;

	//if (predict_stepsize > 1.)
	//	predict_stepsize = 1.;

	hstep_loc = predict_stepsize;
	if(hstep_loc > 1.)
	{
		hstep_loc = 1.;
		predict_stepsize = hstep_loc;
	}

	if (hstep_loc < 0.05)
	{
		hstep_loc = 0.05;
		predict_stepsize = hstep_loc;
	}

	for(i=0; i<10; i++)
	{
		hstep_loc = predict_stepsize;
		if (!RK23_2d(first, second, offset, hstep_loc, hnext, eps, eps_did, type)) return false;
		predict_stepsize = hnext;
		if(eps_did<eps)
			break;
	}
	return true;
}



/*
Use RK45 to do tensor line computation
*/
bool Img_MorseDecomp::get_nextpt_RK45_quad(double first[2], double second[2], double offset[2], bool type)
{

	double t_vec[2] = {0.};
	if (!get_vec_at_regular_grid_image_based2(first, t_vec, type)) return false;
	icVector2 vec (t_vec);;
	if(length(vec) < 5.e-6) return false;

	double eps_did, eps = 1.e-9;
	int i;
	double hstep_loc, hnext;

	/*
	   normalize the predict_stepsize 08/04/2010
	*/

	hstep_loc = predict_stepsize;
	if(hstep_loc > 1.)
	{
		hstep_loc = 1.;
		predict_stepsize = hstep_loc;
	}

	if (hstep_loc < 0.05)
	{
		hstep_loc = 0.05;
		predict_stepsize = hstep_loc;
	}

	for(i=0; i<5; i++)
	{
		hstep_loc = predict_stepsize;
		if (!RK45_2d(first, second, offset, hstep_loc, hnext, eps, eps_did, type)) return false;
		predict_stepsize = hnext;
		if(eps_did<eps)
			break;
	}
	return true;
}


/*we try to implement the RK23 as numerical recipe*/
bool Img_MorseDecomp::RK23_2d(double pre_p[2], double next_p[2], double offset[2], double &hstep_loc, double &hnext,
			 double eps, double &eps_did, bool backward /*void deriv(double cur_p[2], double vec[2])*/
			 )
{
	double dp0[2], dp1[2];
	double temp[2] = {0.};
	double t_vec[2];
	icVector2 vec1, vec2;
	
	/*compute dp0*/
	if (!get_vec_at_regular_grid_image_based2(pre_p, t_vec, backward)) return false;
	// we need to normalize the vector 
	vec1.set(t_vec);
	normalize(vec1);
	dp0[0] = hstep_loc*/*t_vec*/vec1.entry[0];
	dp0[1] = hstep_loc*/*t_vec*/vec1.entry[1];

	/*compute dp1*/
	temp[0]=pre_p[0]+dp0[0];
	temp[1]=pre_p[1]+dp0[1];
	if (!get_vec_at_regular_grid_image_based2(temp, t_vec, backward)) return false;
	vec2.set(t_vec);
	normalize(vec2);
	dp1[0] = hstep_loc*/*t_vec*/vec2.entry[0];
	dp1[1] = hstep_loc*/*t_vec*/vec2.entry[1];
	

	/*compute the next position using dp0, dp1, dp2 and dp3*/
	offset[0] = dp0[0]/2+dp1[0]/2;
	offset[1] = dp0[1]/2+dp1[1]/2;
	next_p[0]=pre_p[0]+offset[0]/*dp0[0]/2+dp1[0]/2*/;
	next_p[1]=pre_p[1]+offset[1]/*dp0[1]/2+dp1[1]/2*/;

	/*evaluate the error*/
	if (!get_vec_at_regular_grid_image_based2(next_p, t_vec, backward)) return false;
	vec1.set(t_vec);
	normalize(vec1);
	
	icVector2 ep;
	ep.entry[0]=(dp1[0]-hstep_loc*/*t_vec*/vec1.entry[0])/3.;
	ep.entry[1]=(dp1[1]-hstep_loc*/*t_vec*/vec1.entry[1])/3.;

	double error = length(ep);

	/*adjust the step size accordingly*/
	hnext = hstep_loc;
	if(error<eps/10.) hnext = 2.*hstep_loc;
	if(error>eps*5) hnext = hstep_loc/2;

	eps_did = error;

	return true;
}



/*we try to implement the RK23 as numerical recipe*/
bool Img_MorseDecomp::RK45_2d(double pre_p[2], double next_p[2], double offset[2], double &hstep_loc, double &hnext,
			 double eps, double &eps_did, bool backward //, bool deriv(double xy[2], int xdim, int ydim, double vxy[2], bool backward) /*void deriv(double cur_p[2], double vec[2])*/
			 )
{
	double dp0[2], dp1[2], dp2[2], dp3[2];
	double temp[2] = {0.};
	double t_vec[2];
	icVector2 vec1, vec2;
	
	/*compute dp0*/
	if (!get_vec_at_regular_grid_image_based2(pre_p, t_vec, backward)) return false;
	// we need to normalize the vector 
	vec1.set(t_vec);
	normalize(vec1);
	dp0[0] = hstep_loc*/*t_vec*/vec1.entry[0];
	dp0[1] = hstep_loc*/*t_vec*/vec1.entry[1];

	/*compute dp1*/
	temp[0]=pre_p[0]+dp0[0]/2;
	temp[1]=pre_p[1]+dp0[1]/2;
	if (!get_vec_at_regular_grid_image_based2(temp, t_vec, backward)) return false;
	vec2.set(t_vec);
	normalize(vec2);
	dp1[0] = hstep_loc*/*t_vec*/vec2.entry[0];
	dp1[1] = hstep_loc*/*t_vec*/vec2.entry[1];


	/*compute dp2*/
	temp[0]=pre_p[0]+dp1[0]/2;
	temp[1]=pre_p[1]+dp1[1]/2;
	if (!get_vec_at_regular_grid_image_based2(temp, t_vec, backward)) return false;
	vec2.set(t_vec);
	normalize(vec2);
	dp2[0] = hstep_loc*vec2.entry[0];
	dp2[1] = hstep_loc*vec2.entry[1];

	/*compute dp3*/
	temp[0]=pre_p[0]+dp2[0];
	temp[1]=pre_p[1]+dp2[1];
	if (!get_vec_at_regular_grid_image_based2(temp, t_vec, backward)) return false;
	vec2.set(t_vec);
	normalize(vec2);
	dp3[0] = hstep_loc*vec2.entry[0];
	dp3[1] = hstep_loc*vec2.entry[1];

	/*compute the next position using dp0, dp1, dp2 and dp3*/
	next_p[0]=pre_p[0]+dp0[0]/6+dp1[0]/3+dp2[0]/3+dp3[0]/6;
	next_p[1]=pre_p[1]+dp0[1]/6+dp1[1]/3+dp2[1]/3+dp3[1]/6;


	/*compute the next position using dp0, dp1, dp2 and dp3*/
	offset[0] = (dp0[0]+2*dp1[0]+2*dp2[0]+dp3[0])/6.;
	offset[1] = (dp0[1]+2*dp1[1]+2*dp2[1]+dp3[1])/6.;
	next_p[0]=pre_p[0]+offset[0]/*dp0[0]/2+dp1[0]/2*/;
	next_p[1]=pre_p[1]+offset[1]/*dp0[1]/2+dp1[1]/2*/;

	/*evaluate the error*/
	if (!get_vec_at_regular_grid_image_based2(next_p, t_vec, backward)) return false;
	vec1.set(t_vec);
	normalize(vec1);
	
	icVector2 ep;
	ep.entry[0]=(dp3[0]-hstep_loc*vec1.entry[0])/6;
	ep.entry[1]=(dp3[1]-hstep_loc*vec1.entry[1])/6;


	double error = length(ep);

	/*adjust the step size accordingly*/
	hnext = hstep_loc;
	if(error<eps/10.) hnext = 2*hstep_loc;
	if(error>eps*5) hnext = hstep_loc/2;

	eps_did = error;

	return true;
}


//void Img_MorseDecomp::set_vec_img(unsigned char ***in_vec_img, float xmin, float xmax, float ymin, float ymax)
void Img_MorseDecomp::set_vec_img(unsigned char *in_vec_img, float xmin, float xmax, float ymin, float ymax)
{
	int i, j;

	vec_img = new MyColor * [YDIM];

//#pragma omp parallel for  
	for (i=0; i<YDIM; i++)
	{
		vec_img[i] = new MyColor[XDIM];
	}

#pragma omp parallel for 
	for (i=0; i<YDIM; i++)
	{
		for (j=0; j<XDIM; j++)
		{
			//vec_img[i][j].val[0] = in_vec_img[i][j][0];
			//vec_img[i][j].val[1] = in_vec_img[i][j][1];
			//vec_img[i][j].val[2] = in_vec_img[i][j][2];
			vec_img[i][j].val[0] = in_vec_img[3*(i*XDIM+j)];
			vec_img[i][j].val[1] = in_vec_img[3*(i*XDIM+j)+1];
			vec_img[i][j].val[2] = in_vec_img[3*(i*XDIM+j)+2];
		}
	}

	min_vx = xmin;
	max_vx = xmax;
	min_vy = ymin;
	max_vy = ymax;
}


//void Img_MorseDecomp::set_fixedPt_mask(bool **in_fixedPt)
void Img_MorseDecomp::set_fixedPt_mask(bool *in_fixedPt)
{
	if (fixedPt_pixels == NULL)
	{
		fixedPt_pixels = new bool *[YDIM];

		int i;
////#pragma omp parallel for 
		for (i=0; i<YDIM; i++)
			fixedPt_pixels[i] = new bool[XDIM];
	}

	int i, j;

#pragma omp parallel for  
	for (i=0; i<YDIM; i++)
	{
		for (j=0; j<XDIM; j++)
		{
			//fixedPt_pixels[i][j] = in_fixedPt[i][j];
			fixedPt_pixels[i][j] = in_fixedPt[i*XDIM+j];
		}
	}
}


void Img_MorseDecomp::build_one_edge(int node1, int node2)
{
	//#pragma omp critical
	dg->add_to_edgelist(node1, node2, dg->elist->nedges);
	////add the edge to the nodes
	dg->add_edge_to_node(node1, dg->elist->nedges-1);
	dg->add_edge_to_node(node2, dg->elist->nedges-1);
}


void Img_MorseDecomp::build_edges()
{
	int i, j;

//#pragma omp parallel for 
	for (i=0; i<XDIM*YDIM; i++)
	{
		for (j=0; j<end_pixels[i].size(); j++)
		{
			build_one_edge(i, end_pixels[i][j]);
		}
		end_pixels[i].clear();
	}
		
	end_pixels.clear();

}


//void Img_MorseDecomp::init_all_edges()  // This is used for triangular mesh only
//{
//}


void Img_MorseDecomp::cal_Img_MorseDecomp(int nsteps)
{
	init_graph();

	printf("building directed graph ...\n");
	build_DG_tau (nsteps);

	printf("finding SCC...\n");
	dg->find_SCCS();

	printf ("building SCC list...\n");
	build_SCCElemList();

	printf("marking valid SCCs...\n");
	mark_all_valid_SCCS();
}


/*
build the SCC list according to the obtained SCC labels of the nodes
*/

void Img_MorseDecomp::build_SCCElemList()
{
	int i;
	////first, we need to search the sub-trees to figure how many scc in the graph and
	////what are they
	//for(i = 0; i < dg->nlist->ndirnodes; i++)
	//{
	//	dg->nlist->dirnodes[i]->visited = 0;
	//}

	/* To get the elements in each component */
	int cur_sccomp = 0;
	int num_morsesets = 0;

	//initialization

	scclist.clear();

	// allocate space for each SCC
	for (i=0; i < dg->num_sccomps; i++)
	{
		SCComponent scc;
		scclist.push_back(scc);
	}


	//record the scc components

//#pragma omp parallel for nowait  (share variables!)
	//#pragma omp parallel for
	for(i = 0; i < dg->nlist->ndirnodes /*the number of total pixels*/; i++)
	{
		int row = i/XDIM;
		int col = i%XDIM;

		pair<int, int> pixel;
		pixel.first = row;
		pixel.second = col;

		#pragma omp critical
		scclist[dg->nlist->dirnodes[i]->sscomp_index].pixels.push_back(pixel);

		/* count the number of fixed points in the SCC 02/21/07 */
		if(fixedPt_pixels[row][col])
		{
			pixel.first = row;
			pixel.second = col;

			#pragma omp critical
			scclist[dg->nlist->dirnodes[i]->sscomp_index].singular_pixels.push_back(pixel);
		}
	}


	// we can leave some room for the later local refinement process

}



void Img_MorseDecomp::mark_all_valid_SCCS()
{

//#pragma omp parallel for  
	for(int i = 0; i < scclist.size(); i++)
		scclist[i].valid = true;

	valid_MS = 0;
	

//#pragma omp parallel for  
	for(int i = 0; i < scclist.size(); i++)
	{
		scclist[i].classification = 0;  // trivial Morse sets

		if(scclist[i].pixels.size() < 2)
		{
			if(scclist[i].singular_pixels.size() == 0)
			{
				scclist[i].valid = false;
				continue;
			}
			continue;
	
		}

		// get the type of the valid Morse sets
		scclist[i].classification = get_Type_MorseSet(i);

		valid_MS++;

		//else if(!is_valid_SCC(i))  // this is test the shape of the Morse set regions (comment out now)
		//{
		//	scclist[i].valid = false;
		//	continue;
		//}

		//// The following test the connectivity of the Morse sets, 
		//// We comment it out at this moment and see what we can get
		
		//int *temp_tris = new int [scclist->scccomponents[i]->nnodes];
		//for (int j=0; j<scclist->scccomponents[i]->nnodes; j++)
		//{
		//	temp_tris[j]=scclist->scccomponents[i]->nodes[j];
		//}

		///*
		//    Added by Guoning 07/01/2010
		//*/

		//if (RemoveDisconnMSOn&&!is_connected_reg(temp_tris, scclist->scccomponents[i]->nnodes))
		//{
		//	scclist->scccomponents[i]->valid = false;
		//	continue;
		//}

		//delete [] temp_tris;

		/*
		   How to remove the boundary Morse sets?? 03/01/2010
		*/
		//else
			//scclist->scccomponents[i]->valid = true;   //this is a valid SCC

	}
}



int Img_MorseDecomp::get_Type_MorseSet(int scc_id)
{
	SCComponent &sccomp = scclist[scc_id];
	int i, j;
	int noutgoing = 0;
	int nincoming = 0;

	for (i=0; i<sccomp.pixels.size(); i++)
	{
		//search the edges associated with each node in this SCC component

		int node_id = sccomp.pixels[i].first*XDIM + sccomp.pixels[i].second;
		DirGraph_Node *cur_n = dg->nlist->dirnodes[node_id];

		for (j=0; j<cur_n->nedges; j++)
		{
			//count the # of outgoing and incoming edges

			Graph_Edge *cur_e = dg->elist->edges[cur_n->edges[j]];

			DirGraph_Node *n1 = dg->nlist->dirnodes[cur_e->node_index1];
			DirGraph_Node *n2 = dg->nlist->dirnodes[cur_e->node_index2];

			if (n1->sscomp_index == n2->sscomp_index) // this is an inner edge
				continue; 

			if (n1 == cur_n) //this is an outgoing edge
				noutgoing ++;

			else 
				nincoming ++;

		}
	}

	if (nincoming == 0)
		return 1; // this is a source
	else if (noutgoing == 0)
		return 2; // this is a sink
	else
		return 3; // it is a saddle
}


void Img_MorseDecomp::init_graph()
{
	int i;
	
	//FILE *fp;
	//fp = fopen("detect_porbit_cooling.txt", "a");
	//fprintf(fp, "Start initializing Directed Graph...\n");
	//fclose(fp);

	if(dg != NULL)
		delete dg;
	
	//fp = fopen("detect_porbit_cooling.txt", "a");
	//fprintf(fp, "Start initializing dg...\n");
	//fprintf(fp, "nodes %d, edges %d...\n", object->tlist.ntris, object->elist.nedges);
	//fclose(fp);

	/*this initialization is suitable for Morse decomposition without using \tau */
	//dg = new DirGraph(XDIM*YDIM, XDIM*YDIM*nsamples*2);
	dg = new DirGraph(XDIM*YDIM, XDIM*YDIM*min(nsamples, 10)*2);

	if(dg->elist == NULL || dg->nlist->dirnodes == NULL || dg->elist->edges == NULL)
	{
		char rout[256], var[256];
		sprintf(rout, "%s", "MorseDecomp::init_graph");
		if(dg->nlist == NULL)
			sprintf(var, "%s", "dg->nlist");
		else
			sprintf(var, "%s", "dg->elist");

		exit(-1);
	}
	

	/*allocate the real space for the node list and edge list*/
//#pragma omp parallel for  
	for(i = 0; i < dg->nlist->curMaxNumENodes; i++)
	{
		dg->nlist->dirnodes[i] = new DirGraph_Node();
		//dg->nlist->dirnodes[i]->edges = NULL;
		//dg->nlist->dirnodes[i]->nedges = 0;
	}
	
	//fp = fopen("detect_porbit_cooling.txt", "a");
	//fprintf(fp, "Start initializing elist...\n");
	//fclose(fp);

//#pragma omp parallel for  
	for(i = 0; i < dg->elist->curMaxNumGedges; i++)
	{
		dg->elist->edges[i] = new Graph_Edge();
	}

	/*initialize the SCC list as well here*/

	//if(scclist != NULL)
	//	delete scclist;

	//scclist = new SCCList();  /*initial size is 1000 now*/

}

extern void	HsvRgb( float[3], float [3] );

//
void 
get_random_Color(int val, float rgb[3])
{
	int c = val%10;

	float hsv[3];
	hsv[0] = 36 * c;
	hsv[1] = 1.-(float)rand()/(2.0*RAND_MAX);
	hsv[2] = 1.;

	HsvRgb(hsv, rgb);
}

void Img_MorseDecomp::mark_valid_SCC_pixels (unsigned char *tex, unsigned char *out) // this is a test code for visualizing the obtained Morse sets
{
	int i, j;


	if (dg == NULL) return;

	if (scclist.empty()) return;

//#pragma omp parallel for 
	for (i=0; i<YDIM; i++)
	{
		for (j=0; j<XDIM; j++)
		{
			out[3*(i*XDIM+j)] = tex[4*(i*XDIM+j)]; 
			out[3*(i*XDIM+j)+1] = tex[4*(i*XDIM+j)+1]; 
			out[3*(i*XDIM+j)+2] = tex[4*(i*XDIM+j)+2]; 
		}
	}
	

//#pragma omp parallel for 
// OK, let us assign different colors
	float rgb[3];
	for (i=0; i<scclist.size(); i++)
	{
		if (!scclist[i].valid) continue;

	    //get_random_Color(i, rgb);

		if (scclist[i].classification == 1)
		{
			rgb[0] = 0; rgb[1] = 1; rgb[2] = 0;
		}
		else if (scclist[i].classification == 2)
		{
			rgb[0] = 1; rgb[1] = rgb[2] = 0.;
		}
		else if (scclist[i].classification == 3)
		{
			rgb[0] = rgb[1] = 0.; rgb[2] = 1;
		}
		else
		{
			rgb[0] = rgb[1] = 0.7; rgb[2] = 0.2;
		}

		for (j=0; j<scclist[i].pixels.size(); j++)
		{
			int row = scclist[i].pixels[j].first;
			int col = scclist[i].pixels[j].second;

			//tex[4*(row*XDIM+col)] = 0.6 * 255 + 0.4*tex[4*(row*XDIM+col)];
			//tex[4*(row*XDIM+col)+2] = 0.6 * 200 + 0.4*tex[4*(row*XDIM+col)+2];
			tex[4*(row*XDIM+col)] = 0.6 * rgb[0]*255 + 0.4*tex[4*(row*XDIM+col)];
			tex[4*(row*XDIM+col)+1] = 0.6 * rgb[1]*255 + 0.4*tex[4*(row*XDIM+col)+1];
			tex[4*(row*XDIM+col)+2] = 0.6 * rgb[2]*255 + 0.4*tex[4*(row*XDIM+col)+2];

			out[3*(row*XDIM+col)] = tex[4*(row*XDIM+col)]; 
			out[3*(row*XDIM+col)+1] = tex[4*(row*XDIM+col)+1]; 
			out[3*(row*XDIM+col)+2] = tex[4*(row*XDIM+col)+2]; 

		}
	}
}