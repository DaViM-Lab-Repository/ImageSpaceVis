#include <stdio.h>
	// yes, I know stdio.h is not good C++, but I like the *printf()
#include <stdlib.h>
#include <ctype.h>

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef WIN32
#include <windows.h>
#pragma warning(disable:4996)
#endif


#define SHOWSECONDTHIRDWINS

// You need to adjust the location of these header files according to your configuration

#include <omp.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include "glut.h"
#include "glui.h"

#include "Skeleton.h"

#include "Img_MorseDecomp.h"

//#include "glui.h"

//
//
//	This is a sample OpenGL / GLUT / GLUI program
//
//	The objective is to draw a 3d object and change the color of the axes
//		with radio buttons
//
//	The left mouse button allows rotation
//	The middle mouse button allows scaling
//	The glui window allows:
//		1. The 3d object to be transformed
//		2. The projection to be changed
//		3. The color of the axes to be changed
//		4. The axes to be turned on and off
//		5. The transformations to be reset
//		6. The program to quit
//
//	Author: Joe Graphics
//


//
// constants:
//
// NOTE: There are a bunch of good reasons to use const variables instead
// of #define's.  However, Visual C++ does not allow a const variable
// to be used as an array size or as the case in a switch() statement.  So in
// the following, all constants are const variables except those which need to
// be array sizes or cases in switch() statements.  Those are #defines.
//
//


// This source code has been modified by Guoning Chen since its release


// title of these windows:

const char *WINDOWTITLE = { "Flow Rotation Operator Explorer" };
const char *GLUITITLE   = { "User Interface Window" };


// what the glui package defines as true and false:

const int GLUITRUE  = { true  };
const int GLUIFALSE = { false };


// the escape key:

#define ESCAPE		0x1b


// initial window size:

const int INIT_WINDOW_SIZE = { /*512*/1024/*512*/ };


// size of the box:

const float BOXSIZE = { 2.f };



// multiplication factors for input interaction:
//  (these are known from previous experience)

const float ANGFACT = { 1. };
const float SCLFACT = { 0.005f };


// able to use the left mouse for either rotation or scaling,
// in case have only a 2-button mouse:

enum LeftButton
{
	ROTATE,
	SCALE
};


// minimum allowable scale factor:

const float MINSCALE = { 0.05f };


// active mouse buttons (or them together):

const int LEFT   = { 4 };
const int MIDDLE = { 2 };
const int RIGHT  = { 1 };


// which projection:

enum Projections
{
	ORTHO,
	PERSP
};


// which button:

enum ButtonVals
{
	RESET,
	QUIT,
	SAVEIMAGE
};


// window background color (rgba):

const float BACKCOLOR[] = { 0., 0., 0., 0. };


// line width for the axes:

const GLfloat AXES_WIDTH   = { 3. };


// the color numbers:
// this order must match the radio button order

enum Colors
{
	RED,
	YELLOW,
	GREEN,
	CYAN,
	BLUE,
	MAGENTA,
	WHITE,
	BLACK
};


// the object numbers:
// 
enum MODELS
{
	BUNNY,
	FELINE,
	DRAGON,
	HAPPY,
	SPHERE,
	TORUS,
};


enum Slidertype{
	TEMP,
	XRANGE,
	YRANGE,
	ZRANGE,
	SRANGE
};

// the color definitions:
// this order must match the radio button order

const GLfloat Colors[8][3] = 
{
	{ 1., 0., 0. },		// red
	{ 1., 1., 0. },		// yellow
	{ 0., 1., 0. },		// green
	{ 0., 1., 1. },		// cyan
	{ 0., 0., 1. },		// blue
	{ 1., 0., 1. },		// magenta
	{ 1., 1., 1. },		// white
	{ 0., 0., 0. },		// black
};


// fog parameters:

const GLfloat FOGCOLOR[4] = { .0, .0, .0, 1. };
const GLenum  FOGMODE     = { GL_LINEAR };
const GLfloat FOGDENSITY  = { 0.30f };
const GLfloat FOGSTART    = { 1.5 };
const GLfloat FOGEND      = { 4. };



//
// non-constant global variables:
//

int	ActiveButton;		// current button that is down
GLuint	AxesList;		// list to hold the axes
int	AxesOn;			// != 0 means to draw the axes
int	DebugOn;			// != 0 means to print debugging info
int	DepthCueOn;		// != 0 means to use intensity depth cueing
GLUI *	Glui;			// instance of glui window
int	GluiWindow;		// the glut id for the glui window
int	LeftButton;		// either ROTATE or SCALE
GLuint	BoxList;		// object display list
int	MainWindow;		// window id for main graphics window
GLfloat	RotMatrix[4][4];	// set by glui rotation widget
float	Scale, Scale2;		// scaling factors
int	WhichColor;		// index into Colors[]
int	WhichProjection;	// ORTHO or PERSP
int	Xmouse, Ymouse;		// mouse values
float	Xrot, Yrot;		// rotation angles in degrees
float	TransXYZ[3];		// set by glui translation widgets

int ArrowsOn = 0;
int RotSumDiffOn = 0;


int SecondWindow;  // the handler for the second windowm
int ThirdWindow;

int FreezeOn = 0;
int ShowRotSignOn = 0;
int ShowLocalRotSignOn = 0;  // we need additional information to enable this visualization
int ShowGroupColorsOn = 0;
int WhichWindow = 0;

int MainWindowImgNum = 0;
int SecondWindowImgNum = 0;
int ThirdWindowImgNum = 0;

//
// function prototypes:
//

void	Animate( void );
void	Buttons( int );
void	Display( void );
void	DoRasterString( float, float, float, unsigned char * );
void	DoStrokeString( float, float, float, float, char * );
float	ElapsedSeconds( void );
void	InitGlui( void );
void	InitGraphics( void );
void	InitLists( void );
void	Keyboard( unsigned char, int, int );
void	MouseButton( int, int, int, int );
void	MouseMotion( int, int );
void	Reset( void );
void	Resize( int, int );
void	Visibility( int );

void	Arrow( float [3], float [3] );
void	Cross( float [3], float [3], float [3] );
float	Dot( float [3], float [3] );
float	Unit( float [3], float [3] );
void	Axes( float );
void	HsvRgb( float[3], float [3] );

void    Display_Model(void);
void    set_view(GLenum mode, Polyhedron *poly);
void    set_scene(GLenum mode, Polyhedron *poly);
void    display_shape(GLenum mode, Polyhedron *this_poly);
void    Choose_Object();


int     display_3D = 0;

void    display_func();
void	Display_3D( void );
void    Sliders( int id );

//// Function for the second windows
void	Keyboard_win2( unsigned char, int, int );
void	MouseButton_win2( int, int, int, int );
void	MouseMotion_win2( int, int );
void    Display_func_win2();
int     win2_size = 512;
void    draw_arrow_head(float head[2], float direct[2]);


void    Display_func_win3();

////////////////////////////////////////////////////////////////////////////////////////////////////\
///  Image-space method

//// For LIC (10/02/2012)
#include <vector>
std::vector<double> weights;   // the list for the weightsm
std::vector<std::pair<int, int>> passed_pixels;
const int img_res = 512/*256*//*512*//*1024*/;
unsigned char noise_tex[img_res][img_res][3];
unsigned char vec_img[img_res][img_res][3];   // render the vector field into an image
float vec_img2[img_res][img_res][3];   // render the vector field into an image
unsigned char LIC_tex[img_res][img_res][4], out_LIC[img_res][img_res][3];
unsigned char LIC_tex2[img_res][img_res][4];
unsigned char LIC_tex_comb[img_res][img_res][3];
unsigned int forward_counter[img_res][img_res];
unsigned int backward_counter[img_res][img_res];

unsigned char err_img[img_res][img_res][3];


/*GLfloat*/unsigned char forward_counter_color[img_res][img_res][3];
/*GLfloat*/unsigned char backward_counter_color[img_res][img_res][3];
/*GLfloat*/unsigned char comb_counter_color[img_res][img_res][3];
float total_rotation[img_res][img_res];
unsigned char rot_sum_color[img_res][img_res][3];
float total_rotation_diff[img_res][img_res][2];
float total_rotation_diff_mag[img_res][img_res];
float total_rotation_diff_Hessian[img_res][img_res][4];
bool visited_pixels[img_res][img_res];
bool fixedPt_pixels[img_res][img_res];

const int secondwin_size = 512;
const int thirdwin_size = 512;
unsigned char secondwin_img[secondwin_size][secondwin_size][3];
unsigned char thirdwin_img[thirdwin_size][thirdwin_size][3];


float vec_mag[img_res][img_res];

float total_rotation_diff_Hessian_directional_mag[img_res][img_res];

unsigned char total_rotation_diff_color[img_res][img_res][3];

int sel_pos_i = 0, sel_pos_j = 0;

//int L=50000/*3000*/;  // 2000 steps is sufficient for an integrator with constant step size but not enough for an adaptive integrator
int L= 2000; //3000;

double L_percentage = 1.;

// Here we generate the vector image instead
float max_vx, min_vx, max_vy, min_vy;
float vec_2d[img_res][img_res][2];

std::vector<std::pair<int, int>> one_streamline, forward_streamline, backward_streamline;
std::vector<float> current_total_rot, forward_rot, backward_rot;

//// Parameter for double gyre
double Amp = 0.1;
double EPSILON = 0.25;
//double OMEGA = 0.628/*0.25*/;
double OMEGA = 0.5/*0.25*/;

double T_window = 5/*15.*/;


//  For record the end position of each particle starting at the center of a pixel
float start_pos[img_res][img_res][2];
float end_positions_forward[img_res][img_res][2];
float spatial_gradient_forward[img_res][img_res][4];   //// the 2x2 spatial gradient matrix
float FTLE_forward[img_res][img_res];
unsigned char FTLE_forward_color[img_res][img_res][3];




void    normalize_Field();

void    update_kernel(int size);
void    gen_test_vecfld ();
void    gen_noise_tex ();
void    replace_noise_tex (unsigned char LIC_tex[img_res][img_res][4]);
void    est_passed_pixels_from (int i, int j, int L);
void    comp_streamline_forward (int i, int j, int L, std::vector<std::pair<int, int>> &passed_pixels);
void    comp_streamline_backward (int i, int j, int L, std::vector<std::pair<int, int>> &passed_pixels);
void    compose_color (unsigned char color_val[3]);
void    compose_color (unsigned char color_val[3], std::vector<std::pair<int, int>> passed_pixels);
void    comp_LIC (unsigned char LIC_tex[img_res][img_res][4]);
void    render_vec_img( Polyhedron *this_poly);
void    jitter_vectorFld (Polyhedron *this_poly, double prob, double perc);
void    interpolate_error_bilinear(Polyhedron *this_poly);

void    render_LIC_img();
void    rotate_vec( Polyhedron *this_poly, double degree);
void    combine_LIC();
void    reflect_vec( Polyhedron *this_poly);

void    comp_LIC_counter(unsigned int forward_[img_res][img_res], unsigned int backward_[img_res][img_res]);
void    init_counters();
void    get_color_map_for_counters();
void    get_color_map_for_counters_forward();

void    comp_LIC_total_rotation(float total_rotation[img_res][img_res]);
void    comp_LIC_total_rotation_opt(float total_rotation[img_res][img_res]);
void    get_color_map_for_rot_sum();
void    init_rotation_sum();

void    trace_streamline_from_image_based(int i, int j, int L);
void    merge_current_forward_back_streamlines();
void    merge_current_forward_back_streamlines(std::vector<std::pair<int, int>> &one_streamline, 
				std::vector<std::pair<int, int>> forward_streamline, 
				std::vector<std::pair<int, int>> backward_streamline,
				std::vector<float>  &current_total_rot, 
				std::vector<float>	forward_rot, 
				std::vector<float>	backward_rot
				);

void    comp_total_rotation_diff(float total_rotation[img_res][img_res], float total_rotation_diff[img_res][img_res][2]);
void    get_color_map_for_rot_sum_diff();

void    comp_total_rotation_diff_Hessian();
void    get_color_map_for_rot_sum_diff_Hessian();

///// Generate a test example (static double gyre)
void    gen_static_double_gyre();
void    cal_vec_double_gyre_static(double x, double y, double &vx, double &vy);

void    cal_SpatialGradient();
void    compute_FTLE_forward();

void    init_particle_list();
void    insert_new_particles();  // insert new particles at the center of each pixel
void    update_all_particles_streaklines (double cur_t, double dt);
void    get_rot_sum(int start_k);

double    sum_rotation_backward_with_streamline (int i, int j, int L, bool &closedloop, bool store, std::vector<std::pair<int, int>> &backward_streamline, std::vector<float> &backward_rot);
double    sum_rotation_forward_with_streamline (int i, int j, int L, bool &closedloop, bool store, std::vector<std::pair<int, int>> &forward_streamline, std::vector<float> &forward_rot);

double    sum_rotation_backward_with_streamline_RK23 (int i, int j, int L, bool &closedloop, bool store, std::vector<std::pair<int, int>> &backward_streamline, std::vector<float> &backward_rot);
double    sum_rotation_forward_with_streamline_RK23 (int i, int j, int L, bool &closedloop, bool stor, std::vector<std::pair<int, int>> &forward_streamline, std::vector<float> &forward_rote);


void 
insert_new_particles_for_streakline_test();
void 
init_particle_list_test();

//// Test the accuracy of the advection
const int PDIM = 32;
float particle_pos[PDIM][PDIM][2];
double cur_ptime = 0.;
double max_time = 1.;
double d_time = max_time/60.;

//// Now, for streakline computation
std::vector<std::pair<float, float>> particles;   // the list of particles
std::vector<bool> particle_active;
std::vector<std::vector<int>> pixels_to_particles; // which particles in the particle list are released at certain pixels


void    init_par_pos();
void    advect_par(double &cur_time, double dt);

////////////////////////////////////////////////////////////////////////////////////////////////
//// Later we need to generalize it to the object-space and consider parallel implementation!


///// For time-dependent vector field (time-varying double gyre)
void    cal_vec_double_gyre(double x, double y, double t, double &vx, double &vy);
void    comp_LIC_total_rotation_timeForward(float total_rotation[img_res][img_res], double start_time);
void    gen_double_gyre_at_time(double time);

/////
void    cal_vec_ibfv_ex(double x, double y, double t, double &vx, double &vy);
void    gen_ibfv_at_time(double time);

void    cal_gyre_saddle(double x, double y, double t, double &vx, double &vy);
void    gen_gyre_saddle_at_time(double time);

void    cal_stuart_vortices(double x, double y, double t, double &vx, double &vy);
void    gen_stuart_vortices_at_time(double time);


void    cal_Oseen_vortices(double x, double y, double t, double &vx, double &vy);
void    gen_Oseen_vortices_at_time(double time);


//// Test 3D steady flow
const int GRIDSIZE = 32/*256*/;
icVector3 data_center;

//// Data domain for ABC flow
//double MIN_COORD = 0, MAX_COORD = 2*PI;
//double MIN_X=0, MIN_Y=0, MIN_Z=0;
//double MAX_X=2*PI, MAX_Y=2*PI, MAX_Z=2*PI;


//// Data domain for Lorenz equation
double MIN_COORD = -50, MAX_COORD = 50;
double MIN_X=-50., MIN_Y=-50, MIN_Z=-10;
double MAX_X=50, MAX_Y=50, MAX_Z=90;



//// Data domain for Rossler equation
//double MIN_COORD = -20, MAX_COORD = 30;
//double MIN_X=-20., MIN_Y=-20, MIN_Z=-20;
//double MAX_X=30, MAX_Y=30, MAX_Z=30;


//// Data domain for general equation
//double MIN_COORD = -2, MAX_COORD = 2;
//double MIN_X=-5., MIN_Y=-5, MIN_Z=-5;
//double MAX_X=5, MAX_Y=5, MAX_Z=5;


double INI_STEP_SIZE = (MAX_X-MIN_X)/(2*GRIDSIZE);

float rot_sum_3d[GRIDSIZE][GRIDSIZE][GRIDSIZE];
float rot_sum_rgb[GRIDSIZE][GRIDSIZE][GRIDSIZE][3];

float rot_sum_3d_diff[GRIDSIZE][GRIDSIZE][GRIDSIZE][3];
float rot_sum_3d_diff_mag[GRIDSIZE][GRIDSIZE][GRIDSIZE];
float rot_sum_diff_rgb[GRIDSIZE][GRIDSIZE][GRIDSIZE][3];
int cell_counters_3d[GRIDSIZE][GRIDSIZE][GRIDSIZE];

// For volume rendering
unsigned char   TextureXY[GRIDSIZE][GRIDSIZE][GRIDSIZE][4]; //XY plane
unsigned char   TextureXZ[GRIDSIZE][GRIDSIZE][GRIDSIZE][4]; //XZ plane
unsigned char   TextureYZ[GRIDSIZE][GRIDSIZE][GRIDSIZE][4]; //YZ plane

/* which way is a surface facing:                    */

const int MINUS = { 0 };
const int PLUS  = { 1 };

float SRange[2];
float MaxAlpha = 0.3;
int    Major;            /* X, Y, or Z                */
int    Xside, Yside, Zside;    /* which side is visible, PLUS or MINUS    */
int Bilinear;

int   Show_ROTSUM_or_DIFF = 0;
int   Show_OnlyLIC = 0;
int   ShowSingularities = 0;
int   ShowSelStreamline = 0;
int   UseRK23 = 1;
int   Disable3DTexture = 0;

////Global Glui element variables
GLUI_HSlider *ScaleSlider, *AlphaSlider;
GLUI_StaticText *TempLabel, *XRangeLabel, *YRangeLabel, *ZRangeLabel, *SRangeLabel;

GLUI_Spinner *Arrow_Scale_Spinner;

float arrow_scale = 0.05;

double  get_rot_ang_3D(double p1[3], double p2[3], double p[3], icVector3 &pre_rotAxis);  // the order of the points matters!!!!!
void    comp_rot_sum_3D_steady(int L);
void    get_rot_sum_3d_color_map();
void    output_rot_sum_result(char *filename);
void    output_rot_sum_diff_result(char *filename);
void    load_rot_sum_result(char *filename);
void    FillXY( void );
void    FillXZ( void );
void    FillYZ( void );
void    UpdateXY( void );
void    UpdateXZ( void );
void    UpdateYZ( void );
void    DetermineVisibility();
void    Drawing();

void    comp_rot_sum_diff_3d();
void    get_color_map_for_rot_sum_diff_3d();
void    get_color_map_for_rot_sum_diff_3d_from_file();


// 1. ABC flow
double AA = 1.;
double BB = sqrt(2./3.);
double CC = 1./sqrt(3.);
void    get_vec_ABC_flow(double x, double y, double z, double &vx, double &vy, double &vz);

// 2. Lorentz equation
double Lorenz_A = 10.;
double Lorenz_B = 28.;
double Lorenz_C = 8./3.;
void    get_vec_Lorenz_eq(double x, double y, double z, double &vx, double &vy, double &vz);

// 3. Rossler equations
double Rossler_A = 0.2;
double Rossler_B = 0.2;
double Rossler_C = 5.7;
//double Rossler_A = 0.1;
//double Rossler_B = 0.1;
//double Rossler_C = 14.;
void    get_vec_Rossler_eq(double x, double y, double z, double &vx, double &vy, double &vz);

void    get_vec_general_3D( double x, double y, double z, double &vxp, double &vyp, double &vzp );


// 4. Chen and Lee system
double Chen_Lee_A = 5;
double Chen_Lee_B = -10;
double Chen_Lee_C = -0.38;
double one_third = 1./3.;
void    get_vec_Chen_Lee_eq(double x, double y, double z, double &vx, double &vy, double &vz);


// 5. Rabin_Fab equation
// RF_A = 0.1, RF_B = 0.98 chaotic, RF_B = 0.14 stable limit cycle
double RF_A = 0.87;  
double RF_B = 1.1;
double RF_C;
void    get_vec_RF_eq(double x, double y, double z, double &vx, double &vy, double &vz);

void    write_ppm(char *filename, unsigned char *img, int dimx, int dimy);
void    write_ppm_flippedY(char *filename, unsigned char *img, int dimx, int dimy);


// For the visualization of 3D streamlines
typedef struct Point3D
{
	float xyz[3];
}Point3D;

typedef struct Streamline3D
{
	int i, j, k;  // index of the starting voxel
	std::vector<Point3D> one_line;
}Streamline3D;

std::vector<Streamline3D> list_streamlines_3d;

void  update_list_streamlines_3d();

void  comp_streamline_3D_steady_at(int i, int j, int k, int L, std::vector<Point3D> &one_line);


double radius_factor = 0.9;
int display_mode = 0; 

int ObjectId = 0;
char object_name[128]= "bunny";

Polyhedron *poly = NULL;


typedef struct MS_Label
{
	unsigned char labels[img_res][img_res]; //-1 -- not in any Morse sets  0 -- source, 1 -- sink, 2 -- saddle
};

int MS_hits[img_res][img_res][2];

std::vector<MS_Label> MS_Labels;

void
store_MS_Label(Img_MorseDecomp &ismd);

void 
count_MS_hits();

void
blend_color_uncertainty(int Ntrials);

//#define MSCOMPARE

int n_I = 0;
int n_II = 0;

unsigned char *tri_labels = NULL;

unsigned char tri_ID_img[img_res][img_res][3];
unsigned char pixel_tri_labels[img_res][img_res];
unsigned char pixel_labels[img_res][img_res];

unsigned char n_I_img[img_res][img_res][3];
unsigned char n_II_img[img_res][img_res][3];
unsigned char combined_n_I_n_II_img[img_res][img_res][3];
unsigned char gray_MS_img[img_res][img_res][3];


void 
render_tri_ID_img(Polyhedron *this_poly);

void 
assign_tri_pixel_labels(Img_MorseDecomp &ismd);

void compute_n_I_and_n_II();

void
get_gray_MS_img(Img_MorseDecomp &ismd);


//
// main program:
//

int
main( int argc, char *argv[] )
{
	// turn on the glut package:
	// (do this before checking argc and argv since it might
	// pull some command line arguments out)

	glutInit( &argc, argv );

	// Load the model and data here
	//FILE *this_file = fopen("../models/bunny.ply", "r");
	//FILE *this_file = fopen("../models/cnoise.ply", "r");
	//FILE *this_file = fopen("../models/dipole.ply", "r");
	//FILE *this_file = fopen("../models/complex_fld1.ply", "r");
	//FILE *this_file = fopen("../models/two_saddles1.ply", "r");
	//FILE *this_file = fopen("../models/four_centers.ply", "r");
	//FILE *this_file = fopen("../models/bnoise.ply", "r");
	//FILE *this_file = fopen("../models/vnoise.ply", "r");
	//FILE *this_file = fopen("../models/many_centers.ply", "r");
	//FILE *this_file = fopen("../models/one_general_fld.ply", "r");  // very good field
	FILE *this_file = fopen("../models/one_repelling_orbit.ply", "r");  // constrained by the accuracy of the streamline tracing
	//FILE *this_file = fopen("../models/one_saddle.ply", "r");
	//FILE *this_file = fopen("../models/one_cw.ply", "r");
	//FILE *this_file = fopen("../models/one_ccw.ply", "r");
	//FILE *this_file = fopen("../models/attractfocus.ply", "r");
	//FILE *this_file = fopen("../models/repelfocus.ply", "r");
	//FILE *this_file = fopen("../models/constant_fld.ply", "r");
	//FILE *this_file = fopen("../models/attachment_elem.ply", "r");
	//FILE *this_file = fopen("../models/diesel_m0830_circle_s2000.ply", "r");
	//FILE *this_file = fopen("../models/diesel_m0630_circle_s2000.ply", "r");
	//FILE *this_file = fopen("../models/diesel_m0730_circle_s2000.ply", "r");
	//FILE *this_file = fopen("../models/ex_field1.ply", "r");
	//FILE *this_file = fopen("../models/abc_slice.ply", "r");
  	//FILE *this_file = fopen("../models/abcslice2.ply", "r");

	//FILE *this_file = fopen("../models/iceland_current_field.ply", "r");
	//FILE *this_file = fopen("../models/fig9ex.ply", "r");
	//FILE *this_file = fopen("../models/multicycles_perver.ply", "r");
	//FILE *this_file = fopen("../models/vis07_exfield_perver.ply", "r");
	//FILE *this_file = fopen("../models/20deg5_0mps_highpass_0001.ply", "r");

	//FILE *this_file = fopen("../models/tile1_00020314.ply", "r");
	//FILE *this_file = fopen("../models/satlantic_00020609.ply", "r");
	//FILE *this_file = fopen("../models/hcci_00137.ply", "r");

	//FILE *this_file = fopen("../models/windfield_houston.ply", "r");
	//FILE *this_file = fopen("../models/windfield_us2.ply", "r");
	//FILE *this_file = fopen("../models/smallscale_swirl.ply", "r");
	//FILE *this_file = fopen("../models/smallscale_swirl2.ply", "r");

	poly = new Polyhedron (this_file);
	fclose(this_file);
	//mat_ident( rotmat );	

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	
	// normalize the domain if it is not the unit square [0,1]x[0,1]
	poly->normalized_domain();
	poly->calc_bounding_sphere();

	poly->calc_face_normals_and_area();
	poly->average_normals();

	normalize_Field();
	capture_singular_tris();

	// setup all the graphics stuff:

	InitGraphics();

#ifdef SHOWSECONDTHIRDWINS
	// initialize the second window
	glutInitWindowPosition (1200, 20);
	glutInitWindowSize (secondwin_size, secondwin_size); 
	SecondWindow = glutCreateWindow ("Cartesian rotation plot");
	glutKeyboardFunc (Keyboard_win2);
	glutDisplayFunc(Display_func_win2); 
	glutMotionFunc ( MouseMotion_win2);
	glutMouseFunc (MouseButton_win2);
	
	// initialize the third window
	glutInitWindowPosition (1200, 590);
	glutInitWindowSize (thirdwin_size, thirdwin_size); 
	ThirdWindow = glutCreateWindow ("Polar Coordinate View");
	//glutKeyboardFunc (Keyboard_win3);
	glutDisplayFunc(Display_func_win3); 
	//glutMotionFunc ( MouseMotion_win3);
	//glutMouseFunc (MouseButton_win3);
#endif


	// create the display structures that will not change:

	InitLists();


	// init all the global variables used by Display():
	// this will also post a redisplay
	// it is important to call this before InitGlui()
	// so that the variables that glui will control are correct
	// when each glui widget is created

	Reset();


	glutSetWindow(MainWindow);


	// Test LIC here
	gen_noise_tex ();

	////gen_test_vecfld();
	//normalize_Field();

	//// We can try two ways to jitter the original vector fields:
	// 1. jitter the vectors defined at the triangular mesh  (may work better!)
	// 2. jitter the vector field that is already projected in the image plane 
	
	//jitter_vectorFld (poly, 0.8, 0.2);

	render_vec_img(poly);

	// Here we can make a test to evaluate the amount of error introduced using the image-space representation
	interpolate_error_bilinear(poly);

	omp_set_dynamic (0);
	omp_set_num_threads (100);
	
	clock_t start, finish; //Used to show the time assumed
	start = clock(); //Get the first time

	comp_LIC(LIC_tex);

	// replace the white noise with the LIC image and do another iteration of LIC


	// SHOW 3D field
	///* DEBUG RK23

	for (int i=0; i<3; i++)
	{
		replace_noise_tex (LIC_tex);
		comp_LIC(LIC_tex);
	}

	


	///* The following achieve an image-space Morse decomposition
/*
	Img_MorseDecomp temp_decomp (img_res, img_res);
	temp_decomp.set_vec_img ((unsigned char*) vec_img, min_vx, max_vx, min_vy, max_vy);
	//temp_decomp.set_fixedPt_mask((bool **) fixedPt_pixels);
	temp_decomp.set_fixedPt_mask((bool *) fixedPt_pixels);
	//temp_decomp.cal_Img_MorseDecomp(20); //20 is the default value 
	temp_decomp.cal_Img_MorseDecomp(5);  

	get_gray_MS_img(temp_decomp);
		
	printf("Highlighting valid pixels...\n");
	temp_decomp.mark_valid_SCC_pixels((unsigned char *)LIC_tex, (unsigned char *)out_LIC);

	write_ppm_flippedY("MorseSets_color.ppm", (unsigned char*)out_LIC, img_res, img_res);

	*/
/*
#ifdef MSCOMPARE

	//load the file with the triangle label information
	FILE *tfp = fopen("tri_MS_labels.txt","r");
	tri_labels = new unsigned char[poly->ntris];
	for (int i=0; i<poly->ntris; i++)
	{
		int label;
		fscanf(tfp, "%d\n", &label);
		tri_labels[i] = (unsigned char) label;
	}
	fclose(tfp);

	//render the triangle ID image
	render_tri_ID_img(poly);

	//assign type values to the individual pixels for both ISMD and the traditional MD
	assign_tri_pixel_labels(temp_decomp);

	//compute two difference images for n_I and n_II and one combined image
	compute_n_I_and_n_II();

	//output the values of n_I and n_II
	printf("n_I = %d,  n_II = %d\n", n_I, n_II);

	write_ppm_flippedY("n_I_img.ppm", (unsigned char*)n_I_img, img_res, img_res);
	write_ppm_flippedY("n_II_img.ppm", (unsigned char*)n_II_img, img_res, img_res);
	write_ppm_flippedY("combined_n_I_n_II_img.ppm", (unsigned char*)combined_n_I_n_II_img, img_res, img_res);

#endif

	printf("There are %d edges in the directed graph.\n", temp_decomp.dg->elist->nedges);
	printf("found %d Morse sets.\n", temp_decomp.valid_MS);

	
	


	/////* The following computes the ensemble of the image-space Morse decomposition
	int Ntrials = 20;
	srand( time( NULL ) );
	for (int trials = 0; trials < Ntrials; trials ++)
	{
		jitter_vectorFld (poly, 0.5, 0.1);
		render_vec_img(poly);
		Img_MorseDecomp temp_decomp (img_res, img_res);
		temp_decomp.set_vec_img ((unsigned char*) vec_img, min_vx, max_vx, min_vy, max_vy);
		//temp_decomp.set_fixedPt_mask((bool **) fixedPt_pixels);
		temp_decomp.set_fixedPt_mask((bool *) fixedPt_pixels);
		//temp_decomp.cal_Img_MorseDecomp(20); //20 is the default value 
		temp_decomp.cal_Img_MorseDecomp(5);
		store_MS_Label(temp_decomp);
		temp_decomp.finalize();
	}
	printf("counting MS hits for each pixel ...\n");
	count_MS_hits();
	printf("blending colors ...\n");
	blend_color_uncertainty(Ntrials);
	write_ppm_flippedY("MorseSets_color.ppm", (unsigned char*)out_LIC, img_res, img_res);
*/

/*
	////// Compute the ensemble of the ISMDs with different integration parameters
	int Ntrials = 20;
	srand( time( NULL ) );

	for (int trials = 0; trials < Ntrials; trials++)
	{
		Img_MorseDecomp temp_decomp (img_res, img_res);
		temp_decomp.set_vec_img ((unsigned char*) vec_img, min_vx, max_vx, min_vy, max_vy);
		temp_decomp.set_fixedPt_mask((bool *) fixedPt_pixels);

		temp_decomp.rand_integrator = rand()%3;
		temp_decomp.rand_stepsize = rand()%10+1;

		printf("rand_integrator = %d,  rand_stepsize = %d. \n", temp_decomp.rand_integrator, temp_decomp.rand_stepsize);
		temp_decomp.cal_Img_MorseDecomp(5); //10
		store_MS_Label(temp_decomp);
		temp_decomp.finalize();
	}
	printf("counting MS hits for each pixel ...\n");
	count_MS_hits();
	printf("blending colors ...\n");
	blend_color_uncertainty(Ntrials);
	write_ppm_flippedY("MorseSets_color_blended.ppm", (unsigned char*)out_LIC, img_res, img_res);
*/		


	//// For steady vector fields
/*	init_rotation_sum();
	//comp_LIC_total_rotation(total_rotation);   // when using a small L, the accurate method should be used.
	printf("computing streamlines and accumulating rotation...\n");
	comp_LIC_total_rotation_opt(total_rotation);
	printf("assigning colors for rotation field...\n");
	get_color_map_for_rot_sum();
	//get_color_map_for_counters_forward();
	printf("computing gradient of the rotation field...\n");
	comp_total_rotation_diff(total_rotation, total_rotation_diff);
	printf("assigning colors to gradient field ..\n.");
	get_color_map_for_rot_sum_diff();
		
	//write_ppm_flippedY("appro_rot_sum_color.ppm", (unsigned char*)rot_sum_color, img_res, img_res);
	//write_ppm_flippedY("appro_rot_sum_diff_color.ppm", (unsigned char*)total_rotation_diff_color, img_res, img_res);
	write_ppm_flippedY("rot_sum_color.ppm", (unsigned char*)rot_sum_color, img_res, img_res);
	write_ppm_flippedY("rot_sum_diff_color.ppm", (unsigned char*)total_rotation_diff_color, img_res, img_res);

*/
	//comp_total_rotation_diff_Hessian();
	//get_color_map_for_rot_sum_diff_Hessian();

	

/***************************************************************************

	// Compute the LIC image for a rotated vector field
	// rotate the original vector field
	gen_noise_tex ();

	rotate_vec(poly, 3.1415926/2.);
	//reflect_vec (poly);
	
	render_vec_img(poly);

	comp_LIC(LIC_tex2);

	// replace the white noise with the LIC image and do another iteration of LIC

	for (int i=0; i<1; i++)
	{
	replace_noise_tex (LIC_tex2);

	comp_LIC(LIC_tex2);
	}

	combine_LIC();

*/

	/**********************************************************************************************/
	// count the number of streamlines that passed through each pixel to find out the convergent/divergent regions
	//init_counters();
	//comp_LIC_counter (forward_counter, backward_counter);
	//get_color_map_for_counters();


	
/*	int iter = 0;
	int nFrames = 400;
	double d_time = 15./(double)(nFrames-1);  // T = 15 is the default.

	char str[512];

	for (iter = 0; iter < nFrames; iter++)
	{
		//
		double cur_st = iter * d_time;

		clock_t loc_start, loc_finish; //Used to show the time assumed
		loc_start = clock();
		gen_double_gyre_at_time(cur_st);
		//gen_ibfv_at_time(cur_st);
		//gen_gyre_saddle_at_time(cur_st);
	
		comp_LIC(LIC_tex);

		init_rotation_sum();
		init_counters();
		comp_LIC_total_rotation_timeForward(total_rotation, cur_st); 
		get_color_map_for_rot_sum();
		//get_color_map_for_counters_forward();
		comp_total_rotation_diff(total_rotation, total_rotation_diff);
		get_color_map_for_rot_sum_diff();


		compute_FTLE_forward();

		//sprintf(str, "saddle_gyre_forward_FTLE_RK2_%3i.ppm", iter);
		sprintf(str, "quadric_gyre_forward_FTLE_RK2_%3i.ppm", iter);
		write_ppm(str, (unsigned char*)FTLE_forward_color, img_res, img_res);


		//sprintf(str, "saddle_gyre_forward_rot_sum_RK2_%3i.ppm", iter);
		sprintf(str, "quadric_gyre_forward_rot_sum_RK2_%3i.ppm", iter);
		write_ppm(str, (unsigned char*)rot_sum_color, img_res, img_res);

		//sprintf(str, "saddle_gyre_forward_rot_sum_diff_RK2_%3i.ppm", iter);
		sprintf(str, "quadric_gyre_forward_rot_sum_diff_RK2_%3i.ppm", iter);
		write_ppm(str, (unsigned char*)total_rotation_diff_color, img_res, img_res);

		//sprintf(str, "double_gyre_forward_counter_RK2_%3i.ppm", iter);
		//write_ppm(str, (unsigned char*)comb_counter_color, img_res, img_res);
		loc_finish = clock(); //Get the current time after finished
		double t = (double)(loc_finish - loc_start)/CLOCKS_PER_SEC;
		printf("Time to compute this iteration is %f seconds\n", t);
	}
	
*/

/*	// The following is to show the total rotation of the streaklines starting from the centers of all pixels
	nFrames = 400;
	d_time = 10./(double)(nFrames-1);  // T = 15 is the default.

	int steps = 50;
	double it_time = .5/(double)(steps-1);

	for (iter = 0; iter < nFrames; iter ++)
	{
		double cur_start = iter * d_time;

		init_particle_list();
		insert_new_particles();
		for (int i = 0; i < steps; i++)
		{
			//double cur_t = cur_start + 2*i*it_time;
			double cur_t = cur_start + i*it_time;
			update_all_particles_streaklines (cur_t, it_time);	
			insert_new_particles();
			//update_all_particles_streaklines (cur_t+it_time, it_time);	
		}

		gen_double_gyre_at_time(cur_start);
		comp_LIC(LIC_tex);

		// compute the rotational sum
		get_rot_sum(0);
		get_color_map_for_rot_sum();
		comp_total_rotation_diff(total_rotation, total_rotation_diff);
		get_color_map_for_rot_sum_diff();
		
		sprintf(str, "double_gyre_forward_streak_rot_sum_RK2_%3i.ppm", iter);
		write_ppm(str, (unsigned char*)rot_sum_color, img_res, img_res);

		sprintf(str, "double_gyre_forward_streak_rot_sum_diff_RK2_%3i.ppm", iter);
		write_ppm(str, (unsigned char*)total_rotation_diff_color, img_res, img_res);
	}


	// we need to compute the backward tracing separately...

	finish = clock(); //Get the current time after finished
	double t = (double)(finish - start)/CLOCKS_PER_SEC;
	
	printf("Time to compute this image is %f seconds\n", t);

	////init_par_pos();
	init_particle_list_test();
	insert_new_particles_for_streakline_test();


	
    // compute the rotational sum of a 3D vector field

	comp_rot_sum_3D_steady (3000);
	get_rot_sum_3d_color_map();
    comp_rot_sum_diff_3d(); 
    get_color_map_for_rot_sum_diff_3d();
	
	//output_rot_sum_diff_result("ABC_flow_rot_sum_diff.txt");
	char filename[512];
	//sprintf(filename, "Lorenz_eq_rot_sum_%d.txt", GRIDSIZE);
	//sprintf(filename, "Rossler_eq_rot_sum_%d.txt", GRIDSIZE);
	//sprintf(filename, "General_eq1_rot_sum_%d.txt", GRIDSIZE);
	//sprintf(filename, "Chen_Lee_eq_rot_sum_%d.txt", GRIDSIZE);
	sprintf(filename, "RF_eq_rot_sum_%d.txt", GRIDSIZE);
	output_rot_sum_result(filename);

	//sprintf(filename, "Lorenz_eq_rot_sum_diff_%d.txt", GRIDSIZE);
	//sprintf(filename, "Rossler_eq_rot_sum_diff_%d.txt", GRIDSIZE);
	//sprintf(filename, "General_eq1_rot_sum_diff_%d.txt", GRIDSIZE);
	//sprintf(filename, "Chen_Lee_eq_rot_sum_diff_%d.txt", GRIDSIZE);
	sprintf(filename, "RF_eq_rot_sum_diff_%d.txt", GRIDSIZE);
	output_rot_sum_diff_result(filename);
*/

    //// Load the pre-saved 3D fields
    //load_rot_sum_result("ABC_flow_rot_sum_diff_32.txt");
    ////load_rot_sum_result("Lorenz_eq_rot_sum_diff_256.txt");
    //get_color_map_for_rot_sum_diff_3d_from_file();

	finish = clock(); //Get the current time after finished
	double t = (double)(finish - start)/CLOCKS_PER_SEC;
	printf("Time to compute these images is %f seconds\n", t);

   // setup all the user interface stuff:

	InitGlui();



	// draw the scene once and wait for some interaction:
	// (will never return)

	glutMainLoop();

	// finalize the object if loaded

	if (poly != NULL)
		poly->finalize();

	// this is here to make the compiler happy:

	return 0;
}



//
// this is where one would put code that is to be called
// everytime the glut main loop has nothing to do
//
// this is typically where animation parameters are set
//
// do not call Display() from here -- let glutMainLoop() do it
//

void
Animate( void )
{
	// put animation stuff in here -- change some global variables
	// for Display() to find:

	//advect_par(cur_ptime, d_time);
/*	update_all_particles_streaklines(cur_ptime, d_time);
	insert_new_particles_for_streakline_test();

	cur_ptime += d_time;

	if (cur_ptime >= max_time)
	{
		//init_par_pos();
		init_particle_list_test();
		cur_ptime = 0;
	}
*/

	// force a call to Display() next time it is convenient:

	glutSetWindow( MainWindow );
	glutPostRedisplay();
}




//
// glui buttons callback:
//

void
Buttons( int id )
{
	GLint viewport[4];
	char name[128];
	unsigned char *imgs = NULL;
	int width, height;

	switch( id )
	{
		case RESET:
			Reset();
			Glui->sync_live();
			glutSetWindow( MainWindow );
			glutPostRedisplay();
			break;

		case QUIT:
			// gracefully close the glui window:
			// gracefully close out the graphics:
			// gracefully close the graphics window:
			// gracefully exit the program:

			Glui->close();
			glutSetWindow( MainWindow );
			glFinish();
			glutDestroyWindow( MainWindow );
			exit( 0 );
			break;

		case SAVEIMAGE:
			if (WhichWindow == 0)   // save the main window
			{
				sprintf(name, "main_window_%0.3i.ppm", MainWindowImgNum);
				MainWindowImgNum ++;
				glutSetWindow( MainWindow );
				glGetIntegerv(GL_VIEWPORT, viewport);
				width = viewport[2];
				height = viewport[3];
				imgs = new unsigned char[3*width*height];
				glReadBuffer(GL_BACK);
				glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, imgs);
				write_ppm_flippedY(name, imgs, width, height);
			}
			else if (WhichWindow == 1)  //save the Cartesian plot
			{
				sprintf(name, "CartesianPlot_%0.3i.ppm", SecondWindowImgNum);
				SecondWindowImgNum ++;
				//glutSetWindow( SecondWindow );
				//glGetIntegerv(GL_VIEWPORT, viewport);
				//width = viewport[2];
				//height = viewport[3];
				//imgs = new unsigned char[3*width*height];
				//glReadBuffer(GL_BACK);
				//glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, imgs);
				//write_ppm_flippedY(name, imgs, width, height);
				write_ppm_flippedY(name, (unsigned char *) secondwin_img, secondwin_size, secondwin_size);
			}
			else                        //save the Polar plot
			{
				sprintf(name, "PolarPlot_%0.3i.ppm", ThirdWindowImgNum);
				ThirdWindowImgNum ++;
				//glutSetWindow( ThirdWindow );
				//glGetIntegerv(GL_VIEWPORT, viewport);
				//width = viewport[2];
				//height = viewport[3];
				//imgs = new unsigned char[3*width*height];
				//glReadBuffer(GL_BACK);
				//glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, imgs);
				//write_ppm_flippedY(name, imgs, width, height);
				write_ppm_flippedY(name, (unsigned char *) thirdwin_img, thirdwin_size, thirdwin_size);
			}
			break;

		default:
			fprintf( stderr, "Don't know what to do with Button ID %d\n", id );
	}

}



//
// draw the complete scene:
//

void
Display( void )
{
	GLsizei vx, vy, v;		// viewport dimensions
	GLint xl, yb;		// lower-left corner of viewport
	GLfloat scale2;		// real glui scale factor

	if( DebugOn != 0 )
	{
		fprintf( stderr, "Display\n" );
	}


	// set which window we want to do the graphics into:

	glutSetWindow( MainWindow );


	// erase the background:

	glDrawBuffer( GL_BACK );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable( GL_DEPTH_TEST );


	// specify shading to be flat:

	glShadeModel( GL_FLAT );


	// set the viewport to a square centered in the window:

	vx = glutGet( GLUT_WINDOW_WIDTH );
	vy = glutGet( GLUT_WINDOW_HEIGHT );
	v = vx < vy ? vx : vy;			// minimum dimension
	xl = ( vx - v ) / 2;
	yb = ( vy - v ) / 2;
	glViewport( xl, yb,  v, v );


	// set the viewing volume:
	// remember that the Z clipping  values are actually
	// given as DISTANCES IN FRONT OF THE EYE
	// USE gluOrtho2D() IF YOU ARE DOING 2D !

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	if( WhichProjection == ORTHO )
		glOrtho( -3., 3.,     -3., 3.,     0.1, 1000. );
	else
		gluPerspective( 90., 1.,	0.1, 1000. );


	// place the objects into the scene:

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();


	// set the eye position, look-at position, and up-vector:
	// IF DOING 2D, REMOVE THIS -- OTHERWISE ALL YOUR 2D WILL DISAPPEAR !

	gluLookAt( 0., 0., 3.,     0., 0., 0.,     0., 1., 0. );


	// translate the objects in the scene:
	// note the minus sign on the z value
	// this is to make the appearance of the glui z translate
	// widget more intuitively match the translate behavior
	// DO NOT TRANSLATE IN Z IF YOU ARE DOING 2D !

	glTranslatef( (GLfloat)TransXYZ[0], (GLfloat)TransXYZ[1], -(GLfloat)TransXYZ[2] );


	// rotate the scene:
	// DO NOT ROTATE (EXCEPT ABOUT Z) IF YOU ARE DOING 2D !

	glRotatef( (GLfloat)Yrot, 0., 1., 0. );
	glRotatef( (GLfloat)Xrot, 1., 0., 0. );
	glMultMatrixf( (const GLfloat *) RotMatrix );


	// uniformly scale the scene:

	glScalef( (GLfloat)Scale, (GLfloat)Scale, (GLfloat)Scale );
	scale2 = 1. + Scale2;		// because glui translation starts at 0.
	if( scale2 < MINSCALE )
		scale2 = MINSCALE;
	glScalef( (GLfloat)scale2, (GLfloat)scale2, (GLfloat)scale2 );


	// set the fog parameters:
	// DON'T NEED THIS IF DOING 2D !

	if( DepthCueOn != 0 )
	{
		glFogi( GL_FOG_MODE, FOGMODE );
		glFogfv( GL_FOG_COLOR, FOGCOLOR );
		glFogf( GL_FOG_DENSITY, FOGDENSITY );
		glFogf( GL_FOG_START, FOGSTART );
		glFogf( GL_FOG_END, FOGEND );
		glEnable( GL_FOG );
	}
	else
	{
		glDisable( GL_FOG );
	}

	// Let us disable lighting right now
	glDisable(GL_LIGHTING);


	// possibly draw the axes:

	if( AxesOn != 0 )
	{
		glColor3fv( &Colors[WhichColor][0] );
		glCallList( AxesList );
	}


	// set the color of the object:

	glColor3fv( Colors[WhichColor] );


	// draw the current object:

	// glCallList( BoxList );

	// Render the texture
	
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, 
					GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);

	// Test noise texture
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, GL_UNSIGNED_BYTE, noise_tex);

	// Test vector field image
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, GL_UNSIGNED_BYTE, vec_img);

	// Test vector field image
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
		GL_RGB, GL_UNSIGNED_BYTE, LIC_tex);

	glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex3f(-1.0, -1.0, 0.);
		glTexCoord2f(0.0,  1.0);  glVertex3f(1.0, -1.0, 0.);
		glTexCoord2f(1.0,  0.0);  glVertex3f(-1.0, 1.0, 0.);
		glTexCoord2f(1.0,  1.0);  glVertex3f(1.0, 1.0, 0.);
	glEnd();
	glDisable(GL_TEXTURE_2D);



	// draw some gratuitous text that just rotates on top of the scene:

	glDisable( GL_DEPTH_TEST );
	glColor3f( 0., 1., 1. );
	//DoRasterString( 0., 1., 0., "Text That Moves" );


	// draw some gratuitous text that is fixed on the screen:
	// the projection matrix is reset to define a scene whose
	// world coordinate system goes from 0-100 in each axis
	// this is called "percent units", and is just a convenience
	// the modelview matrix is reset to identity as we don't
	// want to transform these coordinates

	glDisable( GL_DEPTH_TEST );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluOrtho2D( 0., 100.,     0., 100. );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glColor3f( 1., 1., 1. );
	//DoRasterString( 5., 5., 0., "Text That Doesn't" );


	// Render the loaded object
	//set_view(GL_RENDER, poly);
	//
	//glTranslatef(0.0, 0.0, -3.0);

	//glTranslatef( (GLfloat)TransXYZ[0], (GLfloat)TransXYZ[1], -(GLfloat)TransXYZ[2] );
	//
	//glTranslatef(poly->center.entry[0], poly->center.entry[1], poly->center.entry[2]);

	//glScalef( (GLfloat)scale2, (GLfloat)scale2, (GLfloat)scale2 );
	//
	//glRotatef( (GLfloat)Yrot, 0., 1., 0. );
	//glRotatef( (GLfloat)Xrot, 1., 0., 0. );
	//glMultMatrixf( (const GLfloat *) RotMatrix );
	//

	//glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	//glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);

	//display_shape(GL_RENDER, poly);


	// Render the texture
	//glTranslatef( (GLfloat)TransXYZ[0], (GLfloat)TransXYZ[1], -(GLfloat)TransXYZ[2] );
	//glScalef( (GLfloat)scale2, (GLfloat)scale2, (GLfloat)scale2 );
	//
	//glRotatef( (GLfloat)Yrot, 0., 1., 0. );
	//glRotatef( (GLfloat)Xrot, 1., 0., 0. );
	//glMultMatrixf( (const GLfloat *) RotMatrix );

	// swap the double-buffered framebuffers:

	glutSwapBuffers();


	// be sure the graphics buffer has been sent:
	// note: be sure to use glFlush() here, not glFinish() !

	glFlush();
}


void  Display_3D( void )
{
	GLsizei vx, vy, v;		// viewport dimensions
	GLint xl, yb;		// lower-left corner of viewport
	GLfloat scale2;		// real glui scale factor


	// set which window we want to do the graphics into:

	glutSetWindow( MainWindow );


	// erase the background:
	glClearColor( 0., 0., 0., 1. );

	glDrawBuffer( GL_BACK );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable( GL_DEPTH_TEST );


	// specify shading to be flat:

	glShadeModel( GL_FLAT );


	// set the viewport to a square centered in the window:

	vx = glutGet( GLUT_WINDOW_WIDTH );
	vy = glutGet( GLUT_WINDOW_HEIGHT );
	v = vx < vy ? vx : vy;			// minimum dimension
	xl = ( vx - v ) / 2;
	yb = ( vy - v ) / 2;
	glViewport( xl, yb,  v, v );


	// set the viewing volume:
	// remember that the Z clipping  values are actually
	// given as DISTANCES IN FRONT OF THE EYE
	// ONLY USE gluOrtho2D() IF YOU ARE DOING 2D !

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();

	if( WhichProjection == ORTHO )
		glOrtho(  -15., 15.,     -15., 15.,     0.1, 1000.  );
	else
		gluPerspective( 90., 1.,	0.1, 1000. );


	// place the objects into the scene:

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();


	// set the eye position, look-at position, and up-vector:
	// IF DOING 2D, REMOVE THIS -- OTHERWISE ALL YOUR 2D WILL DISAPPEAR !

	gluLookAt( 0., 0., 30.,     0., 0., 0.,     0., 1., 0. );


	glTranslatef( (GLfloat)TransXYZ[0], (GLfloat)TransXYZ[1], -(GLfloat)TransXYZ[2] );



	glRotatef( (GLfloat)Yrot, 0., 1., 0. );
	glRotatef( (GLfloat)Xrot, 1., 0., 0. );
//	glMultMatrixf( (const GLfloat *) RotMatrix );


	// scale the scene:

	glScalef( (GLfloat)Scale, (GLfloat)Scale, (GLfloat)Scale );
	scale2 = 1. + Scale2;		// because glui translation starts at 0.
	if( scale2 < MINSCALE )
		scale2 = MINSCALE;
	glScalef( (GLfloat)scale2, (GLfloat)scale2, (GLfloat)scale2 );


	// possibly draw the axes:

	//if( AxesOn )
	//	glCallList( AxesList );

	/////////////////////////////////////////////////
	////Project 8 drawing

	////Determine which sides of the cube are facing the user

	if (Disable3DTexture==0)
	{
    DetermineVisibility();
	FillXY();
	FillXZ();
	FillYZ();

	Drawing();
	}

	// Disable texture and draw 3D streamlines
	if (ShowSelStreamline==1)
	{
		int i, j;
		for (i=0; i<list_streamlines_3d.size(); i++)
		{
			int ii=list_streamlines_3d[i].i;
			int jj=list_streamlines_3d[i].j;
			int kk=list_streamlines_3d[i].k;

			//rot_sum_diff_3d_color
			glColor3fv(rot_sum_diff_rgb[ii][jj][kk]);
			Streamline3D &cur_line = list_streamlines_3d[i];
			glBegin(GL_LINE_STRIP);
			for (j=0; j<cur_line.one_line.size(); j++)
			{
				float pos[3]={cur_line.one_line[j].xyz[0], cur_line.one_line[j].xyz[1], cur_line.one_line[j].xyz[2]};

				//pos[0] = ((pos[0]-data_center.entry[0])-MIN_X)/(MAX_X-MIN_X)*10-10;
				//pos[1] = ((pos[1]-data_center.entry[1])-MIN_Y)/(MAX_Y-MIN_Y)*10-10;
				//pos[2] = ((pos[2]-data_center.entry[2])-MIN_Z)/(MAX_Z-MIN_Z)*10-10;
				pos[0] = ((pos[0]-data_center.entry[0])-MIN_X)/(MAX_X-MIN_X)*20-10;
				pos[1] = ((pos[1]-data_center.entry[1])-MIN_Y)/(MAX_Y-MIN_Y)*20-10;
				pos[2] = ((pos[2]-data_center.entry[2])-MIN_Z)/(MAX_Z-MIN_Z)*20-10;

				glVertex3fv (pos);
				//glVertex3fv (cur_line.one_line[j].xyz);
			}
			glEnd();
		}
	}

	// swap the double-buffered framebuffers:
	glutSwapBuffers();

	// be sure the graphics buffer has been sent:
	glFlush();
}



void 
display_func()
{
	if (display_3D == 0)
		render_LIC_img();

	else
		Display_3D();
}

//
// use glut to display a string of characters using a raster font:
//

void
DoRasterString( float x, float y, float z, /*unsigned*/ char *s )
{
	char c;			// one character to print

	glRasterPos3f( (GLfloat)x, (GLfloat)y, (GLfloat)z );
	for( ; ( c = *s ) != '\0'; s++ )
	{
		glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_24, c );
	}
}



//
// use glut to display a string of characters using a stroke font:
//

void
DoStrokeString( float x, float y, float z, float ht, /*unsigned*/ char *s )
{
	char c;			// one character to print
	float sf;		// the scale factor

	glPushMatrix();
		glTranslatef( (GLfloat)x, (GLfloat)y, (GLfloat)z );
		sf = ht / ( 119.05 + 33.33 );
		glScalef( (GLfloat)sf, (GLfloat)sf, (GLfloat)sf );
		for( ; ( c = *s ) != '\0'; s++ )
		{
			glutStrokeCharacter( GLUT_STROKE_ROMAN, c );
		}
	glPopMatrix();
}



//
// return the number of seconds since the start of the program:
//

float
ElapsedSeconds( void )
{
	// get # of milliseconds since the start of the program:

	int ms = glutGet( GLUT_ELAPSED_TIME );

	// convert it to seconds:

	return (float)ms / 1000.;
}

void Sliders( int id )
{
        char str[2345];

        switch( id )
        {

            //case TEMP:
            //     sprintf( str, "T: the current value is from %5.2f to %5.2f \n", TempLowHigh[0], TempLowHigh[1] );
            //     TempLabel->set_text( str );
            //     break;
			
			case XRANGE:
				 sprintf( str, "Maximun opacity:  %5.2f\n ", MaxAlpha );
                 XRangeLabel->set_text( str );
                 break;
			
			//case YRANGE:
			//	 sprintf( str, "Y: the current value is from %5.2f to %5.2f \n", YLowHigh[0], YLowHigh[1] );
   //              YRangeLabel->set_text( str );
   //              break;
			//
			//case ZRANGE:
			//	 sprintf( str, "Z: the current value is from %5.2f to %5.2f \n", ZLowHigh[0], ZLowHigh[1] );
   //              ZRangeLabel->set_text( str );
   //              break;
			
			case SRANGE:
				 sprintf( str, "SRange: the current value is from %5.2f to %5.2f \n", SRange[0], SRange[1] );
				 // compute a set of streamlines
				 if (ShowSelStreamline == 1)
				 {
					 update_list_streamlines_3d();
				 }
                 SRangeLabel->set_text( str );
                 break;

			default:
				fprintf(stderr," Wrong slider operation!\n");
		}

        glutSetWindow( MainWindow );
        glutPostRedisplay();
}

//
// initialize the glui window:
//

void
InitGlui( void )
{
	GLUI_Panel *panel;
	GLUI_RadioGroup *group;
	GLUI_Rotation *rot;
	GLUI_Translation *trans, *scale;


	// setup the glui window:

	glutInitWindowPosition( INIT_WINDOW_SIZE + 50, 0 );
	Glui = GLUI_Master.create_glui( (char *) GLUITITLE );


	Glui->add_statictext( (char *) GLUITITLE );
	Glui->add_separator();
	
	Glui->add_checkbox( "Display 3D", &display_3D );
	Glui->add_checkbox( "Axes", &AxesOn );
	Glui->add_checkbox( "Perspective", &WhichProjection );
	Glui->add_checkbox( "Intensity Depth Cue", &DepthCueOn );

	// Add a rollout for the axes color
	//GLUI_Rollout *rollout = Glui->add_rollout(" Axes Color ", 0);

	//panel = Glui->add_panel(  "Axes Color" );
	//GLUI_Rollout *rollout = Glui->add_rollout_to_panel(panel,  "Axes Color", 1 );
		//group = Glui->add_radiogroup_to_panel( panel, &WhichColor );

		//group = Glui->add_radiogroup_to_panel( rollout, &WhichColor );
		//	Glui->add_radiobutton_to_group( group, "Red" );
		//	Glui->add_radiobutton_to_group( group, "Yellow" );
		//	Glui->add_radiobutton_to_group( group, "Green" );
		//	Glui->add_radiobutton_to_group( group, "Cyan" );
		//	Glui->add_radiobutton_to_group( group, "Blue" );
		//	Glui->add_radiobutton_to_group( group, "Magenta" );
		//	Glui->add_radiobutton_to_group( group, "White" );
		//	Glui->add_radiobutton_to_group( group, "Black" );


	// Add a list for the different models
	//rollout = Glui->add_rollout(" Models ", 0);
	//panel = Glui->add_panel(  "Choose object to open " );
	//		GLUI_Listbox *obj_list = Glui->add_listbox_to_panel(panel, "Objects", &ObjectId, -1, ( GLUI_Update_CB) Choose_Object);
	//		obj_list->add_item (0, "bunny");
	//		obj_list->add_item (1, "feline");
	//		obj_list->add_item (2, "dragon");
	//		obj_list->add_item (3, "happy");
	//		obj_list->add_item (4, "sphere");
	//		obj_list->add_item (5, "torus");
	
	Glui->add_checkbox( "Disable 3D texture", &Disable3DTexture );
	Glui->add_checkbox( "Show the difference", &Show_ROTSUM_or_DIFF );
	Glui->add_checkbox( "Show only LIC", &Show_OnlyLIC );
	Glui->add_checkbox( "Show singularities", &ShowSingularities);
	Glui->add_checkbox( "Show selected integral curve", &ShowSelStreamline);
	Glui->add_checkbox( "Use RK23", &UseRK23);
	
	Glui->add_separator();

	char str[512];

	ScaleSlider =  Glui->add_slider( true, GLUI_HSLIDER_FLOAT, SRange,
		SRANGE, (GLUI_Update_CB) Sliders );
	ScaleSlider->set_float_limits( -10.,10. ); //the range should be changed here
	ScaleSlider->set_w( 300 );
	sprintf( str, "SRange: the current value is from %5.2f to %5.2f \n", SRange[0], SRange[1] );
	SRangeLabel = Glui->add_statictext( str );
    
	AlphaSlider= Glui->add_slider( false, GLUI_HSLIDER_FLOAT, &MaxAlpha,
		XRANGE, (GLUI_Update_CB) Sliders );
	AlphaSlider->set_float_limits( 0.,1. ); //the range should be changed here
	AlphaSlider->set_w( 300 );
	sprintf( str, "Maximun opacity: %5.2f \n", MaxAlpha );
	XRangeLabel = Glui->add_statictext( str );
	Glui->add_separator();

	Glui->add_checkbox( "Show arrow plot", &ArrowsOn );
	Arrow_Scale_Spinner = Glui->add_spinner("Arrow scale", GLUI_SPINNER_FLOAT, &arrow_scale);
	Arrow_Scale_Spinner->set_float_val(0.05);
	Arrow_Scale_Spinner->set_speed(0.01);
	Glui->add_separator();

	
	Glui->add_checkbox( "Show rotation difference", &RotSumDiffOn );
	Glui->add_separator();

	panel = Glui->add_panel( "Object Transformation" );

		rot = Glui->add_rotation_to_panel( panel, "Rotation", (float *) RotMatrix );

		// allow the object to be spun via the glui rotation widget:

		rot->set_spin( 1.0 );


		Glui->add_column_to_panel( panel, GLUIFALSE );
		scale = Glui->add_translation_to_panel( panel, "Scale",  GLUI_TRANSLATION_Y , &Scale2 );
		scale->set_speed( 0.005f );

		Glui->add_column_to_panel( panel, GLUIFALSE );
		trans = Glui->add_translation_to_panel( panel, "Trans XY", GLUI_TRANSLATION_XY, &TransXYZ[0] );
		trans->set_speed( 0.05f );

		Glui->add_column_to_panel( panel, GLUIFALSE );
		trans = Glui->add_translation_to_panel( panel, "Trans Z",  GLUI_TRANSLATION_Z , &TransXYZ[2] );
		trans->set_speed( 0.05f );

	Glui->add_checkbox( "Debug", &DebugOn );

	Glui->add_checkbox( "Freeze previous plot", &FreezeOn );
	
	Glui->add_checkbox( "Show rotation sign", &ShowRotSignOn );
	
	Glui->add_checkbox( "Show local rotation sign", &ShowLocalRotSignOn );

	Glui->add_checkbox( "Show rotation sign", &ShowGroupColorsOn );

	panel = Glui->add_panel( "", GLUIFALSE );

	Glui->add_button_to_panel(panel, "Save images", SAVEIMAGE, (GLUI_Update_CB) Buttons );
	group = Glui->add_radiogroup_to_panel( panel, &WhichWindow );
			Glui->add_radiobutton_to_group( group, "Main window" );
			Glui->add_radiobutton_to_group( group, "Cartesian plot" );
			Glui->add_radiobutton_to_group( group, "Polar plot" );

	Glui->add_separator();

	panel = Glui->add_panel( "", GLUIFALSE );

	Glui->add_button_to_panel( panel, "Reset", RESET, (GLUI_Update_CB) Buttons );

	Glui->add_column_to_panel( panel, GLUIFALSE );

	Glui->add_button_to_panel( panel, "Quit", QUIT, (GLUI_Update_CB) Buttons );


	// tell glui what graphics window it needs to post a redisplay to:

	Glui->set_main_gfx_window( MainWindow );


	// set the graphics window's idle function:

	GLUI_Master.set_glutIdleFunc( NULL );
	//GLUI_Master.set_glutIdleFunc( Animate );
}



//
// initialize the glut and OpenGL libraries:
//	also setup display lists and callback functions
//

void
InitGraphics( void )
{
	// setup the display mode:
	// ( *must* be done before call to glutCreateWindow() )
	// ask for color, double-buffering, and z-buffering:

	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );


	// set the initial window configuration:

	glutInitWindowPosition( 0, 0 );
	glutInitWindowSize( img_res/*INIT_WINDOW_SIZE*/, img_res/*INIT_WINDOW_SIZE*/ );


	// open the window and set its title:

	MainWindow = glutCreateWindow( WINDOWTITLE );
	glutSetWindowTitle( WINDOWTITLE );


	// setup the clear values:

	glClearColor( BACKCOLOR[0], BACKCOLOR[1], BACKCOLOR[2], BACKCOLOR[3] );


	// setup the callback routines:


	// DisplayFunc -- redraw the window
	// ReshapeFunc -- handle the user resizing the window
	// KeyboardFunc -- handle a keyboard input
	// MouseFunc -- handle the mouse button going down or up
	// MotionFunc -- handle the mouse moving with a button down
	// PassiveMotionFunc -- handle the mouse moving with a button up
	// VisibilityFunc -- handle a change in window visibility
	// EntryFunc	-- handle the cursor entering or leaving the window
	// SpecialFunc -- handle special keys on the keyboard
	// SpaceballMotionFunc -- handle spaceball translation
	// SpaceballRotateFunc -- handle spaceball rotation
	// SpaceballButtonFunc -- handle spaceball button hits
	// ButtonBoxFunc -- handle button box hits
	// DialsFunc -- handle dial rotations
	// TabletMotionFunc -- handle digitizing tablet motion
	// TabletButtonFunc -- handle digitizing tablet button hits
	// MenuStateFunc -- declare when a pop-up menu is in use
	// TimerFunc -- trigger something to happen a certain time from now
	// IdleFunc -- what to do when nothing else is going on

	glutSetWindow( MainWindow );
	//glutDisplayFunc( Display );
	glutDisplayFunc( display_func /*render_LIC_img*/ );
	//glutDisplayFunc( Display_Model );render_LIC_img()
	glutReshapeFunc( Resize );
	glutKeyboardFunc( Keyboard );
	glutMouseFunc( MouseButton );
	glutMotionFunc( MouseMotion );
	glutPassiveMotionFunc( NULL );
	glutVisibilityFunc( Visibility );
	glutEntryFunc( NULL );
	glutSpecialFunc( NULL );
	glutSpaceballMotionFunc( NULL );
	glutSpaceballRotateFunc( NULL );
	glutSpaceballButtonFunc( NULL );
	glutButtonBoxFunc( NULL );
	glutDialsFunc( NULL );
	glutTabletMotionFunc( NULL );
	glutTabletButtonFunc( NULL );
	glutMenuStateFunc( NULL );
	glutTimerFunc( 0, NULL, 0 );

	// DO NOT SET THE GLUT IDLE FUNCTION HERE !!
	// glutIdleFunc( NULL );
	// let glui take care of it in InitGlui()
}




//
// initialize the display lists that will not change:
//

void
InitLists( void )
{
	float dx = BOXSIZE / 2.;
	float dy = BOXSIZE / 2.;
	float dz = BOXSIZE / 2.;

	// create the object:

	BoxList = glGenLists( 1 );
	glNewList( BoxList, GL_COMPILE );

		glBegin( GL_QUADS );

			glColor3f( 0., 0., 1. );
			glNormal3f( 0., 0.,  1. );
				glVertex3f( -dx, -dy,  dz );
				glVertex3f(  dx, -dy,  dz );
				glVertex3f(  dx,  dy,  dz );
				glVertex3f( -dx,  dy,  dz );

			glNormal3f( 0., 0., -1. );
				glTexCoord2f( 0., 0. );
				glVertex3f( -dx, -dy, -dz );
				glTexCoord2f( 0., 1. );
				glVertex3f( -dx,  dy, -dz );
				glTexCoord2f( 1., 1. );
				glVertex3f(  dx,  dy, -dz );
				glTexCoord2f( 1., 0. );
				glVertex3f(  dx, -dy, -dz );

			glColor3f( 1., 0., 0. );
			glNormal3f(  1., 0., 0. );
				glVertex3f(  dx, -dy,  dz );
				glVertex3f(  dx, -dy, -dz );
				glVertex3f(  dx,  dy, -dz );
				glVertex3f(  dx,  dy,  dz );

			glNormal3f( -1., 0., 0. );
				glVertex3f( -dx, -dy,  dz );
				glVertex3f( -dx,  dy,  dz );
				glVertex3f( -dx,  dy, -dz );
				glVertex3f( -dx, -dy, -dz );

			glColor3f( 0., 1., 0. );
			glNormal3f( 0.,  1., 0. );
				glVertex3f( -dx,  dy,  dz );
				glVertex3f(  dx,  dy,  dz );
				glVertex3f(  dx,  dy, -dz );
				glVertex3f( -dx,  dy, -dz );

			glNormal3f( 0., -1., 0. );
				glVertex3f( -dx, -dy,  dz );
				glVertex3f( -dx, -dy, -dz );
				glVertex3f(  dx, -dy, -dz );
				glVertex3f(  dx, -dy,  dz );

		glEnd();

	glEndList();


	// create the axes:

	AxesList = glGenLists( 1 );
	glNewList( AxesList, GL_COMPILE );
		glLineWidth( AXES_WIDTH );
			Axes( 1.5 );
		glLineWidth( 1. );
	glEndList();
}



//
// the keyboard callback:
//

void
Keyboard( unsigned char c, int x, int y )
{
	if( DebugOn != 0 )
		fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );

	switch( c )
	{
		case 'o':
		case 'O':
			WhichProjection = ORTHO;
			break;

		case 'p':
		case 'P':
			WhichProjection = PERSP;
			break;

		case 'q':
		case 'Q':
		case ESCAPE:
			Buttons( QUIT );	// will not return here
			break;			// happy compiler

		case 'r':
		case 'R':
			LeftButton = ROTATE;
			break;

		case 's':
		case 'S':
			LeftButton = SCALE;
			break;

		default:
			fprintf( stderr, "Don't know what to do with keyboard hit: '%c' (0x%0x)\n", c, c );
	}


	// synchronize the GLUI display with the variables:

	Glui->sync_live();


	// force a call to Display():

	glutSetWindow( MainWindow );
	glutPostRedisplay();
}



//
// called when the mouse button transitions down or up:
//

void
MouseButton( int button, int state, int x, int y )
{
	int b;			// LEFT, MIDDLE, or RIGHT

	if( DebugOn != 0 )
		fprintf( stderr, "MouseButton: %d, %d, %d, %d\n", button, state, x, y );

	
	// get the proper button bit mask:

	switch( button )
	{
		case GLUT_LEFT_BUTTON:
			b = LEFT;		break;

		case GLUT_MIDDLE_BUTTON:
			b = MIDDLE;		break;

		case GLUT_RIGHT_BUTTON:
			b = RIGHT;		break;

		default:
			b = 0;
			fprintf( stderr, "Unknown mouse button: %d\n", button );
	}


	// button down sets the bit, up clears the bit:

	if( state == GLUT_DOWN )
	{
		Xmouse = x;
		Ymouse = y;
		ActiveButton |= b;		// set the proper bit

		// trace out a streamline
		trace_streamline_from_image_based(img_res-y, x, L);
		merge_current_forward_back_streamlines();

		//Display_func_win2();

		sel_pos_i = img_res-y;
		sel_pos_j = x;
	}
	else
	{
		ActiveButton &= ~b;		// clear the proper bit
	}

	glutSetWindow( MainWindow );
	glutPostRedisplay();

#ifdef SHOWSECONDTHIRDWINS
	glutSetWindow( SecondWindow );
	glutPostRedisplay();

	glutSetWindow( ThirdWindow );
	glutPostRedisplay();
#endif
}



//
// called when the mouse moves while a button is down:
//

void
MouseMotion( int x, int y )
{
	int dx, dy;		// change in mouse coordinates

	if( DebugOn != 0 )
		fprintf( stderr, "MouseMotion: %d, %d\n", x, y );


	dx = x - Xmouse;		// change in mouse coords
	dy = y - Ymouse;

	if( ( ActiveButton & LEFT ) != 0 )
	{
		switch( LeftButton )
		{
			case ROTATE:
				Xrot += ( ANGFACT*dy );
				Yrot += ( ANGFACT*dx );
				break;

			case SCALE:
				Scale += SCLFACT * (float) ( dx - dy );
				if( Scale < MINSCALE )
					Scale = MINSCALE;
				break;
		}
	}


	if( ( ActiveButton & MIDDLE ) != 0 )
	{
		Scale += SCLFACT * (float) ( dx - dy );

		// keep object from turning inside-out or disappearing:

		if( Scale < MINSCALE )
			Scale = MINSCALE;
	}

	Xmouse = x;			// new current position
	Ymouse = y;

	glutSetWindow( MainWindow );
	glutPostRedisplay();
	
	glutSetWindow( SecondWindow );
	glutPostRedisplay();
}



//
// reset the transformations and the colors:
//
// this only sets the global variables --
// the glut main loop is responsible for redrawing the scene
//

void
Reset( void )
{
	ActiveButton = 0;
	AxesOn = GLUITRUE;
	DebugOn = GLUIFALSE;
	DepthCueOn = GLUIFALSE;
	LeftButton = ROTATE;
	Scale  = 1.0;
	Scale2 = 0.0;		// because we add 1. to it in Display()
	WhichColor = WHITE;
	WhichProjection = PERSP;
	Xrot = Yrot = 0.;
	TransXYZ[0] = TransXYZ[1] = TransXYZ[2] = 0.;

	                  RotMatrix[0][1] = RotMatrix[0][2] = RotMatrix[0][3] = 0.;
	RotMatrix[1][0]                   = RotMatrix[1][2] = RotMatrix[1][3] = 0.;
	RotMatrix[2][0] = RotMatrix[2][1]                   = RotMatrix[2][3] = 0.;
	RotMatrix[3][0] = RotMatrix[3][1] = RotMatrix[3][3]                   = 0.;
	RotMatrix[0][0] = RotMatrix[1][1] = RotMatrix[2][2] = RotMatrix[3][3] = 1.;

	RotSumDiffOn = 0;
	ArrowsOn = 0;
}



//
// called when user resizes the window:
//

void
Resize( int width, int height )
{
	if( DebugOn != 0 )
		fprintf( stderr, "ReSize: %d, %d\n", width, height );

	// don't really need to do anything since window size is
	// checked each time in Display():

	glutSetWindow( MainWindow );
	glutPostRedisplay();
}


//
// handle a change to the window's visibility:
//

void
Visibility ( int state )
{
	if( DebugOn != 0 )
		fprintf( stderr, "Visibility: %d\n", state );

	if( state == GLUT_VISIBLE )
	{
		glutSetWindow( MainWindow );
		glutPostRedisplay();
	}
	else
	{
		// could optimize by keeping track of the fact
		// that the window is not visible and avoid
		// animating or redrawing it ...
	}
}


////////////////////////////////////////////////////////////
/// Keyboard and mouse functions for the second window

void
Keyboard_win2( unsigned char c, int x, int y )
{
	if( DebugOn != 0 )
		fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );

	switch( c )
	{
		case 'o':
		case 'O':
			WhichProjection = ORTHO;
			break;

		case 'p':
		case 'P':
			WhichProjection = PERSP;
			break;

		case 'q':
		case 'Q':
		case ESCAPE:
			Buttons( QUIT );	// will not return here
			break;			// happy compiler

		case 'r':
		case 'R':
			LeftButton = ROTATE;
			break;

		case 's':
		case 'S':
			LeftButton = SCALE;
			break;

		default:
			fprintf( stderr, "Don't know what to do with keyboard hit: '%c' (0x%0x)\n", c, c );
	}


	// synchronize the GLUI display with the variables:

	Glui->sync_live();


	// force a call to Display():

	glutSetWindow( SecondWindow );
	glutPostRedisplay();
}



//
// called when the mouse button transitions down or up:
//

void
MouseButton_win2( int button, int state, int x, int y )
{
	int b;			// LEFT, MIDDLE, or RIGHT

	if( DebugOn != 0 )
		fprintf( stderr, "MouseButton: %d, %d, %d, %d\n", button, state, x, y );

	
	// get the proper button bit mask:

	switch( button )
	{
		case GLUT_LEFT_BUTTON:
			b = LEFT;		break;

		case GLUT_MIDDLE_BUTTON:
			b = MIDDLE;		break;

		case GLUT_RIGHT_BUTTON:
			b = RIGHT;		break;

		default:
			b = 0;
			fprintf( stderr, "Unknown mouse button: %d\n", button );
	}


	// button down sets the bit, up clears the bit:

	if( state == GLUT_DOWN )
	{
		Xmouse = x;
		Ymouse = y;
		ActiveButton |= b;		// set the proper bit

		// trace out a streamline
		trace_streamline_from_image_based(img_res-y, x, 1000);

		sel_pos_i = img_res-y;
		sel_pos_j = x;
	}
	else
	{
		ActiveButton &= ~b;		// clear the proper bit
	}

	glutSetWindow( SecondWindow );
	glutPostRedisplay();
}



//
// called when the mouse moves while a button is down:
//

void
MouseMotion_win2( int x, int y )
{
	int dx, dy;		// change in mouse coordinates

	if( DebugOn != 0 )
		fprintf( stderr, "MouseMotion: %d, %d\n", x, y );


	dx = x - Xmouse;		// change in mouse coords
	dy = y - Ymouse;

	if( ( ActiveButton & LEFT ) != 0 )
	{
		switch( LeftButton )
		{
			case ROTATE:
				Xrot += ( ANGFACT*dy );
				Yrot += ( ANGFACT*dx );
				break;

			case SCALE:
				Scale += SCLFACT * (float) ( dx - dy );
				if( Scale < MINSCALE )
					Scale = MINSCALE;
				break;
		}
	}


	if( ( ActiveButton & MIDDLE ) != 0 )
	{
		Scale += SCLFACT * (float) ( dx - dy );

		// keep object from turning inside-out or disappearing:

		if( Scale < MINSCALE )
			Scale = MINSCALE;
	}

	Xmouse = x;			// new current position
	Ymouse = y;

	glutSetWindow( SecondWindow );
	glutPostRedisplay();
}


void 
draw_line_segment(float start[2], float end[2])
{
	glBegin(GL_LINES);
	glVertex2fv (start);
	glVertex2fv (end);
	glEnd();
}


void 
Display_func_win2()
{
	// generate a 2D rotation plot for each selected streamline
	// later, a Polar Coordinate can be shown as well

	// First, draw the X and Y axes of the plot\

	if (FreezeOn == 0)
	{
	glClearColor( 1., 1., 1., 1. );
	glViewport(0, 0, (GLsizei) secondwin_size, (GLsizei) secondwin_size);
	glClear(GL_COLOR_BUFFER_BIT);
	}
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	//gluOrtho2D(-1., 1, -1., 1);
	gluOrtho2D(-0.2, 1.3, -1.2, 1.2);
	glEnable( GL_LINE_SMOOTH );
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// First, cast some shadow
	float start[2] = {0., 0.};
	float end[2] = {1.2, 0.};
	float offset = 0.003;
	float start2[2] = {start[0]+offset, start[1]-offset};
	float end2[2]={end[0]+offset, end[1]-offset};
	float direct[2] = {1., 0.};

	if (FreezeOn == 0)
	{
	glLineWidth(5.);
	glColor3f (0.8, 0.8, 0.8);
	draw_line_segment(start2, end2);
	draw_arrow_head(end2, direct); // cast shadow for X
	
	start2[1] = -1.1-offset;
	end[0] = 0.;
	end[1] = 1.1;
	end2[0] = end[0]+offset;
	end2[1] = end[1]-offset;
	draw_line_segment(start2, end2); 
	
	direct[0] = 0.;
	direct[1] = 1.;
	draw_arrow_head(end2, direct); // cast shadow for Y

	glLineWidth(3.);
	glColor3f (0, 0, 0);
	// Draw X axis
	end[0] = 1.2;
	end[1] = 0.;
	draw_line_segment(start, end);
	//Draw arrow head
	direct[0] = 1.;
	direct[1] = 0.;
	draw_arrow_head(end, direct);

	// Draw Y axis
	start[1] = -1.1;
	end[0] = 0.;
	end[1] = 1.1;
	glColor3f (0, 0, 0);
	draw_line_segment(start, end);
	//Draw arrow head
	direct[0] = 0.;
	direct[1] = 1.;
	draw_arrow_head(end, direct);


	// Let us draw a line to correspond to the value=2pi

	start[0] = -0.02;
	start[1] = 0.5;

	end[0] = 1.2;
	end[1] = 0.5;
	glColor3f (0.9, 0.3, 0.5);
	glLineWidth(1.);
	draw_line_segment(start, end);

	start[0] = -0.02;
	start[1] = -0.5;

	end[0] = 1.2;
	end[1] = -0.5;
	glColor3f (0.3, 0.9, 0.5);
	glLineWidth(1.);
	draw_line_segment(start, end);

	// Draw a sequence of lines corresponding to mark different arc lengths
	glColor3f (0.7, 0.7, 0.7);
	int Nlines = 5;
	double dx = 1.1/(Nlines-1);
	start[0] = end[0] = 1.1;
	start[1] = -1.;
	end[1] = 1.;

	for (int i=1; i<Nlines; i++)
	{
		start[0] = end[0] = i*dx;
		draw_line_segment(start, end);
	}

	/*unsigned*/ 
	glColor3f (0, 0, 0);
	char str[10];
	str[0] = '0';
	//DoRasterString( -0.04, -0.1, 0., str );
	DoStrokeString(-0.04, -0.06, 0., 0.06, str);
	str[0] = '1';
	//DoRasterString( 1.1, -0.1, 0., str );
	DoStrokeString(1.1, -0.06, 0., 0.06, str);
	
	float legend_dx = 1./(Nlines-1);
	for (int i=1; i<Nlines-1; i++)
	{
		float sx = i*dx;
		sprintf(str, "%.2f", legend_dx*i);
		DoStrokeString(sx, -0.06, 0, 0.06, str);
	}
	
	str[0]='2';
	str[1]='?';
	
	//DoRasterString( -0.15, 0.5, 0., "360" );
	DoStrokeString(-0.12, 0.5, 0., 0.06, "360");
	//sprintf(str, "-2");
	str[0]='-';
	str[1]='2';
	str[2]=227;
	//DoRasterString( -0.15, -0.5, 0., "-360" );
	DoStrokeString(-0.15, -0.5, 0., 0.06, "-360");
	}

	// Now display the trend of the rotation along a streamline

	if (current_total_rot.size()>1)
	{
		double interval = 1.1/(current_total_rot.size()-1);
		int i;
		// Note that the 2pi is now a horizontal line with y=0.8
		glLineWidth(5.);
		glColor3f (0.8, 0.8, 0.8);
		offset = 0.001;
		glBegin(GL_LINE_STRIP);
		for (i=0; i<current_total_rot.size(); i++)
		{
			float x = i*interval+offset;
			float y = 0.5*current_total_rot[i]/(2*PI)-offset;
			glVertex2f (x, y);
		}
		glEnd();

		glLineWidth(3.);
		glColor3f (0.3, 0.5, 0.99);
		glBegin(GL_LINE_STRIP);
		for (i=0; i<current_total_rot.size(); i++)
		{
			float x = i*interval;
			float y = 0.5*current_total_rot[i]/(2*PI);

			if (ShowRotSignOn == 1)
			{
				if (y>0)
					glColor3f (0.9, 0.3, 0.5);
				else
					glColor3f (0.3, 0.9, 0.5);


				if (ShowLocalRotSignOn && i>0)
				{
					//
					float rot_orient = current_total_rot[i]-current_total_rot[i-1];

					if (rot_orient>0)
						glColor3f (0.9, 0.3, 0.5);
					else
						glColor3f (0.3, 0.9, 0.5);
				}
			}
			glVertex2f (x, y);
		}
		glEnd();
	}

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, secondwin_size, secondwin_size, GL_RGB, GL_UNSIGNED_BYTE, secondwin_img);
	//merge_current_forward_back_streamlines();

	glutSwapBuffers();

}



////////////////////////////////////////////////////////////////////////////////
// Display the rotation plot in the Polar coordinate system

void draw_hollow_circle2(double cx, double cy, double R)
{
	int i;
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

void 
draw_radial_lines (double cx, double cy, double R, double Nlines)
{
	int i;
	double theta, deta;
	deta = 2* PI/(Nlines-1);
	double x, y;
	theta = 0.;

	glColor3f (0.7, 0.7, 0.7);
	glBegin(GL_LINES);
	for(i=0; i<Nlines; i++, theta += deta)
	{
		x = cx + R * cos(theta);
		y = cy + R * sin(theta);

		glVertex2f (cx, cy);
		glVertex2f (x, y);

	}
	glEnd();

	glColor3f (0, 0, 0);
	theta = 0.;
	for(i=0; i<Nlines; i++, theta += deta)
	{
		x = cx + R * cos(theta);
		y = cy + R * sin(theta);

		char str[20];
		sprintf(str, "%.1f", theta*180./PI);

		if (abs (theta*180./PI-360.)<1e-5)
		{
			sprintf(str, "(%.1f)", theta*180./PI);
			DoStrokeString( x+0.07, y-0.015, 0., 0.035, str );
		}

		else if (abs(theta-0.0)<1.e-5)
			DoStrokeString( x+0.02, y-0.015, 0., 0.035, str );
		else if (abs(theta*180./PI-90.0)<1.e-5)
			DoStrokeString( x-0.02, y+0.015, 0., 0.035, str );
		else if (abs (theta*180./PI-180.)<1e-5)
			DoStrokeString( x-0.08, y-0.015, 0., 0.035, str );
		else if (abs (theta*180./PI-270.)<1e-5)
			DoStrokeString( x-0.03, y-0.035, 0., 0.035, str );

		else if (x >0 && y >0)
			DoStrokeString( x+0.02, y+0.01, 0., 0.035, str );
		else if (x < 0  && y > 0)
			DoStrokeString( x-0.08, y+0.01, 0., 0.035, str );
		else if (x < 0 && y < 0)
			DoStrokeString( x-0.08, y-0.025, 0., 0.035, str );
		else
			DoStrokeString( x+0.02, y-0.025, 0., 0.035, str );
	}
}

void 
Display_func_win3()
{
	if (FreezeOn == 0)
	{
	glClearColor( 1., 1., 1., 1. );
	glViewport(0, 0, (GLsizei) thirdwin_size, (GLsizei) thirdwin_size);
	glClear(GL_COLOR_BUFFER_BIT);
	}
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(-0.6, 0.7, -0.6, 0.7);
	glEnable( GL_LINE_SMOOTH );
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Draw a number of concentric circles

	if (FreezeOn == 0)
	{
	glLineWidth(2.);
	glColor3f (0.5, 0.4, 0.7);
	draw_hollow_circle2(0, 0, 0.1);
	draw_hollow_circle2(0, 0, 0.2);
	draw_hollow_circle2(0, 0, 0.3);
	draw_hollow_circle2(0, 0, 0.4);
	draw_hollow_circle2(0, 0, 0.5);

	DoStrokeString( 0.103, -0.033, 0., 0.035, "0.2" );
	DoStrokeString( 0.203, -0.033, 0., 0.035, "0.4" );
	DoStrokeString( 0.303, -0.033, 0., 0.035, "0.6" );
	DoStrokeString( 0.403, -0.033, 0., 0.035, "0.8" );
	DoStrokeString( 0.483, -0.033, 0., 0.035, "1." );

	glColor3f (0.7, 0.7, 0.7);
	draw_radial_lines (0, 0, 0.5, 17);
	}

	// map the rotational data into the Polar coordinate
	if (current_total_rot.size()>1)
	{
		double interval = 0.5/(current_total_rot.size()-1);
		double offset = 0.001;
		
		int i;
		// Note that the 2pi is now a horizontal line with y=0.8
		glLineWidth(4.);
		glColor3f (0.8, 0.8, 0.8);
		glBegin(GL_LINE_STRIP);
		for (i=0; i<current_total_rot.size(); i++)
		{
			float R = i*interval+offset;
			float x = R * cos (current_total_rot[i]);
			float y = R * sin (current_total_rot[i]);
			glVertex2f (x+offset, y-offset);
		}
		glEnd();

		glLineWidth(3.);
		glColor3f (0.3, 0.5, 0.99);
		glBegin(GL_LINE_STRIP);
		for (i=0; i<current_total_rot.size(); i++)
		{
			float R = i*interval+offset;
			float x = R * cos (current_total_rot[i]);
			float y = R * sin (current_total_rot[i]);
			if (ShowRotSignOn == 1)
			{
				if (current_total_rot[i]>0)
					glColor3f (0.9, 0.3, 0.5);
				else
					glColor3f (0.3, 0.9, 0.5);
			}
			glVertex2f (x, y);
		}
		glEnd();
	}
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, thirdwin_size, thirdwin_size, GL_RGB, GL_UNSIGNED_BYTE, thirdwin_img);

	glutSwapBuffers();
}


//////////////////////////////////////////  EXTRA HANDY UTILITIES:  /////////////////////////////

// size of wings as fraction of length:

#define WINGS	0.10


// axes:

#define X	1
#define Y	2
#define Z	3


// x, y, z, axes:

static float axx[3] = { 1., 0., 0. };
static float ayy[3] = { 0., 1., 0. };
static float azz[3] = { 0., 0., 1. };


void
Arrow( float tail[3], float head[3] )
{
	float u[3], v[3], w[3];		// arrow coordinate system
	float d;			// wing distance
	float x, y, z;			// point to plot
	float mag;			// magnitude of major direction
	float f;			// fabs of magnitude
	int axis;			// which axis is the major


	// set w direction in u-v-w coordinate system:

	w[0] = head[0] - tail[0];
	w[1] = head[1] - tail[1];
	w[2] = head[2] - tail[2];


	// determine major direction:

	axis = X;
	mag = fabs( w[0] );
	if( (f=fabs(w[1]))  > mag )
	{
		axis = Y;
		mag = f;
	}
	if( (f=fabs(w[2]))  > mag )
	{
		axis = Z;
		mag = f;
	}


	// set size of wings and turn w into a Unit vector:

	d = WINGS * Unit( w, w );


	// draw the shaft of the arrow:

	glBegin( GL_LINE_STRIP );
		glVertex3fv( tail );
		glVertex3fv( head );
	glEnd();

	// draw two sets of wings in the non-major directions:

	if( axis != X )
	{
		Cross( w, axx, v );
		(void) Unit( v, v );
		Cross( v, w, u  );
		x = head[0] + d * ( u[0] - w[0] );
		y = head[1] + d * ( u[1] - w[1] );
		z = head[2] + d * ( u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd();
		x = head[0] + d * ( -u[0] - w[0] );
		y = head[1] + d * ( -u[1] - w[1] );
		z = head[2] + d * ( -u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd();
	}


	if( axis != Y )
	{
		Cross( w, ayy, v );
		(void) Unit( v, v );
		Cross( v, w, u  );
		x = head[0] + d * ( u[0] - w[0] );
		y = head[1] + d * ( u[1] - w[1] );
		z = head[2] + d * ( u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd();
		x = head[0] + d * ( -u[0] - w[0] );
		y = head[1] + d * ( -u[1] - w[1] );
		z = head[2] + d * ( -u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd();
	}



	if( axis != Z )
	{
		Cross( w, azz, v );
		(void) Unit( v, v );
		Cross( v, w, u  );
		x = head[0] + d * ( u[0] - w[0] );
		y = head[1] + d * ( u[1] - w[1] );
		z = head[2] + d * ( u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd();
		x = head[0] + d * ( -u[0] - w[0] );
		y = head[1] + d * ( -u[1] - w[1] );
		z = head[2] + d * ( -u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd();
	}
}



float
Dot( float v1[3], float v2[3] )
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}



void
Cross( float v1[3], float v2[3], float vout[3] )
{
	float tmp[3];

	tmp[0] = v1[1]*v2[2] - v2[1]*v1[2];
	tmp[1] = v2[0]*v1[2] - v1[0]*v2[2];
	tmp[2] = v1[0]*v2[1] - v2[0]*v1[1];

	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}



float
Unit( float vin[3], float vout[3] )
{
	float dist, f ;

	dist = vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2];

	if( dist > 0.0 )
	{
		dist = sqrt( dist );
		f = 1. / dist;
		vout[0] = f * vin[0];
		vout[1] = f * vin[1];
		vout[2] = f * vin[2];
	}
	else
	{
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}

	return dist;
}



// the stroke characters 'X' 'Y' 'Z' :

static float xx[] = {
		0.f, 1.f, 0.f, 1.f
	      };

static float xy[] = {
		-.5f, .5f, .5f, -.5f
	      };

static int xorder[] = {
		1, 2, -3, 4
		};


static float yx[] = {
		0.f, 0.f, -.5f, .5f
	      };

static float yy[] = {
		0.f, .6f, 1.f, 1.f
	      };

static int yorder[] = {
		1, 2, 3, -2, 4
		};


static float zx[] = {
		1.f, 0.f, 1.f, 0.f, .25f, .75f
	      };

static float zy[] = {
		.5f, .5f, -.5f, -.5f, 0.f, 0.f
	      };

static int zorder[] = {
		1, 2, 3, 4, -5, 6
		};


// fraction of the length to use as height of the characters:

const float LENFRAC = 0.10f;


// fraction of length to use as start location of the characters:

const float BASEFRAC = 1.10f;


//
//	Draw a set of 3D axes:
//	(length is the axis length in world coordinates)
//

void
Axes( float length )
{
	int i, j;			// counters
	float fact;			// character scale factor
	float base;			// character start location


	glBegin( GL_LINE_STRIP );
		glVertex3f( length, 0., 0. );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., length, 0. );
	glEnd();
	glBegin( GL_LINE_STRIP );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., 0., length );
	glEnd();

	fact = LENFRAC * length;
	base = BASEFRAC * length;

	glBegin( GL_LINE_STRIP );
		for( i = 0; i < 4; i++ )
		{
			j = xorder[i];
			if( j < 0 )
			{
				
				glEnd();
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( base + fact*xx[j], fact*xy[j], 0.0 );
		}
	glEnd();

	glBegin( GL_LINE_STRIP );
		for( i = 0; i < 5; i++ )
		{
			j = yorder[i];
			if( j < 0 )
			{
				
				glEnd();
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( fact*yx[j], base + fact*yy[j], 0.0 );
		}
	glEnd();

	glBegin( GL_LINE_STRIP );
		for( i = 0; i < 6; i++ )
		{
			j = zorder[i];
			if( j < 0 )
			{
				
				glEnd();
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( 0.0, fact*zy[j], base + fact*zx[j] );
		}
	glEnd();

}




//
// routine to convert HSV to RGB
//
// Reference:  Foley, van Dam, Feiner, Hughes,
//		"Computer Graphics Principles and Practices,"
//		Additon-Wesley, 1990, pp592-593.


void
HsvRgb( float hsv[3], float rgb[3] )
{
	float h, s, v;			// hue, sat, value
	float r, g, b;			// red, green, blue
	float i, f, p, q, t;		// interim values


	// guarantee valid input:

	h = hsv[0] / 60.;
	while( h >= 6. )	h -= 6.;
	while( h <  0. ) 	h += 6.;

	s = hsv[1];
	if( s < 0. )
		s = 0.;
	if( s > 1. )
		s = 1.;

	v = hsv[2];
	if( v < 0. )
		v = 0.;
	if( v > 1. )
		v = 1.;


	// if sat==0, then is a gray:

	if( s == 0.0 )
	{
		rgb[0] = rgb[1] = rgb[2] = v;
		return;
	}


	// get an rgb from the hue itself:
	
	i = floor( h );
	f = h - i;
	p = v * ( 1. - s );
	q = v * ( 1. - s*f );
	t = v * ( 1. - ( s * (1.-f) ) );

	switch( (int) i )
	{
		case 0:
			r = v;	g = t;	b = p;
			break;
	
		case 1:
			r = q;	g = v;	b = p;
			break;
	
		case 2:
			r = p;	g = v;	b = t;
			break;
	
		case 3:
			r = p;	g = q;	b = v;
			break;
	
		case 4:
			r = t;	g = p;	b = v;
			break;
	
		case 5:
			r = v;	g = p;	b = q;
			break;
	}


	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
}


void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


  glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (WhichProjection == ORTHO)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);

	glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}


void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

  glEnable (GL_POLYGON_OFFSET_FILL);
  glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	for (i=0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
      glLoadName(i+1);

		Triangle *temp_t=this_poly->tlist[i];

		switch (display_mode) {
		case 0:
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			} else {
				mat_diffuse[0] = 0.6;
				mat_diffuse[1] = 0.8;
				mat_diffuse[2] = 0.7;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i==this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 6:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 10:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
		
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}
	}
}


void Display_Model(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv (GL_VIEWPORT, viewport);
 
	set_view(GL_RENDER, poly);
	set_scene(GL_RENDER, poly);
	display_shape(GL_RENDER, poly);
	glFlush();
	glutSwapBuffers();
	glFinish();
}


void Choose_Object()
{
	int w, h;
	switch(ObjectId){
		case BUNNY:
			strcpy(object_name, "bunny");
			break;

		case FELINE:
			strcpy(object_name, "feline");
			break;

		case DRAGON:
			strcpy(object_name, "dragon");
			break;

		case HAPPY:
			strcpy(object_name, "happy");
			break;

		case TORUS:
			strcpy(object_name, "torus");
			break;

		case SPHERE:
			strcpy(object_name, "sphere");
			break;

	//Glui->add_radiobutton_to_group(group,"lucy");
	//Glui->add_radiobutton_to_group(group,"lion");
	//Glui->add_radiobutton_to_group(group,"heptoroid");
	//Glui->add_radiobutton_to_group(group,"igea");

		//case ICOSAHEDRON:
		//	//strcpy(object_name, "icosahedron");
		//	strcpy(object_name, "Armadillo");
		//	break;

		//case OCTAHEDRON:
		//	strcpy(object_name, "lucy");
		//	break;

		//case HEXAHEDRON:
		//	strcpy(object_name, "heptoroid");
		//	break;

		//case DODECAHEDRON:
		//	strcpy(object_name, "igea");
		//	break;

		//case TETRAHEDRON:
		//	strcpy(object_name, "brain");
		//	break;

	}

	poly->finalize();

    Reset();

	char tmp_str[512];

	sprintf (tmp_str, "../models/%s.ply", object_name);

	FILE *this_file = fopen(tmp_str, "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);

    ////Following codes build the edge information
	clock_t start, finish; //Used to show the time assumed
	start = clock(); //Get the first time

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	finish = clock(); //Get the current time after finished
	double t = (double)(finish - start)/CLOCKS_PER_SEC;

	printf("\n");
	printf("The number of the edges of the object %s is %d \n",object_name, poly->nedges);
	printf("The Euler Characteristics of the object %s is %d \n",object_name, (poly->nverts - poly->nedges + poly->ntris));

	printf("Time to building the edge link is %f seconds\n", t);
				
	Glui->sync_live();
	glutSetWindow( MainWindow );
	glutPostRedisplay();
}


//// For LIC (10/02/2012)
//#include <vector>
//std::vector<double> weights;   // the list for the weights
//std::vector<std::pair<int, int>> passed_pixels;
//const int img_res = 512;
//unsigned char noise_tex[img_res][img_res][4];
//unsigned char vec_img[img_res][img_res][2];   // render the vector field into an image
//unsigned char LIC_tex[img_res][img_res][4];
//void    update_kernel(int size);
//void    gen_test_vecfld ();
//void    gen_noise_tex ();
//void    comp_LIC ();
//void    est_passed_pixels_from (int i, int j, int L);
//void    comp_streamline_forward (int i, int j, int L);
//void    comp_streamline_backward (int i, int j, int L);



void 
update_kernel(int size)
{
}

// normalize vector

void 
normalize_vec (float vec[2])
{
	float len = sqrt((double)(vec[0]*vec[0] + vec[1]*vec[1]));

	if (len < 1.e-6) return;

	vec[0] /= len;
	vec[1] /= len;
}

double 
get_vec_len (float vec[2])
{
	return (sqrt((double)(vec[0]*vec[0]+vec[1]*vec[1])));
}


void 
normalize_vec (double vec[2])
{
	float len = sqrt((double)(vec[0]*vec[0] + vec[1]*vec[1]));

	if (len < 1.e-6) return;

	vec[0] /= len;
	vec[1] /= len;
}

template <class T>
T get_vec3_len(T vec[3])
{
	return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

template <class T>
void normalize_vec3 (T vec[3])
{
	T len = get_vec3_len (vec);

	if (len < 1.e-6) return;

	vec[0] /= len;
	vec[1] /= len;
	vec[2] /= len;
}

double 
get_vec_len (double vec[2])
{
	return (sqrt((double)(vec[0]*vec[0]+vec[1]*vec[1])));
}


double bilinear_interpolate(double a, double b, 
						  double f00, double f01, double f10, double f11)
{
	return (f00*(1-a)*(1-b)+f10*a*(1-b)+f01*(1-a)*b+f11*a*b);
}


// The following routine assume that the vector field has been stored in a 2D square texture vec_img
// 
void get_vec_at_regular_grid_image_based(float x, float y, int xdim, int ydim, float &vx, float &vy)
{
	int lrow = (int) y;
	int lcol = (int) x;

	if (lrow<0 || lrow>=xdim-1 || lcol<0 || lcol>=ydim-1) // reaching the boundary
	{
		vx = vy = 0.;
		return;
	}

	// Do we need to consider the center of the pixel instead, which is (lcol+.5, lrow+.5)
	double a = (x - lcol);
	double b = (y - lrow);

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
	double f00 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol][0]/255.;
	double f01 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol][0]/255.;
	double f10 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol+1][0]/255.;
	double f11 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol+1][0]/255.;
	vx = bilinear_interpolate(a, b, f00, f01, f10, f11);

	f00 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol][1]/255.;
	f01 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol][1]/255.;
	f10 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol+1][1]/255.;
	f11 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol+1][1]/255.;
	vy = bilinear_interpolate(a, b, f00, f01, f10, f11);
}



bool get_vec_at_regular_grid_image_based2(double xy[2], int xdim, int ydim, double vxy[2], bool backward)
{
	double x = xy[0]-.5;
	double y = xy[1]-.5;

	int lrow = (int) y;
	int lcol = (int) x;

	if (lrow<0 || lrow>=xdim-1 || lcol<0 || lcol>=ydim-1) // reaching the boundary
	{
		vxy[0] = vxy[1] = 0.;
		return false;
	}

	// Do we need to consider the center of the pixel instead, which is (lcol+.5, lrow+.5)
	double a = (x - lcol);
	double b = (y - lrow);

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
	//double f00 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol][0]/255.;
	//double f01 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol][0]/255.;
	//double f10 = min_vx + (max_vx - min_vx) * vec_img[lrow][lcol+1][0]/255.;
	//double f11 = min_vx + (max_vx - min_vx) * vec_img[lrow+1][lcol+1][0]/255.;
	double f00 = min_vx + (max_vx - min_vx) * vec_img2[lrow][lcol][0];
	double f01 = min_vx + (max_vx - min_vx) * vec_img2[lrow+1][lcol][0];
	double f10 = min_vx + (max_vx - min_vx) * vec_img2[lrow][lcol+1][0];
	double f11 = min_vx + (max_vx - min_vx) * vec_img2[lrow+1][lcol+1][0];
	vxy[0] = bilinear_interpolate(a, b, f00, f01, f10, f11);


	//f00 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol][1]/255.;
	//f01 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol][1]/255.;
	//f10 = min_vy + (max_vy - min_vy) * vec_img[lrow][lcol+1][1]/255.;
	//f11 = min_vy + (max_vy - min_vy) * vec_img[lrow+1][lcol+1][1]/255.;
	f00 = min_vy + (max_vy - min_vy) * vec_img2[lrow][lcol][1];
	f01 = min_vy + (max_vy - min_vy) * vec_img2[lrow+1][lcol][1];
	f10 = min_vy + (max_vy - min_vy) * vec_img2[lrow][lcol+1][1];
	f11 = min_vy + (max_vy - min_vy) * vec_img2[lrow+1][lcol+1][1];
	vxy[1] = bilinear_interpolate(a, b, f00, f01, f10, f11);

	if (backward) {vxy[0] =-vxy[0]; vxy[1] = -vxy[1];}
	return true;
}

double predict_stepsize = 1.;/* = quadmesh->xinterval;*/


/*we try to implement the RK23 as numerical recipe*/
bool RK23_2d(double pre_p[2], double next_p[2], int xdim, int ydim, double offset[2], double &hstep_loc, double &hnext,
			 double eps, double &eps_did, bool backward, bool deriv(double xy[2], int xdim, int ydim, double vxy[2], bool backward) /*void deriv(double cur_p[2], double vec[2])*/
			 )
{
	double dp0[2], dp1[2];
	double temp[2] = {0.};
	double t_vec[2];
	icVector2 vec1, vec2;
	
	/*compute dp0*/
	if (!deriv(pre_p, xdim, ydim, t_vec, backward)) return false;
	// we need to normalize the vector 
	vec1.set(t_vec);
	normalize(vec1);
	dp0[0] = hstep_loc*/*t_vec*/vec1.entry[0];
	dp0[1] = hstep_loc*/*t_vec*/vec1.entry[1];

	/*compute dp1*/
	temp[0]=pre_p[0]+dp0[0];
	temp[1]=pre_p[1]+dp0[1];
	if (!deriv(temp, xdim, ydim, t_vec, backward)) return false;
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
	if (!deriv(next_p, xdim, ydim, t_vec, backward)) return false;
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
bool RK45_2d(double pre_p[2], double next_p[2], int xdim, int ydim, double offset[2], double &hstep_loc, double &hnext,
			 double eps, double &eps_did, bool backward, bool deriv(double xy[2], int xdim, int ydim, double vxy[2], bool backward) /*void deriv(double cur_p[2], double vec[2])*/
			 )
{
	double dp0[2], dp1[2], dp2[2], dp3[2];
	double temp[2] = {0.};
	double t_vec[2];
	icVector2 vec1, vec2;
	
	/*compute dp0*/
	if (!deriv(pre_p, xdim, ydim, t_vec, backward)) return false;
	// we need to normalize the vector 
	vec1.set(t_vec);
	normalize(vec1);
	dp0[0] = hstep_loc*/*t_vec*/vec1.entry[0];
	dp0[1] = hstep_loc*/*t_vec*/vec1.entry[1];

	/*compute dp1*/
	temp[0]=pre_p[0]+dp0[0]/2;
	temp[1]=pre_p[1]+dp0[1]/2;
	if (!deriv(temp, xdim, ydim, t_vec, backward)) return false;
	vec2.set(t_vec);
	normalize(vec2);
	dp1[0] = hstep_loc*/*t_vec*/vec2.entry[0];
	dp1[1] = hstep_loc*/*t_vec*/vec2.entry[1];


	/*compute dp2*/
	temp[0]=pre_p[0]+dp1[0]/2;
	temp[1]=pre_p[1]+dp1[1]/2;
	if (!deriv(temp, xdim, ydim, t_vec, backward)) return false;
	vec2.set(t_vec);
	normalize(vec2);
	dp2[0] = hstep_loc*vec2.entry[0];
	dp2[1] = hstep_loc*vec2.entry[1];

	/*compute dp3*/
	temp[0]=pre_p[0]+dp2[0];
	temp[1]=pre_p[1]+dp2[1];
	if (!deriv(temp, xdim, ydim, t_vec, backward)) return false;
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
	if (!deriv(next_p, xdim, ydim, t_vec, backward)) return false;
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

/*
Use RK23 to do tensor line computation
*/
bool get_nextpt_RK23_quad(double first[2], double second[2], double offset[2], int xdim, int ydim, bool type)
{

	double t_vec[2] = {0.};
	if (!get_vec_at_regular_grid_image_based2(first, xdim, ydim, t_vec, type)) return false;
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
		if (!RK23_2d(first, second, xdim, ydim, offset, hstep_loc, hnext, eps, eps_did, type, get_vec_at_regular_grid_image_based2)) return false;
		predict_stepsize = hnext;
		if(eps_did<eps)
			break;
	}
	return true;
}


/*
Use RK45 to do tensor line computation
*/
bool get_nextpt_RK45_quad(double first[2], double second[2], double offset[2], int xdim, int ydim, bool type)
{

	double t_vec[2] = {0.};
	if (!get_vec_at_regular_grid_image_based2(first, xdim, ydim, t_vec, type)) return false;
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
		if (!RK45_2d(first, second, xdim, ydim, offset, hstep_loc, hnext, eps, eps_did, type, get_vec_at_regular_grid_image_based2)) return false;
		predict_stepsize = hnext;
		if(eps_did<eps)
			break;
	}
	return true;
}


void
gen_test_vecfld()
{
	// consider a center defined at the pixel (img_res/2, img_res/2)


	int i, j;
	float interval = 1./(img_res-1);
	float center[2] = {img_res/2., img_res/2.};
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			//vec_2d[i][j][0] = -(-0.5+j*interval);
			//vec_2d[i][j][1] = (-0.5+i*interval);
			vec_2d[i][j][0] = i*interval-0.5;
			vec_2d[i][j][1] = j*interval-0.5;
			//vec_2d[i][j][0] = -j*interval;
			//vec_2d[i][j][1] = i*interval;
			//vec_2d[i][j][0] = -(j-center[1]);
			//vec_2d[i][j][1] = (i-center[0]);
			normalize_vec(vec_2d[i][j]);
		}
	}

	// find out the max_vx, min_vx, max_vy, and min_vy
	max_vy = max_vx = -1.e8;
	min_vy = min_vx = 1.e8;
	for (i=0; i<img_res; i++)
		for (j=0; j<img_res; j++)
		{
			if (vec_2d[i][j][0] > max_vx) max_vx = vec_2d[i][j][0];
			if (vec_2d[i][j][0] < min_vx) min_vx = vec_2d[i][j][0];
			if (vec_2d[i][j][1] > max_vy) max_vy = vec_2d[i][j][1];
			if (vec_2d[i][j][1] < min_vy) min_vy = vec_2d[i][j][1];
		}

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			vec_img[i][j][0] = (unsigned char) 255*(vec_2d[i][j][0] - min_vx)/(max_vx - min_vx);
			vec_img[i][j][1] = (unsigned char) 255*(vec_2d[i][j][1] - min_vy)/(max_vy - min_vy);
			vec_img[i][j][2] = 120;
		}
	}
}


// generate random noise
void 
gen_noise_tex ()
{
    for (int x = 0; x < img_res; x++)
    for (int y = 0; y < img_res; y++)
    {
        noise_tex[x][y][0] = 
		noise_tex[x][y][1] = 
		noise_tex[x][y][2] = (unsigned char) 255*(rand() % 32768) / 32768.0;
    }
}


void    
replace_noise_tex (unsigned char LIC_tex[img_res][img_res][4])
{
	unsigned char max_v, min_v;
	max_v = 0, min_v = 255;
    for (int x = 0; x < img_res; x++)
    for (int y = 0; y < img_res; y++)
	{
		if (LIC_tex[x][y][0] > max_v) max_v = LIC_tex[x][y][0];
		if (LIC_tex[x][y][0] < min_v) min_v = LIC_tex[x][y][0];
	}

    for (int x = 0; x < img_res; x++)
    for (int y = 0; y < img_res; y++)
    {
		if (LIC_tex[x][y][0] < 100)
		{
        noise_tex[x][y][0] = 
		noise_tex[x][y][1] = 
		noise_tex[x][y][2] = unsigned char(0.9*LIC_tex[x][y][0])%255;
		}

		//else if (LIC_tex[x][y][0] > 200)
		//{
  //      noise_tex[x][y][0] = 
		//noise_tex[x][y][1] = 
		//noise_tex[x][y][2] = 255;
		//}

		else
		{
        noise_tex[x][y][0] = 
		noise_tex[x][y][1] = 
		noise_tex[x][y][2] = unsigned char(1.1*LIC_tex[x][y][0])%255;
		}
    }
}



bool 
RK2_integrate_multi_steps(double trace_gp[2], int N, bool forward_backward)
{
	double vec[2], vec2[2];
	double offset[2]={0.};

	double step_size = 1./N;

	for (int i=0; i<N; i++)
	{
		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec, forward_backward);
		
		if (get_vec_len(vec) < 1.e-6) return false;

		normalize_vec(vec);
		double tx = trace_gp[0]+step_size/2.*vec[0];
		double ty = trace_gp[1]+step_size/2.*vec[1];

		double tp[2]={tx,ty};
		double vec2[2];
		get_vec_at_regular_grid_image_based2(tp, img_res, img_res, vec2, forward_backward);
		normalize_vec(vec2);
		vec[0] = (vec[0] + vec2[0])/2.;
		vec[1] = (vec[1] + vec2[1])/2.;

		normalize_vec(vec);

		trace_gp[0] =trace_gp[0]+step_size*vec[0];
		trace_gp[1] =trace_gp[1]+step_size*vec[1];
	}

	return true;
}


bool 
Euler_integrate_multi_steps(double trace_gp[2], int N, bool forward_backward)
{
	double vec[2], vec2[2];
	double offset[2]={0.};

	double step_size = 1./N;

	for (int i=0; i<N; i++)
	{
		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec, forward_backward);
		
		if (get_vec_len(vec) < 1.e-6) return false;

		normalize_vec(vec);

		trace_gp[0] =trace_gp[0]+step_size*vec[0];
		trace_gp[1] =trace_gp[1]+step_size*vec[1];
	}

	return true;
}


bool 
RK4_integrate_multi_steps(double trace_gp[2], int N, bool forward_backward)
{
	double vec[2], vec2[2], vec3[2], vec4[2];
	double offset[2]={0.};
	double tmp[2];

	double step_size = 1./N;

	for (int i=0; i<N; i++)
	{
		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec, forward_backward);
		
		if (get_vec_len(vec) < 1.e-6) return false;

		normalize_vec(vec);

		tmp[0] =trace_gp[0]+step_size/2*vec[0];
		tmp[1] =trace_gp[1]+step_size/2*vec[1];
		
		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec2, forward_backward);
		normalize_vec(vec2);

		tmp[0] =trace_gp[0]+step_size/2*vec2[0];
		tmp[1] =trace_gp[1]+step_size/2*vec2[1];

		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec3, forward_backward);
		normalize_vec(vec3);

		tmp[0] =trace_gp[0]+step_size*vec3[0];
		tmp[1] =trace_gp[1]+step_size*vec3[1];
		
		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec4, forward_backward);
		normalize_vec(vec4);

		vec[0] = (vec[0]+2*vec2[0]+2*vec3[0]+vec4[0])/6.;
		vec[1] = (vec[1]+2*vec2[1]+2*vec3[1]+vec4[1])/6.;
		normalize_vec(vec);
		
		trace_gp[0] =trace_gp[0]+step_size*vec[0];
		trace_gp[1] =trace_gp[1]+step_size*vec[1];
	}

	return true;
}

void    
est_passed_pixels_from (int i, int j, int L)
{
	//passed_pixels.clear();
	std::vector<std::pair<int, int>> passed_pixels;
	std::pair<int, int> start_pixel;
	start_pixel.first = i;
	start_pixel.second = j;
	passed_pixels.push_back (start_pixel);

	comp_streamline_forward (i, j, L, passed_pixels);
	comp_streamline_backward (i, j, L, passed_pixels);

	compose_color(LIC_tex[i][j], passed_pixels);
	passed_pixels.clear();
}

// compute forward direction 
void    
comp_streamline_forward (int i, int j, int L, std::vector<std::pair<int, int>> &passed_pixels)
{
	//float x = i+.5;
	//float y = j+.5;
	float x = j+.5;
	float y = i+.5;

	double vec[2];
	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (int step=0; step < L; step++)
	{
	/*	// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec, false);

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 1.e-6) return;

		normalize_vec(vec);
		//x = x+vec[0];
		//y = y+vec[1];
		//get_nextpt_RK23_quad(trace_gp, trace_gp, offset, img_res, img_res, false);
		double tx = trace_gp[0]+0.125*vec[0];
		double ty = trace_gp[1]+0.125*vec[1];

		double tp[2]={tx,ty};
		double vec2[2];
		get_vec_at_regular_grid_image_based2(tp, img_res, img_res, vec2, false);
		normalize_vec(vec2);
		vec[0] += vec2[0];
		vec[1] += vec2[1];

		normalize_vec(vec);

		trace_gp[0] =trace_gp[0]+0.25*vec[0];
		trace_gp[1] =trace_gp[1]+0.25*vec[1];
		*/

		if (!RK2_integrate_multi_steps(trace_gp, 1, false)) return;

		//i = (int) x;
		//j = (int) y;
		//j = (int) x;
		//i = (int) y;
		j = (int) trace_gp[0];
		i = (int) trace_gp[1];

		std::pair<int, int> pixel;
		pixel.first = i;
		pixel.second = j;

		passed_pixels.push_back(pixel);

		if (i<0 || i>=img_res || j<0 || j>=img_res)
			return;
	}
}

void    
comp_streamline_backward (int i, int j, int L, std::vector<std::pair<int, int>> &passed_pixels)
{
	//float x = i+.5;
	//float y = j+.5;
	float x = j+.5;
	float y = i+.5;

	double vec[2];
	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

	for (int step=0; step < L; step++)
	{
	/*	// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
		get_vec_at_regular_grid_image_based2(trace_gp, img_res, img_res, vec, true);

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 1.e-6) return;

		normalize_vec(vec);
		//x = x-vec[0];
		//y = y-vec[1];
		//get_nextpt_RK23_quad(trace_gp, trace_gp, offset, img_res, img_res, false);
		
		double tx = trace_gp[0]+0.125*vec[0];
		double ty = trace_gp[1]+0.125*vec[1];

		double tp[2]={tx,ty};
		double vec2[2];
		get_vec_at_regular_grid_image_based2(tp, img_res, img_res, vec2, true);
		normalize_vec(vec2);
		vec[0] += vec2[0];
		vec[1] += vec2[1];

		normalize_vec(vec);

		trace_gp[0] =trace_gp[0]+0.25*vec[0];
		trace_gp[1] =trace_gp[1]+0.25*vec[1];
		*/

		if (!RK2_integrate_multi_steps(trace_gp, 1, true)) return;

		//i = (int) x;
		//j = (int) y;
		//j = (int) x;
		//i = (int) y;
		j = (int) trace_gp[0];
		i = (int) trace_gp[1];

		std::pair<int, int> pixel;
		pixel.first = i;
		pixel.second = j;

		passed_pixels.push_back(pixel);

		if (i<0 || i>=img_res || j<0 || j>=img_res)
			return;
	}
}


void
compose_color (unsigned char color_val[3])
{
	// currently, let us simply average the colors of the obtained pixels
	int i;

	int sum_color[3] = {0};
	for (i=0; i<passed_pixels.size(); i++)
	{
		int row = passed_pixels[i].first;
		int col = passed_pixels[i].second;
		sum_color[0] += noise_tex[row][col][0];
		sum_color[1] += noise_tex[row][col][1];
		sum_color[2] += noise_tex[row][col][2];
	}

	color_val[0] = sum_color[0]/passed_pixels.size();
	color_val[1] = sum_color[1]/passed_pixels.size();
	color_val[2] = sum_color[2]/passed_pixels.size();
}

void
compose_color (unsigned char color_val[3], std::vector<std::pair<int, int>> passed_pixels)
{
	// currently, let us simply average the colors of the obtained pixels
	int i;

	int sum_color[3] = {0};
	for (i=0; i<passed_pixels.size(); i++)
	{
		int row = passed_pixels[i].first;
		int col = passed_pixels[i].second;
		sum_color[0] += noise_tex[row][col][0];
		sum_color[1] += noise_tex[row][col][1];
		sum_color[2] += noise_tex[row][col][2];
	}

	color_val[0] = sum_color[0]/passed_pixels.size();
	color_val[1] = sum_color[1]/passed_pixels.size();
	color_val[2] = sum_color[2]/passed_pixels.size();
}

void 
comp_LIC(unsigned char LIC_tex[img_res][img_res][4])
{
	int L=10;
#pragma omp parallel for 
	for (int i = 0; i<img_res; i++)
		for (int j=0; j<img_res; j++)
		{
			if (fixedPt_pixels[i][j])
			{
				LIC_tex[i][j][0] = LIC_tex[i][j][1] = LIC_tex[i][j][2] = 200;
				LIC_tex[i][j][3] = 255;
				continue;
			}
			else
			{
			est_passed_pixels_from(i, j, L);
			//compose_color(LIC_tex[i][j]);
			//LIC_tex[i][j][3] = 125;    // for combining with second LIC texture
			LIC_tex[i][j][3] = 255;
			}
		}
}


void draw_arrow_head(float head[2], float direct[2])
{
    glPushMatrix();
    glTranslatef(head[0], head[1], 0);
    glRotatef(atan2(direct[1], direct[0])*360/(2*M_PI), 0, 0, 1);
    ////glScalef(0.03, 0.03, 1);
    glScalef(0.2, 0.2, 1);

    glBegin(GL_TRIANGLES);
        glVertex2f(0, 0);
        glVertex2f(-0.35, 0.12);
        glVertex2f(-0.35, -0.12);
    glEnd();

    glPopMatrix();
}


void 
jitter_vectorFld (Polyhedron *this_poly, double prob, double perc)
{
	// Should we reset the rand seed??

	int i;
	for (i=0; i<this_poly->nverts; i++)
	{
		double cur_prob = (double)rand()/(double)RAND_MAX;

		if (rand()<1-prob) continue;

		double scale = (double)rand()/(double)RAND_MAX*perc;
		double theta = 2*M_PI*(double)rand()/(double)RAND_MAX;

		double vec_len = sqrt(this_poly->vlist[i]->vec.entry[0]*this_poly->vlist[i]->vec.entry[0] 
			+ this_poly->vlist[i]->vec.entry[1]*this_poly->vlist[i]->vec.entry[1]);

		icVector2 offset(vec_len*scale*cos(theta), vec_len*scale*sin(theta));

		this_poly->vlist[i]->vx = this_poly->vlist[i]->vec.entry[0] + offset.entry[0];
		this_poly->vlist[i]->vy = this_poly->vlist[i]->vec.entry[1] + offset.entry[1];
	}
}

void
store_MS_Label(Img_MorseDecomp &ismd)
{
	int i, j;

	MS_Label label1;

	for (i=0; i<img_res; i++)
		for (j=0; j<img_res; j++)
			label1.labels[i][j]=0;


	for (i=0; i<ismd.scclist.size(); i++)
	{
		if (!ismd.scclist[i].valid) continue;

	    //get_random_Color(i, rgb);
		unsigned char type = ismd.scclist[i].classification;

		for (j=0; j<ismd.scclist[i].pixels.size(); j++)
		{
			int row = ismd.scclist[i].pixels[j].first;
			int col = ismd.scclist[i].pixels[j].second;

			label1.labels[row][col] = type;
		}
	}

	MS_Labels.push_back(label1);
}


void 
count_MS_hits()
{
	int iters, i, j;

	for (i=0; i<img_res; i++)
		for (j=0; j<img_res; j++)
		{
			MS_hits[i][j][0] = MS_hits[i][j][1] = 0;
		}

	for (i=0; i<img_res; i++)
		for (j=0; j<img_res; j++)
		{
			int type_count_so = 0, type_count_si = 0, type_count_sa = 0;

			for (iters=0; iters<MS_Labels.size(); iters++)
			{
				MS_Label &cur_labels = MS_Labels[iters];

				if (cur_labels.labels[i][j]==0)  //not within a Morse set
					continue;           

				MS_hits[i][j][0] = MS_hits[i][j][0]+1;

				if (cur_labels.labels[i][j]==1)  //source type
					type_count_so++;
				else if (cur_labels.labels[i][j]==2)  //sink type
					type_count_si++;
				else if (cur_labels.labels[i][j]==3)  //saddle type
					type_count_sa++;
			}

			if (type_count_so>type_count_si && type_count_so>type_count_sa)
				MS_hits[i][j][1] = 1;
			else if (type_count_si>type_count_so && type_count_si>type_count_sa)
				MS_hits[i][j][1] = 2;
			else if (type_count_sa>type_count_si && type_count_sa>type_count_so)
				MS_hits[i][j][1] = 3;
			else
				MS_hits[i][j][1] = 4;  //uncertain type
		}

		MS_Labels.clear();
}


void
blend_color_uncertainty(int Ntrials)
{
	int i, j;
	float hsv[3], rgb[3];
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			if (MS_hits[i][j][0] == 0)
			{
				out_LIC[i][j][0] = LIC_tex[i][j][0];
				out_LIC[i][j][1] = LIC_tex[i][j][1];
				out_LIC[i][j][2] = LIC_tex[i][j][2];
			}

			else
			{
				if (MS_hits[i][j][1] == 1)  // Green for sources
					hsv[0] = 120.;   
				else if (MS_hits[i][j][1] == 2)  // Red for sinks
					hsv[0] = 0.;
				else if (MS_hits[i][j][1] == 3)  // Blue for saddles
					hsv[0] = 240.;
				else
					hsv[0] = 330;

				hsv[1] = (float) MS_hits[i][j][0]/(float)Ntrials;
				//hsv[1] = 1.;
				hsv[2] = hsv[1]; //1.;
				HsvRgb(hsv, rgb);
				out_LIC[i][j][0] = 0.8*rgb[0]*255+0.2*LIC_tex[i][j][0];
				out_LIC[i][j][1] = 0.8*rgb[1]*255+0.2*LIC_tex[i][j][1];
				out_LIC[i][j][2] = 0.8*rgb[2]*255+0.2*LIC_tex[i][j][2];

				LIC_tex[i][j][0] = out_LIC[i][j][0];
				LIC_tex[i][j][1] = out_LIC[i][j][1];
				LIC_tex[i][j][2] = out_LIC[i][j][2];
			}
		}
	}
}

void    
render_vec_img( Polyhedron *this_poly)
{
	//glViewport(0, 0, (GLsizei) 512, (GLsizei) 512);
	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClearColor( 0., 0., 0., 1. );
	glClear(GL_COLOR_BUFFER_BIT);

	glRenderMode(GL_SMOOTH);

	glDrawBuffer(GL_BACK);
	int i, j;

	// first search the max_vx, min_vx, max_vy, min_vy
	max_vy = max_vx = -1.e8;
	min_vy = min_vx = 1.e8;
	for (j=0; j<this_poly->nverts; j++)
	{
		if (this_poly->vlist[j]->vx > max_vx) max_vx = this_poly->vlist[j]->vx;
		if (this_poly->vlist[j]->vx < min_vx) min_vx = this_poly->vlist[j]->vx;
		if (this_poly->vlist[j]->vy > max_vy) max_vy = this_poly->vlist[j]->vy;
		if (this_poly->vlist[j]->vy < min_vy) min_vy = this_poly->vlist[j]->vy;
	}

	// Can we zoom in to certain region
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	//glTranslatef(0.75, 0.15, 0);  // This combination works well for the multicycle data set zoom in
	//glScalef(3, 3, 3);
	//glTranslatef(0.5, 0.5, 0);  // This combination works well for the multicycle data set whole domain
	//glScalef(1.7, 1.7, 1.7);

	////glTranslatef(0.25, 0.13, 0);  // This combination works well for the satlantic data set
	////glScalef(2.2, 2.2, 2.2);

	//glTranslatef(-0.5, -0.5, 0);

	//glTranslatef(0.5, 0.2, 0);
	//glScalef(1.2, 1.2, 1.2);         // for 20deg5_0mps_highpass
	////////glTranslatef(-0., -0.25, 0);   // for 20deg5_0mps_highpass
	//glTranslatef(-0.5, -0.5, 0);

	////////////////////////////////////////////////////

	for (i=0; i<this_poly->ntris; i++) {
		Triangle *temp_t=this_poly->tlist[i];
		float rgb[3];
		rgb[2] = 0.0;
		//glBegin(GL_POLYGON);
		glBegin(GL_TRIANGLES);
		for (j=0; j<3; j++)
		{
			Vertex *v = temp_t->verts[j];
			//rgb[1] =  (v->vx - min_vx)/(max_vx - min_vx);
			//rgb[0] =  (v->vy - min_vy)/(max_vy - min_vy);
			rgb[0] =  (v->vx - min_vx)/(max_vx - min_vx);
			rgb[1] =  (v->vy - min_vy)/(max_vy - min_vy);

			if (length(v->vec) < 1.e-8) rgb[0] = rgb[1] = 0.;
			glColor3fv (rgb);
			glVertex2f (v->x, v->y);
		}
		glEnd();
	}
	glPopMatrix();

	// initialize the vec_img and vec_img2
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			vec_img[i][j][0] = vec_img[i][j][1] = vec_img[i][j][2] = 0;
			vec_img2[i][j][0] = vec_img2[i][j][1] = vec_img2[i][j][2] = 0.;
		}
	}

	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, vec_img);

	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_FLOAT, vec_img2);

	// compute the vector field magnitude value now
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			float vec[2] = {vec_img2[i][j][0], vec_img2[i][j][1]};
			vec_mag[i][j] = get_vec_len(vec);

			fixedPt_pixels[i][j] = false;
		}
	}

	// initialize the vec_img and vec_img2
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			//vec_img[i][j][0] = vec_img[i][j][1] = vec_img[i][j][2] = 0;
			//vec_img2[i][j][0] = vec_img2[i][j][1] = vec_img2[i][j][2] = 0.;
			//double vx = vec_img[i][j][0]/255.*(max_vx-min_vx) + min_vx;
			//double vy = vec_img[i][j][1]/255.*(max_vy-min_vy) + min_vy;
			icVector2 vec(vec_img2[i][j][0], vec_img2[i][j][1]);
			double vf_mag = length(vec);

			if (vf_mag < 1.e-6)
			{
				vec_img[i][j][0] = vec_img[i][j][1] = vec_img[i][j][2] = 0;
				vec_img2[i][j][0] = vec_img2[i][j][1] = vec_img2[i][j][2] = 0.;
				fixedPt_pixels[i][j] = true;
			}
		}
	}
	write_ppm("render_vec_img.ppm", (unsigned char*)vec_img, img_res, img_res);

	// also, compute the winding number for each pixel to locate fixed points explicitly
	                    

	//for (i=1; i<img_res-1; i++)
	//{
	//	for (j=1; j<img_res-1; j++)
	//	{
	//	}
	//}

extern std::vector<Singularity> singularity_list;
	// we project the detected fixed points on the image space
	int row, col, row1, col1; 
	double dx = 1./(double) (img_res-1);
	for (i=0; i<singularity_list.size(); i++)
	{
		row = (int)(singularity_list[i].gpos.entry[1]*img_res);
		col = (int)(singularity_list[i].gpos.entry[0]*img_res);

		//vec_img[row][col][0] = vec_img[row][col][1] = vec_img[row][col][2] = 0;
		//vec_img2[row][col][0] = vec_img2[row][col][1] = vec_img2[row][col][2] = 0.;

		// The range of the fixed point also affects the result a lot!!!
		//for (int ki=-5; ki<6; ki++)
		//{
		//	for (int kj=-5; kj<6; kj++)
		//	{
		//		row1 = row + ki;
		//		if (row1 < 0) continue;

		//		col1 = col + kj;
		//		if (col1 < 0) continue;

		//		//vec_img[row1][col1][0] = vec_img[row1][col1][1] = vec_img[row1][col1][2] = 0;
		//		//vec_img2[row1][col1][0] = vec_img2[row1][col1][1] = vec_img2[row1][col1][2] = 0.;

		//		fixedPt_pixels[row1][col1] = true;

		//	}
		//}
		//for (int ki=-2; ki<3; ki++)
		//{
		//	for (int kj=-2; kj<3; kj++)
		for (int ki=-0; ki<1; ki++)
		{
			for (int kj=-0; kj<1; kj++)
			{
				row1 = row + ki;
				if (row1 < 0) continue;

				col1 = col + kj;
				if (col1 < 0) continue;

				//vec_img[row1][col1][0] = vec_img[row1][col1][1] = vec_img[row1][col1][2] = 0;
				//vec_img2[row1][col1][0] = vec_img2[row1][col1][1] = vec_img2[row1][col1][2] = 0.;

				fixedPt_pixels[row1][col1] = true;

			}
		}
	}
}


void 
render_tri_ID_img(Polyhedron *this_poly)
{
	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClearColor( 0., 0., 0., 1. );
	glClear(GL_COLOR_BUFFER_BIT);

	glRenderMode(GL_SMOOTH);

	glDrawBuffer(GL_BACK);
	int i, j;


	// We can zoom in to certain region
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	//glLoadIdentity();
	//glTranslatef(0.75, 0.15, 0);  // This combination works well for the multicycle data set zoom in
	//glScalef(3, 3, 3);
	//glTranslatef(0.5, 0.5, 0);  // This combination works well for the multicycle data set whole domain
	//glScalef(1.7, 1.7, 1.7);

	////glTranslatef(0.25, 0.13, 0);  // This combination works well for the satlantic data set
	////glScalef(2.2, 2.2, 2.2);

	//glTranslatef(-0.5, -0.5, 0);

	//glTranslatef(0.5, 0.2, 0);
	//glScalef(1.2, 1.2, 1.2);         // for 20deg5_0mps_highpass
	////glTranslatef(-0., -0.25, 0);   // for 20deg5_0mps_highpass
	//glTranslatef(-0.5, -0.5, 0);

	////////////////////////////////////////////////////
	unsigned char rgb[3];

	for (i=0; i<this_poly->ntris; i++) {
		Triangle *temp_t=this_poly->tlist[i];
		rgb[2] = i/(256*256);
		rgb[0] = i%256;
		rgb[1] = (i-rgb[2]*256*256 - rgb[1])/256;

		glBegin(GL_TRIANGLES);
		for (j=0; j<3; j++)
		{
			Vertex *v = temp_t->verts[j];
			glColor3f ((float)rgb[0]/255., (float)rgb[1]/255., (float)rgb[2]/255.);
			glVertex2f (v->x, v->y);
		}
		glEnd();
	}
	glPopMatrix();
	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, tri_ID_img);
	
	write_ppm("tri_ID_img.ppm", (unsigned char*)tri_ID_img, img_res, img_res);
}

void 
assign_tri_pixel_labels(Img_MorseDecomp &ismd)
{
	int i, j;
	for (i=0; i<img_res; i++)
		for (j=0; j<img_res; j++)
		{
			pixel_labels[i][j]=0;
			pixel_tri_labels[i][j]=0;
		}

	// get the labels for the pixels based on the ISMD result
	for (i=0; i<ismd.scclist.size(); i++)
	{
		if (!ismd.scclist[i].valid) continue;

	    //get_random_Color(i, rgb);
		unsigned char type = ismd.scclist[i].classification;

		for (j=0; j<ismd.scclist[i].pixels.size(); j++)
		{
			int row = ismd.scclist[i].pixels[j].first;
			int col = ismd.scclist[i].pixels[j].second;

			pixel_labels[row][col] = type;
		}
	}

	// get the labels for the pixels based on the tranditional Morse decomposition
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			int tri_ID = (int)tri_ID_img[i][j][2]*256*256 + (int)tri_ID_img[i][j][1]*256 + (int)tri_ID_img[i][j][0];
			pixel_tri_labels[i][j] = tri_labels[tri_ID];
		}
	}
}


void 
compute_n_I_and_n_II()
{
	n_I = n_II = 0;

	int i, j;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			combined_n_I_n_II_img[i][j][0] = n_I_img[i][j][0] = n_II_img[i][j][0] = gray_MS_img[i][j][0];
			combined_n_I_n_II_img[i][j][1] = n_I_img[i][j][1] = n_II_img[i][j][1] = gray_MS_img[i][j][1];
			combined_n_I_n_II_img[i][j][2] = n_I_img[i][j][2] = n_II_img[i][j][2] = gray_MS_img[i][j][2];
		}
	}

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			if (pixel_labels[i][j]>0 && pixel_tri_labels[i][j] == 0)
			{
				n_I ++;
				combined_n_I_n_II_img[i][j][0] = n_I_img[i][j][0] = 0.6*255+0.4*n_I_img[i][j][0];
				combined_n_I_n_II_img[i][j][1] = n_I_img[i][j][1] = 0.6*255+0.4*n_I_img[i][j][1];
			}

			else if (pixel_labels[i][j]==0 && pixel_tri_labels[i][j] > 0)
			{
				n_II ++;
				combined_n_I_n_II_img[i][j][1] = n_II_img[i][j][1] = 0.6*255+0.4*n_II_img[i][j][1];
				combined_n_I_n_II_img[i][j][2] = n_II_img[i][j][2] = 0.6*255+0.4*n_II_img[i][j][2];
			}
		}
	}
}


void
get_gray_MS_img(Img_MorseDecomp &ismd)
{
	int i, j;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			gray_MS_img[i][j][0] = LIC_tex[i][j][0]; 
			gray_MS_img[i][j][1] = LIC_tex[i][j][1]; 
			gray_MS_img[i][j][2] = LIC_tex[i][j][2]; 
		}
	}
	for (i=0; i<ismd.scclist.size(); i++)
	{
		if (!ismd.scclist[i].valid) continue;
		for (j=0; j<ismd.scclist[i].pixels.size(); j++)
		{
			int row = ismd.scclist[i].pixels[j].first;
			int col = ismd.scclist[i].pixels[j].second;
			
			gray_MS_img[row][col][0] = 0.7*100 + 0.3 * LIC_tex[row][col][0]; 
			gray_MS_img[row][col][1] = 0.7*100 + 0.3 * LIC_tex[row][col][1]; 
			gray_MS_img[row][col][2] = 0.7*100 + 0.3 * LIC_tex[row][col][2]; 

		}
	}
}

////////////////////////////////////////////////////////////////////////////////////

void    
render_LIC_img()
{
	//glClearColor( BACKCOLOR[0], BACKCOLOR[1], BACKCOLOR[2], BACKCOLOR[3] );
	glClearColor( 1., 1., 1., 1. );
	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	//gluOrtho2D(-1., 1, -1., 1);
	gluOrtho2D(0., 1, 0., 1);
	glClear(GL_COLOR_BUFFER_BIT);

	                    
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, 
					GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);

	glEnable (GL_BLEND);
	// Test noise texture
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, GL_UNSIGNED_BYTE, noise_tex);

	// Test vector field image
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, GL_UNSIGNED_BYTE, vec_img);

	// Render LIC image
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, GL_UNSIGNED_BYTE, LIC_tex);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, GL_UNSIGNED_BYTE, LIC_tex_comb);
	
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_res, img_res, 0,
	//	GL_RGBA, GL_UNSIGNED_BYTE, LIC_tex);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_res, img_res, 0,
	//	GL_RGBA, GL_UNSIGNED_BYTE, LIC_tex2);

	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, /*GL_FLOAT*/GL_UNSIGNED_BYTE, forward_counter_color);

	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, /*GL_FLOAT*/GL_UNSIGNED_BYTE, comb_counter_color); //comb_counter_color
	
	//if (RotSumDiffOn == 0)
	//{
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, /*GL_FLOAT*/GL_UNSIGNED_BYTE, rot_sum_color); //comb_counter_color 
	//}

	//else
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
	//	GL_RGB, /*GL_FLOAT*/GL_UNSIGNED_BYTE, total_rotation_diff_color); //comb_counter_color 

	//glBegin(GL_QUAD_STRIP);
	//	glTexCoord2f(0.0,  0.0);  glVertex2f(0, 0);
	//	glTexCoord2f(0.0,  1.0);  glVertex2f(0, 1.0);
	//	glTexCoord2f(1.0,  0.0);  glVertex2f(1.0, 0.);
	//	glTexCoord2f(1.0,  1.0);  glVertex2f(1.0, 1.0);
	//glEnd();

	///////
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_res, img_res, 0,
	//	GL_RGBA, GL_UNSIGNED_BYTE, LIC_tex2);

	if (Show_OnlyLIC == 1)
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img_res, img_res, 0,
			GL_RGBA, GL_UNSIGNED_BYTE, LIC_tex);
	}

	else if (Show_ROTSUM_or_DIFF ==0 )
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
				GL_RGB, /*GL_FLOAT*/GL_UNSIGNED_BYTE, rot_sum_color);
	}
	else
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_res, img_res, 0,
			GL_RGB, /*GL_FLOAT*/GL_UNSIGNED_BYTE, total_rotation_diff_color);	
	}
	glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0, 0);
		glTexCoord2f(0.0,  1.0);  glVertex2f(0, 1.0);
		glTexCoord2f(1.0,  0.0);  glVertex2f(1.0, 0.);
		glTexCoord2f(1.0,  1.0);  glVertex2f(1.0, 1.0);
	glEnd();	

	
	glDisable(GL_TEXTURE_2D);

	if (ArrowsOn==1)
	{

		//glEnable(GL_COLOR_MATERIAL);
		//// Draw the arrow plot
		//glColor3f (1, 1, 0);
		//for (int i=0; i<poly->nverts; i++) {
		//	Vertex *v=poly->vlist[i];
		//	float head[2] = {v->x+0.05*v->vx, v->y+0.05*v->vy};
		//	float direct[2] = {0.05*v->vx, 0.05*v->vy};
		//	glBegin(GL_LINES);
		//		glVertex2f (v->x, v->y);
		//		glVertex2fv (head);
		//	glEnd();

		//	draw_arrow_head (head, direct);
		//}

		glEnable(GL_COLOR_MATERIAL);
		glColor3f (1, 1, 0);

		for (int i=0; i<img_res; i++)
		{
			for (int j=0; j<img_res; j++)
			{
				if (i%16!=0 || j%16 != 0) continue;
				float x = (float)j/(img_res-1), y = (float)i/(img_res-1); 
				float head[2];
				float direct[2] = {total_rotation_diff[i][j][0], total_rotation_diff[i][j][1]};
				normalize_vec(direct);
				head[0] = x + arrow_scale*direct[0];
				head[1] = y + arrow_scale*direct[1];


				glBegin(GL_LINES);
				glVertex2f (x, y);
				glVertex2fv (head);
				glEnd();

				draw_arrow_head (head, direct);
			}

			// Draw an arrow at the selected pixel
			float x = (float)sel_pos_j/(img_res-1), y = (float)sel_pos_i/(img_res-1); 
			float head[2];
			float direct[2] = {total_rotation_diff[sel_pos_i][sel_pos_j][0], total_rotation_diff[sel_pos_i][sel_pos_j][1]};
			normalize_vec(direct);
			head[0] = x + arrow_scale*direct[0];
			head[1] = y + arrow_scale*direct[1];


			glBegin(GL_LINES);
			glVertex2f (x, y);
			glVertex2fv (head);
			glEnd();

			draw_arrow_head (head, direct);
		}
	}

	//// Display the streamline
	//glColor3f (0, 0, 0); // the color should correspond to the color of the 2D plot
	glDisable(GL_TEXTURE_2D);

	if (ShowSelStreamline)
	{
	glEnable( GL_LINE_SMOOTH );
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	//glEnable(GL_BLEND); 
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor3f (0.1, 0.6, 1.);
	glLineWidth(3.);
	int i;
	//glBegin(GL_POINTS);
	glBegin(GL_LINE_STRIP);
	for (i=1; i<one_streamline.size()-1; i++)
	{
		std::pair<int, int> &one_px = one_streamline[i];
		glVertex2f ((GLfloat)one_px.second/(GLfloat)img_res, (GLfloat)one_px.first/(GLfloat)img_res);
	}
	glEnd();
	}


	if (ShowSingularities)
	   display_singularities();

	////Display the advected particles
	//int i, j;

	//glColor3f (1, 1, 0);
	//glPointSize(2.);
	//glBegin(GL_POINTS);
	//for (i=0; i<PDIM; i++)
	//{
	//	for (j=0; j<PDIM; j++)
	//		glVertex2f (particle_pos[i][j][0], particle_pos[i][j][1]);
	//}
	//glEnd();


	////// Display the streaklines

	//int k;

	//float hsv[3] = {0, 1., 1.};
	//float rgb[3];

	//for (i=0; i<PDIM; i++)
	//{
	//	for (j=0; j<PDIM; j++)
	//	{
	//		std::vector<int> &cur_pixel = pixels_to_particles[i*PDIM+j];

	//		hsv[0] = 360.*(float)(i*PDIM+j)/(float)(PDIM*PDIM);

	//		HsvRgb(hsv, rgb);

	//		glColor3fv(rgb);

	//		glBegin(GL_LINE_STRIP);
	//		//for (k=0; k<cur_pixel.size(); k++)
	//		for (k=cur_pixel.size()-1; k>=0; k--)
	//		{
	//			if (particle_active[cur_pixel[k]])
	//				glVertex2f (particles[cur_pixel[k]].first, particles[cur_pixel[k]].second);
	//		}
	//		glEnd();
	//	}
	//}

	glutSwapBuffers();


	// be sure the graphics buffer has been sent:
	// note: be sure to use glFlush() here, not glFinish() !

	glFlush();
}


/*------------------------------------------*/
//Normalize the whole field
double DistanceThreshold = 1e-5;
double  dmax   = 3/512.;
void normalize_Field()
{
	int i;
    double r;
	Vertex *cur_v;

	for(i = 0; i < poly->nverts; i++)
	{
		cur_v = poly->vlist[i];
		r = cur_v->vx*cur_v->vx + cur_v->vy*cur_v->vy;
		//r *= r;
					
		if (r < DistanceThreshold) 
		{
			r = DistanceThreshold;
			cur_v->vx *= dmax/r; 
			cur_v->vy *= dmax/r; 
		}

		r = cur_v->vx*cur_v->vx + cur_v->vy*cur_v->vy;
		r *= r;

		if (r > dmax*dmax) { 
			r  = sqrt(r); 
			cur_v->vx *= dmax/r; 
			cur_v->vy *= dmax/r; 
		}

		cur_v->vec.set(cur_v->vx, cur_v->vy);
	}
}

void    
rotate_vec( Polyhedron *this_poly, double degree)
{
	int i, j;

	for (i=0; i<this_poly->nverts; i++)
	{
		Vertex *v = this_poly->vlist[i];
		double tmp_v[2];
		tmp_v[0] = v->vx*cos(degree) - v->vy*sin(degree);
		tmp_v[1] = v->vx*sin(degree) + v->vy*cos(degree);

		v->vx = tmp_v[0];
		v->vy = tmp_v[1];
	}
}

void    
reflect_vec( Polyhedron *this_poly)
{
	int i, j;

	for (i=0; i<this_poly->nverts; i++)
	{
		Vertex *v = this_poly->vlist[i];
		double tmp_v[2];
		tmp_v[0] = v->vy;
		tmp_v[1] = v->vx;

		//v->vx = tmp_v[0];
		//v->vy = tmp_v[1];

		v->vx += tmp_v[0];
		v->vy += tmp_v[1];
	}

}

/*  The obtained color values will be over-ranged using this simple summation
*/
void   
combine_LIC()
{
	int i, j;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			LIC_tex_comb[i][j][0] = 0.5*(LIC_tex[i][j][0] + LIC_tex2[i][j][0]);
			LIC_tex_comb[i][j][1] = 0.5*(LIC_tex[i][j][1] + LIC_tex2[i][j][1]);
			LIC_tex_comb[i][j][2] = 0.5*(LIC_tex[i][j][2] + LIC_tex2[i][j][2]);
		}
	}

}


////////////////////////////////////////////////

////////////////////////////////////////////////


// compute forward direction 
void    
count_streamline_forward (int i, int j, int L, unsigned int forward_count[img_res][img_res])
{
	//float jitter = 0.5*(float)rand()/32768.0;
	//float x = i+.5 + jitter;
	//jitter = 0.5*(float)rand()/32768.0;
	//float y = j+.5 + jitter;
	float jitter = 0.5*(float)rand()/32768.0;
	float y = i+.5 + jitter;
	jitter = 0.5*(float)rand()/32768.0;
	float x = j+.5 + jitter;

	float vec[2];

	for (int step=0; step < L; step++)
	{
		// get the vector value at current location (decode the vec_img)
		vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 1.e-6) return;

		normalize_vec(vec);
		x = x+vec[0];
		y = y+vec[1];

		//i = (int) x;
		//j = (int) y;
		j = (int) x;
		i = (int) y;

		forward_count[i][j] += 1;

		if (i<0 || i>=img_res || j<0 || j>=img_res)
			return;
	}
}

void    
count_streamline_backward (int i, int j, int L, unsigned int backward_count[img_res][img_res])
{
	//float jitter = 0.5*(float)rand()/32768.0;
	//float x = i+.5 + jitter;
	//jitter = 0.5*(float)rand()/32768.0;
	//float y = j+.5 + jitter;
	float jitter = 0.5*(float)rand()/32768.0;
	float y = i+.5 + jitter;
	jitter = 0.5*(float)rand()/32768.0;
	float x = j+.5 + jitter;

	float vec[2];

	for (int step=0; step < L; step++)
	{
		// get the vector value at current location (decode the vec_img)
		vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 1.e-6) return;

		normalize_vec(vec);
		x = x-vec[0];
		y = y-vec[1];

		//i = (int) x;
		//j = (int) y;
		j = (int) x;
		i = (int) y;

		backward_count[i][j] += 1;

		if (i<0 || i>=img_res || j<0 || j>=img_res)
			return;
	}
}

void    
count_passed_pixels_from (int i, int j, int L,
		unsigned int forward_count[img_res][img_res], unsigned int backward_count[img_res][img_res])
{
	//passed_pixels.clear();
	//std::pair<int, int> start_pixel;
	//start_pixel.first = i;
	//start_pixel.second = j;
	//passed_pixels.push_back (start_pixel);

	for (int ii=0; ii<1; ii++)
	{
		count_streamline_forward (i, j, L, forward_count);
		count_streamline_backward (i, j, L, backward_count);
	}
}



// update the counters via LIC method
void 
comp_LIC_counter(unsigned int forward_[img_res][img_res], unsigned int backward_[img_res][img_res])
{
	int L=40;
	for (int i = 0; i<img_res; i++)
		for (int j=0; j<img_res; j++)
		{
			count_passed_pixels_from(i, j, L, forward_, backward_);
		}
}

void 
init_counters()
{
	int i, j;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			forward_counter[i][j]=0;
			backward_counter[i][j] = 0;
		}
	}
}


void    
get_color_map_for_counters()
{
	int i, j;
	int max_c = -1, min_c = INT_MAX;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			if (forward_counter[i][j] > max_c) 
				max_c = forward_counter[i][j];
			if (forward_counter[i][j] < min_c)
				min_c = forward_counter[i][j];
		}
	}

	// get the color map for the forward counter
	int diff = max_c - min_c;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			forward_counter_color[i][j][0] = 255 * (float)(forward_counter[i][j] -  min_c)/(float)diff;
			forward_counter_color[i][j][1] = 0;
			forward_counter_color[i][j][2] = 0/*255 * (float)(forward_counter[i][j] -  min_c)/(float)diff*/;
		}
	}

	max_c = -1; min_c = INT_MAX;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			if (backward_counter[i][j] > max_c) 
				max_c = backward_counter[i][j];
			if (backward_counter[i][j] < min_c)
				min_c = backward_counter[i][j];
		}
	}

	// get the color map for the backward counter
	diff = max_c - min_c;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			backward_counter_color[i][j][0] = 0;
			backward_counter_color[i][j][1] = 255 * (float)(backward_counter[i][j] -  min_c)/(float)diff;
			backward_counter_color[i][j][2] = 0;
		}
	}


	// combine 
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			comb_counter_color[i][j][0] = forward_counter_color[i][j][0];
			comb_counter_color[i][j][1] = backward_counter_color[i][j][1];
			comb_counter_color[i][j][2] = 0;
		}
	}

	// combine with LIC

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			comb_counter_color[i][j][0] = 0.7*forward_counter_color[i][j][0] + 0.3 * LIC_tex[i][j][0];
			comb_counter_color[i][j][1] = 0.7*backward_counter_color[i][j][1] + 0.3 * LIC_tex[i][j][1];
			comb_counter_color[i][j][2] = 0.3 * LIC_tex[i][j][2];
		}
	}
}

void    
get_color_map_for_counters_forward()
{
	int i, j;
	int max_c = -1, min_c = INT_MAX;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			if (forward_counter[i][j] > max_c) 
				max_c = forward_counter[i][j];
			if (forward_counter[i][j] < min_c)
				min_c = forward_counter[i][j];
		}
	}

	// get the color map for the forward counter
	int diff = max_c - min_c;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			forward_counter_color[i][j][0] = 255 * (float)(forward_counter[i][j] -  min_c)/(float)diff;
			forward_counter_color[i][j][1] = 0;
			forward_counter_color[i][j][2] = 0/*255 * (float)(forward_counter[i][j] -  min_c)/(float)diff*/;
		}
	}


	// combine 
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			comb_counter_color[i][j][0] = forward_counter_color[i][j][0];
			comb_counter_color[i][j][1] = 0;
			comb_counter_color[i][j][2] = 0;
		}
	}

	// combine with LIC

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			comb_counter_color[i][j][0] = 0.7*forward_counter_color[i][j][0] + 0.3 * LIC_tex[i][j][0];
			comb_counter_color[i][j][1] = 0.3 * LIC_tex[i][j][1];
			comb_counter_color[i][j][2] = 0.3 * LIC_tex[i][j][2];
		}
	}
}

void 
init_rotation_sum()
{
	int i, j;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			total_rotation [i][j] = 0.;
			visited_pixels[i][j] = false;
		}
	}
}


void    
sum_rotation_forward (int i, int j, int L, float total_rotation[img_res][img_res], bool &closedloop)
{
	// with jitter starting position
	//float jitter = 0.5*(float)rand()/32768.0;
	//float x = i+.5 + jitter;
	//jitter = 0.5*(float)rand()/32768.0;
	//float y = j+.5 + jitter;

	// without jittering
	//float x = i + 0.5;
	//float y = j + 0.5;
	float x = j + 0.5;
	float y = i + 0.5;

	// store the start position of each pixel in the world coordinate system for later computation
	start_pos[i][j][0] = x;
	start_pos[i][j][1] = y;

	float vec[2], vec2[2];
	float pre_vec[2];

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
	get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	double pre_ang, cur_ang, ang_diff;

	pre_ang = atan2(vec[1], vec[0]);

	closedloop = false;
	for (int step=0; step < L; step++)
	{
		if (fixedPt_pixels[i][j]) break;
		// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
		get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 5e-6) /*return*/break;

		normalize_vec(vec);

		float tx = x + vec[0];
		float ty = y + vec[1];

		int tj = (int)tx;
		int ti = (int)ty;
		if (fixedPt_pixels[ti][tj]) break;
		//vec2[0] = min_vx + (max_vx - min_vx) * vec_img[ti][tj][0]/255;
		//vec2[1] = min_vy + (max_vy - min_vy) * vec_img[ti][tj][1]/255;
		get_vec_at_regular_grid_image_based (tx, ty, img_res, img_res, vec2[0], vec2[1]);

		normalize_vec(vec2);

		vec[0] = 0.5*(vec[0] + vec2[0]);
		vec[1] = 0.5*(vec[1] + vec2[1]);

		normalize_vec(vec);

		x = x+vec[0];
		y = y+vec[1];

		//i = (int) x;
		//j = (int) y;
		j = (int) x;
		i = (int) y;

		std::pair<int, int> one_px;
		one_px.first = i;
		one_px.second = j;
		one_streamline.push_back(one_px);

		//icVector2 vec1(vec[0], vec[1]), vec2(pre_vec[0], pre_vec[1]);
		//double pre_ang = atan2 (pre_vec[1], pre_vec[0]);
		cur_ang = atan2 (vec[1], vec[0]);

		ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		//pre_vec[0] = vec[0];
		//pre_vec[1] = vec[1];
		pre_ang = cur_ang;
		
		//total_rotation[i][j] += ang_diff;
		rot_sum_tmp += ang_diff;

		if (i<=0 || i>=img_res-1 || j<=0 || j>=img_res-1)
			/*return*/break;

		/*if (step > 2 && i==start_i && j==start_j)*/ 
		if (rot_sum_tmp>=2*PI && abs(i-start_i)<4 && abs(j-start_j)<4) { closedloop=true; break;}  // to counter the closed loops
	}

	total_rotation[start_i][start_j] += rot_sum_tmp;
}

void    
sum_rotation_backward (int i, int j, int L, float total_rotation[img_res][img_res], bool &closedloop)
{
	// with jitter starting position
	//float jitter = 0.5*(float)rand()/32768.0;
	//float x = i+.5 + jitter;
	//jitter = 0.5*(float)rand()/32768.0;
	//float y = j+.5 + jitter;

	// without jittering
	//float x = i + 0.5;
	//float y = j + 0.5;
	float x = j + 0.5;
	float y = i + 0.5;

	float vec[2], vec2[2];
	float pre_vec[2];

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
	get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;
	double pre_ang, cur_ang, ang_diff;
	pre_ang = atan2(vec[1], vec[0]);

	closedloop = false;

	for (int step=0; step < L; step++)
	{
		// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
		if (fixedPt_pixels[i][j]) break;

		get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 5e-6) /*return*/break;

		normalize_vec(vec);
		float tx = x - vec[0];
		float ty = y - vec[1];

		int tj = (int)tx;
		int ti = (int)ty;
		if (fixedPt_pixels[ti][tj]) break;

		//vec2[0] = min_vx + (max_vx - min_vx) * vec_img[ti][tj][0]/255;
		//vec2[1] = min_vy + (max_vy - min_vy) * vec_img[ti][tj][1]/255;
		get_vec_at_regular_grid_image_based (tx, ty, img_res, img_res, vec2[0], vec2[1]);

		normalize_vec(vec2);

		vec[0] = 0.5*(vec[0] + vec2[0]);
		vec[1] = 0.5*(vec[1] + vec2[1]);

		normalize_vec(vec);

		x = x-vec[0];
		y = y-vec[1];

		//i = (int) x;
		//j = (int) y;
		j = (int) x;
		i = (int) y;
		
		std::pair<int, int> one_px;
		one_px.first = i;
		one_px.second = j;
		one_streamline.push_back(one_px);

		//icVector2 vec1(vec[0], vec[1]), vec2(pre_vec[0], pre_vec[1]);
		//double pre_ang = atan2 (pre_vec[1], pre_vec[0]);
		cur_ang = atan2 (vec[1], vec[0]);

		ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		//pre_vec[0] = vec[0];
		//pre_vec[1] = vec[1];

		pre_ang = cur_ang;

		//total_rotation[i][j] -= ang_diff;
		//rot_sum_tmp +=ang_diff;
		rot_sum_tmp -=ang_diff;

		if (i<=0 || i>=img_res-1 || j<=0 || j>=img_res-1)
			/*return*/break;

		if (rot_sum_tmp>=2*PI && abs(i-start_i)<4 && abs(j-start_j)<4) {closedloop = true; break;}  // to counter the closed loops
	}

	total_rotation[start_i][start_j] += rot_sum_tmp;
}


void 
comp_LIC_total_rotation(float total_rotation[img_res][img_res])
{
	//int L=20;
	int i, j;
	//for (i = 0; i<img_res; i++)
	//{
	//	for (j=0; j<img_res; j++)
	//	{
//#pragma omp parallel for 
	for (i=0; i<img_res*img_res; i++)
	{
		int row = i/img_res;
		int col = i%img_res;
			//if (fixedPt_pixels[i][j]) {total_rotation[i][j]=0.; continue;}
			if (fixedPt_pixels[row][col]) {total_rotation[row][col]=0.; continue;}

			bool closedloop = false;

			//one_streamline.clear();
			// we have to define them as the local variables so that each thread will modify
			// only the local ones not the global ones!!!!
			forward_streamline.clear();
			backward_streamline.clear();
			forward_rot.clear();
			backward_rot.clear();

			std::vector<std::pair<int, int>> one_streamline, forward_streamline, backward_streamline;
			std::vector<float> current_total_rot, forward_rot, backward_rot;

			//sum_rotation_forward (i, j, L, total_rotation, closedloop);
			double rotsum =0.;
			if (UseRK23 == 0)
				//rotsum = sum_rotation_forward_with_streamline (i, j, L, closedloop, false);
				rotsum = sum_rotation_forward_with_streamline (row, col, L, closedloop, false, forward_streamline,
				forward_rot);
			else
				//rotsum = sum_rotation_forward_with_streamline_RK23 (i, j, L, closedloop, false);
				rotsum = sum_rotation_forward_with_streamline_RK23 (row, col, L, closedloop, false,
				forward_streamline, forward_rot);

			if (!closedloop)
			{
				if (UseRK23 == 0)
					//sum_rotation_backward(i, j, L, total_rotation, closedloop);
					//rotsum += sum_rotation_backward_with_streamline (i, j, L, closedloop, false); 
					rotsum += sum_rotation_backward_with_streamline (row, col, L, closedloop, false, 
					backward_streamline, backward_rot); 
				else
					//rotsum += sum_rotation_backward_with_streamline_RK23 (i, j, L, closedloop, false); 
					rotsum += sum_rotation_backward_with_streamline_RK23 (row, col, L, closedloop, false,
					backward_streamline, backward_rot); 
			}

			//total_rotation[i][j] = rotsum;
			total_rotation[row][col] = rotsum;
		}
	//}

	// average the total rotation
	
	//for (i = 0; i<img_res; i++)
	//{
	//	for (j=0; j<img_res; j++)
	//		total_rotation[i][j] /= (forward_counter[i][j]+backward_counter[i][j]);
	//}
	
}

void    comp_LIC_total_rotation_opt(float total_rotation[img_res][img_res])
{
	  // default L is 1000
	int i, j;

//#pragma omp parallel for //private(loc_one_streamline, loc_forward_streamline, loc_backward_streamline, loc_current_total_rot, loc_forward_rot, loc_backward_rot)
	//for (i = 0; i<img_res; i++)
	//{
	//	for (j=0; j<img_res; j++)
	//	{
 	for (i=0; i<img_res*img_res; i++)
	{
		int row = i/img_res;
		int col = i%img_res;

			//if (visited_pixels[i][j]) continue;
			//if (fixedPt_pixels[i][j]) {total_rotation[i][j]=0.; continue;}
			
		    if (visited_pixels[row][col]) continue;
			if (fixedPt_pixels[row][col]) {total_rotation[row][col]=0.; continue;}


			//forward_streamline.clear();
			//backward_streamline.clear();
			//forward_rot.clear();
			//backward_rot.clear();

			std::vector<std::pair<int, int>> loc_one_streamline, loc_forward_streamline, loc_backward_streamline;
			std::vector<float> loc_current_total_rot, loc_forward_rot, loc_backward_rot;

			bool closedloop = false;
			//predict_stepsize = 1.;

			//sum_rotation_forward (i, j, L, total_rotation, closedloop);
			//printf ("i=%d\n", i);
			double rotsum;
			if (UseRK23 == 0)
				//rotsum = sum_rotation_forward_with_streamline (i, j, L, closedloop, false);
				rotsum = sum_rotation_forward_with_streamline (row, col, L, closedloop, false,
				loc_forward_streamline, loc_forward_rot);
			else
				//rotsum = sum_rotation_forward_with_streamline_RK23 (i, j, L, closedloop, false);
				rotsum = sum_rotation_forward_with_streamline_RK23 (row, col, L, closedloop, false,
				loc_forward_streamline, loc_forward_rot);

			if (!closedloop)
			{
				if (UseRK23==0)
					//rotsum += sum_rotation_backward_with_streamline (i, j, L, closedloop, false);
					rotsum += sum_rotation_backward_with_streamline (row, col, L, closedloop, false,
					loc_backward_streamline, loc_backward_rot);
				else
					//rotsum += sum_rotation_backward_with_streamline_RK23 (i, j, L, closedloop, false);
					rotsum += sum_rotation_backward_with_streamline_RK23 (row, col, L, closedloop, false,
					loc_backward_streamline, loc_backward_rot);
			}

			merge_current_forward_back_streamlines(loc_one_streamline, loc_forward_streamline, loc_backward_streamline,
				loc_current_total_rot, loc_forward_rot, loc_backward_rot);
			std::pair<int, int> one_px;
			//one_px.first = i;
			//one_px.second = j;
			one_px.first = row;
			one_px.second = col;
			loc_one_streamline.push_back(one_px);

			// set the same rotation value for all the pixels that the current streamline passes
			//#pragma omp critical
			for (int k=0; k<loc_one_streamline.size(); k++)
			{
				// now we need to take a percentage of the rotation along the specified arc length
				int rot_pos = L_percentage * (loc_current_total_rot.size()-1);

				#pragma omp critical
				total_rotation[loc_one_streamline[k].first][loc_one_streamline[k].second] = 
					/*rotsum*/ /*total_rotation[i][j]*/loc_current_total_rot[rot_pos];
				#pragma omp critical
				visited_pixels[loc_one_streamline[k].first][loc_one_streamline[k].second] = true;
			}
		//}
			loc_one_streamline.clear();
			loc_forward_streamline.clear();
			loc_backward_streamline.clear();
			loc_forward_rot.clear();
			loc_backward_rot.clear();
			loc_current_total_rot.clear();
	}
}


double    
sum_rotation_forward_with_streamline (int i, int j, int L, bool &closedloop, bool store, 
std::vector<std::pair<int, int>> &forward_streamline, std::vector<float> &forward_rot)
{
	// with jitter starting position
	float jitter = 0.5*(float)rand()/32768.0;
	float x = j+.5 + jitter;
	jitter = 0.5*(float)rand()/32768.0;
	float y = i+.5 + jitter;

	// without jittering
	////float x = i + 0.5;
	////float y = j + 0.5;
	//float x = j + 0.5;
	//float y = i + 0.5;

	// store the start position of each pixel in the world coordinate system for later computation
	start_pos[i][j][0] = x;
	start_pos[i][j][1] = y;

	float vec[2], vec2[2];
	float pre_vec[2];

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

	//get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);
	get_vec_at_regular_grid_image_based ((x-.5), (y-.5), img_res, img_res, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	double pre_ang, cur_ang, ang_diff;

	pre_ang = atan2(vec[1], vec[0]);

	closedloop = false;
	int step;

#pragma omp single 
	for (step=0; step < L; step++)
	{
		if (fixedPt_pixels[i][j]) break;
		// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
		//get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);
		get_vec_at_regular_grid_image_based ((x-.5), (y-.5), img_res, img_res, vec[0], vec[1]);

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 5.e-6) break;

		normalize_vec(vec);

		///*  DEBUG EULER INTEGRATOR
		//float tx = x + .25*vec[0];
		//float ty = y + .25*vec[1];
		float tx = x + .5*vec[0];
		float ty = y + .5*vec[1];
		//float tx = x + vec[0];
		//float ty = y + vec[1];

		int tj = (int)tx;
		int ti = (int)ty;
		if (fixedPt_pixels[ti][tj]) break;
		//vec2[0] = min_vx + (max_vx - min_vx) * vec_img[ti][tj][0]/255;
		//vec2[1] = min_vy + (max_vy - min_vy) * vec_img[ti][tj][1]/255;
		//get_vec_at_regular_grid_image_based (tx, ty, img_res, img_res, vec2[0], vec2[1]);
		get_vec_at_regular_grid_image_based ((tx-.5), (ty-.5), img_res, img_res, vec2[0], vec2[1]);

		normalize_vec(vec2);

		vec[0] = 0.5*(vec[0] + vec2[0]);
		vec[1] = 0.5*(vec[1] + vec2[1]);

		normalize_vec(vec);
		//*/

		x = x+vec[0];
		y = y+vec[1];
		//x = x+.5*vec[0];   // Note that small step size works fine for vortex structure
		//y = y+.5*vec[1];

		//i = (int) x;
		//j = (int) y;
		j = (int) x;
		i = (int) y;

		//if (store)
		//{
		std::pair<int, int> one_px;
		one_px.first = i;
		one_px.second = j;
		//one_streamline.push_back(one_px);
		forward_streamline.push_back(one_px); // this is not a global variable, change that!
		//}
		//icVector2 vec1(vec[0], vec[1]), vec2(pre_vec[0], pre_vec[1]);
		//double pre_ang = atan2 (pre_vec[1], pre_vec[0]);
		cur_ang = atan2 (vec[1], vec[0]);

		ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		//pre_vec[0] = vec[0];
		//pre_vec[1] = vec[1];
		pre_ang = cur_ang;
		
		//total_rotation[i][j] += ang_diff;
		rot_sum_tmp += ang_diff;

		//current_total_rot.push_back(rot_sum_tmp);
		//if (store)
		forward_rot.push_back(rot_sum_tmp);

		if (i<=0 || i>=img_res-1 || j<=0 || j>=img_res-1)
			break;
		
		/*if (step > 2 && i==start_i && j==start_j)*/ 
		if (abs(rot_sum_tmp)>=2*PI && abs(i-start_i)<3 && abs(j-start_j)<3) {closedloop=true; break;}  // to counter the closed loops
	}

	return rot_sum_tmp;
}


double    
sum_rotation_forward_with_streamline_RK23 (int i, int j, int L, bool &closedloop, bool store,
std::vector<std::pair<int, int>> &forward_streamline, std::vector<float> &forward_rot)
{
	// with jitter starting position
	float jitter = 0.5*(float)rand()/32768.0;
	float x = j+.5 + jitter;
	jitter = 0.5*(float)rand()/32768.0;
	float y = i+.5 + jitter;

	// without jittering
	////float x = i + 0.5;
	////float y = j + 0.5;
	//float x = j + 0.5;
	//float y = i + 0.5;

	// store the start position of each pixel in the world coordinate system for later computation
	start_pos[i][j][0] = x;
	start_pos[i][j][1] = y;

	float vec[2], vec2[2];
	float pre_vec[2];

	get_vec_at_regular_grid_image_based ((x-.5), (y-.5), img_res, img_res, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	double pre_ang, cur_ang, ang_diff;

	pre_ang = atan2(vec[1], vec[0]);

	closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

#pragma omp single 
	for (step=0; step < L; step++)
	{
		if (fixedPt_pixels[i][j]) break;

		get_nextpt_RK23_quad(trace_gp, trace_gp, offset, img_res, img_res, false);
		//get_nextpt_RK45_quad(trace_gp, trace_gp, offset, img_res, img_res, false);

		j = (int) trace_gp[0];
		i = (int) trace_gp[1];

		std::pair<int, int> one_px;
		one_px.first = i;
		one_px.second = j;
		forward_streamline.push_back(one_px);


		//cur_ang = atan2 (vec[1], vec[0]);
		cur_ang = atan2 (offset[1], offset[0]);

		ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		pre_ang = cur_ang;
		
		//total_rotation[i][j] += ang_diff;
		rot_sum_tmp += ang_diff;

		forward_rot.push_back(rot_sum_tmp);

		if (i<=0 || i>=img_res-1 || j<=0 || j>=img_res-1)
			break;
		
		if (abs(rot_sum_tmp)>=2*PI && abs(i-start_i)<3 && abs(j-start_j)<3) {closedloop=true; break;}  // to counter the closed loops
	}

	return rot_sum_tmp;
}

double    
sum_rotation_backward_with_streamline (int i, int j, int L, bool &closedloop, bool store,
std::vector<std::pair<int, int>> &backward_streamline, std::vector<float> &backward_rot)
{
	// with jitter starting position
	float jitter = 0.5*(float)rand()/32768.0;
	float x = j+.5 + jitter;
	jitter = 0.5*(float)rand()/32768.0;
	float y = i+.5 + jitter;

	// without jittering
	//////float x = i + 0.5;
	//////float y = j + 0.5;
	//float x = j + 0.5;
	//float y = i + 0.5;

	float vec[2], vec2[2];
	float pre_vec[2];

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
	//get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);
	get_vec_at_regular_grid_image_based ((x-.5), (y-.5), img_res, img_res, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;
	double pre_ang, cur_ang, ang_diff;
	pre_ang = atan2(vec[1], vec[0]);

	closedloop = false;
	int step;

#pragma omp single 
	for (step=0; step < L; step++)
	{
		if (fixedPt_pixels[i][j]) break;
		// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;
		//get_vec_at_regular_grid_image_based (x, y, img_res, img_res, vec[0], vec[1]);
		get_vec_at_regular_grid_image_based ((x-.5), (y-.5), img_res, img_res, vec[0], vec[1]);

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		if (get_vec_len(vec) < 5.e-6) break;

		normalize_vec(vec);

		///* DEBUG EULER INTEGRATOR
		//float tx = x - .25*vec[0];
		//float ty = y - .25*vec[1];
		float tx = x - .5*vec[0];
		float ty = y - .5*vec[1];
		//float tx = x - vec[0];
		//float ty = y - vec[1];

		int tj = (int)tx;
		int ti = (int)ty;
		if (fixedPt_pixels[ti][tj]) break;
		//vec2[0] = min_vx + (max_vx - min_vx) * vec_img[ti][tj][0]/255;
		//vec2[1] = min_vy + (max_vy - min_vy) * vec_img[ti][tj][1]/255;
		//get_vec_at_regular_grid_image_based (tx, ty, img_res, img_res, vec2[0], vec2[1]);
		get_vec_at_regular_grid_image_based ((tx-.5), (ty-.5), img_res, img_res, vec2[0], vec2[1]);

		normalize_vec(vec2);

		vec[0] = 0.5*(vec[0] + vec2[0]);
		vec[1] = 0.5*(vec[1] + vec2[1]);

		normalize_vec(vec);
		//*/

		x = x-vec[0];
		y = y-vec[1];
		//x = x-.5*vec[0];
		//y = y-.5*vec[1];

		//i = (int) x;
		//j = (int) y;
		j = (int) x;
		i = (int) y;

		//if (store)
		//{
		std::pair<int, int> one_px;
		one_px.first = i;
		one_px.second = j;
		//one_streamline.push_back(one_px);
		backward_streamline.push_back(one_px);
		//}

		//icVector2 vec1(vec[0], vec[1]), vec2(pre_vec[0], pre_vec[1]);
		//double pre_ang = atan2 (pre_vec[1], pre_vec[0]);
		cur_ang = atan2 (vec[1], vec[0]);

		ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		//pre_vec[0] = vec[0];
		//pre_vec[1] = vec[1];

		pre_ang = cur_ang;

		//total_rotation[i][j] -= ang_diff;
		//rot_sum_tmp +=ang_diff;
		rot_sum_tmp -=ang_diff;

		//if (store)
		backward_rot.push_back(rot_sum_tmp);

		if (i<=0 || i>=img_res-1 || j<=0 || j>=img_res-1)
			break;

		/*if (step > 2 && i==start_i && j==start_j)*/
		if (abs(rot_sum_tmp)>=2*PI && abs(i-start_i)<3 && abs(j-start_j)<3){closedloop=true; break;}  // to counter the closed loops
	}

	return rot_sum_tmp;
}

double    
sum_rotation_backward_with_streamline_RK23 (int i, int j, int L, bool &closedloop, bool store,
std::vector<std::pair<int, int>> &backward_streamline, std::vector<float> &backward_rot)
{
	// with jitter starting position
	float jitter = 0.5*(float)rand()/32768.0;
	float x = j+.5 + jitter;
	jitter = 0.5*(float)rand()/32768.0;
	float y = i+.5 + jitter;

	// without jittering
	//////float x = i + 0.5;
	//////float y = j + 0.5;
	//float x = j + 0.5;
	//float y = i + 0.5;

	float vec[2], vec2[2];
	float pre_vec[2];

	get_vec_at_regular_grid_image_based ((x-.5), (y-.5), img_res, img_res, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;
	double pre_ang, cur_ang, ang_diff;
	pre_ang = atan2(vec[1], vec[0]);

	closedloop = false;
	int step;

	double trace_gp[2] = {x, y};
	double offset[2] = {0.};

#pragma omp single 
	for (step=0; step < L; step++)
	{
		if (fixedPt_pixels[i][j]) break;
		get_nextpt_RK23_quad(trace_gp, trace_gp, offset, img_res, img_res, true);
		//get_nextpt_RK45_quad(trace_gp, trace_gp, offset, img_res, img_res, true);

		j = (int) trace_gp[0];
		i = (int) trace_gp[1];

		std::pair<int, int> one_px;
		one_px.first = i;
		one_px.second = j;
		backward_streamline.push_back(one_px);

		cur_ang = atan2 (offset[1], offset[0]);

		ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		pre_ang = cur_ang;

		rot_sum_tmp -=ang_diff;

		backward_rot.push_back(rot_sum_tmp);

		if (i<=0 || i>=img_res-1 || j<=0 || j>=img_res-1)
			break;

		/*if (step > 2 && i==start_i && j==start_j)*/
		if (abs(rot_sum_tmp)>=2*PI && abs(i-start_i)<3 && abs(j-start_j)<3){closedloop=true; break;}  // to counter the closed loops
	}

	return rot_sum_tmp;
}

void    
trace_streamline_from_image_based(int i, int j, int L)
{

	forward_streamline.clear();
	backward_streamline.clear();
	forward_rot.clear();
	backward_rot.clear();

	// Add the first pixel
	std::pair<int, int> one_px;
	one_px.first = i;
	one_px.second = j;
	forward_streamline.push_back(one_px);


	double rot_sum = 0;
	bool closedloop=false;
	if (UseRK23==0)
		rot_sum += sum_rotation_forward_with_streamline (i, j, L, closedloop, true, forward_streamline, forward_rot);
	else
		rot_sum += sum_rotation_forward_with_streamline_RK23 (i, j, L, closedloop, true, forward_streamline, forward_rot);

	if (!closedloop)
	{
		if (UseRK23==0)
			rot_sum += sum_rotation_backward_with_streamline(i, j, L, closedloop, true,
			backward_streamline, backward_rot);
		else
			rot_sum += sum_rotation_backward_with_streamline_RK23(i, j, L, closedloop, true,
			backward_streamline, backward_rot);
	}

	printf("the total rotation of this streamline from pixel (%d, %d) is %f\n",i, j, rot_sum);


	// sort the streamline from start to the end
}


void 
merge_current_forward_back_streamlines(std::vector<std::pair<int, int>> &one_streamline, 
				std::vector<std::pair<int, int>> forward_streamline, 
				std::vector<std::pair<int, int>> backward_streamline,
				std::vector<float>  &current_total_rot, 
				std::vector<float>	forward_rot, 
				std::vector<float>	backward_rot
				)
{
	one_streamline.clear();
	current_total_rot.clear();
	int i;
	float total_rot = 0;
	if( backward_rot.size()>0)
	{
		for (i=/*backward_streamline*/backward_rot.size()-1; i>=1; i--)
		{
			one_streamline.push_back(backward_streamline[i]);
			total_rot += (backward_rot[i]-backward_rot[i-1]);
			current_total_rot.push_back(total_rot);
		}
		one_streamline.push_back(backward_streamline[0]);
		current_total_rot.push_back(total_rot);

		one_streamline.push_back(forward_streamline[0]);
		for (i=1; i</*current_total_rot*/forward_rot.size(); i++)
		{
			one_streamline.push_back(forward_streamline[i]);
			total_rot += (forward_rot[i]-forward_rot[i-1]);
			current_total_rot.push_back(total_rot);
		}
	}
	else
	{
		for (i=0; i<forward_rot.size(); i++)
		{
			one_streamline.push_back(forward_streamline[i]);
			total_rot += (forward_rot[i]-forward_rot[i-1]);
			current_total_rot.push_back(total_rot);
		}
		
	}
}

void 
merge_current_forward_back_streamlines()
{
	one_streamline.clear();
	current_total_rot.clear();
	int i;
	float total_rot = 0;
	for (i=/*backward_streamline*/backward_rot.size()-1; i>=1; i--)
	{
		one_streamline.push_back(backward_streamline[i]);
		total_rot += (backward_rot[i]-backward_rot[i-1]);
		current_total_rot.push_back(total_rot);
	}
	one_streamline.push_back(backward_streamline[0]);
	current_total_rot.push_back(total_rot);

	one_streamline.push_back(forward_streamline[0]);
	for (i=1; i</*current_total_rot*/forward_rot.size(); i++)
	{
		one_streamline.push_back(forward_streamline[i]);
		total_rot += (forward_rot[i]-forward_rot[i-1]);
		current_total_rot.push_back(total_rot);
	}
}

void   
get_color_map_for_rot_sum()
{
	int i, j;
	float max_rot=-1.e8, min_rot=1.e8;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{	
			if (total_rotation[i][j] > max_rot) max_rot = total_rotation[i][j];

			if (total_rotation[i][j] < min_rot) min_rot = total_rotation[i][j];
		}
	}

	// convert to color map
	float hsv[3]={0.}, rgb[3]={0.};
	float range = max_rot - min_rot;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			hsv[2] = 1.;

			hsv[1] = 1. - 2.*(total_rotation[i][j] - min_rot)/range;

			if (hsv[1] < 0)
			{
				hsv[0] = 0.;
				hsv[1] = abs(hsv[1]);
			}
			else
				hsv[0] = 240.;

			//if(total_rotation[i][j]>max_rot/4. || total_rotation[i][j]<min_rot/4.)
			//	hsv[1] = 1;

			HsvRgb (hsv, rgb);
			rot_sum_color[i][j][0] = (unsigned char)(rgb[0]*255.);
			rot_sum_color[i][j][1] = (unsigned char)(rgb[1]*255.);
			rot_sum_color[i][j][2] = (unsigned char)(rgb[2]*255.);
		}
	}

	// combine with LIC

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			//rot_sum_color[i][j][0] = 0.7*rot_sum_color[i][j][0] + 0.3 * LIC_tex[i][j][0];
			//rot_sum_color[i][j][1] = 0.7*rot_sum_color[i][j][1] + 0.3 * LIC_tex[i][j][1];
			//rot_sum_color[i][j][2] = 0.7*rot_sum_color[i][j][2] + 0.3 * LIC_tex[i][j][2];
			rot_sum_color[i][j][0] = 0.8*rot_sum_color[i][j][0] + 0.2 * LIC_tex[i][j][0];
			rot_sum_color[i][j][1] = 0.8*rot_sum_color[i][j][1] + 0.2 * LIC_tex[i][j][1];
			rot_sum_color[i][j][2] = 0.8*rot_sum_color[i][j][2] + 0.2 * LIC_tex[i][j][2];
		}
	}
}




///////////////////////////////////////////////////////////
void   
cal_vec_double_gyre_static(double x, double y, double &vx, double &vy)
{
	vx = -M_PI*sin(M_PI*x)*cos(M_PI*y);
	vy = M_PI*sin(M_PI*y)*cos(M_PI*x);
}


void    
gen_static_double_gyre()
{
	// First, create a small mesh say 64x64 mesh, then compute the vector values at the grid points of the mesh, render it to the vector image
	int i, j;
	double x1, x2, y, vx, vy;
	double dx = 2./63;
	float rgb[3]={0.};
	min_vx = 1.e8, max_vx=-1.e8, min_vy = 1.e8, max_vy=-1.e8;

	// first, find out the maximum and minimum vector value
	for (i=0; i<64; i++)
	{
		x1 = i*dx;

		for (j=0; j<64; j++)
		{
			y = j*dx;
			cal_vec_double_gyre_static(x1, y, vx, vy);

			if (vx < min_vx) min_vx = vx;
			if (vx > max_vx) max_vx = vx;
			if (vy < min_vy) min_vy = vy;
			if (vy > max_vy) max_vy = vy;
		}
	}

	double x_rang = max_vx - min_vx;
	double y_rang = max_vy - min_vy;

	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 2, 0, 2);
	glClear(GL_COLOR_BUFFER_BIT);

	glDrawBuffer(GL_BACK);

	for (i=0; i<63; i++)
	{
		x1 = i*dx;
		x2 = x1 + dx;

		glBegin(GL_QUAD_STRIP);
		for (j=0; j<64; j++)
		{
			y = j*dx;
			cal_vec_double_gyre_static(x1, y, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x1, y);

			cal_vec_double_gyre_static(x2, y, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x2, y);

		}
		glEnd();
	}
	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, vec_img);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_FLOAT, vec_img2);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
}

//void 
//get_vector_field( float x, float y, float z, float t, float &vxp, float &vyp, float &vzp ) 
//{ 
//	//a(t) and b(t)
//	float a_t = EPSILON * sin( OMEGA * PI * t);
//	float b_t = 1.0f - 2.0f * EPSILON * sin( OMEGA * PI * t); 
//
//	//f(x,t)
//	float f_x = a_t * pow(x, 2) + b_t * x;
//
//	//f'(x,t)
//	float df_x = 2.0f * a_t * x + b_t;
//
//	//u, v
//	vxp = -1.0f * PI * AAAA * sin(PI * f_x) * cos ( PI * y);
//	vyp = PI * AAAA * cos( PI * f_x) * sin( PI * y) * df_x;
//	vzp = 0;
//
//} 



void   
cal_vec_double_gyre(double x, double y, double t, double &vx, double &vy)
{
//	//a(t) and b(t)
	double a_t = EPSILON * sin( OMEGA * PI * t);
	double b_t = 1.0f - 2.0f * EPSILON * sin( OMEGA * PI * t); 
//
//	//f(x,t)
	double f_x = a_t * x * x + b_t * x;
//
//	//f'(x,t)
	double df_x = 2.0f * a_t * x + b_t;
//
//	//u, v
	vx = -PI * Amp * sin(PI * f_x) * cos ( PI * y);
	vy = PI * Amp * cos( PI * f_x) * sin( PI * y) * df_x;
}


void    
gen_double_gyre_at_time(double time)
{
	// First, create a small mesh say 64x64 mesh, then compute the vector values at the grid points of the mesh, render it to the vector image
	int i, j;
	double x1, x2, y, vx, vy;
	double dx = 2./63;
	float rgb[3]={0.};
	min_vx = 1.e8, max_vx=-1.e8, min_vy = 1.e8, max_vy=-1.e8;

	// first, find out the maximum and minimum vector value
	for (i=0; i<64; i++)
	{
		x1 = i*dx;

		for (j=0; j<64; j++)
		{
			y = j*dx;
			cal_vec_double_gyre(x1, y, time, vx, vy);

			if (vx < min_vx) min_vx = vx;
			if (vx > max_vx) max_vx = vx;
			if (vy < min_vy) min_vy = vy;
			if (vy > max_vy) max_vy = vy;
		}
	}

	double x_rang = max_vx - min_vx;
	double y_rang = max_vy - min_vy;

	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 2, 0, 2);
	glClear(GL_COLOR_BUFFER_BIT);

	glDrawBuffer(GL_BACK);

	for (i=0; i<63; i++)
	{
		x1 = i*dx;
		x2 = x1 + dx;

		glBegin(GL_QUAD_STRIP);
		for (j=0; j<64; j++)
		{
			y = j*dx;
			cal_vec_double_gyre(x1, y, time, vx, vy);
			// get the color for the current vertex
			//rgb[1] =  (vx - min_vx)/x_rang;
			//rgb[0] =  (vy - min_vy)/y_rang;
			rgb[0] =  (vx - min_vx)/x_rang;
			rgb[1] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x1, y);

			cal_vec_double_gyre(x2, y, time, vx, vy);
			// get the color for the current vertex
			//rgb[1] =  (vx - min_vx)/x_rang;
			//rgb[0] =  (vy - min_vy)/y_rang;
			rgb[0] =  (vx - min_vx)/x_rang;
			rgb[1] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x2, y);

		}
		glEnd();
	}
	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, vec_img);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_FLOAT, vec_img2);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
}

/* We try to compute the rotation of the pathline starting from a given position after a certain time determined by the number of integration steps right now
*/
void    
sum_rotation_forward_timeDep (int i, int j, int L, float total_rotation[img_res][img_res], double start_time)
{
	// with jitter starting position
	//float jitter = 0.5*(double)rand()/32768.0;
	//double x = i+.5 + jitter;
	//jitter = 0.5*(double)rand()/32768.0;
	//double y = j+.5 + jitter;

	// without jittering
	//double x = i + 0.5;
	//double y = j + 0.5;
	double y = i + 0.5;
	double x = j + 0.5;


	x *= 2;
	y *= 2;
	x /= img_res;
	y /= img_res;

	// store the start position of each pixel in the world coordinate system for later computation
	start_pos[i][j][0] = x;
	start_pos[i][j][1] = y;

	double vec[2];
	double pre_vec[2];

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

	cal_vec_double_gyre (x, y, start_time, vec[0], vec[1]);
	//cal_vec_ibfv_ex (x, y, start_time, vec[0], vec[1]);
	//cal_gyre_saddle	(x, y, start_time, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	double pre_ang = atan2 (pre_vec[1], pre_vec[0]);
	double cur_ang = atan2 (vec[1], vec[0]);

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;

	double d_time = T_window/(L-1);
	double half_dt = d_time/2.;

	for (int step=0; step < L; step++)
	{
		// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		// We need a better strategy to interpolate the vector value at (x, y, t) especially with t being considered now!!

		/*
		   We are now testing the flow generated by formula so that we do not need to worry about the interpolation between two sampled time steps
		   01/28/2013
		*/
		double cur_t = start_time+d_time*step;

		cal_vec_double_gyre (x, y, cur_t, vec[0], vec[1]);
		//cal_vec_ibfv_ex (x, y, cur_t, vec[0], vec[1]);
		//cal_gyre_saddle (x, y, cur_t, vec[0], vec[1]);

		if (get_vec_len(vec) < 1.e-6) break;


		// Euler integration
		//x = x+vec[0];
		//y = y+vec[1];
		//// Euler integration
		//normalize_vec(vec);
		//x = x+d_time*vec[0];
		//y = y+d_time*vec[1];

		// we should try to use RK4 for the integration.

		double tx1 = x+half_dt*vec[0];
		double ty1 = y+half_dt*vec[1];
		
		double vec1[2];
		cal_vec_double_gyre (tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
		//cal_vec_double_gyre (tx1, ty1, cur_t, vec1[0], vec1[1]);
		//cal_gyre_saddle (tx1, ty1, cur_t, vec1[0], vec1[1]);

		//cal_vec_ibfv_ex (tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
		//normalize_vec(vec1);
		//tx1 = x+half_dt*vec1[0];
		//ty1 = y+half_dt*vec1[1];

		////// RK2
		double total_vec[2];
		total_vec[0] = (vec[0]+vec1[0])/2.;
		total_vec[1] = (vec[1]+vec1[1])/2.;
		x = x+d_time*total_vec[0];
		y = y+d_time*total_vec[1];

		//double vec2[2];
		//cal_vec_double_gyre (tx1, ty1, cur_t+half_dt, vec2[0], vec2[1]);
		//normalize_vec(vec2);
		//tx1 = x+d_time*vec2[0];
		//ty1 = y+d_time*vec2[1];

		//double vec3[2];
		//cal_vec_double_gyre (tx1, ty1, cur_t+d_time, vec3[0], vec3[1]);
		//normalize_vec(vec3);

		////////// RK4
		//double total_vec[2];
		//total_vec[0] = (vec[0]+2.0*vec1[0]+2.0*vec2[0]+vec3[0])/6.;
		//total_vec[1] = (vec[1]+2.0*vec1[1]+2.0*vec2[1]+vec3[1])/6.;
		//x = x+d_time*total_vec[0];
		//y = y+d_time*total_vec[1];


		//i = (int) x*img_res/2;
		//j = (int) y*img_res/2;
		j = (int) x*img_res/2;
		i = (int) y*img_res/2;

		//icVector2 vec1(vec[0], vec[1]), vec2(pre_vec[0], pre_vec[1]);
		//cur_ang = atan2 (vec[1], vec[0]);
		cur_ang = atan2 (total_vec[1], total_vec[0]);

		double ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		pre_vec[0] = vec[0];
		pre_vec[1] = vec[1];
		pre_ang = cur_ang;
		
		//total_rotation[i][j] += ang_diff;
		rot_sum_tmp += ang_diff;

		//forward_counter[i][j] += 1;

		if (i<0 || i>=img_res || j<0 || j>=img_res)
			break;
	}

	total_rotation[start_i][start_j] += rot_sum_tmp;

	// record the end position of the particle starting at (i, j)
	end_positions_forward[start_i][start_j][0] = x;
	end_positions_forward[start_i][start_j][1] = y;
}

void    
sum_rotation_backward_timeDep (int i, int j, int L, float total_rotation[img_res][img_res], double start_time)
{
	// with jitter starting position
	double jitter = 0.5*(double)rand()/32768.0;
	double y = i+.5 + jitter;
	jitter = 0.5*(double)rand()/32768.0;
	double x = j+.5 + jitter;

	// without jittering
	//float x = i + 0.5;
	//float y = j + 0.5;
	x /= img_res;
	y /= img_res;

	// store the start position of each pixel in the world coordinate system for later computation
	start_pos[i][j][0] = x;
	start_pos[i][j][1] = y;

	double vec[2];
	double pre_vec[2];

	//pre_vec[0] = vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
	//pre_vec[1] = vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

	cal_vec_double_gyre (x, y, start_time, vec[0], vec[1]);
	//cal_vec_ibfv_ex (x, y, start_time, vec[0], vec[1]);
	pre_vec[0] = vec[0];
	pre_vec[1] = vec[1];

	double pre_ang = atan2 (pre_vec[1], pre_vec[0]);
	double cur_ang = atan2 (vec[1], vec[0]);

	float rot_sum_tmp = 0.;
	int start_i = i, start_j=j;
	double d_time = T_window/(L-1);
	double half_dt = d_time/2.;

	for (int step=0; step < L; step++)
	{
		// get the vector value at current location (decode the vec_img)
		//vec[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

		//vec[0] = min_vx + (max_vx - min_vx) * vec_img2[i][j][0];
		//vec[1] = min_vy + (max_vy - min_vy) * vec_img2[i][j][1];

		// We need a better strategy to interpolate the vector value at (x, y, t) especially with t being considered now!!

		/*
		   We are now testing the flow generated by formula so that we do not need to worry about the interpolation between two sampled time steps
		   01/28/2013
		*/
		double cur_t = start_time+d_time*step;

		cal_vec_double_gyre (x, y, cur_t, vec[0], vec[1]);
		//cal_vec_ibfv_ex (x, y, cur_t, vec[0], vec[1]);

		if (get_vec_len(vec) < 1.e-6) break;


		// Euler integration
		//x = x+vec[0];
		//y = y+vec[1];
		//// Euler integration
		//normalize_vec(vec);
		//x = x+d_time*vec[0];
		//y = y+d_time*vec[1];

		// we should try to use RK4 for the integration.

		double tx1 = x-half_dt*vec[0];
		double ty1 = y-half_dt*vec[1];
		
		double vec1[2];
		//cal_vec_double_gyre (tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
		cal_vec_double_gyre (tx1, ty1, cur_t, vec1[0], vec1[1]);

		//cal_vec_ibfv_ex (tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
		//normalize_vec(vec1);
		//tx1 = x+half_dt*vec1[0];
		//ty1 = y+half_dt*vec1[1];

		////// RK2
		double total_vec[2];
		total_vec[0] = (vec[0]+vec1[0])/2.;
		total_vec[1] = (vec[1]+vec1[1])/2.;
		x = x-d_time*total_vec[0];
		y = y-d_time*total_vec[1];

		//double vec2[2];
		//cal_vec_double_gyre (tx1, ty1, cur_t+half_dt, vec2[0], vec2[1]);
		//normalize_vec(vec2);
		//tx1 = x+d_time*vec2[0];
		//ty1 = y+d_time*vec2[1];

		//double vec3[2];
		//cal_vec_double_gyre (tx1, ty1, cur_t+d_time, vec3[0], vec3[1]);
		//normalize_vec(vec3);

		////////// RK4
		//double total_vec[2];
		//total_vec[0] = (vec[0]+2.0*vec1[0]+2.0*vec2[0]+vec3[0])/6.;
		//total_vec[1] = (vec[1]+2.0*vec1[1]+2.0*vec2[1]+vec3[1])/6.;
		//x = x+d_time*total_vec[0];
		//y = y+d_time*total_vec[1];


		//i = (int) x*img_res/2;
		//j = (int) y*img_res/2;
		j = (int) x*img_res/2;
		i = (int) y*img_res/2;

		//icVector2 vec1(vec[0], vec[1]), vec2(pre_vec[0], pre_vec[1]);
		//cur_ang = atan2 (vec[1], vec[0]);
		cur_ang = atan2 (total_vec[1], total_vec[0]);

		normalize_vec(vec);
		x = x-vec[0];
		y = y-vec[1];

		j = (int) x*img_res;
		i = (int) y*img_res;

		//icVector2 vec1(vec[0], vec[1]), vec2(pre_vec[0], pre_vec[1]);

		double ang_diff = cur_ang - pre_ang;

		if (ang_diff < -M_PI)
			ang_diff += 2*M_PI;

		else if (ang_diff > M_PI)
			ang_diff -= 2*M_PI;

		pre_ang = cur_ang;

		//total_rotation[i][j] -= ang_diff;
		//rot_sum_tmp +=ang_diff;
		rot_sum_tmp -=ang_diff;

		if (i<0 || i>=img_res || j<0 || j>=img_res)
			return;
	}

	total_rotation[start_i][start_j] += rot_sum_tmp;

	// record the end position of the particle starting at (i, j)
	end_positions_forward[start_i][start_j][0] = x;
	end_positions_forward[start_i][start_j][1] = y;
}



void 
comp_LIC_total_rotation_timeForward(float total_rotation[img_res][img_res], double start_time)
{
	//int L=4000;
	//int L=500;
	int L=1000;
	int i, j;
	for (i = 0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			sum_rotation_forward_timeDep (i, j, L, total_rotation, start_time);
		}
	}	
}


/////////////////////////////
/*
   After obtaining the total rotation for each pixel, we compute the difference with its neighboring pixels along X and Y directions, respectively
   This will give rise to an estimate gradient field
*/
void    
comp_total_rotation_diff(float total_rotation[img_res][img_res],
						float total_rotation_diff[img_res][img_res][2])
{
	// 
	int i, j;

	// Initialization
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
			total_rotation_diff[i][j][0] = total_rotation_diff[i][j][1] = 0.;
	}

	// First compute the difference along Y direction
	for (i=0; i<img_res; i++)
	{
		if (i==0)
		{
			for (j=0; j<img_res; j++)
			{
				//total_rotation_diff[i][j][1] = 2*(total_rotation[i+1][j]-total_rotation[i][j])/(total_rotation[i+1][j]+total_rotation[i][j]);
				if (!fixedPt_pixels[i][j] && !fixedPt_pixels[i+1][j])
				total_rotation_diff[i][j][1] = total_rotation[i+1][j]-total_rotation[i][j];
			}
		}

		else if (i == img_res-1)
		{
			for (j=0; j<img_res; j++)
			{
				//total_rotation_diff[i][j][1] = 2*(total_rotation[i][j]-total_rotation[i-1][j])/(total_rotation[i][j]+total_rotation[i-1][j]);
				if (!fixedPt_pixels[i][j] && !fixedPt_pixels[i-1][j])
				total_rotation_diff[i][j][1] = total_rotation[i][j]-total_rotation[i-1][j];
			}
		}
		else
		{
			for (j=0; j<img_res; j++)
			{
				//total_rotation_diff[i][j][1] = /*0.5**/(total_rotation[i+1][j]-total_rotation[i-1][j])/(total_rotation[i+1][j]+total_rotation[i-1][j]);
				if (!fixedPt_pixels[i-1][j] && !fixedPt_pixels[i+1][j])
				total_rotation_diff[i][j][1] = 0.5*(total_rotation[i+1][j]-total_rotation[i-1][j]);
			}
		}
	}


	// compute the difference along X direction
	for (j=0; j<img_res; j++)
	{
		if (j==0)
		{
			for (i=0; i<img_res; i++)
			{
				//total_rotation_diff[i][j][0] = 2*(total_rotation[i][j+1]-total_rotation[i][j])/(total_rotation[i][j+1]+total_rotation[i][j]);
				if (!fixedPt_pixels[i][j] && !fixedPt_pixels[i][j+1])
				total_rotation_diff[i][j][0] = total_rotation[i][j+1]-total_rotation[i][j];
			}
		}

		else if (j == img_res-1)
		{
			for (i=0; i<img_res; i++)
			{
				//total_rotation_diff[i][j][0] = 2*(total_rotation[i][j]-total_rotation[i][j-1])/(total_rotation[i][j]+total_rotation[i][j-1]);
				if (!fixedPt_pixels[i][j] && !fixedPt_pixels[i][j-1])
				total_rotation_diff[i][j][0] = total_rotation[i][j]-total_rotation[i][j-1];
			}
		}
		else
		{
			for (i=0; i<img_res; i++)
			{
				//total_rotation_diff[i][j][0] = /*0.5**/(total_rotation[i][j+1]-total_rotation[i][j-1])/(total_rotation[i][j+1]+total_rotation[i][j-1]);
				if (!fixedPt_pixels[i][j+1] && !fixedPt_pixels[i][j-1])
				total_rotation_diff[i][j][0] = 0.5*(total_rotation[i][j+1]-total_rotation[i][j-1]);
			}
		}
	}
}



void    
get_color_map_for_rot_sum_diff()
{
	int i, j;
	float max_len=-1.e8, min_len=1.e8;

	// initialization
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			total_rotation_diff_color[i][j][0]=
			total_rotation_diff_color[i][j][1]=
			total_rotation_diff_color[i][j][2]=0.;
		}
	}

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{	
			//double len = get_vec_len(total_rotation_diff[i][j]);

			//// Note the following directional gradient computation will not work in 3D!!!!!
			icVector2 tmp_vec;
			tmp_vec.entry[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
			tmp_vec.entry[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

			normalize(tmp_vec);
			icVector2 proj_vec;
			proj_vec.entry[0] = -tmp_vec.entry[1];
			proj_vec.entry[1] = tmp_vec.entry[0];

			icVector2 grad_vec(total_rotation_diff[i][j][0], total_rotation_diff[i][j][1]);
			double len = fabs(dot(grad_vec, proj_vec));
				
			total_rotation_diff_mag[i][j] = len/*log(len)/2.*/;  // 

			//if (len > max_len) max_len = len;

			//if (len < min_len) min_len = len;
			
			if (total_rotation_diff_mag[i][j] > max_len) max_len = total_rotation_diff_mag[i][j];
			if (total_rotation_diff_mag[i][j] < min_len) min_len = total_rotation_diff_mag[i][j];
		}
	}


	// convert to color map
	float hsv[3]={0.}, rgb[3]={0.};
	float range = max_len - min_len;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			//double len = get_vec_len(total_rotation_diff[i][j]);
			
			hsv[2] = 1.;

			hsv[1] = 1. - 2.*(total_rotation_diff_mag[i][j] - min_len)/range;

			if (hsv[1] < 0)
			{
				hsv[0] = 0.;
				hsv[1] = abs(hsv[1]);
				//hsv[1] = 0.;
				hsv[2] = 0.;
			}
			else{
				hsv[0] = 240.;
				//hsv[2] = 1.;
				hsv[2] = abs(hsv[1]);
				hsv[1] =0.;
			}

			//if(total_rotation[i][j]>max_rot/4. || total_rotation[i][j]<min_rot/4.)
			//	hsv[1] = 1;

			HsvRgb (hsv, rgb);
			total_rotation_diff_color[i][j][0] = (unsigned char)(rgb[0]*255.);
			total_rotation_diff_color[i][j][1] = (unsigned char)(rgb[1]*255.);
			total_rotation_diff_color[i][j][2] = (unsigned char)(rgb[2]*255.);
		}
	}

	// combine with LIC

	// Probably do not combined with the LIC at this moment?
/*
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			total_rotation_diff_color[i][j][0] = 0.7*total_rotation_diff_color[i][j][0] + 0.3 * LIC_tex[i][j][0];
			total_rotation_diff_color[i][j][1] = 0.7*total_rotation_diff_color[i][j][1] + 0.3 * LIC_tex[i][j][1];
			total_rotation_diff_color[i][j][2] = 0.7*total_rotation_diff_color[i][j][2] + 0.3 * LIC_tex[i][j][2];
		}
	}
	*/
}


///////////////////////////////////////////////////////////////////////////////
void    
comp_total_rotation_diff_Hessian()
{
	// 
	int i, j;

	//// We first need to compute the first derivative of the difference magnitude

	// First compute the difference along Y direction
	for (i=0; i<img_res; i++)
	{
		if (i==0)
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff[i][j][1] = total_rotation_diff_mag[i+1][j]-total_rotation_diff_mag[i][j];
			}
		}

		else if (i == img_res-1)
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff[i][j][1] = total_rotation_diff_mag[i][j]-total_rotation_diff_mag[i-1][j];
			}
		}
		else
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff[i][j][1] = 0.5*(total_rotation_diff_mag[i+1][j]-total_rotation_diff_mag[i-1][j]);
			}
		}
	}


	// compute the difference along X direction
	for (j=0; j<img_res; j++)
	{
		if (j==0)
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff[i][j][0] = total_rotation_diff_mag[i][j+1]-total_rotation_diff_mag[i][j];
			}
		}

		else if (j == img_res-1)
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff[i][j][0] = total_rotation_diff_mag[i][j]-total_rotation_diff_mag[i][j-1];
			}
		}
		else
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff[i][j][0] = 0.5*(total_rotation_diff_mag[i][j+1]-total_rotation_diff_mag[i][j-1]);
			}
		}
	}


	//// Next, we need to compute the second derivative, which is the Hessian
	// compute the difference along X direction for X component
	for (j=0; j<img_res; j++)
	{
		if (j==0)
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff_Hessian[i][j][0] = total_rotation_diff[i][j+1][0]-total_rotation_diff[i][j][0];
			}
		}

		else if (j == img_res-1)
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff_Hessian[i][j][0] = total_rotation_diff[i][j][0]-total_rotation_diff[i][j-1][0];
			}
		}
		else
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff_Hessian[i][j][0] = 0.5*(total_rotation_diff[i][j+1][0]-total_rotation_diff[i][j-1][0]);
			}
		}
	}

	//compute the difference along Y direction for X component
	for (i=0; i<img_res; i++)
	{
		if (i==0)
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff_Hessian[i][j][1] = total_rotation_diff[i+1][j][0]-total_rotation_diff[i][j][0];
			}
		}

		else if (i == img_res-1)
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff_Hessian[i][j][1] = total_rotation_diff[i][j][0]-total_rotation_diff[i-1][j][0];
			}
		}
		else
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff_Hessian[i][j][1] = 0.5*(total_rotation_diff[i+1][j][0]-total_rotation_diff[i-1][j][0]);
			}
		}
	}

	// compute the difference along X direction for Y component
	for (j=0; j<img_res; j++)
	{
		if (j==0)
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff_Hessian[i][j][2] = total_rotation_diff[i][j+1][1]-total_rotation_diff[i][j][1];
			}
		}

		else if (j == img_res-1)
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff_Hessian[i][j][2] = total_rotation_diff[i][j][1]-total_rotation_diff[i][j-1][1];
			}
		}
		else
		{
			for (i=0; i<img_res; i++)
			{
				total_rotation_diff_Hessian[i][j][2] = 0.5*(total_rotation_diff[i][j+1][1]-total_rotation_diff[i][j-1][1]);
			}
		}
	}

	//compute the difference along Y direction for Y component
	for (i=0; i<img_res; i++)
	{
		if (i==0)
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff_Hessian[i][j][3] = total_rotation_diff[i+1][j][1]-total_rotation_diff[i][j][1];
			}
		}

		else if (i == img_res-1)
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff_Hessian[i][j][3] = total_rotation_diff[i][j][1]-total_rotation_diff[i-1][j][1];
			}
		}
		else
		{
			for (j=0; j<img_res; j++)
			{
				total_rotation_diff_Hessian[i][j][3] = 0.5*(total_rotation_diff[i+1][j][1]-total_rotation_diff[i-1][j][1]);
			}
		}
	}
}

float CalDeterminant(float mat[4])
{
	return (mat[0]*mat[3]-mat[1]*mat[2]);
}

void GetEigenVectors(float mat[4], float evalues[2], float ev[2][2])
{
	//first calculate the dominant of the matrix, if it is zero, return
	if(abs(CalDeterminant(mat)) < 1e-30)
		return;

	//ev[0][0] = mat[3] - evalues[0];
	//ev[0][1] = - mat[2];
	//
	//ev[1][0] = mat[3] - evalues[1];
	//ev[1][1] = - mat[2];

	ev[0][0] = -mat[1];
	ev[0][1] = mat[0] - evalues[0];

	ev[1][0] = -mat[1];
	ev[1][1] = mat[0] - evalues[1];

}

void CalEigenForMatrix2x2(float mat[4], float evec[2][2], float ev[2])
{
	double la, lb, lc, ld, A, B, C, delta;

	la = mat[0];
	lb = mat[1];
	lc = mat[2];
	ld = mat[3];

	A = 1;
	B = -(la + ld);
	C = (la * ld - lb * lc);                

	delta =  B*B - 4*A*C;

	if(delta >= 0)
	{
		ev[0] = (-B + sqrt(delta))/2;
		ev[1] = (-B - sqrt(delta))/2;

		////find the largest (absolute) eigen values, and store it to the first element
		//if(abs(evalues[0]) < abs(evalues[1]))
		//{
			//double temp = ev[0];
			//ev[0] = ev[1];
			//ev[1] = temp;
		//}

		//for real eigen values, we calculate the eigen vectors of it
		GetEigenVectors(mat, ev, evec);

		normalize_vec(evec[0]);
		normalize_vec(evec[1]);
	}

	else
	{
		evec[0][0]=evec[0][1] = 0;
		evec[1][0]=evec[1][1] = 0;
	}
}



void    
get_color_map_for_rot_sum_diff_Hessian()
{
	int i, j;
	float max_len=-1.e8, min_len=1.e8;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{	
			double len = get_vec_len(total_rotation_diff[i][j]);

			//// Note the following directional gradient computation will not work in 3D!!!!!
			
			//icVector2 tmp_vec;
			//tmp_vec.entry[0] = min_vx + (max_vx - min_vx) * vec_img[i][j][0]/255;
			//tmp_vec.entry[1] = min_vy + (max_vy - min_vy) * vec_img[i][j][1]/255;

			//normalize(tmp_vec);
			//icVector2 proj_vec;
			//proj_vec.entry[0] = -tmp_vec.entry[1];
			//proj_vec.entry[1] = tmp_vec.entry[0];

			float evec[2][2];
			float ev[2];
			float mat[4] = {total_rotation_diff_Hessian[i][j][0],
				total_rotation_diff_Hessian[i][j][1],
				total_rotation_diff_Hessian[i][j][2],
				total_rotation_diff_Hessian[i][j][3]};
			CalEigenForMatrix2x2(mat, evec, ev);

			//total_rotation_diff[i][j][0] = evec[0][0];
			//total_rotation_diff[i][j][1] = evec[0][1];

			icVector2 grad_vec(total_rotation_diff[i][j][0], total_rotation_diff[i][j][1]);
			icVector2 maj_evec(evec[0][0], evec[0][1]);
			normalize(grad_vec);
			normalize(maj_evec);
			float len1 = fabs(dot(grad_vec, maj_evec));

			//normalize_vec(evec[0]);
			//double len = fabs(tmp_vec.entry[0]*evec[0][0] + tmp_vec.entry[1]*evec[0][1]);
		
			//if (ev[0] < 0 && len < 1.e-4) len = 1;
			//else
			//	len = 0;
				
			total_rotation_diff_mag[i][j] = len/*log(len)/2.*/;  // 

			//if (len > max_len) max_len = len;

			//if (len < min_len) min_len = len;
			
			if (total_rotation_diff_mag[i][j] > max_len) max_len = total_rotation_diff_mag[i][j];
			if (total_rotation_diff_mag[i][j] < min_len) min_len = total_rotation_diff_mag[i][j];
		}
	}


	// convert to color map
	float hsv[3]={0.}, rgb[3]={0.};
	float range = max_len - min_len;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			//double len = get_vec_len(total_rotation_diff[i][j]);
			
			hsv[2] = 1.;

			hsv[1] = 1. - 2.*(total_rotation_diff_mag[i][j] - min_len)/range;

			if (hsv[1] < 0)
			{
				hsv[0] = 0.;
				hsv[1] = abs(hsv[1]);
			}
			else
				hsv[0] = 240.;

			//if(total_rotation[i][j]>max_rot/4. || total_rotation[i][j]<min_rot/4.)
			//	hsv[1] = 1;

			HsvRgb (hsv, rgb);
			total_rotation_diff_color[i][j][0] = (unsigned char)(rgb[0]*255.);
			total_rotation_diff_color[i][j][1] = (unsigned char)(rgb[1]*255.);
			total_rotation_diff_color[i][j][2] = (unsigned char)(rgb[2]*255.);
		}
	}

	// combine with LIC

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			total_rotation_diff_color[i][j][0] = 0.7*total_rotation_diff_color[i][j][0] + 0.3 * LIC_tex[i][j][0];
			total_rotation_diff_color[i][j][1] = 0.7*total_rotation_diff_color[i][j][1] + 0.3 * LIC_tex[i][j][1];
			total_rotation_diff_color[i][j][2] = 0.7*total_rotation_diff_color[i][j][2] + 0.3 * LIC_tex[i][j][2];
		}
	}
}


///////////////////////////////////////////////////////////////////////////////

void 
cal_SpatialGradient()
{
	int i, j;

	// initialization

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			spatial_gradient_forward[i][j][0]=
			spatial_gradient_forward[i][j][1]=
			spatial_gradient_forward[i][j][2]=
			spatial_gradient_forward[i][j][3]= 0.;
		}
	}

	// a temporary array to store the world coordinates of the starting positions of those particles which are at the centers of pixels
	//float start_pos[img_res][img_res][2];

	//for (i=0; i<img_res; i++)
	//{
	//	for (j=0; j<img_res; j++)
	//	{
	//		start_pos[i][j][0] = 2*i
	//	}
	//}

	// compute the gradient
	for (i=1; i<img_res-1; i++)
	{
		for (j=1; j<img_res-1; j++)
		{
			//spatial_gradient_forward[i][j][0] = (end_positions_forward[i+1][j][0] - end_positions_forward[i-1][j][0]) / (start_pos[i+1][j][0] - start_pos[i-1][j][0]);
			//spatial_gradient_forward[i][j][1] = (end_positions_forward[i][j+1][0] - end_positions_forward[i][j-1][0]) / (start_pos[i][j+1][1] - start_pos[i][j-1][1]);
			//spatial_gradient_forward[i][j][2] = (end_positions_forward[i+1][j][1] - end_positions_forward[i-1][j][1]) / (start_pos[i+1][j][0] - start_pos[i-1][j][0]);
			//spatial_gradient_forward[i][j][3] = (end_positions_forward[i][j+1][1] - end_positions_forward[i][j-1][1]) / (start_pos[i][j+1][1] - start_pos[i][j-1][1]);

			spatial_gradient_forward[i][j][0] = (end_positions_forward[i][j+1][0] - end_positions_forward[i][j-1][0]) / (start_pos[i][j+1][0] - start_pos[i][j-1][0]);
			spatial_gradient_forward[i][j][1] = (end_positions_forward[i+1][j][0] - end_positions_forward[i-1][j][0]) / (start_pos[i+1][j][1] - start_pos[i-1][j][1]);
			spatial_gradient_forward[i][j][2] = (end_positions_forward[i][j+1][1] - end_positions_forward[i][j-1][1]) / (start_pos[i][j+1][0] - start_pos[i][j-1][0]);
			spatial_gradient_forward[i][j][3] = (end_positions_forward[i+1][j][1] - end_positions_forward[i-1][j][1]) / (start_pos[i+1][j][1] - start_pos[i-1][j][1]);
		}
	}

	// handle the gradient at the boundary 
}


float** 
cal_CauchyGreen(float inputMatrix[][2])
{
	float** CauchyGreenDistortion;
	CauchyGreenDistortion = new float*[2];	
	CauchyGreenDistortion[0] = new float[2];
	CauchyGreenDistortion[1] = new float[2];

	//Phi' * Phi
	CauchyGreenDistortion[0][0] = inputMatrix[0][0]*inputMatrix[0][0]+inputMatrix[1][0]*inputMatrix[1][0];
	CauchyGreenDistortion[0][1] = inputMatrix[0][0]*inputMatrix[0][1]+inputMatrix[1][0]*inputMatrix[1][1];
	CauchyGreenDistortion[1][0] = inputMatrix[0][1]*inputMatrix[0][0]+inputMatrix[1][1]*inputMatrix[1][0];
	CauchyGreenDistortion[1][1] = inputMatrix[0][1]*inputMatrix[0][1]+inputMatrix[1][1]*inputMatrix[1][1];

	return CauchyGreenDistortion;

}

// Return the largest eigen value
float 
cal_EigenVector2x2(float** inputMatrix)
{
	float a = 1.;
	float b = -(inputMatrix[0][0] + inputMatrix[1][1]);
	float c = inputMatrix[0][0]*inputMatrix[1][1] - inputMatrix[0][1]*inputMatrix[1][0];

	float det = pow(b,2) - 4*a*c;
	float lambda1, lambda2;
	if(det >= 0){
		lambda1 = (-1.*b + sqrt(det))/(2.*a);
		lambda2 = (-1.*b - sqrt(det))/(2.*a);
	}
	else{
		return 0.0f;
		//return (-1.*b)/(2.*a);
	}
	if (lambda1>lambda2) return lambda1;
	
	return lambda2;

}


void
cal_FTLE_forward()
{
	int i, j;
	float M[2][2];
	float largest_eigenvalue;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			// first compute the Cauchy Green Tensor
			M[0][0] = spatial_gradient_forward[i][j][0];
			M[0][1] = spatial_gradient_forward[i][j][1];
			M[1][0] = spatial_gradient_forward[i][j][2];
			M[1][1] = spatial_gradient_forward[i][j][3];
			float** CachyGreenDistrotionTensor= cal_CauchyGreen(M);

			largest_eigenvalue = cal_EigenVector2x2(CachyGreenDistrotionTensor);

			if (largest_eigenvalue > 0.)
				FTLE_forward[i][j] = 10. * log (largest_eigenvalue)/T_window;
			else
				FTLE_forward[i][j] = 0.;
		}
	}
}


void 
get_FTLE_forward_color()
{
	int i, j;

	float max_=-1.e8, min_=1.e8;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			if (FTLE_forward[i][j] < min_) min_ = FTLE_forward[i][j];

			if (FTLE_forward[i][j] > max_) max_ = FTLE_forward[i][j];
		}
	}

	// compute the color using a simple rainbow color scheme at this moment
	float hsv[3] = {0, 1, 1};
	float rgb[3];
	float rang = max_ - min_;
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			hsv[0] = 240. - 240 * (FTLE_forward[i][j] - min_)/rang;

			HsvRgb(hsv, rgb);

			FTLE_forward_color[i][j][0] = 255 * rgb[0];
			FTLE_forward_color[i][j][1] = 255 * rgb[1];
			FTLE_forward_color[i][j][2] = 255 * rgb[2];
		}
	}
}

void 
compute_FTLE_forward()
{
	cal_SpatialGradient();
	cal_FTLE_forward();
	get_FTLE_forward_color();
}

/*  need to handle boundary

void 
cal_SpatialGradients(){
	for (int t = 1; t < INTEGRATIONTIME*10; t++){
		for (int i = 1; i < GRIDSIZEX - 1; i++){
			for (int j = 1; j < GRIDSIZEY - 1; j++){

				SpatailGradientFlowX_dx[t][i][j] = (flowMapX[t][i+1][j] - flowMapX[t][i-1][j])/(flowMapX[0][i+1][j] - flowMapX[0][i-1][j]);
				SpatailGradientFlowX_dy[t][i][j] = (flowMapX[t][i][j+1] - flowMapX[t][i][j-1])/(flowMapY[0][i][j+1] - flowMapY[0][i][j-1]);
				SpatailGradientFlowY_dx[t][i][j] = (flowMapY[t][i+1][j] - flowMapY[t][i-1][j])/(flowMapX[0][i+1][j] - flowMapX[0][i-1][j]);
				SpatailGradientFlowY_dy[t][i][j] = (flowMapY[t][i][j+1] - flowMapY[t][i][j-1])/(flowMapY[0][i][j+1] - flowMapY[0][i][j-1]);							
				//printf("%f, %f\n", SpatailGradientFlowX_dy[t][i][j], SpatailGradientFlowY_dx[t][i][j]);
				
			}
		}
	}
}

void 
cal_FTLEField()
{

	float Phi[2][2];float** CachyGreenDistrotionTensor;
	float eigenVector;
	for (int t = 1; t < INTEGRATIONTIME*10; t++){
		for (int i = 1; i < GRIDSIZEX - 1; i++){
			for (int j = 1; j < GRIDSIZEY - 1; j++){
				Phi[0][0] = SpatailGradientFlowX_dy[t][i][j];								
				Phi[0][1] = SpatailGradientFlowX_dx[t][i][j];								
				Phi[1][0] = SpatailGradientFlowY_dy[t][i][j];								
				Phi[1][1] = SpatailGradientFlowY_dx[t][i][j];		
				CachyGreenDistrotionTensor = CalculateCauchyGreen(Phi);

				eigenVector = CalcualteEigenVector(CachyGreenDistrotionTensor);
				
				if(eigenVector > 0.){
					FTLEfield[t][i][j] = 10.0f * log(eigenVector)/(float)t;
				}
				else{
					FTLEfield[t][i][j] = 0.;
				}
			}
		}
	}
}

*/

////////////////////////////////////////////////////////////////////////////////////////////
// File output section

void 
write_ppm(char *filename, unsigned char *img, int dimx, int dimy)
{
FILE *fp = fopen(filename, "wb"); /* b - binary mode */
  (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
  int i, j;
  for (j = 0; j < dimy; ++j)
  {
    for (i = 0; i < dimx; ++i)
    {
      static unsigned char color[3];
      color[0] = img[3*(j*dimx+i)];  /* red */
      color[1] = img[3*(j*dimx+i)+1];  /* green */
      color[2] = img[3*(j*dimx+i)+2];  /* blue */
      (void) fwrite(color, 1, 3, fp);
    }
  }
  (void) fclose(fp);
}


void 
write_ppm_flippedY(char *filename, unsigned char *img, int dimx, int dimy)
{
FILE *fp = fopen(filename, "wb"); /* b - binary mode */
  (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
  int i, j;
  for (j = dimy-1; j >=0; j--)
  {
    for (i = 0; i < dimx; ++i)
    {
      static unsigned char color[3];
      color[0] = img[3*(j*dimx+i)];  /* red */
      color[1] = img[3*(j*dimx+i)+1];  /* green */
      color[2] = img[3*(j*dimx+i)+2];  /* blue */
      (void) fwrite(color, 1, 3, fp);
    }
  }
  (void) fclose(fp);
}


void    
init_par_pos()
{
	int i, j;

	//float p_interval = 1./(PDIM-1);
	//float x_interval = 4./(PDIM-1);
	float p_interval = 1.5/(PDIM-1);
	float x_interval = 1.5/(PDIM-1);

	for (i=0; i<PDIM; i++)
	{
		//float py = -0.5+i*p_interval;
		float py = -0.75+i*p_interval;
		for (j=0; j<PDIM; j++)
		{
			//float px = -2.+j*x_interval/*p_interval*/;
			float px = -0.75+j*x_interval/*p_interval*/;

			particle_pos[i][j][0] = px;
			particle_pos[i][j][1] = py;
		}
	}
}


void   
advect_par(double &cur_t, double dt)
{
	int i, j;

	double half_dt = dt/2.;

	double vec[2], vec1[2];

	for (i=0; i<PDIM; i++)
	{
		for (j=0; j<PDIM; j++)
		{
			double x, y;
			x = particle_pos[i][j][0];
			y = particle_pos[i][j][1];

			//cal_vec_double_gyre (x, y, cur_t, vec[0], vec[1]);
			//cal_vec_ibfv_ex(x, y, cur_t, vec[0], vec[1]);
			//cal_gyre_saddle(x, y, cur_t, vec[0], vec[1]);
			//cal_stuart_vortices(x, y, cur_t, vec[0], vec[1]);
			cal_Oseen_vortices(x, y, cur_t, vec[0], vec[1]);

			if (get_vec_len(vec) < 1.e-6) continue;


			normalize_vec(vec);

			double tx1 = x+half_dt*vec[0];
			double ty1 = y+half_dt*vec[1];
		
			//cal_vec_double_gyre (tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
			//cal_vec_ibfv_ex (tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
			//cal_gyre_saddle(tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
			//cal_stuart_vortices(tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
			cal_Oseen_vortices(tx1, ty1, cur_t+half_dt, vec1[0], vec1[1]);
			
			normalize_vec(vec1);
			tx1 = x+half_dt*vec1[0];
			ty1 = y+half_dt*vec1[1];

			////// RK2
			double total_vec[2];
			total_vec[0] = (vec[0]+vec1[0])/2.;
			total_vec[1] = (vec[1]+vec1[1])/2.;
			normalize_vec(total_vec);
			particle_pos[i][j][0] = x+dt*total_vec[0];
			particle_pos[i][j][1] = y+dt*total_vec[1];
		}
	}
}

///////////////////////////////////////////////////////////////////////

void   
cal_vec_ibfv_ex(double x, double y, double t, double &vx, double &vy)
{
   double dx, dy, r;
   double sa;

   //sa = 0.001 * cos (t*2.*PI);
   //sa = 0.01 * cos (t*1.*PI);
   sa = 0.001 * cos (t*1.*PI);

   //dx = x - 0.5;         
   //dy = y - 0.5; 
   //dx = x - 0.9;         
   //dy = y - 0.7; 
   dx = x - 1.;         
   dy = y - 1.; 
   r  = dx*dx + dy*dy; 
   if (r < 1.e-6) r = 1.e-6;
   //vx = sa*dx/r + 0.0002 * cos (t*1.*PI);  
   ////vx = sa*dx/r + 0.0002 * sin (t*1.*PI);  
   ////vy = sa*dy/r - 0.0002 * sin (t*10.*PI);
   ////vx = sa*dx/r + 0.0002 * cos (t*2.*PI);  
   //vy = sa*dy/r - 0.002 * sin (t*2.*PI);
   ////vy = sa*dy/r + 0.0002 * sin (t*2.*PI);
  
   vx = sa*dx/r + 0.002 * dy* cos (t*0.5*PI);  
   //vy = sa*dy/r - 0.002 * dx* cos (t*.4*PI);
   vy = sa*dy/r - 0.002 * dx* sin (t*.4*PI);
}


void    
gen_ibfv_at_time(double time)
{
	// First, create a small mesh say 64x64 mesh, then compute the vector values at the grid points of the mesh, render it to the vector image
	int i, j;
	double x1, x2, y, vx, vy;
	double dx = 4./63;
	float rgb[3]={0.};
	min_vx = 1.e8, max_vx=-1.e8, min_vy = 1.e8, max_vy=-1.e8;

	// first, find out the maximum and minimum vector value
	for (i=0; i<64; i++)
	{
		x1 = -2.+i*dx;

		for (j=0; j<64; j++)
		{
			y = -2.+j*dx;
			cal_vec_ibfv_ex(x1, y, time, vx, vy);

			if (vx < min_vx) min_vx = vx;
			if (vx > max_vx) max_vx = vx;
			if (vy < min_vy) min_vy = vy;
			if (vy > max_vy) max_vy = vy;
		}
	}

	double x_rang = max_vx - min_vx;
	double y_rang = max_vy - min_vy;

	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 2, 0, 2);
	glClear(GL_COLOR_BUFFER_BIT);

	glDrawBuffer(GL_BACK);

	for (i=0; i<63; i++)
	{
		x1 = -2.+i*dx;
		x2 = x1 + dx;

		glBegin(GL_QUAD_STRIP);
		for (j=0; j<64; j++)
		{
			y = -2.+j*dx;
			cal_vec_ibfv_ex(x1, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x1, y);

			cal_vec_ibfv_ex(x2, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x2, y);

		}
		glEnd();
	}
	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, vec_img);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_FLOAT, vec_img2);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
}


double 
c_fun(double x)
{
	if ( x < -PI/2. || x > PI/2.)
		return 0;
	else
		return cos(x);
}

void    
cal_gyre_saddle(double x, double y, double t, double &vx, double &vy)
{
	double a_t = sin(2*PI*t)/3.;

	double k, l;

	k = (y >= fabs(x)? 1:0) - (y <= -fabs(x)? 1:0);
	l = (x >= fabs(y)? 1:0) - (x <= -fabs(y)? 1:0);

	if (x >= -.5 && x <= .5 && y >= -.5 && y <= .5)
	{
		vx = -sin(PI*x)*cos(PI*y) + a_t*sin(PI*y)*cos(PI*x);
		vy = sin(PI*y)*cos(PI*x) - a_t*sin(PI*x)*cos(PI*y);
	}
	else
	{
		vx = k*a_t*c_fun(PI*x - a_t*PI*(y-k/2.)) - l*c_fun(PI*y - a_t*PI*(x-l/2.));
		vy = k*c_fun(PI*x-a_t*PI*(y-k/2.)) - l*a_t*c_fun(PI*y-a_t*PI*(x-l/2.));
	}
}

void    
gen_gyre_saddle_at_time(double time)
{
	// First, create a small mesh say 64x64 mesh, then compute the vector values at the grid points of the mesh, render it to the vector image
	int i, j;
	double x1, x2, y, vx, vy;
	double dx = 2./63;
	float rgb[3]={0.};
	min_vx = 1.e8, max_vx=-1.e8, min_vy = 1.e8, max_vy=-1.e8;

	// first, find out the maximum and minimum vector value
	for (i=0; i<64; i++)
	{
		x1 = -1.+i*dx;

		for (j=0; j<64; j++)
		{
			y = -1.+j*dx;
			cal_gyre_saddle(x1, y, time, vx, vy);

			if (vx < min_vx) min_vx = vx;
			if (vx > max_vx) max_vx = vx;
			if (vy < min_vy) min_vy = vy;
			if (vy > max_vy) max_vy = vy;
		}
	}

	double x_rang = max_vx - min_vx;
	double y_rang = max_vy - min_vy;

	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(-1., 1., -1., 1.);
	glClear(GL_COLOR_BUFFER_BIT);

	glDrawBuffer(GL_BACK);

	for (i=0; i<63; i++)
	{
		x1 = -1.+i*dx;
		x2 = x1 + dx;

		glBegin(GL_QUAD_STRIP);
		for (j=0; j<64; j++)
		{
			y = -1.+j*dx;
			cal_gyre_saddle(x1, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x1, y);

			cal_gyre_saddle(x2, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x2, y);

		}
		glEnd();
	}

	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, vec_img);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_FLOAT, vec_img2);

	// compute the vector field magnitude value now
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			float vec[2] = {vec_img2[i][j][0], vec_img2[i][j][1]};
			vec_mag[i][j] = get_vec_len(vec);

			fixedPt_pixels[i][j] = false;
		}
	}

	// initialize the vec_img and vec_img2
	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			//vec_img[i][j][0] = vec_img[i][j][1] = vec_img[i][j][2] = 0;
			//vec_img2[i][j][0] = vec_img2[i][j][1] = vec_img2[i][j][2] = 0.;
			//double vx = vec_img[i][j][0]/255.*(max_vx-min_vx) + min_vx;
			//double vy = vec_img[i][j][1]/255.*(max_vy-min_vy) + min_vy;
			icVector2 vec(vec_img2[i][j][0], vec_img2[i][j][1]);
			double vf_mag = length(vec);

			if (vf_mag < 1.e-6)
			{
				vec_img[i][j][0] = vec_img[i][j][1] = vec_img[i][j][2] = 0;
				vec_img2[i][j][0] = vec_img2[i][j][1] = vec_img2[i][j][2] = 0.;
				fixedPt_pixels[i][j] = true;
			}
		}
	}

	write_ppm("render_vec_img.ppm", (unsigned char*)vec_img, img_res, img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
}


void    
cal_stuart_vortices(double x, double y, double t, double &vx, double &vy)
{
	double dominator = cosh(2*y)-.25*cos(2*(x-t));
	vx = sinh(2*y)/dominator + 1;
	vy = -.25*sin(2*(x-t))/dominator;
}

void
gen_stuart_vortices_at_time(double time)
{
	// First, create a small mesh say 64x64 mesh, then compute the vector values at the grid points of the mesh, render it to the vector image
	int i, j;
	double x1, x2, y, vx, vy;
	double dx = 2./63;
	float rgb[3]={0.};
	min_vx = 1.e8, max_vx=-1.e8, min_vy = 1.e8, max_vy=-1.e8;

	// first, find out the maximum and minimum vector value
	for (i=0; i<64; i++)
	{
		x1 = -1.+i*dx;

		for (j=0; j<64; j++)
		{
			y = -1.+j*dx;
			cal_stuart_vortices(x1, y, time, vx, vy);

			if (vx < min_vx) min_vx = vx;
			if (vx > max_vx) max_vx = vx;
			if (vy < min_vy) min_vy = vy;
			if (vy > max_vy) max_vy = vy;
		}
	}

	double x_rang = max_vx - min_vx;
	double y_rang = max_vy - min_vy;

	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(-1., 1., -1., 1.);
	glClear(GL_COLOR_BUFFER_BIT);

	glDrawBuffer(GL_BACK);

	for (i=0; i<63; i++)
	{
		x1 = -1.+i*dx;
		x2 = x1 + dx;

		glBegin(GL_QUAD_STRIP);
		for (j=0; j<64; j++)
		{
			y = -1.+j*dx;
			cal_stuart_vortices(x1, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x1, y);

			cal_stuart_vortices(x2, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x2, y);

		}
		glEnd();
	}
	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, vec_img);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_FLOAT, vec_img2);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
}


double Oseen_vortex_centers[3][3]={ {-0.4, -0.2, 0.1},
{-0.1, 0.5, 0.4}, 
{0.5,0, 0.3}
}; 

void    
cal_Oseen_vortices(double x, double y, double t, double &vx, double &vy)
{
	// Determine the centers of these vertices
	int i;
	double rc;
	double alpha = 1.25643;

	double r, rsquare;

	vx = vy = 0;

	for (i=0; i<3; i++)
	{
		double dx, dy;
		
		if (i == 0)
		{
			dx = x - (Oseen_vortex_centers[i][0]+cos(2*PI*t));
			dy = y - (Oseen_vortex_centers[i][1]+sin(2*PI*t));
			rc = sqrt(4*1.1*t);
		}

		else if (i==1)
		{
			dx = x - (Oseen_vortex_centers[i][0]+cos(2.1*PI*t));
			dy = y - (Oseen_vortex_centers[i][1]+sin(1.6*PI*t));
			rc = sqrt(4*1.*t);
		}
		else
		{
			dx = x - (Oseen_vortex_centers[i][0]-cos(PI*t));
			dy = y - (Oseen_vortex_centers[i][1]-sin(PI*t));
			rc = sqrt(4*1.9*t);
		}

		rsquare = dx*dx + dy*dy;
		r = sqrt (rsquare);

		double ang = atan2(dy, dx);

		double coeff = (1+.5/alpha)*rc/r*(1-exp(-alpha*rsquare/(rc*rc)));

		double theta;

		if (i==0)
		{
			theta = 2*t;
			coeff = (0.5+cos(ang-theta)) * coeff;
		}
		else if (i==1)
		{
			theta = 5*t;
			coeff = (0.8+cos(ang-theta)) * 0.8 * coeff;
		}
		else
		{
			theta = -10*t;
			coeff = (0.3+cos(ang-theta)) * coeff;
		}

		vx += -(coeff*dy);
		vy += (coeff*dx);
	}

	vx /= 3;
	vy /= 3;
}


void    
gen_Oseen_vortices_at_time(double time)
{
	// First, create a small mesh say 64x64 mesh, then compute the vector values at the grid points of the mesh, render it to the vector image
	int i, j;
	double x1, x2, y, vx, vy;
	double dx = 2./63;
	float rgb[3]={0.};
	min_vx = 1.e8, max_vx=-1.e8, min_vy = 1.e8, max_vy=-1.e8;

	// first, find out the maximum and minimum vector value
	for (i=0; i<64; i++)
	{
		x1 = -1.+i*dx;

		for (j=0; j<64; j++)
		{
			y = -1.+j*dx;
			cal_Oseen_vortices(x1, y, time, vx, vy);

			if (vx < min_vx) min_vx = vx;
			if (vx > max_vx) max_vx = vx;
			if (vy < min_vy) min_vy = vy;
			if (vy > max_vy) max_vy = vy;
		}
	}

	double x_rang = max_vx - min_vx;
	double y_rang = max_vy - min_vy;

	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(-1., 1., -1., 1.);
	glClear(GL_COLOR_BUFFER_BIT);

	glDrawBuffer(GL_BACK);

	for (i=0; i<63; i++)
	{
		x1 = -1.+i*dx;
		x2 = x1 + dx;

		glBegin(GL_QUAD_STRIP);
		for (j=0; j<64; j++)
		{
			y = -1.+j*dx;
			cal_Oseen_vortices(x1, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x1, y);

			cal_Oseen_vortices(x2, y, time, vx, vy);
			// get the color for the current vertex
			rgb[1] =  (vx - min_vx)/x_rang;
			rgb[0] =  (vy - min_vy)/y_rang;
			glColor3fv (rgb);
			glVertex2f (x2, y);

		}
		glEnd();
	}
	// save the rendered image into the vec_img
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, vec_img);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_FLOAT, vec_img2);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
}


void 
init_particle_list()
{
	particles.clear();
	pixels_to_particles.clear();
	pixels_to_particles.resize(img_res*img_res);
	particle_active.clear();
}

void    
insert_new_particles()
{
	int i, j;

	double interval = 2./(img_res-1);

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			// add a new particle
			std::pair<float, float> par;
			par.first = 2*(j+.5)/img_res;
			par.second = 2*(i+.5)/img_res;

			particles.push_back(par);
			particle_active.push_back(true);

			// add to the pixels_to_particles pointers
			pixels_to_particles[i*img_res+j].push_back(particles.size()-1);
		}
	}
}

// Update the positions of the particles. Move one step forward for each particle
void   
update_all_particles_streaklines (double cur_t, double dt)
{
	int i;

	for (i=0; i<particles.size(); i++)
	{
		if (!particle_active[i]) continue;

		double vec[2];
		std::pair<float, float> &par = particles[i];
		cal_vec_double_gyre ((double)par.first, (double)par.second, cur_t, vec[0], vec[1]);
		if (get_vec_len(vec) < 1.e-6) { particle_active[i]=false; continue;}

		normalize_vec(vec);

		// to verify the framework of streakline computation, let us use Euler integration at this moment
		par.first += dt*vec[0];
		par.second += dt*vec[1];
		
		//double tx = par.first+dt/2.*vec[0];
		//double ty = par.second+dt/2.*vec[1];
		//
		//double vec1[2];
		//cal_vec_double_gyre (tx, ty, cur_t+dt/2., vec1[0], vec1[1]);
		//normalize_vec(vec1);

		//par.first = par.first + 0.5*dt*(vec[0]+vec1[0]);
		//par.second = par.second + 0.5*dt*(vec[1]+vec1[1]);


		int row = (int) par.first*img_res/2;
		int col = (int) par.second*img_res/2;
		
		if (row<0 || row>img_res-1 || col<0 || col>img_res-1)
			particle_active[i] = false;
	}
}


void    
get_rot_sum(int start_k)
{
	int i, j, k;
	icVector2 pre_pos, cur_pos, next_pos, edge1, edge2;

	for (i=0; i<img_res; i++)
	{
		for (j=0; j<img_res; j++)
		{
			std::vector<int> &par_list = pixels_to_particles[i*img_res+j];

			float rot_sum = 0;
			//for (k=start_k; k<par_list.size()-2; k++)
			//{
			//	if (particle_active[par_list[k]] || particle_active[par_list[k]] || particle_active[par_list[k]])
			//		break;
			//	pre_pos.set(particles[par_list[k]].first, particles[par_list[k]].second);
			//	cur_pos.set(particles[par_list[k+1]].first, particles[par_list[k+1]].second);
			//	next_pos.set(particles[par_list[k+2]].first, particles[par_list[k+2]].second);
			//	edge1 = cur_pos - pre_pos;
			//	edge2 = next_pos - cur_pos;
			//	float ang1 = atan2(edge1.entry[1], edge1.entry[0]);
			//	float ang2 = atan2(edge2.entry[1], edge2.entry[0]);
			//	float ang_diff = ang2 -ang1;

			//	if (ang_diff > PI)
			//		ang_diff -= 2*PI;
			//	else if (ang_diff < -PI)
			//		ang_diff += 2*PI;

			//	rot_sum += ang_diff;
			//}

			for (k=par_list.size()-1; k>=2; k--)
			{
				if (!particle_active[par_list[k]] || !particle_active[par_list[k-1]] || !particle_active[par_list[k-2]])
					break;
				pre_pos.set(particles[par_list[k]].first, particles[par_list[k]].second);
				cur_pos.set(particles[par_list[k-1]].first, particles[par_list[k-1]].second);
				next_pos.set(particles[par_list[k-2]].first, particles[par_list[k-2]].second);
				edge1 = cur_pos - pre_pos;
				edge2 = next_pos - cur_pos;
				float ang1 = atan2(edge1.entry[1], edge1.entry[0]);
				float ang2 = atan2(edge2.entry[1], edge2.entry[0]);
				float ang_diff = ang2 -ang1;

				if (ang_diff > PI)
					ang_diff -= 2*PI;
				else if (ang_diff < -PI)
					ang_diff += 2*PI;

				rot_sum += ang_diff;
			}
			total_rotation[i][j] = rot_sum;
		}
	}
}

//void
void 
init_particle_list_test()
{
	particles.clear();
	pixels_to_particles.clear();
	pixels_to_particles.resize(PDIM*PDIM);
	particle_active.clear();
}

void 
insert_new_particles_for_streakline_test()
{
	int i, j;
	
	// for double gyre 

	double interval = 2./(PDIM-1); 

	for (i=0; i<PDIM; i++)
	{
		for (j=0; j<PDIM; j++)
		{
			// add a new particle
			std::pair<float, float> par;
			par.first = 2*(i+.5)/PDIM;
			par.second = 2*(j+.5)/PDIM;

			particles.push_back(par);
			particle_active.push_back(true);

			// add to the pixels_to_particles pointers
			pixels_to_particles[i*PDIM+j].push_back(particles.size()-1);
		}
	}
}


//// Test 3D steady flow

// the order of the points matters!!!!!
// We may also need to know the preceding rotation axis
double 
get_rot_ang_3D(double p1[3], double p2[3], double p3[3], icVector3 &pre_rotAxis)  
{
	icVector3 evec1(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]);
	icVector3 evec2(p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]);

	icVector3 cur_rotAxis = cross (evec1, evec2); 

	icVector3 B = cross (cur_rotAxis, evec1); // here we use evec1 as X, cur_rotAxis as Z to compute Y, in order to form a local coordinate system

	double a = dot(evec2, evec1);
	double b = dot(evec2, B);

	double ang = atan2(b, a);
	pre_rotAxis = cur_rotAxis;
	return ang;
}

// 1. ABC flow
//double AA = 1.;
//double BB = sqrt(2./3.);
//double CC = 1./sqrt(3.);


void 
get_vec_ABC_flow(double x, double y, double z, double &vx, double &vy, double &vz)
{
	vx = AA*sin(z) + CC*cos(y);
	vy = BB*sin(x) + AA*cos(z);
	vz = CC*sin(y) + BB*cos(x);
}


void    
get_vec_Lorenz_eq(double x, double y, double z, double &vx, double &vy, double &vz)
{
	vx = Lorenz_A*(y-x);
	vy = Lorenz_B*x-y-x*z;
	vz = x*y-Lorenz_C*z;
}

void    
get_vec_Rossler_eq(double x, double y, double z, double &vx, double &vy, double &vz)
{
	vx = -y-z;
	vy = x+Rossler_A*y;
	vz = Rossler_B+z*(x-Rossler_C);
}


void
get_vec_general_3D( double x, double y, double z, double &vxp, double &vyp, double &vzp )
{
	//*vxp = -3 + 6.*x - 4.*x*(y+1.) - 4.*z;
	//*vyp = 12.*x - 4.*x*x - 12.*z + 4.*z*z;
	//*vzp = 3. + 4.*x - 4.*x*(y+1.) - 6.*z + 4.*(y+1.)*z;
	//float a = -4., b = 1.;
	//*vxp = a*x+y*z;
	//*vyp = b*y+x*z;
	//*vzp = z - x*y;
	vxp = 2./3.*y*z;
	vyp = 2./3.*x*z;
	vzp = - 4./3.*x*y;
}


void 
get_vec_Chen_Lee_eq(double x, double y, double z, double &vx, double &vy, double &vz)
{
	vx = -z*y + Chen_Lee_A*x;
	vy = x*z + Chen_Lee_B*y;
	vz = one_third*y*x + Chen_Lee_C*z;
}


void    
get_vec_RF_eq(double x, double y, double z, double &vx, double &vy, double &vz)
{
	vx = y*(z-1+x*x) + RF_A*x;
	vy = x*(3*z+1-x*x) + RF_A*y;
	vz = -2*z*(RF_B+x*y);
}



// Accumulate the rotation for the streamline starting at a pixel center
void 
comp_rot_sum_3D_steady_at(int i, int j, int k, int L)
{
	double x = i+.5;
	double y = j+.5;
	double z = k+.5;

	double startx = x = x/GRIDSIZE*(MAX_X-MIN_X)+MIN_X;
	double starty = y = y/GRIDSIZE*(MAX_Y-MIN_Y)+MIN_Y;
	double startz = z = z/GRIDSIZE*(MAX_Z-MIN_Z)+MIN_Z;

	std::vector<double> streamline;
	streamline.push_back(x);
	streamline.push_back(y);
	streamline.push_back(z);

	double step_size = INI_STEP_SIZE/10.;
	// forward tracing
	double total_rot = 0;
	double vec[3];
	for (int i=0; i<L; i++)
	{
		//get_vec_ABC_flow(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Lorenz_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Rossler_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_general_3D(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Chen_Lee_eq(x, y, z, vec[0], vec[1], vec[2]);
		get_vec_RF_eq(x, y, z, vec[0], vec[1], vec[2]);

		double len = get_vec3_len <double> (vec);

		if (len < 1.e-6) break;

		// normalize the vector
		vec[0] /= len;
		vec[1] /= len;
		vec[2] /= len;

		x += (step_size*vec[0]);
		y += (step_size*vec[1]);
		z += (step_size*vec[2]);

		streamline.push_back(x);
		streamline.push_back(y);
		streamline.push_back(z);

		if (x<MIN_X || x>MAX_X || y<MIN_Y || y>MAX_Y || z<MIN_Z || z>MAX_Z)
			break;		
	}
	// compute the rotation sum
	icVector3 pre_axis, cur_axis;
	double pre_ang, cur_ang = 0;
	
	for (int ii=0; ii<streamline.size()-6; ii+=3)
	{
		double p1[3], p2[3], p3[3];

		p1[0] = streamline[ii];
		p1[1] = streamline[ii+1];
		p1[2] = streamline[ii+2];

		p2[0] = streamline[ii+3];
		p2[1] = streamline[ii+4];
		p2[2] = streamline[ii+5];

		p3[0] = streamline[ii+6];
		p3[1] = streamline[ii+7];
		p3[2] = streamline[ii+8];

		cur_ang = get_rot_ang_3D(p1, p2, p3, cur_axis);

		//if (i>0)
		//{
			total_rot += cur_ang;  // it seems that the function "get_rot_ang_3D" should take care of the SIGN of the rotation (POTENTIAL BUG)
		//}

		pre_ang = cur_ang;
		pre_axis = cur_axis;
	}
	rot_sum_3d[i][j][k] = total_rot;

	// Now compute the backward tracing
	streamline.clear();

	x = startx; 
	y = starty;
	z = startz;

	streamline.push_back(x);
	streamline.push_back(y);
	streamline.push_back(z);

	for (int i=0; i<L; i++)
	{
		//get_vec_ABC_flow(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Lorenz_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Rossler_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_general_3D(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Chen_Lee_eq(x, y, z, vec[0], vec[1], vec[2]);
		get_vec_RF_eq(x, y, z, vec[0], vec[1], vec[2]);

		double len = get_vec3_len <double> (vec);

		if (len < 1.e-6) break;

		// normalize the vector
		vec[0] /= len;
		vec[1] /= len;
		vec[2] /= len;

		x -= (step_size*vec[0]);
		y -= (step_size*vec[1]);
		z -= (step_size*vec[2]);

		streamline.push_back(x);
		streamline.push_back(y);
		streamline.push_back(z);

		if (x<MIN_X || x>MAX_X || y<MIN_Y || y>MAX_Y || z<MIN_Z || z>MAX_Z)
			break;		
	}

	// compute the rotation sum
	total_rot = 0.;
	for (int j=0; j<streamline.size()-6; j+=3)
	{
		double p1[3], p2[3], p3[3];

		p1[0] = streamline[j];
		p1[1] = streamline[j+1];
		p1[2] = streamline[j+2];

		p2[0] = streamline[j+3];
		p2[1] = streamline[j+4];
		p2[2] = streamline[j+5];

		p3[0] = streamline[j+6];
		p3[1] = streamline[j+7];
		p3[2] = streamline[j+8];

		cur_ang = get_rot_ang_3D(p1, p2, p3, cur_axis);

		//if (i>0)
		//{
			total_rot -= cur_ang;  // it seems that the function "get_rot_ang_3D" should take care of the SIGN of the rotation (POTENTIAL BUG)
		//}

		pre_ang = cur_ang;
		pre_axis = cur_axis;
	}
	rot_sum_3d[i][j][k] += total_rot;
}

void 
comp_rot_sum_3D_steady(int L)
{
	int i, j, k;
	for (i = 0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				comp_rot_sum_3D_steady_at(i, j, k, L);
				printf("finish (%d, %d, %d)\n", i, j, k);
			}
		}
	}	

	int stp=0;
}

//float rot_sum_3d[GRIDSIZE][GRIDSIZE][GRIDSIZE];
//float rot_sum_rgb[GRIDSIZE][GRIDSIZE][GRIDSIZE][3];

void    
get_rot_sum_3d_color_map()
{
	int i, j, k;
	
	float hsv[3], rgb[3];

	double max_rot = -1.e8, min_rot = 1.e8;
	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				if (rot_sum_3d[i][j][k] > max_rot) max_rot = rot_sum_3d[i][j][k];
				if (rot_sum_3d[i][j][k] < min_rot) min_rot = rot_sum_3d[i][j][k];
			}
		}
	}

	double range = max_rot - min_rot;

	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				hsv[2] = 1.;

				hsv[1] = 1. - 2.*(rot_sum_3d[i][j][k] - min_rot)/range;

				if (hsv[1] < 0)
				{
					hsv[0] = 0.;
					hsv[1] = abs(hsv[1]);
				}
				else
					hsv[0] = 240.;

				HsvRgb (hsv, rgb);

				rot_sum_rgb[i][j][k][0] = rgb[0];
				rot_sum_rgb[i][j][k][1] = rgb[1];
				rot_sum_rgb[i][j][k][2] = rgb[2];
			}
		}
	}
}


/////FillXY
void
FillXY( void )
{
    int x, y, z, zz;
    float alpha;            /* opacity at this voxel    */
    float r, g, b;            /* running color composite    */

    for( x = 0; x < GRIDSIZE; x++ )
    {
        for( y = 0; y < GRIDSIZE; y++ )
        {
            r = g = b = 0.;
            for( zz = 0; zz < GRIDSIZE; zz++ )
            {
                /* which direction to fill:    */

                if( Zside == PLUS )
                    z = zz;
                else
                    z = ( GRIDSIZE-1 ) - zz;


				if (Show_ROTSUM_or_DIFF ==0 )
				{
					if( rot_sum_3d[x][y][z] < SRange[0] || rot_sum_3d[x][y][z] > SRange[1])
					{
						r = g = b = 0.;
						alpha = 0.;
					}
					else
					{
						r = rot_sum_rgb[x][y][z][0];
						g = rot_sum_rgb[x][y][z][1];
						b = rot_sum_rgb[x][y][z][2];
						alpha = MaxAlpha;
					}
				}

				else
				{
					//float len = get_vec3_len <float> (rot_sum_3d_diff[x][y][z]);
					if( rot_sum_3d_diff_mag[x][y][z] < SRange[0] || rot_sum_3d_diff_mag[x][y][z] > SRange[1])
					{
						r = g = b = 0.;
						alpha = 0.;
					}
					else
					{
						r = rot_sum_diff_rgb[x][y][z][0];
						g = rot_sum_diff_rgb[x][y][z][1];
						b = rot_sum_diff_rgb[x][y][z][2];
						alpha = MaxAlpha;
					}
				}

				TextureXY[zz][y][x][0] = (unsigned char) ( 255.*r + .5 );
                TextureXY[zz][y][x][1] = (unsigned char) ( 255.*g + .5 );
                TextureXY[zz][y][x][2] = (unsigned char) ( 255.*b + .5 );
                TextureXY[zz][y][x][3] = (unsigned char) ( 255.*alpha + .5 );
            }
        }
    }
}


/////FillXZ
void
FillXZ( void )
{
    int x, y, z, yy;
    float alpha;            /* opacity at this voxel    */
    float r, g, b;            /* running color composite    */
    for( x = 0; x < GRIDSIZE; x++ )
    {
        for( z = 0; z < GRIDSIZE; z++ )
        {
            r = g = b = 0.;
            for( yy = 0; yy < GRIDSIZE; yy++ )
            {
                // which direction to fill:    

                if( Yside == PLUS )
                    y = yy;
                else
                    y = ( GRIDSIZE-1 ) - yy;


				if (Show_ROTSUM_or_DIFF ==0 )
				{
					if( rot_sum_3d[x][y][z] < SRange[0] || rot_sum_3d[x][y][z] > SRange[1])
					{
						r = g = b = 0.;
						alpha = 0.;
					}
					else
					{
						r = rot_sum_rgb[x][y][z][0];
						g = rot_sum_rgb[x][y][z][1];
						b = rot_sum_rgb[x][y][z][2];
						alpha = MaxAlpha;
					}
				}

				else
				{
					//float len = get_vec3_len <float> (rot_sum_3d_diff[x][y][z]);
					if( rot_sum_3d_diff_mag[x][y][z] < SRange[0] || rot_sum_3d_diff_mag[x][y][z] > SRange[1])
					{
						r = g = b = 0.;
						alpha = 0.;
					}
					else
					{
						r = rot_sum_diff_rgb[x][y][z][0];
						g = rot_sum_diff_rgb[x][y][z][1];
						b = rot_sum_diff_rgb[x][y][z][2];
						alpha = MaxAlpha;
					}
				}

                TextureXZ[yy][z][x][0] = (unsigned char) ( 255.*r + .5 );
                TextureXZ[yy][z][x][1] = (unsigned char) ( 255.*g + .5 );
                TextureXZ[yy][z][x][2] = (unsigned char) ( 255.*b + .5 );
                TextureXZ[yy][z][x][3] = (unsigned char) ( 255.*alpha + .5 );
            }
        }
    }


}

/////FillYZ
void
FillYZ( void )
{
    int x, y, z, xx;
    float alpha;            /* opacity at this voxel    */
    float r, g, b;            /* running color composite    */

    for( y = 0; y < GRIDSIZE; y++ )
    {
        for( z = 0; z < GRIDSIZE; z++ )
        {
            r = g = b = 0.;
            for( xx = 0; xx < GRIDSIZE; xx++ )
            {
                /* which direction to fill:    */

                if( Xside == PLUS )
                    x = xx;
                else
                    x = ( GRIDSIZE-1 ) - xx;


				if (Show_ROTSUM_or_DIFF ==0 )
				{
					if( rot_sum_3d[x][y][z] < SRange[0] || rot_sum_3d[x][y][z] > SRange[1])
					{
						r = g = b = 0.;
						alpha = 0.;
					}
					else
					{
						r = rot_sum_rgb[x][y][z][0];
						g = rot_sum_rgb[x][y][z][1];
						b = rot_sum_rgb[x][y][z][2];
						alpha = MaxAlpha;
					}
				}
				else
				{
					//float len = get_vec3_len <float> (rot_sum_3d_diff[x][y][z]);
					if( rot_sum_3d_diff_mag[x][y][z] < SRange[0] || rot_sum_3d_diff_mag[x][y][z] > SRange[1])
					{
						r = g = b = 0.;
						alpha = 0.;
					}
					else
					{
						r = rot_sum_diff_rgb[x][y][z][0];
						g = rot_sum_diff_rgb[x][y][z][1];
						b = rot_sum_diff_rgb[x][y][z][2];
						alpha = MaxAlpha;
					}
				}

                TextureYZ[xx][z][y][0] = (unsigned char) ( 255.*r + .5 );
                TextureYZ[xx][z][y][1] = (unsigned char) ( 255.*g + .5 );
                TextureYZ[xx][z][y][2] = (unsigned char) ( 255.*b + .5 );
                TextureYZ[xx][z][y][3] = (unsigned char) ( 255.*alpha + .5 );
            }
        }
    }
}


void    UpdateXY( void )
{
}
void    UpdateXZ( void )
{
}

void    UpdateYZ( void )
{
}

/**
 ** determine which sides of the cube are facing the user
 **
 **/
void
DetermineVisibility()
{
    float xr, yr;
    float cx, sx;
    float cy, sy;
    float nzx, nzy, nzz;    /* z component of normal for x side, y side, and z side    */

    xr = Xrot * ( M_PI/180. );
    yr = Yrot * ( M_PI/180. );

    cx = cos( xr );
    sx = sin( xr );
    cy = cos( yr );
    sy = sin( yr );

    nzx = -sy;
    nzy =  sx * cy;
    nzz =  cx * cy;


    /* which sides of the cube are showing:                */
    /* the Xside being shown to the user is MINUS or PLUS        */

    Xside = ( nzx > 0.  ?  PLUS  :  MINUS );
    Yside = ( nzy > 0.  ?  PLUS  :  MINUS );
    Zside = ( nzz > 0.  ?  PLUS  :  MINUS );


    /* which direction needs to be composited:            */

    if( fabs(nzx) > fabs(nzy)  &&  fabs(nzx) > fabs(nzz)  )
        Major = X;
    else if( fabs(nzy) > fabs(nzx)  &&  fabs(nzy) > fabs(nzz)  )
        Major = Y;
    else
        Major = Z;
}



void    
Drawing()
{
	float x0, dx, xcoord;
	float y0, dy, ycoord;
	float z0, dz, zcoord;
	int x, y, z;

//	glDisable( GL_COLOR_MATERIAL);

    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );

    int filter = GL_NEAREST;
    if( Bilinear )
        filter = GL_LINEAR;

    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter );
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
    glEnable( GL_TEXTURE_2D );


    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glEnable( GL_BLEND );

	/////////////
	///Test code here
    //Major = Z; ////11/12 11:20PM

    if( Major == Z )
    {
        if( Zside == PLUS )
        {
            z0 = -10.;
            dz = 20. / (float)( GRIDSIZE - 1 );
        }
        else
        {
            z0 = 10.;
            dz = -20. / (float)( GRIDSIZE - 1 );
        }

        for( z = 0, zcoord = z0; z < GRIDSIZE; z++, zcoord += dz )
        {
            glTexImage2D( GL_TEXTURE_2D, 0, 4, GRIDSIZE, GRIDSIZE, 0, GL_RGBA, GL_UNSIGNED_BYTE, &TextureXY[z][0][0][0] );

            glBegin( GL_QUADS );
            
                glTexCoord2f( 0.f, 0.f );
                glVertex3f( -10.f, -10.f, zcoord );
    
                glTexCoord2f( 1.f, 0.f );
                glVertex3f( 10.f, -10.f, zcoord );
    
                glTexCoord2f( 1.f, 1.f );
                glVertex3f( 10.f, 10.f, zcoord );
    
                glTexCoord2f( 0.f, 1.f );
                glVertex3f( -10.f, 10.f, zcoord );

            glEnd();
        }
    }

	else if( Major == Y )
    {
        if( Yside == PLUS )
        {
            y0 = -10.;
            dy = 20. / (float)( GRIDSIZE - 1 );
        }
        else
        {
            y0 = 10.;
            dy = -20. / (float)( GRIDSIZE - 1 );
        }

        for( y = 0, ycoord = y0; y < GRIDSIZE; y++, ycoord += dy )
        {
            glTexImage2D( GL_TEXTURE_2D, 0, 4, GRIDSIZE, GRIDSIZE, 0, GL_RGBA, GL_UNSIGNED_BYTE, &TextureXZ[y][0][0][0] );

            glBegin( GL_QUADS );
            
                glTexCoord2f( 0.f, 0.f );
                glVertex3f( -10.f, ycoord, -10.f );
    
                glTexCoord2f( 1.f, 0.f );
                glVertex3f( 10.f, ycoord, -10.f );
    
                glTexCoord2f( 1.f, 1.f );
                glVertex3f( 10.f, ycoord, 10.f );
    
                glTexCoord2f( 0.f, 1.f );
                glVertex3f( -10.f, ycoord,  10.f);

            glEnd();
        }
    }

	else if( Major == X )
    {
        if( Xside == PLUS )
        {
            x0 = -10.;
            dx = 20. / (float)( GRIDSIZE - 1 );
        }
        else
        {
            x0 = 10.;
            dx = -20. / (float)( GRIDSIZE - 1 );
        }

        for( x = 0, xcoord = x0; x < GRIDSIZE; x++, xcoord += dx )
        {
            glTexImage2D( GL_TEXTURE_2D, 0, 4, GRIDSIZE, GRIDSIZE, 0, GL_RGBA, GL_UNSIGNED_BYTE, &TextureYZ[x][0][0][0] );

            glBegin( GL_QUADS );
            
                glTexCoord2f( 0.f, 0.f );
                glVertex3f( xcoord, -10.f, -10.f );
    
                glTexCoord2f( 1.f, 0.f );
                glVertex3f( xcoord, 10.f, -10.f );
    
                glTexCoord2f( 1.f, 1.f );
                glVertex3f( xcoord, 10.f,  10.f );
    
                glTexCoord2f( 0.f, 1.f );
                glVertex3f( xcoord, -10.f, 10.f);

            glEnd();
        }
    }

    glDisable( GL_BLEND );


	glDisable( GL_TEXTURE_2D );


}




//float rot_sum_3d_diff[GRIDSIZE][GRIDSIZE][GRIDSIZE];
//float rot_sum_diff_rgb[GRIDSIZE][GRIDSIZE][GRIDSIZE][3];

void    comp_rot_sum_diff_3d()
{
	// 
	int i, j, k;

	// First compute the difference along X direction
	for (i=0; i<GRIDSIZE; i++)
	{
		if (i==0)
		{
			for (j=0; j<GRIDSIZE; j++)
			{
				for (k=0; k<GRIDSIZE; k++)
				{
					rot_sum_3d_diff[i][j][k][2] = rot_sum_3d[i+1][j][k]-rot_sum_3d[i][j][k];
				}
			}
		}

		else if (i == GRIDSIZE-1)
		{
			for (j=0; j<GRIDSIZE; j++)
			{
				for (k=0; k<GRIDSIZE; k++)
				{
					rot_sum_3d_diff[i][j][k][2] = rot_sum_3d[i][j][k]-rot_sum_3d[i-1][j][k];
				}
			}
		}
		else
		{
			for (j=0; j<GRIDSIZE; j++)
			{
				for (k=0; k<GRIDSIZE; k++)
				{
					rot_sum_3d_diff[i][j][k][2] = 0.5*(rot_sum_3d[i+1][j][k]-rot_sum_3d[i-1][j][k]);
				}
			}
		}
	}


	// compute the difference along Y direction
	for (j=0; j<GRIDSIZE; j++)
	{
		if (j==0)
		{
			for (i=0; i<GRIDSIZE; i++)
			{
				for (k=0; k<GRIDSIZE; k++)
				{
					rot_sum_3d_diff[i][j][k][1] = rot_sum_3d[i][j+1][k]-rot_sum_3d[i][j][k];
				}
			}
		}

		else if (j == GRIDSIZE-1)
		{
			for (i=0; i<GRIDSIZE; i++)
			{
				for (k=0; k<GRIDSIZE; k++)
				{
					rot_sum_3d_diff[i][j][k][1] = rot_sum_3d[i][j][k]-rot_sum_3d[i][j-1][k];
				}
			}
		}
		else
		{
			for (i=0; i<GRIDSIZE; i++)
			{
				for (k=0; k<GRIDSIZE; k++)
				{
					rot_sum_3d_diff[i][j][k][1] = 0.5*(rot_sum_3d[i][j+1][k]-rot_sum_3d[i][j-1][k]);
				}
			}
		}
	}

	for (k=0; k<GRIDSIZE; k++)
	{
		if (k==0)
		{
			for (i=0; i<GRIDSIZE; i++)
			{
				for (j=0; j<GRIDSIZE; j++)
				{
					rot_sum_3d_diff[i][j][k][0] = rot_sum_3d[i][j][k+1]-rot_sum_3d[i][j][k];
				}
			}
		}
		else if (k==GRIDSIZE-1)
		{
			for (i=0; i<GRIDSIZE; i++)
			{
				for (j=0; j<GRIDSIZE; j++)
				{
					rot_sum_3d_diff[i][j][k][0] = rot_sum_3d[i][j][k]-rot_sum_3d[i][j][k-1];
				}
			}
		}
		else
		{
			for (i=0; i<GRIDSIZE; i++)
			{
				for (j=0; j<GRIDSIZE; j++)
				{
					rot_sum_3d_diff[i][j][k][0] = 0.5*(rot_sum_3d[i][j][k+1]-rot_sum_3d[i][j][k-1]);
				}
			}
		}
	}

	// we need to compute the logorithm
	//for (i=0; i<GRIDSIZE; i++)
	//{
	//	for (j=0; j<GRIDSIZE; j++)
	//	{
	//		for (k=0; k<GRIDSIZE; k++)
	//		{
	//			double len = get_vec3_len <float> (rot_sum_3d_diff[i][j][k]);
	//		}
	//	}
	//}
}

void    
get_color_map_for_rot_sum_diff_3d()
{
	int i, j, k;
	
	float hsv[3], rgb[3];

	double max_rot = -1.e8, min_rot = 1.e8;
	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				double len = get_vec3_len <float> (rot_sum_3d_diff[i][j][k]);

				// we need to compute the logorithm
				rot_sum_3d_diff_mag[i][j][k] = log(len)/10.;  // this is specific for ABC flow
				if (rot_sum_3d_diff_mag[i][j][k] > max_rot) max_rot = rot_sum_3d_diff_mag[i][j][k];
				if (rot_sum_3d_diff_mag[i][j][k] < min_rot) min_rot = rot_sum_3d_diff_mag[i][j][k];
			}
		}
	}

	double mean = 0.5*(max_rot - min_rot); 

	double range = max_rot - min_rot;

	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				//double len = get_vec3_len <float> (rot_sum_3d_diff[i][j][k]);
				
				hsv[2] = 1.;

				hsv[1] = 1. - 2.*(rot_sum_3d_diff_mag[i][j][k] - min_rot)/range;

				if (hsv[1] < 0)
				{
					hsv[0] = 0.;
					hsv[1] = abs(hsv[1]);
				}
				else
					hsv[0] = 240.;

				HsvRgb (hsv, rgb);

				rot_sum_diff_rgb[i][j][k][0] = rgb[0];
				rot_sum_diff_rgb[i][j][k][1] = rgb[1];
				rot_sum_diff_rgb[i][j][k][2] = rgb[2];
			}
		}
	}
}


void    
get_color_map_for_rot_sum_diff_3d_from_file()
{
	int i, j, k;
	
	float hsv[3], rgb[3];

	double max_rot = -1.e8, min_rot = 1.e8;
	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				if (rot_sum_3d_diff_mag[i][j][k] > max_rot) max_rot = rot_sum_3d_diff_mag[i][j][k];
				if (rot_sum_3d_diff_mag[i][j][k] < min_rot) min_rot = rot_sum_3d_diff_mag[i][j][k];
			}
		}
	}

	double mean = 0.5*(max_rot - min_rot); 

	double range = max_rot - min_rot;

	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				//double len = get_vec3_len <float> (rot_sum_3d_diff[i][j][k]);
				
				hsv[2] = 1.;

				hsv[1] = 1. - 2.*(rot_sum_3d_diff_mag[i][j][k] - min_rot)/range;

				if (hsv[1] < 0)
				{
					hsv[0] = 0.;
					hsv[1] = abs(hsv[1]);
				}
				else
					hsv[0] = 240.;

				HsvRgb (hsv, rgb);

				rot_sum_diff_rgb[i][j][k][0] = rgb[0];
				rot_sum_diff_rgb[i][j][k][1] = rgb[1];
				rot_sum_diff_rgb[i][j][k][2] = rgb[2];
			}
		}
	}
}


void    output_rot_sum_result(char *filename)
{
	// output the file into ASCII format at this moment
	// HEADER: NX, NY, NZ  (resolution along each dimension)
	// HEADER: MINX, MINY, MINZ, MAXX, MAXY, MAXZ
	FILE *fp = fopen(filename, "w");
	int i, j, k;  // k - X dimension, j - Y dimension, i - Z dimension
	fprintf(fp, "%d %d %d\n", GRIDSIZE, GRIDSIZE, GRIDSIZE);
	fprintf(fp, "%f %f %f %f %f %f\n", MIN_X, MIN_Y, MIN_Z, MAX_X, MAX_Y, MAX_Z);

	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				fprintf(fp, "%f ", rot_sum_3d[i][j][k]);
			}
		}
	}
	fclose(fp);
}

void    output_rot_sum_diff_result(char *filename)
{
	// output the file into ASCII format at this moment
	// HEADER: NX, NY, NZ  (resolution along each dimension)
	// HEADER: MINX, MINY, MINZ, MAXX, MAXY, MAXZ
	FILE *fp = fopen(filename, "w");
	int i, j, k;  // k - X dimension, j - Y dimension, i - Z dimension
	fprintf(fp, "%d %d %d\n", GRIDSIZE, GRIDSIZE, GRIDSIZE);
	fprintf(fp, "%f %f %f %f %f %f\n", MIN_X, MIN_Y, MIN_Z, MAX_X, MAX_Y, MAX_Z);

	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				//fprintf(fp, "%f %f %f\n", rot_sum_3d_diff[i][j][k][0],rot_sum_3d_diff[i][j][k][1],rot_sum_3d_diff[i][j][k][2]);
				fprintf(fp, "%f ", rot_sum_3d_diff_mag[i][j][k]);
			}
		}
	}
	fclose(fp);
}

void    load_rot_sum_result(char *filename)
{
	FILE *fp = fopen(filename, "r");
	int i, j, k;
	int gridsize_x, gridsize_y, gridsize_z;
	float minx, miny, minz, maxx, maxy, maxz;
	fscanf(fp, "%d %d %d\n", &gridsize_x, &gridsize_y, &gridsize_z);
	fscanf(fp, "%f %f %f %f %f %f\n", &minx, &miny, &minz, &maxx, &maxy, &maxz);

	MIN_X = minx, MIN_Y = miny, MIN_Z = minz, MAX_X = maxx, MAX_Y = maxy, MAX_Z = maxz;

	for (i=0; i<gridsize_x; i++)
	{
		for (j=0; j<gridsize_y; j++)
		{
			for (k=0; k<gridsize_z; k++)
			{
				fscanf(fp, "%f ", &rot_sum_3d_diff_mag[i][j][k]);
			}
		}
	}
	fclose(fp);

	data_center.entry[0] = (maxx-minx)/2;
	data_center.entry[1] = (maxy-miny)/2;
	data_center.entry[2] = (maxz-minz)/2;

	printf ("data range (%f, %f) x (%f, %f) x (%f, %f), center at (%f, %f, %f)\n", minx, maxx, miny, maxy, minz, maxz,
		data_center.entry[0], data_center.entry[1], data_center.entry[2]);
}


void 
comp_streamline_3D_steady_at(int i, int j, int k, int L, std::vector<Point3D> &one_line)
{
	double x = i+.5;
	double y = j+.5;
	double z = k+.5;

	double startx = x = x/GRIDSIZE*(MAX_X-MIN_X)+MIN_X;
	double starty = y = y/GRIDSIZE*(MAX_Y-MIN_Y)+MIN_Y;
	double startz = z = z/GRIDSIZE*(MAX_Z-MIN_Z)+MIN_Z;

	//std::vector<double> streamline;
	//streamline.push_back(x);
	//streamline.push_back(y);
	//streamline.push_back(z);

	std::vector<Point3D> tmp_forward;

	Point3D pt;
	pt.xyz[0] = x;
	pt.xyz[1] = y;
	pt.xyz[2] = z;
	tmp_forward.push_back(pt);

	double step_size = INI_STEP_SIZE/10.;
	// forward tracing
	double total_rot = 0;
	double vec[3];
	for (int i=0; i<L; i++)
	{
		get_vec_ABC_flow(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Lorenz_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Rossler_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_general_3D(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Chen_Lee_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_RF_eq(x, y, z, vec[0], vec[1], vec[2]);

		double len = get_vec3_len <double> (vec);

		if (len < 1.e-6) break;

		// normalize the vector
		vec[0] /= len;
		vec[1] /= len;
		vec[2] /= len;

		x += (step_size*vec[0]);
		y += (step_size*vec[1]);
		z += (step_size*vec[2]);

		//streamline.push_back(x);
		//streamline.push_back(y);
		//streamline.push_back(z);

		pt.xyz[0] = x;
		pt.xyz[1] = y;
		pt.xyz[2] = z;
		tmp_forward.push_back(pt);

		if (x<MIN_X || x>MAX_X || y<MIN_Y || y>MAX_Y || z<MIN_Z || z>MAX_Z)
			break;		
	}

	//// compute the rotation sum
	//icVector3 pre_axis, cur_axis;
	//double pre_ang, cur_ang = 0;
	//
	//for (int ii=0; ii<streamline.size()-6; ii+=3)
	//{
	//	double p1[3], p2[3], p3[3];

	//	p1[0] = streamline[ii];
	//	p1[1] = streamline[ii+1];
	//	p1[2] = streamline[ii+2];

	//	p2[0] = streamline[ii+3];
	//	p2[1] = streamline[ii+4];
	//	p2[2] = streamline[ii+5];

	//	p3[0] = streamline[ii+6];
	//	p3[1] = streamline[ii+7];
	//	p3[2] = streamline[ii+8];

	//	cur_ang = get_rot_ang_3D(p1, p2, p3, cur_axis);

	//	//if (i>0)
	//	//{
	//		total_rot += cur_ang;  // it seems that the function "get_rot_ang_3D" should take care of the SIGN of the rotation (POTENTIAL BUG)
	//	//}

	//	pre_ang = cur_ang;
	//	pre_axis = cur_axis;
	//}
	//rot_sum_3d[i][j][k] = total_rot;

	//// Now compute the backward tracing
	//streamline.clear();

	x = startx; 
	y = starty;
	z = startz;

	//streamline.push_back(x);
	//streamline.push_back(y);
	//streamline.push_back(z);

	std::vector<Point3D> tmp_backward;

	for (int i=0; i<L; i++)
	{
		get_vec_ABC_flow(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Lorenz_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Rossler_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_general_3D(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_Chen_Lee_eq(x, y, z, vec[0], vec[1], vec[2]);
		//get_vec_RF_eq(x, y, z, vec[0], vec[1], vec[2]);

		double len = get_vec3_len <double> (vec);

		if (len < 1.e-6) break;

		// normalize the vector
		vec[0] /= len;
		vec[1] /= len;
		vec[2] /= len;

		x -= (step_size*vec[0]);
		y -= (step_size*vec[1]);
		z -= (step_size*vec[2]);

		//streamline.push_back(x);
		//streamline.push_back(y);
		//streamline.push_back(z);
		pt.xyz[0] = x;
		pt.xyz[1] = y;
		pt.xyz[2] = z;
		tmp_backward.push_back(pt);

		if (x<MIN_X || x>MAX_X || y<MIN_Y || y>MAX_Y || z<MIN_Z || z>MAX_Z)
			break;		
	}

	//// compute the rotation sum
	//total_rot = 0.;
	//for (int j=0; j<streamline.size()-6; j+=3)
	//{
	//	double p1[3], p2[3], p3[3];

	//	p1[0] = streamline[j];
	//	p1[1] = streamline[j+1];
	//	p1[2] = streamline[j+2];

	//	p2[0] = streamline[j+3];
	//	p2[1] = streamline[j+4];
	//	p2[2] = streamline[j+5];

	//	p3[0] = streamline[j+6];
	//	p3[1] = streamline[j+7];
	//	p3[2] = streamline[j+8];

	//	cur_ang = get_rot_ang_3D(p1, p2, p3, cur_axis);

	//	//if (i>0)
	//	//{
	//		total_rot -= cur_ang;  // it seems that the function "get_rot_ang_3D" should take care of the SIGN of the rotation (POTENTIAL BUG)
	//	//}

	//	pre_ang = cur_ang;
	//	pre_axis = cur_axis;
	//}
	//rot_sum_3d[i][j][k] += total_rot;

	// Merge the backward and forward tracing

	for (i=tmp_backward.size()-1; i>=0; i--)
		one_line.push_back(tmp_backward[i]);
	for (i=0; i<tmp_forward.size(); i++)
		one_line.push_back(tmp_forward[i]);
	tmp_backward.clear();
	tmp_forward.clear();
}

void 
update_list_streamlines_3d()
{
	// According to the user specified threshold, compute streamlines from the centers of the voxels which are in the range
	int i, j, k;

	for (i=0; i<GRIDSIZE; i++)
		for (j=0; j<GRIDSIZE; j++)
			for (k=0; k<GRIDSIZE; k++)
				cell_counters_3d[i][j][k]=0;
	
	for (i=0; i<list_streamlines_3d.size(); i++)
		list_streamlines_3d[i].one_line.clear();

	list_streamlines_3d.clear();
	double dx = (MAX_X-MIN_X)/(GRIDSIZE-1);
	double dy = (MAX_Y-MIN_Y)/(GRIDSIZE-1);
	double dz = (MAX_Z-MIN_Z)/(GRIDSIZE-1);
	double xrang = MAX_X-MIN_X;
	double yrang = MAX_Y-MIN_Y;
	double zrang = MAX_Z-MIN_Z;

	for (i=0; i<GRIDSIZE; i++)
	{
		for (j=0; j<GRIDSIZE; j++)
		{
			for (k=0; k<GRIDSIZE; k++)
			{
				if( rot_sum_3d_diff_mag[i][j][k] < SRange[0] || rot_sum_3d_diff_mag[i][j][k] > SRange[1])
					continue;

				if (cell_counters_3d[i][j][k]>1) continue;

				Streamline3D line;
				line.i = i;
				line.j = j;
				line.k = k;

				// start a streamline at the center of this voxel
				comp_streamline_3D_steady_at(i, j, k, L, line.one_line);
				list_streamlines_3d.push_back(line);

				// increase the counters of the voxels that the streamline passes through
				int pre_i=0, pre_j=0, pre_k=0, cur_i=0, cur_j=0, cur_k=0;

				for (int ii=0; ii<line.one_line.size(); ii++)
				{
					cur_i = (line.one_line[ii].xyz[0]-MIN_X)/xrang*GRIDSIZE;
					cur_j = (line.one_line[ii].xyz[1]-MIN_Y)/yrang*GRIDSIZE;
					cur_k = (line.one_line[ii].xyz[2]-MIN_Z)/zrang*GRIDSIZE;

					if (cur_i<0||cur_j<0||cur_k<0||cur_i>=GRIDSIZE||cur_j>=GRIDSIZE||cur_k>=GRIDSIZE)
						continue;

					if (ii==0)
					{
						pre_i=cur_i; pre_j=cur_j; pre_k=cur_k;
					}

					if ( ii > 0 && cur_i==pre_i && cur_j==pre_j && cur_k==pre_k) continue;

					cell_counters_3d[cur_i][cur_j][cur_k]++;

					pre_i=cur_i; pre_j=cur_j; pre_k=cur_k;
				}
			}
		}
	}
}



void   
interpolate_error_bilinear(Polyhedron *this_poly)
{
	std::vector<float> err_est;
	double err_max = 0, err_min = 1.e8;

	int i, j;
	
	// estimate the error by directly comparing the inteporlated vector with the original one at each vertex
	for (i=0; i<this_poly->nverts; i++)
	{
		double xy[2], vxy[2];
		xy[0] = this_poly->vlist[i]->x * img_res;
		xy[1] = this_poly->vlist[i]->y * img_res;

		if (this_poly->vlist[i]->boundary
			|| fabs(this_poly->vlist[i]->vx)<1e-6 && fabs(this_poly->vlist[i]->vy)<1e-6)
		{
			err_est.push_back(0.);
			continue;
		}

		//estimate the position in the image space
		if (!get_vec_at_regular_grid_image_based2(xy, img_res, img_res, vxy, false))
		{
			err_est.push_back(0.);
			continue;
		}


		double dx = vxy[0] - this_poly->vlist[i]->vx;
		double dy = vxy[1] - this_poly->vlist[i]->vy;

		double err = sqrt(dx*dx+dy*dy);
		err_est.push_back(err);
	}

	for (i=0; i<err_est.size(); i++)
	{
		// We need to skip the boundary vertices!

		if (err_est[i]>err_max) err_max = err_est[i];
		if (err_est[i]<err_min) err_min = err_est[i];
	}

	// generate a color plot
	glViewport(0, 0, (GLsizei) img_res, (GLsizei) img_res);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	glClearColor( 0., 0., 0., 1. );
	glClear(GL_COLOR_BUFFER_BIT);

	glRenderMode(GL_SMOOTH);

	glDrawBuffer(GL_BACK);

	for (i=0; i<this_poly->ntris; i++) {
		Triangle *temp_t=this_poly->tlist[i];
		float rgb[3];
		rgb[2] = 0.0;
		glBegin(GL_TRIANGLES);
		for (j=0; j<3; j++)
		{
			Vertex *v = temp_t->verts[j];
			rgb[0] =  (err_est[temp_t->verts[j]->index] - err_min)/(err_max - err_min);
			rgb[1] = rgb[2] = rgb[0];
			glColor3fv (rgb);
			glVertex2f (v->x, v->y);
		}
		glEnd();
	}
	glPopMatrix();
	
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, img_res, img_res, GL_RGB, GL_UNSIGNED_BYTE, err_img);

	write_ppm_flippedY("err_img.ppm", (unsigned char*)err_img, img_res, img_res);

	printf("max_err = %f, min_err = %f\n", err_max, err_min);
}
