#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define __USE_XOPEN_EXTENDED /**< defined to access usleep with the C99 standard on Linux, see unistd.h */
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/freeglut.h>
#include <sys/stat.h>
#include <omp.h>
#include <stdarg.h>
#include <fenv.h>
/* #define NDEBUG */
#include <assert.h>

#define ISLINUX 0u
#define ISWIN32 1u
#define ISWIN64 2u

#ifdef __unix__
 #define DIRMOD ,0775)
 #include <errno.h>
 #define OSID ISLINUX
#elif defined(__WIN32__) || defined(_WIN32) || defined(__MSDOS__)
 #include <io.h>
 #define DIRMOD )
 #ifdef _WIN64
  #define OSID ISWIN64
 #else
  #define OSID ISWIN32
 #endif
#else
 #error "operating system unknown"
#endif

#define VERSION_INFO "4th version"
#define VERSION_DATE "xxx xx 2013"
#define SUBVERSION_INFO "0"
#define SUBVERSION_DATE "xxx xx xxxx"

/**
 * progress:
 * 2nd step Apr 21 2012
 * - worked out the main loop in parallel
 * 3rd step May 20 2012
 * - added a function which imports files of 'particle matrices' to facilitate the handling of particle systems
 * - added a gnuplot handler
 * - corrected and improved the sellmeier equations
 *   1st subversion May 25 2012
 * -- improved the overwriting of output files and plots
 * 4th step xxx xx 2013
 * - using restricted pointers
 * - added an assembler sincos function
 * - improved lots of functions
 * - debugged the orthogonalization function (amongst others)
 * - found a better name for the program
 * - looking for proper 'const-ness'
 * - added comments and doxygen descriptions for all functions
 * - checking code with cppcheck
 * - added assertions about the orthogonality of the polarization
 * - cleaner structure of the threshold values
 * - added the 'uchar' type to resemble 'bit values' (switches and that sort of thing)
 * - started to use rcs
 */

#define RELEASE_BUILD 0 /**< hides or prints certain information important for a public usage */
#if RELEASE_BUILD
 #define MEMORY_COUNTER 1 /**< controls the sumation and printing of allocated memory */
#else
 #define MEMORY_COUNTER 0
#endif
#define DEBUG_MOTION 0 /**< controls the printing of additional information in the viewer */
#define PRINT_PROG_OF_RAYS 1 /**< controls if the progress of the ray tracing is printed to the screen */
#define PARALLEL_PROCESSING 1 /**< controls the usage of OpenMP */
#define use_restrict 1 /**< controls if restricted pointers are used */
#define STOP_AT_ERR 0 /**< controls if 'Enter' is required at an error_msg call or not */

#if use_restrict
 #define res_pt __restrict
#else
 #define res_pt
#endif

#define MAXOF(a,b) ((a)>(b)?(a):(b))
#define MINOF(a,b) ((a)<(b)?(a):(b))
#define POW2(x) ((x)!=0.?(x)*(x):0.)
#define PRINT_TOKEN_INT(tok) fprintf(stdout,#tok " is %i\n",tok)
#define PRINT_TOKEN_DOUBLE(tok) fprintf(stdout,#tok " is %g\n",tok)
#define FP_ERR_CLEAR feclearexcept(FE_ALL_EXCEPT)
#define FP_ERR_CHECK fp_error(__FILE__,__LINE__,__func__)
#define ERR_ARG __FILE__,__LINE__,__func__

#ifdef M_PI
 #if RELEASE_BUILD
  #warning "overwriting 'M_PI' with (double)3.1415926535897932"
 #endif
#endif
#undef M_PI
#define M_PI    (double)3.1415926535897932
#ifdef M_SQRT2
 #if RELEASE_BUILD
  #warning "overwriting 'M_SQRT2' with (double)1.4142135623730951"
 #endif
#endif
#undef M_SQRT2
#define M_SQRT2 (double)1.4142135623730951
#define M_PIh   (double)1.5707963267948966
#define M_PI3h  (double)4.7123889803846898
#define M_PI2   (double)6.2831853071795864
#define M_EULER (double)2.7182818284590452
#define M_C0    (double)299792458.

#define DBL_EPSILON2 (DBL_EPSILON*2.) /**< a threshold value equaling the maximum absolute error in a dot3 function */
#define DBL_EPSILON4 (DBL_EPSILON*4.) /**< a threshold value used in different cirumstances instead of the machine epsilon */
#define DBL_EPSILON6 (DBL_EPSILON*6.) /**< a threshold value used in different cirumstances instead of the machine epsilon */
#define DBL_EPSILON10 (DBL_EPSILON*1.e1) /**< another threshold value used in different cirumstances instead of the machine epsilon */
#define DBL_EPSILON20 (DBL_EPSILON*8.e1) /**< another threshold value used in different cirumstances instead of the machine epsilon */
#define MINI_AMPLITUDE 1.e-10 /**< the minimum amplitude of a ray after which the ray is said to be exhausted */
#define MINI_INTENSITY 1.e-15 /**< the minimum intensity of a ray after which the ray is said to be exhausted */
#define MINI_TRANSMISSION_COEF 1.e-15 /**< the minimum transmission coefficient. if this value is undershoot, the transmitted ray is not traced. */
#define MINI_REFLECTION_COEF 1.e-15 /**< the minimum reflection coefficient. if this value is undershoot, the reflected ray is not traced. */
#define MAX_CHILD 50u /**< the maximum number of child rays. if this value is overshoot, the next child ray will not be traced. */
#define EPS_RAY_PROP_LEN 1.e-10 /**< the geometric distance which is added to a ray after a hit to propagate it from the point of intersection */

#define MAX_NUMBER_OF_RAYS (unsigned int)UINT_MAX
#define INFOLENGTH 64u
#define L_INFOLENGTH 256u
#define FILENAME_MAX1 (uint)(FILENAME_MAX+1)

#define DEFAULT_ZOOM 200.
#define DEFAULT_BIN_SPHERE3_ALPHA .1

/**
 * ISO C standard for complex arithmetics
 * http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/complex.h.html
 */
#ifndef complex
 #define complex _Complex
#endif

typedef enum ray_end_state_enum
{
	EMPTY_STATE,
	FIRSTTIME_HIT_STATE,
	REGULAR_HIT_STATE,
	EXHAUSTED_RAY_STATE,
	LOST_RAY_STATE,
	EXCEPTION_STATE
}ray_end_state; /**< describes what happened to a ray at the end of tracing */

typedef enum bin_hit_print_type_enum
{
	E_FIELD_AMPLITUDE,
	MOD_E_FIELD_AMPLITUDE,
	INTENSITY,
	POLARISATION_DENSITY
}bin_hit_print_type; /**< defines what should be printed with the OpenGL viewer */

typedef enum prog_ray_order_enum
{
	INIT_PROG,
	COUNT_PROG,
	PRINT_PROG,
	CLEAR_OUT
}prog_ray_order; /**< an option identifier of what shall be printed when calling prog_of_rays */

typedef enum coordinate_sys_enum
{
	CARTESIAN_CS,
	SPHERICAL_CS,
	MIXED_CS
}coordinate_sys; /**< an identifier of the coordinate system */

typedef enum incidence_enum
{
	NONE,
	OBLIQUE,
	VERTICAL,
	RIGHT_ANGLED
}incidence; /**< an identifier of the type of incidence */

typedef enum draw_opt_enum
{
	ALLOC_DATA,
	COPY_RAYS,
	GEN_LISTS,
	DRAW_EM_ALL,
	DRAW_SOME,
	FREE_EM_ALL
}draw_opt; /**< an option identifier to be used for the drawing function to specify what has to be done */

typedef enum colorfun_enum
{
	SINCOS_MAP,
	VISIBLE_LIGHT
}colorfun; /**< an identifier of the function used to color the output */

typedef enum trackball_mode_enum
{
	RESET=0,
	MOUSEMOTION,
	KNOWMOUSEBUTTON
}trackball_mode;

typedef unsigned int uint;
typedef unsigned long int ulong;
typedef unsigned char uchar;
typedef complex double cdoub;

typedef struct index2_struct
{
	uint ia,
		 ib;
}index2; /**< used to address the bins of the detections sphere easily */

typedef struct point3_struct
{
	double x[3];
}point3; /**< a point in R^3 */

typedef struct line3_struct
{
	double o[3], /**< origin */
		   r[3], /**< direction */
		   l; /**< scalar multiplier */
}line3; /**< a vector */

typedef struct plane3_struct
{
	double o[3], /**< a point in the plane */
		   n[3]; /**< the normalized normal vector */
}plane3; /**< a plane in R^3 */

typedef struct sphere3_struct
{
	double o[3], /**< the center of the sphere */
		   r; /**< the radius */
}sphere3; /**< a sphere in R^3 */

typedef struct intrsec_struct
{
	double p[3], /**< the point of the intersection */
		   normal[3], /**< the (normalized) normal vector with respect to the surface, pointing outwards */
		   angl, /**< the angle between the surface normal and the ray vector */
		   cangl, /**< the cosine of angl */
		   mu_f; /**< the permeability of the medium the ray is going to enter */
	cdoub n_f; /**< the refractive index of the medium the ray is going to enter */
	incidence incdnc; /**< identifier of the incidence */
}intrsec; /**< a structure containing information about an intersection between a ray and a particle in R^3 */

typedef struct hit_screen_struct
{
	double p[3], /**< the point of the intersection */
		   lam, /**< the wavelength of the ray */
		   cos_incdnc, /**< the cosine of the angle of incidence */
		   opol[3], /**< the 'o' (normalized) polarization-vector of the ray */
		   ppol[3]; /**< the 'p' (normalized) polarization-vector of the ray */
	cdoub coamp, /**< the complex amplitude of the 'o' polarization */
		  cpamp, /**< the complex amplitude of the 'p' polarization */
		  oint, /**< the intensity of the 'o' polarization */
		  pint; /**< the intensity of the 'p' polarization */
	uint tir; /**< the number of tir events of the ray */
	char info[INFOLENGTH + 1]; /**< an information about this event */
	ray_end_state state; /**< the precise end state of the ray */
}hit_screen; /**< a structure describe the last resort of the ray (not necessarily the screen) */

typedef struct gen_ray_info_struct
{
	ulong count_gen, /**< the number of generated rays */
		  count_hit, /**< the number of rays that hit the screen */
		  count_lost, /**< the number of rays which have been lost, i.e. found outside of the screen */
		  count_exhstd; /**< the number of rays which have been exhausted during their flight */
	char info[L_INFOLENGTH + 1]; /**< a general information */
}gen_ray_info; /**< a general information about the rays traced through the system */

typedef struct bin_hit_screen_struct
{
	double lam, /**< the wavelength of the rays */
		   rad, /**< the radius of the sphere */
		   **amp, /**< the amplitude at each bin of the screen */
		   **mod_amp, /**< the modulus of the amplitude at each bin of the screen */
		   **ray_int, /**< the intensity at each bin of the screen */
		   **pol_dist, /**< the polarization distribution at each bin of the screen */
		   **pol_dens, /**< the polarization density at each bin of the screen */
		   amp_max, /**< the maximum amplitude detected */
		   amp_sum, /**< the sum of the amplitudes detected */
		   mod_amp_max, /**< the maximum modulus of the amplitude detected */
		   mod_amp_sum, /**< the sum of the modulus of the amplitudes detected */
		   ray_int_max, /**< the maximum of the intensities detected */
		   ray_int_sum, /**< the sum of the intensities detected */
		   pol_max, /**< the maximum polarization ???density?? detected */
		   pol_sum; /**< the sum of the polarization ???densities?? detected */
	point3 *tir_coor, /**< the coordinate of a tir event in cartesian coordinates */
		   *exhstd_coor, /**< the coordinate of an exhausted ray event in cartesian coordinates */
		   *lost_coor; /**< the coordinate of a lost ray event in cartesian coordinates */
	cdoub **camp; /**< the complex amplitude at each bin of the screen */
	uint res_polar, /**< the polar resolution of the binarized sphere/screen */
		 res_azim, /**< the azimuthal resolution of the binarized sphere/screen */
		 nbin; /**< the number of the bins */
	index2 *idx; /**< this array stores the correct indexing of the bins to address them easily */
	ulong tir; /**< the count of tir events */
	coordinate_sys screen_hit_coor_sys; /**< the coordinate system being used to store the screen hit events (and not the other events, e.g. tir) */
	gen_ray_info global_info; /**< some informations about the rays */
}bin_hit_screen; /**< a structure of the binarized sphere which acts as the detection screen */

typedef struct ray_struct
{
	line3 v; /**< the vector which is followed by the ray */
	double lam, /**< the wavenlength of the ray */
		   oamp, /**< the 'o' amplitude */
		   pamp, /**< the 'p' amplitude */
		   ophase, /**< the 'o' phase */
		   pphase, /**< the 'p' phase */
		   oint, /**< the 'o' (or TE) intensity */
		   pint, /**< the 'p' intensity */
		   opol[3], /**< the 'o' (normalized) polarization-vector of the ray */
		   ppol[3], /**< the 'p' (normalized) polarization-vector of the ray. it is always orthogonal to 'o' and the directional vector in a right-handed sense. */
		   mu_i, /**< the permeability of the medium the ray is in */
		   travel; /**< the optical path traveled by the ray */
	cdoub n_i; /**< the refractive index of the medium the ray is in */
	uint tir, /**< the number of tir events of this ray */
		 trans_child, /**< the number of transmission which led to the creation of this ray. warning: when used this in a loop, it has to be used *with* the last number. */
		 hits; /**< the number of times a particle was hit by this child */
	char info[INFOLENGTH+1]; /**< an information about this ray */
}ray; /**< the structure of the ray being traced through the system */

typedef struct glray_s_struct
{
	line3 v; /**< the vector which is initially followed by the ray */
	point3 *trace; /**< an array of points used for ray tracing */
	double oamp, /**< the 'o' (or TE) amplitude */
		   pamp, /**< the 'p' (or TM) amplitude */
		   ophase, /**< the 'o' phase */
		   pphase, /**< the 'p' phase */
		   oint, /**< the 'o' (or TE) intensity */
		   pint, /**< the 'p' intensity */
		   opol[3], /**< the initial 'o' (normalized) polarization-vector of the ray */
		   ppol[3]; /**< the initial 'p' (normalized) polarization-vector of the ray */
	const double *rgba[4]; /**< the color to be used for drawing the ray */
	uint n_trace, /**< number of trace points */
		 n_child, /**< number of ray childs a.k.a. segments */
		 trace_len, /**< allocated space for the trace array */
		 child_len, /**< allocated space for the child array */
		 *child; /**< array to store the number of childs of one ray */
}glray_s; /**< a single ray used for the OpenGL viewer */

typedef struct glray_struct
{
	glray_s *glrs; /**< an array of glray_s variables */
	uint n_glrs; /**< the number of rays inside this structure */
}glray; /**< a container of glray_s used for the OpenGL viewer */

typedef struct sphrcl_prtcl_struct
{
	sphere3 s; /**< the location and size of the particle */
	cdoub n; /**< the refractive index of the particle */
	double mu; /**< the permeability of the particle */
	uint no; /**< the number of reference of the particle */
	ulong hits, /**< the time this particle has been hit */
		  exhstds; /**< the number of rays which have been died inside of the particle */
	char info[L_INFOLENGTH + 1]; /**< an information */
}sphrcl_prtcl; /**< a structure containing the data of a spherical particle */

typedef struct vertex3_struct
{
	double n[3],
		   x[3];
}vertex3; /**< a structure containing the information about a vertex in R^3 to be used in the OpenGL viewer */

typedef struct patch3_struct
{
	vertex3 *vt; /**< an array of vertices */
	uint n_vt, /**< the number of vertices */
		 gl_primitive; /**< the render method to be used, e.g. triangle stripes */
	double rgba[4]; /**< the normal color of all patches */
}patch3; /**< a structure to describe a patch of an object used in the OpenGL viewer, e.g. the bins of a sphere */

typedef struct colorval_struct
{
	double rgba[4],
		   val;
}colorval; /**< a structure which contains a double value and its map into rgba space */

typedef struct colorbox_struct
{
	colorval *cval; /**< an array of colorval variables to store the lines that fill up the colorbox */
	double max, /**< the maximum double value of the colorscale */
		   min; /**< the minimum double value of the colorscale */
	uint ncval; /**< the number of the lines in the cval array */
}colorbox; /**< a structure to describe the colorbox in the OpenGL viewer */

typedef struct boundingbox_struct
{
	point3 rect[8]; /**< the vertex points of the rectangular box */
	char xwidth[INFOLENGTH + 1], /**< the width to be printed at the box */
		 yheight[INFOLENGTH + 1], /**< the height to be printed at the box */
		 zdepth[INFOLENGTH + 1]; /**< the depth to be printed at the box */
}boundingbox; /**< a structure to be used in the OpenGL viewer which is used to circumscribe all the particles */
