#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define __USE_XOPEN_EXTENDED /**< Defined to access usleep with the C99 standard on Linux, see unistd.h */
#include <unistd.h>
#include <tgmath.h>
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
 * - corrected and improved the Sellmeier equations
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
 * August 2015
 * - using git, stopping rcs
 * May 2017
 * - cleaned and re-formatted code and comments for better readability
 * June 2018
 * - updated const correctness and static code checks
 */

#define RELEASE_BUILD 0 /**< Hides or prints certain information important for a public usage */
#if RELEASE_BUILD
 #define MEMORY_COUNTER 1 /**< Controls the summation and printing of allocated memory */
#else
 #define MEMORY_COUNTER 0
#endif
#define DEBUG_MOTION 0 /**< Controls the printing of additional information in the viewer */
#define PRINT_PROG_OF_RAYS 1 /**< Controls if the progress of the ray tracing is printed to the screen */
#define PARALLEL_PROCESSING 0 /**< Controls the usage of OpenMP */
#define use_restrict 1 /**< Controls if restricted pointers are used */
#define STOP_AT_ERR 0 /**< Controls if 'Enter' is required at an error_msg call or not */

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

#define DBL_EPSILON2 (DBL_EPSILON * 2.) /**< A threshold value equalling the maximum absolute error in a dot3 function */
#define DBL_EPSILON4 (DBL_EPSILON * 4.) /**< A threshold value used in different circumstances instead of the machine epsilon */
#define DBL_EPSILON6 (DBL_EPSILON * 6.) /**< A threshold value used in different circumstances instead of the machine epsilon */
#define DBL_EPSILON10 (DBL_EPSILON * 1.e1) /**< Another threshold value used in different circumstances instead of the machine epsilon */
#define DBL_EPSILON20 (DBL_EPSILON * 8.e1) /**< Another threshold value used in different circumstances instead of the machine epsilon */
#define MINI_AMPLITUDE 1.e-10 /**< The minimum amplitude of a ray after which the ray is said to be exhausted */
#define MINI_INTENSITY 1.e-15 /**< The minimum intensity of a ray after which the ray is said to be exhausted */
#define MINI_TRANSMISSION_COEF 1.e-15 /**< The minimum transmission coefficient. If this value is undershoot, the transmitted ray is not traced. */
#define MINI_REFLECTION_COEF 1.e-15 /**< The minimum reflection coefficient. If this value is undershoot, the reflected ray is not traced. */
#define MAX_CHILD 50u /**< The maximum number of child rays. If this value is overshoot, the next child ray will not be traced. */
#define EPS_RAY_PROP_LEN 1.e-10 /**< The geometric distance which is added to a ray after a hit to propagate it from the point of intersection */

#define MAX_NUMBER_OF_RAYS (unsigned int)UINT_MAX
#define INFOLENGTH 64u
#define L_INFOLENGTH 256u
#define FILENAME_MAX1 (uint)(FILENAME_MAX + 1)

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
}ray_end_state; /**< Describes what happened to a ray at the end of tracing */

typedef enum bin_hit_print_type_enum
{
    E_FIELD_AMPLITUDE,
    MOD_E_FIELD_AMPLITUDE,
    INTENSITY,
    POLARISATION_DENSITY
}bin_hit_print_type; /**< Defines what should be printed with the OpenGL viewer */

typedef enum prog_ray_order_enum
{
    INIT_PROG,
    COUNT_PROG,
    PRINT_PROG,
    CLEAR_OUT
}prog_ray_order; /**< An option identifier of what shall be printed when calling prog_of_rays */

typedef enum coordinate_sys_enum
{
    CARTESIAN_CS,
    SPHERICAL_CS,
    MIXED_CS
}coordinate_sys; /**< An identifier of the coordinate system */

typedef enum incidence_enum
{
    NONE,
    OBLIQUE,
    VERTICAL,
    RIGHT_ANGLED
}incidence; /**< An identifier of the type of incidence */

typedef enum draw_opt_enum
{
    ALLOC_DATA,
    COPY_RAYS,
    GEN_LISTS,
    DRAW_EM_ALL,
    DRAW_SOME,
    FREE_EM_ALL
}draw_opt; /**< An option identifier to be used for the drawing function to specify what has to be done */

typedef enum colorfun_enum
{
    SINCOS_MAP,
    VISIBLE_LIGHT
}colorfun; /**< An identifier of the function used to colour the output */

typedef enum trackball_mode_enum
{
    RESET = 0,
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
}index2; /**< Used to address the bins of the detections sphere easily */

typedef struct point3_struct
{
    double x[3];
}point3; /**< A point in R^3 */

typedef struct line3_struct
{
    double o[3], /**< Origin */
           r[3], /**< Direction */
           l; /**< Scalar multiplier */
}line3; /**< A vector */

typedef struct plane3_struct
{
    double o[3], /**< A point in the plane */
           n[3]; /**< The normalized normal vector */
}plane3; /**< A plane in R^3 */

typedef struct sphere3_struct
{
    double o[3], /**< The center of the sphere */
           r; /**< The radius */
}sphere3; /**< A sphere in R^3 */

typedef struct intrsec_struct
{
    double p[3], /**< The point of the intersection */
           normal[3], /**< The (normalized) normal vector with respect to the surface, pointing outwards */
           angl, /**< The angle between the surface normal and the ray vector */
           cangl, /**< The cosine of angl */
           mu_f; /**< The permeability of the medium the ray is going to enter */
    cdoub n_f; /**< The refractive index of the medium the ray is going to enter */
    incidence incdnc; /**< Identifier of the incidence */
}intrsec; /**< A structure containing information about an intersection between a ray and a particle in R^3 */

typedef struct hit_screen_struct
{
    double p[3], /**< The point of the intersection */
           lam, /**< The wavelength of the ray */
           cos_incdnc, /**< The cosine of the angle of incidence */
           opol[3], /**< The 'o' (normalized) polarization-vector of the ray */
           ppol[3]; /**< The 'p' (normalized) polarization-vector of the ray */
    cdoub coamp, /**< The complex amplitude of the 'o' polarization */
          cpamp, /**< The complex amplitude of the 'p' polarization */
          oint, /**< The intensity of the 'o' polarization */
          pint; /**< The intensity of the 'p' polarization */
    uint tir; /**< The number of TIR events of the ray */
    char info[INFOLENGTH + 1]; /**< An information about this event */
    ray_end_state state; /**< The precise end state of the ray */
}hit_screen; /**< A structure describe the last resort of the ray (not necessarily the screen) */

typedef struct gen_ray_info_struct
{
    ulong count_gen, /**< The number of generated rays */
          count_hit, /**< The number of rays that hit the screen */
          count_lost, /**< The number of rays which have been lost, i.e. found outside of the screen */
          count_exhstd; /**< The number of rays which have been exhausted during their flight */
    char info[L_INFOLENGTH + 1]; /**< A general information */
}gen_ray_info; /**< A general information about the rays traced through the system */

typedef struct bin_hit_screen_struct
{
    double lam, /**< The wavelength of the rays */
           rad, /**< The radius of the sphere */
           **amp, /**< The amplitude at each bin of the screen */
           **mod_amp, /**< The modulus of the amplitude at each bin of the screen */
           **ray_int, /**< The intensity at each bin of the screen */
           **pol_dist, /**< The polarization distribution at each bin of the screen */
           **pol_dens, /**< The polarization density at each bin of the screen */
           amp_max, /**< The maximum amplitude detected */
           amp_sum, /**< The sum of the amplitudes detected */
           mod_amp_max, /**< The maximum modulus of the amplitude detected */
           mod_amp_sum, /**< The sum of the modulus of the amplitudes detected */
           ray_int_max, /**< The maximum of the intensities detected */
           ray_int_sum, /**< The sum of the intensities detected */
           pol_max, /**< The maximum polarization ???density?? detected */
           pol_sum; /**< The sum of the polarization ???densities?? detected */
    point3 *tir_coor, /**< The coordinate of a TIR event in Cartesian coordinates */
           *exhstd_coor, /**< The coordinate of an exhausted ray event in Cartesian coordinates */
           *lost_coor; /**< The coordinate of a lost ray event in Cartesian coordinates */
    cdoub **camp; /**< The complex amplitude at each bin of the screen */
    uint res_polar, /**< The polar resolution of the binarised sphere/screen */
         res_azim, /**< The azimuthal resolution of the binarised sphere/screen */
         nbin; /**< The number of the bins */
    index2 *idx; /**< This array stores the correct indexing of the bins to address them easily */
    ulong tir; /**< The count of TIR events */
    coordinate_sys screen_hit_coor_sys; /**< The coordinate system being used to store the screen hit events (and not the other events, e.g. TIR) */
    gen_ray_info global_info; /**< Some informations about the rays */
}bin_hit_screen; /**< A structure of the binarised sphere which acts as the detection screen */

typedef struct ray_struct
{
    line3 v; /**< The vector which is followed by the ray */
    double lam, /**< The wavelength of the ray */
           oamp, /**< The 'o' amplitude */
           pamp, /**< The 'p' amplitude */
           ophase, /**< The 'o' phase */
           pphase, /**< The 'p' phase */
           oint, /**< The 'o' (or TE) intensity */
           pint, /**< The 'p' intensity */
           opol[3], /**< The 'o' (normalized) polarization-vector of the ray */
           ppol[3], /**< The 'p' (normalized) polarization-vector of the ray. It is always orthogonal to 'o' and the directional vector in a right-handed sense. */
           mu_i, /**< The permeability of the medium the ray is in */
           travel; /**< The optical path travelled by the ray */
    cdoub n_i; /**< The refractive index of the medium the ray is in */
    uint tir, /**< The number of TIR events of this ray */
         trans_child, /**< The number of transmission which led to the creation of this ray. Warning: when using this in a loop, it has to be used *with* the last number. */
         hits; /**< The number of times a particle was hit by this child */
    char info[INFOLENGTH+1]; /**< An information about this ray */
}ray; /**< The structure of the ray being traced through the system */

typedef struct glray_s_struct
{
    line3 v; /**< The vector which is initially followed by the ray */
    point3 *trace; /**< An array of points used for ray tracing */
    double oamp, /**< The 'o' (or TE) amplitude */
           pamp, /**< The 'p' (or TM) amplitude */
           ophase, /**< The 'o' phase */
           pphase, /**< The 'p' phase */
           oint, /**< The 'o' (or TE) intensity */
           pint, /**< The 'p' intensity */
           opol[3], /**< The initial 'o' (normalized) polarization-vector of the ray */
           ppol[3]; /**< The initial 'p' (normalized) polarization-vector of the ray */
    const double *rgba[4]; /**< The colour to be used for drawing the ray */
    uint n_trace, /**< Number of trace points */
         n_child, /**< Number of ray children a.k.a. segments */
         trace_len, /**< Allocated space for the trace array */
         child_len, /**< Allocated space for the child array */
         *child; /**< Array to store the number of children of one ray */
}glray_s; /**< A single ray used for the OpenGL viewer */

typedef struct glray_struct
{
    glray_s *glrs; /**< An array of glray_s variables */
    uint n_glrs; /**< The number of rays inside this structure */
}glray; /**< A container of glray_s used for the OpenGL viewer */

typedef struct sphrcl_prtcl_struct
{
    sphere3 s; /**< The location and size of the particle */
    cdoub n; /**< The refractive index of the particle */
    double mu; /**< The permeability of the particle */
    uint no; /**< The number of reference of the particle */
    ulong hits, /**< The time this particle has been hit */
          exhstds; /**< The number of rays which have been died inside of the particle */
    char info[L_INFOLENGTH + 1]; /**< An information */
}sphrcl_prtcl; /**< A structure containing the data of a spherical particle */

typedef struct vertex3_struct
{
    double n[3],
           x[3];
}vertex3; /**< A structure containing the information about a vertex in R^3 to be used in the OpenGL viewer */

typedef struct patch3_struct
{
    vertex3 *vt; /**< An array of vertices */
    uint n_vt, /**< The number of vertices */
         gl_primitive; /**< The render method to be used, e.g. triangle stripes */
    double rgba[4]; /**< The normal colour of all patches */
}patch3; /**< A structure to describe a patch of an object used in the OpenGL viewer, e.g. the bins of a sphere */

typedef struct colorval_struct
{
    double rgba[4],
           val;
}colorval; /**< A structure which contains a double value and its map into rgba space */

typedef struct colorbox_struct
{
    colorval *cval; /**< An array of colorval variables to store the lines that fill up the colorbox */
    double max, /**< The maximum double value of the colorscale */
           min; /**< The minimum double value of the colorscale */
    uint ncval; /**< The number of the lines in the cval array */
}colorbox; /**< A structure to describe the colorbox in the OpenGL viewer */

typedef struct boundingbox_struct
{
    point3 rect[8]; /**< The vertex points of the rectangular box */
    char xwidth[INFOLENGTH + 1], /**< The width to be printed at the box */
         yheight[INFOLENGTH + 1], /**< The height to be printed at the box */
         zdepth[INFOLENGTH + 1]; /**< The depth to be printed at the box */
}boundingbox; /**< A structure to be used in the OpenGL viewer which is used to circumscribe all the particles */
