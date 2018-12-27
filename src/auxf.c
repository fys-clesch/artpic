#include "artpic.h"
#include "lina.h"
#include "alloc.h"
#include "ray_aux.h"
#include "color.h"
#include "rot.h"
#include "msg.h"
#include "auxf.h"

/** \brief Gets the maximum number of digits of a double.
 *
 * \param x double The input double.
 * \return uint The maximum number of digits.
 *
 */
uint get_nodd(double x)
{
    char strf[128],
         stre[128];
    snprintf(strf, 128, "%f", x);
    snprintf(stre, 128, "%e", x);
    uint f = strlen(strf),
         e = strlen(stre);
    return MAXOF(f, e);
}

/** \brief Gets the maximum number of digits of an integer.
 *
 * \param x int The input int.
 * \return uint The maximum number of digits.
 *
 */
uint get_nodi(int x)
{
    char str[128];
    snprintf(str, 128, "%i", x);
    return strlen(str);
}

/** \brief Gets the maximum number of digits of an unsigned integer.
 *
 * \param x uint The input uint.
 * \return uint The maximum number of digits.
 *
 */
uint get_nodui(uint x)
{
    char str[128];
    snprintf(str, 128, "%u", x);
    return strlen(str);
}

/** \brief Copy a 3*double array.
 *
 * \param dest double *res_pt Pointer to the destination.
 * \param source const double *res_pt Pointer to the source.
 * \return void
 *
 * Source and destination are not allowed to overlap - this might cause problems when called in other functions.
 */
void cp3(double *res_pt dest, const double *res_pt source)
{
    static const uint i = 3 * sizeof(double);
    memcpy(dest, source, i);
}

/** \brief An integer exponentiation function.
 *
 * \param x int The base.
 * \param n uint The power.
 * \return int The result.
 *
 */
int int_pow(int x, uint n)
{
    int pow = 1;
    if(n)
        for(;;)
        {
            if(n & 1) pow *= x; /**< n is odd. */
            if(n >>= 1) x *= x;
            else break;
        }
    return pow;
}

/** \brief Returns the absolute of the real part of a complex number.
 *
 * \param z const cdoub A complex number.
 * \return double The absolute real.
 *
 */
double cabs_real(const cdoub z)
{
    return fabs(creal(z));
}

/** \brief Returns the absolute of the imaginary part of a complex number.
 *
 * \param z const cdoub A complex number.
 * \return double The absolute imaginary.
 *
 */
double cabs_imag(const cdoub z)
{
    return fabs(cimag(z));
}

/** \brief Computes the case-sensitive arctan function in the interval [0,2*pi).
 *
 * \param x double The numerator.
 * \param y const double The denominator.
 * \return double The return value.
 *
 */
double atan2_up(double x, const double y)
{
    x = atan2(x, y);
    if(x < 0.)
        x += M_PI2;
    return x;
}

/** \brief Compares two pointers to an int.
 *
 * \param p1 const void* The first pointer to an int.
 * \param p2 const void* The second pointer to an int.
 * \return int States whether the pointers are equal or unequal.
 *
 * In case that the first int is greater than the second one, return 1, if both are equal, return 0, otherwise -1.
 */
int comp_int(const void *p1, const void *p2)
{
    if(*(int *)p1 < * (int *)p2) return -1;
    else if(*(int *)p1 > *(int *)p2) return 1;
    return 0;
}

/** \brief Compares two pointers to a double.
 *
 * \param p1 const void* The first pointer to a double.
 * \param p2 const void* The second pointer to a double.
 * \return int States whether the pointers are equal or unequal.
 *
 * In case that the first double is greater than the second one, return 1, if both are equal, return 0, otherwise -1.
 */
int comp_double(const void *p1, const void *p2)
{
    if(*(double *)p1 < * (double *)p2) return -1;
    else if(*(double *)p1 > *(double *)p2) return 1;
    return 0;
}

/** \brief Finds the maximum value in a double array.
 *
 * \param m const double *res_pt The array.
 * \param xm uint *res_pt The first index where the maximum value is found.
 * \param ym uint *res_pt The second index where the maximum value is found.
 * \param row const uint The range of the first index.
 * \param col const uint The range of the second index.
 * \return double The maximum value.
 *
 * Only the first occurrence of the maximum value is returned.
 * row and col are used to help indexing.
 */
double find_max(const double *res_pt m, uint *res_pt xm, uint *res_pt ym, const uint row, const uint col)
{
    double t1 = -DBL_MAX;
    uint x;
    for(x = 0; x < row; x++)
    {
        uint t2 = x * col, y;
        for(y = 0; y < col; y++)
            if(m[t2 + y] > t1)
            {
                t1 = m[t2 + y];
                *xm = x;
                *ym = y;
            }
    }
    return t1;
}

/** \brief Finds the minimum value in a double array
 *
 * \param m const double *res_pt The array.
 * \param xm uint *res_pt The first index where the minimum value is found.
 * \param ym uint *res_pt The first index where the minimum value is found.
 * \param row const uint The range of the first index.
 * \param col const uint The range of the first index.
 * \return double The minimum value.
 *
 * Only the first occurrence of the minimum value is returned.
 * Row and col are used to help indexing.
 */
double find_min(const double *res_pt m, uint *res_pt xm, uint *res_pt ym, const uint row, const uint col)
{
    double t1 = DBL_MAX;
    uint x;
    for(x = 0; x < row; x++)
    {
        uint t2 = x * col, y;
        for(y = 0; y < col; y++)
            if(m[t2 + y] < t1)
            {
                t1 = m[t2 + y];
                *xm = x;
                *ym = y;
            }
    }
    return t1;
}

/** \brief Solve ax^2+bx+c=0 for x.
 *
 * \param a const double The parameter of x^2.
 * \param b double The parameter of x^1.
 * \param c const double The parameter of x^0.
 * \param r[] double The solution.
 * \param warn_complex const uchar If(i) then print a warning to the screen in case of complex solution.
 * \return uchar If two solutions are found, return 2, if one solution is found, return 1, if a complex solution is found, return 0.
 *
 */
uchar solve_pq(const double a, double b, const double c, double r[2], const uchar warn_complex)
{
    double t = b * b - 4. * a * c;
    if(isgreater(t, 0.))
    {
#ifdef FP_FAST_FMA
        b = -.5 * fma(sign(b), sqrt(t), b);
#else
        b = -.5 * (b + sign(b) * sqrt(t));
#endif
        r[0] = b / a;
        r[1] = c / b;
        return 2;
    }
    else if(isless(t, 0.))
    {
        if(warn_complex)
            fprintf(stdout, "! solution in solve_pq found to be complex\n");
        r[0] = -.5 * sign(b) * sqrt(t);
        r[1] = -.5 * b;
        return 0;
    }
    else
    {
        r[0] = r[1] = -2. * c / b;
        return 1;
    }
}

/** \brief Gets the sign of a double.
 *
 * \param x const double The float to get the sign from.
 * \return double The sign in front of unity.
 *
 * Returns 0. if x is 0.
 */
double sign(const double x)
{
    return (x > 0.) ? 1. : ((x < 0.) ? -1. : 0.);
}

/** \brief Takes an array of hit_screen variables and sorts it into a bin_hit_screen variable.
 *
 * \param hso const hit_screen *res_pt The hit_screen array.
 * \param no const uint The length of the hit_screen array.
 * \param bhs bin_hit_screen *res_pt The variable used to store the hits.
 * \param sumup const uchar A flag used to control the summation of the electric field in each bin.
 * \return void
 *
 * All the detections in the array have to be of the same wavelength (tested for security reasons) and a Cartesian coordinate system (not tested).
 * Also updates the variable global_info of 'bhs'.
 */
void sort_detection(const hit_screen *res_pt hso, const uint no, bin_hit_screen *res_pt bhs, const uchar sumup)
{
    hit_screen *hsc = alloc_hit_screen(no); /**< Copy of hs due to coordinate- and some sign changes. */
    memcpy(hsc, hso, no * sizeof(hit_screen));
    const double ra = M_PI2 / ((double)(*bhs).res_azim),
                 rp = M_PI / ((double)(*bhs).res_polar);
    static double c_lam = 0.;
    static uint i_exhstd = 0,
                i_lost = 0;
    uint i;
    (*bhs).tir += hsc[0].tir;
    if(sumup) /**< In this case, a hit_screen state may be equal to EMPTY_STATE, because not all rays allocated had been traced... */
    {
        double t1 = 0xDEAD;
        for(i = 0; i < no; i++)
            if(hsc[i].state != EMPTY_STATE) /**< Find a valid wavelength: */
            {
                t1 = hsc[i].lam;
                if(c_lam == 0.)
                    c_lam = t1; /**< In this case, all rays have hit the screen but no particle, because this function shouldn't be called before. */
                break;
            }
        for(i = 1; i < no; i++) /**< Test if there's another wavelength in the array: */
            if(hsc[i].state != EMPTY_STATE && t1 != hsc[i].lam)
            {
                fprintf(stderr, "hsc[%u].lam = %g\n", i, hsc[i].lam);
                error_msg("mixed up wavelengths. good bye.", ERR_ARG);
                free(hsc);
                exit(EXIT_FAILURE);
            }
            else
                (*bhs).tir += hsc[i].tir;
        (*bhs).tir += hsc[0].tir;
        (*bhs).lam = c_lam;
        c_lam = 0.; /**< Reset the constant wavelength for a new run, e.g. for a second bundle of rays. */
        memcpy(&(*bhs).global_info, &GLOBAL_RAY_INFO, sizeof(gen_ray_info));
    }
    else /**< ...here not, because this function is called after a certain threshold of hits was reached. */
    {
        if(c_lam == 0.)
            c_lam = hsc[0].lam;
        for(i = 1; i < no; i++)
            if(c_lam != hsc[i].lam)
            {
                fprintf(stderr, "hsc[%u].lam = %g\n", i, hsc[i].lam);
                error_msg("mixed up wavelengths. good bye.", ERR_ARG);
                free(hsc);
                exit(EXIT_FAILURE);
            }
            else
                (*bhs).tir += hsc[i].tir;
        (*bhs).tir += hsc[0].tir;
    }
    for(i = 0; i < no; i++) /**< Scan through the different rays and store them according to their state: */
        if(hsc[i].state == REGULAR_HIT_STATE || hsc[i].state == FIRSTTIME_HIT_STATE)
            c_to_s3_ip(hsc[i].p); /**< Transform the coordinate system to a spherical one. */
        else if(hsc[i].state == EXCEPTION_STATE || hsc[i].state == LOST_RAY_STATE)
        {
            if(i_lost == ALLOCATED_LOST_EVENT_MEMORY) /**< Check if there are already too many lost rays stored. */
            {
                (*bhs).lost_coor = realloc_point3((*bhs).lost_coor, i_lost, ALLOCATED_LOST_EVENT_MEMORY);
                ALLOCATED_LOST_EVENT_MEMORY += ALLOCATED_LOST_EVENT_MEMORY;
            }
            cp3((*bhs).lost_coor[i_lost].x, hsc[i].p);
            i_lost++;
        }
        else if(hsc[i].state == EXHAUSTED_RAY_STATE)
        {
            if(i_exhstd == ALLOCATED_EXHAUSTED_EVENT_MEMORY) /**< Check if there are already too many exhausted rays stored. */
            {
                (*bhs).exhstd_coor = realloc_point3((*bhs).exhstd_coor, i_exhstd, ALLOCATED_EXHAUSTED_EVENT_MEMORY);
                ALLOCATED_EXHAUSTED_EVENT_MEMORY += ALLOCATED_EXHAUSTED_EVENT_MEMORY;
            }
            cp3((*bhs).exhstd_coor[i_exhstd].x, hsc[i].p);
            i_exhstd++;
        }
    (*bhs).screen_hit_coor_sys = SPHERICAL_CS;
    for(i = 0; i < no; i++) /**< Sort the screen hits into the bins of the detection screen/sphere: */
        if(hsc[i].state == REGULAR_HIT_STATE || hsc[i].state == FIRSTTIME_HIT_STATE)
        {
            double tp = hsc[i].p[1] / rp,
                   ta = hsc[i].p[2] / ra;
            ulong tpi = lrint(floor(tp)),
                  tai = lrint(floor(ta));
            if(tpi == 0 || tpi == (*bhs).res_polar - 1 || tai == (*bhs).res_azim)
                tai = 0;
            assert(tai < (*bhs).res_azim && tpi < (*bhs).res_polar);
            uint j,
                 tai3 = tai * 3; /**< There are 3 entries in each bin, to step from one to another, multiply the address by 3. */
            /**< Compute and sum the complex e-field at a bin, weighted by the polarisation. thus Cartesian coordinates have to be used: */
            for(j = 0; j < 3; j++)
                (*bhs).camp[tpi][tai3 + j] += hso[i].coamp * hso[i].opol[j] +
                                              hso[i].cpamp * hso[i].ppol[j];
            /**< Compute and sum the intensity at each bin: */
            assert(hso[i].cos_incdnc != 0xDEAD && hso[i].cos_incdnc < 0.);
            (*bhs).ray_int[tpi][tai] += ((hso[i].oint + hso[i].pint) * fabs(hso[i].cos_incdnc));
            /**< Compute distribution of polarisation vectors, weighted by the e-field amplitude: */
            //        if(hsc[i].opol[0]<0.) hsc[i].opol[0]=-hsc[i].opol[0]; /**< Map opol[0] to positive half sphere */
            //        if(hsc[i].opol[0]<DBL_EPSILON) hsc[i].opol[1]=-hsc[i].opol[1];
            //        if(hsc[i].ppol[0]<0.) hsc[i].ppol[0]=-hsc[i].ppol[0]; /**< Map ppol[0] ... */
            //        if(hsc[i].ppol[0]<DBL_EPSILON) hsc[i].ppol[1]=-hsc[i].ppol[1];
            for(j = 0; j < 3; j++)
                hsc[i].opol[j] = fabs(hsc[i].opol[j]);
            for(j = 0; j < 3; j++)
                hsc[i].ppol[j] = fabs(hsc[i].ppol[j]);
            double t1 = cabs(hsc[i].coamp);
            for(j = 0; j < 3; j++)
                (*bhs).pol_dist[tpi][tai3 + j] += t1 * hsc[i].opol[j];
            t1 = cabs(hsc[i].cpamp);
            for(j = 0; j < 3; j++)
                (*bhs).pol_dist[tpi][tai3 + j] += t1 * hsc[i].ppol[j];
        }
    if(sumup)
    {
        double polnorm,
               s_a = 0., s_p = 0., s_ma = 0., s_i = 0.,
               m_a = 0., m_p = 0., m_ma = 0., m_i = 0.;
        for(i = 0; i < (*bhs).nbin; i++)
        {
            uint j = (*bhs).idx[i].ia,
                 k = (*bhs).idx[i].ib,
                 l0 = 3 * k,
                 l1 = l0 + 1,
                 l2 = l0 + 2;
            /**< Compute the modulus of the amplitude by taking the modulus of the complex vector: */
            (*bhs).mod_amp[j][k] = POW2(creal((*bhs).camp[j][l0])) + POW2(cimag((*bhs).camp[j][l0])) +
                                   POW2(creal((*bhs).camp[j][l1])) + POW2(cimag((*bhs).camp[j][l1])) +
                                   POW2(creal((*bhs).camp[j][l2])) + POW2(cimag((*bhs).camp[j][l2]));
            if((*bhs).mod_amp[j][k] > m_ma)
                m_ma = (*bhs).mod_amp[j][k];
            s_ma += (*bhs).mod_amp[j][k];
            /**< Compute the scalar amplitude of the e-field: */
            (*bhs).amp[j][k] = sqrt(POW2(creal((*bhs).camp[j][l0])) +
                                    POW2(creal((*bhs).camp[j][l1])) +
                                    POW2(creal((*bhs).camp[j][l2])));
            if((*bhs).amp[j][k] > m_a)
                m_a = (*bhs).amp[j][k];
            s_a += (*bhs).amp[j][k];
            /**< Compute the sum and maximum of the intensity: */
            if((*bhs).ray_int[j][k] > m_i)
                m_i = (*bhs).ray_int[j][k];
            s_i += (*bhs).ray_int[j][k];
            /**< Norm the distribution of the polarisation in one bin: */
            double t1 = MAXOF(fabs((*bhs).pol_dist[j][l0]), MAXOF(fabs((*bhs).pol_dist[j][l1]), fabs((*bhs).pol_dist[j][l2])));
            if(t1 != 0.)
            {
                (*bhs).pol_dist[j][l0] /= t1;
                (*bhs).pol_dist[j][l1] /= t1;
                (*bhs).pol_dist[j][l2] /= t1;
                /**< Compute the mean deviation of the polarisation: */
                polnorm = POW2((*bhs).pol_dist[j][l0]) + POW2((*bhs).pol_dist[j][l1]) + POW2((*bhs).pol_dist[j][l2]);
                (*bhs).pol_dens[j][k] = (polnorm > DBL_MIN ? sqrt(polnorm) : 0.);
                if((*bhs).pol_dens[j][k] > m_p)
                    m_p = (*bhs).pol_dens[j][k];
                s_p += (*bhs).pol_dens[j][k];
            }
        }
        (*bhs).amp_max = m_a;
        (*bhs).amp_sum = s_a;
        (*bhs).mod_amp_max = m_ma;
        (*bhs).mod_amp_sum = s_ma;
        (*bhs).ray_int_max = m_i;
        (*bhs).ray_int_sum = s_i;
        (*bhs).pol_max = m_p;
        (*bhs).pol_sum = s_p;
        i_exhstd = i_lost = 0; /**< Reset these values for a new run. */
    }
    free(hsc);
}

/** \brief Sets a point3 variable.
 *
 * \param p point3* The pointer to a point3 variable.
 * \param x0 const double The first coordinate.
 * \param x1 const double The second coordinate.
 * \param x2 const double The third coordinate.
 * \return void
 *
 */
void set_point3(point3 *p, const double x0, const double x1, const double x2)
{
    (*p) = (point3){.x = {x0, x1, x2}};
}

/** \brief Sets a sphere3 variable.
 *
 * \param p sphere3* The pointer to a sphere3 variable.
 * \param x0 const double The first coordinate.
 * \param x1 const double The second coordinate.
 * \param x2 const double The third coordinate.
 * \param r const double The radius.
 * \return void
 *
 */
void set_sphere3(sphere3 *p, const double x0, const double x1, const double x2, const double r)
{
    (*p) = (sphere3){.o = {x0, x1, x2}, .r = r};
}

/** \brief Sets a sphrcl_prtcl variable.
 *
 * \param p sphrcl_prtcl* The pointer to a sphrcl_prtcl variable.
 * \param n const cdoub The refractive index of the particle.
 * \param mu const double The permeability of the particle.
 * \param x0 const double The first coordinate.
 * \param x1 const double The second coordinate.
 * \param x2 const double The third coordinate.
 * \param r const double The radius of the particle.
 * \return void
 *
 */
void set_sphrcl_prtcl(sphrcl_prtcl *p, const cdoub n, const double mu, const double x0, const double x1, const double x2, const double r)
{
    (*p) = (sphrcl_prtcl){
        .s = (sphere3) {.o = {x0, x1, x2}, .r = r}, .n = n, .mu = mu};
    snprintf((*p).info, L_INFOLENGTH, "no %u from (%g, %g, %g). n = %g + %gi. mu = %g. rad = %g.",
             (*p).no, x0, x1, x2, creal(n), cimag(n), mu, r);
}

/** \brief Sets an index2 variable for a proper indexing of a sphere.
 *
 * \param indx index2* The output variable via pointer.
 * \param pol const uint The polar resolution of the sphere.
 * \param azi const uint The azimuthal resolution of the sphere.
 * \return void
 *
 */
void set_index2_bin_sphere3(index2 *indx, const uint pol, const uint azi)
{
    uint i, j, k;
    for(i = 0, k = 0; i < pol; i++)
    {
        if(i != 0 && i != pol - 1)
            for(j = 0; j < azi; j++, k++)
            {
                indx[k].ia = i;
                indx[k].ib = j;
            }
        else
        {
            indx[k].ia = i;
            indx[k].ib = 0;
            k++;
        }
    }
    assert(k == (pol - 2)*azi + 2);
}

/** \brief Copy an array of rays to a glray_s array to display them in the OpenGL viewer.
 *
 * \param rs const ray *const res_pt The array of rays used in the calculation.
 * \param glrs glray_s *res_pt The arrays of rays used in the viewer.
 * \param nors const uint The number of rays.
 * \return void
 *
 * Just to be used for the initial state of rays, i.e. before the actual tracing.
 */
void copy_ray_to_glray_s(const ray *const res_pt rs, glray_s *res_pt glrs, const uint nors)
{
    uint i;
    for(i = 0; i < nors; i++)
    {
        glrs[i] = (glray_s){
            .v = rs[i].v,
            .oamp = rs[i].oamp, .pamp = rs[i].pamp,
            .ophase = rs[i].ophase, .pphase = rs[i].pphase,
            .oint = rs[i].oint, .pint = rs[i].pint,
            .opol = {rs[i].opol[0], rs[i].opol[1], rs[i].opol[2]},
            .ppol = {rs[i].ppol[0], rs[i].ppol[1], rs[i].ppol[2]}};
        add3(glrs[i].v.o, glrs[i].v.r, glrs[i].v.r);
        ulong t1 = lrint(floor(rs[i].lam * 1e6)); /**< Using a uint cast brings at least valgrind into trouble. */
        assert(t1 >= 380);
        t1 -= 380;
        assert(t1 < 401);
        uint j;
        for(j = 0; j < 3; j++)
            glrs[i].rgba[j] = &LIGHTTOCOLOR[t1][j];
    }
}

/** \brief Writes the data from the current ray into a glray variable for visualization.
 *
 * \param const ray *const res_pt r The current ray to be traced.
 * \param glr glray *res_pt The glray-variable where the data is stored in.
 * \param n_glr const uint The number of the specific variable in glr.
 * \return void
 *
 * Pointers are used to facilitate the addressing of the counters.
 * Additional scoping is used for security and a better overview.
 */
void make_trace(const ray *const res_pt r, glray *res_pt glr, const uint n_glr)
{
    {
        /**< First, make the trace. */
        static const uint noi = 64;
        uint *i = &((*glr).glrs[n_glr].n_trace);
        if(!(*i))
        {
            assert(!(*glr).glrs[n_glr].trace_len);
            (*glr).glrs[n_glr].trace = alloc_point3(noi);
            (*glr).glrs[n_glr].trace_len = noi;
        }
        else if(!(*i % noi))
        {
            (*glr).glrs[n_glr].trace = realloc_point3((*glr).glrs[n_glr].trace, (*glr).glrs[n_glr].trace_len, noi);
            (*glr).glrs[n_glr].trace_len += noi;
        }
        cp3((*glr).glrs[n_glr].trace[*i].x, (*r).v.o);
        (*i)++;
        assert(*i == (*glr).glrs[n_glr].n_trace);
    }
    {
        /**< Then, care for proper indexing: */
        static const uint noj = 32;
        uint *j = &((*glr).glrs[n_glr].n_child);
        if(!(*j) && !(*glr).glrs[n_glr].child_len)
        {
            (*glr).glrs[n_glr].child = alloc_uint(noj);
            (*glr).glrs[n_glr].child_len = noj;
        }
        else if(!(*j % noj))
        {
            (*glr).glrs[n_glr].child = realloc_uint((*glr).glrs[n_glr].child, (*glr).glrs[n_glr].child_len, noj);
            (*glr).glrs[n_glr].child_len += noj;
        }
        if(*j != (*r).trans_child)
            (*j)++; /**< Add a new child. */
        (*glr).glrs[n_glr].child[*j]++; /**< Count the rays of this child. */
    }
}

/** \brief Reads only one char and discards the rest until a newline or EOF occurred.
 *
 * \param void
 * \return int Returns the char like getchar in terms of an int.
 *
 */
int getchar_dscrd_rmng(void)
{
    const int first = getchar(); /**< Very pedantic. */
    int c = first;
    while(c != '\n' && c != EOF)
        c = getchar(); /**< Keep getting more. */
    return first;
}

/** \brief Pause the program for some seconds.
 *
 * \param sec uint The pause time in seconds.
 * \return void
 *
 */
void wait(uint sec)
{
    clock_t endwait = clock() + sec * CLOCKS_PER_SEC;
    while(clock() < endwait) {}
}

/** \brief Compute the sine and cosine of an angle alpha in one step.
 *
 * \param alpha const double The angle in radian.
 * \param cosd double *res_pt The cosine of alpha.
 * \param sind double *res_pt The sine of alpha.
 * \return void
 *
 * If available, uses assembler code to compute the sine and cosine in one step.
 */
void sincosd(const double alpha, double *res_pt cosd, double *res_pt sind)
{
#if defined(__i386__)&&!defined(NO_ASM)
#if defined __GNUC__
#define ASM_SINCOS
    __asm__ __volatile__("fsincos" : "=t"(*cosd), "=u"(*sind) : "0"(alpha));
#elif defined _MSC_VER
#define ASM_SINCOS
#warning "assembler functionality not yet tested"
    __asm fld alpha
    __asm fsincos
    __asm fstp *cosd
    __asm fstp *sind
#endif
#endif
#ifndef ASM_SINCOS /**< Fall-back */
#warning "no assembler sincos function available"
    sincos_sqrt(alpha, cosd, sind);
#endif
}

/** \brief Fall-back function for sincosd.
 *
 * \param alpha double The angle in radians.
 * \param cosd double *res_pt The cosine of alpha.
 * \param sind double *res_pt The sine of alpha.
 * \return void
 *
 * Safes about 25 percent compared to a standard call of cos and sin.
 */
void sincos_sqrt(double alpha, double *res_pt cosd, double *res_pt sind)
{
    *cosd = cos(alpha);
    double t1 = sqrt(1. - (*cosd) * (*cosd));
    alpha = fabs(alpha);
    if(alpha < M_PI)
        *sind = t1;
    else if(alpha < M_PI2)
        *sind = -t1;
    else
    {
        alpha = fmod(alpha, M_PI2);
        if(alpha < M_PI)
            *sind = t1;
        else
            *sind = -t1;
    }
}

/** \brief Initialises a plane from three distinct points in R^3.
 *
 * \param p1 const double *res_pt The pointer to the first point3 variable.
 * \param p2 const double *res_pt The pointer to the second point3 variable.
 * \param p3 const double *res_pt The pointer to the third point3 variable.
 * \param r plane3 *res_pt The pointer to the output variable.
 * \return void
 *
 * The arrays are not allowed to overlap.
 * The points have to given in Cartesian coordinates.
 */
void init_plane3_pform(const double *res_pt p1, const double *res_pt p2, const double *res_pt p3, plane3 *res_pt r)
{
    double t1[3],
           t2[3];
    sub3(p2, p1, t1);
    sub3(p3, p1, t2);
    normedcross3(t1, t2, (*r).n);
    cp3((*r).o, p1);
}

/** \brief Initialises a plane from a given normal vector and a point within the plane in R^3.
 *
 * \param n const point3 *res_pt The pointer to the normal vector.
 * \param o const point3 *res_pt The pointer to the point within the plane.
 * \param r plane3 *res_pt The pointer to the output variable.
 * \return void
 *
 * The arrays are not allowed to overlap.
 */
void init_plane3(const point3 *res_pt n, const point3 *res_pt o, plane3 *res_pt r)
{
    double t = len_squ3((*n).x);
    cp3((*r).o, (*o).x);
    cp3((*r).n, (*n).x);
    if(fabs(t - 1.) > DBL_EPSILON)
        nvec3_ip((*r).n, sqrt(t));
}
