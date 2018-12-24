#include "artpic.h"
#include "lina.h"
#include "ray.h"
#include "refr_data.h"
#include "msg.h"
#include "auxf.h"
#include "tests.h"
#include "rot.h"

gen_ray_info GLOBAL_RAY_INFO;

/** \brief Compares two media.
 *
 * \param n1 const cdoub A complex refractive index.
 * \param mu1 const double A permeability.
 * \param n2 const cdoub A complex refractive index.
 * \param mu2 const double A permeability.
 * \return uchar Returns 1 if the media are equal, otherwise 0.
 *
 */
uchar compare_media(const cdoub n1, const double mu1, const cdoub n2, const double mu2)
{
    return ((mu1 == mu2 && n1 == n2) ? 1 : 0);
}

/** \brief Generates a 'cylindrical' bundle of rays.
 *
 * \param rs ray *res_pt The output.
 * \param lam const double The wavelength of the rays.
 * \param no const uint The number of rays, i.e. the length of the array rs.
 * \param res_azim const uint The azimuthal resolution of the bundle, i.e. the number of rays at a constant radius.
 * \param res_rad const uint The radial resolution of the bundle, i.e. the number of rays at a constant (azimuthal) angle.
 * \param max_rad const double The maximum radial circumference.
 * \param i const point3 *res_pt The starting point of the bundle.
 * \param f const point3 *res_pt The end point of the bundle.
 * \param pol const point3 *res_pt The polarisation of the bundle.
 * \return void
 *
 * Generates a bundle of rays which point from i to f, i.e. they start from one plane (with point i) and end up on in parallel one (with point f).
 * The polarisation is homogeneous and linear and has to be entered in the xy plane.
 * The rays generated rays are coherent.
 */
void gen_startrays_straight(ray *res_pt rs, const double lam,
                            const uint no,
                            const uint res_azim, const uint res_rad, const double max_rad,
                            const point3 *res_pt i, const point3 *res_pt f, const point3 *res_pt pol)
{
    uint l, j, k;
    const double z_axis[3] = {0., 0., 1.}, /**< This shall not be modified, it is the initial direction of the bundle. */
                             ra = M_PI2 / res_azim,
                             rr = max_rad / res_rad;
    double axy, ayz, axz,
           alpha, beta,
           t1[3],
           t2,
           cpol[3], normalpol[3];
    uchar transform_pol = 1;
    sub3((*f).x, (*i).x, t1);
    envec3_ip(t1); /**< t1 is a normed vector pointing from i to f. */
    cp3(cpol, (*pol).x);
    envec3_ip(cpol);
    normedcross3(z_axis, cpol, normalpol);
    snprintf(GLOBAL_RAY_INFO.info, L_INFOLENGTH,
             "polarisation at origin is (x %g, y %g) and declared as 'o' polarized."
             " rays starting in (%g, %g, %g) direction.",
             cpol[0], cpol[1],
             z_axis[0], z_axis[1], z_axis[2]);
    for(l = 0, j = 1; j <= res_rad; j++)
        for(k = 0, t2 = rr * j; k < res_azim; k++, l++)
        {
            rs[l] = (ray){.lam = lam,
                          .v = (line3){.o = {t2, M_PIh, ra * k},
                                       .l = 0.,
                                       .r = {z_axis[0], z_axis[1], z_axis[2]}},
                          .oamp = 1.,
                          .pamp = 0.,
                          .ophase = 0.,
                          .pphase = 0.,
                          .oint = 1.,
                          .pint = 0.,
                          .mu_i = MU_VAC,
                          .n_i = N_VAC,
                          .opol = {cpol[0], cpol[1], cpol[2]},
                          .ppol = {normalpol[0], normalpol[1], normalpol[2]}};
            s_to_c3_ip(rs[l].v.o);
        }
    assert(l == no);
    GLOBAL_RAY_INFO.count_gen = l;
    GLOBAL_RAY_INFO.count_hit =
    GLOBAL_RAY_INFO.count_lost =
    GLOBAL_RAY_INFO.count_exhstd = 0;
    axy = angl3_xy(t1);
    ayz = angl3_yz(t1);
    axz = angl3_xz(t1);
    alpha = atan2(t1[0], -t1[1]);
    beta = acos(t1[2]);
    if(axy == 0xDEAD && ayz == 0xDEAD)
        transform_pol = 0; /**< t1 is at the origin. */
    else if(axy == 0xDEAD && ayz == M_PIh) /**< t1 is parallel to the z axis. */
    {
        transform_pol = 0;
        for(l = 0; l < no; l++)
            add3(rs[l].v.o, (*i).x, rs[l].v.o);
    }
    else if(axy == 0xDEAD && ayz == M_PI3h) /**< t1 is anti-parallel ... */
    {
        transform_pol = 0;
        for(l = 0; l < no; l++)
        {
            rev3_ip(rs[l].v.r);
            rev3_ip(rs[l].ppol);
            add3(rs[l].v.o, (*i).x, rs[l].v.o);
        }
    }
    else if(axy == 0. && ayz == 0xDEAD) /**< t1 is parallel to the x axis. */
        for(l = 0; l < no; l++)
        {
            rot3_y(rs[l].v.r, M_PIh);
            rottrans3_y(rs[l].v.o, M_PIh, (*i).x);
        }
    else if(axy == M_PI && ayz == 0xDEAD) /**< t1 is anti-parallel ... */
        for(l = 0; l < no; l++)
        {
            rot3_y(rs[l].v.r, M_PI3h);
            rottrans3_y(rs[l].v.o, M_PI3h, (*i).x);
        }
    else if(axz == 0xDEAD && ayz == 0.) /**< t1 is parallel to the y axis. */
        for(l = 0; l < no; l++)
        {
            rot3_x(rs[l].v.r, M_PI3h);
            rottrans3_x(rs[l].v.o, M_PI3h, (*i).x);
        }
    else if(axz == 0xDEAD && ayz == M_PI) /**< t1 is anti-parallel ... */
        for(l = 0; l < no; l++)
        {
            rot3_x(rs[l].v.r, M_PIh);
            rottrans3_x(rs[l].v.o, M_PIh, (*i).x);
        }
    else
        for(l = 0; l < no; l++)
        {
            rot3_xpz(rs[l].v.r, beta, alpha);
            rottrans3_xpz(rs[l].v.o, beta, alpha, (*i).x);
        }
    if(transform_pol)
        for(l = 0; l < no; l++)
        {
            rot3_xpz(rs[l].opol, beta, alpha);
            normedcross3(rs[l].v.r, rs[l].opol, rs[l].ppol);
        }
    for(l = 0; l < no; l++)
    {
        snprintf(rs[l].info, INFOLENGTH, "no. %u from (x %g,y %g,z %g)", l, rs[l].v.o[0], rs[l].v.o[1], rs[l].v.o[2]);
        assert_ray(rs, "init", ERR_ARG);
    }
}

/** \brief Computes the maximum amplitude of a ray independent on its phase.
 *
 * \param r const ray* A pointer to a ray.
 * \return double The output.
 *
 * This function does not consider a non-linear state of polarisation.
 */
double get_ray_max_amp(const ray *r)
{
    return sqrt(POW2((*r).pamp) + POW2((*r).oamp));
}

/** \brief Computes the amplitude of a ray which is dependent on the phase.
 *
 * \param r const ray* A pointer to a ray.
 * \return double The output.
 *
 */
double get_ray_phased_amp(const ray *r)
{
    return sqrt(POW2((*r).pamp * cos((*r).pphase)) + POW2((*r).oamp * cos((*r).ophase))); /**< Check if this is coherent with the use of the other phase. Maybe use a sin. */
}

/** \brief Computes the sum of the o- and p-intensity of a ray.
 *
 * \param r const ray* A pointer to a ray.
 * \return double The output.
 *
 */
double get_ray_int(const ray *r)
{
    return (*r).oint + (*r).pint;
}

/** \brief Propagates a ray according to the latest event.
 *
 * \param r ray* A pointer to a ray.
 * \return void
 *
 * The propagation length is stored in the scalar multiplier of the ray's line3 variable.
 * The phase will be mapped to [0, 2 pi].
 */
void propagate_ray(ray *r)
{
    double t1[3],
           t2 = (*r).v.l;
    addnmul3((*r).v.o, t2, (*r).v.r, t1); /**< Get the end point. */
    if(fabs(dist3((*r).v.o, t1) - t2) > DBL_EPSILON6)
        printf("The directional vector does not seem to be correctly normed: %g\n",
               fabs(dist3((*r).v.o, t1) - t2));
    assert(fabs(dist3((*r).v.o, t1) - t2) <= DBL_EPSILON6); /**< Checks if the directional vector was correctly normed. */
    cp3((*r).v.o, t1);
    (*r).travel += t2;
    (*r).v.l = 0.;
    const double t3 = M_PI2 * t2 / (*r).lam;
    t2 = t3 * creal((*r).n_i);
    (*r).ophase = fmod((*r).ophase + t2, M_PI2);
    (*r).pphase = fmod((*r).pphase + t2, M_PI2);
    t2 = exp(-t3 * cimag((*r).n_i));
    (*r).oamp *= t2;
    (*r).pamp *= t2;
}

/** \brief Propagates a ray a tiny distance so that it gets away from the latest event.
 *
 * \param r ray* A pointer to a ray.
 * \return void
 *
 */
void propagate_ray_eps(ray *r)
{
    double t1[3],
           t2 = EPS_RAY_PROP_LEN;
    assert(fabs(len_squ3((*r).v.r) - 1.) <= DBL_EPSILON6);
    addnmul3((*r).v.o, t2, (*r).v.r, t1);
    cp3((*r).v.o, t1);
    (*r).travel += t2;
    assert((*r).v.l == 0.);
    const double t3 = M_PI2 * t2 / (*r).lam;
    t2 = t3 * creal((*r).n_i);
    (*r).ophase = fmod((*r).ophase + t2, M_PI2);
    (*r).pphase = fmod((*r).pphase + t2, M_PI2);
    t2 = exp(-t3 * cimag((*r).n_i));
    (*r).oamp *= t2;
    (*r).pamp *= t2;
}

/** \brief Sets the correct refractive index and the permeability of the intersection and checks which index is behind the surface the ray hits.
 *
 * \param ri const ray *res_pt The ingoing ray.
 * \param i intrsec *res_pt The intersection.
 * \param n_prtcl const cdoub The refractive index of the particle the ray hits.
 * \param mu_prtcl const double The permeability of the particle the ray hits.
 * \return uchar tells whether The media are different or not.
 *
 * If ri.v.r and i.normal point in the same half sphere, the ray is travelling into the object.
 * If they point in different half spheres, the ray is travelling out of it.
 * If the index of the travelling ray matches the upcoming index, return 0.
 */
uchar set_refr_index(const ray *res_pt ri, intrsec *res_pt i, const cdoub n_prtcl, const double mu_prtcl)
{
    assert((*i).incdnc == OBLIQUE || (*i).incdnc == VERTICAL);
    if((*i).cangl > 0.)
    {
        assert(dot3((*ri).v.r, (*i).normal) > 0.);
        (*i).n_f = n_prtcl;
        (*i).mu_f = mu_prtcl;
    }
    else
    {
        assert(dot3((*ri).v.r, (*i).normal) < 0.);
        (*i).n_f = N_VAC;
        (*i).mu_f = MU_VAC;
    }
    return !compare_media((*ri).n_i, (*ri).mu_i, (*i).n_f, (*i).mu_f);
}

/** \brief Registers an event into a hit_screen variable.
 *
 * \param r const ray *res_pt A ray that ended somewhere.
 * \param hs hit_screen *res_pt The pointer to the input variable.
 * \param cangl_i const double The cosine of the angle between the ray and the screen's normal.
 * \param state const ray_end_state A variable describing the end state of the ray.
 * \return void
 *
 */
void reg_hit(const ray *res_pt r, hit_screen *res_pt hs, const double cangl_i, const ray_end_state state)
{
    *hs = (hit_screen){
        .p = {(*r).v.o[0], (*r).v.o[1], (*r).v.o[2]},
        .lam = (*r).lam,
        .cos_incdnc = cangl_i,
        .opol = {(*r).opol[0], (*r).opol[1], (*r).opol[2]},
        .ppol = {(*r).ppol[0], (*r).ppol[1], (*r).ppol[2]},
        .coamp = (*r).oamp * (cos((*r).ophase) + sin((*r).ophase) * 1.i),
        .cpamp = (*r).pamp * (cos((*r).pphase) + sin((*r).pphase) * 1.i),
        .oint = (*r).oint,
        .pint = (*r).pint,
        .tir = (*r).tir,
        .state = state};
    strncpy((*hs).info, (*r).info, INFOLENGTH);
}
