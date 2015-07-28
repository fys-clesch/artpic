#include "artpic.h"
#include "lina.h"
#include "msg.h"
#include "auxf.h"
#include "rot.h"
#include "refr_data.h"
#include "fresnel.h"
#include "ray_aux.h"
#include "tests.h"
#include "ray.h"

/**
 * definition of the surface normal according to (amongst others) E. Hecht's book 'Optics' on page 103 (4th ed.) figure 4.24
 */

/** \brief computes the cosine of the angle between a normal vector and the ray's vector
 *
 * \param r const ray *res_pt a pointer to a ray
 * \param nv const double *res_pt a pointer to a normal vector
 * \return double the output
 *
 * both vectors have to be normalized which is not checked
 */
double cangl_in(const ray *res_pt r, const double *res_pt nv)
{
	return (*r).v.r[0] * nv[0] + (*r).v.r[1] * nv[1] + (*r).v.r[2] * nv[2];
}

/** \brief computes the complex cosine of the angle of the transmitted beam
 *
 * \param ni const cdoub the refractive index of the medium the incident ray is in
 * \param nt const cdoub the refractive index of the medium the transmitted ray is in
 * \param cangl_i const double the cosine of the angle that the incident ray makes with the surface normal
 * \param tir uchar* a return value to check if tir occured or not
 * \return cdoub the cosine of a probably complex argument
 *
 * if the real part of t1 is less or equal than zero, tir occured
 * according to IEEE Std 1003.1, 2004 Edition the csqurt-function behaves as follow:
 * These functions shall return the complex square root value, in the range of the right half-plane (including the imaginary axis).
 */
cdoub ccangl_out(const cdoub ni, const cdoub nt, const double cangl_i, uchar *tir)
{
	const double t1 = 1. - cangl_i * cangl_i;
	const cdoub t2 = 1. - ni * ni / (nt * nt) * t1;
	const double t3 = creal(t2),
	             sign = cangl_i < 0. ? -1. : 1.; /**< this is used to set the proper orientation of the angle */
	if(t3 > 0.)
	{
		*tir = 0;
		return sign * csqrt(t2);
	}
	else if(t3 == 0.)
	{
		*tir = 1;
		return sign * csqrt(t2);
	}
	else
	{
		*tir = 1;
		return sign * (-csqrt(t2)); /**< see the csqrt definition */
	}
}

/** \brief computes the directional vector of the reflected ray
 *
 * \param in const double *res_pt the directional vector of the incident ray
 * \param ref double *res_pt the directional vector of the reflected ray
 * \param nv const double *res_pt the surface normal
 * \param cangl_i const double the cosine of the angle between the surface normal and the incident ray
 * \return void
 *
 */
void ray_refl(const double *res_pt in, double *res_pt ref, const double *res_pt nv, const double cangl_i)
{
	double t1 = -2.*cangl_i;
	addnmul3(in, t1, nv, ref);
	if(fabs((t1 = len_squ3(ref)) - 1.) > DBL_EPSILON)
		nvec3_ip(ref, sqrt(t1)); /**< normalizes the vector */
}

/** \brief computes the transmitted ray
 *
 * \param in const double *res_pt the incident ray
 * \param trans double *res_pt the output as the transmitted ray
 * \param nv const double *res_pt the surface normal
 * \param ni const cdoub the refractive index of the medium the incident ray is in
 * \param nt const cdoub the refractive index of the medium the transmitted ray is in
 * \param cangl_i const double the cosine of the angle of the incident ray with the surface normal
 * \param cangl_o const cdoub the complex cosine of the angle of the transmitted ray with the surface normal
 * \return void
 *
 * nv points inwards object
 * multiply t2 by -1. if the definition of nv is different from E. Hecht's
 */
void ray_trans(const double *res_pt in, double *res_pt trans, const double *res_pt nv, const cdoub ni, const cdoub nt, const double cangl_i, const cdoub cangl_o)
{
	const cdoub tn = ni / nt,
	            t1 = cangl_o - tn * cangl_i;
	double t2;
	lincomb3(creal(tn), in, creal(t1), nv, trans);
	if(fabs((t2 = len_squ3(trans)) - 1.) > DBL_EPSILON)
		nvec3_ip(trans, sqrt(t2));
#ifndef NDEBUG
	if(fabs(dot3(nv, trans) - creal(cangl_o)) > DBL_EPSILON4)
	{
		fprintf(stderr, "fabs(dot3(nv,trans)-creal(cangl_o)): %g\n"
		        "incoming vector: (%15g, %15g, %15g)\n"
		        "normal         : (%15g, %15g, %15g)\n"
		        "outgoing       : (%15g, %15g, %15g)\n"
		        "dot3(incoming, normal)    : %15g\n"
		        "ingoing cos(angle)        : %15g\n"
		        "dot3(outgoing, normal)    : %15g\n"
		        "outgoing creal(cos(angle)): %15g\n"
		        "%%%% - %%                    : %15g\n"
		        "outgoing cabs(cos(angle)) : %15g\n"
		        "%%%%%%%% - %%                  : %15g\n"
		        "outgoing cimag(cos(angle)): %15g",
		        fabs(dot3(nv, trans) - creal(cangl_o)), in[0], in[1], in[2], nv[0], nv[1], nv[2], trans[0], trans[1], trans[2],
		        dot3(in, nv), cangl_i, dot3(trans, nv), creal(cangl_o),
		        dot3(trans, nv) - creal(cangl_o), cabs(cangl_o), dot3(trans, nv) - cabs(cangl_o), cimag(cangl_o));
		error_msg("this is probably due to a not detected tir event or due to an intersection of particles with different refractive indices.", ERR_ARG);
	}
#endif
}

/** \brief handles the reflection of a ray at an intersection
 *
 * \param ri const ray *res_pt a pointer to a ray
 * \param rrefl ray *res_pt the output
 * \param isec const intrsec *res_pt a pointer to an intersection
 * \param ccangl_f const cdoub the complex cosine of the angle between the surface normal and the transmitted ray
 * \return uchar if the fresnel coefficients for both polarisations are smaller than MINI_REFLECTION_COEF, return 0
 *
 */
uchar handle_refl(const ray *res_pt ri, ray *res_pt rrefl, const intrsec *res_pt isec, const cdoub ccangl_f)
{
	double t1 = fabs((*isec).cangl), t2;
	const cdoub cabs_cangl_f = cabs_real(ccangl_f),
	            cro = cfres_refl_opol((*ri).n_i, (*isec).n_f, (*ri).mu_i, (*isec).mu_f, t1, cabs_cangl_f), /**< TE component */
	            crp = cfres_refl_ppol((*ri).n_i, (*isec).n_f, (*ri).mu_i, (*isec).mu_f, t1, cabs_cangl_f); /**< TM component */
	t1 = cabs(cro);
	if(t1 < MINI_REFLECTION_COEF) t1 = 0.;
	t2 = cabs(crp);
	if(t2 < MINI_REFLECTION_COEF) t2 = 0.;
	if(t1 > (1. + DBL_EPSILON) && t2 > (1. + DBL_EPSILON))
	{
		fprintf(stderr, "\novershoot of the reflection coefficient\n"
		        "- orthogonal (TE): %g\n- parallel (TM)  : %g", t1 - 1., t2 - 1.);
		error_msg("this is probably due to the use of complex refractive indices and a tir event", ERR_ARG);
	}
	if(t1 == 0. && t2 == 0.) return 0;
	else
	{
		memcpy(rrefl, ri, sizeof(ray));
		if((*rrefl).oamp > 0.)(*rrefl).oamp *= t1;
		(*rrefl).ophase += atan2_up(cimag(cro), creal(cro));
		if((*rrefl).pamp > 0.)(*rrefl).pamp *= t2;
		(*rrefl).pphase += atan2_up(cimag(crp), creal(crp));
		(*rrefl).oint *= (t1 * t1);
		(*rrefl).pint *= (t2 * t2);
		if((*isec).incdnc == OBLIQUE) ray_refl((*ri).v.r, (*rrefl).v.r, (*isec).normal, (*isec).cangl);
		else rev3((*ri).v.r, (*rrefl).v.r);
		normedcross3((*rrefl).v.r, (*rrefl).opol, (*rrefl).ppol);
		return 1;
	}
}

/** \brief handles the transmission of a ray at an intersection
 *
 * \param ri const ray *res_pt a pointer to a ray
 * \param rtrans ray *res_pt the output
 * \param isec const intrsec *res_pt a pointer to an intersection
 * \param ccangl_f const cdoub the complex cosine of the angle between the surface normal and the transmitted ray
 * \return uchar refer to the description:
 *
 * if the amplitude of the fresnel coefficients are smaller than MINI_TRANSMISSION_COEF, a fp-error may have occured
 * in this case, return 0
 * if a polarization of the ray is not transmitted, ('half-tir') return 1
 * if both polarizations are transmitted, return 2
 */
uchar handle_trans(const ray *res_pt ri, ray *res_pt rtrans, const intrsec *res_pt isec, const cdoub ccangl_f)
{
	const double fabs_cangl_i = fabs((*isec).cangl);
	double t1, t2;
	const cdoub cabs_cangl_f = cabs_real(ccangl_f),
	            cto = cfres_trans_opol((*ri).n_i, (*isec).n_f, (*ri).mu_i, (*isec).mu_f, fabs_cangl_i, cabs_cangl_f), /**< TE component */
	            ctp = cfres_trans_ppol((*ri).n_i, (*isec).n_f, (*ri).mu_i, (*isec).mu_f, fabs_cangl_i, cabs_cangl_f); /**< TM component */
	t1 = cabs(cto);
	if(t1 < MINI_TRANSMISSION_COEF) t1 = 0.;
	t2 = cabs(ctp);
	if(t2 < MINI_TRANSMISSION_COEF) t2 = 0.;
	if(t1 > (2. + DBL_EPSILON))
	{
		fprintf(stderr, "\novershoot of the transmission coefficient\n"
		        "- orthogonal (TE): %g\n"
		        "- parallel (TM)  : %g\n"
		        "- orthogonal (TE) - 2.: %g",
		        t1, t2, t1 - 2.);
		error_msg("this is probably due to the use of complex refractive indices and a tir event", ERR_ARG);
	}
	if(t1 == 0. && t2 == 0.)
	{
		error_msg("unexpected result. this event should have been handled by the detection of tir before.", ERR_ARG);
		return 0;
	}
	else
	{
		memcpy(rtrans, ri, sizeof(ray));
		if((*rtrans).oamp > 0.)(*rtrans).oamp *= t1;
		(*rtrans).ophase += atan2_up(cimag(cto), creal(cto));
		if((*rtrans).pamp > 0.)(*rtrans).pamp *= t2;
		(*rtrans).pphase += atan2_up(cimag(ctp), creal(ctp));
		const double t3 = cabs((*isec).n_f * cabs_cangl_f / ((*ri).n_i * fabs_cangl_i));
		(*rtrans).oint *= (t1 * t1 * t3);
		(*rtrans).pint *= (t2 * t2 * t3);
		if((*isec).incdnc == OBLIQUE)
			ray_trans((*ri).v.r, (*rtrans).v.r, (*isec).normal, (*ri).n_i, (*isec).n_f, (*isec).cangl, ccangl_f);
		else
			cp3((*rtrans).v.r, (*ri).v.r);
		assert(dot3((*rtrans).v.r, (*ri).v.r) > 0.); /**< the transmitted ray should point in the same half-sphere as the primary one */
		(*rtrans).n_i = (*isec).n_f;
		(*rtrans).mu_i = (*isec).mu_f;
		normedcross3((*rtrans).v.r, (*rtrans).opol, (*rtrans).ppol);
		if(t2 == 0. || t1 == 0.) return 1; /**< 'half-tir' */
		else return 2;
	}
}

/** \brief handles both the transmission and reflection of a ray at an intersection
 *
 * \param ri ray *res_pt the incoming ray. it will be used as the reflected ray.
 * \param rsec ray *res_pt the secondary ray. it will be used as the transmitted ray in case of no tir.
 * \param isec const intrsec *res_pt a pointer to an intersection
 * \return uchar returns 1 if the secondary ray is created, 0 otherwise
 *
 */
uchar handle_reflntrans(ray *res_pt ri, ray *res_pt rsec, const intrsec *res_pt isec)
{
	assert((*ri).n_i != (*isec).n_f);
	assert_ray(ri, "input ray", ERR_ARG);
	uchar tir;
	const cdoub cangl_f = ccangl_out((*ri).n_i, (*isec).n_f, (*isec).cangl, &tir);
	ray rswap;
	if(!tir)
	{
		uchar ctrans = handle_trans(ri, rsec, isec, cangl_f);
		if(ctrans != 2)
			(*ri).tir++; /**< one polarisation was not transmitted ('half-tir') */
		uchar crefl = handle_refl(ri, &rswap, isec, cangl_f);
		if(crefl) memcpy(ri, &rswap, sizeof(ray));
		if(!ctrans && !crefl) /**< no transmission and no reflection: */
		{
			error_msg("fp error in handle_trans, handle_refl or ccangl_out", ERR_ARG);
			(*ri).oamp = (*ri).pamp = 0.;
			goto sec_dead;
		}
		else if(!ctrans && crefl) /**< no transmission: */
		{
			error_msg("unexpected case in handle_reflntrans, most likely due to a fp error", ERR_ARG);
			goto sec_dead;
		}
		else if(ctrans && !crefl) /**< no reflection: */
		{
			memcpy(ri, rsec, sizeof(ray)); /**< overwrite the ingoing ray with the secondary */
			goto sec_dead;
		}
		else /**< both rays are active. get away from the interface: */
		{
			propagate_ray_eps(ri);
			propagate_ray_eps(rsec);
			(*rsec).trans_child++;
			uchar tc = assert_ray(ri, "reflected", ERR_ARG);
			tc &= assert_ray(rsec, "secondary", ERR_ARG);
			if(!tc) print_intrsec(isec, "standard");
			return 1;
		}
sec_dead:
		propagate_ray_eps(ri);
		(*rsec).lam = 0xDEAD; /**< the ray is not occupied yet and can be used for another intersection */
		assert_ray(ri, "transmitted, no reflection", ERR_ARG);
		return 0;
	}
	else /**< tir occured: */
	{
		printf("TIRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR\n");
		getchar();
		if(!handle_refl(ri, &rswap, isec, cangl_f))
			error_msg("tir occured, but reflection coefficients are near 0. "
			          "check ccangl_out in handle_refl (or before) for rounding errors.", ERR_ARG);
		memcpy(ri, &rswap, sizeof(ray));
		propagate_ray_eps(ri);
		(*ri).tir++;
		assert_ray(ri, "reflected, tir", ERR_ARG);
		return 0;
	}
}

void map_pol_to_lc(ray *r, const double alpha)
{
	const double Ax = (*r).oamp,
	             ax2 = POW2(Ax),
	             Ay = (*r).pamp,
	             ay2 = POW2(Ay),
	             txp = (*r).ophase,
	             typ = (*r).pphase,
	             salpha = sin(alpha),
	             calpha = cos(alpha),
	             spx = sin(txp),
	             cpx = cos(txp),
	             spy = sin(typ),
	             cpy = cos(typ);
	(*r).oamp = sqrt(ax2 + ay2 + (Ax - Ay) * (Ax + Ay) * cos(2.*alpha) - 2.*Ax * Ay * cos(txp - typ) * sin(2.*alpha)) / M_SQRT2;
	(*r).pamp = sqrt(ax2 + ay2 + (-ax2 + ay2) * cos(2.*alpha) + 2.*Ax * Ay * cos(txp - typ) * sin(2.*alpha)) / M_SQRT2;
	(*r).ophase = atan2_up(Ax * calpha * spx - Ay * salpha * spy, Ax * calpha * cpx - Ay * cpy * salpha);
	(*r).pphase = atan2_up(Ax * salpha * spx + Ay * calpha * spy, Ay * calpha * cpy + Ax * cpx * salpha);
	double m[9];
	gen_mat_rot3_n(m, (*r).v.r, -alpha); /**< the rotation of the polarisation vector has to compensate for the rotation of the amplitude and phase */
	inner_product3_mat_vec_ip(m, (*r).opol);
	inner_product3_mat_vec_ip(m, (*r).ppol);
}

/** \brief rotates a ray's polarisation according to a new angle alpha to match the local coordinate system at an intersection
 *
 * \param r ray* a pointer to the ray variable
 * \param alpha const double the angle between a tangential vector of the intersection and the orthogonal polarisation
 * \return void
 *
 * this is a unitary transform and does not change the polarisation state of the ray
 * this function should not be called when the angle is 0 or +/- pi
 * the angle is counted positively in a counter-clockwise sense
 */
void map_pol_to_lc_opt(ray *r, const double alpha)
{
	const double ax = (*r).oamp,
	             ax2 = POW2(ax),
	             ay = (*r).pamp,
	             ay2 = POW2(ay),
	             salpha = sin(alpha),
	             calpha = cos(alpha),
	             s2alpha = sin(2.*alpha),
	             c2alpha = cos(2.*alpha),
	             cpdiff = cos((*r).ophase - (*r).pphase),
	             spx = sin((*r).ophase),
	             cpx = cos((*r).ophase),
	             spy = sin((*r).pphase),
	             cpy = cos((*r).pphase),
	             t1 = 2.*ax * ay * cpdiff * s2alpha,
	             t2 = ax2 + ay2,
	             t3 = ax * spx,
	             t4 = ay * spy,
	             t5 = ay * cpy,
	             t6 = ax * cpx;
	(*r).oamp = sqrt(t2 + (ax - ay) * (ax + ay) * c2alpha - t1) / M_SQRT2;
	(*r).pamp = sqrt(t2 + (ay2 - ax2) * c2alpha + t1) / M_SQRT2;
	(*r).ophase = atan2_up(t3 * calpha - t4 * salpha, t6 * calpha - t5 * salpha);
	(*r).pphase = atan2_up(t3 * salpha + t4 * calpha, t5 * calpha + t6 * salpha);
	double m[9];
	gen_mat_rot3_n(m, (*r).v.r, -alpha); /**< the rotation of the polarisation vector has to compensate for the rotation of the amplitude and phase */
	inner_product3_mat_vec_ip(m, (*r).opol);
	inner_product3_mat_vec_ip(m, (*r).ppol);
}

/** \brief maps an ingoing ray into the orthogonal and the parallel component with respect to the plane spanned by the ray vector and the surface normal
 *
 * \param ri ray *res_pt the ingoing ray
 * \param isec const intrsec *res_pt the intersection
 * \return uchar returns 0 if the components are already aligned with the local coordinates and 1 if not
 *
 */
uchar orthogo_ray_at_intersec(ray *res_pt ri, const intrsec *res_pt isec)
{
	assert(fabs((*isec).cangl) > DBL_EPSILON);
	assert((*isec).incdnc != VERTICAL);
	double tang[3], tswap[3];
#ifndef NDEBUG
	if(!assert_ray(ri, "before orthogonalization", ERR_ARG))
		print_intrsec(isec, "orthogo-start");
	const double c1 = get_ray_phased_amp(ri);
#endif
	normedcross3((*ri).v.r, (*isec).normal, tang); /**< tang is tangential to the surface */
	const double tcangl = save_dot3(tang, (*ri).opol); /**< the cosine of the angle between tang and the orthogonal polarisation */
	tswap[0] = fabs(tcangl);
	if(fabs(tswap[0] - 1.) <= DBL_EPSILON2) /**< the orthogonal polarisation is (anti-)parallel to tang */
		return 0;
	else if(tswap[0] <= DBL_EPSILON2) /**< the parallel polarisation is (anti-)parallel to tang */
	{
		cp3(tswap, (*ri).opol);
		rev3_ip(tswap); /**< reverse to keep to the correct orientation between p, o and the directional vector */
		(*ri).ophase += M_PI; /**< compensate the prior step with pi */
		/**< swap the polarisation vectors: */
		cp3((*ri).opol, (*ri).ppol);
		cp3((*ri).ppol, tswap);
		/**< swap the amplitudes: */
		tswap[0] = (*ri).pamp;
		(*ri).pamp = (*ri).oamp;
		(*ri).oamp = tswap[0];
		/**< swap the phases: */
		tswap[0] = (*ri).pphase;
		(*ri).pphase = (*ri).ophase;
		(*ri).ophase = tswap[0];
		uchar tc = assert_orthogo(tang, (*ri).ppol, "tang", "p", ERR_ARG);
		tc &= assert_orthogo(tang, (*ri).v.r, "tang", "r", ERR_ARG);
		tc &= assert_parallel(tang, (*ri).opol, "tang", "o", ERR_ARG);
		if(!tc) print_intrsec(isec, "orthogo-special");
	}
	else /**< the general case */
	{
		tswap[0] = dot3(tang, (*ri).ppol);
		if(tswap[0] > 0.)
			tswap[1] = -acos(tcangl);
		else
			tswap[1] = acos(tcangl);
		map_pol_to_lc(ri, tswap[1]);
		uchar tc = assert_orthogo(tang, (*ri).ppol, "tang", "p", ERR_ARG);
		tc &= assert_orthogo(tang, (*ri).v.r, "tang", "r", ERR_ARG);
		tc &= assert_direct_parallel(tang, (*ri).opol, "tang", "o", ERR_ARG);
		if(!tc)
		{
			fprintf(stdout, "*** angle between tang and p': %g\n"
			        "*** angle between tang and o': %g\n", tswap[0], tcangl);
			print_intrsec(isec, "orthogo-general");
		}
	}
#ifndef NDEBUG
	assert_ray(ri, "after orthogonalization", ERR_ARG);
	const double c2 = get_ray_phased_amp(ri);
	if(fabs(c1 - c2) >= DBL_EPSILON6)
	{
		error_msg("numerical deviation. press enter to check it.", ERR_ARG);
		fprintf(stderr, "%g - %g = %g\n", c1, c2, fabs(c1 - c2));
		if(STOP_AT_ERR) getchar();
	}
#endif
	return 1;
}
