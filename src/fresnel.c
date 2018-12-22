#include "artpic.h"
#include "auxf.h"
#include "fresnel.h"

/**
 * A brief explanation of the following formulas can be found at http://de.wikipedia.org/wiki/Fresnelsche_Formeln
 * The 'incident' angle (cangl_i) has to be positive up to the definition
 */

/** \brief Calculates the Fresnel coefficient of the transmitted component orthogonal to the plane of incidence (a.k.a. TE)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return cdoub The output
 *
 * The cosine of the angles have to be positive
 */
cdoub cfres_trans_opol(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    assert(cangl_i > 0. && cangl_i <= 1. + DBL_EPSILON);
    cdoub t1 = (nt * cangl_o) / (ni * cangl_i), ret;
    if(mui != mut)
        t1 *= mui / mut;
    ret = 2. / (1. + t1);
    return ret;
}

/** \brief Calculates the Fresnel coefficient of the transmitted component parallel to the plane of incidence (a.k.a. TM)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return cdoub The output
 *
 * The cosine of the angles have to be positive
 */
cdoub cfres_trans_ppol(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    assert(cangl_i > 0. && cangl_i <= 1. + DBL_EPSILON);
    cdoub t1 = nt / ni, t2 = cangl_o / cangl_i, ret;
    if(mui != mut)
        t1 *= mui / mut;
    ret = 2. / (t1 + t2);
    return ret;
}

/** \brief Calculates the relative intensity of the transmitted ray orthogonal to the plane of incidence (a.k.a. TE)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return double The output
 *
 */
double fres_trans_opol_int(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    cdoub t1 = cfres_trans_opol(ni, nt, mui, mut, cangl_i, cangl_o);
    double t2;
    t1 *= conj(t1);
    t2 = cabs(nt * cangl_o / (ni * cangl_i)) * creal(t1);
    return t2;
}

/** \brief Calculates the relative intensity of the transmitted ray parallel to the plane of incidence (a.k.a. TM)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return double The output
 *
 */
double fres_trans_ppol_int(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    cdoub t1 = cfres_trans_ppol(ni, nt, mui, mut, cangl_i, cangl_o);
    double t2;
    t1 *= conj(t1);
    t2 = cabs(nt * cangl_o / (ni * cangl_i)) * creal(t1);
    return t2;
}

/** \brief Calculates the Fresnel coefficient of the reflected component orthogonal to the plane of incidence (a.k.a. TE)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return cdoub The output
 *
 * The cosine of the angles have to be positive
 */
cdoub cfres_refl_opol(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    assert(cangl_i > 0. && cangl_i <= 1. + DBL_EPSILON);
    cdoub t1 = (nt * cangl_o) / (ni * cangl_i), ret;
    if(mui != mut)
        t1 *= mui / mut;
    ret = (1. - t1) / (1. + t1);
    return ret;
}

/** \brief Calculates the Fresnel coefficient of the reflected component parallel to the plane of incidence (a.k.a. TM)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return cdoub The output
 *
 * The cosine of the angles have to be positive
 */
cdoub cfres_refl_ppol(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    assert(cangl_i > 0. && cangl_i <= 1. + DBL_EPSILON);
    cdoub t1 = nt / ni, t2 = cangl_o / cangl_i, ret;
    if(mui != mut)
        t1 *= mui / mut;
    ret = (t1 - t2) / (t1 + t2);
    return ret;
}

/** \brief Calculates the relative intensity of the reflected ray orthogonal to the plane of incidence (a.k.a. TE)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return double The output
 *
 */
double fres_refl_opol_int(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    cdoub t1 = cfres_refl_opol(ni, nt, mui, mut, cangl_i, cangl_o);
    t1 *= conj(t1);
    return creal(t1);
}

/** \brief Calculates the relative intensity of the reflected ray parallel to the plane of incidence (a.k.a. TM)
 *
 * \param ni const cdoub The refractive index of the medium the incident ray is in
 * \param nt const cdoub The refractive index of the medium the transmitted ray is in
 * \param mui const double The permeability of the medium the incident ray is in
 * \param mut const double The permeability of the medium the transmitted ray is in
 * \param cangl_i const double The cosine of the angle that the incident ray makes with the surface normal
 * \param cangl_o const cdoub The cosine of the angle that the outgoing ray makes with the surface normal
 * \return double The output
 *
 */
double fres_refl_ppol_int(const cdoub ni, const cdoub nt, const double mui, const double mut, const double cangl_i, const cdoub cangl_o)
{
    cdoub t1 = cfres_refl_ppol(ni, nt, mui, mut, cangl_i, cangl_o);
    t1 *= conj(t1);
    return creal(t1);
}
