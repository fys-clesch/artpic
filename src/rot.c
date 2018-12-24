#include "artpic.h"
#include "auxf.h"
#include "rot.h"

/** \brief Computes a matrix to rotate a point in R^3 around the origin.
 *
 * \param m double* The output.
 * \param ea const double* An array of three angles of rotation.
 * \return void
 *
 * ea contains three Euler angles, whereas the line of nodes is the line perpendicular to both z and z' axis.
 * The order is xyz, right handed, see http://en.wikipedia.org/wiki/Euler_angles.
 */
void gen_mat_rot3_euler_xyz(double *m, const double *ea)
{
    double c1, c2, c3, s1, s2, s3, t1, t2;
    sincosd(ea[0], &c1, &s1);
    sincosd(ea[1], &c2, &s2);
    sincosd(ea[2], &c3, &s3);
    t1 = c1 * s3;
    t2 = s2 * c3;
    m[0] = c2 * c3;
    m[1] = c2 * s3;
    m[2] = -s2;
    m[3] = s1 * t2 - t1;
    m[4] = c1 * c3 + s1 * s2 * s3;
    m[5] = c2 * s1;
    m[6] = s1 * s3 + c1 * t2;
    m[7] = -c3 * s1 + s2 * t1;
    m[8] = c1 * c2;
}

/** \brief Computes a matrix to rotate a point in R^3 around the origin.
 *
 * \param m double* The output.
 * \param ea const double* An array of two angles of rotation.
 * \return void
 *
 * ea contains two Euler angles, whereas the line of nodes is the line perpendicular to both z and z' axis.
 * The order is zx, right handed.
 */
void gen_mat_rot3_euler_zx(double *m, const double *ea)
{
    double c1, c2, s1, s2;
    sincosd(ea[0], &c1, &s1);
    sincosd(ea[1], &c2, &s2);
    m[0] = c1;
    m[1] = -c2 * s1;
    m[2] = s1 * s2;
    m[3] = s1;
    m[4] = c1 * c2;
    m[5] = -c1 * s2;
    m[6] = 0.;
    m[7] = s2;
    m[8] = c2;
}

/** \brief Computes a matrix to rotate a point in R^3 around the origin with Euler angles.
 *
 * \param m double* The output.
 * \param gamma const double The third Euler angle.
 * \param beta const double The second Euler angle.
 * \param alpha const double The first Euler angle.
 * \return void
 *
 * The angles are Euler angles
 * First rotates the point by alpha around z, then rotates it by beta around the x' and finally rotates it by gamma around z''.
 * Assumes a right handed system.
 */
void gen_mat_rot3_euler_zpxpz(double *m, const double gamma, const double beta, const double alpha)
{
    double sa, ca, sb, cb, sg, cg, t1;
    sincosd(alpha, &ca, &sa);
    sincosd(beta, &cb, &sb);
    sincosd(gamma, &cg, &sg);
    t1 = 1. - cg;
    m[0] = cg + sa * sa * sb * sb * t1;
    m[1] = -sa * sb * sb * ca * t1 - cb * sg;
    m[2] = sa * sb * cb * t1 - ca * sb * sg;
    m[3] = -sa * sb * sb * ca * t1 + cb * sg;
    m[4] = cg + ca * ca * sb * sb * t1;
    m[5] = -ca * sb * cb * t1 - sa * sb * sg;
    m[6] = sa * sb * cb * t1 + ca * sb * sg;
    m[7] = -ca * sb * cb * t1 + sa * sb * sg;
    m[8] = cg + cb * cb * t1;
}

/** \brief Computes a matrix to rotate a point in R^3 around a given vector.
 *
 * \param m double *res_pt The output.
 * \param n const double *res_pt The vector to be rotated around.
 * \param alpha const double The angle.
 * \return void
 *
 * Assumes a right handed system.
 */
void gen_mat_rot3_n(double *res_pt m, const double *res_pt n, const double alpha)
{
    double ca = cos(alpha),
           sa = sin(alpha),
           cam1 = 1. - ca,
           t1 = n[0] * n[1] * cam1,
           t2 = n[0] * n[2] * cam1,
           t3 = n[1] * n[2] * cam1;
    m[0] = n[0] * n[0] * cam1 + ca;
    m[1] = t1 - n[2] * sa;
    m[2] = t2 + n[1] * sa;
    m[3] = t1 + n[2] * sa;
    m[4] = n[1] * n[1] * cam1 + ca;
    m[5] = t3 - n[0] * sa;
    m[6] = t2 - n[1] * sa;
    m[7] = t3 + n[0] * sa;
    m[8] = n[2] * n[2] * cam1 + ca;
}

/** \brief Rotates a point in R^3 around the origin.
 *
 * \param p double* The point to be rotated.
 * \param axy const double The angle in the xy plane.
 * \param ayz const double The angle in the yz plane.
 * \return void
 *
 * Rotates the point first around the x axis by ayz, then around the z axis by axy.
 * Assumes a right handed system.
 */
void rot3_zx(double *p, const double axy, const double ayz)
{
    double sxy, syz, cxy, cyz, t1, t2, t3[3];
    cp3(t3, p);
    if(axy != 0.)
        sincosd(axy, &cxy, &sxy);
    else
    {
        sxy = 0.;
        cxy = 1.;
    }
    if(ayz != 0.)
        sincosd(ayz, &cyz, &syz);
    else
    {
        syz = 0.;
        cyz = 1.;
    }
    t1 = cyz * t3[1];
    t2 = syz * t3[2];
    p[0] = cxy * t3[0] - sxy * t1 + sxy * t2;
    p[1] = sxy * t3[0] + cxy * t1 - cxy * t2;
    p[2] = syz * t3[1] + cyz * t3[2];
}

/** \brief Rotates and translates a point in R^3 around the origin.
 *
 * \param p double* The point to be transformed.
 * \param axy const double The angle in the xy plane.
 * \param ayz const double The angle in the yz plane.
 * \param t const double* The vector used for the translation.
 * \return void
 *
 * Rotates the point first around the x axis by ayz, then around the z axis by axy and finally translates it.
 * Assumes a right handed system.
 */
void rottrans3_zx(double *p, const double axy, const double ayz, const double *t)
{
    rot3_zx(p, axy, ayz);
    p[0] += t[0];
    p[1] += t[1];
    p[2] += t[2];
}

/** \brief Rotates a point in R^3 around the origin.
 *
 * \param p double* The point to be rotated.
 * \param ayz const double The angle in the yz plane.
 * \param axy const double The angle in the xy plane.
 * \return void
 *
 * Rotates the point first around the z axis by axy, then around the x axis by ayz.
 * Assumes a right handed system.
 */
void rot3_xz(double *p, const double ayz, const double axy)
{
    double sxy, syz, cxy, cyz, t1, t2, t3[3];
    cp3(t3, p);
    if(axy != 0.)
        sincosd(axy, &cxy, &sxy);
    else
    {
        sxy = 0.;
        cxy = 1.;
    }
    if(ayz != 0.)
        sincosd(ayz, &cyz, &syz);
    else
    {
        syz = 0.;
        cyz = 1.;
    }
    t1 = sxy * t3[0];
    t2 = cxy * t3[1];
    p[0] = cxy * t3[0] - sxy * t3[1];
    p[1] = cyz * t1 + cyz * t2 - syz * t3[2];
    p[2] = syz * t1 + syz * t2 + cyz * t3[2];
}

/** \brief Rotates and translates a point in R^3 around the origin.
 *
 * \param p double* The point to be transformed.
 * \param ayz const double The angle in the yz plane.
 * \param axy const double The angle in the xy plane.
 * \param t const double* The vector used for the translation.
 * \return void
 *
 * Rotates the point first around the z axis by axy, then around the x axis by ayz and finally translates it.
 * Assumes a right handed system.
 */
void rottrans3_xz(double *p, const double ayz, const double axy, const double *t)
{
    rot3_xz(p, ayz, axy);
    p[0] += t[0];
    p[1] += t[1];
    p[2] += t[2];
}

/** \brief Rotates a point in R^3 around the origin.
 *
 * \param p double* The point to be rotated.
 * \param beta const double The pitch angle.
 * \param alpha const double The angle for the first rotation around the z axis.
 * \return void
 *
 * Rotates the point first around the z axis by alpha, then around the x' axis by beta (pitch).
 * Assumes a right handed system.
 */
void rot3_xpz(double *p, const double beta, const double alpha)
{
    double sa, sb, ca, cb, t1, t2, t3[3], t4, t5, t6;
    cp3(t3, p);
    if(alpha != 0.)
        sincosd(alpha, &ca, &sa);
    else
    {
        sa = 0.;
        ca = 1.;
    }
    if(beta != 0.)
        sincosd(beta, &cb, &sb);
    else
    {
        sb = 0.;
        cb = 1.;
    }
    t1 = 1. - cb;
    t2 = sa * t1;
    t4 = sa * sb;
    t5 = ca * sb;
    t6 = ca * t2;
    p[0] = (cb + ca * ca * t1) * t3[0] + t6 * t3[1] + t4 * t3[2];
    p[1] = t6 * t3[0] + (cb + sa * t2) * t3[1] - t5 * t3[2];
    p[2] = -t4 * t3[0] + t5 * t3[1] + cb * t3[2];
}

/** \brief Rotates and translates a point in R^3 around the origin.
 *
 * \param p double* The point to be transformed.
 * \param beta const double The pitch angle.
 * \param alpha const double The angle for the first rotation around the z axis.
 * \param t const double* The vector used for the translation.
 * \return void
 *
 * Rotates the point first around the z axis by alpha, then around the x' axis by beta (pitch), then translates it.
 * Assumes a right handed system.
 */
void rottrans3_xpz(double *p, const double beta, const double alpha, const double *t)
{
    rot3_xpz(p, beta, alpha);
    p[0] += t[0];
    p[1] += t[1];
    p[2] += t[2];
}

/** \brief Rotates a point in R^3 around the z axis.
 *
 * \param p double* The point to be rotated.
 * \param axy const double The angle in the xy plane.
 * \return void
 *
 * Rotates the point around the z axis by axy.
 * Assumes a right handed system.
 */
void rot3_z(double *p, const double axy)
{
    double t1[3], sxy, cxy;
    cp3(t1, p);
    if(axy != 0.)
        sincosd(axy, &cxy, &sxy);
    else
    {
        sxy = 0.;
        cxy = 1.;
    }
    p[0] = cxy * t1[0] - sxy * t1[1];
    p[1] = sxy * t1[0] + cxy * t1[1];
}

/** \brief Rotates and translates point in R^3 around the z axis.
 *
 * \param p double* The point to be transformed.
 * \param axy const double The angle in the xy plane.
 * \param t const double* The vector used for the translation.
 * \return void
 *
 * Rotates the point around the z axis by axy, then translates it.
 * Assumes a right handed system.
 */
void rottrans3_z(double *p, const double axy, const double *t)
{
    rot3_z(p, axy);
    p[0] += t[0];
    p[1] += t[1];
    p[2] += t[2];
}

/** \brief Rotates a point in R^3 around the x axis.
 *
 * \param p double* The point to be rotated.
 * \param ayz const double The angle in the yz plane.
 * \return void
 *
 * Rotates the point around the x axis by ayz.
 * Assumes a right handed system.
 */
void rot3_x(double *p, const double ayz)
{
    double t1[3], syz, cyz;
    cp3(t1, p);
    if(ayz != 0.)
        sincosd(ayz, &cyz, &syz);
    else
    {
        syz = 0.;
        cyz = 1.;
    }
    p[1] = cyz * t1[1] - syz * t1[2];
    p[2] = syz * t1[1] + cyz * t1[2];
}

/** \brief Rotates and translates point in R^3 around the x axis.
 *
 * \param p double* The point to be transformed.
 * \param ayz const double The angle in the yz plane.
 * \param t const double* The vector used for the translation.
 * \return void
 *
 * Rotates the point around the x axis by ayz, then translates it.
 * Assumes a right handed system.
 */
void rottrans3_x(double *p, const double ayz, const double *t)
{
    rot3_x(p, ayz);
    p[0] += t[0];
    p[1] += t[1];
    p[2] += t[2];
}

/** \brief Rotates a point in R^3 around the y axis.
 *
 * \param p double* The point to be rotated.
 * \param azx const double The angle in the zx plane.
 * \return void
 *
 * Rotates the point around the y axis by azx.
 * Assumes a right handed system.
 */
void rot3_y(double *p, const double azx)
{
    double t1[3], szx, czx;
    cp3(t1, p);
    if(azx != 0.)
        sincosd(azx, &czx, &szx);
    else
    {
        szx = 0.;
        czx = 1.;
    }
    p[0] = czx * t1[0] + szx * t1[2];
    p[2] = -szx * t1[0] + czx * t1[2];
}

/** \brief Rotates and translates point in R^3 around the y axis.
 *
 * \param p double* The point to be transformed.
 * \param azx const double The angle in the yz plane.
 * \param t const double* The vector used for the translation.
 * \return void
 *
 * Rotates the point around the y axis by azx, then translates it.
 * Assumes a right handed system.
 */
void rottrans3_y(double *p, const double azx, const double *t)
{
    rot3_y(p, azx);
    p[0] += t[0];
    p[1] += t[1];
    p[2] += t[2];
}

/** \brief Rotates a point in R^3 around the origin with Euler angles.
 *
 * \param p double* The point to be rotated.
 * \param gamma const double The third Euler angle.
 * \param beta const double The second Euler angle.
 * \param alpha const double The first Euler angle.
 * \return void
 *
 * Rotates the point first around z axis by alpha, then around the x' axis by beta (pitch) and finally around z' by gamma.
 * Assumes a right handed system.
 */
void rot3_zpxpz(double *p, const double gamma, const double beta, const double alpha)
{
    double sa, sb, ca, cb, sg, cg,
           t1, t2, t3[3], t4, t5, t6, t7;
    cp3(t3, p);
    if(alpha != 0.)
        sincosd(alpha, &ca, &sa);
    else
    {
        sa = 0.;
        ca = 1.;
    }
    if(beta != 0.)
        sincosd(beta, &cb, &sb);
    else
    {
        sb = 0.;
        cb = 1.;
    }
    if(gamma != 0.)
        sincosd(gamma, &cg, &sg);
    else
    {
        sg = 0.;
        cg = 1.;
    }
    t1 = 1. - cg;
    t2 = sa * sb;
    t4 = ca * sb;
    t5 = t2 * sb * ca * t1;
    t7 = cb * t1;
    t6 = t2 * t7;
    p[0] = (cg + sa * t2 * sb * t1) * t3[0] +
           (-t5 - cb * sg) * t3[1] +
           (t6 - t4 * sg) * t3[2];
    p[1] = (-t5 + cb * sg) * t3[0] +
           (cg + ca * t4 * sb * t1) * t3[1] +
           (-t4 * t7 - t2 * sg) * t3[2];
    p[2] = (t6 + t4 * sg) * t3[0] +
           (-t4 * t7 + t2 * sg) * t3[1] +
           (cg + cb * t7) * t3[2];
}

/** \brief Rotates and translates a point in R^3 around the origin with Euler angles.
 *
 * \param p double* The point to be transformed.
 * \param gamma const double The third Euler angle.
 * \param beta const double The second Euler angle.
 * \param alpha const double The first Euler angle.
 * \param t const double* The vector used for the translation.
 * \return void
 *
 * Rotates the point first around z axis by alpha, then around the x' axis by beta (pitch), then around z' by gamma and finally translates it.
 * Assumes a right handed system.
 */
void rottrans3_zpxpz(double *p, const double gamma, const double beta, const double alpha, const double *t)
{
    rot3_zpxpz(p, gamma, beta, alpha);
    p[0] += t[0];
    p[1] += t[1];
    p[2] += t[2];
}

/** \brief Computes the angle between the x axis and the point mapped to the xy plane in R^3.
 *
 * \param p const double* The point.
 * \return double The angle in radians spanning [0, M_PI2].
 *
 * In case that the point coincides with the z axis, return 0xDEAD.
 * Returns the angle between the x axis and the point mapped to xy plane, which is the rotation of it around the z axis.
 */
double angl3_xy(const double *p)
{
    double t1 = p[0], t2 = p[1];
    if(t1 != 0.)
    {
        if(t2 != 0.)
        {
            double t3 = sqrt(t1 * t1 + t2 * t2);
            assert(fabs(t1 / t3) <= 1.);
            if(t2 < 0.)
                return acos(t1 / t3) + M_PI;
            else
                return acos(t1 / t3);
        }
        else
            return t1 > 0. ? 0. : M_PI;
    }
    else if(t2 != 0.)
    {
        if(t2 > 0.)
            return M_PIh;
        else
            return M_PI3h;
    }
    else
        return 0xDEAD;
}

/** \brief Computes the angle between the y axis and the point mapped to the yz plane in R^3.
 *
 * \param p const double* The point.
 * \return double The angle in radians spanning [0, M_PI2].
 *
 * In case that the point coincides with the x axis, return 0xDEAD.
 * Returns the angle between the y axis and the point mapped to yz plane, which is the rotation of it around the x axis.
 */
double angl3_yz(const double *p)
{
    double t1 = p[1], t2 = p[2];
    if(t1 != 0.)
    {
        if(t2 != 0.)
        {
            double t3 = sqrt(t1 * t1 + t2 * t2);
            assert(fabs(t1 / t3) <= 1.);
            if(t2 < 0.)
                return acos(t1 / t3) + M_PI;
            else
                return acos(t1 / t3);
        }
        else
            return t1 > 0. ? 0. : M_PI;
    }
    else if(t2 != 0.)
    {
        if(t2 > 0.)
            return M_PIh;
        else
            return M_PI3h;
    }
    else
        return 0xDEAD;
}

/** \brief Computes the angle between the x axis and the point mapped to the xz plane in R^3.
 *
 * \param p const double* The point.
 * \return double The angle in radians spanning [0, M_PI2].
 *
 * In case that the point coincides with the y axis, return 0xDEAD.
 * Returns the angle between the x axis and the point mapped to xz plane, which is the rotation of it around the y axis.
 */
double angl3_xz(const double *p)
{
    double t1 = p[0], t2 = p[2];
    if(t1 != 0.)
    {
        if(t2 != 0.)
        {
            double t3 = sqrt(t1 * t1 + t2 * t2);
            assert(fabs(t1 / t3) <= 1.);
            if(t2 < 0.)
                return acos(t1 / t3) + M_PI;
            else
                return acos(t1 / t3);
        }
        else return t1 > 0. ? 0. : M_PI;
    }
    else if(t2 != 0.)
    {
        if(t2 > 0.)
            return M_PIh;
        else
            return M_PI3h;
    }
    else
        return 0xDEAD;
}
