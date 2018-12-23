#include "artpic.h"
#include "msg.h"
#include "auxf.h"
#include "lina.h"

/**
 * If not otherwise mentioned, the input has to be in Cartesian coordinates.
 */

/** \brief Computes the Euclidean norm of a vector in R^3.
 *
 * \param x const double* The array.
 * \return double The norm.
 *
 */
double len3(const double *x)
{
    return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

/** \brief Computes the square of the Euclidean norm of a vector in R^3.
 *
 * \param x const double* The array.
 * \return double The norm.
 *
 */
double len_squ3(const double *x)
{
    return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}

/** \brief Computes the inner product of two vectors in R^3.
 *
 * \param x const double *res_pt The first array.
 * \param y const double *res_pt The second array.
 * \return double The scalar output.
 *
 * The input arrays are not allowed to overlap.
 */
double dot3(const double *res_pt x, const double *res_pt y)
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

/** \brief Computes the inner product of two normed vectors in R^3.
 *
 * \param x const double *res_pt The first array.
 * \param y const double *res_pt The second array.
 * \return double The scalar output.
 *
 * The input arrays are not allowed to overlap.
 * Compensates for rounding errors.
 * x and y have to be normed to a unity Euclidean norm.
 */
double save_dot3(const double *res_pt x, const double *res_pt y)
{
    const double t1 = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    if(t1 > 1.)
    {
        if(t1 > 1. + DBL_EPSILON2)
        {
            fprintf(stderr, "deviation: %g\n", t1 - 1.);
            error_msg("input has to be normed", ERR_ARG);
        }
        return 1.;
    }
    else if(t1 < -1.)
    {
        if(t1 < -1. - DBL_EPSILON2)
        {
            fprintf(stderr, "deviation: %g\n", t1 + 1.);
            error_msg("input has to be normed", ERR_ARG);
        }
        return -1.;
    }
    return t1;
}

/** \brief Norms a vector in R^3 with respect to the Euclidean norm.
 *
 * \param x const double *res_pt The input array.
 * \param y double *res_pt The output array.
 * \return void
 *
 * The arrays are not allowed to overlap.
 */
void nvec3(const double *res_pt x, double *res_pt y)
{
    const double len = 1. / sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    y[0] = x[0] * len;
    y[1] = x[1] * len;
    y[2] = x[2] * len;
}

/** \brief Norms a vector in R^3 with respect to a given norm.
 *
 * \param x double* The input array.
 * \param len const double The norm to be used.
 * \return void
 *
 * Inplace operation.
 */
void nvec3_ip(double *x, const double len)
{
    const double t = 1. / len;
    x[0] *= t;
    x[1] *= t;
    x[2] *= t;
}

/** \brief Norms a vector in R^3 with respect to the Euclidean norm.
 *
 * \param x double* The array.
 * \return void
 *
 * Inplace operation.
 */
void envec3_ip(double *x)
{
    const double t = 1. / len3(x);
    x[0] *= t;
    x[1] *= t;
    x[2] *= t;
}

/** \brief Computes the Euclidean distance between two vectors in R^3.
 *
 * \param x const double *res_pt The first array.
 * \param y const double *res_pt The second array.
 * \return double The distance.
 *
 * The arrays are not allowed to overlap.
 */
double dist3(const double *res_pt x, const double *res_pt y)
{
    double z[3] = {x[0] - y[0], x[1] - y[1], x[2] - y[2]},
           t = len_squ3(z);
    if(t > 0.) return sqrt(t);
    else return 0.;
}

/** \brief Computes the squared Euclidean distance between two vectors in R^3.
 *
 * \param x const double *res_pt The first array.
 * \param y const double *res_pt The second array.
 * \return double The distance.
 *
 * The arrays are not allowed to overlap.
 */
double dist_squ3(const double *res_pt x, const double *res_pt y)
{
    double z[3] = {x[0] - y[0], x[1] - y[1], x[2] - y[2]};
    return len_squ3(z);
}

/** \brief Adds two vectors in R^3 and saves the result in a third one.
 *
 * \param a const double *res_pt The first input array.
 * \param b const double *res_pt The second input array.
 * \param c double* The output array.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void add3(const double *res_pt a, const double *res_pt b, double *c)
{
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

/** \brief Subtracts two vectors in R^3 and saves the result in a third one.
 *
 * \param a const double *res_pt The first input array.
 * \param b const double *res_pt The second input array.
 * \param c double* The output array.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void sub3(const double *res_pt a, const double *res_pt b, double *c)
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

/** \brief Multiplies two vectors in R^3 entry-wise and saves the result in a third one.
 *
 * \param a const double *res_pt The first input array.
 * \param b const double *res_pt The second input array.
 * \param c double* The output array.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void mul3(const double *res_pt a, const double *res_pt b, double *c)
{
    c[0] = a[0] * b[0];
    c[1] = a[1] * b[1];
    c[2] = a[2] * b[2];
}

/** \brief Multiplies a scalar with a vector in R^3 and saves the result in a second one.
 *
 * \param a const double The scalar multiplier.
 * \param b const double *res_pt The input array.
 * \param c double* The output array.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void linmul3(const double a, const double *res_pt b, double *c)
{
    c[0] = a * b[0];
    c[1] = a * b[1];
    c[2] = a * b[2];
}

/** \brief Multiplies a scalar with a vector in R^3.
 *
 * \param a const double The scalar multiplier.
 * \param b double* The input array.
 * \return void
 *
 * Inplace operation
 */
void linmul3_ip(const double a, double *b)
{
    b[0] *= a;
    b[1] *= a;
    b[2] *= a;
}

/** \brief Adds a linearly scaled vector in R^3 to a second one and stores the result in a third one.
 *
 * \param a const double *res_pt The additive array.
 * \param l const double The scalar multiplier.
 * \param b const double *res_pt The linearly scaled array.
 * \param c double* The output array.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void addnmul3(const double *res_pt a, const double l, const double *res_pt b, double *c)
{
#ifdef FP_FAST_FMA
    c[0] = fma(l, b[0], a[0]);
    c[1] = fma(l, b[1], a[1]);
    c[2] = fma(l, b[2], a[2]);
#else
    c[0] = a[0] + l * b[0];
    c[1] = a[1] + l * b[1];
    c[2] = a[2] + l * b[2];
#endif
}

/** \brief Adds two linearly scaled vectors in R^3 and stores the result in a third one.
 *
 * \param la const double The first multiplier.
 * \param a const double *res_pt The first linearly scaled array.
 * \param lb const double The second multiplier.
 * \param b const double *res_pt The second linearly scaled array.
 * \param c double* The output array.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void lincomb3(const double la, const double *res_pt a, const double lb, const double *res_pt b, double *c)
{
    c[0] = la * a[0] + lb * b[0];
    c[1] = la * a[1] + lb * b[1];
    c[2] = la * a[2] + lb * b[2];
}

/** \brief Changes the sign of a vector in R^3.
 *
 * \param x const double *res_pt The input array.
 * \param y double *res_pt The output array.
 * \return void
 *
 */
void rev3(const double *res_pt x, double *res_pt y)
{
    y[0] = -x[0];
    y[1] = -x[1];
    y[2] = -x[2];
}

/** \brief Changes the sign of a vector in R^3.
 *
 * \param x double* The array.
 * \return void
 *
 * Inplace operation
 */
void rev3_ip(double *x)
{
    x[0] = -x[0];
    x[1] = -x[1];
    x[2] = -x[2];
}

/** \brief Computes the crossproduct of two vectors in R^3 and stores the result in a third one.
 *
 * \param a const double* The first array.
 * \param b const double* The second array.
 * \param c double *res_pt The output array.
 * \return void
 *
 * The output array is not allowed to overlap with the input arrays.
 */
void cross3(const double *a, const double *b, double *res_pt c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/** \brief Computes the crossproduct of two vectors in R^3.
 *
 * \param a const double *res_pt The first array.
 * \param b double *res_pt The second array which will be used to store the result.
 * \return void
 *
 * Inplace operation.
 * The arrays are not allowed to overlap.
 */
void cross3_ip(const double *res_pt a, double *res_pt b)
{
    double t[3] = {a[1] *b[2] - a[2] *b[1],
                   a[2] *b[0] - a[0] *b[2],
                   a[0] *b[1] - a[1] *b[0]
                  };
    cp3(b, t);
}

/** \brief Computes the normed crossproduct of vectors in R^3 and stores the result in a third one.
 *
 * \param a const double *res_pt The first array.
 * \param b const double *res_pt The second array.
 * \param c double *res_pt The output array.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void normedcross3(const double *res_pt a, const double *res_pt b, double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    double len = len_squ3(c);
    if(fabs(len - 1.) > DBL_EPSILON)
    {
        len = 1. / sqrt(len);
        c[0] *= len;
        c[1] *= len;
        c[2] *= len;
    }
}

/** \brief Computes the normed crossproduct of two vectors in R^3.
 *
 * \param a const double* The first array.
 * \param b double* The second array which will be used to store the result.
 * \return void
 *
 * Inplace operation.
 * The arrays are not allowed to overlap.
 */
void normedcross3_ip(const double *res_pt a, double *res_pt b)
{
    double t[3];
    normedcross3(a, b, t);
    cp3(b, t);
}

/** \brief Computes the triple product of three vectors in R^3.
 *
 * \param a const double* The first array.
 * \param b const double* The second array.
 * \param c const double* The third array.
 * \return double The scalar output.
 *
 * If all vectors are normed and orthogonal, the result will be unity.
 */
double tripprod3(const double *a, const double *b, const double *c)
{
    return a[0] * b[1] * c[2] + a[1] * b[2] * c[0] + a[2] * b[0] * c[1] -
           a[0] * b[2] * c[1] - a[1] * b[0] * c[2] - a[2] * b[1] * c[0];
}

/** \brief Computes the inner product of a matrix and a vector in R^3.
 *
 * \param m const double *res_pt A 3*3 matrix.
 * \param v const double *res_pt A column vector.
 * \param out double *res_pt The output vector.
 * \return void
 *
 * None of the pointers are allowed to overlap.
 */
void inner_product3_mat_vec(const double *res_pt m, const double *res_pt v, double *res_pt out)
{
    uint i, j;
    for(i = 0; i < 3; i++)
        for(out[i] = 0., j = 0; j < 3; j++)
            out[i] += m[j + i * 3] * v[j];
}

/** \brief Computes the inner product of a matrix and a vector in R^3.
 *
 * \param m const double *res_pt A 3*3 matrix.
 * \param v double *res_pt A column vector also used for the output.
 * \return void
 *
 * The input vector will be overwritten.
 */
void inner_product3_mat_vec_ip(const double *res_pt m, double *res_pt v)
{
    double t[3] = {0., 0., 0.};
    uint i, j;
    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
            t[i] += m[j + i * 3] * v[j];
    cp3(v, t);
}

/** \brief Computes the cosine of the angle between two vectors in R^3.
 *
 * \param a const double *res_pt The first array.
 * \param b const double *res_pt The second array.
 * \return double The scalar output.
 *
 */
double cangl3(const double *res_pt a, const double *res_pt b)
{
    double t1 = len_squ3(a), t2 = len_squ3(b);
    if(fabs(t1 - 1.) > DBL_EPSILON || fabs(t2 - 1.) > DBL_EPSILON)
        return dot3(a, b) / sqrt(t1 * t2);
    else
        return dot3(a, b);
}

/** \brief Computes the angle between two vectors in R^3.
 *
 * \param a const double *res_pt The first array.
 * \param b const double *res_pt The second array.
 * \return double The scalar output in radians.
 *
 */
double angl3(const double *res_pt a, const double *res_pt b)
{
    return acos(cangl3(a, b));
}

/** \brief Converts spherical into Cartesian coordinates.
 *
 * \param s const double *res_pt The spherical coordinates.
 * \param c double *res_pt The Cartesian coordinates.
 * \return void
 *
 * The arrays are not allowed to overlap.
 * The spherical coordinates are s = (rad, theta_polar, phi_azimut).
 * The Cartesian coordinates are c = (x, y, z).
 */
void s_to_c3(const double *res_pt s, double *res_pt c)
{
    const double t1 = s[0] * sin(s[1]);
    c[0] = t1 * cos(s[2]);
    c[1] = t1 * sin(s[2]);
    c[2] = s[0] * cos(s[1]);
}

/** \brief Converts spherical into Cartesian coordinates.
 *
 * \param s double* The spherical coordinates.
 * \return void
 *
 * Inplace operation.
 * The spherical coordinates are s = (rad, theta_polar, phi_azimut).
 */
void s_to_c3_ip(double *s)
{
    const double t1 = s[0] * sin(s[1]),
                 t2[3] = {t1 * cos(s[2]), t1 * sin(s[2]), s[0] *cos(s[1])};
    cp3(s, t2);
}

/** \brief Converts Cartesian into spherical coordinates.
 *
 * \param c const double *res_pt The Cartesian coordinates.
 * \param s double *res_pt The spherical coordinates.
 * \return void
 *
 * The arrays are not allowed to overlap.
 * The Cartesian coordinates are c = (x, y, z).
 * The spherical coordinates are s = (rad, theta_polar, phi_azimut).
 */
void c_to_s3(const double *res_pt c, double *res_pt s)
{
    s[0] = len3(c);
    s[1] = acos(c[2] / s[0]);
    s[2] = atan2(c[1], c[0]);
    if(s[2] < 0.)
        s[2] += M_PI2;
}

/** \brief Converts Cartesian into spherical coordinates.
 *
 * \param c const double *res_pt The Cartesian coordinates.
 * \return void
 *
 * Inplace operation.
 * The Cartesian coordinates are c = (x, y, z).
 * The spherical coordinates are s = (rad, theta_polar, phi_azimut).
 */
void c_to_s3_ip(double *c)
{
    double t1[3];
    t1[0] = len3(c);
    t1[1] = acos(c[2] / t1[0]);
    t1[2] = atan2(c[1], c[0]);
    if(t1[2] < 0.)
        t1[2] += M_PI2;
    cp3(c, t1);
}

/** \brief Computes the minimal distance between a point and a vector in R^3.
 *
 * \param p const point3 *res_pt The pointer to the point3 variable.
 * \param l const line3 *res_pt The pointer to the line3 variable.
 * \return double The scalar distance.
 *
 * The input arrays are not allowed to overlap.
 * The distance is computed according to the Euclidean norm.
 */
double dist_poi_lin(const point3 *res_pt p, const line3 *res_pt l)
{
    double t1[3];
    sub3((*p).x, (*l).o, t1);
    cross3_ip((*l).r, t1);
    const double t2 = len3(t1);
    if(t2 == 0.)
        return 0.;
    else
        return t2 / len3((*l).r);
}

/** \brief Computes the distance between two parallel vectors in R^3.
 *
 * \param l1 const line3 *res_pt The pointer to the first line3 variable.
 * \param l2 const line3 *res_pt The pointer to the second line3 variable.
 * \return double The scalar distance.
 *
 * The input arrays are not allowed to overlap.
 * The distance is computed according to the Euclidean norm.
 */
double dist_lin_lin_parallel(const line3 *res_pt l1, const line3 *res_pt l2)
{
    double t1[3];
    sub3((*l2).o, (*l1).o, t1);
    cross3_ip((*l1).r, t1);
    const double t2 = len3(t1);
    if(t2 == 0.)
        return 0.;
    else
        return t2 / len3((*l1).r);
}

/** \brief Computes the minimal distance between two skew vectors in R^3.
 *
 * \param l1 const line3 *res_pt The pointer to the first line3 variable.
 * \param l2 const line3 *res_pt The pointer to the second line3 variable.
 * \return double The scalar distance.
 *
 * The input arrays are not allowed to overlap.
 * The distance is computed according to the Euclidean norm.
 */
double dist_lin_lin_skew(const line3 *res_pt l1, const line3 *res_pt l2)
{
    double t1[3], t2[3];
    sub3((*l2).o, (*l1).o, t1);
    normedcross3((*l1).r, (*l2).r, t2);
    const double t3 = dot3(t2, t1);
    if(t3 > 0.)
        return t3;
    else
        return -t3;
}

/** \brief Computes the minimal distance between a point and a plane in R^3.
 *
 * \param p const point3 *res_pt The pointer to the point3 variable.
 * \param e const plane3 *res_pt The pointer to the plane3 variable.
 * \return double The scalar distance.
 *
 * The input arrays are not allowed to overlap.
 * The distance is computed according to the Euclidean norm.
 */
double dist_poi_pla(const point3 *res_pt p, const plane3 *res_pt e)
{
    double t1[3];
    sub3((*p).x, (*e).o, t1);
    const double t2 = dot3((*e).n, t1);
    if(t2 > 0.)
        return t2 / len3((*e).n);
    else if(t2 < 0.)
        return -t2 / len3((*e).n);
    else
        return 0.;
}

/** \brief Computes the distance between a vector and a plane in R^3.
 *
 * \param l const line3 *res_pt The pointer to the line3 variable.
 * \param e const plane3 *res_pt The pointer to the plane3 variable.
 * \return double The scalar distance.
 *
 * The input arrays are not allowed to overlap.
 * The distance is computed according to the Euclidean norm.
 */
double dist_lin_pla(const line3 *res_pt l, const plane3 *res_pt e)
{
    double t1[3];
    sub3((*l).o, (*e).o, t1);
    const double t2 = dot3((*e).n, t1);
    if(t2 > 0.)
        return t2 / len3((*e).n);
    else if(t2 < 0.)
        return -t2 / len3((*e).n);
    else
        return 0.;
}

/** \brief Computes the distance two planes in R^3.
 *
 * \param e1 const plane3 *res_pt The pointer to the first plane3 variable.
 * \param e2 const plane3 *res_pt The pointer to the second plane3 variable.
 * \return double The scalar distance.
 *
 * The input arrays are not allowed to overlap.
 * The distance is computed according to the Euclidean norm.
 */
double dist_pla_pla(const plane3 *res_pt e1, const plane3 *res_pt e2)
{
    double t1[3];
    sub3((*e2).o, (*e1).o, t1);
    const double t2 = dot3((*e1).n, t1);
    if(t2 > 0.)
        return t2 / len3((*e1).n);
    else if(t2 < 0.)
        return -t2 / len3((*e1).n);
    else
        return 0.;
}
