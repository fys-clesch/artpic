#include "artpic.h"
#include "auxf.h"
#include "lina.h"
#include "intersec.h"

/** \brief Computes the intersection of two vector in R^3.
 *
 * \param l1 const line3 *res_pt The pointer to the first line3 variable.
 * \param l2 line3 *res_pt The pointer to the second line3 variable.
 * \param p point3 *res_pt The output.
 * \param cangl double *res_pt The cosine of the angle between the vectors.
 * \return uchar indicates if The vectors intersect or not.
 *
 * The arrays are not allowed to overlap.
 * The linear multiplier of the second line3 variable is changed such that it points to the intersection.
 * If no intersection is found, this function returns 0, otherwise 1.
 */
uchar intersec_lin_lin(const line3 *res_pt l1, line3 *res_pt l2, point3 *res_pt p, double *res_pt cangl)
{
    double t1[3], t2[3], t3, t4;
    cross3((*l1).r, (*l2).r, t1);
    if(len_squ3(t1) == 0.)
        return 0;
    sub3((*l2).o, (*l1).o, t2);
    if(dot3(t1, t2) != 0.)
        return 0;
    t3 = (*l1).r[1] / (*l1).r[0];
    t4 = (*l1).o[1] - (*l2).o[1] + t3 * ((*l2).o[0] - (*l1).o[0]);
    t4 /= ((*l2).r[1] - t3 * (*l2).r[0]);
    (*l2).l = t4;
    addnmul3((*l2).o, (*l2).l, (*l2).r, (*p).x);
    *cangl = cangl3((*l1).r, (*l2).r);
    return 1;
}

/** \brief Computes the intersection between a vector and a plane in R^3.
 *
 * \param e const plane3 *res_pt The pointer to the plane3 variable.
 * \param l line3 *res_pt The pointer to the line3 variable.
 * \param p point3 *res_pt The output.
 * \return uchar Indicates if an intersection is found or not.
 *
 * The arrays are not allowed to overlap.
 * The linear multiplier of the line3 variable is changed such that it points to the intersection.
 * If no intersection is found, this function returns 0, otherwise 1.
 */
uchar intersec_lin_pla(const plane3 *res_pt e, line3 *res_pt l, point3 *res_pt p)
{
    double t1[3], t2, t3;
    if(fabs(t3 = dot3((*e).n, (*l).r)) < DBL_EPSILON)
        return 0;
    sub3((*e).o, (*l).o, t1);
    t2 = dot3((*e).n, t1);
    t2 /= t3;
    (*l).l = t2;
    addnmul3((*l).o, (*l).l, (*l).r, (*p).x);
    return 1;
}

/** \brief Computes the intersection between a line and a triangle.
 *
 * \param tria const point3 *res_pt An array of vertices to describe the triangle.
 * \param pin const point3 *res_pt A point which is inside of the triangle.
 * \param l line3 *res_pt The pointer to the line3 variable.
 * \param ret const uchar Controls if the function only checks for an intersection or also returns an output.
 * \param is intrsec *res_pt The output.
 * \return uchar indicates If an intersection is found or not.
 *
 * NOT TESTED YET
 */
uchar intersec_lin_tria(const point3 *res_pt tria, const point3 *res_pt pin, line3 *res_pt l, const uchar ret, intrsec *res_pt is)
{
    plane3 pt;
    double t1[3], t2, t3, t4[3], t5[3];
    /**
     * First, generate a plane from the vertices of the triangle and test if the line crosses it.
     * Then check if the crossing lies within the boundaries of the triangle.
     */
    init_plane3_pform(tria[0].x, tria[1].x, tria[2].x, &pt);
    if(dist_poi_pla(pin, &pt) < 0.)
        rev3_ip(pt.n);
    if(fabs(t3 = dot3(pt.n, (*l).r)) < DBL_EPSILON)
        return 0;
    sub3(pt.o, (*l).o, t1);
    t2 = dot3(pt.n, t1) / t3;
    if(t2 < 0.)
        return 0;
    addnmul3((*l).o, t2, (*l).r, t1);
    sub3(tria[1].x, tria[0].x, t4);
    sub3(tria[2].x, tria[0].x, t5);
    /**
     * Now, two ways are possible: check for an intersection between a line and a plane or a point and a plane.
     * The latter is unsafer but faster: not using it.
     */
    double r1x, r2y, r1y, r2x,
           goy, gry, oy, gox, grx, ox,
           lam1, lam2;
    uint i, j;
    if(t4[0] > t4[1])
        if(t4[0] > t4[2])
            i = 0;
        else
            i = 2;
    else if(t4[1] > t4[2])
        i = 1;
    else
        i = 2;
    if(t5[(i + 1) % 3] > t5[(i + 2) % 3])
        j = (i + 1) % 3;
    else
        j = (i + 2) % 3;
    r1x = t4[i];
    r2x = t5[i];
    gox = (*l).o[i];
    grx = (*l).r[i];
    ox = tria[0].x[i];
    r1y = t4[j];
    r2y = t5[j];
    goy = (*l).o[j];
    gry = (*l).r[j];
    oy = tria[0].x[j];
    lam2 = r1x / (r1x * r2y - r1y * r2x) * (goy + t2 * gry - oy - r1y / r1x * (gox - t2 * grx + ox));
    lam1 = (gox + t2 * grx - ox - lam2 * r2x) / r1x;
    if(lam2 < 0. || lam1 < 0. || lam1 + lam2 > 1.)
        return 0;
    if(ret)
    {
        (*l).l = t2;
        cp3((*is).p, t1);
        cp3((*is).normal, pt.n);
        (*is).cangl = cangl3((*l).r, pt.n);
        assert(fabs((*is).cangl) <= 1.);
        (*is).angl = acos((*is).cangl);
    }
    return 1;
}

/** \brief Computes the intersection between a line and a sphere in R^3.
 *
 * \param e const sphere3 *res_pt The pointer to the sphere3 variable.
 * \param l line3 *res_pt The pointer to the line3 variable.
 * \param ret const uchar Controls if the function only checks for an intersection or also returns an output.
 * \param is intrsec3 *res_pt The output.
 * \return uchar Indicates if an intersection is found or not.
 *
 * The arrays are not allowed to overlap.
 * If ret is set to 1, then return the point p and cosine of the angle and the normal vector the intersection.
 * The variable is will be overwritten in any case.
 * The multiplier to the vector to reach point p is always returned.
 * If there is no intersection or the multiplier is negative, return 0, otherwise 1.
 */
uchar intersec_lin_sph(const sphere3 *res_pt e, line3 *res_pt l, const uchar ret, intrsec *res_pt is)
{
    double t1[3], t2[2];
    sub3((*l).o, (*e).o, (*is).p);
    t1[0] = (*l).r[0] * (*l).r[0] + (*l).r[1] * (*l).r[1] + (*l).r[2] * (*l).r[2];
    t1[1] = (*is).p[0] * (*l).r[0] + (*is).p[1] * (*l).r[1] + (*is).p[2] * (*l).r[2];
    t1[2] = (*is).p[0] * (*is).p[0] + (*is).p[1] * (*is).p[1] + (*is).p[2] * (*is).p[2];
    t1[1] *= 2.;
    t1[2] -= ((*e).r * (*e).r);
    uchar pq = solve_pq(t1[0], t1[1], t1[2], t2, 0);
    if(pq != 2)
    {
        if(pq)
            (*is).incdnc = VERTICAL;
        else
            (*is).incdnc = NONE;
        return 0;
    }
    if(t2[0] < 0. || t2[1] < 0.) /**< Check if the intersection is ahead of l. */
    {
        if(t2[0] < 0. && t2[1] < 0.)
        {
            (*is).incdnc = NONE;
            return 0;
        }
        else if(t2[0] < 0.) t2[0] = t2[1];
    }
    else if(t2[0] > t2[1])
        t2[0] = t2[1];
    (*l).l = t2[0];
    if(ret)
    {
        addnmul3((*l).o, (*l).l, (*l).r, (*is).p); /**< Point of intersection. */
        /**< Attention here: */
        sub3((*e).o, (*is).p, t1); /**< Normal vector points to centre of sphere. */
        double t3 = len_squ3(t1);
        if(fabs(t3 - 1.) > DBL_EPSILON)
            nvec3_ip(t1, sqrt(t3));
        cp3((*is).normal, t1);
        (*is).cangl = cangl3((*l).r, t1); /**< Cosine of angle between l and normal vector. */
        assert(fabs((*is).cangl) <= 1.);
        (*is).angl = acos((*is).cangl);
        if(fabs((*is).cangl) < DBL_EPSILON6)
            (*is).incdnc = VERTICAL;
        else
            (*is).incdnc = OBLIQUE;
    }
    return 1;
}
