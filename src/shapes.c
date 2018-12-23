#include "artpic.h"
#include "alloc.h"
#include "auxf.h"
#include "lina.h"
#include "shapes.h"

/** \brief Generates a table of sine and cosine values.
 *
 * \param s_t double *res_pt The output array of sine values.
 * \param c_t double *res_pt The output array of cosine values.
 * \param res const uint The number of steps from 0 to full_angl.
 * \param full_angl const double The maximum angle to be used for the computation.
 * \return void
 *
 * The input arrays are not allowed to overlap.
 */
void gen_spherical_table(double *res_pt s_t, double *res_pt c_t, const uint res, const double full_angl)
{
    uint i;
    double angl = full_angl / (double)res, t1 = angl;
    s_t = alloc_double(res + 1); /**< Duplicate of the first entry at the end. */
    c_t = alloc_double(res + 1);
    s_t[0] = s_t[res] = 0.;
    c_t[0] = c_t[res] = 1.;
    for(i = 1; i < res; i++, t1 += angl)
        sincosd(t1, &c_t[i], &s_t[i]);
}

/** \brief Generates vertices to render a sphere.
 *
 * \param vsph vertex3* The output array of vertices.
 * \param res_polar const uint The polar resolution, i.e. the steps taken from 0 to pi.
 * \param res_azim const uint The azimuthal resolution, i.e. the steps taken from 0 to 2*pi.
 * \param r const double The radius of the sphere.
 * \return void
 *
 * This is called by gen_patch3_bin_sphere3.
 */
void gen_sphere3_vertex3(vertex3 *vsph, const uint res_polar, const uint res_azim, const double r)
{
    uint i, j, k;
    double *st_a, *ct_a,
           t1, t2, z, tp, tlen;
    st_a = alloc_double(res_azim);
    ct_a = alloc_double(res_azim);
    st_a[0] = 0.;
    ct_a[0] = 1.;
    t2 = t1 = M_PI2 / (double)res_azim;
    for(i = 1; i < res_azim; i++, t2 += t1)
        sincosd(t2, &ct_a[i], &st_a[i]);
    t2 = t1 = M_PI / (double)res_polar;
    for(i = 1, k = 0; i < res_polar; i++, t2 += t1)
    {
        sincosd(t2, &z, &tp); /**< From top to bottom. */
        double zr = z * r,
               tpr = tp * r;
        for(j = 0; j < res_azim; j++, k++)
        {
            vsph[k].n[0] = ct_a[j] * tp;
            vsph[k].n[1] = st_a[j] * tp;
            vsph[k].n[2] = z;
            vsph[k].x[0] = ct_a[j] * tpr;
            vsph[k].x[1] = st_a[j] * tpr;
            vsph[k].x[2] = zr;
            if(fabs((tlen = len_squ3(vsph[k].n)) - 1.) > DBL_EPSILON)
                nvec3_ip(vsph[k].n, sqrt(tlen));
        }
    }
    free(st_a);
    free(ct_a);
}

/** \brief Generates vertices to render a sphere.
 *
 * \param vsph vertex3* The output array of vertices.
 * \param res_polar const uint The polar resolution, i.e. the steps taken from 0 to pi.
 * \param res_azim const uint The azimuthal resolution, i.e. the steps taken from 0 to 2*pi.
 * \param r const double The radius of the sphere.
 * \return void
 *
 * This is called by ddraw_sphere3.
 * This is similar to gen_sphere3_vertex3 but to be used with direct drawing of spheres.
 */
void fgen_sphere3_vertex3(vertex3 *vsph, const uint res_polar, const uint res_azim, const double r)
{
    uint i, j, k;
    double *st_a, *ct_a,
           t1, t2, z, tp;
    st_a = (double *)malloc((res_azim + 1) * sizeof(double));
    ct_a = (double *)malloc((res_azim + 1) * sizeof(double));
    st_a[0] = st_a[res_azim] = 0.;
    ct_a[0] = ct_a[res_azim] = 1.;
    t2 = t1 = M_PI2 / (double)res_azim;
    for(i = 1; i < res_azim; i++, t2 += t1)
        sincosd(t2, &ct_a[i], &st_a[i]);
    vsph[0].n[0] = vsph[0].n[1] = 0.;
    vsph[0].n[2] = 1.;
    vsph[0].x[0] = vsph[0].x[1] = 0.;
    vsph[0].x[2] = r;
    t2 = t1 = M_PI / (double)res_polar;
    for(i = k = 1; i < res_polar; i++, k++, t2 += t1)
    {
        sincosd(t2, &z, &tp);
        double zr = z * r,
               tpr = tp * r;
        for(j = 0; j < res_azim; j++, k++)
        {
            vsph[k].n[0] = ct_a[j] * tp;
            vsph[k].n[1] = st_a[j] * tp;
            vsph[k].n[2] = z;
            vsph[k].x[0] = ct_a[j] * tpr;
            vsph[k].x[1] = st_a[j] * tpr;
            vsph[k].x[2] = zr;
        }
        memcpy(&vsph[k], &vsph[k - res_azim], sizeof(vertex3));
    }
    vsph[k].n[0] = vsph[k].n[1] = 0.;
    vsph[k].n[2] = -1.;
    vsph[k].x[0] = vsph[k].x[1] = 0.;
    vsph[k].x[2] = -r;
    free(st_a);
    free(ct_a);
}

/** \brief Generates an array of patch3s to render a sphere.
 *
 * \param res_polar const uint The polar resolution, i.e. the steps taken from 0 to pi.
 * \param res_azim const uint The azimuthal resolution, i.e. the steps taken from 0 to 2*pi.
 * \param r const double The radius of the sphere.
 * \param nptc uint* The total number of patches.
 * \return patch3* The resulting array of patches.
 *
 * No gl_quad_strip is used to make further addressing of the single patches easier.
 * The resulting array shall be used in gl_lists.
 */
patch3 *gen_patch3_bin_sphere3(const uint res_polar, const uint res_azim, const double r, uint *nptc)
{
    uint i, j, k, l, t1;
    patch3 *ptc;
    vertex3 *vtc;
    vtc = alloc_vertex3((res_polar - 1) * res_azim);
    *nptc = (res_polar - 2) * res_azim + 2;
    ptc = alloc_patch3(*nptc);
    gen_sphere3_vertex3(vtc, res_polar, res_azim, r);
    /**< North pole patch. */
    ptc[0].gl_primitive = GL_TRIANGLE_FAN;
    ptc[0].vt = alloc_vertex3(ptc[0].n_vt = (res_azim + 2));
    ptc[0].vt[0].n[0] = ptc[0].vt[0].n[1] = 0.;
    ptc[0].vt[0].n[2] = 1.;
    ptc[0].vt[0].x[0] = ptc[0].vt[0].x[1] = 0.;
    ptc[0].vt[0].x[2] = r;
    for(i = 1, j = 0; i <= res_azim; i++, j++)
        memcpy(&ptc[0].vt[i], &vtc[j], sizeof(vertex3));
    memcpy(&ptc[0].vt[i], &vtc[0], sizeof(vertex3));
    /**< Main patches. */
    for(i = k = 1, l = 0; i < res_polar - 1; i++)
    {
        for(j = 0; j < res_azim - 1; j++, k++, l++)
        {
            ptc[k].gl_primitive = GL_QUADS;
            ptc[k].vt = alloc_vertex3(ptc[k].n_vt = 4);
            t1 = l + res_azim;
            memcpy(&ptc[k].vt[0], &vtc[l], sizeof(vertex3));
            memcpy(&ptc[k].vt[1], &vtc[t1], sizeof(vertex3));
            memcpy(&ptc[k].vt[2], &vtc[t1 + 1], sizeof(vertex3));
            memcpy(&ptc[k].vt[3], &vtc[l + 1], sizeof(vertex3));
        }
        ptc[k].gl_primitive = GL_QUADS;
        ptc[k].vt = alloc_vertex3(ptc[k].n_vt = 4);
        t1 = l + res_azim;
        memcpy(&ptc[k].vt[0], &vtc[l], sizeof(vertex3));
        memcpy(&ptc[k].vt[1], &vtc[t1], sizeof(vertex3));
        memcpy(&ptc[k].vt[2], &vtc[t1 + 1 - res_azim], sizeof(vertex3));
        memcpy(&ptc[k].vt[3], &vtc[l + 1 - res_azim], sizeof(vertex3));
        k++;
        l++;
    }
    assert(l == (res_polar - 2)*res_azim);
    /**< South pole patch. */
    i = l = (res_polar - 1) * res_azim - 1;
    ptc[k].gl_primitive = GL_TRIANGLE_FAN;
    ptc[k].vt = alloc_vertex3(ptc[k].n_vt = (res_azim + 2));
    ptc[k].vt[0].n[0] = ptc[k].vt[0].n[1] = 0.;
    ptc[k].vt[0].n[2] = -1.;
    ptc[k].vt[0].x[0] = ptc[k].vt[0].x[1] = 0.;
    ptc[k].vt[0].x[2] = -r;
    for(j = 1; j <= res_azim; j++, l--)
        memcpy(&ptc[k].vt[j], &vtc[l], sizeof(vertex3));
    memcpy(&ptc[k].vt[j], &vtc[i], sizeof(vertex3));
    free(vtc);
    return ptc;
}

/** \brief Generates an array of patch3s to render a sphere.
 *
 * \param res_polar const uint The polar resolution, i.e. the steps taken from 0 to pi.
 * \param res_azim const uint The azimuthal resolution, i.e. the steps taken from 0 to 2*pi.
 * \param r const double The radius of the sphere.
 * \param nptc uint* The total number of patches.
 * \return patch3* The resulting array of patches.
 *
 * In contrast to gen_patch3_bin_sphere3 the output is optimized for size and render speed.
 * The resulting array shall be used in gl_lists.
 */
patch3 *gen_patch3_sphere3(const uint res_polar, const uint res_azim, const double r, uint *nptc)
{
    uint i, j, k, l;
    patch3 *ptc;
    vertex3 *vtc;
    vtc = alloc_vertex3((res_polar - 1) * res_azim);
    ptc = alloc_patch3(*nptc = res_polar);
    gen_sphere3_vertex3(vtc, res_polar, res_azim, r);
    /**< North pole patch. */
    ptc[0].gl_primitive = GL_TRIANGLE_FAN;
    ptc[0].vt = alloc_vertex3(ptc[0].n_vt = (res_azim + 2));
    ptc[0].vt[0].n[0] = ptc[0].vt[0].n[1] = 0.;
    ptc[0].vt[0].n[2] = 1.;
    ptc[0].vt[0].x[0] = ptc[0].vt[0].x[1] = 0.;
    ptc[0].vt[0].x[2] = r;
    for(i = 1, j = 0; i <= res_azim; i++, j++)
        memcpy(&ptc[0].vt[i], &vtc[j], sizeof(vertex3));
    memcpy(&ptc[0].vt[i], &vtc[0], sizeof(vertex3));
    /**< Main patch. */
    for(i = 1, l = 0; i < res_polar - 1; i++)
    {
        ptc[i].gl_primitive = GL_QUAD_STRIP;
        ptc[i].vt = alloc_vertex3(ptc[i].n_vt = ((res_azim + 1) << 1));
        for(j = k = 0; j < res_azim; j++, l++)
        {
            memcpy(&ptc[i].vt[k++], &vtc[l], sizeof(vertex3));
            memcpy(&ptc[i].vt[k++], &vtc[l + res_azim], sizeof(vertex3));
        }
        memcpy(&ptc[i].vt[k++], &vtc[l - res_azim], sizeof(vertex3));
        memcpy(&ptc[i].vt[k], &vtc[l], sizeof(vertex3));
    }
    assert(l == (res_polar - 2)*res_azim && i == res_polar - 1);
    /**< South pole patch. */
    k = l = (res_polar - 1) * res_azim - 1;
    ptc[i].gl_primitive = GL_TRIANGLE_FAN;
    ptc[i].vt = alloc_vertex3(ptc[i].n_vt = (res_azim + 2));
    ptc[i].vt[0].n[0] = ptc[i].vt[0].n[1] = 0.;
    ptc[i].vt[0].n[2] = -1.;
    ptc[i].vt[0].x[0] = ptc[i].vt[0].x[1] = 0.;
    ptc[i].vt[0].x[2] = -r;
    for(j = 1; j <= res_azim; j++, k--)
        memcpy(&ptc[i].vt[j], &vtc[k], sizeof(vertex3));
    memcpy(&ptc[i].vt[j], &vtc[l], sizeof(vertex3));
    free(vtc);
    return ptc;
}

/** \brief Directly renders a sphere.
 *
 * \param res_polar const uint The polar resolution, i.e. the steps taken from 0 to pi.
 * \param res_azim const uint The azimuthal resolution, i.e. the steps taken from 0 to 2*pi.
 * \param r const double The radius of the sphere.
 * \return void
 *
 * After each call of this function, the sphere has to be recalculated for another render step.
 */
void ddraw_sphere3(const uint res_polar, const uint res_azim, const double r)
{
    uint i, j, l, e;
    e = (res_azim + 1) * (res_polar - 1) + 2;
    vertex3 *vsph = (vertex3 *)malloc(e * sizeof(vertex3));
    fgen_sphere3_vertex3(vsph, res_polar, res_azim, r);
    /**< North pole patch. */
    glBegin(GL_TRIANGLE_FAN);
    for(i = 0; i <= res_azim; i++)
    {
        glNormal3dv(vsph[i].n);
        glVertex3dv(vsph[i].x);
    }
    glNormal3dv(vsph[1].n);
    glVertex3dv(vsph[1].x);
    glEnd();
    /**< Main patch. */
    for(i = l = 1; i < res_polar - 1; i++)
    {
        glBegin(GL_QUAD_STRIP);
        for(j = 0; j <= res_azim; j++, l++)
        {
            glNormal3dv(vsph[l].n);
            glVertex3dv(vsph[l].x);
            glNormal3dv(vsph[l + res_azim + 1].n);
            glVertex3dv(vsph[l + res_azim + 1].x);
        }
        glEnd();
    }
    /**< South pole patch. */
    j = e - 1;
    glBegin(GL_TRIANGLE_FAN);
    for(i = 0; i <= res_azim; i++, j--)
    {
        glNormal3dv(vsph[j].n);
        glVertex3dv(vsph[j].x);
    }
    glNormal3dv(vsph[e - 2].n);
    glVertex3dv(vsph[e - 2].x);
    glEnd();
    free(vsph);
}

/** \brief Directly draws a Sierpinski sponge.
 *
 * \param n_lvls uint The number of recursions used in the computation to hollow the sponge.
 * \param offset const double* The origin of the sponge.
 * \param scl double A scaling factor to control the size.
 * \return void
 *
 * Taken from freeglut.
 * vertices of the triangle: r0 = (1, 0, 0), r1 = (-1/3, 2*sqrt(2)/3, 0), r2 = (-1/3, -sqrt(2)/3, sqrt(6)/3), r3 = (-1/3, -sqrt(2)/3, -sqrt(6)/3), |r0| = |r1| = |r2| = |r3| = 1.
 * The distance between any two points is 2*sqrt(6)/3.
 */
void sierpinskisponge(uint n_lvls, const double *offset, double scl)
{
    static double tet_r[4][3] = {{1., 0., 0.},
        { -1. / 3., .9428090415820634, 0.},
        { -1. / 3., -.4714045207910316, .8164965809277259},
        { -1. / 3., -.4714045207910316, -.8164965809277259}
    };
    static uint num_tet_fac = 4;
    uint i;
    if(!n_lvls)
    {
        double v[3];
        static uint vert_tet_i[4][3] = {{1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}}; /**< Vertex indices */
        for(i = 0; i < num_tet_fac; i++)
        {
            glBegin(GL_LINE_LOOP);
            glNormal3d(-tet_r[i][0], -tet_r[i][1], -tet_r[i][2]); /**< The unit normals are the negative of the coordinates of the point. */
            uint j;
            for(j = 0; j < 3; j++)
            {
                v[0] = offset[0] + scl * tet_r[vert_tet_i[i][j]][0];
                v[1] = offset[1] + scl * tet_r[vert_tet_i[i][j]][1];
                v[2] = offset[2] + scl * tet_r[vert_tet_i[i][j]][2];
                glVertex3dv(v);
            }
            glEnd();
        }
    }
    else
    {
        double local_offset[3]; /**< Use a local variable to avoid buildup of roundoff errors. */
        n_lvls--;
        scl *= .5;
        for(i = 0; i < num_tet_fac; i++)
        {
            local_offset[0] = offset[0] + scl * tet_r[i][0];
            local_offset[1] = offset[1] + scl * tet_r[i][1];
            local_offset[2] = offset[2] + scl * tet_r[i][2];
            sierpinskisponge(n_lvls, local_offset, scl);
        }
    }
}

/** \brief Draws a 3d like line with an arrow tip.
 *
 * \param x1 const double *res_pt The starting point.
 * \param x2 const double *res_pt The end point.
 * \param len const double The length of the tip, measured in terms of the length of the arrow.
 * \param width const double The width or diameter of the tip, measured in terms of the length of the arrow.
 * \param nsegs uint The number of segments of the tip.
 * \return void
 *
 */
void draw_arrowv(const double *res_pt x1, const double *res_pt x2, const double len, const double width, uint nsegs)
{
#define MAX_NSEGS 20u
    if(!nsegs || nsegs > MAX_NSEGS)
    {
        fprintf(stdout, "! arrow segments greater %u or equal 0, setting to 10\n", MAX_NSEGS);
        nsegs = 10;
    }
    uint i;
    double t1, t2 = len / dist3(x1, x2),
               x3[3],
               v1[3], v2[3], v3[3],
               vt[3] = {1., 0., 0.},
               vn[nsegs + 1][3];
    sub3(x1, x2, x3);
    addnmul3(x2, t2, x3, x3);
    sub3(x2, x3, v1);
    if(fabs(v1[0]) >= fabs(v1[1]))
    {
        vt[0] = 0.;
        vt[1] = 1.;
    }
    cross3(v1, vt, v2);
    t1 = len3(v2);
    if(t1 > 1e-6)
    {
        t2 = width / t1;
        for(i = 0; i < 3; i++)
            v2[i] *= t2;
    }
    cross3(v1, v2, v3);
    t1 = len3(v3);
    if(t1 > 1e-6)
    {
        t2 = width / t1;
        for(i = 0; i < 3; i++)
            v3[i] *= t2;
    }
    static double pifrac_s[(MAX_NSEGS * (1 + MAX_NSEGS)) >> 1],
                  pifrac_c[(MAX_NSEGS * (1 + MAX_NSEGS)) >> 1];
    static char init = 0;
    if(!init)
    {
        uint t = 0;
        for(i = 1; i <= MAX_NSEGS; i++)
        {
            t1 = M_PI2 / i;
            uint j;
            for(j = 0; j < i; j++, t++)
                sincosd(j * t1, &pifrac_c[t], &pifrac_s[t]);
        }
        assert(t == (MAX_NSEGS * (1 + MAX_NSEGS)) >> 1);
        init = 1;
    }
    uint lbound = (((nsegs - 1) * nsegs) >> 1);
    for(i = 0; i < nsegs; i++)
    {
        uint j;
        for(j = 0; j < 3; j++)
#ifdef FP_FAST_FMA
            vn[i][j] = fma(pifrac_c[lbound + i], v2[j], x3[j]) + pifrac_s[lbound + i] * v3[j];
#else
            vn[i][j] = x3[j] + pifrac_c[lbound + i] * v2[j] + pifrac_s[lbound + i] * v3[j];
#endif
    }
    cp3(vn[nsegs], vn[0]);
    glBegin(GL_LINES);
    glVertex3dv(x1);
    glVertex3dv(x3);
    for(i = 0; i < nsegs; i++)
    {
        glVertex3dv(vn[i]);
        glVertex3dv(x2);
        glVertex3dv(vn[i]);
        glVertex3dv(vn[i + 1]);
    }
    glEnd();
#undef MAX_NSEGS
}

/** \brief A 'single point' version of draw_arrowv.
 *
 * \param x1 const double The first entry of the first point.
 * \param x2 const double The second entry of the first point.
 * \param x3 const double The third entry of the first point.
 * \param y1 const double The first entry of the second point.
 * \param y2 const double The second entry of the second point.
 * \param y3 const double The third entry of the second point.
 * \param len const double The length of the tip, measured in terms of the length of the arrow.
 * \param width const double The width or diameter of the tip, measured in terms of the length of the arrow.
 * \param nsegs uint The number of segments of the tip.
 * \return void
 *
 */
void draw_arrow(const double x1, const double x2, const double x3, const double y1, const double y2, const double y3, const double len, const double width, uint nsegs)
{
    const double x[3] = {x1, x2, x3}, y[3] = {y1, y2, y3};
    draw_arrowv(x, y, len, width, nsegs);
}
