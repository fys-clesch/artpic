#include "artpic.h"
#include "msg.h"
#include "alloc.h"
#include "shapes.h"
#include "color.h"
#include "lina.h"
#include "viewer.h"
#include "draw.h"

extern double const rotY,
                    norY,
                    rotX;
extern double const bin_sphere3_alpha;
extern const uchar use_light;
extern uint draw_ray_n;

/** \brief Draw a character string.
 *
 * \param s const char* The string.
 * \return void
 *
 * glutBitmapString() can be used instead but it requires a typecast.
 */
void drawstring(const char *s)
{
    uint len, i;
    len = strlen(s);
    for(i = 0; i < len; i++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, s[i]);
}

/** \brief Draw the minimum and maximum value at the colourbox.
 *
 * \param vmin const double* The maximum value mapped to the colourbox.
 * \param vmax const double* The minimum value mapped to the colourbox.
 * \param slen const uint Maximum string length.
 * \return void
 *
 */
void drawstring_cb(const double *vmin, const double *vmax, const uint slen)
{
    uint i;
    char *c = alloc_char(slen + 1);
    snprintf(c, slen, "%-6.3f", *vmin);
    glRasterPos2d(-1., 0.);
    for(i = 0; i < slen; i++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, c[i]);
    snprintf(c, slen, "%-6.3f", *vmax);
    glRasterPos2d(-1., 1.);
    for(i = 0; i < slen; i++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, c[i]);
    free(c);
}

/** \brief Draw a character string in block setting.
 *
 * \param s const char* The character string.
 * \param linewidth const uint The width of the line in characters.
 * \param sx double The initial x position to start drawing.
 * \param sy double The initial y position to start drawing.
 * \param ix const double The x increment after a linebreak.
 * \param iy const double The y increment after a linebreak.
 * \return void
 *
 */
void drawblockstring2d(const char *s, const uint linewidth, double sx, double sy, const double ix, const double iy)
{
    char *buf = alloc_char(linewidth + 1);
    uint len = strlen(s),
         i, j;
    for(i = 0; i < len;)
    {
        strncpy(buf, &s[i], linewidth);
        if(isspace(buf[0]))
        {
            memmove(buf, &buf[1], linewidth - 1); /**< Clear spaces at the first position. */
            if(i + linewidth < len)
            {
                buf[linewidth - 1] = s[i + linewidth]; /**< Get one additional character due to removing the first position. */
                i++;
            }
        }
        if(isalpha(buf[linewidth - 1]) && isalpha(s[i + linewidth]) &&
                i + linewidth < len)
        {
            buf[linewidth - 1] = '-';
            i += (linewidth - 1);
        }
        else
            i += linewidth;
        glRasterPos2d(sx, sy);
        sx += ix;
        sy += iy;
        for(j = 0; j < strlen(buf); j++)
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13, buf[j]);
    }
    free(buf);
}

/** \brief Draws a xyz coordinate system in the lower left corner.
 *
 * \param void
 * \return void
 *
 */
void draw_coord_ov(void)
{
    const double origin[3] = {0., 0., 0.},
                             x[3] = {.5, 0., 0.},
                                    y[3] = {0., .5, 0.},
                                           z[3] = {0., 0., .5};
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0., 10., 0., 10., -2., 2.);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(1.1, 1.1, 0.);
    glRotated(rotX, 1., 0., 0.);
    glRotated(rotY, norY, 1., 0.);
    draw_arrowv(origin, x, .1, .05, 6);
    draw_arrowv(origin, y, .1, .05, 6);
    draw_arrowv(origin, z, .1, .05, 6);
    glRasterPos3d(.7, 0., 0.);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'x');
    glRasterPos3d(0., .7, 0.);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'y');
    glRasterPos3d(0., 0., .7);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'z');
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}

/** \brief draws A colourbox in the upper right corner.
 *
 * \param free_cb const uchar Frees the allocated memory for the colourbox if set to 1.
 * \return void
 *
 */
void draw_cb(const uchar free_cb)
{
    static colourval *cb_stripe;
    static uchar allocd = 0;
    if(free_cb && allocd)
        free(cb_stripe);
    else if(!free_cb)
    {
        const uint res = 200;
        static uchar chg = 1;
        uint i;
        if(!allocd)
        {
            cb_stripe = alloc_colourval(res);
            allocd = 1;
        }
        if(chg)
        {
            for(i = 0; i < res; i++)
                map_cfun(cb_stripe[i].rgba, i, 0, res, VISIBLE_LIGHT);
            chg = 0;
        }
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        gluOrtho2D(0., 10., 0., 10.);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslated(9.25, 7.25, 0.);
        const double cbheight = 2.5, cbwidth = .5,
                     s = cbheight / res;
        double y = s;
        glBegin(GL_LINE_STRIP);
        for(i = 0; i < res; i++, y += s)
        {
            glColor4dv(cb_stripe[i].rgba);
            glVertex2d(0., y);
            glVertex2d(cbwidth, y);
        }
        glEnd();
        /**< Draw a boundary: */
        glColor3dv(CWHITE);
        glBegin(GL_LINE_LOOP);
        glVertex2d(0., cbheight);
        glVertex2d(0., 0.);
        glVertex2d(cbwidth, 0.);
        glVertex2d(cbwidth, cbheight);
        glEnd();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
    }
}

/** \brief Handles the drawing of rays traced through the system.
 *
 * \param rs const glray* The array of rays which are used in the viewer.
 * \param bundles const uint The number of bundles of rays, e.g. for different wavelength.
 * \param opt const draw_opt A flag which tells this function what to do.
 * \param n_rs const uint The number of rays to draw.
 * \param ... A list of integers to identify the rays.
 * \return void
 *
 */
void handle_glray(const glray *rs, const uint bundles, const draw_opt opt, const uint n_rs, ...)
{
#define MAX_BUNDLES 20u
    static const glray *crs[MAX_BUNDLES];
    static uint allocdbundles = 0;
    if(opt == COPY_RAYS)
    {
        crs[allocdbundles++] = rs;
        if(allocdbundles > bundles || allocdbundles > MAX_BUNDLES)
        {
            error_msg("too many ray bundles shall be stored. setting to maximum number of bundles to avoid crash.", ERR_ARG);
            allocdbundles = MAX_BUNDLES;
        }
    }
    else if(opt == DRAW_SOME)
    {
        va_list vl;
        uint i, crs_i[n_rs];
        va_start(vl, n_rs);
        for(i = 0; i < n_rs; i++)
            crs_i[i] = va_arg(vl, uint); /**< Get the ray children to draw. */
        va_end(vl);
        for(i = 0; i < allocdbundles; i++)
        {
            uint j;
            for(j = 0; j < n_rs; j++)
            {
                uint jr = crs_i[j], /**< jr is the specific ray child to draw. */
                     k, ks = 0, kt = 0; /**< ks counts the number of segments in a child, kt the total number of segments. */
                if(jr >= (*crs[i]).n_glrs)
                {
                    jr %= (*crs[i]).n_glrs;
                    uint ji;
                    uchar skip = 0;
                    for(ji = 0; ji < n_rs; ji++)
                        if(jr == crs_i[ji]) /**< Check if this ray was drawn before. */
                        {
                            skip = 1;
                            break;
                        }
                    if(skip)
                        continue; /**< Do not draw a ray twice. */
                    else
                    {
                        crs_i[j] = jr; /**< Insert this ray in the list. */
                        draw_ray_n = crs_i[j];
                    }
                }
                for(k = 0; k <= (*crs[i]).glrs[jr].n_child; k++)
                {
                    uint l;
                    glBegin(GL_LINE_STRIP);/**< Draw a strip for each ray child. */
                    glColor3dv(*((*crs[i]).glrs[jr].rgba));
                    for(l = 0; l < (*crs[i]).glrs[jr].child[ks]; l++)
                        glVertex3dv((*crs[i]).glrs[jr].trace[kt++].x);
                    glEnd();
                    ks++;
                }
                if(kt != (*crs[i]).glrs[jr].n_trace)
                {
                    print_trace_child_count(&(*crs[i]).glrs[jr]);
                    error_msg("indexing is wrong", ERR_ARG);
                }
            }
        }
    }
    else if(opt == DRAW_EM_ALL)
    {
        uint i;
        for(i = 0; i < allocdbundles; i++)
        {
            uint j;
            for(j = 0; j < (*crs[i]).n_glrs; j++) /**< Run through the number of rays in a bundle: */
            {
                uint k, ks = 0, kt = 0;
                for(k = 0; k <= (*crs[i]).glrs[j].n_child; k++) /**< Run through the total number of child rays: */
                {
                    uint l;
                    glBegin(GL_LINE_STRIP);
                    glColor3dv(*((*crs[i]).glrs[j].rgba));
                    for(l = 0; l < (*crs[i]).glrs[j].child[ks]; l++) /**< Run through the children: */
                        glVertex3dv((*crs[i]).glrs[j].trace[kt++].x);
                    glEnd();
                    ks++; /**< Go to the next index of children. */
                }
                if(kt != (*crs[i]).glrs[j].n_trace)
                {
                    print_trace_child_count(&(*crs[i]).glrs[j]);
                    error_msg("indexing is wrong", ERR_ARG);
                }
            }
        }
    }
    else if(opt == FREE_EM_ALL)
        while(allocdbundles)
            crs[allocdbundles--] = NULL;
    else
        error_msg("wrong specifier", ERR_ARG);
#undef MAX_BUNDLES
}

/** \brief This handles the drawing of a simple sphere.
 *
 * \param res_polar const uint The polar resolution of the sphere.
 * \param res_azim const uint The azimuthal resolution of the sphere.
 * \param rad const double The radius of the sphere.
 * \param rgba const double* The colour which should be used to paint the sphere.
 * \param opt const draw_opt a flag which tells this function what to do.
 * \return void
 *
 * Draws and creates a sphere by a list functionality of OpenGL.
 * At the moment, only for testing purposes. shouldn't be called more than once.
 */
void handle_sphere3(const uint res_polar, const uint res_azim, const double rad, const double *rgba, const draw_opt opt)
{
    static uchar initdlist = 0;
    static uint listid, ptc_count;
    if(opt == DRAW_EM_ALL && !initdlist)
    {
        patch3 *ptc = gen_patch3_sphere3(res_polar, res_azim, rad, &ptc_count);
        listid = glGenLists(1);
        if(!listid)
            error_msg("no empty display list", ERR_ARG);
        uint i;
        glNewList(listid, GL_COMPILE);
        glColor4dv(rgba);
        for(i = 0; i < ptc_count; i++)
        {
            uint j;
            glBegin(ptc[i].gl_primitive);
            for(j = 0; j < ptc[i].n_vt; j++)
            {
                glNormal3dv(ptc[i].vt[j].n);
                glVertex3dv(ptc[i].vt[j].x);
            }
            glEnd();
        }
        glEndList();
        free_patch3(ptc, ptc_count);
        initdlist = 1;
    }
    if(opt == DRAW_EM_ALL)
    {
        if(initdlist)
            glCallList(listid);
        else
            return;
    }
    else if(opt == FREE_EM_ALL)
    {
        if(initdlist)
        {
            glDeleteLists(listid, 1);
            initdlist = 0;
        }
    }
    else
        error_msg("wrong specifier", ERR_ARG);
}

/** \brief Handles the drawing of a bin_hit_screen variable.
 *
 * \param *bhs const bin_hit_screenconst The input, i.e. the results of the ray tracing.
 * \param opt const draw_opt A flag which tells this function what to do.
 * \param cfun const colourfun The function which will be used to map the scalar results into RGB space.
 * \param ptype const bin_hit_print_type The physical quantity which should be printed on the screen/sphere.
 * \return void
 *
 */
void handle_bin_sphere3(const bin_hit_screen *const bhs, const draw_opt opt, const colourfun cfun, const bin_hit_print_type ptype)
{
    static patch3 *ptc;
    static bin_hit_screen const *cbhs;
    static uchar allocd = 0, initdlist = 0, update = 0;
    static uint listid, ptc_count;
    static double t_alpha = DEFAULT_BIN_SPHERE3_ALPHA;
    if(t_alpha != bin_sphere3_alpha)
        update = 1;
    if(opt == ALLOC_DATA && !allocd)
    {
        cbhs = bhs;
        ptc = gen_patch3_bin_sphere3((*bhs).res_polar, (*bhs).res_azim, (*bhs).rad, &ptc_count);
        if((*bhs).nbin < ptc_count)
            error_msg("the patching is incorrect. colouring the sphere will be... i do not know, let us see.", ERR_ARG);
        else if((*bhs).nbin > ptc_count)
        {
            error_msg("there are more bins than patches. exiting.", ERR_ARG);
            exit(EXIT_FAILURE);
        }
        colour_bin_patch3(ptc, ptc_count, bhs, cfun, bin_sphere3_alpha, ptype);
        allocd = 1;
    }
    else if(opt == GEN_LISTS && !initdlist)
    {
gen_list_clause:
        listid = glGenLists(1);
        if(!listid)
            error_msg("no empty display list", ERR_ARG);
        uint i;
        glNewList(listid, GL_COMPILE);
        for(i = 0; i < ptc_count; i++)
        {
            uint j;
            glBegin(ptc[i].gl_primitive);
            for(j = 0; j < ptc[i].n_vt; j++)
            {
                glColor4dv(ptc[i].rgba);
                glNormal3dv(ptc[i].vt[j].n);
                glVertex3dv(ptc[i].vt[j].x);
            }
            glEnd();
        }
        glEndList();
        initdlist = 1;
        if(update)
        {
            glCallList(listid);
            update = 0;
        }
    }
    else if(opt == DRAW_EM_ALL)
    {
        if(!update)
        {
            if(allocd)
            {
                if(!initdlist)
                    goto gen_list_clause;
                else if(initdlist)
                    glCallList(listid);
            }
            else
                return;
        }
        else
        {
            assert(allocd);
            glDeleteLists(listid, 1);
            colour_bin_patch3(ptc, ptc_count, cbhs, cfun, bin_sphere3_alpha, ptype);
            t_alpha = bin_sphere3_alpha;
            goto gen_list_clause;
        }
    }
    else if(opt == FREE_EM_ALL)
    {
        if(allocd)
        {
            free_patch3(ptc, ptc_count);
            allocd = 0;
        }
        if(initdlist)
        {
            glDeleteLists(listid, 1);
            initdlist = 0;
        }
    }
    else
        error_msg("wrong specifier", ERR_ARG);
}

/** \brief Handles the drawing of spherical particles.
 *
 * \param spp const sphrcl_prtcl *const res_pt An array of particles.
 * \param bbox const boundingbox *res_pt The box which encloses all particles.
 * \param np const uint The number of particles.
 * \param opt const draw_opt A flag which tells this function what to do.
 * \return void
 *
 * When called the first time, the address of spp and bbox are copied for the later drawing procedure. No memory has to be allocated.
 */
void handle_prtcls_boxed(const sphrcl_prtcl *const res_pt spp, const boundingbox *res_pt bbox, const uint np, const draw_opt opt)
{
    static const sphrcl_prtcl *cspp;
    static const boundingbox *cbbox;
    static uchar initdlist = 0;
    static uint listid, pcount = 0;
    if(opt == ALLOC_DATA)
    {
        cspp = spp;
        cbbox = bbox;
        pcount = np;
    }
    else if(opt == DRAW_EM_ALL)
    {
        if(!pcount)
            return;
        if(!initdlist)
        {
            uint i;
            double xstr[3], ystr[3], zstr[3];
            sub3((*cbbox).rect[1].x, (*cbbox).rect[0].x, xstr);
            sub3((*cbbox).rect[2].x, (*cbbox).rect[1].x, ystr);
            sub3((*cbbox).rect[4].x, (*cbbox).rect[0].x, zstr);
            addnmul3((*cbbox).rect[0].x, .5, xstr, xstr);
            addnmul3((*cbbox).rect[1].x, .5, ystr, ystr);
            addnmul3((*cbbox).rect[0].x, .5, zstr, zstr);
            listid = glGenLists(1);
            if(!listid)
                error_msg("no empty display list", ERR_ARG);
            glNewList(listid, GL_COMPILE);
            glColor4dv(CGREY);
            glMatrixMode(GL_MODELVIEW);
            for(i = 0; i < pcount; i++)
            {
                glPushMatrix();
                glTranslated(cspp[i].s.o[0], cspp[i].s.o[1], cspp[i].s.o[2]);
                ddraw_sphere3(30, 30, cspp[i].s.r);
                glPopMatrix();
            }
            if(use_light)
                glDisable(GL_LIGHTING);
            glColor3dv(CLIGHTGREY);
            glBegin(GL_LINE_LOOP);
            for(i = 0; i < 4; i++)
                glVertex3dv((*cbbox).rect[i].x);
            glEnd();
            glBegin(GL_LINE_LOOP);
            for(i = 4; i < 8; i++)
                glVertex3dv((*cbbox).rect[i].x);
            glEnd();
            glBegin(GL_LINES);
            for(i = 0; i < 4; i++)
            {
                glVertex3dv((*cbbox).rect[i].x);
                glVertex3dv((*cbbox).rect[i + 4].x);
            }
            glEnd();
            glRasterPos3dv(xstr);
            drawstring((*cbbox).xwidth);
            glRasterPos3dv(ystr);
            drawstring((*cbbox).yheight);
            glRasterPos3dv(zstr);
            drawstring((*cbbox).zdepth);
            glEndList();
            initdlist = 1;
        }
        glCallList(listid);
    }
    else if(opt == FREE_EM_ALL)
    {
        if(initdlist)
        {
            glDeleteLists(listid, 1);
            initdlist = 0;
            cbbox = NULL;
            cspp = NULL;
        }
    }
    else
        error_msg("wrong specifier", ERR_ARG);
}

/** \brief Handles the drawing of spherical particles.
 *
 * \param const sphrcl_prtcl *const res_pt spp The single particle.
 * \param opt const draw_opt A flag which tells this function what to do.
 * \return void
 *
 * Can be called multiple times to draw multiple particles.
 */
void handle_prtcls(const sphrcl_prtcl *const res_pt spp, const draw_opt opt)
{
#define MAX_SINGLE_PARTICLES 20u
    static const sphrcl_prtcl *cspp[MAX_SINGLE_PARTICLES];
    static uchar initdlist = 0;
    static uint listid, pcount = 0;
    if(opt == ALLOC_DATA)
    {
        if(pcount < MAX_SINGLE_PARTICLES)
            cspp[pcount++] = spp;
    }
    else if(opt == DRAW_EM_ALL)
    {
        if(!pcount)
            return;
        if(!initdlist)
        {
            uint i;
            listid = glGenLists(1);
            if(!listid)
                error_msg("no empty display list", ERR_ARG);
            glNewList(listid, GL_COMPILE);
            glColor4dv(CGREY);
            glMatrixMode(GL_MODELVIEW);
            for(i = 0; i < pcount; i++)
            {
                glPushMatrix();
                glTranslated((*cspp[i]).s.o[0], (*cspp[i]).s.o[1], (*cspp[i]).s.o[2]);
                ddraw_sphere3(30, 30, (*cspp[i]).s.r);
                glPopMatrix();
            }
            glEndList();
            initdlist = 1;
        }
        glCallList(listid);
    }
    else if(opt == FREE_EM_ALL)
    {
        if(initdlist)
        {
            glDeleteLists(listid, 1);
            initdlist = 0;
        }
        if(pcount)
        {
            uint i;
            for(i = 0; i < pcount; i++)
                cspp[i] = NULL;
        }
    }
    else
        error_msg("wrong specifier", ERR_ARG);
#undef MAX_SINGLE_PARTICLES
}
