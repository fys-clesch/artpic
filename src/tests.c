#include "artpic.h"
#include "lina.h"
#include "msg.h"
#include "ray.h"
#include "ray_aux.h"
#include "fresnel.h"
#include "refr_data.h"
#include "tests.h"

#define CHECK_OCCURANCE \
        char id[128];\
        snprintf(id,128,"%.121s_%.5i",fname,line);\
        static char firsttime[1024][128];\
        static uint tested[1024],count=0;\
        uint i,j=1;\
        if(!count) for(i=0;i<1024;i++) tested[i]=0;\
        for(i=0;i<count;i++) /**< Check if the function was tested before */\
            if(!(j=strncmp(id,firsttime[i],128))) {tested[i]++;break;}\
        if(j) strncpy(firsttime[count++],id,128); /**< If not, then make a new entry with the id of the function */
#define PRINT_OCCURANCE \
        if(!j) fprintf(stderr,"anyways, this function was tested %u times before",tested[i]);\
        else fprintf(stderr,"anyways, this function was not tested before");

/** \brief Checks a ray for proper (orthogonal) polarisation
 *
 * \param r const ray *res_pt The ray to be checked
 * \param info const char *res_pt Something to identify the ray
 * \param file const char *res_pt The file in which this function was called
 * \param line int The line where this functions was called
 * \param fname const char *res_pt The name of the function where this functions was called
 * \return uchar Returns 0 if something went wrong
 *
 */
uchar assert_ray(const ray *res_pt r, const char *res_pt info, const char *res_pt file, int line, const char *res_pt fname)
{
    double t1 = fabs(dot3((*r).v.r, (*r).opol)),
           t2 = fabs(dot3((*r).v.r, (*r).ppol)),
           t3 = fabs(dot3((*r).ppol, (*r).opol)),
           t4 = fabs(tripprod3((*r).v.r, (*r).opol, (*r).ppol) - 1.),
           t5 = fabs(len3((*r).v.r) - 1.),
           t6 = fabs(len3((*r).opol) - 1.),
           t7 = fabs(len3((*r).ppol) - 1.);
    static const double cd1 = DBL_EPSILON20,
                        cd2 = DBL_EPSILON20,
                        cd3 = DBL_EPSILON20,
                        cd4 = DBL_EPSILON20,
                        cd5 = DBL_EPSILON20,
                        cd6 = DBL_EPSILON20,
                        cd7 = DBL_EPSILON20;
    uchar c1 = t1 > cd1,
          c2 = t2 > cd2,
          c3 = t3 > cd3,
          c4 = t4 > cd4,
          c5 = t5 > cd5,
          c6 = t6 > cd6,
          c7 = t7 > cd7;
    CHECK_OCCURANCE
    if(c1 | c2 | c3 | c4 | c5 | c6 | c7) fprintf(stderr, "news from ray '%s':\n", info);
    if(c1) fprintf(stderr, "- expected less than %g, got: %g\n-> 'r' and 'o' polarisation are not orthogonal\n", cd1, t1);
    if(c2) fprintf(stderr, "- expected less than %g, got: %g\n-> 'r' and 'p' polarisation are not orthogonal\n", cd2, t2);
    if(c3) fprintf(stderr, "- expected less than %g, got: %g\n-> 'p' and 'o' polarisation are not orthogonal\n", cd3, t3);
    if(c4) fprintf(stderr, "- expected less than %g, got: %g\n-> triple product is not unity\n", cd4, t4);
    if(c5) fprintf(stderr, "- expected less than %g, got: %g\n-> 'r' is not normed\n", cd5, t5);
    if(c6) fprintf(stderr, "- expected less than %g, got: %g\n-> 'o' is not normed\n", cd6, t6);
    if(c7) fprintf(stderr, "- expected less than %g, got: %g\n-> 'p' is not normed\n", cd7, t7);
    if(c1 | c2 | c3 | c4 | c5 | c6 | c7)
    {
        PRINT_OCCURANCE
        error_msg("something went wrong here", file, line, fname);
        return 0;
    }
    return 1;
}

/** \brief Checks two vectors for proper orthogonality and their unity norm
 *
 * \param x const double *res_pt The first vector
 * \param y const double *res_pt The second vector
 * \param infox const char *res_pt Something to identify the test
 * \param infoy const char *res_pt Something to identify the test
 * \param file const char *res_pt The file in which this function was called
 * \param line int The line where this functions was called
 * \param fname const char *res_pt The name of the function where this functions was called
 * \return uchar Returns 0 if something went wrong
 *
 */
uchar assert_orthogo(const double *res_pt x, const double *res_pt y, const char *res_pt infox, const char *res_pt infoy, const char *res_pt file, int line, const char *res_pt fname)
{
    double t1 = dot3(x, y),
           t2 = fabs(t1),
           t3 = fabs(len3(x) - 1.),
           t4 = fabs(len3(y) - 1.);
    static const double cd2 = DBL_EPSILON20,
                        cd3 = DBL_EPSILON20,
                        cd4 = DBL_EPSILON20;
    uchar c1 = t2 > cd2,
          c2 = t3 > cd3,
          c3 = t4 > cd4;
    CHECK_OCCURANCE
    if(c1 | c2 | c3) fprintf(stderr, "%s: news from vectors '%s' and '%s'\n", __func__, infox, infoy);
    if(c1) fprintf(stderr, "- expected less than %g, got: (%c) %g\n-> '%s' and '%s' are not orthogonal\n", cd2, t1 < 0. ? '-' : '+', t2, infox, infoy);
    if(c2) fprintf(stderr, "- expected less than %g, got: %g\n-> '%s' is not normalized\n", cd3, t3, infox);
    if(c3) fprintf(stderr, "- expected less than %g, got: %g\n-> '%s' is not normalized\n", cd4, t4, infoy);
    if(c1 | c2 | c3)
    {
        PRINT_OCCURANCE
        error_msg("something went wrong here", file, line, fname);
        return 0;
    }
    return 1;
}

/** \brief Checks two vectors for proper (anti-)parallelity and their unity norm
 *
 * \param x const double *res_pt The first vector
 * \param y const double *res_pt The second vector
 * \param infox const char *res_pt Something to identify the test
 * \param infoy const char *res_pt Something to identify the test
 * \param file const char *res_pt The file in which this function was called
 * \param line int The line where this functions was called
 * \param fname const char *res_pt The name of the function where this functions was called
 * \return uchar Returns 0 if something went wrong
 *
 */
uchar assert_parallel(const double *res_pt x, const double *res_pt y, const char *res_pt infox, const char *res_pt infoy, const char *res_pt file, int line, const char *res_pt fname)
{
    double t1 = dot3(x, y),
           t2 = fabs(fabs(t1) - 1.),
           t3 = fabs(len3(x) - 1.),
           t4 = fabs(len3(y) - 1.);
    static const double cd2 = DBL_EPSILON20,
                        cd3 = DBL_EPSILON20,
                        cd4 = DBL_EPSILON20;
    uchar c1 = t2 > cd2,
          c2 = t3 > cd3,
          c3 = t4 > cd4;
    CHECK_OCCURANCE
    if(c1 | c2 | c3) fprintf(stderr, "%s: news from vectors '%s' and '%s'\n", __func__, infox, infoy);
    if(c1) fprintf(stderr, "- expected less than %g, got: (%c) %g\n-> '%s' and '%s' are not parallel\n", cd2, t1 < 0. ? '-' : '+', t2, infox, infoy);
    if(c2) fprintf(stderr, "- expected less than %g, got: %g\n-> '%s' is not normalized\n", cd3, t3, infox);
    if(c3) fprintf(stderr, "- expected less than %g, got: %g\n-> '%s' is not normalized\n", cd4, t4, infoy);
    if(c1 | c2 | c3)
    {
        PRINT_OCCURANCE
        error_msg("something went wrong here", file, line, fname);
        return 0;
    }
    return 1;
}

/** \brief Checks two vectors for proper directional parallelity and their unity norm
 *
 * \param x const double *res_pt The first vector
 * \param y const double *res_pt The second vector
 * \param infox const char *res_pt Something to identify the test
 * \param infoy const char *res_pt Something to identify the test
 * \param file const char *res_pt The file in which this function was called
 * \param line int The line where this functions was called
 * \param fname const char *res_pt The name of the function where this functions was called
 * \return uchar Returns 0 if something went wrong
 *
 */
uchar assert_direct_parallel(const double *res_pt x, const double *res_pt y, const char *res_pt infox, const char *res_pt infoy, const char *res_pt file, int line, const char *res_pt fname)
{
    double t1 = dot3(x, y),
           t2 = fabs(t1 - 1.),
           t3 = fabs(len3(x) - 1.),
           t4 = fabs(len3(y) - 1.);
    static const double cd2 = DBL_EPSILON20,
                        cd3 = DBL_EPSILON20,
                        cd4 = DBL_EPSILON20;
    uchar c1 = t2 > cd2,
          c2 = t3 > cd3,
          c3 = t4 > cd4;
    CHECK_OCCURANCE
    if(c1 | c2 | c3) fprintf(stderr, "%s: news from vectors '%s' and '%s'\n", __func__, infox, infoy);
    if(c1) fprintf(stderr, "- expected less than %g, got: (%c) %g\n-> '%s' and '%s' are not parallel\n", cd2, t1 < 0. ? '-' : '+', t2, infox, infoy);
    if(c2) fprintf(stderr, "- expected less than %g, got: %g\n-> '%s' is not normalized\n", cd3, t3, infox);
    if(c3) fprintf(stderr, "- expected less than %g, got: %g\n-> '%s' is not normalized\n", cd4, t4, infoy);
    if(c1 | c2 | c3)
    {
        PRINT_OCCURANCE
        error_msg("something went wrong here", file, line, fname);
        return 0;
    }
    return 1;
}

/** \brief Testing various functions
 *
 * \param void
 * \return void
 *
 */
void testing(void)
{
    ray r, ref, trans;
    uchar tir;
    intrsec i;
    r.v.r[0] = 1.;
    r.v.r[1] = 0.;
    r.v.r[2] = 0.;
    r.opol[0] = 0.;
    r.opol[1] = 0.;
    r.opol[2] = 1.;
    r.oamp = 0.;
    r.pamp = 0.;
    r.lam = 1.;
    r.mu_i = 1.;
    r.trans_child = 0;
    r.travel = 0.;
    r.tir = 0;
    i.normal[0] = 0.;
    i.normal[1] = 1.;
    i.normal[2] = 0.;
    i.mu_f = 1.;
    i.incdnc = OBLIQUE;
    cdoub c_out;
    //r.n_i=1.;
    //i.n_f=1.5+1.1e-7*1.i;
    //i.n_f=1.5+.5*1.i;
    //i.n_f=1.5+0.*1.i;
    i.n_f = 1.;
    r.n_i = 1.5 + 1.1e-7 * 1.i;
    r.n_i = 1.5 + .5 * 1.i;
    r.n_i = 1.5 + 0.*1.i;
    for(r.v.r[0] -= .2, r.v.r[1] = .2;
        r.v.r[1] < 1. && r.v.r[0] > 0.;
        r.v.r[1] += .06, r.v.r[0] -= .06)
    {
        FP_ERR_CLEAR;
        envec3_ip(r.v.r);
        normedcross3(r.v.r, r.opol, r.ppol);
        i.cangl = cangl_in(&r, i.normal);
        c_out = ccangl_out(r.n_i, i.n_f, i.cangl, &tir);
        handle_refl(&r, &ref, &i, c_out);
        handle_trans(&r, &trans, &i, c_out);
        double inangle = acos(i.cangl),
               abss = asin(cabs(r.n_i / i.n_f) * sin(inangle)) * 180 / M_PI,
               reals = asin(creal(r.n_i / i.n_f) * sin(inangle)) * 180 / M_PI,
               outdeg = acos(creal(c_out)) * 180 / M_PI,
               outdot = acos(dot3(i.normal, trans.v.r)) * 180 / M_PI;
        printf("c_in      : %12g c_out      : %12g  %s\n", i.cangl, creal(c_out), (tir ? "TIR" : ""));
        printf("c_in (deg): %12g c_out (deg): %12g\n", inangle * 180 / M_PI, outdeg);
        printf("%12s         abs snell says : %12g\n", "", abss);
        printf("%12s        real snell says : %12g\n", "", reals);
        printf("%12s            outdot says : %12g\n", "", outdot);
        printf("%12s      c_out - abs snell : %12g\n", "", outdeg - abss);
        printf("%12s     c_out - real snell : %12g\n", "", outdeg - reals);
        printf("%12s         outdot - c_out : %12g\n", "", outdot - outdeg);
        printf("%12s     outdot - abs snell : %12g\n", "", outdot - abss);
        printf("%12s            +           : %12g\n", "", outdot + abss);
        printf("%12s     outdot - real snell: %12g\n", "", outdot - reals);
        printf("%12s            +           : %12g\n", "", outdot + reals);
        printf("in[0]     : %12g in[1]      : %12g\n", r.v.r[0], r.v.r[1]);
        printf("ref[0]    : %12g ref[1]     : %12g\n", ref.v.r[0], ref.v.r[1]);
        printf("trans[0]  : %12g trans[1]   : %12g\n\n", trans.v.r[0], trans.v.r[1]);
        if(!(fabs(dot3(ref.v.r, ref.opol)) < DBL_EPSILON) || !(fabs(dot3(ref.v.r, ref.ppol)) < DBL_EPSILON) || !(fabs(dot3(ref.ppol, ref.opol)) < DBL_EPSILON))
        {
            printf("ref-dot: %g\n", dot3(ref.v.r, ref.opol));
            printf("ref-dot: %g\n", dot3(ref.v.r, ref.ppol));
            printf("ref-dot: %g\n", dot3(ref.ppol, ref.opol));
        }
        assert_ray(&trans, "trans", ERR_ARG);
        assert_ray(&ref, "ref", ERR_ARG);
        getchar();
        FP_ERR_CHECK;
    }
    /*cdoub ni=1.5+(0.*1.i),
          no=1.+(0.*1.i);
    //cdoub no=1.5+(0.*1.i),
    //      ni=1.+(0.*1.i);
    double c_in,t1,t2,tep,teo,abs_c_out;
    uchar tir;

    printf("   %12s %12s %12s %12s\n"
           "   %12s %12s %12s %12s\n"
           ,"c_in","ref","trans","sum","----","---","-----","---");
    printf("OPOL: \n");
    for(c_in=.745;c_in<.76;c_in+=.0005)
    {
        cdoub c_out=ccangl_out(ni,no,c_in,&tir);
        abs_c_out=cabs_real(c_out);
        if(tir!=1)
        {
            t1=fres_refl_opol_int(ni,no,1.,1.,c_in,abs_c_out);
            t2=fres_trans_opol_int(ni,no,1.,1.,c_in,abs_c_out);
            teo=cabs(cfres_trans_opol(ni,no,1.,1.,c_in,abs_c_out));
            printf("   %12g %12g %12g %12g %g\n",acos(c_in)*180./M_PI,t1,t2,t1+t2,c_in);
            printf("   %12g %12s %12g\n",abs_c_out*180./M_PI,"",teo);
        }
        else
        {
            t1=fres_refl_opol_int(ni,no,1.,1.,c_in,abs_c_out);
            printf("   %12g %12g %12s %12s %g\n",acos(c_in)*180./M_PI,t1,"TIR","TIR",c_in);
        }
    }
    getchar();
    printf("PPOL: \n");
    for(c_in=.745;c_in<.76;c_in+=.0005)
    {
        cdoub c_out=ccangl_out(ni,no,c_in,&tir);
        abs_c_out=cabs_real(c_out);
        if(tir!=1)
        {
            t1=fres_refl_ppol_int(ni,no,1.,1.,c_in,c_out);
            t2=fres_trans_ppol_int(ni,no,1.,1.,c_in,c_out);
            tep=cabs(cfres_trans_ppol(ni,no,1.,1.,c_in,abs_c_out));
            printf("   %12g %12g %12g %12g %g\n",acos(c_in),t1,t2,t1+t2,c_in);
            printf("   %12g %12s %12g\n",abs_c_out*180./M_PI,"",tep);
        }
        else
        {
            t1=fres_refl_opol_int(ni,no,1.,1.,c_in,c_out);
            printf("   %12g %12g %12s %12s %g\n",acos(c_in),t1,"TIR","TIR",c_in);
        }
    }*/
    /*line3 l1;
    //line3 l2;
    plane3 e;
    //plane3 e2;
    sphere3 s;
    //double t1;
    intrsec isec;

    point3 p;
    p.x[0]=-2.;
    p.x[1]=1.;
    p.x[2]=3.;

    l1.o[0]=1.;
    l1.o[1]=1.;
    l1.o[2]=0.;
    l1.r[0]=3.;l1.r[1]=2.;l1.r[2]=5.;
    l1.l=1.;

    e.n[0]=1/5.;
    e.n[1]=1/3.;
    e.n[2]=1/2.;
    e.o[0]=5.;e.o[1]=0;e.o[2]=0;

    s.o[0]=2.;s.o[1]=2.;s.o[2]=2.;
    s.r=1.;

    //t1=dist_pla_pla(&e,&e2);
    intersec_lin_sph(&s,&l1,1,&isec);
    intersec_lin_pla(&e,&l1,&p);

    print_point3(&p,2);
    print_line3(&l1,2);
    print_plane3(&e,2);*/
    /*cdoub cp1=5+4*_Complex_I,cp2=0.5+8*_Complex_I;
    cp1/=cp2;
    printf("%g + I * %g\n",creal(cp1),cimag(cp1));*/
    /*line3 l1;
    l1.o[0]=1.;
    l1.o[1]=1.;
    l1.o[2]=0.;
    l1.r[0]=3.;
    l1.r[1]=2.;
    l1.r[2]=5.;
    l1.l=1.;
    uint n;
    boundingbox bndngbox;
    load_prtcls("E:/Programming/Mie/artpic/mat.dat",rfrct_indx_h2o,632.8e-6,20.,20.,&l1,&n,&bndngbox,1);*/
    puts("done testing");
    /** Direct draw of some test objects, insert into maindisplay:
     * lhandle_sphere3(20,40,1.,CRED,DRAW_EM_ALL);
     * ddraw_sphere3(20,40,1.);
     */
}
