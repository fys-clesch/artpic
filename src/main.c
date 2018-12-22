/**
 * artpic - a ray tracing program in C
 * Copyright (C) 2018  Clemens Schaefermeier
 * clemens(at)fh-muenster.de
 *
 * This program is free software: you can
 * redistribute it and/or modify it under the
 * terms of the GNU General Public License as
 * published by the Free Software Foundation,
 * either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that
 * it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * give credit where credit is due
 */

#include "artpic.h"
#include "viewer.h"
#include "draw.h"
#include "intersec.h"
#include "auxf.h"
#include "msg.h"
#include "ray.h"
#include "ray_aux.h"
#include "refr_data.h"
#include "tofig.h"
#include "tofile.h"
#include "fromfile.h"
#include "tests.h"
#include "alloc.h"

/**
 ** long-term TODO:
 * substitute assertions with error messages (only run time errors)
 * look for a clean printing structure
 *
 ** maybe TODO:
 * record where TIR occurred (in the particle variable)
 * record how much energy is lost during the propagation in the particles
 * handle the lists in OpenGL in an explicit function for a better overview
 * polarization density: change to compute the stokes vector
 *
 ** TODO:
 * check multiple scattering, as some error was detected in the main loop with uninitialized rays
 * problem with the tracing visualization: the start of a transmitted ray isn't stored properly
 * -> add make trace when updating i1
 */

int main(int argc, char **argv)
{
//testing(); getchar(); exit(EXIT_SUCCESS);
    FP_ERR_CLEAR;
    {
        int i;
        for(i = 0; i < argc; i++)
            fprintf(stdout, "arg[%d]: %s\n", i, argv[i]);
    }
#if RELEASE_BUILD
    print_intro();
#endif
    init_allocations();
    const uchar use_ogl = 1,
                output_results = 0;
    uint nors,
         ui,
         res_azim, res_rad,
         nprtcls,
         res_polar_bin, res_azim_bin;
    double max_rad,
           lambda,
           mu_prtcl;
    point3 init,
           fin,
           polarisation;
    sphere3 screen;
    sphrcl_prtcl *mat_prtcls;
    /**< Set the detection screen: */
    res_polar_bin = 91;
    res_azim_bin = 360;
    bin_hit_screen *b_htscrn = alloc_bin_hit_screen(1, 4., res_polar_bin, res_azim_bin);
    /**< Set the rays and some detection containers: */
    lambda = 632.8e-6;
    max_rad = 3.95e-1;
    res_rad = 2; /* 80 */
    res_azim = 2;
    if(res_azim < MAX_NUMBER_OF_RAYS / res_rad)
        nors = res_azim * res_rad;
    else
        error_msg("there are too many rays to be created - change the data type.", ERR_ARG);
    ray *rs = alloc_ray(nors);
    hit_screen *htscrn = alloc_hit_screen(nors);
    /**< Set the particles: */
    nprtcls = 3;
    mu_prtcl = 1.;
    sphrcl_prtcl *prtcls = alloc_sphrcl_prtcl(nprtcls);
//  set_sphrcl_prtcl(&prtcls[0],rfrct_indx_h2o(lambda,"mm"),mu_prtcl,0.,0.,0.,.4);
    set_sphrcl_prtcl(&prtcls[0], 1.4 + 0.i, mu_prtcl, 0., 0., 0., .4);
    set_sphrcl_prtcl(&prtcls[1], n_vac, mu_prtcl, 0., 1., 0., .42);
    set_sphrcl_prtcl(&prtcls[2], n_vac, mu_prtcl, 0., 0., 1., .41);
    set_sphere3(&screen, 0., 0., 0., b_htscrn[0].rad);
    /**< Initiate the bundle of rays: */
    set_point3(&polarisation, 1., 0., 0.);
    set_point3(&init, -1., 0., 0.);
    set_point3(&fin, 0., 0., 0.);
    gen_startrays_straight(rs, lambda, nors, res_azim, res_rad, max_rad, &init, &fin, &polarisation);
    /**< Initiate the drawing of rays: */
    glray *glrays;
    if(use_ogl)
    {
        glrays = alloc_glray(1, nors);
        copy_ray_to_glray_s(rs, glrays[0].glrs, nors);
        for(ui = 0; ui < nors; ui++)
            make_trace(&rs[ui], &glrays[0], ui);
    }
    FP_ERR_CHECK;
    prog_of_rays(0, INIT_PROG);
#if PARALLEL_PROCESSING
#if OSID == ISWIN32 & RELEASE_BUILD
#warning "on windows xp 32 bit and mingw32 with gcc 4.6.2" \
"there have been found deviations when using the built - in openmp 3.0." \
"deviations have been detected in about 20 percent of the results at the 4th decimal place."
#endif
    /**< Parallel start ==> */
    #pragma omp parallel for default(shared) lastprivate(ui) schedule(guided,2)
#endif
    for(ui = 0; ui < nors; ui++)
    {
        intrsec isec;
        double min_l;
        uint min_l_ind = 0,
             n_subrs = 256,
             n_subhtscrn = 512,
             uj,
             i1, i2, i3, i4,
             j1, j2,
             n_move,
             max_subhtscrn = 1 << 15;
        ray first,
            second;
        memcpy(&first, &rs[ui], sizeof(ray));
        /**
        * First run through the system. The code can be easily shortened by
        * beginning with the inner loop, but if the ray hits no particle,
        * this procedure would cost unnecessary memory allocations.
        */
        for(uj = 0, min_l = DBL_MAX; uj < nprtcls; uj++) /**< Check for the nearest intersection with a particle */
        {
            if(intersec_lin_sph(&prtcls[uj].s, &first.v, 0, &isec))
                if(first.v.l < min_l)
                {
                    min_l = first.v.l;
                    min_l_ind = uj;
                }
        }
        if(intersec_lin_sph(&screen, &first.v, 1, &isec)) /**< Check for the intersection with the screen */
        {
            if(first.v.l < min_l)
            {
                propagate_ray(&first);
                if(use_ogl)
                    make_trace(&first, &glrays[0], ui);
                #pragma omp critical(reg_hit)
                {
                    reg_hit(&first, &htscrn[ui], isec.cangl, FIRSTTIME_HIT_STATE);
                }
                #pragma omp critical(prog_of_rays)
                {
                    prog_of_rays(1, COUNT_PROG);
                }
                #pragma omp atomic
                global_ray_info.count_hit++;
                continue; /**< If the screen was hit, proceed with the next ray */
            }
        }
        else
        {
            #pragma omp critical(reg_hit)
            {
                reg_hit(&first, &htscrn[ui], 0xDEAD, EXCEPTION_STATE);
            }
            #pragma omp atomic
            global_ray_info.count_lost++;
            fprintf(stderr, "ray %u:\n", ui);
            error_msg("direct loss of a ray", ERR_ARG);
            continue;
        }
        /**< If a particle has been hit by the ray before the screen, allocate memory and start the inner loop: */
        ray *subrs = alloc_ray(n_subrs);
        hit_screen *subhtscrn = alloc_hit_screen(n_subhtscrn);
        /**< Handle the first intersection: */
        intersec_lin_sph(&prtcls[min_l_ind].s, &first.v, 1, &isec);
        #pragma omp atomic
        prtcls[min_l_ind].hits++;
        propagate_ray(&first);
        if(set_refr_index(&first, &isec, prtcls[min_l_ind].n, prtcls[min_l_ind].mu))
        {
            if(isec.incdnc == OBLIQUE)
                orthogo_ray_at_intersec(&first, &isec);
            i2 = handle_reflntrans(&first, &second, &isec);
            if(i2)
                memcpy(&subrs[1], &second, sizeof(ray));
            first.hits++;
        }
        else
        {
            i2 = 0;
            propagate_ray_eps(&first);
        }
        if(use_ogl)
            make_trace(&first, &glrays[0], ui);
        memcpy(&subrs[0], &first, sizeof(ray));
        i4 = i2 + 1;
        subrs[i4].lam = 0xDEAD;
        /**< Now, we have at maximum 2 rays from the initial one */
        for(i1 = 0, i3 = n_subrs, j1 = 0, j2 = n_subhtscrn, n_move = 0; i1 <= i2;) /**< --> Start inner loop */
        {
            /**
            * i1: current ray index
            * i2: total rays to work out minus 1
            * i3: allocated rays
            * i4: rays to work out and the next ray to be saved by handle_reflntrans
            * n_move: times of moving subrs by n_subrs
            * j1: current hit_screen index
            * j2: allocated hit_screen variables
            */
            assert(subrs[i1].lam != 0.);
            if(i1 == n_subrs) /**< Save precious memory by overwriting the already handled rays */
            {
                assert(i1 <= i4 && (i1 + i4) <= i3);
                memmove(&subrs[0], &subrs[i1], i4 * sizeof(ray));
                i1 = 0;
                i2 -= n_subrs;
                i4 -= n_subrs;
                n_move++;
            }
            if(i1 + i4 == i3) /**< Check for memory re-allocation of rays */
            {
                subrs = realloc_ray(subrs, 0., i3, n_subrs);
                i3 += n_subrs;
            }
            if(j1 == j2) /**< Check for memory re-allocation of hit detection */
            {
                if(j2 > max_subhtscrn - n_subhtscrn)
                {
                    #pragma omp critical(sort_detection)
                    {
                        sort_detection(subhtscrn, j1, &b_htscrn[0], 0);
                    }
                    j1 = 0;
                }
                else
                {
                    subhtscrn = realloc_hit_screen(subhtscrn, j2, n_subhtscrn);
                    j2 += n_subhtscrn;
                }
            }
            if(get_ray_int(&subrs[i1]) < MINI_INTENSITY) /**< Check boundary condition for intensity drop-out */
            {
                if(use_ogl)
                    make_trace(&subrs[i1], &glrays[0], ui);
                reg_hit(&subrs[i1], &subhtscrn[j1], 0xDEAD, EXHAUSTED_RAY_STATE);
                i1++;
                j1++;
                if(use_ogl)
                    make_trace(&subrs[i1], &glrays[0], ui);
                #pragma omp atomic
                global_ray_info.count_exhstd++;
                continue;
            }
            /**< Check for the nearest intersection with the particles: */
            for(uj = 0, min_l = DBL_MAX; uj < nprtcls; uj++)
            {
                if(intersec_lin_sph(&prtcls[uj].s, &subrs[i1].v, 0, &isec))
                    if(subrs[i1].v.l < min_l)
                    {
                        min_l = subrs[i1].v.l;
                        min_l_ind = uj;
                    }
            }
            /**< Check for a hit on the detector: */
            if(intersec_lin_sph(&screen, &subrs[i1].v, 1, &isec))
                if(subrs[i1].v.l < min_l)
                {
                    propagate_ray(&subrs[i1]);
                    if(use_ogl)
                        make_trace(&subrs[i1], &glrays[0], ui);
                    reg_hit(&subrs[i1], &subhtscrn[j1], isec.cangl, REGULAR_HIT_STATE);
                    i1++;
                    j1++;
                    if(use_ogl)
                        make_trace(&subrs[i1], &glrays[0], ui);
                    #pragma omp atomic
                    global_ray_info.count_hit++;
                    continue;
                }
            if(min_l == DBL_MAX)
            {
                if(use_ogl)
                    make_trace(&subrs[i1], &glrays[0], ui);
                reg_hit(&subrs[i1], &subhtscrn[j1], 0xDEAD, EXCEPTION_STATE);
                i1++;
                j1++;
                if(use_ogl)
                    make_trace(&subrs[i1], &glrays[0], ui);
                #pragma omp atomic
                global_ray_info.count_lost++;
                continue;
            }
            /**< Set the intersection: */
            intersec_lin_sph(&prtcls[min_l_ind].s, &subrs[i1].v, 1, &isec);
            #pragma omp atomic
            prtcls[min_l_ind].hits++;
            /**< Propagate the ray to the intersection: */
            propagate_ray(&subrs[i1]);
            /**< Check if ray is travelling into another medium: */
            if(set_refr_index(&subrs[i1], &isec, prtcls[min_l_ind].n, prtcls[min_l_ind].mu))
            {
                if(isec.incdnc == OBLIQUE)
                    orthogo_ray_at_intersec(&subrs[i1], &isec);
                assert(i1 + i4 < i3);
                assert(subrs[i4].lam == 0xDEAD);
                if(handle_reflntrans(&subrs[i1], &subrs[i4], &isec))
                {
                    if(get_ray_int(&subrs[i4]) < MINI_INTENSITY) /**< This is to increment the event-counter of the particle */
                    {
                        #pragma omp atomic
                        prtcls[min_l_ind].exhstds++;
                    }
                    i2++;
                    i4++;
                    assert(subrs[i4].lam == 0.);
                    subrs[i4].lam = 0xDEAD;
                }
                subrs[i1].hits++;
            }
            else
                propagate_ray_eps(&subrs[i1]);
            if(use_ogl)
                make_trace(&subrs[i1], &glrays[0], ui);
            if(get_ray_int(&subrs[i1]) < MINI_INTENSITY) /**< This is to increment the event-counter of the particle */
            {
                #pragma omp atomic
                prtcls[min_l_ind].exhstds++;
            }
        } /**< End inner loop, all child rays from one starting rays are worked out <-- */
        #pragma omp critical(sort_detection)
        {
            sort_detection(subhtscrn, j1, &b_htscrn[0], 0);
        }
        #pragma omp critical(count_gen_update) /**< This can be substituted by atomic capture with OMP 3.1 */
        {
            global_ray_info.count_gen += (i2 + n_move * n_subrs);
        }
        #pragma omp critical(prog_of_rays)
        {
            prog_of_rays(i2 + 1 + n_move * n_subrs, PRINT_PROG); /**< +1 due to the first ray */
        }
        free(subrs);
        free(subhtscrn);
    } /**< <== Parallel end */
    sort_detection(htscrn, nors, &b_htscrn[0], 1);
    prog_of_rays(0, CLEAR_OUT);
    FP_ERR_CHECK;
    if(output_results)
    {
        print_bin_hit_screen("tenet / opera", &b_htscrn[0], 10, 0, INTENSITY, "testing", 0);
        gnuplot_bin_hit_screen("plot / test", &b_htscrn[0], E_FIELD_AMPLITUDE, 0, 0, 1, (res_polar_bin >> 1) + 1);
        gnuplot_bin_hit_screen("plot / test", &b_htscrn[0], MOD_E_FIELD_AMPLITUDE, 0, 0, 1, (res_polar_bin >> 1) + 1);
        gnuplot_bin_hit_screen("plot / test", &b_htscrn[0], INTENSITY, 0, 0, 1, (res_polar_bin >> 1) + 1);
        gnuplot_bin_hit_screen("plot / test", &b_htscrn[0], POLARISATION_DENSITY, 0, 0, 1, (res_polar_bin >> 1) + 1);
    }
    boundingbox bbox;
    line3 trans_prtcls;
    trans_prtcls = (line3)
    {
        .o = {0., 0., 0.},
        .r = {1., 1., 0.},
        .l = 0.
    };
    mat_prtcls = load_prtcls("mat.dat", rfrct_indx_h2o, lambda, 1., 1., &trans_prtcls, &nprtcls, &bbox, 1);
    if(use_ogl)
    {
        handle_glray(glrays, 1, COPY_RAYS, 0);
        handle_bin_sphere3(&b_htscrn[0], ALLOC_DATA, SINCOS_MAP, INTENSITY);
        handle_prtcls_boxed(mat_prtcls, &bbox, nprtcls, ALLOC_DATA);
        handle_prtcls(&prtcls[0], ALLOC_DATA);
        handle_prtcls(&prtcls[1], ALLOC_DATA);
        go_freeglut(argc, argv);
    }
    if(use_ogl)
        free_glray(glrays, 1);
    free(mat_prtcls);
    free(rs);
    free(htscrn);
    free(prtcls);
    free_bin_hit_screen(b_htscrn, 1);
    print_used_mem();
    return EXIT_SUCCESS;
}
