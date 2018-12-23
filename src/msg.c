#include "artpic.h"
#include "alloc.h"
#include "ray_aux.h"
#include "auxf.h"
#include "msg.h"

/** \brief Prints an information about the program to the screen-
 *
 * \param void
 * \return void
 *
 */
void print_intro(void)
{
    char *str = alloc_char(1024);
    snprintf(str, 1024, "\
             artpic - a ray tracing program in C\n\
             Copyright (C) 2018  Clemens Schaefermeier\n\
             clemens(at)fh-muenster.de\n\n\
             This program is free software: you can\n\
             redistribute it and/or modify it under the\n\
             terms of the GNU General Public License as\n\
             published by the Free Software Foundation,\n\
             either version 3 of the License, or (at your\n\
             option) any later version.\n\n\
             This program is distributed in the hope that\n\
             it will be useful, but WITHOUT ANY WARRANTY;\n\
             without even the implied warranty of\n\
             MERCHANTABILITY or FITNESS FOR A PARTICULAR\n\
             PURPOSE.  See the GNU General Public License\n\
             for more details.\n\n\
             %s built on %s\n\
             gcc version %s\n\n",
             VERSION_INFO, VERSION_DATE, __VERSION__);
    fncyprint(str, 5);
    wait(2);
    free(str);
}

/** \brief Prints an error message to stderr.
 *
 * \param msg const char *res_pt The message.
 * \param file const char *res_pt The name of the file where the error occurred.
 * \param line int The line number where the error occurred.
 * \param fname const char *res_pt The name of the function where the error occurred.
 * \return void
 *
 */
void error_msg(const char *res_pt msg, const char *res_pt file, int line, const char *res_pt fname)
{
    fprintf(stderr, "\a\ndire straits in file %s line %i from '%s':\n%s\n", file, line, fname, msg);
    perror("interpreted as");
    errno = 0;
#if STOP_AT_ERR
    getchar();
#endif
}

/** \brief Prints floating point exceptions.
 *
 * \param file const char* The name of the file which called this function.
 * \param line int The line where this function was called.
 * \param fname const char* The name of the function which called this function.
 * \return void
 *
 */
void fp_error(const char *file, int line, const char *fname)
{
    int set_excepts = fetestexcept(FE_ALL_EXCEPT & (!FE_INVALID));
    if(set_excepts & FE_DIVBYZERO) error_msg("there was a division by zero", file, line, fname);
    if(set_excepts & FE_INEXACT) error_msg("sometimes double is not enough", file, line, fname);
    if(set_excepts & FE_INVALID) error_msg("there was an invalid operation", file, line, fname);
    if(set_excepts & FE_OVERFLOW) error_msg("overflow fp error. whatever that means.", file, line, fname);
    if(set_excepts & FE_UNDERFLOW) error_msg("underflow fp error. whatever that means.", file, line, fname);
#if __GNUC__>=4&&__GNUC_MINOR__>=8&&(OSID==ISWIN32||OSID==ISWIN64)
    if(set_excepts & FE_DENORMAL) error_msg("denormal error. whatever that means.", file, line, fname);
#endif
    FP_ERR_CLEAR;
}

/** \brief Prints a ray variable to the screen.
 *
 * \param r const ray *res_pt A pointer to a ray variable.
 * \param c const char *res_pt Some additional string to identify the ray.
 * \return void
 *
 */
void print_ray(const ray *res_pt r, const char *res_pt c)
{
    fprintf(stdout, "ray %s:\nwavelength: %g\n"
            "orthogonal component: %g exp(i %g), oint: %g\n"
            "parallel            : %g exp(i %g), pint: %g\n"
            "medium: n = %g + i %g, mu = %g\n"
            "travel: %g\ntir: %u\nhits: %u\ninfo: %s\nchild:\n",
            c, (*r).lam, (*r).oamp, (*r).ophase, (*r).oint, (*r).pamp, (*r).pphase, (*r).pint,
            creal((*r).n_i), cimag((*r).n_i), (*r).mu_i,
            (*r).travel, (*r).tir, (*r).hits, (*r).info);
    print_line3(&(*r).v, (*r).trans_child, 'r');
    print_doub3((*r).opol, (*r).trans_child, 'o');
    print_doub3((*r).ppol, (*r).trans_child, 'p');
}

/** \brief Prints an intrsec variable to the screen.
 *
 * \param i const intrsec *res_pt A pointer to an intrsec variable.
 * \param c const char *res_pt Some additional string to identify the intrsec.
 * \return void
 *
 */
void print_intrsec(const intrsec *res_pt i, const char *res_pt c)
{
    fprintf(stdout, "intrsec %s:\n"
            "angle     : %g\ncos(angle): %g\n"
            "mu_f      : %g\nn_f       : %g + i %g\n"
            "type      : %s\n\n",
            c, (*i).angl * 180. / M_PI, (*i).cangl, (*i).mu_f, creal((*i).n_f), cimag((*i).n_f),
            (*i).incdnc == NONE ? "none" : (*i).incdnc == OBLIQUE ? "oblique" : (*i).incdnc == VERTICAL ? "vertical" : "right-angled");
}

/** \brief Prints a double matrix to the screen.
 *
 * \param m const double* The matrix.
 * \param row uint The range of the first index.
 * \param col uint The range of the second index.
 * \param c const char* A name to identify the matrix.
 * \return void
 *
 */
void print_doub_mat(const double *m, uint row, uint col, const char *c)
{
    uint i, ii, j;
    fprintf(stdout, "%s:\n", c);
    for(i = 0; i < row; i++)
    {
        for(j = 0, ii = i * col; j < col; j++)
            fprintf(stdout, "%15.4g ", m[ii + j]);
        fprintf(stdout, "\n");
    }
}

/** \brief Prints a double matrix to the screen.
 *
 * \param m const uint* The matrix.
 * \param row uint The range of the first index.
 * \param col uint The range of the second index.
 * \param c const char* A name to identify the matrix.
 * \return void
 *
 */
void print_uint_mat(const uint *m, uint row, uint col, const char *c)
{
    uint i, ii, j;
    fprintf(stdout, "%s:\n", c);
    for(i = 0; i < row; i++)
    {
        for(j = 0, ii = i * col; j < col; j++) fprintf(stdout, "%8u ", m[ii + j]);
        fprintf(stdout, "\n");
    }
}

/** \brief Prints a double[3] variable to the screen.
 *
 * \param p const double* The pointer to the double array.
 * \param i const int An integer to identify the variable.
 * \param c const char A character to identify the variable.
 * \return void
 *
 */
void print_doub3(const double *p, const int i, const char c)
{
    fprintf(stdout, "%8s%15g\n%c[%3i]: %15g\n%8s%15g\n",
            "", p[0], c, i, p[1], "", p[2]);
}

/** \brief Prints a point3 variable to the screen.
 *
 * \param p const point3* The pointer to the point3 variable.
 * \param i const int An integer to identify the variable.
 * \param c const char A character to identify the variable.
 * \return void
 *
 */
void print_point3(const point3 *p, const int i, const char c)
{
    fprintf(stdout, "%8s%15g\n%c[%3i]: %15g\n%8s%15g\n",
            "", (*p).x[0], c, i, (*p).x[1], "", (*p).x[2]);
}

/** \brief Prints a line3 variable to the screen.
 *
 * \param l const line3* The pointer to the line3 variable.
 * \param i const int An integer to identify the variable.
 * \param c const char A character to identify the variable.
 * \return void
 *
 */
void print_line3(const line3 *l, const int i, const char c)
{
    fprintf(stdout, "%8s%15g%3s%15g\n%c[%3i]: %15g + %15g * %g\n%8s%15g%3s%15g\n",
            "", (*l).o[0], "", (*l).r[0],
            c, i, (*l).o[1], (*l).r[1], (*l).l,
            "", (*l).o[2], "", (*l).r[2]);
}

/** \brief Prints a plane3 variable to the screen.
 *
 * \param p const plane3* The pointer to the plane3 variable.
 * \param i const int An integer to identify the variable.
 * \return void
 *
 */
void print_plane3(const plane3 *p, const int i)
{
    fprintf(stdout, "%8s%15g%9s%15g\no[%3i]: %15g n[%3i]: %15g\n%8s%15g%9s%15g\n",
            "", (*p).o[0], "", (*p).n[0],
            i, (*p).o[1], i, (*p).n[1],
            "", (*p).o[2], "", (*p).n[2]);
}

/** \brief Prints an unsigned char in a binary format.
 *
 * \param i uchar The char to be converted.
 * \return void
 *
 */
void gimme_bin(uchar i)
{
    int count = 7;
    char out[8] = "ABCDEFGH";
    while(count >= 0)
    {
        out[count--] = (i & 1) ? '1' : '0';
        i >>= 1;
    }
    fprintf(stdout, "dec: %4i -> bin: %s\n", i, out);
}

/** \brief Handles the printing of the calculation progress.
 *
 * \param n const uint The number of the rays which have been traced (un)successfully at the time of calling this function.
 * \param order const prog_ray_order The precise order to this function.
 * \return void
 *
 * The first call to this function has to be done with INIT_PROG, the last one with CLEAR_OUT.
 */
void prog_of_rays(const uint n, const prog_ray_order order)
{
    static char bar[100], printbar[100];
    static ulong sn = 0;
    static double start_t;
    if(order == INIT_PROG)
    {
#if PARALLEL_PROCESSING
        fprintf(stdout, "\nparallel processing info:"
                " there are %i processors available to this program\n",
                omp_get_num_procs());
#endif
        start_t = omp_get_wtime();
        memset(bar, '-', 100);
#if PRINT_PROG_OF_RAYS
        fprintf(stdout, "\nprogress of driving %lu rays through your setup:\n|%-100.100s|",
                GLOBAL_RAY_INFO.count_gen, printbar); /**< GLOBAL_RAY_INFO.count_gen is a global variable. */
        fflush(stdout);
#endif
    }
    else if(order != CLEAR_OUT)
    {
        sn += n;
        if(order == PRINT_PROG)
        {
            memcpy(printbar, bar, sn * 100 / GLOBAL_RAY_INFO.count_gen);
#if PRINT_PROG_OF_RAYS
            fprintf(stdout, "%c|%-100.100s|", 13, printbar);
            fflush(stdout);
#endif
        }
    }
    else
    {
        assert(sn * 100 / GLOBAL_RAY_INFO.count_gen <= 100);
        memcpy(printbar, bar, sn * 100 / GLOBAL_RAY_INFO.count_gen);
#if PRINT_PROG_OF_RAYS
        fprintf(stdout, "%c|%-100.100s|", 13, printbar);
#endif
        fprintf(stdout,
                "\nhere are the statistics:\n"
                "execution time %.2g s\n"
                "number of generated rays %10lu\n"
                "          absorbed       %10lu\n"
                "          lost           %10lu\n"
                "          detected       %10lu (%lu %%)\n",
                fabs(omp_get_wtime() - start_t),
                GLOBAL_RAY_INFO.count_gen, GLOBAL_RAY_INFO.count_exhstd,
                GLOBAL_RAY_INFO.count_lost, GLOBAL_RAY_INFO.count_hit,
                100 * GLOBAL_RAY_INFO.count_hit / GLOBAL_RAY_INFO.count_gen);
        sn = 0;
    }
}

/** \brief Print a char array in a fancy style.
 *
 * \param str const char* The message.
 * \param d const uint The time in microseconds to sleep between each letter-animation.
 * \return void
 *
 */
void fncyprint(const char *str, const uint d)
{
    int i;
    uint j, len = strlen(str);
    for(j = 0; j < len; j++)
    {
        if(iscntrl(str[j])) i = str[j];
        else if(d <= 100 && str[j] < 122) i = str[j] + 6;
        else if(str[j] > 96) i = 127;
        else if(str[j] <= 96 && str[j] > 64) i = 96;
        else if(str[j] <= 64 && str[j] > 47) i = 64;
        else i = 47;
        while(str[j] != i)
        {
            fprintf(stdout, "%c", i--);
            fflush(stdout);
            usleep(d);
            fprintf(stdout, "%c", 8);
        }
        fprintf(stdout, "%c", str[j]);
    }
    fflush(stdout);
}

/** \brief Prints a glray_s variable to check for proper initialisation.
 *
 * \param const glray_s *const res_pt glr The ray to be checked.
 * \return void
 *
 */
void print_trace_child_count(const glray_s *const res_pt glr)
{
    uint i, count = 0;
    for(i = 0; i <= (*glr).n_child; i++) count += (*glr).child[i];
    fprintf(stdout,
            "total trace count is:                %u\n"
            "in %4u child rays there are traces: %u\n", (*glr).n_trace, (*glr).n_child, count);
}
