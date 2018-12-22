#include "artpic.h"
#include "auxf.h"
#include "alloc.h"
#include "msg.h"
#include "tofig.h"

/** \brief Calls gnuplot to print a given polar slice to a ps file
 *
 * \param pname const char *res_pt The name of the ps file
 * \param bhs const bin_hit_screen *res_pt The pointer to the variable which stores the results
 * \param type const bin_hit_print_type The type of the information which should be printed to the ps file
 * \param showplot const char If set to 1, the system tries to open the plot file
 * \param unattndd const char If set to 1, no warning will be displayed in case 'fname' already exists. The file will be overwritten.
 * \param slcs const uint The count of slices which are printed to the ps file
 * \param ... The number of the slices, e.g. 5,10,15
 * \return void
 *
 * The first polar slice is index with 1 and not with 0
 */
void gnuplot_bin_hit_screen(const char *res_pt pname, const bin_hit_screen *res_pt bhs, const bin_hit_print_type type, const uchar showplot, const uchar unattndd, const uint slcs, ...)
{
    static uint run = 1;
    uint i, j,
         *az_slcs = alloc_uint(slcs),
          *pind = alloc_uint(slcs);
    const double ra = M_PI2 / ((double)(*bhs).res_azim);
    double ta;
    struct stat fbuf;
    va_list vl;
    va_start(vl, slcs);
    for(i = 0; i < slcs; i++)
    {
        az_slcs[i] = va_arg(vl, uint);
        if(az_slcs[i] > (*bhs).res_polar)
        {
            error_msg("your polar slice number is not available. choosing highest polar slice.", ERR_ARG);
            az_slcs[i] = (*bhs).res_polar;
        }
        else if(az_slcs[i] == 0)
        {
            error_msg("your polar slice number is not available."
                      " remember: in this function, the first element is indexed with '1'."
                      " choosing polar slice 1.", ERR_ARG);
            az_slcs[i] = 1;
        }
    }
    va_end(vl);
    assert(i == slcs);
    /**< Sort the slices and facilitate the addressing of the them: */
    qsort(az_slcs, slcs, sizeof(uint), comp_int);
    for(j = 0, i = 0; j < (*bhs).nbin; j++)
        if((*bhs).idx[j].ia == az_slcs[i] - 1)
        {
            pind[i++] = j; /**< Get the right index for each polar slice */
            if(i == slcs) break;
        }
    /**< Prepare for proper file creation: */
    if((i = strcspn(pname, "/")) < strlen(pname)) /**< Create a directory to store the file if it doesn't exist */
    {
        uint k = i++;
        for(; i < strlen(pname); i++)
            if(pname[i] == '/')
                k = i;
        char dirname[FILENAME_MAX] = ""; /**< Due to a compiler problem on Windows in snprintf */
        snprintf(dirname, ++k, "%s", pname);
        mkdir(dirname DIRMOD;
    }
    char typeid[12], modfname[FILENAME_MAX1];
    snprintf(typeid, 12, "%s",
             type == E_FIELD_AMPLITUDE ? "amp" :
             type == MOD_E_FIELD_AMPLITUDE ? "modamp" :
             type == INTENSITY ? "int" :
             type == POLARISATION_DENSITY ? "pol" : "err");
fnaming:
    i = (uint)snprintf(modfname, FILENAME_MAX1, "%s_%s_%.0f_%u_%u.ps",
                       pname, typeid, (*bhs).lam * 1e6, slcs, run);
    if(i > FILENAME_MAX)
    {
        fprintf(stderr, "produced filename: '%s'\n", modfname);
        error_msg("maximum filename length exceeded: another file could be overwritten, clear identification unlikely. "
                  "pressing 'Enter' will continue the program.", ERR_ARG);
    }
    if(!stat(modfname, &fbuf))
    {
        if(!unattndd)
        {
            fprintf(stdout,
                    "\nthe file '%s' already exists. shall the new file be renamed? then enter 'y': ",
                    modfname);
            fflush(stdout);
            int yn = getchar_dscrd_rmng();
            if(yn == 'y' || yn == 'Y')
            {
                run++;
                goto fnaming;
            }
            else if(yn == '\n' || yn == 'n' || yn == 'N') /**< This maybe elaborated for more security */
            {}
            else
            {
                error_msg("something went wrong here. renaming", ERR_ARG);
                run++;
                goto fnaming;
            }
        }
        else goto fnaming;
    }
    /**< Create a temporary file to store the data: */
    char tempbuf[] = "temp_data_file";
    FILE *fout = fopen(tempbuf, "wb");
    if(NULL == fout)
    {
        error_msg("can not open temporary file to write data", ERR_ARG);
        free(az_slcs);
        free(pind);
        return;
    }
    else if((*bhs).screen_hit_coor_sys != SPHERICAL_CS)
    {
        error_msg("data has to be in spherical coordinates", ERR_ARG);
        fclose(fout);
        remove(tempbuf);
        free(az_slcs);
        free(pind);
        return;
    }
    /**< Address the quantity which is chosen to be printed: */
    double **pnt;
    char typestr[64], tlong[32];
    switch(type)
    {
    case E_FIELD_AMPLITUDE:
        snprintf(typestr, 64, "e-field amplitude / a.u.");
        snprintf(tlong, 32, "amplitude");
        pnt = (*bhs).amp;
        break;
    case MOD_E_FIELD_AMPLITUDE:
        snprintf(typestr, 64, "modulus of e-field amplitude");
        snprintf(tlong, 32, "mod amplitude / a.u.");
        pnt = (*bhs).mod_amp;
        break;
    case INTENSITY:
        snprintf(typestr, 64, "ray intensity");
        snprintf(tlong, 32, "intensity / a.u.");
        pnt = (*bhs).ray_int;
        break;
    case POLARISATION_DENSITY:
        snprintf(typestr, 64, "polarisation distribution");
        snprintf(tlong, 32, "polarisation / %%");
        pnt = (*bhs).pol_dens;
        break;
    default:
        error_msg("wrong output identifier", ERR_ARG);
        snprintf(typestr, 64, "e-field amplitude / a.u.");
        snprintf(tlong, 32, "amplitude");
        pnt = (*bhs).amp;
    }
    /**< Plot the data in a temporary file (see 'tempbuf') */
    fprintf(stdout, "\nplotting to '%s'... ", modfname);
    fflush(stdout);
    uint k;
    for(k = 0; k < slcs; k++)
    {
        for(j = 0, i = pind[k], ta = 0.; (*bhs).idx[i].ia == (*bhs).idx[i + j].ia; j++)
        {
            fprintf(fout, "%g\t%g\n", ta, pnt[(*bhs).idx[i].ia][(*bhs).idx[i + j].ib]);
            ta += ra;
        }
        if(k < slcs - 1) fprintf(fout, "\n\n");
        else fprintf(fout, "\n");
    }
    fclose(fout);
    /**< Write the gnuplot command file: */
    char tempgnu[] = "gnuplot_artpic_cmd_file";
    FILE *gnuout = fopen(tempgnu, "wb");
    if(NULL == gnuout)
    {
        error_msg("can not open temporary file to gnuplot commands", ERR_ARG);
        remove(tempbuf);
        free(az_slcs);
        free(pind);
        return;
    }
    fprintf(gnuout,
            "set terminal postscript color dashed enhanced lw 1.\n"
            "set timestamp 'artpic - %s'\n"
            "set key top Left title 'slice' nobox\n"
            "set grid back\n"
            "set logscale y\n"
            "set xrange [0:%g]\n"
            "set yrange [1e-15:*]\n"
            "set ylabel '%s'\n"
            "set xlabel 'azimuthal degree'\n"
            "set xtics (0,\
                 '{/Symbol p}/4' pi/4,\
                 '{/Symbol p}/2' pi/2,\
                 '3/4{/Symbol p}' pi*3/4,\
                 '{/Symbol p}' pi,\
                 '5/4{/Symbol p}' pi*5/4,\
                 '3/2{/Symbol p}' pi*3/2,\
                 '7/4{/Symbol p}' pi*7/4,\
                 '2{/Symbol p}' pi*2)\n"
            "set format y '%%.0l{/Symbol \\327}10^{%%L}'\n"
            "set output '%s'\n"
            "plot ",
            VERSION_INFO, M_PI2, tlong, modfname);
    for(i = 0; i < slcs; i++)
    {
        fprintf(gnuout, "'%s' index %u title '%u' with linespoints", tempbuf, i, az_slcs[i]);
        if(i < slcs - 1) fprintf(gnuout, ", ");
        else fprintf(gnuout, "\n");
    }
    fprintf(gnuout, "unset output\n");
    fclose(gnuout);
    i = strlen(modfname) + 12;
    char openbuf[i];
    if(showplot)
    {
        if(OSID == ISLINUX) snprintf(openbuf, i, "evince %s", modfname);
        else
        {
            if(OSID == ISWIN64) snprintf(openbuf, i, "gsview64 %s", modfname);
            else if(OSID == ISWIN32) snprintf(openbuf, i, "gsview32 %s", modfname);
        }
    }
    strncpy(modfname, "gnuplot ", FILENAME_MAX);
    strncat(modfname, tempgnu, FILENAME_MAX);
    /**< Call gnuplot: */
    if(system(modfname)) error_msg("there was some error while calling gnuplot", ERR_ARG);
    else
    {
        fprintf(stdout, "done\n");
        if(showplot) system(openbuf);
    }
    remove(tempgnu);
    remove(tempbuf);
    free(pind);
    free(az_slcs);
}
