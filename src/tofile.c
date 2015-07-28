#include "artpic.h"
#include "auxf.h"
#include "alloc.h"
#include "msg.h"
#include "tofile.h"

/** \brief saves a double array to the file f
 *
 * \param m const double* the double array
 * \param row const uint the range of the first index
 * \param col const uint the range of the second index
 * \param f const char* the name of the file
 * \return void
 *
 * the format of the ouput is such that it can be used to include the file as a const double array
 * row and col are used to help indexing and formating the output
 */
void mem_to_f(const double *m, const uint row, const uint col, const char *f)
{
	FILE *wfile = fopen(f, "w");
	if(NULL == wfile)
	{
		fprintf(stderr, "\a\ncan not open file: %s\n", f);
		getchar();
		exit(EXIT_FAILURE);
	}
	uint i = 0, x, y;
	const uint z = row * col;
	fprintf(wfile, "const double m2c[%u]={", z);
	for(x = 0; x < row; x++)
	{
		uint t = x * col;
		for(y = 0; y < col; y++)
		{
			fprintf(wfile, "%g", m[t + y]);
			i++;
			if(i < z) fprintf(wfile, ",");
		}
	}
	fprintf(wfile, "};");
	assert(i == z);
	fclose(wfile);
}

/** \brief prints results stored in a bin_hit_screen variable in a file
 *
 * \param fname const char *res_pt the name of the file
 * \param bhs const bin_hit_screen *res_pt the bin_hit_screen variable which stores the results
 * \param pprec const uint a flag to specify the printing precision in terms of decimal places
 * \param polar_slice uint the number of the polar slice to be printed. if specified with 0, all are printed.
 * \param type const bin_hit_print_type specifies which information shall be printed to file
 * \param add_info const char *res_pt some additional information which will be printed to the file
 * \param unattndd const uchar if set to 1, no warning will be displayed in case 'fname' already exists. the file will be overwritten.
 * \return void
 *
 */
void print_bin_hit_screen(const char *res_pt fname, const bin_hit_screen *res_pt bhs, const uint pprec, uint polar_slice, const bin_hit_print_type type, const char *res_pt add_info, const uchar unattndd)
{
	if(polar_slice > (*bhs).res_polar)
	{
		error_msg("your polar slice number is not available. choosing highest polar slice.", ERR_ARG);
		polar_slice = (*bhs).res_polar;
	}
	static uint run = 1;
	uint i,
	     fsize = 0,
	     pind = 0;
	const double ra = M_PI2 / ((double)(*bhs).res_azim), rp = M_PI / ((double)(*bhs).res_polar),
	             lam = (*bhs).lam * 1e6;
	double ta, tp;
	if(polar_slice != 0) /**< determine the bins in the chosen azimuthal slice and the size of the output */
	{
		for(i = 0; i < (*bhs).nbin; i++)
			if((*bhs).idx[i].ia == polar_slice - 1) break;
		pind = i;
		for(i = 0; (*bhs).idx[pind + i].ia == (*bhs).idx[pind].ia; i++)
			fsize += sizeof(double); /**< fsize can be used for binary output */
	}
	else fsize = (*bhs).nbin * sizeof(double);
	if((i = strcspn(fname, "/")) < strlen(fname)) /**< create a directory to store the file if it doesn't exist */
	{
		uint k = i++;
		for(; i < strlen(fname); i++)
			if(fname[i] == '/') k = i;
		char dirname[FILENAME_MAX] = ""; /**< due to a compiler problem on Windows in snprintf */
		snprintf(dirname, ++k, "%s", fname);
		mkdir(dirname DIRMOD;
	}
	char typeid[12], modfname[FILENAME_MAX1];
	struct stat fbuf;
	snprintf(typeid, 12, "%s",
	         type == E_FIELD_AMPLITUDE ? "amp" :
	         type == MOD_E_FIELD_AMPLITUDE ? "modamp" :
	         type == INTENSITY ? "int" :
	         type == POLARISATION_DENSITY ? "pol" : "err");
fnaming:
	i = (uint)snprintf(modfname, FILENAME_MAX1, "%s_%s_%.0f_%u.txt", fname, typeid, lam, run);
	if(i > FILENAME_MAX)
	{
		fprintf(stderr, "produced filename: '%s'\n", modfname);
		error_msg("maximum filename length exceeded: another file could be overwritten, clear identification unlikely. "
		          "pressing 'Enter' will continue the program.", ERR_ARG);
	}
	if(!stat(modfname, &fbuf)) /**< check if the file already exists */
	{
		if(!unattndd)
		{
			fprintf(stdout, "\nthe file '%s' already exists. shall the new file be renamed? then enter 'y': ", modfname);
			fflush(stdout);
			int yn = getchar_dscrd_rmng();
			if(yn == 'y' || yn == 'Y')
			{
				run++;
				goto fnaming;
			}
			else if(yn == '\n' || yn == 'n' || yn == 'N') /**< this maybe elaborated for more security */
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
	FILE *fout = fopen(modfname, "w");
	if(NULL == fout)
	{
		error_msg("can not open file to write", ERR_ARG);
		return;
	}
	else if((*bhs).screen_hit_coor_sys != SPHERICAL_CS)
	{
		error_msg("data has to be in spherical coordinates", ERR_ARG);
		fprintf(fout, "dire straits: data in 'print_bin_hit_screen' has to be in polar coordinates");
		fclose(fout);
	}
	else
	{
		double **pnt, tmax = 0., tsum = 0.;
		char typestr[64], tlong[32], tshort[16];
		switch(type)
		{
		case E_FIELD_AMPLITUDE:
			snprintf(typestr, 64, "e-field amplitude");
			snprintf(tlong, 32, "amplitude");
			snprintf(tshort, 16, "amp");
			pnt = (*bhs).amp;
			tmax = (*bhs).amp_max;
			tsum = (*bhs).amp_sum;
			break;
		case MOD_E_FIELD_AMPLITUDE:
			snprintf(typestr, 64, "modulus of e-field amplitude");
			snprintf(tlong, 32, "modulus amplitude");
			snprintf(tshort, 16, "mod_amp");
			pnt = (*bhs).mod_amp;
			tmax = (*bhs).mod_amp_max;
			tsum = (*bhs).mod_amp_sum;
			break;
		case INTENSITY:
			snprintf(typestr, 64, "ray intensity");
			snprintf(tlong, 32, "intensity");
			snprintf(tshort, 16, "ray_int");
			pnt = (*bhs).ray_int;
			tmax = (*bhs).ray_int_max;
			tsum = (*bhs).ray_int_sum;
			break;
		case POLARISATION_DENSITY:
			snprintf(typestr, 64, "polarisation distribution");
			snprintf(tlong, 32, "polarisation");
			snprintf(tshort, 16, "pol");
			pnt = (*bhs).pol_dens;
			tmax = (*bhs).pol_max;
			tsum = (*bhs).pol_sum;
			break;
		default:
			error_msg("wrong output identifier", ERR_ARG);
			snprintf(typestr, 64, "e-field amplitude");
			snprintf(tlong, 32, "amplitude");
			snprintf(tshort, 16, "amp");
			pnt = (*bhs).amp;
			tmax = (*bhs).amp_max;
			tsum = (*bhs).amp_sum;
		}
		char buffer[16];
		time_t tt;
		time(&tt);
		fprintf(fout, "this is artpic %s from %s\nthank you for evoking %s on %s\n",
		        VERSION_INFO, VERSION_DATE, __func__, ctime(&tt));
		snprintf(buffer, 16, "%u", polar_slice);
		if(polar_slice != 0) /**< find the maximum and the sum of the chosen quantity in the chosen slice: */
			for(i = 0; (*bhs).idx[pind].ia == (*bhs).idx[pind + i].ia; i++)
			{
				tsum += pnt[(*bhs).idx[pind].ia][(*bhs).idx[pind + i].ib];
				if(pnt[(*bhs).idx[pind].ia][(*bhs).idx[pind + i].ib] > tmax)
					tmax = pnt[(*bhs).idx[pind].ia][(*bhs).idx[pind + i].ib];
			}
		/**< print the file header: */
		fprintf(fout, "wavelength                     : %g nm\n"
		        "polar resolution of detector   : %u\n"
		        "azimutal resolution of detector: %u\n"
		        "chosen polar slice             : %s\n"
		        "output file number             : %u\n"
		        "output type                    : %s\n"
		        "information                    : %s\n"
		        "additional information         : %s\n"
		        "number of generated rays       : %10lu\n"
		        "          absorbed             : %10lu\n"
		        "          lost                 : %10lu\n"
		        "          detected             : %10lu (%lu %%)\n"
		        "total count of TIR             : %lu\n"
		        "maximum %-23.23s: %g\n"
		        "total %-25.25s: %g\n\n",
		        lam, (*bhs).res_polar, (*bhs).res_azim, (polar_slice != 0) ? buffer : "all",
		        run, typestr, (*bhs).global_info.info, add_info,
		        (*bhs).global_info.count_gen, (*bhs).global_info.count_exhstd,
		        (*bhs).global_info.count_lost, (*bhs).global_info.count_hit,
		        100 * (*bhs).global_info.count_hit / (*bhs).global_info.count_gen,
		        (*bhs).tir,
		        tlong, tmax, tlong, tsum);
		fprintf(stdout, "\nprinting %u byte of data to '%s'... ", fsize, modfname);
		fflush(stdout);
		fprintf(fout, "%*s\t%*s\t%*s\n", pprec, tshort, pprec, "polar", pprec, "azimut");
		/**< print the information, finally: */
		uchar ch = 0;
		uint j;
		if(!polar_slice)
			for(i = j = 0, tp = ta = 0.; i < (*bhs).nbin; i++)
			{
				if(j != (*bhs).idx[i].ia)
				{
					tp += rp;
					ch = 1;
				}
				j = (*bhs).idx[i].ia;
				if((*bhs).idx[i].ib == 0) ta = 0.;
				if(ch)
				{
					fprintf(fout, "\n");
					ch = 0;
				}
				fprintf(fout, "\n%*g\t%*g\t%*g", pprec, pnt[j][(*bhs).idx[i].ib], pprec, tp, pprec, ta);
				ta += ra;
			}
		else
			for(j = 0, ta = 0., tp = (double)((polar_slice - 1) * rp); (*bhs).idx[pind].ia == (*bhs).idx[pind + j].ia; j++)
			{
				fprintf(fout, "\n%*g\t%*g\t%*g", pprec, pnt[(*bhs).idx[pind].ia][(*bhs).idx[pind + j].ib], pprec, tp, pprec, ta);
				ta += ra;
			}
		fclose(fout);
		fprintf(stdout, "done\n");
		run++;
		stat(modfname, &fbuf);
		fsize = (uint)fbuf.st_size;
		fprintf(stdout, "total size of '%s': %u kbyte\n", modfname, fsize >> 10);
	}
}
