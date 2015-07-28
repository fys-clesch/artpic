#include "artpic.h"
#include "msg.h"
#include "alloc.h"
#include "auxf.h"
#include "lina.h"
#include "rot.h"

/** \brief reads an int data matrix file and stores it into memory
 *
 * \param fname const char* the name of the file to be read
 * \param m double* the array to store the data
 * \param row const uint the count of the first index of the matrix
 * \param col const uint the count of the second index of the matrix
 * \return void
 *
 */
void intf_to_mem(const char *fname, double *m, const uint row, const uint col)
{
	FILE *rfile = fopen(fname, "r");
	if(NULL == rfile)
	{
		fprintf(stderr, "\a\ncan not open file: %s\n", fname);
		getchar();
		exit(EXIT_FAILURE);
	}
	uint x = 0, y = 0, i = 0,
	     maxlen = get_nodi(INT_MAX);
	const uint o = row * col;
	int t;
	char format[10];
	snprintf(format, 10, "%%%ud", maxlen);
	while(fscanf(rfile, format, &t) != EOF)
	{
		double z = (double)t;
		m[x * col + y] = z;
		if(y < (col - 1)) y++;
		else
		{
			y = 0;
			x++;
		}
		i++;
	}
	fclose(rfile);
	if(i > o) error_msg("memory access violation", ERR_ARG);
	else if(i < o) error_msg("file smaller than expected", ERR_ARG);
}

/** \brief reads a double data matrix file and stores it into memory
 *
 * \param fname const char* the name of the file to be read
 * \param m double* the array to store the data
 * \param row const uint the count of the first index of the matrix
 * \param col const uint the count of the second index of the matrix
 * \return void
 *
 */
void doubf_to_mem(const char *fname, double *m, const uint row, const uint col)
{
	FILE *rfile = fopen(fname, "r");
	if(NULL == rfile)
	{
		fprintf(stderr, "\a\ncan not open file: %s\n", fname);
		getchar();
		exit(EXIT_FAILURE);
	}
	uint x = 0, y = 0, i = 0,
	     maxlen = get_nodd(DBL_MAX);
	const uint o = row * col;
	double z;
	char format[10];
	snprintf(format, 10, "%%%ulg", maxlen);
	while(fscanf(rfile, format, &z) != EOF)
	{
		m[x * col + y] = z;
		if(y < (col - 1)) y++;
		else
		{
			y = 0;
			x++;
		}
		i++;
	}
	fclose(rfile);
	if(i > o) error_msg("memory access violation", ERR_ARG);
	else if(i < o) error_msg("file smaller than expected", ERR_ARG);
}

/** \brief generate particles from a file
 *
 * \param fname const char *res_pt the name of the file
 * \param (*reffun)(double, const char *) cdoub the function to evaluate the refractive index of the medium. first input is the wavelength, second the length unit, e.g. 'mm'.
 * \param wavlen const double the wavelength in milimeters
 * \param xspread const double the center to center distance in the x direction
 * \param yspread const double the center to center distance in the y direction
 * \param tvec const line3 *res_pt the vector along which the particles are displaced. the particles thereby placed in a plane orthogonal to this vector, starting from the x-y plane.
 * \param prtcls uint *res_pt the number of particles which have been loaded from the file
 * \param bbox boundingbox *res_pt the boundingbox variable which is written in the cause of this function
 * \param chatty const uchar a flag to control the output printed on screen during the execution
 * \return sphrcl_prtcl* an array of particles which is created due to this function
 *
 * the file must contain a matrix whose entries correspond to binarized particles
 * a 0 in the file is a particle
 * the particles are initialized with a refractive index according to the given wavelength
 * the plane on which the particles are generated will be transformed along the vector
 */
sphrcl_prtcl *load_prtcls(const char *res_pt fname,
						cdoub(*reffun)(double, const char *),
						const double wavlen,
						const double xspread, const double yspread,
						const line3 *res_pt tvec,
						uint *res_pt prtcls, boundingbox *res_pt bbox,
						const uchar chatty)
{
	FILE *rfile = fopen(fname, "rb");
	if(NULL == rfile)
	{
		fprintf(stderr, "\a\ncan not open file: %s\n", fname);
		getchar();
		exit(EXIT_FAILURE);
	}
	uint i = 0, j = 0, k = 0, o = 0,
	     np, ti, col = 0, row, *mat;
	char tc;
	if(chatty) fprintf(stdout, "\nloading particle matrix from '%s'... ", fname);
	fflush(stdout);
	while(fscanf(rfile, "%c", &tc) != EOF) /**< determine the length of a column */
		if(tc == '\n')
		{
			if(j) break; /**< this happens in case of an additional linebreak at the end of the file */
			else if(!o) col = o = i; /**< first linebreak */
			else if((o += col) != i) error_msg("reading error, columns are not constant", ERR_ARG);
			j = 1;
		}
		else if(isdigit(tc))
		{
			i++; /**< the counter */
			j = 0;
		}
		else if(tc == '\t') j = 0;
	if(!j) o += col; /**< adds the last column to the total number */
	mat = alloc_uint(o);
	rewind(rfile);
	/**< now really extract data from the file: */
	uint maxlen = get_nodui(UINT_MAX);
	char format[10];
	snprintf(format, 10, "%%%uu", maxlen);
	np = i = j = 0;
	while(fscanf(rfile, format, &ti) != EOF)
	{
		if(!ti) np++;
		mat[i * col + j] = ti;
		if(j < (col - 1)) j++;
		else
		{
			j = 0;
			i++;
		}
	}
	if(i * col + j != o) error_msg("read the wrong number of entries. check the input file for linefeeds", ERR_ARG);
	row = o / col;
	fclose(rfile);
	if(chatty)
	{
		fprintf(stdout, "smoothly acquired a %u by %u matrix\n", col, row);
		print_uint_mat(mat, o / col, col, "particle matrix");
	}
	*prtcls = np;
	sphrcl_prtcl *nprtcls = alloc_sphrcl_prtcl(np);
	cdoub ri = reffun(wavlen, "mm");
	double mu_prtcl = 1.,
	       alpha, beta,
	       startx, starty, startz;
	const double psize = (xspread / col < yspread / row) ? (xspread / (2.*col)) : (yspread / (2.*row)); /**< this determines the particle's radius */
	startx = -xspread / 2. + ((col & 1) ? 0. : (xspread / (2.*col)));
	starty = yspread / 2. - ((row & 1) ? 0. : (yspread / (2.*row)));
	startz = 0.;
	/**< get the point to construct the box: */
	double rp1[3] = {startx, starty, startz},
	                rp2[3] = {startx + xspread / col, starty, startz},
	                         rp3[3] = {startx, starty - yspread / row, startz};
	/**< initilize the edge vertices of the box: */
	(*bbox).rect[0].x[0] = startx - psize;
	(*bbox).rect[0].x[1] = starty + psize;
	(*bbox).rect[1].x[0] = (*bbox).rect[0].x[0] + xspread;
	(*bbox).rect[1].x[1] = (*bbox).rect[0].x[1];
	(*bbox).rect[3].x[0] = (*bbox).rect[0].x[0];
	(*bbox).rect[3].x[1] = (*bbox).rect[0].x[1] - yspread;
	(*bbox).rect[2].x[0] = (*bbox).rect[1].x[0];
	(*bbox).rect[2].x[1] = (*bbox).rect[3].x[1];
	for(i = 0; i < 4; i++)
		memcpy(&(*bbox).rect[i + 4], &(*bbox).rect[i], sizeof(point3));
	(*bbox).rect[0].x[2] = (*bbox).rect[1].x[2] =
	                           (*bbox).rect[2].x[2] = (*bbox).rect[3].x[2] = -psize;
	(*bbox).rect[4].x[2] = (*bbox).rect[5].x[2] =
	                           (*bbox).rect[6].x[2] = (*bbox).rect[7].x[2] = psize;
	snprintf((*bbox).xwidth, INFOLENGTH, "%g mm", fabs((*bbox).rect[0].x[0] - (*bbox).rect[1].x[0]));
	snprintf((*bbox).yheight, INFOLENGTH, "%g mm", fabs((*bbox).rect[0].x[1] - (*bbox).rect[3].x[1]));
	snprintf((*bbox).zdepth, INFOLENGTH, "%g mm", fabs((*bbox).rect[0].x[2] - (*bbox).rect[4].x[2]));
	/**< apply the rotation and translation according to the vector */
	if((*tvec).r[0] == 0. && (*tvec).r[1] == 0. && (*tvec).r[2] == 0.) alpha = beta = 0.;
	else
	{
		alpha = atan2((*tvec).r[0], -(*tvec).r[1]);
		beta = acos((*tvec).r[2]);
	}
	rottrans3_xpz(rp1, beta, alpha, (*tvec).o);
	rottrans3_xpz(rp2, beta, alpha, (*tvec).o);
	rottrans3_xpz(rp3, beta, alpha, (*tvec).o);
	for(i = 0; i < 8; i++) rottrans3_xpz((*bbox).rect[i].x, beta, alpha, (*tvec).o);
	sub3(rp2, rp1, rp2);
	sub3(rp3, rp1, rp3);
	for(i = 0; i < row; i++)
		for(j = 0; j < col; j++)
			if(!mat[i * col + j])
				set_sphrcl_prtcl(&nprtcls[k++], ri, mu_prtcl,
									rp1[0] + j * rp2[0] + i * rp3[0],
									rp1[1] + j * rp2[1] + i * rp3[1],
									rp1[2] + j * rp2[2] + i * rp3[2],
									psize);
	assert(k == np);
	free(mat);
	return nprtcls;
}
