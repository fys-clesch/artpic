#include "artpic.h"
#include "msg.h"
#include "auxf.h"
#include "refr_data.h"

uint allocated_tir_event_memory, /**< counts the allocated memory for tir. considered for reallocation. */
     allocated_exhausted_event_memory, /**< counts the allocated memory for exhausted rays. considered for reallocation. */
     allocated_lost_event_memory; /**< counts the allocated memory for lost rays. considered for reallocation. */

#if MEMORY_COUNTER
/** \brief sums the allocated memory and stores it in memcounter
 *
 * \param i const uint the allocated memory measured in byte
 * \return void
 *
 */
double memcounter;
#define MEMCOUNTER_FLAG_LEN 16u
char memcounter_flag[MEMCOUNTER_FLAG_LEN];
void count_mem_(const uint i)
{
	static uint j = 0;
	static double k = 1.;
	memcounter += (((double)i) * k);
	if(memcounter > 1e3)
	{
		memcounter /= 1e3;
		k /= 1e3;
		j += 3;
	}
	snprintf(memcounter_flag, MEMCOUNTER_FLAG_LEN, "10^%u byte", j);
}
#define count_mem(x) count_mem_(x)
#else
#define count_mem(x) ((void)0)
#endif

/** \brief initializes the memory counter
 *
 * \param void
 * \return void
 *
 */
void init_allocations(void)
{
#if MEMORY_COUNTER
	memcounter = 0.;
	strncpy(memcounter_flag, "10^0 byte", MEMCOUNTER_FLAG_LEN);
#else
	fprintf(stdout, "\nmemory counter not active\n");
#endif
}
#if MEMORY_COUNTER
#undef MEMCOUNTER_FLAG_LEN
#endif

/** \brief prints the allocated memory to the screen
 *
 * \param void
 * \return void
 *
 */
void print_used_mem(void)
{
#if MEMORY_COUNTER
	fprintf(stdout,
			"\nallocated %.3f * %s (decimal interpretation) by functions in file %s during the execution. "
			"1 byte equals %i bit on the machine this version was built.\n",
			memcounter, memcounter_flag, __FILE__, CHAR_BIT);
#endif
}

/** \brief allocates memory for a char array
 *
 * \param len const uint the length of the array
 * \return char* the array
 *
 */
char *alloc_char(const uint len)
{
	char *m = (char *)calloc(len, sizeof(char));
	count_mem(len * sizeof(char));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	return m;
}

/** \brief allocates memory for an unsigned char array
 *
 * \param len const uint the length of the array
 * \return uchar* the array
 *
 */
uchar *alloc_uchar(const uint len)
{
	uchar *m = (uchar *)calloc(len, sizeof(uchar));
	count_mem(len * sizeof(uchar));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	return m;
}

/** \brief allocates memory for a complex double array
 *
 * \param len const uint the length of the array
 * \return cdoub* the array
 *
 */
cdoub *alloc_cdoub(const uint len)
{
	uint i;
	cdoub *m = (cdoub *)malloc(len * sizeof(cdoub));
	count_mem(len * sizeof(cdoub));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < len; i++) m[i] = 0. + 0.i;
	return m;
}

/** \brief allocates memory for a double array
 *
 * \param len const uint the length of the array
 * \return double* the array
 *
 */
double *alloc_double(const uint len)
{
	uint i;
	double *m = (double *)malloc(len * sizeof(double));
	count_mem(len * sizeof(double));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < len; i++) m[i] = 0.;
	return m;
}

/** \brief allocates memory for an unsigned int array
 *
 * \param len const uint the length of the array
 * \return uint* the array
 *
 */
uint *alloc_uint(const uint len)
{
	uint i;
	uint *m = (uint *)malloc(len * sizeof(uint));
	count_mem(len * sizeof(uint));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < len; i++) m[i] = 0;
	return m;
}

/** \brief reallocates memory for a uint array
 *
 * \param m uint* the existing array
 * \param n const uint the length of the existing array
 * \param plus const uint the additional length of the array
 * \return uint* the array
 *
 */
uint *realloc_uint(uint *m, const uint n, const uint plus)
{
	uint np = n + plus;
	if(np < n)
		error_msg("integer overflow. continuing will likely be erroneous, at least data will be lost", ERR_ARG);
	uint *back = m;
	m = (uint *)realloc(m, np * sizeof(uint));
	count_mem(plus * sizeof(uint));
	if(NULL == m)
	{
		error_msg("memory reallocation failed", ERR_ARG);
		m = back;
	}
	else
	{
		uint i;
		for(i = n; i < np; i++) m[i] = 0;
	}
	return m;
}

/** \brief allocates memory for an index2 array
 *
 * \param len const uint the length of the array
 * \return index2* the array
 *
 */
index2 *alloc_index2(const uint len)
{
	uint i;
	index2 *m = (index2 *)malloc(len * sizeof(index2));
	count_mem(len * sizeof(index2));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < len; i++) m[i] = (index2){.ia = 0, .ib = 0};
	return m;
}

/** \brief allocates memory for a ray array
 *
 * \param n const uint the length of the array
 * \return ray* the array
 *
 */
ray *alloc_ray(const uint n)
{
	uint i;
	static const uint size = sizeof(ray);
	ray *m = (ray *)malloc(n * size);
	static ray r0;
	static uchar init = 0;
	count_mem(n * size);
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	if(!init)
	{
		r0 = (ray){
			.v = (line3) {.o = {0., 0., 0.}, .r = {0., 0., 0.}, .l = 0.},
			.lam = 0., .oamp = 0., .pamp = 0., .ophase = 0., .pphase = 0., .oint = 0., .pint = 0.,
			.opol = {1., 0., 0.}, .ppol = {0., 1., 0.},
			.mu_i = mu_vac, .travel = 0., .n_i = n_vac, .tir = 0, .trans_child = 0, .hits = 0};
		memset(r0.info, 0, INFOLENGTH + 1);
		init = 1;
	}
	for(i = 0; i < n; i++) memcpy(&m[i], &r0, size);
	return m;
}

/** \brief reallocates memory for a ray array
 *
 * \param m ray* the existing array
 * \param lam const double the wavelength of the existing ray array
 * \param n const uint the length of the existing array
 * \param plus const uint the additional length of the array
 * \return ray* the array
 *
 * the new ray variables in the arrays will be initialized with the given wavelegth
 */
ray *realloc_ray(ray *m, const double lam, const uint n, const uint plus)
{
	uint np = n + plus;
	if(np < n)
		error_msg("integer overflow. continuing will likely be erroneous, at least data will be lost", ERR_ARG);
	ray *back = m, r0;
	r0 = (ray){
		.v = (line3) {.o = {0., 0., 0.}, .r = {0., 0., 0.}, .l = 0.},
		.lam = lam, .oamp = 0., .pamp = 0., .ophase = 0., .pphase = 0., .oint = 0., .pint = 0.,
		.opol = {1., 0., 0.}, .ppol = {0., 1., 0.},
		.mu_i = mu_vac, .travel = 0., .n_i = n_vac, .tir = 0, .trans_child = 0, .hits = 0};
	memset(r0.info, 0, INFOLENGTH + 1);
	m = (ray *)realloc(m, np * sizeof(ray));
	count_mem(plus * sizeof(ray));
	if(NULL == m)
	{
		error_msg("memory reallocation failed", ERR_ARG);
		m = back;
	}
	else
	{
		uint i;
		for(i = n; i < np; i++) memcpy(&m[i], &r0, sizeof(ray));
	}
	return m;
}

/** \brief allocates memory for a glray_s array
 *
 * \param n const uint the length of the array
 * \return glray_s* the array
 *
 */
glray_s *alloc_glray_s(const uint n)
{
	uint i;
	glray_s *m = (glray_s *)malloc(n * sizeof(glray_s));
	count_mem(n * sizeof(glray_s));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
		m[i] = (glray_s){
			.v = (line3){.o = {0., 0., 0.}, .r = {0., 0., 0.}, .l = 0.},
			.oamp = 0., .pamp = 0., .ophase = 0., .pphase = 0., .oint = 0., .pint = 0.,
			.opol = {1., 0., 0.}, .ppol = {0., 1., 0.},
			.rgba = {NULL, NULL, NULL, NULL},
			.n_trace = 0,
		    .n_child = 0,
		    .trace_len = 0,
		    .child_len = 0};
	return m;
}

/** \brief allocates memory for a glray array
 *
 * \param n const uint the length of the array
 * \param rs const uint the rays inside of this array
 * \return glray* the array
 *
 * is to be used as a container for multiple rays in OpenGL
 */
glray *alloc_glray(const uint n, const uint rs)
{
	uint i;
	glray *m = (glray *)malloc(n * sizeof(glray));
	count_mem(n * sizeof(glray));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
	{
		m[i].glrs = alloc_glray_s(rs);
		m[i].n_glrs = rs;
	}
	return m;
}

/** \brief frees an allocated glray array
 *
 * \param m glray* the array
 * \param n const uint the number of rays inside the array
 * \return void
 *
 */
void free_glray(glray *m, const uint n)
{
	uint i;
	for(i = 0; i < n; i++)
	{
		uint j;
		for(j = 0; j < m[i].n_glrs; j++)
			if(m[i].glrs[j].n_trace)
			{
				free(m[i].glrs[j].trace);
				free(m[i].glrs[j].child);
			}
		free(m[i].glrs);
	}
	free(m);
}

/** \brief allocates memory for a gen_ray_info variable
 *
 * \param void
 * \return gen_ray_info* the pointer to the variable
 *
 */
gen_ray_info *alloc_gen_ray_info(void)
{
	gen_ray_info *m = (gen_ray_info *)malloc(sizeof(gen_ray_info));
	count_mem(sizeof(gen_ray_info));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	(*m).count_gen =
	(*m).count_hit =
	(*m).count_lost =
	(*m).count_exhstd = 0;
	memset((*m).info, 0, L_INFOLENGTH + 1);
	return m;
}

/** \brief allocates memory for a point3 array
 *
 * \param n const uint the length of the array
 * \return point3* the array
 *
 */
point3 *alloc_point3(const uint n)
{
	uint i;
	point3 *m = (point3 *)malloc(n * sizeof(point3));
	count_mem(n * sizeof(point3));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++) m[i] = (point3){.x = {0., 0., 0.}};
	return m;
}

/** \brief reallocates memory for a point3 array
 *
 * \param m point3* the existing array
 * \param n const uint the length of the existing array
 * \param plus const uint the additional length of the array
 * \return point3* the array
 *
 */
point3 *realloc_point3(point3 *m, const uint n, const uint plus)
{
	uint np = n + plus;
	if(np < n)
		error_msg("integer overflow. continuing will likely be erroneous, at least data will be lost", ERR_ARG);
	point3 *back = m;
	m = (point3 *)realloc(m, np * sizeof(point3));
	count_mem(plus * sizeof(point3));
	if(NULL == m)
	{
		error_msg("memory reallocation failed", ERR_ARG);
		m = back;
	}
	else
	{
		uint i;
		for(i = n; i < np; i++) m[i] = (point3){.x = {0., 0., 0.}};
	}
	return m;
}

/** \brief allocates memory for a line3 array
 *
 * \param n const uint the length of the array
 * \return line3* the array
 *
 */
line3 *alloc_line3(const uint n)
{
	uint i;
	line3 *m = (line3 *)malloc(n * sizeof(line3));
	count_mem(n * sizeof(line3));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
		m[i] = (line3){.o = {0., 0., 0.}, .r = {0., 0., 0.}, .l = 0.};
	return m;
}

/** \brief allocates memory for a plane3 array
 *
 * \param n const uint the length of the array
 * \return plane3* the array
 *
 */
plane3 *alloc_plane3(const uint n)
{
	uint i;
	plane3 *m = (plane3 *)malloc(n * sizeof(plane3));
	count_mem(n * sizeof(plane3));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
		m[i] = (plane3){.o = {0., 0., 0.}, .n = {0., 0., 0.}};
	return m;
}

/** \brief allocates memory for a sphere3 array
 *
 * \param n const uint the length of the array
 * \return sphere3* the array
 *
 */
sphere3 *alloc_sphere3(const uint n)
{
	uint i;
	sphere3 *m = (sphere3 *)malloc(n * sizeof(sphere3));
	count_mem(n * sizeof(sphere3));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++) m[i] = (sphere3){.o = {0., 0., 0.}, .r = 0.};
	return m;
}

/** \brief allocates memory for an intrsec array
 *
 * \param n const uint the length of the array
 * \return intrsec* the array
 *
 */
intrsec *alloc_intrsec(const uint n)
{
	uint i;
	intrsec *m = (intrsec *)malloc(n * sizeof(intrsec));
	count_mem(n * sizeof(intrsec));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
		m[i] = (intrsec){
			.p = {0., 0., 0.}, .normal = {0., 0., 0.},
			.angl = 0., .cangl = 0., .mu_f = mu_vac, .n_f = n_vac, .incdnc = NONE};
	return m;
}

/** \brief reallocates memory for a intrsec array
 *
 * \param m intrsec* the existing array
 * \param n const uint the length of the existing array
 * \param plus const uint the additional length of the array
 * \return intrsec* the array
 *
 */
intrsec *realloc_intrsec(intrsec *m, const uint n, const uint plus)
{
	uint np = n + plus;
	if(np < n)
		error_msg("integer overflow. continuing will likely be erroneous, at least data will be lost", ERR_ARG);
	intrsec *back = m;
	m = (intrsec *)realloc(m, np * sizeof(intrsec));
	count_mem(plus * sizeof(intrsec));
	if(NULL == m)
	{
		error_msg("memory reallocation failed", ERR_ARG);
		m = back;
	}
	else
	{
		uint i;
		for(i = n; i < np; i++)
			m[i] = (intrsec){
				.p = {0., 0., 0.}, .normal = {0., 0., 0.},
				.angl = 0., .cangl = 0., .mu_f = mu_vac, .n_f = n_vac, .incdnc = NONE};
	}
	return m;
}

/** \brief frees an allocated double matrix
 *
 * \param m double** the matrix
 * \param row const uint the range of the first index
 * \return void
 *
 */
void free_omatrix(double **m, const uint row)
{
	uint i;
	for(i = 0; i < row; i++) free(m[i]);
	free(m);
}

/** \brief allocates memory for a double matrix
 *
 * \param row const uint the length of the first index
 * \param col const uint the length of the second index
 * \return double** the matrix
 *
 */
double **alloc_omatrix(const uint row, const uint col)
{
	double **m;
	uint i, j;
	m = (double **)malloc(row * sizeof(double *));
	count_mem(row * sizeof(double *));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < row; i++)
	{
		m[i] = (double *)malloc(col * sizeof(double));
		count_mem(col * sizeof(double));
		if(NULL == m[i])
		{
			error_msg("memory allocation failed. good bye.", ERR_ARG);
			free_omatrix(m, i);
			exit(EXIT_FAILURE);
		}
		for(j = 0; j < col; j++) m[i][j] = 0.;
	}
	return m;
}

/** \brief allocates memory for a double matrix used for a bin_sphere3
 *
 * \param pol const uint the polar slices of the bin_sphere3 variable
 * \param azi const uint the azimuthal slices of the bin_sphere3 variable
 * \param field const uint the number of azimuthal elements in the 'most polar' slices
 * \return double** the matrix
 *
 * same functionallity as alloc_omatrix, but specially designed for the usage in alloc_bin_hit_screen
 */
double **alloc_omatrix_for_bin_sphere3(const uint pol, const uint azi, const uint field)
{
	double **m;
	uint i, j, t1 = field * azi;
	if(t1 < azi) error_msg("integer overflow", ERR_ARG);
	m = (double **)malloc(pol * sizeof(double *));
	count_mem(pol * sizeof(double *));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < pol; i++)
	{
		uint k;
		(i == 0 || i == pol - 1) ? (k = field) : (k = t1);
		m[i] = (double *)malloc(k * sizeof(double));
		count_mem(k * sizeof(double));
		if(NULL == m[i])
		{
			error_msg("memory allocation failed. good bye.", ERR_ARG);
			free_omatrix(m, i);
			exit(EXIT_FAILURE);
		}
		for(j = 0; j < k; j++) m[i][j] = 0.;
	}
	return m;
}

/** \brief frees an allocated complex double matrix
 *
 * \param m cdoub** the matrix
 * \param row const uint the range of the first index
 * \return void
 *
 */
void free_ocmatrix(cdoub **m, const uint row)
{
	uint i;
	for(i = 0; i < row; i++) free(m[i]);
	free(m);
}

/** \brief allocates memory for a double matrix
 *
 * \param row const uint the length of the first index
 * \param col const uint the length of the second index
 * \return cdoub** the matrix
 *
 */
cdoub **alloc_ocmatrix(const uint row, const uint col)
{
	cdoub **m;
	uint i, j;
	m = (cdoub **)malloc(row * sizeof(cdoub *));
	count_mem(row * sizeof(cdoub *));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < row; i++)
	{
		m[i] = (cdoub *)malloc(col * sizeof(cdoub));
		count_mem(col * sizeof(cdoub));
		if(NULL == m[i])
		{
			error_msg("memory allocation failed. good bye.", ERR_ARG);
			free_ocmatrix(m, i);
			exit(EXIT_FAILURE);
		}
		for(j = 0; j < col; j++) m[i][j] = 0. + 0.i;
	}
	return m;
}

/** \brief allocates memory for a complex double matrix used for a bin_sphere3
 *
 * \param pol const uint the polar slices of the bin_sphere3 variable
 * \param azi const uint the azimuthal slices of the bin_sphere3 variable
 * \param field const uint the number of azimuthal elements in the 'most polar' slices
 * \return cdoub** the matrix
 *
 * same functionallity as alloc_omatrix, but specially designed for the usage in alloc_bin_hit_screen
 */
cdoub **alloc_ocmatrix_for_bin_sphere3(const uint pol, const uint azi, const uint field)
{
	cdoub **m;
	uint i, j, t1 = field * azi;
	if(t1 < azi) error_msg("integer overflow", ERR_ARG);
	m = (cdoub **)malloc(pol * sizeof(cdoub *));
	count_mem(pol * sizeof(cdoub *));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < pol; i++)
	{
		uint k;
		(i == 0 || i == pol - 1) ? (k = field) : (k = t1);
		m[i] = (cdoub *)malloc(k * sizeof(cdoub));
		count_mem(k * sizeof(cdoub));
		if(NULL == m[i])
		{
			error_msg("memory allocation failed. good bye.", ERR_ARG);
			free_ocmatrix(m, i);
			exit(EXIT_FAILURE);
		}
		for(j = 0; j < k; j++) m[i][j] = 0.;
	}
	return m;
}

/** \brief allocates memory for a hit_screen array
 *
 * \param n const uint the length of the array
 * \return hit_screen* the array
 *
 */
hit_screen *alloc_hit_screen(const uint n)
{
	uint i;
	static const uint size = sizeof(hit_screen);
	hit_screen *m = (hit_screen *)malloc(n * size);
	static hit_screen h0;
	static uchar init = 0;
	count_mem(n * size);
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	if(!init)
	{
		h0 = (hit_screen){
			.p = {0., 0., 0.}, .lam = 0., .cos_incdnc = 0xDEAD, .opol = {0., 0., 0.}, .ppol = {0., 0., 0.},
			.coamp = 0. + 0.i, .cpamp = 0. + 0.i, .oint = 0., .pint = 0., .tir = 0, .state = EMPTY_STATE};
		memset(h0.info, 0, INFOLENGTH + 1);
		init = 1;
	}
	for(i = 0; i < n; i++) memcpy(&m[i], &h0, size);
	return m;
}

/** \brief reallocates memory for a hit_screen array
 *
 * \param m hit_screen* the existing array
 * \param n const uint the length of the existing array
 * \param plus const uint the additional length of the array
 * \return hit_screen* the array
 *
 */
hit_screen *realloc_hit_screen(hit_screen *m, const uint n, const uint plus)
{
	uint np = n + plus;
	if(np < n)
		error_msg("integer overflow. continuing will likely be erroneous, at least data will be lost", ERR_ARG);
	hit_screen *back = m, h0 = (hit_screen){
		.p = {0., 0., 0.}, .lam = 0., .cos_incdnc = 0xDEAD, .opol = {0., 0., 0.}, .ppol = {0., 0., 0.},
		.coamp = 0. + 0.i, .cpamp = 0. + 0.i, .oint = 0., .pint = 0., .tir = 0, .state = EMPTY_STATE};
	memset(h0.info, 0, INFOLENGTH + 1);
	m = (hit_screen *)realloc(m, np * sizeof(hit_screen));
	count_mem(plus * sizeof(hit_screen));
	if(NULL == m)
	{
		error_msg("memory reallocation failed", ERR_ARG);
		m = back;
	}
	else
	{
		uint i;
		for(i = n; i < np; i++) memcpy(&m[i], &h0, sizeof(hit_screen));
	}
	return m;
}

/** \brief allocates memory for a bin_hit_screen array
 *
 * \param n const uint the length of the array
 * \param rad const double the radius of the sphere
 * \param res_polar const uint the polar resolution of the sphere
 * \param res_azim const uint the azimuthal resolution of the sphere
 * \return bin_hit_screen* the array
 *
 * also sets the initial value for
 * - allocated_tir_event_memory
 * - allocated_lost_event_memory
 * - allocated_exhausted_event_memory
 */
bin_hit_screen *alloc_bin_hit_screen(const uint n, const double rad, const uint res_polar, const uint res_azim)
{
	if(res_polar < 3)
	{
		error_msg("the polar resolution is smaller than 3. restart with higher resolution.", ERR_ARG);
		exit(EXIT_FAILURE);
	}
	uint i;
	bin_hit_screen *m = (bin_hit_screen *)malloc(n * sizeof(bin_hit_screen));
	count_mem(n * sizeof(bin_hit_screen));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	allocated_tir_event_memory =
	allocated_lost_event_memory =
	allocated_exhausted_event_memory = 1000;
	for(i = 0; i < n; i++)
	{
		m[i] = (bin_hit_screen){
			.lam = 0., .rad = rad,
			.amp_max = 0., .amp_sum = 0., .mod_amp_max = 0., .mod_amp_sum = 0.,
			.ray_int_max = 0., .ray_int_sum = 0., .pol_max = 0., .pol_sum = 0.,
			.res_polar = res_polar, .res_azim = res_azim, .tir = 0,
			.screen_hit_coor_sys = CARTESIAN_CS,
			.global_info = (gen_ray_info){
				.count_gen = 0, .count_hit = 0, .count_lost = 0, .count_exhstd = 0}};
		m[i].amp = alloc_omatrix_for_bin_sphere3(res_polar, res_azim, 1);
		m[i].mod_amp = alloc_omatrix_for_bin_sphere3(res_polar, res_azim, 1);
		m[i].ray_int = alloc_omatrix_for_bin_sphere3(res_polar, res_azim, 1);
		m[i].pol_dist = alloc_omatrix_for_bin_sphere3(res_polar, res_azim, 3); /**< at each coordinate there are 3 entries to account for the polarization in x, y and z */
		m[i].pol_dens = alloc_omatrix_for_bin_sphere3(res_polar, res_azim, 1);
		m[i].tir_coor = alloc_point3(allocated_tir_event_memory);
		m[i].exhstd_coor = alloc_point3(allocated_exhausted_event_memory);
		m[i].lost_coor = alloc_point3(allocated_lost_event_memory);
		m[i].camp = alloc_ocmatrix_for_bin_sphere3(res_polar, res_azim, 3);
		m[i].nbin = (res_polar - 2) * res_azim + 2; /**< the poles do not have an azimutal spacing */
		m[i].idx = alloc_index2(m[i].nbin);
		set_index2_bin_sphere3(m[i].idx, m[i].res_polar, m[i].res_azim);
		memset(m[i].global_info.info, 0, L_INFOLENGTH + 1);
	}
	return m;
}

/** \brief frees an allocated bin_hit_screen array
 *
 * \param m bin_hit_screen* the array
 * \param n const uint the length of the array
 * \return void
 *
 */
void free_bin_hit_screen(bin_hit_screen *m, const uint n)
{
	uint i, j = m[0].res_polar;
	for(i = 0; i < n; i++)
	{
		free_omatrix(m[i].amp, j);
		free_omatrix(m[i].mod_amp, j);
		free_omatrix(m[i].ray_int, j);
		free_omatrix(m[i].pol_dist, j);
		free_omatrix(m[i].pol_dens, j);
		free(m[i].tir_coor);
		free(m[i].exhstd_coor);
		free(m[i].lost_coor);
		free_ocmatrix(m[i].camp, j);
		free(m[i].idx);
	}
	free(m);
}

/** \brief allocates memory and initializes a sphere3 array
 *
 * \param ss sphere3* a sphere3 pointer
 * \param ps const point3* a point3 array to set the center of each sphere3
 * \param rs const double* a double array to set the radius of each sphere3
 * \param n const uint the length of the array
 * \return void
 *
 */
void alloc_nsphere3(sphere3 *ss, const point3 *ps, const double *rs, const uint n)
{
	uint i;
	ss = alloc_sphere3(n);
	for(i = 0; i < n; i++)
	{
		cp3(ss[i].o, ps[i].x);
		ss[i].r = rs[i];
	}
}

/** \brief allocates memory for a sphrcl_prtcl array
 *
 * \param n const uint the length of the array
 * \return sphrcl_prtcl* the array
 *
 */
sphrcl_prtcl *alloc_sphrcl_prtcl(const uint n)
{
	uint i;
	static uint j = 0;
	sphrcl_prtcl *m = (sphrcl_prtcl *)malloc(n * sizeof(sphrcl_prtcl));
	count_mem(n * sizeof(sphrcl_prtcl));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
	{
		m[i] = (sphrcl_prtcl){
			.s = (sphere3){.o = {0., 0., 0.}, .r = 0.}, .n = n_vac, .mu = mu_vac, .hits = 0, .exhstds = 0, .no = j++};
		memset(m[i].info, 0, L_INFOLENGTH + 1);
	}
	return m;
}

/** \brief allocates memory for a vertex3 array
 *
 * \param n const uint the length of the array
 * \return vertex3* the array
 *
 */
vertex3 *alloc_vertex3(const uint n)
{
	uint i;
	vertex3 *m = (vertex3 *)malloc(n * sizeof(vertex3));
	count_mem(n * sizeof(vertex3));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
		m[i] = (vertex3){.n = {0., 0., 0.}, .x = {0., 0., 0.}};
	return m;
}

/** \brief allocates memory for a patch3 array
 *
 * \param n const uint the length of the array
 * \return patch3* the array
 *
 * mind that the member .vt still has to be allocated
 */
patch3 *alloc_patch3(const uint n)
{
	uint i;
	patch3 *m = (patch3 *)malloc(n * sizeof(patch3));
	count_mem(n * sizeof(patch3));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
		m[i] = (patch3){.n_vt = 0, .gl_primitive = 0, .rgba = {0., 0., 0., 0.}};
	return m;
}

/** \brief frees an allocated patch3 array
 *
 * \param m patch3* the array
 * \param n const uint the length of the array
 * \return void
 *
 */
void free_patch3(patch3 *m, const uint n)
{
	uint i;
	for(i = 0; i < n; i++) free(m[i].vt);
	free(m);
}

/** \brief allocates memory for a colorval array
 *
 * \param n const uint the length of the array
 * \return colorval* the array
 *
 */
colorval *alloc_colorval(const uint n)
{
	uint i;
	colorval *m = (colorval *)malloc(n * sizeof(colorval));
	count_mem(n * sizeof(colorval));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
		m[i] = (colorval){.rgba = {0xDEAD, 0xDEAD, 0xDEAD, 0xDEAD}, .val = 0xDEAD};
	return m;
}

/** \brief allocates memory for a colorbox array
 *
 * \param n const uint the length of the array
 * \param ncval const uint the number of colorvals inside the colorbox
 * \return colorbox* the array
 *
 */
colorbox *alloc_colorbox(const uint n, const uint ncval)
{
	uint i;
	colorbox *m = (colorbox *)malloc(n * sizeof(colorbox));
	count_mem(n * sizeof(colorbox));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
	{
		m[i] = (colorbox){.max = 0., .min = 0., .ncval = ncval};
		m[i].cval = alloc_colorval(ncval);
	}
	return m;
}

/** \brief allocates memory for a boundingbox array
 *
 * \param n const uint the length of the array
 * \return boundingbox* the array
 *
 */
boundingbox *alloc_boundingbox(const uint n)
{
	uint i, j;
	const uint k = 8;
	boundingbox *m = (boundingbox *)malloc(n * sizeof(boundingbox));
	count_mem(n * sizeof(boundingbox));
	if(NULL == m)
	{
		error_msg("memory allocation failed. good bye.", ERR_ARG);
		free(m);
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < k; j++)
			m[i].rect[j] = (point3){.x = {0., 0., 0.}};
		memset(m[i].xwidth, 0, INFOLENGTH + 1);
		memset(m[i].yheight, 0, INFOLENGTH + 1);
		memset(m[i].zdepth, 0, INFOLENGTH + 1);
	}
	return m;
}
