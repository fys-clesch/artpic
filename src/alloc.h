void init_allocations(void);
void print_used_mem(void);
char *alloc_char(const uint len);
uchar *alloc_uchar(const uint len);
cdoub *alloc_cdoub(const uint len);
double *alloc_double(const uint len);
uint *alloc_uint(const uint len);
uint *realloc_uint(uint *m, const uint n, const uint plus);
index2 *alloc_index2(const uint len);
ray *alloc_ray(const uint n);
ray *realloc_ray(ray *m, const double lam, const uint n, const uint plus);
glray_s *alloc_glray_s(const uint n);
glray *alloc_glray(const uint n, const uint rs);
void free_glray(glray *m, const uint n);
gen_ray_info *alloc_gen_ray_info(void);
point3 *alloc_point3(const uint n);
point3 *realloc_point3(point3 *m, const uint n, const uint plus);
line3 *alloc_line3(const uint n);
plane3 *alloc_plane3(const uint n);
sphere3 *alloc_sphere3(const uint n);
intrsec *alloc_intrsec(const uint n);
intrsec *realloc_intrsec(intrsec *m, const uint n, const uint plus);
void free_omatrix(double **m, const uint row);
double **alloc_omatrix(const uint row, const uint col);
double **alloc_omatrix_for_bin_sphere3(const uint pol, const uint azi, const uint field);
void free_ocmatrix(cdoub **m, const uint row);
cdoub **alloc_ocmatrix(const uint row, const uint col);
cdoub **alloc_ocmatrix_for_bin_sphere3(const uint pol, const uint azi, const uint field);
hit_screen *alloc_hit_screen(const uint n);
hit_screen *realloc_hit_screen(hit_screen *m, const uint n, const uint plus);
bin_hit_screen *alloc_bin_hit_screen(const uint n, const double rad, const uint res_polar, const uint res_azim);
void free_bin_hit_screen(bin_hit_screen *m, const uint n);
sphere3 *alloc_nsphere3(const point3 *ps, const double *rs, const uint n);
sphrcl_prtcl *alloc_sphrcl_prtcl(const uint n);
vertex3 *alloc_vertex3(const uint n);
patch3 *alloc_patch3(const uint n);
void free_patch3(patch3 *m, const uint n);
colorval *alloc_colorval(const uint n);
colorbox *alloc_colorbox(const uint n, const uint ncval);
boundingbox *alloc_boundingbox(const uint n);

extern uint allocated_tir_event_memory,
            allocated_exhausted_event_memory,
            allocated_lost_event_memory;
