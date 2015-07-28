uint get_nodd(double x);
uint get_nodi(int x);
uint get_nodui(uint x);
void cp3(double *res_pt dest, const double *res_pt source);
int int_pow(int x, uint n);
double cabs_real(const cdoub z);
double cabs_imag(const cdoub z);
double atan2_up(double x, const double y);
int comp_int(const void *p1, const void *p2);
int comp_double(const void *p1, const void *p2);
double find_max(const double *res_pt m, uint *res_pt xm, uint *res_pt ym, const uint row, const uint col);
double find_min(const double *res_pt m, uint *res_pt xm, uint *res_pt ym, const uint row, const uint col);
uchar solve_pq(const double a, double b, const double c, double r[2], const uchar warn_complex);
double sign(const double x);
void sort_detection(const hit_screen *res_pt hso, const uint no, bin_hit_screen *res_pt bhs, const uchar sumup);
void set_point3(point3 *p, const double x0, const double x1, const double x2);
void set_sphere3(sphere3 *p, const double x0, const double x1, const double x2, const double r);
void set_sphrcl_prtcl(sphrcl_prtcl *p, const cdoub n, const double mu, const double x0, const double x1, const double x2, const double r);
void set_index2_bin_sphere3(index2 *indx, const uint pol, const uint azi);
void copy_ray_to_glray_s(const ray const *res_pt rs, glray_s *res_pt glrs, const uint nors);
void make_trace(const ray const *res_pt r, glray *res_pt glr, const uint n_glr);
int getchar_dscrd_rmng(void);
void wait(uint sec);
void sincosd(const double alpha, double *res_pt cosd, double *res_pt sind) __attribute__((optimize(0)));
void sincos_sqrt(double alpha, double *res_pt cosd, double *res_pt sind);
void init_plane3_pform(const double *res_pt p1, const double *res_pt p2, const double *res_pt p3, plane3 *res_pt r);
void init_plane3(const point3 *res_pt n, const point3 *res_pt o, plane3 *res_pt r);