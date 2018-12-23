uchar compare_media(const cdoub n1, const double mu1, const cdoub n2, const double mu2);
void gen_startrays_straight(ray *res_pt rs, const double lam, const uint no, const uint res_azim, const uint res_rad, const double max_rad, const point3 *res_pt i, const point3 *res_pt f, const point3 *res_pt pol);
double get_ray_max_amp(const ray *r);
double get_ray_phased_amp(const ray *r);
double get_ray_int(const ray *r);
void propagate_ray(ray *r);
void propagate_ray_eps(ray *r);
uchar set_refr_index(const ray *res_pt ri, intrsec *res_pt i, const cdoub n_prtcl, const double mu_prtcl);
void reg_hit(const ray *res_pt r, hit_screen *res_pt hs, const double cangl_i, const ray_end_state state);

extern gen_ray_info GLOBAL_RAY_INFO;
