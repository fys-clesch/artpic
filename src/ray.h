double cangl_in(const ray *res_pt r, const double *res_pt nv);
cdoub ccangl_out(const cdoub ni, const cdoub nt, const double cangl_i, uchar *tir);
void ray_refl(const double *res_pt in, double *res_pt ref, const double *res_pt nv, const double cangl_i);
void ray_trans(const double *res_pt in, double *res_pt trans, const double *res_pt nv, const cdoub ni, const cdoub nt, const double cangl_i, const cdoub cangl_o);
uchar handle_refl(const ray *res_pt ri, ray *res_pt rrefl, const intrsec *res_pt isec, const cdoub ccangl_f);
uchar handle_trans(const ray *res_pt ri, ray *res_pt rtrans, const intrsec *res_pt isec, const cdoub ccangl_f);
uchar handle_reflntrans(ray *res_pt ri, ray *res_pt rsec, const intrsec *res_pt isec);
void map_pol_to_lc(ray *r, const double alpha);
void map_pol_to_lc_opt(ray *r, const double alpha);
uchar orthogo_ray_at_intersec(ray *res_pt ri, const intrsec *res_pt isec);