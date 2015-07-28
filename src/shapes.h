void gen_spherical_table(double *res_pt s_t, double *res_pt c_t, const uint res, const double full_angl);
void gen_sphere3_vertex3(vertex3 *vsph, const uint res_polar, const uint res_azim, const double r);
void fgen_sphere3_vertex3(vertex3 *vsph, const uint res_polar, const uint res_azim, const double r);
patch3 *gen_patch3_bin_sphere3(const uint res_polar, const uint res_azim, const double r, uint *nptc);
patch3 *gen_patch3_sphere3(const uint res_polar, const uint res_azim, const double r, uint *nptc);
void ddraw_sphere3(const uint res_polar, const uint res_azim, const double r);
void sierpinskisponge(uint n_lvls, const double *offset, double scl);
void draw_arrowv(const double *res_pt x1, const double *res_pt x2, const double len, const double width, uint nsegs);
void draw_arrow(const double x1, const double x2, const double x3, const double y1, const double y2, const double y3, const double len, const double width, uint nsegs);
