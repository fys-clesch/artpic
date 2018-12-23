void drawstring(const char *s);
void drawstring_cb(const double *vmin, const double *vmax, const uint slen);
void drawblockstring2d(const char *s, const uint linewidth, double sx, double sy, const double ix, const double iy);
void draw_coord_ov(void);
void draw_cb(const uchar free_cb);
void handle_glray(const glray *rs, const uint bundles, const draw_opt opt, const uint n_rs, ...);
void handle_sphere3(const uint res_polar, const uint res_azim, const double rad, const double *rgba, const draw_opt opt);
void handle_bin_sphere3(const bin_hit_screen *const bhs, const draw_opt opt, const colourfun cfun, const bin_hit_print_type ptype);
void handle_prtcls_boxed(const sphrcl_prtcl *const res_pt spp, const boundingbox *res_pt bbox, const uint np, const draw_opt opt);
void handle_prtcls(const sphrcl_prtcl *const res_pt spp, const draw_opt opt);
