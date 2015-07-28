void norm_rgb(double *rgb);
void constrain_rgb(double *rgb);
void rgbtohsv(double *rgb, double *hsv);
void hsvtorgb(double *hsv, double *rgb);
void color_bin_patch3(patch3 *ptc, const uint nptc, const bin_hit_screen const *bhs, const colorfun cfun, const double alpha, const bin_hit_print_type type);
void map_cfun(double *rgba, double nval, const double min, const double max, const colorfun cfun);

extern const double lighttocolor[401][3],
                    cwhite[4],
                    cblack[4],
                    ctrans[4],
                    cgrey[4],
                    clightgrey[4],
                    cred[4],
                    cgreen[4],
                    cblue[4];
extern const float light_ambient[4],
                   light_diffuse[4],
                   light_specular[4],
                   light_position[4],
                   mat_ambient[4],
                   mat_diffuse[4],
                   mat_specular[4],
                   high_shininess[1];
