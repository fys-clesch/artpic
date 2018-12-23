void norm_rgb(double *rgb);
void constrain_rgb(double *rgb);
void rgbtohsv(double *rgb, double *hsv);
void hsvtorgb(double *hsv, double *rgb);
void colour_bin_patch3(patch3 *ptc, const uint nptc, const bin_hit_screen *const bhs, const colourfun cfun, const double alpha, const bin_hit_print_type type);
void map_cfun(double *rgba, double nval, const double min, const double max, const colourfun cfun);

extern const double LIGHTTOCOLOR[401][3],
                    CWHITE[4],
                    CBLACK[4],
                    CTRANS[4],
                    CGREY[4],
                    CLIGHTGREY[4],
                    CRED[4],
                    CGREEN[4],
                    CBLUE[4];
extern const float LIGHT_AMBIENT[4],
                   LIGHT_DIFFUSE[4],
                   LIGHT_SPECULAR[4],
                   LIGHT_POSITION[4],
                   MAT_AMBIENT[4],
                   MAT_DIFFUSE[4],
                   MAT_SPECULAR[4],
                   HIGH_SHININESS[1];
