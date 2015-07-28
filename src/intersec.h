uchar intersec_lin_lin(const line3 *res_pt l1, line3 *res_pt l2, point3 *res_pt p, double *res_pt cangl);
uchar intersec_lin_pla(const plane3 *res_pt e, line3 *res_pt l, point3 *res_pt p);
uchar intersec_lin_tria(const point3 *res_pt tria, const point3 *res_pt pin, line3 *res_pt l, const uchar ret, intrsec *res_pt is);
uchar intersec_lin_sph(const sphere3 *res_pt e, line3 *res_pt l, const uchar ret, intrsec *res_pt is);
