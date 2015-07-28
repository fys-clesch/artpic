void intf_to_mem(const char *fname, double *m, const uint row, const uint col);
void doubf_to_mem(const char *fname, double *m, const uint row, const uint col);
sphrcl_prtcl *load_prtcls(const char *res_pt fname, cdoub (*reffun)(double, const char *), const double wavlen, const double xspread, const double yspread, const line3 *res_pt tvec, uint *res_pt prtcls, boundingbox *res_pt bbox, const uchar chatty);
