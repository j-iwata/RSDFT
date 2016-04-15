void libtetrabz_mp_libtetrabz_occ_(int *ltetra0, double *bvec, int *nb0, int *nge, 
                                  double *eig, int *ngw, double *wght0);
void libtetrabz_mp_libtetrabz_fermieng_(int *ltetra0, double *bvec, int *nb0, int *nge, 
                                       double *eig, int *ngw, double *wght0, double *ef, double *nelec);
void libtetrabz_mp_libtetrabz_dos_(int *ltetra0, double *bvec, int *nb0, int *nge, double *eig, 
                                  int *ngw, double *wght0, int *ne0, double *e0);
void libtetrabz_mp_libtetrabz_doubledelta_(int *ltetra0, double *bvec, int *nb0, int *nge, 
                                          double *eig1, double *eig2, int *ngw, double *wght0);
void libtetrabz_mp_libtetrabz_occstep_(int *ltetra0, double *bvec, int *nb0, int *nge, double *eig1, 
                                      double *eig2, int *ngw, double *wght0);
void libtetrabz_mp_libtetrabz_polstat_(int *ltetra0, double *bvec, int *nb0, int *nge, double *eig1, 
                                      double *eig2, int *ngw, double *wght0);
void libtetrabz_mp_libtetrabz_fermigr_(int *ltetra0, double *bvec, int *nb0, int *nge, double *eig1, 
                                      double *eig2, int *ngw, double *wght0, int *ne0, double *e0);
void libtetrabz_mp_libtetrabz_polimg_(int *ltetra0, double *bvec, int *nb0, int *nge, double *eig1, 
                                     double *eig2, int *ngw, double *wght0, int *ne0, double *e0);
