#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libtetrabz.h"
/**/
/*use libtetrabz, only : libtetrabz_fermieng, libtetrabz_dos*/
/**/
int main(int argc, char *argv[])
{
/**/
  int ltetra, nb, ng, nge[3], ngw[3], ii, i0, i1, i2, ik, ne, ie, nke, nkw;
  double bvec[9], ef, nelec, kvec[3], pi, sumdos;
  double *eig, *wght, *e0, *wght_dos;
  FILE *fp;
  /**/
  printf("Which tetrahedron method ?(1 = Linear, 2 = Optimized): ");
  scanf("%d", &ltetra);
  printf("\n k-point mesh ?: ");
  scanf("%d", &ng);
  /**/
  pi = acos(-1.0);
  nb = 1;
  for(ii = 0; ii < 3; ii++) nge[ii] = ng;
  for(ii = 0; ii < 3; ii++) ngw[ii] = ng;
  nelec = 0.5;
  bvec[0] = 1.0; bvec[1] = 0.0; bvec[2] = 0.0;
  bvec[3] = 0.0; bvec[4] = 1.0; bvec[5] = 0.0;
  bvec[6] = 0.0; bvec[7] = 0.0; bvec[8] = 1.0;
  nke = nge[0] * nge[1] * nge[2];
  nkw = ngw[0] * ngw[1] * ngw[2];
  /**/
  eig = (double*)malloc(nb * nke * sizeof(double));
  wght = (double*)malloc(nb * nkw * sizeof(double));
  /**/
  ik = 0;
  for(i2 = 0; i2 < nge[2]; i2++){
    for(i1 = 0; i1 < nge[1]; i1++){
      for(i0 = 0; i0 < nge[0]; i0++){
        /**/
        kvec[0] = 2.0 * pi * (double)(i0 - nge[0] / 2) / (double)nge[0];
        kvec[1] = 2.0 * pi * (double)(i1 - nge[1] / 2) / (double)nge[1];
        kvec[2] = 2.0 * pi * (double)(i2 - nge[2] / 2) / (double)nge[2];
        /**/
        eig[nb * ik] = - cos(kvec[0]) - cos(kvec[1]) - cos(kvec[2]);
        /**/
        ik += 1;
      }
    }
  }
  /**/
  libtetrabz_mp_libtetrabz_fermieng_(&ltetra,bvec,&nb,nge,eig,ngw,wght,&ef,&nelec);
  /**/
  printf("  E_F = %25.15e\n", ef);
  /**/
  ne = 100;
  wght_dos = (double*)malloc(ne * nb * nkw * sizeof(double));
  e0 = (double*)malloc(ne * sizeof(double));
  /**/
  for(ie = 0; ie < ne; ie++){
    e0[ie] = 6.0 / (double)(ne - 1) * (double)ie - 3.0;
  }
  /**/
  libtetrabz_mp_libtetrabz_dos_(&ltetra,bvec,&nb,nge,eig,ngw,wght_dos,&ne,e0);
  /**/
  fp = fopen("dos.dat", "w");
  for(ie = 0; ie < ne; ie++){
    /**/
    sumdos = 0.0;
    for(ik = 0; ik < nb * nkw; ik++) sumdos += wght_dos[ie + ik * ne];
    /**/
    fprintf(fp,"%25.15e  %25.15e\n", e0[ie], sumdos);
  }
  /**/
  fclose(fp);
  /**/
}
