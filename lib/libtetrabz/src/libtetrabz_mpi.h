#pragma once
#include "mpi.h"

void libtetrabz_mpi_mp_libtetrabz_mpi_occ_(int *ltetra0, MPI_Comm *comm0, double *bvec, int *nb0, 
                                          int *nge, double *eig, int *ngw, double *wght0);
void libtetrabz_mpi_mp_libtetrabz_mpi_fermieng_(int *ltetra0, MPI_Comm *comm0, double *bvec, int *nb0, 
                                               int *nge, double *eig, int *ngw, double *wght0, 
                                               double *ef, double *nelec);
void libtetrabz_mpi_mp_libtetrabz_mpi_dos_(int *ltetra0, MPI_Comm *comm0, double *bvec, 
                                          int *nb0, int *nge, double *eig, int *ngw, double *wght0, 
                                          int *ne0, double *e0);
void libtetrabz_mpi_mp_libtetrabz_mpi_doubledelta_(int *ltetra0, MPI_Comm *comm0, double *bvec, 
                                                  int *nb0, int *nge, double *eig1, double *eig2, 
                                                  int *ngw, double *wght0);
void libtetrabz_mpi_mp_libtetrabz_mpi_occstep_(int *ltetra0, MPI_Comm *comm0, double *bvec, int *nb0, 
                                              int *nge, double *eig1, double *eig2, int *ngw, double *wght0);
void libtetrabz_mpi_mp_libtetrabz_mpi_polstat_(int *ltetra0, MPI_Comm *comm0, double *bvec, int *nb0, 
                                              int *nge, double *eig1, double *eig2, int *ngw, double *wght0);
void libtetrabz_mpi_mp_libtetrabz_mpi_fermigr_(int *ltetra0, MPI_Comm *comm0, double *bvec, int *nb0, 
                                              int *nge, double *eig1, double *eig2, int *ngw, double *wght0, 
                                              int *ne0, double *e0);
void libtetrabz_mpi_mp_libtetrabz_mpi_polimg_(int *ltetra0, MPI_Comm *comm0, double *bvec, int *nb0, 
                                             int *nge, double *eig1, double *eig2, int *ngw, double *wght0, 
                                             int *ne0, double *e0);
