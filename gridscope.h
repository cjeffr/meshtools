/*
 * GRID-SCOPE - C Translation from Fortran
 * 
 * OVER THE LIMIT FOR MEMORY USAGE OF SMS (32bit)
 * Phase 0. Checking Grid
 * Phase 1. Pull out the sub-grid from the global(original) grid.
 * Phase 2. Rezone the nodal attributes on the sub-grid.
 * Phase 3. Recombine edited sub-grid on global grid.
 *
 * Original Fortran: Copyleft 2008- by Seizo Tanaka,
 *                   and C.H.L. at University of Notre Dame
 * C Translation: 2026
 */

#ifndef GRIDSCOPE_H
#define GRIDSCOPE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MAX_STRING_LENGTH 256
#define INITIAL_NEAN_MAX 10
#define INITIAL_NBN_MAX 20
#define EPS_DEFAULT 1.0e-12
#define EPS_SEARCH 1.0e-08
#define OUTFLAG -9999999.0

/* Grid data structures */
typedef struct {
    int ne;                     /* Number of elements */
    int np;                     /* Number of nodes */
    double **xyd;               /* Node coordinates [3][np] - x, y, depth */
    int **nm;                   /* Element connectivity [ne][3] */
} Grid;

/* Boundary data structures */
typedef struct {
    int nope;                   /* Number of elevation boundary segments */
    int neta;                   /* Total elevation boundary nodes */
    int nvdl_max;               /* Max nodes per elevation boundary segment */
    int *nvdll;                 /* Nodes per elevation boundary [nope] */
    int **nbdv;                 /* Elevation boundary nodes [nope][nvdl_max] */
} ElevationBoundary;

typedef struct {
    int nbou;                   /* Number of flow boundary segments */
    int nvel;                   /* Total flow boundary nodes */
    int nvel_max;               /* Max nodes per flow boundary segment */
    int *nvell;                 /* Nodes per flow boundary [nbou] */
    int *ibtype;                /* Boundary type [nbou] */
    int **nbvv;                 /* Flow boundary nodes [nbou][nvel_max] */
    int **ibconn;               /* Connected nodes [nbou][nvel_max] */
    double ***bar;              /* Barrier heights [3][nbou][nvel_max] */
} FlowBoundary;

/* Complete mesh structure */
typedef struct {
    Grid grid;
    ElevationBoundary elev_bc;
    FlowBoundary flow_bc;
    int nodemax;
    int *nsequencer;            /* Node sequencer array */
    char agrid[MAX_STRING_LENGTH];
} Mesh;

/* Nodal attributes (fort.13) structures */
typedef struct {
    int nattr;                  /* Number of attributes */
    int max_valpernode;         /* Maximum values per node */
    char **attrname;            /* Attribute names */
    char **units;               /* Attribute units */
    int *valuespernode;         /* Values per node for each attribute */
    double **defaultattrval;    /* Default values [nattr][max_valpernode] */
} NodalAttributes;

/* Function prototypes - Main phases */
int grid_check(int ifile, const char *fort14, const char *fort13);
int pull_out_subgrid(int ifile, const char *fort14, const char *project_name);
int pull_out_nodals(const char *fort14, const char *project_name, const char *fort13);
int merge_domain(int ifile, const char *fort14, const char *project_name, const char *fort13);

/* Memory allocation functions */
int read14_alloc(FILE *fp, int *ne, int *np, int *nope, int *nbou, 
                 int *nvdl_max, int *nvel_max, int *nodemax);
Mesh* allocate_mesh(int ne, int np, int nope, int nbou, 
                    int nvdl_max, int nvel_max, int nodemax);
void free_mesh(Mesh *mesh);

/* Grid I/O functions */
int read14(FILE *fp, Mesh *mesh);
int write_grid(const char *filename, const Mesh *mesh);

/* Grid selection functions */
int select_circle(int np, double **xyd, int *nprop);
int select_rectan(int np, double **xyd, int *nprop);
int select_arbitr(int np, double **xyd, int *nprop);

/* Grid processing functions */
int mkeline(int np, int ne, int **nm, int *netable, int nbn_max, 
            int *nbl, int **nblnc);
void eline(int i, int *n1, int *n2, int **nm, int m);
int readetab(int ne, int *netable);

/* Pullout functions */
int pullout(int ne, int np, int nope, int nbou, int nvdl_max, int nvel_max,
            char *agrid, double **xyd, int **nm, 
            int *nvdll, int **nbdv, int *nvell, int *ibtype, 
            int **nbvv, int **ibconn, double ***bar,
            const char *project_name);

/* Merge functions */
int merge_(int nope, int nbou, int nvdl_max, int nvel_max,
          char *agrid, int *nvdll, int **nbdv, int *nvell, int *ibtype,
          int **nbvv, int **ibconn, double ***bar,
          int nopes, int nbous, int nvdls_max, int nvels_max,
          int *nvdlls, int **nbdvs, int *nvells, int *ibtypes,
          int **nbvvs, int **ibconns, double ***bars,
          int ifile, const char *project_name, const char *fort14,
          const char *fort13, int nodemax, int *nsequencer,
          Grid *global_grid, Grid *sub_grid);

/* Check functions */
int checkgrd(int ne, int np, int nope, int nbou, int nvdl_max, int nvel_max,
             int nodemax, int *nsequencer, int ifile, const char *fort13,
             double **xyd, int **nm, int *neta, int *nvdll, int **nbdv,
             int *nvel, int *nvell, int *ibtype, int **nbvv, 
             int **ibconn, double ***bar);

int eliminate_node(int ne, int np, int *nen, int *npn, 
                   int nope, int nbou, int nvdl_max, int nvel_max, int *nwork,
                   double **xyd, int **nm, int *neta, int *nvdll, int **nbdv,
                   int *nvel, int *nvell, int *ibtype, int **nbvv,
                   int **ibconn, double ***bar);

int overlap_e(int ne, int np, int *nen, int *npn, double **xyd, int **nm,
              int nean_max, int *numean, int **nean);

void element_area(int np, double **xyd, int n1, int n2, int n3, double *earea);

/* Nodal attributes functions */
int read13_alloc(FILE *fp, int *nattr, int *max_valpernode);
int pullout13(FILE *fp13, int ne, int np, double **xyd, int **nm,
              int nes, int nps, double **xyds, int **nms,
              const char *gridname, int *netable, int *nptable,
              int neos, int *nemap, int npos, int *list_ehasn,
              int nattr, int max_valpernode, int nodemax, int *nsequencer);

/* Search functions */
int count_list(int np, int ne, int **nm, double **xyd, int neos, int *nemap,
               double *xmax, double *xmin, int idiv, int **ne_piece_add,
               int **ne_piece, int *ne_piece_sum);

int mkelocation(int np, int ne, int **nm, double **xyd, int neos, int *nemap,
                int idiv, double *xmax, double *xmin, int ne_in_piece_sum,
                int **ne_piece_add, int **ne_piece, int *ne_piece_list);

int search_inelm(int ne, int np, double **xyd, int **nm, int ilist, int *list_gelem,
                 int *m, double xc, double yc, int idiv, double *xmax, double *xmin,
                 int ne_in_piece_sum, int **ne_piece_add, int **ne_piece,
                 int *ne_piece_list);

void mkareacood(int ne, int np, double **xyd, int **nm, int m, 
                double xc, double yc, double *sl);

/* Utility functions */
double** allocate_2d_double(int rows, int cols);
int** allocate_2d_int(int rows, int cols);
double*** allocate_3d_double(int dim1, int dim2, int dim3);
void free_2d_double(double **arr, int rows);
void free_2d_int(int **arr, int rows);
void free_3d_double(double ***arr, int dim1, int dim2);

/* Global grid data (for merge operations) */
extern Grid *global_grid_data;
extern Grid *sub_grid_data;

#endif /* GRIDSCOPE_H */
