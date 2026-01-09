/*
 * gridscope.c - Core implementation
 * C translation of gridscope.f90
 */

#include "gridscope.h"

/* Global variables for module data (equivalent to Fortran modules) */
Grid *global_grid_data = NULL;
Grid *sub_grid_data = NULL;

/* ===================================================================
 * UTILITY FUNCTIONS - Memory Allocation
 * =================================================================== */

double** allocate_2d_double(int rows, int cols) {
    double **arr = (double **)malloc(rows * sizeof(double *));
    if (!arr) return NULL;
    
    for (int i = 0; i < rows; i++) {
        arr[i] = (double *)calloc(cols, sizeof(double));
        if (!arr[i]) {
            for (int j = 0; j < i; j++) free(arr[j]);
            free(arr);
            return NULL;
        }
    }
    return arr;
}

int** allocate_2d_int(int rows, int cols) {
    int **arr = (int **)malloc(rows * sizeof(int *));
    if (!arr) return NULL;
    
    for (int i = 0; i < rows; i++) {
        arr[i] = (int *)calloc(cols, sizeof(int));
        if (!arr[i]) {
            for (int j = 0; j < i; j++) free(arr[j]);
            free(arr);
            return NULL;
        }
    }
    return arr;
}

double*** allocate_3d_double(int dim1, int dim2, int dim3) {
    double ***arr = (double ***)malloc(dim1 * sizeof(double **));
    if (!arr) return NULL;
    
    for (int i = 0; i < dim1; i++) {
        arr[i] = allocate_2d_double(dim2, dim3);
        if (!arr[i]) {
            for (int j = 0; j < i; j++) free_2d_double(arr[j], dim2);
            free(arr);
            return NULL;
        }
    }
    return arr;
}

void free_2d_double(double **arr, int rows) {
    if (!arr) return;
    for (int i = 0; i < rows; i++) {
        if (arr[i]) free(arr[i]);
    }
    free(arr);
}

void free_2d_int(int **arr, int rows) {
    if (!arr) return;
    for (int i = 0; i < rows; i++) {
        if (arr[i]) free(arr[i]);
    }
    free(arr);
}

void free_3d_double(double ***arr, int dim1, int dim2) {
    if (!arr) return;
    for (int i = 0; i < dim1; i++) {
        if (arr[i]) free_2d_double(arr[i], dim2);
    }
    free(arr);
}

Mesh* allocate_mesh(int ne, int np, int nope, int nbou, 
                    int nvdl_max, int nvel_max, int nodemax) {
    Mesh *mesh = (Mesh *)calloc(1, sizeof(Mesh));
    if (!mesh) return NULL;
    
    /* Allocate grid data */
    mesh->grid.ne = ne;
    mesh->grid.np = np;
    mesh->grid.xyd = allocate_2d_double(3, np);
    mesh->grid.nm = allocate_2d_int(ne, 3);
    
    /* Allocate elevation boundaries */
    mesh->elev_bc.nope = nope;
    mesh->elev_bc.nvdl_max = nvdl_max;
    mesh->elev_bc.nvdll = (int *)calloc(nope, sizeof(int));
    mesh->elev_bc.nbdv = allocate_2d_int(nope, nvdl_max);
    
    /* Allocate flow boundaries */
    mesh->flow_bc.nbou = nbou;
    mesh->flow_bc.nvel_max = nvel_max;
    mesh->flow_bc.nvell = (int *)calloc(nbou, sizeof(int));
    mesh->flow_bc.ibtype = (int *)calloc(nbou, sizeof(int));
    mesh->flow_bc.nbvv = allocate_2d_int(nbou, nvel_max);
    mesh->flow_bc.ibconn = allocate_2d_int(nbou, nvel_max);
    mesh->flow_bc.bar = allocate_3d_double(3, nbou, nvel_max);
    
    /* Allocate node sequencer */
    mesh->nodemax = nodemax;
    mesh->nsequencer = (int *)calloc(nodemax, sizeof(int));
    
    if (!mesh->grid.xyd || !mesh->grid.nm || !mesh->elev_bc.nvdll || 
        !mesh->elev_bc.nbdv || !mesh->flow_bc.nvell || !mesh->flow_bc.ibtype ||
        !mesh->flow_bc.nbvv || !mesh->flow_bc.ibconn || !mesh->flow_bc.bar ||
        !mesh->nsequencer) {
        free_mesh(mesh);
        return NULL;
    }
    
    return mesh;
}

void free_mesh(Mesh *mesh) {
    if (!mesh) return;
    
    free_2d_double(mesh->grid.xyd, 3);
    free_2d_int(mesh->grid.nm, mesh->grid.ne);
    
    if (mesh->elev_bc.nvdll) free(mesh->elev_bc.nvdll);
    free_2d_int(mesh->elev_bc.nbdv, mesh->elev_bc.nope);
    
    if (mesh->flow_bc.nvell) free(mesh->flow_bc.nvell);
    if (mesh->flow_bc.ibtype) free(mesh->flow_bc.ibtype);
    free_2d_int(mesh->flow_bc.nbvv, mesh->flow_bc.nbou);
    free_2d_int(mesh->flow_bc.ibconn, mesh->flow_bc.nbou);
    free_3d_double(mesh->flow_bc.bar, 3, mesh->flow_bc.nbou);
    
    if (mesh->nsequencer) free(mesh->nsequencer);
    
    free(mesh);
}

/* ===================================================================
 * FILE I/O FUNCTIONS
 * =================================================================== */

int read14_alloc(FILE *fp, int *ne, int *np, int *nope, int *nbou, 
                 int *nvdl_max, int *nvel_max, int *nodemax) {
    char line[MAX_STRING_LENGTH];
    int i, j, k, n, nhy;
    
    *nvdl_max = 0;
    *nvel_max = 0;
    *nodemax = 0;
    
    /* Read header */
    if (!fgets(line, sizeof(line), fp)) return -1;  /* agrid */
    if (fscanf(fp, "%d %d", ne, np) != 2) return -1;
    
    /* Read nodes to find nodemax */
    for (k = 0; k < *np; k++) {
        if (fscanf(fp, "%d", &i) != 1) return -1;
        if (i > *nodemax) *nodemax = i;
        /* Skip coordinates */
        if (!fgets(line, sizeof(line), fp)) return -1;
    }
    
    printf("  |\n");
    
    /* Read elements */
    for (k = 0; k < *ne; k++) {
        if (!fgets(line, sizeof(line), fp)) return -1;
    }
    
    /* Read elevation boundaries */
    if (fscanf(fp, "%d", nope) != 1) return -1;
    if (!fgets(line, sizeof(line), fp)) return -1;  /* Skip rest of line */
    if (!fgets(line, sizeof(line), fp)) return -1;  /* Total elevation boundary nodes */
    
    for (k = 0; k < *nope; k++) {
        if (fscanf(fp, "%d", &i) != 1) return -1;
        if (i >= *nvdl_max) *nvdl_max = i;
        for (j = 0; j < i; j++) {
            if (!fgets(line, sizeof(line), fp)) return -1;
        }
    }
    
    /* Read flow boundaries */
    if (fscanf(fp, "%d", nbou) != 1) return -1;
    if (!fgets(line, sizeof(line), fp)) return -1;
    if (!fgets(line, sizeof(line), fp)) return -1;
    
    for (k = 0; k < *nbou; k++) {
        if (fscanf(fp, "%d", &i) != 1) return -1;
        if (i >= *nvel_max) *nvel_max = i;
        for (j = 0; j < i; j++) {
            if (!fgets(line, sizeof(line), fp)) return -1;
        }
    }
    
    return 0;
}

int read14(FILE *fp, Mesh *mesh) {
    char line[MAX_STRING_LENGTH];
    int i, j, k, jn, je, nhy;
    int ne = mesh->grid.ne;
    int np = mesh->grid.np;
    int nope = mesh->elev_bc.nope;
    int nbou = mesh->flow_bc.nbou;
    int nvdl_max = mesh->elev_bc.nvdl_max;
    int nvel_max = mesh->flow_bc.nvel_max;
    
    /* Initialize */
    memset(mesh->nsequencer, 0, mesh->nodemax * sizeof(int));
    
    /* Read agrid */
    if (!fgets(mesh->agrid, sizeof(mesh->agrid), fp)) return -1;
    
    /* Read ne, np */
    if (fscanf(fp, "%d %d", &ne, &np) != 2) return -1;
    
    /* Read nodes */
    for (k = 0; k < np; k++) {
        if (fscanf(fp, "%d %lf %lf %lf", &jn,
                   &mesh->grid.xyd[0][k],
                   &mesh->grid.xyd[1][k],
                   &mesh->grid.xyd[2][k]) != 4) return -1;
        mesh->nsequencer[jn - 1] = k;  /* Fortran is 1-indexed */
    }
    
    printf("  + \n");
    
    /* Read elements and convert to 0-based node indices */
    for (k = 0; k < ne; k++) {
        if (fscanf(fp, "%d %d %d %d %d", &je, &nhy,
                   &mesh->grid.nm[k][0],
                   &mesh->grid.nm[k][1],
                   &mesh->grid.nm[k][2]) != 5) return -1;
        
        for (j = 0; j < 3; j++) {
            if (mesh->grid.nm[k][j] <= 0) {
                printf("%d %d %d\n", k, j, mesh->grid.nm[k][j]);
            }
            mesh->grid.nm[k][j] = mesh->nsequencer[mesh->grid.nm[k][j] - 1];
        }
    }
    
    /* Read elevation boundaries */
    if (fscanf(fp, "%d", &mesh->elev_bc.nope) != 1) return -1;
    if (fscanf(fp, "%d", &mesh->elev_bc.neta) != 1) return -1;
    
    for (k = 0; k < nope; k++) {
        if (fscanf(fp, "%d", &mesh->elev_bc.nvdll[k]) != 1) return -1;
        for (j = 0; j < mesh->elev_bc.nvdll[k]; j++) {
            if (fscanf(fp, "%d", &mesh->elev_bc.nbdv[k][j]) != 1) return -1;
            mesh->elev_bc.nbdv[k][j] = mesh->nsequencer[mesh->elev_bc.nbdv[k][j] - 1];
        }
    }
    
    /* Read flow boundaries */
    if (fscanf(fp, "%d", &mesh->flow_bc.nbou) != 1) return -1;
    if (fscanf(fp, "%d", &mesh->flow_bc.nvel) != 1) return -1;
    
    for (k = 0; k < nbou; k++) {
        if (fscanf(fp, "%d %d", &mesh->flow_bc.nvell[k], 
                   &mesh->flow_bc.ibtype[k]) != 2) return -1;
        
        int ibtype = mesh->flow_bc.ibtype[k];
        
        /* Handle different boundary types */
        if (ibtype == 0 || ibtype == 1 || ibtype == 2 || ibtype == 10 ||
            ibtype == 11 || ibtype == 12 || ibtype == 20 || ibtype == 21 ||
            ibtype == 22 || ibtype == 30 || ibtype == 52) {
            /* Simple boundary */
            for (j = 0; j < mesh->flow_bc.nvell[k]; j++) {
                if (fscanf(fp, "%d", &mesh->flow_bc.nbvv[k][j]) != 1) return -1;
                mesh->flow_bc.nbvv[k][j] = mesh->nsequencer[mesh->flow_bc.nbvv[k][j] - 1];
            }
        } else if (ibtype == 3 || ibtype == 13 || ibtype == 23) {
            /* Barrier boundary */
            for (j = 0; j < mesh->flow_bc.nvell[k]; j++) {
                if (fscanf(fp, "%d %lf %lf", &mesh->flow_bc.nbvv[k][j],
                           &mesh->flow_bc.bar[0][k][j],
                           &mesh->flow_bc.bar[1][k][j]) != 3) return -1;
                mesh->flow_bc.nbvv[k][j] = mesh->nsequencer[mesh->flow_bc.nbvv[k][j] - 1];
            }
        } else if (ibtype == 4 || ibtype == 24) {
            /* Connected boundary */
            for (j = 0; j < mesh->flow_bc.nvell[k]; j++) {
                if (fscanf(fp, "%d %d %lf %lf %lf", 
                           &mesh->flow_bc.nbvv[k][j],
                           &mesh->flow_bc.ibconn[k][j],
                           &mesh->flow_bc.bar[0][k][j],
                           &mesh->flow_bc.bar[1][k][j],
                           &mesh->flow_bc.bar[2][k][j]) != 5) return -1;
                mesh->flow_bc.nbvv[k][j] = mesh->nsequencer[mesh->flow_bc.nbvv[k][j] - 1];
                mesh->flow_bc.ibconn[k][j] = mesh->nsequencer[mesh->flow_bc.ibconn[k][j] - 1];
            }
        } else {
            printf(" IBTYPE is not allowed %d\n", ibtype);
            printf("\n**** Hit the Enter-Key to stop ****\n");
            getchar();
            return -1;
        }
    }
    
    return 0;
}

/* ===================================================================
 * GRID PROCESSING FUNCTIONS
 * =================================================================== */

void eline(int i, int *n1, int *n2, int **nm, int m) {
    switch (i) {
        case 2:  /* case 3 in Fortran (1-indexed) */
            *n1 = nm[m][0];
            *n2 = nm[m][1];
            break;
        case 0:  /* case 1 in Fortran */
            *n1 = nm[m][1];
            *n2 = nm[m][2];
            break;
        case 1:  /* case 2 in Fortran */
            *n1 = nm[m][2];
            *n2 = nm[m][0];
            break;
    }
}

void element_area(int np, double **xyd, int n1, int n2, int n3, double *earea) {
    double x1 = xyd[0][n1];
    double x2 = xyd[0][n2];
    double x3 = xyd[0][n3];
    double y1 = xyd[1][n1];
    double y2 = xyd[1][n2];
    double y3 = xyd[1][n3];
    
    *earea = (x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2);
}

int mkeline(int np, int ne, int **nm, int *netable, int nbn_max, 
            int *nbl, int **nblnc) {
    int *numean = (int *)calloc(np, sizeof(int));
    int **nean = NULL;
    int nean_max = INITIAL_NEAN_MAX;
    int icheck = 0;
    int n, m, i, j, n1, n2, m1, m2, i1, i2;
    
    /* Search for elements around a node */
    while (1) {
        nean = allocate_2d_int(np, nean_max);
        if (!nean) {
            free(numean);
            return -1;
        }
        
        memset(numean, 0, np * sizeof(int));
        int restart = 0;
        
        for (m = 0; m < ne; m++) {
            for (i = 0; i < 3; i++) {
                n = nm[m][i];
                numean[n]++;
                if (numean[n] > nean_max) {
                    icheck++;
                    if (icheck == 1) {
                        printf("      Salvere000: System require much more nean_max\n");
                    }
                    nean_max += 5;
                    printf("                  System try to re-calculate with nean_max=%d\n", nean_max);
                    restart = 1;
                    break;
                }
                nean[n][numean[n] - 1] = m;
            }
            if (restart) break;
        }
        
        if (!restart) break;
        free_2d_int(nean, np);
    }
    
    if (icheck != 0) {
        printf("                  Success!\n\n\n");
    }
    
    /* Search the border of global grid and sub-grid */
    *nbl = 0;
    for (m = 0; m < ne; m++) {
        if (netable[m] == 1) {
            for (i = 0; i < 3; i++) {
                eline(i, &n1, &n2, nm, m);
                icheck = 0;
                
                for (i1 = 0; i1 < numean[n1]; i1++) {
                    m1 = nean[n1][i1];
                    for (i2 = 0; i2 < numean[n2]; i2++) {
                        m2 = nean[n2][i2];
                        if (m1 == m2) {
                            if (m1 != m) {
                                if (netable[m1] == 1) {
                                    icheck++;
                                    break;
                                }
                            }
                        }
                    }
                    if (icheck != 0) break;
                }
                
                if (icheck == 0) {
                    nblnc[0][*nbl] = n1;
                    nblnc[1][*nbl] = n2;
                    (*nbl)++;
                }
            }
        }
    }
    
    free(numean);
    free_2d_int(nean, np);
    return 0;
}

int readetab(int ne, int *netable) {
    FILE *fp = fopen(".tmp01", "r");
    if (!fp) return -1;
    
    int icdigit, n1 = 0, i, m, n, j;
    
    if (fscanf(fp, "%d", &icdigit) != 1) {
        fclose(fp);
        return -1;
    }
    
    int *nw = (int *)malloc(icdigit * sizeof(int));
    if (!nw) {
        fclose(fp);
        return -1;
    }
    
    memset(netable, 0, ne * sizeof(int));
    
    m = (ne % icdigit == 0) ? ne / icdigit : ne / icdigit + 1;
    
    for (i = 0; i < m; i++) {
        if (fscanf(fp, "%d", &n) != 1) break;
        
        for (j = 0; j < icdigit; j++) {
            nw[j] = n % 2;
            n = n / 2;
        }
        
        for (j = icdigit - 1; j >= 0; j--) {
            if (n1 >= ne) break;
            netable[n1] = nw[j];
            n1++;
        }
    }
    
    free(nw);
    fclose(fp);
    return 0;
}

void mkareacood(int ne, int np, double **xyd, int **nm, int m, 
                double xc, double yc, double *sl) {
    sl[0] = sl[1] = sl[2] = 0.0;
    
    double x1 = xyd[0][nm[m][0]];
    double x2 = xyd[0][nm[m][1]];
    double x3 = xyd[0][nm[m][2]];
    double y1 = xyd[1][nm[m][0]];
    double y2 = xyd[1][nm[m][1]];
    double y3 = xyd[1][nm[m][2]];
    
    double a0 = (x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2);
    double a1 = ((xc - x2) * (yc - y3) - (xc - x3) * (yc - y2)) / a0;
    double a2 = ((x1 - xc) * (y1 - y3) - (x1 - x3) * (y1 - yc)) / a0;
    double a3 = ((x1 - x2) * (y1 - yc) - (x1 - xc) * (y1 - y2)) / a0;
    
    sl[0] = a1;
    sl[1] = a2;
    sl[2] = a3;
}
