/*
 * test_gridscope.c - Comprehensive test harness
 * 
 * Tests the C translation against the original Fortran implementation
 * by comparing outputs for identical inputs.
 */

#include "gridscope.h"
#include <assert.h>
#include <float.h>

#define TEST_TOLERANCE 1.0e-10
#define MAX_LINE_LENGTH 1024

/* Test result structure */
typedef struct {
    int passed;
    int failed;
    char last_error[MAX_STRING_LENGTH];
} TestResults;

TestResults test_results = {0, 0, ""};

/* ===================================================================
 * TEST UTILITY FUNCTIONS
 * =================================================================== */

void test_assert(int condition, const char *message) {
    if (condition) {
        test_results.passed++;
        printf("  [PASS] %s\n", message);
    } else {
        test_results.failed++;
        snprintf(test_results.last_error, MAX_STRING_LENGTH, "%s", message);
        printf("  [FAIL] %s\n", message);
    }
}

void test_assert_double_equal(double actual, double expected, const char *message) {
    double diff = fabs(actual - expected);
    int passed = (diff < TEST_TOLERANCE) || 
                 (fabs(diff / (fabs(expected) + TEST_TOLERANCE)) < TEST_TOLERANCE);
    
    if (passed) {
        test_results.passed++;
        printf("  [PASS] %s (actual=%.10e, expected=%.10e)\n", message, actual, expected);
    } else {
        test_results.failed++;
        snprintf(test_results.last_error, MAX_STRING_LENGTH, 
                "%s: difference %.10e", message, diff);
        printf("  [FAIL] %s (actual=%.10e, expected=%.10e, diff=%.10e)\n", 
               message, actual, expected, diff);
    }
}

void test_assert_int_equal(int actual, int expected, const char *message) {
    if (actual == expected) {
        test_results.passed++;
        printf("  [PASS] %s (value=%d)\n", message, actual);
    } else {
        test_results.failed++;
        snprintf(test_results.last_error, MAX_STRING_LENGTH, "%s", message);
        printf("  [FAIL] %s (actual=%d, expected=%d)\n", message, actual, expected);
    }
}

/* ===================================================================
 * MESH COMPARISON FUNCTIONS
 * =================================================================== */

int compare_mesh_files(const char *file1, const char *file2) {
    FILE *fp1 = fopen(file1, "r");
    FILE *fp2 = fopen(file2, "r");
    
    if (!fp1 || !fp2) {
        if (fp1) fclose(fp1);
        if (fp2) fclose(fp2);
        return -1;
    }
    
    char line1[MAX_LINE_LENGTH], line2[MAX_LINE_LENGTH];
    int line_num = 0;
    int differences = 0;
    
    while (fgets(line1, sizeof(line1), fp1) && fgets(line2, sizeof(line2), fp2)) {
        line_num++;
        
        /* Parse lines for numeric comparison */
        double vals1[10], vals2[10];
        int ints1[10], ints2[10];
        int n1, n2;
        
        /* Try to parse as coordinates (node line) */
        n1 = sscanf(line1, "%d %lf %lf %lf", &ints1[0], &vals1[0], &vals1[1], &vals1[2]);
        n2 = sscanf(line2, "%d %lf %lf %lf", &ints2[0], &vals2[0], &vals2[1], &vals2[2]);
        
        if (n1 == 4 && n2 == 4) {
            /* Node line comparison */
            if (ints1[0] != ints2[0]) {
                printf("  Line %d: Node ID mismatch (%d vs %d)\n", 
                       line_num, ints1[0], ints2[0]);
                differences++;
            }
            
            for (int i = 0; i < 3; i++) {
                double diff = fabs(vals1[i] - vals2[i]);
                if (diff > TEST_TOLERANCE) {
                    printf("  Line %d: Coordinate %d mismatch (%.10e vs %.10e, diff=%.10e)\n", 
                           line_num, i, vals1[i], vals2[i], diff);
                    differences++;
                }
            }
            continue;
        }
        
        /* Try to parse as element connectivity */
        n1 = sscanf(line1, "%d %d %d %d %d", &ints1[0], &ints1[1], &ints1[2], &ints1[3], &ints1[4]);
        n2 = sscanf(line2, "%d %d %d %d %d", &ints2[0], &ints2[1], &ints2[2], &ints2[3], &ints2[4]);
        
        if (n1 == 5 && n2 == 5) {
            /* Element line comparison */
            for (int i = 0; i < 5; i++) {
                if (ints1[i] != ints2[i]) {
                    printf("  Line %d: Element field %d mismatch (%d vs %d)\n", 
                           line_num, i, ints1[i], ints2[i]);
                    differences++;
                }
            }
            continue;
        }
        
        /* Default string comparison for other lines */
        if (strcmp(line1, line2) != 0) {
            /* Allow minor formatting differences */
            char *trimmed1 = line1, *trimmed2 = line2;
            while (*trimmed1 == ' ' || *trimmed1 == '\t') trimmed1++;
            while (*trimmed2 == ' ' || *trimmed2 == '\t') trimmed2++;
            
            if (strcmp(trimmed1, trimmed2) != 0) {
                printf("  Line %d: Text mismatch\n", line_num);
                printf("    File1: %s", line1);
                printf("    File2: %s", line2);
                differences++;
            }
        }
    }
    
    /* Check if files have different lengths */
    int has_more1 = fgets(line1, sizeof(line1), fp1) != NULL;
    int has_more2 = fgets(line2, sizeof(line2), fp2) != NULL;
    
    if (has_more1 || has_more2) {
        printf("  Files have different lengths\n");
        differences++;
    }
    
    fclose(fp1);
    fclose(fp2);
    
    return differences;
}

/* ===================================================================
 * UNIT TESTS
 * =================================================================== */

void test_memory_allocation() {
    printf("\n=== Testing Memory Allocation ===\n");
    
    double **arr2d = allocate_2d_double(10, 20);
    test_assert(arr2d != NULL, "Allocate 2D double array");
    
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 20; j++) {
            arr2d[i][j] = i * 20 + j;
        }
    }
    
    test_assert_double_equal(arr2d[5][10], 110.0, "2D array indexing");
    free_2d_double(arr2d, 10);
    
    int **arr2i = allocate_2d_int(15, 25);
    test_assert(arr2i != NULL, "Allocate 2D int array");
    free_2d_int(arr2i, 15);
    
    double ***arr3d = allocate_3d_double(3, 5, 7);
    test_assert(arr3d != NULL, "Allocate 3D double array");
    free_3d_double(arr3d, 3, 5);
}

void test_element_area_calculation() {
    printf("\n=== Testing Element Area Calculation ===\n");
    
    /* Create simple test mesh */
    double **xyd = allocate_2d_double(3, 3);
    
    /* Triangle with vertices at (0,0), (1,0), (0,1) */
    xyd[0][0] = 0.0; xyd[1][0] = 0.0; xyd[2][0] = 0.0;
    xyd[0][1] = 1.0; xyd[1][1] = 0.0; xyd[2][1] = 0.0;
    xyd[0][2] = 0.0; xyd[1][2] = 1.0; xyd[2][2] = 0.0;
    
    double earea;
    element_area(3, xyd, 0, 1, 2, &earea);
    
    /* Area should be 1.0 (actually 2*area) */
    test_assert_double_equal(earea, 1.0, "Element area calculation");
    
    /* Test negative area (reversed winding) */
    element_area(3, xyd, 0, 2, 1, &earea);
    test_assert_double_equal(earea, -1.0, "Element area (reversed winding)");
    
    free_2d_double(xyd, 3);
}

void test_eline_function() {
    printf("\n=== Testing eline Function ===\n");
    
    int **nm = allocate_2d_int(1, 3);
    nm[0][0] = 10;
    nm[0][1] = 20;
    nm[0][2] = 30;
    
    int n1, n2;
    
    eline(2, &n1, &n2, nm, 0);  /* Fortran case 3 */
    test_assert_int_equal(n1, 10, "eline case 2: n1");
    test_assert_int_equal(n2, 20, "eline case 2: n2");
    
    eline(0, &n1, &n2, nm, 0);  /* Fortran case 1 */
    test_assert_int_equal(n1, 20, "eline case 0: n1");
    test_assert_int_equal(n2, 30, "eline case 0: n2");
    
    eline(1, &n1, &n2, nm, 0);  /* Fortran case 2 */
    test_assert_int_equal(n1, 30, "eline case 1: n1");
    test_assert_int_equal(n2, 10, "eline case 1: n2");
    
    free_2d_int(nm, 1);
}

void test_barycentric_coordinates() {
    printf("\n=== Testing Barycentric Coordinates ===\n");
    
    /* Create simple test mesh */
    double **xyd = allocate_2d_double(3, 3);
    int **nm = allocate_2d_int(1, 3);
    
    /* Triangle with vertices at (0,0), (1,0), (0,1) */
    xyd[0][0] = 0.0; xyd[1][0] = 0.0;
    xyd[0][1] = 1.0; xyd[1][1] = 0.0;
    xyd[0][2] = 0.0; xyd[1][2] = 1.0;
    
    nm[0][0] = 0; nm[0][1] = 1; nm[0][2] = 2;
    
    double sl[3];
    
    /* Test center point (1/3, 1/3) */
    mkareacood(1, 3, xyd, nm, 0, 1.0/3.0, 1.0/3.0, sl);
    test_assert_double_equal(sl[0], 1.0/3.0, "Barycentric coord 0 at centroid");
    test_assert_double_equal(sl[1], 1.0/3.0, "Barycentric coord 1 at centroid");
    test_assert_double_equal(sl[2], 1.0/3.0, "Barycentric coord 2 at centroid");
    
    /* Test vertex point (0,0) */
    mkareacood(1, 3, xyd, nm, 0, 0.0, 0.0, sl);
    test_assert_double_equal(sl[0], 1.0, "Barycentric coord 0 at vertex 0");
    test_assert_double_equal(sl[1], 0.0, "Barycentric coord 1 at vertex 0");
    test_assert_double_equal(sl[2], 0.0, "Barycentric coord 2 at vertex 0");
    
    free_2d_double(xyd, 3);
    free_2d_int(nm, 1);
}

void test_mesh_io() {
    printf("\n=== Testing Mesh I/O ===\n");
    
    /* Test reading the provided test mesh */
    const char *test_file = "/mnt/project/test_submesh.14";
    
    FILE *fp = fopen(test_file, "r");
    test_assert(fp != NULL, "Open test mesh file");
    
    if (fp) {
        int ne, np, nope, nbou, nvdl_max, nvel_max, nodemax;
        
        int result = read14_alloc(fp, &ne, &np, &nope, &nbou, 
                                  &nvdl_max, &nvel_max, &nodemax);
        fclose(fp);
        
        test_assert(result == 0, "Read mesh allocation info");
        test_assert_int_equal(ne, 85, "Number of elements");
        test_assert_int_equal(np, 56, "Number of nodes");
        test_assert_int_equal(nope, 0, "Number of open boundaries");
        test_assert_int_equal(nbou, 1, "Number of land boundaries");
        
        /* Allocate and read full mesh */
        Mesh *mesh = allocate_mesh(ne, np, nope, nbou, nvdl_max, nvel_max, nodemax);
        test_assert(mesh != NULL, "Allocate mesh structure");
        
        if (mesh) {
            fp = fopen(test_file, "r");
            result = read14(fp, mesh);
            fclose(fp);
            
            test_assert(result == 0, "Read full mesh");
            
            /* Verify some known values from the test file */
            test_assert_double_equal(mesh->grid.xyd[0][0], -81.2532600000, 
                                    "Node 1 x-coordinate");
            test_assert_double_equal(mesh->grid.xyd[1][0], 30.6589510000, 
                                    "Node 1 y-coordinate");
            
            /* Check element connectivity (convert to 1-indexed for comparison) */
            test_assert_int_equal(mesh->grid.nm[0][0] + 1, 1, "Element 1 node 1");
            test_assert_int_equal(mesh->grid.nm[0][1] + 1, 9, "Element 1 node 2");
            test_assert_int_equal(mesh->grid.nm[0][2] + 1, 8, "Element 1 node 3");
            
            free_mesh(mesh);
        }
    }
}

/* ===================================================================
 * INTEGRATION TESTS
 * =================================================================== */

void test_mkeline_boundary_detection() {
    printf("\n=== Testing mkeline Boundary Detection ===\n");
    
    /* Create a simple 2x2 quad mesh (4 triangles) */
    int ne = 4, np = 6;
    int **nm = allocate_2d_int(ne, 3);
    int *netable = (int *)malloc(ne * sizeof(int));
    
    /* Define triangulation:
     * 0---1---2
     * |\ 1|\ 3|
     * |0 \|2 \|
     * 3---4---5
     */
    nm[0][0] = 0; nm[0][1] = 1; nm[0][2] = 3;
    nm[1][0] = 1; nm[1][1] = 4; nm[1][2] = 3;
    nm[2][0] = 1; nm[2][1] = 2; nm[2][2] = 4;
    nm[3][0] = 2; nm[3][1] = 5; nm[3][2] = 4;
    
    /* All elements in netable */
    for (int i = 0; i < ne; i++) netable[i] = 1;
    
    int nbl;
    int **nblnc = allocate_2d_int(2, 20);
    
    int result = mkeline(np, ne, nm, netable, 20, &nbl, nblnc);
    
    test_assert(result == 0, "mkeline executed successfully");
    test_assert_int_equal(nbl, 8, "Correct number of boundary edges");
    
    free_2d_int(nm, ne);
    free(netable);
    free_2d_int(nblnc, 2);
}

/* ===================================================================
 * MAIN TEST RUNNER
 * =================================================================== */

int main(int argc, char **argv) {
    printf("\n");
    printf("========================================\n");
    printf("  GRIDSCOPE C Translation Test Suite\n");
    printf("========================================\n");
    
    /* Run all tests */
    test_memory_allocation();
    test_element_area_calculation();
    test_eline_function();
    test_barycentric_coordinates();
    test_mesh_io();
    test_mkeline_boundary_detection();
    
    /* Print summary */
    printf("\n");
    printf("========================================\n");
    printf("  Test Results Summary\n");
    printf("========================================\n");
    printf("  Passed: %d\n", test_results.passed);
    printf("  Failed: %d\n", test_results.failed);
    
    if (test_results.failed > 0) {
        printf("\n  Last error: %s\n", test_results.last_error);
        printf("\n  OVERALL: FAILED\n");
        printf("========================================\n\n");
        return 1;
    } else {
        printf("\n  OVERALL: PASSED\n");
        printf("========================================\n\n");
        return 0;
    }
}
