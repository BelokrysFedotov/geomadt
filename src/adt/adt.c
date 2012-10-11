/*****************************************************************************/
/*                                                                           */
/* adt.c                                                                     */
/*                                                                           */
/* List of functions that operate with Alternating Digital Tree (ADT):       */
/*                                                                           */
/*      adt_construct   - construct ADT (with no mask)                       */
/*      adt_init        - construct ADT (with mask)                          */
/*      adt_search      - search in ADT                                      */
/*      adt_free        - free ADT                                           */
/*                                                                           */
/* Version 14.06.01                                                          */
/*                                                                           */
/* Problem to be solved:                                                     */
/*      find all parallelipipeds from the tree                               */
/*      that intersect the specified one                                     */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*
 * To insert lots of self-checks for internal errors as well as self-test,
 * define the SELF_CHECK symbol.
 * It is best to define the symbol using the -DSELF_CHECK compiler switch,
 * but you could write "#define SELF_CHECK" below.
 */

/* #define SELF_CHECK */

/*---------------------------------------------------------------------------*
 * To inline local static functions with GNU C compiler
 * comment the definition of inline symbol below.
 * For compatibility with strict ANSI sintax (to prevent inlining)
 * inline symbol should be defined: "#define inline".
 */

#define inline

/*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "adt.h"

#ifdef UNIX
  #define STDERR stderr
#else
  #define STDERR stdout
#endif


/*
 * Explicit declaration of functions malloc() and free()
 * because some non-ANSI C compilers lack stdlib.h
 */

#ifndef _STDLIB_H_
extern void *malloc();
extern void free();
#endif

/*---------------------------------------------------------------------------*/

/*===========================================================================*/

struct adtnode                  /* node of Alternating Digital Tree (ADT) */
{
    struct adtnode *left;       /* left child */
    struct adtnode *center;     /* center child */
    struct adtnode *right;      /* right child */
    REAL *region;               /* pointer to region coordinates */
    int index;                  /* index in the user's list of objects */
};                              /* Remark: index can be computed from pointers
                                   or from above hystory of computations */

struct adtree                   /* Alternating Digital Tree (ADT) */
{
    struct adtnode *root;       /* pointer to the root of the tree:
                                   usualy equal to that of nodes */
    int dim;                    /* dimension of the user's space:
                                   2 or 3 */
    REAL *domain;               /* pointer to initial domain region coordinates:
                                   dimension [2][dim] */
                        /* ------- additional information: */
    int nnodes;                 /* number of nodes in the tree */
    struct adtnode *nodes;      /* pointer to array of adtnode structures:
                                   dimension [nnodes] */
    REAL *regions;              /* pointer to array of region coordinates:
                                   dimension [nnodes][1+2*dim] */
    int height;                 /* actual height of the tree */
    int height_opt;             /* optimal height of the tree: log(nnodes) */
                        /* ------- the other working storage: */
    int *indices;               /* indices found in search: dimension[nnodes] */
    struct adtnode **stack;     /* stack of levels: dimension [height+1] */
    REAL *stack_coord;          /* stack of coordinates: dimension [height+1] */
    REAL *subdomain;            /* pointer to local subdomain coordinates:
                                   dimension [2][dim] */
    REAL *subregion;            /* pointer to local region coordinates:
                                   dimension [2][dim] */
};

/*---------------------------------------------------------------------------*
 *      Add a node to ADT.
 */
inline static void adt_add_node(const int i, struct adtree *adt)
{
    struct adtnode *ptr;
    register int dim, j;
    REAL *reg, *dom, center;

    dim = adt->dim;
    ptr = adt->root;
    dom = adt->subdomain;
    reg = (adt->nodes[i]).region;
    for (j = 0; j < 2*dim; j++)
        dom[j] = adt->domain[j];

    if (ptr == NULL)
    {
        adt->root = &adt->nodes[i];
        j = 0;
    }
    else
    {
        for (j = 0; ; j = (++j < dim) ? j : 0)          /* j = {0,...,dim-1} */
        {
            center = *ptr->region;
            if (reg[1+j+dim] < center)
            {
                /* go to left */
                if (ptr->left != NULL)
                {
                    dom[j+dim] = center;
                    ptr = ptr->left;
                }
                else
                {
                    ptr->left = &adt->nodes[i];
                    break;
                }
            }
            else if (reg[1+j] > center)
            {
                /* go to right */
                if (ptr->right != NULL)
                {
                    dom[j] = center;
                    ptr = ptr->right;
                }
                else
                {
                    ptr->right = &adt->nodes[i];
                    break;
                }
            }
            else
            {
                /* go to center */
                if (ptr->center != NULL)
                {
                    ptr = ptr->center;
                }
                else
                {
                    ptr->center = &adt->nodes[i];
                    break;
                }
            }
        }
        j = (++j < dim) ? j : 0;
    }
    reg[0] = (dom[j] + dom[j+dim]) * 0.5e0;
    return;
}

/*---------------------------------------------------------------------------*
 *      Set min and max coordinates of the region for ADT.
 */
inline static void adt_minmax
        (const int dim, const REAL point1[], const REAL point2[], REAL region[])
{
    REAL t1, t2;
    int j;

    for (j = 0; j < dim; j++)
    {
        t1 = point1[j];
        t2 = point2[j];
        if (t1 < t2)
        {
            region[j] = t1;
            region[j+dim] = t2;
        }
        else
        {
            region[j] = t2;
            region[j+dim] = t1;
        }
    }
}

/*---------------------------------------------------------------------------*
 *      Compute height of ADT branch.
 */
inline static int adt_height(struct adtnode *node)
{
    int il, ic, ir, ih;
/*
 *  if (node == NULL)
 *      return 0;
 *  il = adt_height(node->left);
 *  ic = adt_height(node->center);
 *  ir = adt_height(node->right);
 *  ih = (il > ir) ? il : ir;
 *  ih = (ih > ic) ? ih : ic;
 *  return ++ih;
 */
    return (node == NULL) ? 0 :
                (ih = (il = adt_height(node->left)) >
                      (ir = adt_height(node->right)) ? il : ir) >
                      (ic = adt_height(node->center)) ? ++ih : ++ic;
}

/*---------------------------------------------------------------------------*
 *      Two regions are intersected, is not it?
 */
inline static int adt_intersect(register const int dim,
        const REAL reg1[], const REAL reg2[])
{
    register int j;

    for (j = dim; --j >= 0; )           /* fast version */
/*  for (j = 0; j < dim; j++) */        /* slow version */
    {
        if (reg1[j] > reg2[j+dim] ||
            reg2[j] > reg1[j+dim])
        {
                return 0;
        }
    }
    return 1;
}

/*****************************************************************************/
/*      External functions                                                   */
/*****************************************************************************/

/*---------------------------------------------------------------------------*
 *      Initialize ADT (with mask).
 */
void *adt_init
        (int nnodes, const int dim, int *mask, void *user_data,
    void (*adt_region)(int i, REAL reg_min[], REAL reg_max[], void *user_data))
{
    int i, j, nnodes_total;
    REAL *dom, *reg;
    struct adtnode *node;
    struct adtree *adt;

    /* Find actual number of nodes in the tree */
    nnodes_total = nnodes;
    if (mask)
        for (nnodes = 0, i = 0; i < nnodes_total; i++)
            if (mask[i])
                nnodes++;

    /* Allocate the tree and initial setup */
    adt        = (struct adtree *)  malloc(sizeof(struct adtree));
    if (adt == NULL)
        return NULL;
    adt->nodes = (struct adtnode *) malloc(sizeof(struct adtnode) * nnodes);
    adt->domain    = (REAL *)       malloc(sizeof(REAL) *    2*dim);
    adt->regions   = (REAL *)       malloc(sizeof(REAL) * (1+2*dim) * nnodes);
    adt->subdomain = (REAL *)       malloc(sizeof(REAL) *    2*dim);
    adt->subregion = (REAL *)       malloc(sizeof(REAL) *    2*dim);
    adt->dim       = dim;
    adt->nnodes    = nnodes;
    adt->root      = NULL;
    if (adt->nodes     == NULL ||
        adt->domain    == NULL ||
        adt->regions   == NULL ||
        adt->subdomain == NULL ||
        adt->subregion == NULL)
            return NULL;

    /* Compute the regions for the tree nodes and initial setup */
    dom  = adt->domain;
    reg  = adt->regions;
    node = adt->nodes;
    (*adt_region)(0, dom+0, dom+dim, user_data);
    adt_minmax(dim, dom, dom+dim, dom);
    for (i = 0; i < nnodes_total; i++)
    {
        if (mask && !mask[i])
            continue;

        node->index  = i;
        node->left   = NULL;
        node->center = NULL;
        node->right  = NULL;
        node->region = reg;
        (*adt_region)(i, reg+1, reg+1+dim, user_data);

        /* Find min and max coordinates:
         * the first half containing dim coondinates specify min
         * and the second half of coordinates specify max
         */
        adt_minmax(dim, reg+1, reg+1+dim, reg+1);
        for (j = 0; j < dim; j++)
            dom[j] = (reg[1+j] < dom[j]) ? reg[1+j] : dom[j];   /* min */
        for (j = dim; j < 2*dim; j++)
            dom[j] = (reg[1+j] > dom[j]) ? reg[1+j] : dom[j];   /* max */
        reg += 1 + 2*dim;
        node++;
    }

    /* Constuct the tree */
    for (i = 0; i < nnodes; i++)
    {
        adt_add_node(i, adt);
    }

    /* Compute height of the tree */
    /* log2(nnodes) = ln nnodes / ln 2 */
    adt->height     = adt_height(adt->root);
    adt->height_opt = (int)ceil(log((float)nnodes+1) / log(2.0));

    /* Allocate stack */
    adt->indices     = (int *) malloc(sizeof(int) * nnodes);
    adt->stack       = (struct adtnode **) malloc(sizeof(struct adtnode *)
                                                    * (adt->height+1));
    adt->stack_coord = (REAL *) malloc(sizeof(REAL) * (adt->height+1));
    if (adt->indices     == NULL ||
        adt->stack       == NULL ||
        adt->stack_coord == NULL)
            return NULL;
    return (void *)adt;
}

/*---------------------------------------------------------------------------*
 *      Construct ADT.
 */
void *adt_construct
        (const int nnodes, const int dim, void *user_data,
    void (*adt_region)(int i, REAL reg_min[], REAL reg_max[], void *user_data))
{
    int i, j;
    REAL *dom, *reg;
    struct adtnode *node;
    struct adtree *adt;

    /* Allocate the tree and initial setup */
    adt        = (struct adtree *)  malloc(sizeof(struct adtree));
    if (adt == NULL)
        return NULL;
    adt->nodes = (struct adtnode *) malloc(sizeof(struct adtnode) * nnodes);
    adt->domain    = (REAL *)       malloc(sizeof(REAL) *    2*dim);
    adt->regions   = (REAL *)       malloc(sizeof(REAL) * (1+2*dim) * nnodes);
    adt->subdomain = (REAL *)       malloc(sizeof(REAL) *    2*dim);
    adt->subregion = (REAL *)       malloc(sizeof(REAL) *    2*dim);
    adt->dim       = dim;
    adt->nnodes    = nnodes;
    adt->root      = NULL;
    if (adt->nodes     == NULL ||
        adt->domain    == NULL ||
        adt->regions   == NULL ||
        adt->subdomain == NULL ||
        adt->subregion == NULL)
            return NULL;

    /* Compute the regions for the tree nodes and initial setup */
    dom  = adt->domain;
    reg  = adt->regions;
    node = adt->nodes;
    (*adt_region)(0, dom+0, dom+dim, user_data);
    adt_minmax(dim, dom, dom+dim, dom);
    for (i = 0; i < nnodes; i++)
    {
        node->index  = i;
        node->left   = NULL;
        node->center = NULL;
        node->right  = NULL;
        node->region = reg;
        (*adt_region)(i, reg+1, reg+1+dim, user_data);

        /* Find min and max coordinates:
         * the first half containing dim coondinates specify min
         * and the second half of coordinates specify max
         */
        adt_minmax(dim, reg+1, reg+1+dim, reg+1);
        for (j = 0; j < dim; j++)
            dom[j] = (reg[1+j] < dom[j]) ? reg[1+j] : dom[j];   /* min */
        for (j = dim; j < 2*dim; j++)
            dom[j] = (reg[1+j] > dom[j]) ? reg[1+j] : dom[j];   /* max */
        reg += 1 + 2*dim;
        node++;
    }

    /* Constuct the tree */
    for (i = 0; i < nnodes; i++)
    {
        adt_add_node(i, adt);
    }

    /* Compute height of the tree */
    /* log2(nnodes) = ln nnodes / ln 2 */
    adt->height     = adt_height(adt->root);
    adt->height_opt = (int)ceil(log((float)nnodes+1) / log(2.0));

    /* Allocate stack */
    adt->indices     = (int *) malloc(sizeof(int) * nnodes);
    adt->stack       = (struct adtnode **) malloc(sizeof(struct adtnode *)
                                                    * (adt->height+1));
    adt->stack_coord = (REAL *) malloc(sizeof(REAL) * (adt->height+1));
    if (adt->indices     == NULL ||
        adt->stack       == NULL ||
        adt->stack_coord == NULL)
            return NULL;
    return (void *)adt;
}

/*---------------------------------------------------------------------------*
 *      Search in ADT.
 */
int adt_search(int **indices, const REAL reg_min[], const REAL reg_max[],
               const void *v_adt)
{
    struct adtree *adt = (struct adtree *)v_adt;
    register int i, j, dim;
    int n;
    REAL *reg, *dom, center;
    struct adtnode *ptr;

    if (adt == NULL)
        return -1;
    dim = adt->dim;
    ptr = adt->root;
    dom = adt->subdomain;
    reg = adt->subregion;
    for (j = 0; j < 2*dim; j++)
        dom[j] = adt->domain[j];
    adt_minmax(dim, reg_min, reg_max, reg);
    for (j = 0; j <= adt->height; j++)
        adt->stack[j] = adt->root;
    n = i = j = 0; /* j = {0,...,dim-1} */
    while (1)
    {
        center = *ptr->region;
        if (ptr->left   != NULL            &&
            ptr->left   != adt->stack[i+1] &&
            ptr->right  != adt->stack[i+1] &&
            ptr->center != adt->stack[i+1] &&
            reg[j]      <  center)
        {
            /* goto left */
            adt->stack_coord[i] = dom[j+dim];
            dom[j+dim] = center;
            adt->stack[i] = ptr;
            ptr = ptr->left;
            i++;
            j = (++j < dim) ? j : 0;
        }
        else if (ptr->right  != NULL            &&
                 ptr->right  != adt->stack[i+1] &&
                 ptr->center != adt->stack[i+1] &&
                 reg[j+dim]  >  center)
        {
            /* goto right */
            adt->stack_coord[i] = dom[j];
            dom[j] = center;
            adt->stack[i] = ptr;
            ptr = ptr->right;
            i++;
            j = (++j < dim) ? j : 0;
        }
        else if (ptr->center != NULL            &&
                 ptr->center != adt->stack[i+1])
        {
            /* goto center */
            adt->stack[i] = ptr;
            ptr = ptr->center;
            i++;
            j = (++j < dim) ? j : 0;
        }
        else
        {
            /* check intersection */
            if (adt_intersect(dim, reg, ptr->region+1))
            {
                adt->indices[n++] = ptr->index;
            }
            /* return to up */
            i--;
            j = (--j < 0) ? dim-1 : j;
            if (i < 0) /* we are at the root */
            {
                break;
            }
            else if (adt->stack[i]->left == ptr) /* return from left */
            {
                dom[j+dim] = adt->stack_coord[i];
            }
            else if (adt->stack[i]->right == ptr) /* return from right */
            {
                dom[j] = adt->stack_coord[i];
            }
            adt->stack[i+1] = ptr;
            ptr = adt->stack[i];
        }
    }
    *indices = adt->indices;
    return n;
}

/*---------------------------------------------------------------------------*
 *      Free ADT.
 *      Free all arrays of ADT in reverse order.
 *      Note: the pointer to ADT is not modified (is not set to NULL).
 */
void adt_free(const void *v_adt)
{
    

    struct adtree *adt = (struct adtree *)v_adt;
    
    if (v_adt == NULL) return;

    free(adt->stack_coord);
    free(adt->stack);
    free(adt->indices);
    free(adt->subregion);
    free(adt->subdomain);
    free(adt->regions);
    free(adt->domain);
    free(adt->nodes);
    free(adt);
}

/*****************************************************************************/
/*      Self check                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*
 *      Full search in list of ADT nodes.
 */
#if defined SELF_CHECK || defined ADT_FULL_SEARCH
int adt_full_search(int **indices, const REAL reg_min[], const REAL reg_max[],
                    const void *v_adt)
{
    struct adtree *adt = (struct adtree *)v_adt;
    int dim, i, n = 0;
    struct adtnode *ptr;

    if (adt == NULL)
        return -1;
    dim = adt->dim;
    adt_minmax(dim, reg_min, reg_max, adt->subregion);

    for (i = 0; i < adt->nnodes; i++)
    {
        ptr = adt->nodes+i;
        if (adt_intersect(dim, adt->subregion, ptr->region+1))
            adt->indices[n++] = ptr->index;
    }
    *indices = adt->indices;
    return n;
}
#endif /* SELF_CHECK || ADT_FULL_SEARCH */

/*===========================================================================*
 *      Example of test programs
 */
#ifdef SELF_CHECK

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "adt.h"

#ifndef REAL
#define REAL double
#endif

#ifndef REAL2
typedef REAL REAL2[2];
#endif

#ifndef REAL3
typedef REAL REAL3[3];
#endif

void main_t1();
void main_t2();
void main_t3();
void main_t4();
void main_test();

int main() { main_t4(); return 0;}

void main_t1()
{
    int dim = 2;
    REAL points[][2] = {           {4,4},                                  /*1*/
                 {2,4},                             {6,4},                 /*3*/
        {2,2},           {2,6},            {6,2},            {6,6},        /*7*/
   {1,2},   {3,2},   {1,6},   {3,6},   {5,2},   {7,2},   {5,6},   {7,6},  /*15*/
/*1,1,1,3, 3,1,3,3, 1,5,1,7, 3,5,3,7, 5,1,5,3, 7,1,7,3, 5,5,5,7, 7,5,7,7*//*31*/
                        };

                /* x2 ^                                                    */
                /*  8 +-----------------------+           points           */
                /*  7 | 20    22    28    30  |             0         j=0  */
                /*  6 |  9  4 10    13  6 14  |           /   \            */
                /*  5 | 19    21    27    29  |         1       2     j=1  */
                /*  4 |     1     0     2     |        / \     / \         */
                /*  3 | 16    18    24    26  |       3   4   5   6   j=0  */
                /*  2 |  7  3  8    11  5 12  |      / \ / \ / \ / \       */
                /*  1 | 15    17    23    25  |      . . . . . . . .       */
                /*  0 +-----------------------+ ->                         */
                /*    0 01 02 03 04 05 06 07 08 x1                         */

    int nnodes;

    nnodes = sizeof(points) / sizeof(REAL) / dim;
    main_test(dim, nnodes, points);
}

void main_t2()
{
    int dim = 2;
    int m = 30;
    int nnodes = m * m;
    REAL2 *points;
    int i, j;

    points = (REAL2 *) malloc(sizeof(REAL2) * nnodes);

    for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
    {
        points[i+j*m][0] = i;
        points[i+j*m][1] = j;
    }
    main_test(dim, nnodes, points);
}

void main_t3()
{
    int dim = 3;
    int m = 10;
    int nnodes = m * m * m;
    REAL3 *points;
    int i, j, k;

    points = (REAL3 *) malloc(sizeof(REAL3) * nnodes);

    for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
    for (k = 0; k < m; k++)
    {
        points[i+j*m+k*m*m][0] = i;
        points[i+j*m+k*m*m][1] = j;
        points[i+j*m+k*m*m][2] = k;
    }
    main_test(dim, nnodes, points);
}

void main_t4()
{
    int dim = 3;
    int nnodes = 1023;
    REAL3 *points;
    int i;

    points = (REAL3 *) malloc(sizeof(REAL3) * nnodes);

    printf("Input seed for srand: ");
    scanf("%d",&i);
    srand(i);
    for (i = 0; i < nnodes; i++)
    {
        points[i][0] = (float)rand() / RAND_MAX;
        points[i][1] = (float)rand() / RAND_MAX;
        points[i][2] = (float)rand() / RAND_MAX;
    }
    main_test(dim, nnodes, points);
}

void main_test(int dim, int nnodes, REAL points[])
{
    void adt_region2(int i, REAL reg_min[2], REAL reg_max[2], void *user_data);
    void adt_region3(int i, REAL reg_min[3], REAL reg_max[3], void *user_data);
    void *adt;
    int i, ii, n, nf, *indices;
    double t;

    printf("main_test: nnodes = %i\n",nnodes);
    t = clock();
    if (dim == 2)
        adt = adt_construct(nnodes, dim, (void *)points, adt_region2);
    else if (dim == 3)
        adt = adt_construct(nnodes, dim, (void *)points, adt_region3);
    t = (clock() - t) / CLOCKS_PER_SEC;
    printf("adt_construct time = %g sec\n",t);
    if (adt == NULL)
    {
        printf("Error in adt_construct\n");
        exit(1);
    }

    /* Check adt_search for random values */
    if (1)
    {
        for (i = 0; i < nnodes; i++)
        {
            for (ii = 0; ii < 2*dim; ii++)
                points[ii] = (float)rand() / RAND_MAX;
            n = adt_search(&indices, points, points+dim, adt);
            nf = adt_full_search(&indices, points, points+dim, adt);
            if (n != nf)
            {
                fprintf(STDERR, "!!! Error: n=%i != nf=%i\n", n, nf);
                exit(-1);
            }
        }
    }

    /* Check adt_search */
    if (0)
    {
        for (i = 0; i < nnodes; i++)
        {
            n = adt_search(&indices, points+i*dim, points+i*dim, adt);
            if (n > 0 && 0)
            {
                printf("i=%i n=%i :", i, n);
                for (ii = 0; ii < n; ii++)
                    printf(" %i", indices[ii]);
                printf("\n");
            }
            nf = adt_full_search(&indices, points+i*dim, points+i*dim, adt);
            if (nf > 0 && 0)
            {
                printf("i=%i n=%i :", i, nf);
                for (ii = 0; ii < nf; ii++)
                    printf(" %i", indices[ii]);
                printf("\n");
            }
            if (n != nf)
            {
                fprintf(STDERR, "!!! Error: n=%i != nf=%i\n", n, nf);
                exit(-1);
            }
        }
    }

    /* Timing of adt_search */
    if (1)
    {
        t = clock();
        for (i = 0; i < nnodes; i++)
            n = adt_search(&indices, points+i*dim, points+i*dim, adt);
        t = (clock() - t) / CLOCKS_PER_SEC;
        printf("adt_search time = %g sec\n",t);
    }

    /* Timing of adt_full_search */
    if (0)
    {
        t = clock();
        for (i = 0; i < nnodes; i++)
            nf = adt_full_search(&indices, points+i*dim, points+i*dim, adt);
        t = (clock() - t) / CLOCKS_PER_SEC;
        printf("adt_full_search time = %g sec\n",t);
    }

    adt_free(adt);
}

/*---------------------------------------------------------------------------*
 *      Compute the region of the object, 
 *      i.e. coordinates of min and max points of the rectangle.
 */
REAL eps = 0.001;

void adt_region2(int i, REAL reg_min[2], REAL reg_max[2], void *user_data)
{
    REAL2 *points = (REAL2 *)user_data;

    reg_min[0] = points[i][0] - eps;
    reg_min[1] = points[i][1] - eps;
    reg_max[0] = points[i][0] + eps;
    reg_max[1] = points[i][1] + eps;
    return;
}

void adt_region3(int i, REAL reg_min[3], REAL reg_max[3], void *user_data)
{
    REAL3 *points = (REAL3 *)user_data;

    reg_min[0] = points[i][0] - eps;
    reg_min[1] = points[i][1] - eps;
    reg_min[2] = points[i][2] - eps;
    reg_max[0] = points[i][0] + eps;
    reg_max[1] = points[i][1] + eps;
    reg_max[2] = points[i][2] + eps;
    return;
}

#endif /* SELF_CHECK */

/*****************************************************************************/
/*                                                                           */
/*      Timing results of ADT search in 2D and 3D                            */
/*      for the set of points (x,y) and (x,y,z), respectively,               */
/*      where x,y,z = 1,...,M.                                               */
/*                                                                           */
/*      Pentium II 266 MHz workstation under Linux: cc -O3                   */
/*                                                                           */
/*        mesh    |  number   |    ADT    |   ADT  |   full  |               */
/*         MxM    | of nodes  | construct | search |  search | full/ADT      */
/*      ----------+-----------+-----------+--------+---------+----------     */
/*        25x25   |      625  |    0.00   |   0.01 |    0.06 |       6       */
/*        50x50   |     2500  |    0.01   |   0.03 |    1.02 |      34       */
/*       100x100  |    10000  |    0.03   |   0.11 |   34.23 |     311       */
/*       200x200  |    40000  |    0.16   |   0.49 |  766.24 |    1563       */
/*       400x400  |   160000  |    0.68   |   2.18 |    --   |   ~8000       */
/*       800x800  |   640000  |    2.82   |  10.31 |    --   |  ~40000       */
/*      1600x1600 |  2560000  |   13.61   |  81.42 |    --   | ~200000  swap */
/*                                                                           */
/*       mesh     |  number   |    ADT    |   ADT  |   full  |               */
/*      M x M     | of nodes  | construct | search |  search | full/ADT      */
/*    ------------+-----------+-----------+--------+---------+----------     */
/*       64x64    |     4096  |    0.02   |   0.04 |    3.69 |      92       */
/*      216x216   |    46656  |    0.19   |   0.57 | 1042.95 |    1829       */
/*      512x512   |   262144  |    1.13   |   3.68 |    --   |  ~10000       */
/*     1000x1000  |  1000000  |    4.45   |  15.17 |    --   | ~100000       */
/*                                                                           */
/*       mesh     |  number   |    ADT    |   ADT  |   full  |               */
/*     M x M x M  | of nodes  | construct | search |  search | full/ADT      */
/*    ------------+-----------+-----------+--------+---------+----------     */
/*     16x16x16   |     4096  |    0.02   |   0.04 |    3.22 |      80       */
/*     36x36x36   |    46656  |    0.21   |   0.61 | 1055.06 |    1729       */
/*     64x64x64   |   262144  |    1.31   |   4.03 |    --   |  ~10000       */
/*    100x100x100 |  1000000  |    5.39   |  16.69 |    --   | ~100000       */
/*                                                                           */
/*****************************************************************************/
