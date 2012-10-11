/*****************************************************************************/
/*                                                                           */
/* adt.h                                                                     */
/*                                                                           */
/* Include file for programs that operate with Alternating Digital Tree      */
/*                                                                           */
/* Version 14.06.01                                                          */
/*                                                                           */
/* List of functions that operate with Alternating Digital Tree (ADT):       */
/*                                                                           */
/*      adt_construct   - construct ADT (with no mask)                       */
/*      adt_init        - construct ADT (with mask)                          */
/*      adt_search      - search in ADT                                      */
/*      adt_free        - free ADT                                           */
/*                                                                           */
/*****************************************************************************/

#ifndef ADT_H_
#define ADT_H_

#ifndef REAL
#define REAL double
#endif

void *adt_construct             /* construct ADT */
        (const int nnodes,      /* number of nodes in the tree */
         const int dim,         /* dimension of the user's space: 2 or 3 */
         void *user_data,       /* pointer to the user's data */
         void (*adt_region)     /* interface to user's function */
                (int i,                 /* index of the node (on entry) */
                 REAL reg_min[],        /* min corner coordinates (on return) */
                 REAL reg_max[],        /* max corner coordinates (on return) */
                 void *user_data));     /* user's data (on entry) */

void *adt_init                  /* construct ADT */
        (int nnodes,            /* total number of nodes */
         const int dim,         /* dimension of the user's space: 2 or 3 */
         int *mask,             /* pointer to mask array or NULL */
         void *user_data,       /* pointer to the user's data */
         void (*adt_region)     /* interface to user's function */
                (int i,                 /* index of the node (on entry) */
                 REAL reg_min[],        /* min corner coordinates (on return) */
                 REAL reg_max[],        /* max corner coordinates (on return) */
                 void *user_data));     /* user's data (on entry) */

int adt_search                  /* seach in ADT */
        (int **indices,         /* indices list for nodes found (on return) */
         const REAL reg_min[],  /* min corner coordinates (on entry) */
         const REAL reg_max[],  /* max corner coordinates (on entry) */
         const void *adt);      /* pointer to ADT */

void adt_free                   /* free ADT */
        (const void *adt);      /* pointer to ADT */

#endif /* ADT_H_ */
