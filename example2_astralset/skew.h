/* code from http://graphics.lcs.mit.edu/~seth/geomlib/linelinecp.c */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define EPS (1.0E-7F)
#define MINLEN (EPS)
#define MAXRECIPLEN (1.0F / MINLEN)

/* Global Variables */

typedef struct pointstruct
{
   float px;
   float py;
   float pz;
} POINT;

/* plane ax + by + cz + d == 0  */
typedef struct planestruct
{
   float a;
   float b;
   float c;
   float d;
} PLANE;
 
typedef struct vecstruct
{
   float dx;
   float dy;
   float dz;
} VECTOR;


/* Prototypes */

float edist(POINT *a,POINT *b);
POINT* pplusv(POINT *a,POINT *p,VECTOR *v);
VECTOR* vcross(VECTOR *c,VECTOR *a,VECTOR *b);
float vdot(VECTOR *a,VECTOR *b);
float vmag(VECTOR *a);
float recip_vmag(VECTOR *a);
VECTOR* vnorm(VECTOR *a);
float point_plane_dist(POINT *a,PLANE *P);
POINT* plerp(POINT *p,POINT *a, POINT *b,float t);
POINT* intersect_line_plane(POINT *p,POINT *a, POINT *b,PLANE *M);
POINT* intersect_dline_plane(POINT *p,POINT *a, VECTOR *adir,PLANE *M);
PLANE* plane_from_two_vectors_and_point(PLANE *M,VECTOR *u,VECTOR *v,POINT *p);
int line_line_closest_points3d(POINT *pA, POINT *pB,POINT *a, VECTOR *adir,POINT *b, VECTOR *bdir);

/* Functions */

/* Function to get distance bewteen two points */
float edist(POINT *a,POINT *b)
{
   float rval;

   rval=sqrt((b->px - a->px)*(b->px - a->px)+(b->py - a->py)*(b->py - a->py)+(b->pz - a->pz)*(b->pz - a->pz));
   return (rval);
}

/* Function to get second point from a vector and first point (a = p + v) */
POINT* pplusv(POINT *a,POINT *p,VECTOR *v)
{
    a->px = p->px + v->dx;
    a->py = p->py + v->dy;
    a->pz = p->pz + v->dz;
    return (a);
}

/*  Function to get the cross-product of two vectors (c = a x b) */
VECTOR* vcross(VECTOR *c,VECTOR *a,VECTOR *b)
{
    c->dx = a->dy * b->dz - a->dz * b->dy;
    c->dy = a->dz * b->dx - a->dx * b->dz;
    c->dz = a->dx * b->dy - a->dy * b->dx;
    return (c);
}

/* Function to get the dot-product of two vectors */
float vdot(VECTOR *a,VECTOR *b)
{
    return (a->dx * b->dx + a->dy * b->dy + a->dz * b->dz);
}

/* Function to get the magnitude of a vector */
float vmag(VECTOR *a)
{
    return (sqrt(vdot(a, a)));
}

/* Function to get the reciprocal of the magnitude of a vector */
float recip_vmag(VECTOR *a)
{
    return (1.0 / sqrt(vdot(a, a)));
}

/* Function to normalise a vector (a /= ||a||) */
VECTOR* vnorm(VECTOR *a)
{
   float d;

   if ((d = recip_vmag(a)) > MAXRECIPLEN)
   {
     /* printf ("vector at 0x%x: %g %g %g?",(int)a, a->dx, a->dy, a->dz);*/

      a->dx = a->dy = a->dz = 0.0;

      return ('\0');
   }
   else
   {
      a->dx *= d;
      a->dy *= d;
      a->dz *= d;
      return (a);
   }
}

/* Function to get the distance between a point and a plane */
float point_plane_dist(POINT *a,PLANE *P)
{
   float rval;

   rval = (a->px * P->a + a->py * P->b + a->pz * P->c + P->d);

   return (rval);
}

/* p = (1-t)a + tb */
POINT* plerp(POINT *p,POINT *a, POINT *b,float t)
{
    p->px = t * b->px + (1-t) * a->px;
    p->py = t * b->py + (1-t) * a->py;
    p->pz = t * b->pz + (1-t) * a->pz;
    return (p);
}

/* Function to intersect a line and a plane and return the intersection point */
POINT* intersect_line_plane(POINT *p,POINT *a, POINT *b,PLANE *M)
{
   float Mdota, Mdotb, denom, t;
   int ip,iq;

   Mdota = point_plane_dist (a, M);
   Mdotb = point_plane_dist (b, M);

   denom = Mdota - Mdotb;

   if(fabs(denom / (fabs(Mdota) + fabs(Mdotb))) < EPS)
   {
      printf("int_line_plane(): no intersection?\n");
      p->px = p->py = p->pz = 0.0;
      return (NULL);
   }
   else
   {
      t = Mdota / denom;
      plerp (p, a, b, t);
      return (p);
   }
}

/* <point,direction> line form, return the intersection point from an intersected line and the plane in question */
POINT* intersect_dline_plane(POINT *p,POINT *a, VECTOR *adir,PLANE *M)               
{
   static POINT B, *b = &B;

   pplusv (b, a, adir);
   return (intersect_line_plane (p, a, b, M));
}

PLANE* plane_from_two_vectors_and_point(PLANE *M,VECTOR *u,VECTOR *v,POINT *p)
{
   vnorm (vcross ((VECTOR *)M, u, v));

   /* plane contains p */
   M->d = -(M->a * p->px + M->b * p->py + M->c * p->pz);
   return (M);
}

/* Function that computes a point on line B (pB) closest to line A  */
/* and then computes a point (pA) on line A that is closest to line B */
/* Returns: 0 if parallel; 1 if lines intersect (coincident); 2 if lines are skew (as expected) */
int line_line_closest_points3d(POINT *pA, POINT *pB,POINT *a, VECTOR *adir,POINT *b, VECTOR *bdir)
{
   static VECTOR Cdir, *cdir = &Cdir;
   static PLANE Ac, *ac = &Ac, Bc, *bc = &Bc;

   /* connecting line is perpendicular to both */
   vcross (cdir, adir, bdir);

   /* lines are near-parallel -- all points are closest */
   if(!vnorm(cdir))
   {
      *pA = *a;
      *pB = *b;
      return (0);
   }

   /* form plane containing line A, parallel to cdir  */
   plane_from_two_vectors_and_point (ac, cdir, adir, a);

   /* form plane containing line B, parallel to cdir  */
   plane_from_two_vectors_and_point (bc, cdir, bdir, b);

   /* closest point on A is line A ^ bc */
   intersect_dline_plane (pA, a, adir, bc);

   /* closest point on B is line B ^ ac  */
   intersect_dline_plane (pB, b, bdir, ac);

   if(edist (pA, pB) < EPS) return (1);     /* coincident */
   else return (2);                          /* distinct */
}



