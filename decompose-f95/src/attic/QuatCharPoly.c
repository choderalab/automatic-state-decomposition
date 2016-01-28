/*******************************************************************************
 * -/_|:|_|_\-
 *
 * File: QuatCharPoly.c
 *
 * Function: Rapid calculation of RMSD using a quaternion-based
 * characteristic polynomial
 *
 * Author(s): Douglas Theobald
 * Department of Chemistry and Biochemistry
 * UCB 215
 * University of Colorado at Boulder
 * Boulder, CO 80309-0215
 *
 * theobal@colorado.edu
 * dtheobald@hotmail.com
 *
 * Copyright: Copyright (c) 2005 Douglas L. Theobald
 *
 * QuatCharPoly.c is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * QuatCharPoly.c is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with theseus.c in the file 'COPYING'; if not, write to the:
 *
 * Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307 USA
 *
 * Source: started anew.
 *
 * Change History:
 * 5/5/05 6:21 PM Started source
 * 7/18/05 11:22 AM Fixed errors in calls to CoordsInnerProd()
 * Removed undefined mysquare() call in CalcQuarticCoeffs()
 * (thanks Chris Pettitt!)
 * ******************************************************************************/ 

/* This method is highly derived from: Horn, B. K. P. (1987). 
   "Closed-form solution of absolute orientation using unit quaternions." 
   J Opt Soc Am A 4(4):629Ð642. 
   
   Here the superposition problem is solved using a simple and very numerically stable 
   Newton-Raphson procedure. The minimizing rotation is ignored, only the minimum RMSD 
   is calculated. The main function is QuatCharPoly(), which takes as arguments: 

   const double **coords1 
     A 3 x N array of structure coordinates. 
   const double **coords2 
     The target structure coordinates. 
   const int len 
     The length N of the coords (# of atoms) 
   double *coeff 
     A pointer to an array of 3 doubles, to hold the last three quartic coefficients 

   QuatCharPoly() returns the sum of squared deviations at the least squares minimum 
   for coords1 and coords2. It *does not* return the RMSD. For the example provided 
   below, the minimum least-squares RMSD for the two 7-atom fragments should be 0.719106 A. 
   
   NB #1: QuatCharPoly() returns the sum of squared deviations (sumdev^2) between the 
   two structures in coords1 and coords2. RMSD = sqrt(sumdev^2 / N). 

   NB #2: If you are doing a full superposition (the usual least squares way), you MUST 
   center each structure first. That is, you must translate each structure so that its 
   centroid is at the origin. You can use CenterCoords() for this. 

   NB #3: Please note how I store structure coordinates in the double **coords arrays. 
   They are 3xN arrays, not Nx3 arrays as is also commonly used (where the x, y, z axes 
   are interleaved). The difference is something like this for storage of a structure 
   with 8 atoms: 

   Nx3: xyzxyzxyzxyzxyzxyzxyzxyz 
   3xN: xxxxxxxxyyyyyyyyzzzzzzzz 

   The functions can be easily modified, however, to accomodate any data format preference. 
   I chose this format because it is readily used in vectorized functions (SIMD, Altivec, 
   MMX, SSE2, etc.). 

   If you use this QCP RMSD calculation method in a publication, please reference: 

   Douglas L. Theobald (2005) 
   "Rapid calculation of RMSD using a quaternion-based characteristic polynomial." 
   Acta Crystallographica A 61(4):478-480. 
*/ 

/* gcc -O3 -Wall -Werror -ansi -o QuatCharPoly QuatCharPoly.c */ 

#include 
#include 
#include 

double QuatCharPoly(const double **coords1, const double **coords2, const int len, double *coeff); 
static void CalcQuarticCoeffs(const double **coords1, const double **coords2, const int len, double *coeff); 
static double QCProot(double *coeff, double guess, const double delta); 
static double eval_horn_NR_corrxn(const double *c, const double x); 
static double CoordsInnerProd(const double **coords, const int len); 
static void CenterCoords(double **coords, const int len); 
double **MatInit(const int rows, const int cols); 
void MatDestroy(double **matrix); 
static void PrintCoords(const double **coords, const int len); 

int main() { 
  double coeff[3], sumdev2, rmsd; 
  double **frag_a, **frag_b; 
  int len = 7; 

  frag_a = MatInit(3, len); 
  frag_b = MatInit(3, len); 

  frag_a[0][0] = -2.803; frag_a[1][0] = -15.373; frag_a[2][0] = 24.556; frag_a[0][1] = 0.893; frag_a[1][1] = -16.062; frag_a[2][1] = 25.147; frag_a[0][2] = 1.368; frag_a[1][2] = -12.371; frag_a[2][2] = 25.885; frag_a[0][3] = -1.651; frag_a[1][3] = -12.153; frag_a[2][3] = 28.177; frag_a[0][4] = -0.440; frag_a[1][4] = -15.218; frag_a[2][4] = 30.068; frag_a[0][5] = 2.551; frag_a[1][5] = -13.273; frag_a[2][5] = 31.372; frag_a[0][6] = 0.105; frag_a[1][6] = -11.330; frag_a[2][6] = 33.567; frag_b[0][0] = -14.739; frag_b[1][0] = -18.673; frag_b[2][0] = 15.040; frag_b[0][1] = -12.473; frag_b[1][1] = -15.810; frag_b[2][1] = 16.074; frag_b[0][2] = -14.802; frag_b[1][2] = -13.307; frag_b[2][2] = 14.408; frag_b[0][3] = -17.782; frag_b[1][3] = -14.852; frag_b[2][3] = 16.171; frag_b[0][4] = -16.124; frag_b[1][4] = -14.617; frag_b[2][4] = 19.584; frag_b[0][5] = -15.029; frag_b[1][5] = -11.037; frag_b[2][5] = 18.902; frag_b[0][6] = -18.577; frag_b[1][6] = -10.001; frag_b[2][6] = 17.996; 

  PrintCoords((const double **) frag_a, len); 
  PrintCoords((const double **) frag_b, len); 

  CenterCoords(frag_a, len); 
  CenterCoords(frag_b, len); 

  PrintCoords((const double **) frag_a, len); 
  PrintCoords((const double **) frag_b, len); 

  sumdev2 = QuatCharPoly((const double **) frag_a, (const double **) frag_b, len, coeff); 

  rmsd = sqrt(sumdev2 / len); 

  printf("\nsumdev^2 = %f\nrmsd = %f\n", sumdev2, rmsd); 

  MatDestroy(frag_a); 
  MatDestroy(frag_b); 

  exit(EXIT_SUCCESS); 
} 

/* returns the sum of the squared deviations (sumdev^2) between the two structures rmsd = sqrt(sumdev^2 / atom_num) */ 
double QuatCharPoly(const double **coords1, const double **coords2, const int len, double *coeff) { 
  double innerprod; 
  double lambdamax; 

  innerprod = CoordsInnerProd(coords1, len) + CoordsInnerProd(coords2, len); 
  CalcQuarticCoeffs(coords1, coords2, len, coeff); 
  lambdamax = QCProot(coeff, 0.5 * innerprod, 1e-6); 
  return (innerprod - (2.0 * lambdamax)); 
} 

/* A lot of register variables, but this sort of thing scales very well with new and improved processors */ 
static void CalcQuarticCoeffs(const double **coords1, const double **coords2, const int len, double *coeff)  { 
  double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz; 
  double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2, SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2, SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy, SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy; double x1, x2, y1, y2, z1, z2; 
  const double *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2], *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2]; 
  int i; 
  
  Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0; 

  for (i = 0; i < len; ++i) { 
    x1 = fx1[i]; y1 = fy1[i]; z1 = fz1[i]; 
    x2 = fx2[i]; y2 = fy2[i]; z2 = fz2[i]; 

    Sxx += (x1 * x2); 
    Sxy += (x1 * y2); 
    Sxz += (x1 * z2); 
    Syx += (y1 * x2); 
    Syy += (y1 * y2); 
    Syz += (y1 * z2); 
    Szx += (z1 * x2); 
    Szy += (z1 * y2); 
    Szz += (z1 * z2); 
  } 

  Sxx2 = Sxx * Sxx; 
  Syy2 = Syy * Syy; 
  Szz2 = Szz * Szz; 
  Sxy2 = Sxy * Sxy; 
  Syz2 = Syz * Syz; 
  Sxz2 = Sxz * Sxz; 
  Syx2 = Syx * Syx; 
  Szy2 = Szy * Szy; 
  Szx2 = Szx * Szx; 

  SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz); 
  Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2; 

  /* coeff[4] = 1.0; */ 
  /* coeff[3] = 0.0; */ 
  coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2); 
  coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz); 
  SxzpSzx = Sxz+Szx; 
  SyzpSzy = Syz+Szy; 
  SxypSyx = Sxy+Syx; 
  SyzmSzy = Syz-Szy; 
  SxzmSzx = Sxz-Szx; 
  SxymSyx = Sxy-Syx; 
  SxxpSyy = Sxx+Syy; 
  SxxmSyy = Sxx-Syy; 

  Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2; 

  coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz)) + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz)) + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz)) + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz)); 
} 

/* Newton-Raphson root finding */ 
static double QCProot(double *coeff, double guess, const double delta) { 
  int i; 
  double oldg; 

  for (i = 0; i < 50; ++i) { 
    oldg = guess; 
    /* guess -= (eval_horn_quart(coeff, guess) / eval_horn_quart_deriv(coeff, guess)); */ 
    guess -= eval_horn_NR_corrxn(coeff, guess); 
    
    if (fabs(guess - oldg) < fabs(delta*guess)) return(guess); 
  } 

  fprintf(stderr, "\n\n ERROR21: Newton-Raphson root-finding in \'QCProot()\' did not converge \n"); 
  exit(EXIT_FAILURE); 
} 

/* Evaluates the Newton-Raphson correction for the Horn quartic. only 11 FLOPs */ 
static double eval_horn_NR_corrxn(const double *c, const double x) { 
  double x2 = x*x; 
  double b = (x2 + c[2])*x; 
  double a = b + c[1]; 
  return((a*x + c[0])/(2.0*x2*x + b + a)); 
} 

/* Evaluates the Horn quartic for coefficients c and given x. */ 
double eval_horn_quart(const double *c, const double x) { 
  return(((x*x + c[2])*x + c[1])*x + c[0]); 
} 

/* Evaluates the derivative of the Horn quartic for coefficients c and given x. */ double eval_horn_quart_deriv(const double *c, const double x) { return(2.0*(2.0*x*x + c[2])*x + c[1]); } 

/* Calculate the inner product of some coordinates. This is the same as the squared radius of gyration without normalization for the number of atoms. */ 
static double CoordsInnerProd(const double **coords, const int len) { 
  int i; 
  double sum, tmpx, tmpy, tmpz; 
  const double *x = coords[0], *y = coords[1], *z = coords[2]; 

  sum = 0.0; 

  for (i = 0; i < len; ++i) { 
    tmpx = x[i]; tmpy = y[i]; tmpz = z[i]; 

    sum += (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz); 
  } 

  return(sum); 
} 

static void CenterCoords(double **coords, const int len) { 
  int i; 
  double xsum, ysum, zsum; 
  double *x = coords[0], *y = coords[1], *z = coords[2]; 

  xsum = ysum = zsum = 0.0; 
  for (i = 0; i < len; ++i) { 
    xsum += x[i]; ysum += y[i]; zsum += z[i]; 
  } 

  xsum /= len; ysum /= len; zsum /= len; 

  for (i = 0; i < len; ++i) { 
    x[i] -= xsum; y[i] -= ysum; z[i] -= zsum; 
  } 
} 

double **MatInit(const int rows, const int cols) { 
  int i; 
  double **matrix = NULL; 
  double *matspace = NULL; 
  matspace = (double *) calloc((rows * cols), sizeof(double)); 
  if (matspace == NULL) { 
    perror("\n ERROR"); 
    puts("\n ERROR: Failure to allocate room for pointers"); 
    exit(EXIT_FAILURE); 
  } 

  /* allocate room for the pointers to the rows */ 
  matrix = (double **) malloc(rows * sizeof(double *)); 
  if (matrix == NULL) { 
    perror("\n ERROR"); 
    puts("\n ERROR: Failure to allocate room for pointers"); 
    exit(EXIT_FAILURE); 
  } 

  /* now 'point' the pointers */ 
  for (i = 0; i < rows; i++) 
    matrix[i] = matspace + (i * cols); 
  return(matrix); 
} 

void MatDestroy(double **matrix) { 
  if (matrix[0] != NULL) 
    free(matrix[0]); 

  if (matrix != NULL) 
    free(matrix); 

  matrix = NULL; 
} 

static void PrintCoords(const double **coords, const int len) { 
  int i; 
  for (i = 0; i < len; ++i) printf("\n % 8.3f % 8.3f % 8.3f", coords[0][i], coords[1][i], coords[2][i]); putchar('\n'); 
}

