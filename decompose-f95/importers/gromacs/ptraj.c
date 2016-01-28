/*  _______________________________________________________________________
 *
 *                        RDPARM/PTRAJ: 2006
 *  _______________________________________________________________________
 *
 *  This file is part of rdparm/ptraj.
 *
 *  rdparm/ptraj is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  rdparm/ptraj is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You can receive a copy of the GNU General Public License from
 *  http://www.gnu.org or by writing to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *  ________________________________________________________________________
 *
 *  AMBER CVS tracking:
 *
 *  $Header: /share_kubo/Storage/cvsroot/decompose-f95/importers/gromacs/ptraj.c,v 1.1 2006/09/27 18:50:21 jchodera Exp $
 *
 *  Revision: $Revision: 1.1 $
 *  Date: $Date: 2006/09/27 18:50:21 $
 *  Last checked in by $Author: jchodera $
 *  ________________________________________________________________________
 *
 *
 *  CONTACT INFO: To learn who the code developers are, who to contact for
 *  more information, and to know what version of the code this is, refer
 *  to the CVS information and the include files (contributors.h && version.h)
 *
 */

#include "contributors.h"
#include "version.h"

/*  ________________________________________________________________________
 */


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define PTRAJ_MODULE
#include "ptraj.h"


/*
 *  This is the main routine for the ptraj functionality and in general this
 *  code is likely not to be modified by general users/developers, except:
 *
 *  (*) Adding new "actions" or routines to process/analyze coordinate sets.
 *      In this case, ptrajSetup() will likely need modification for any
 *      "action" on the coordinates or ptrajSetupAnalyze() for any routine that
 *      proprocesses data accumulated by the actions.
 *  (*) Adding new coordinate formats.  Modify checkCoordinates(), 
 *      ptrajSetupIO(), printCoordinateInfo(), ptrajPreProcessInputCoordinates(),
 *      and ptrajProcessInputCoordinates()
 *  (*) Modifying the code to parse "mask" strings: all the *mask* routines
 *
 *  For more information on adding new actions, see the details comments in
 *  actions.c.
 */


   int
atomToResidue(int atom, int residues, int *ipres)
{
  int i;

  if (ipres == NULL) return -1;

  for (i = 0; i < residues; i++)
    if (atom >= ipres[i] && atom < ipres[i+1])
      return (i+1);

  return -1;
}



   int
atomToSolventMolecule(int atom, int molecules, int *start, int *stop)
{
  int i;

  if (start == NULL || stop == NULL) return -1;

  for (i = 0; i < molecules; i++)
    if (atom <= start[i])
      return -1;
    else if (atom > start[i] && atom <= stop[i])
      return (i+1);

  return -1;
}



   int
atomToMolecule(int atom, int molecules, int *mol)
{
  int i, a;

  if (mol == NULL) return -1;

  a = 0;
  for (i = 0; i < molecules; i++) {
    a += mol[i];
    if (atom <= a) 
      return (i+1);
  }
  return -1;
}






   int
isActiveDetailed(int atom, int residue, int *mask, int atoms, int residues,
		 Name *atomName, Name *residueName, int *ipres)
{
  int i;
 
  if (residue >= residues || residue < 0) {
    printf("WARNING: residue out of range in isActiveDetailed, res %i (total %i)\n",
	   residue, residues);
    return 0;
  }
  for (i = ipres[residue]-1; i < ipres[residue+1]-1; i++)
    if ( mask[i] && strcmp(atomName[i], atomName[atom]) == 0 ) 
      return 1;

  return 0;

}


   int
isActive(int atom, int residue, int *mask, ptrajState *state)
{
  return( isActiveDetailed(atom, residue, mask,
			   state->atoms, 
			   state->residues,
			   state->atomName,
			   state->residueName,
			   state->ipres) );
}



   int
isActiveResidueDetailed(int residue, int *mask, int atoms, int residues,
			Name *atomName, Name *residueName, int *ipres)
{
  int i, total;

  total = 0;
  if (residue >= residues || residue < 0) {
    printf("WARNING: residue out of range in isActiveResidueDetailed(), res %i (total %i)\n",
	   residue, residues);
    return 0;
  }
  for (i = ipres[residue]-1; i < ipres[residue+1]-1; i++)
    if ( ! mask[i] )
      return 0;
    else
      total++;

  return total;

}


   int
isActiveResidue(int residue, int *mask, ptrajState *state)
{
  return( isActiveResidueDetailed(residue, mask,
				  state->atoms, 
				  state->residues,
				  state->atomName,
				  state->residueName,
				  state->ipres) );
}



   void
printAtomMaskDetailed(int *mask, int atoms, int residues,
		      Name *atomName, Name *residueName, int *ipres)
{
  int i, j, curres;
  char tmpatom[20];
  int *resactive, *ressimilar;
  int printed, numactive, numresactive, numressimilar;
  int incurres, innextres;

  printed = 0;
  numactive = 0;
  numresactive = 0;
  numressimilar = 0;

  /*
   *  This routine is kind of junky since in general I want to avoid printing
   *  as much detail as possible.  Therefore, we check to see is certain ranges 
   *  residues have all atoms active, etc. to avoid printing each atom in a residue.
   *  This makes it ugly and obtuse.
   */


  if (mask == NULL) {
    printf("[No atoms are selected]\n");
    return;
  }

  j=0;
  for (i=0; i < atoms; i++)
    if (mask[i]) j++;

  if (j == 0) {
    printf("[No atoms are selected]\n");
    return;
  }

     /*
      *  check if all atoms are active and if so print an asterisk
      */

  j = 0;
  for (i=0; i < atoms; i++) {
    j += mask[i];
  }
  if ( j == atoms ) {
    fprintf(stdout, "  * (All atoms are selected)\n");
    return;
  }
  numactive = j;

     /*
      *  determine which residues have all the atoms in that residue active
      */

  resactive = (int *) safe_malloc(sizeof(int) * residues);
  for (i=0; i < residues; i++) {
    resactive[i] = 0.0;
  }

  curres = 0;
  j = 0;
  for (i=0; i < atoms; i++) {
    if (i == ipres[curres+1]-1) {
      if (j == ipres[curres+1] - ipres[curres]) {
	resactive[curres] = 1.0;
	numresactive++;
      }
      j = 0;
      curres++;
    }
    if (mask[i])
      j++;
  }
  if (j == ipres[curres+1] - ipres[curres]) {
    resactive[curres] = 1.0;
  }

     /*
      *  determine the range over which the residues are fully active
      */

  for (curres = residues-2; curres >= 0; curres--) {
    if (resactive[curres]) {
      resactive[curres] += resactive[curres+1];
      numresactive--;
    }
  }


  /*
   *  determine ranges over which residues have the same atoms active
   *  as the next residue
   */
  ressimilar = (int *) safe_malloc(sizeof(int) * residues);
  for (i=0; i < residues; i++) {
    ressimilar[i] = 0.0;
  }

  for (curres = residues-2; curres >=0; curres--) {

    incurres = 0;
    innextres = 0;
    for (i = ipres[curres]-1; i < ipres[curres+1]-1; i++) {
      if ( mask[i] ) {
	incurres++;
      }

      if (isActiveDetailed(i, curres+1, mask,
			   atoms, residues, atomName, residueName, ipres))
	innextres++;
    }
    if (incurres && innextres == incurres) {
      ressimilar[curres] = ressimilar[curres+1] + 1;
    } else {
      numressimilar++;
    }

  }

   
     /*
      *  do the actual printing
      */

  j = 0;
  for (curres = 0; curres < residues; curres++) {

    if (resactive[curres] ) {

      /*
       *  If all of the atoms are active in this residue print either the
       *  residue number or range as appropriate
       */

      if (resactive[curres] > 2) {
	if (j!=0 && j%10 != 0) fprintf(stdout, ",");
	fprintf(stdout, ":%i-%i", curres+1, curres+resactive[curres]);
	curres += resactive[curres]-1;
      } else {
	if (j!=0 && j%10 != 0) fprintf(stdout, ",");
	fprintf(stdout, ":%i", curres+1);
      }
      j++;
      if (j != 0 && j % 10 == 0) {
	fprintf(stdout, "\n    ");
	j = 0;
      }
    } else if (ressimilar[curres]) {

      /*
       *  If there is a set of residues with a similar atom selection...
       */

      if (ressimilar[curres] > 2) {
	if (j!=0 && j%10 != 0) fprintf(stdout, ",");
	fprintf(stdout, ":%i-%i", curres+1, curres+ressimilar[curres]+1);
	curres += ressimilar[curres];
      } else {
	if (j!=0 && j%10 != 0) fprintf(stdout, ",");
	fprintf(stdout, ":%i", curres+1);
      }

      for (i = ipres[curres]-1; i < ipres[curres+1]-1; i++) {
	if ( mask[i] ) {
	  if (printed) 
	    fprintf(stdout, ",");
	  else {
	    fprintf(stdout, "@");
	    printed = 1;
	  }
	  strcpy(tmpatom, atomName[i]);
	  fprintf(stdout, "%s", strtok(tmpatom, " "));
	}
      }
      j++;
      if (j != 0 && j % 10 == 0) {
	fprintf(stdout, "\n    ");
	j = 0;
      }


    } else {

      /*
       *  Print individual atoms
       */

      if (numactive > 10 && numressimilar > 10 && numresactive > 10) {
	fprintf(stdout, "\n    ");
	numactive = 0;
      }
      for (i = ipres[curres]-1; i < ipres[curres+1]-1; i++) {
	if ( mask[i] ) {
	  if (j!=0 && j%10 != 0) fprintf(stdout, ",");
	  
	  strcpy(tmpatom, atomName[i]);
	  fprintf(stdout, ":%i@%s", curres+1, strtok(tmpatom, " "));
	  j++;
	}
	if (j != 0 && j % 10 == 0) {
	  fprintf(stdout, "\n    ");
	  j = 0;
	}
      }
    }
  }
  if (j!=0 && j%10 != 0) fprintf(stdout, "\n");
  safe_free(resactive);
  safe_free(ressimilar);
}


   void
printAtom(FILE *fpout, int atom, ptrajState *state)
{
  int curres;
  char buffer[50];


  curres = atomToResidue(atom+1, state->residues, state->ipres)-1;
  sprintf(buffer, ":%i               ", curres+1);
  sprintf(buffer+6, "%s %s", state->atomName[atom], state->residueName[curres]);
  fprintf(fpout, "%s", buffer);
}


   void
printAtomCompact(FILE *fpout, int atom, ptrajState *state)
{
  int curres;
  char buffer[50], *bufferp;

  curres = atomToResidue(atom+1, state->residues, state->ipres)-1;

  sprintf(buffer, ":%i", curres+1);
  bufferp = buffer+strlen(buffer);
  sprintf(bufferp, "@%s", state->atomName[atom]);

  fprintf(fpout, "%10s", buffer);


}







   parseEntry *
parseToken(char **textp, int operator)
{
  char *text;
  int i, j;
  parseEntry *p;

  text = *textp;

  for (i=0;; i++) {
    if (text[i] == (char) 0 ||
	text[i] == ',' ||
	text[i] == ':' ||
	text[i] == '@' ||
	(i>0 && text[i] == '-' && !isalpha(text[i-1])) ||
	isspace(text[i]))
      break;
  }

  p = (parseEntry *) safe_malloc(sizeof(parseEntry));
  p->isnumber = 0;
  
  if ( i > 0 ) 
    p->token = (char *) safe_malloc(sizeof(char) * (i+1));

  for (j=0; j < i; j++) {
    p->token[j] = text[j];
  }
  p->token[i] = (char) 0;


  p->operator = PARSE_NOOP;
  switch ( text[i] ) {

  case '-':
    p->operator = PARSE_CONTINUATION;

  case ',':

    text++;

  }

  if ( isdigit( p->token[0] ) ) {
    p->isnumber = 1;
    for (j=1; j < i; j++) {
      if ( isdigit(p->token[j]) == 0 ) {
	p->isnumber = 0;
      }
    }
  }

  text = text+i;
  skipWhitespace(text);
  /*
   *  extra special check to handle extra spacing between continuation
   *  operators
   */
  if (text[0] == '-') {
    p->operator = PARSE_CONTINUATION;
    text += 1;
    skipWhitespace(text);
  } else if (text[0] == ',') {
    text += 1;
    skipWhitespace(text);
  }

  *textp = text;
  return(p);


}


   int
isMatch(char *s1, char *s2) 
{
  int i;

  /*
   *  straight match
   */

  if ( strcmp(s1, s2) == 0 ) return 1;

  /*
   *  fuzzy match: NOTE this will break if s1 > s2
   */

  if ( strchr(s1, '?') != NULL || strchr(s1, '*') != NULL ) {

    /*
     *  look over the minimal map between the two strings
     */
    for (i=0; i < strlen(s1) && i < strlen(s2); i++) {

      if ( s1[i] != s2[i] ) {
	switch( s1[i] ) {

	case '*':           /* wild card, multiple characters */
	  return 1;
	case '?':           /* wild card, single character    */
	  if ( isspace(s2[i]) ) return 0;
	  break;
	default:            /* mismatch                       */
	  return 0;
	}
      }
    }
    return 1;
  }
  return 0;
}



   void
sortIndex(double *values, int *index, int size) 
{

#ifdef SHELL_SORT
  double *sortValue;
  double approx_ln_to_log2, ftmp;
  int logsplit, itmp, psort;
  int i, j, k;

  approx_ln_to_log2 = 1.0 / log( (double) 2.0 ) + 0.000001;

  sortValue = (double *) safe_malloc(sizeof(double) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  logsplit = ( log( (double) size ) * approx_ln_to_log2 );
  i = size;

  for (psort = 1; psort <= logsplit; psort++) {
    i >>= 1;
    for (j = i+1; j <= size; j++) {
      k = j - i;
      ftmp = sortValue[j-1];
      itmp = index[j-1];
      while (k >= 1 && sortValue[k-1] > ftmp ) {
	sortValue[k+i-1] = sortValue[k-1];
	index[k+i-1] = index[k-1];
	k -= i;
      }
      sortValue[k+i-1] = ftmp;
      index[k+i-1] = itmp;
    }
  }

  safe_free(sortValue);
#else
  double *sortValue;
  double ftmp;
  int itmp;
  int i,j;
    
  sortValue = (double *) safe_malloc(sizeof(double) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  for (i=1; i < size; i++) {
    ftmp = sortValue[i];
    itmp = index[i];
    j = i-1;
    while (j >= 0 && sortValue[j] > ftmp) {
      sortValue[j+1] = sortValue[j];
      index[j+1] = index[j];
      j--;
    }
    sortValue[j+1] = ftmp;
    index[j+1] = itmp;
  }
  safe_free(sortValue);
#endif
}


   void
sortIndexFloat(float *values, int *index, int size) 
{

#ifdef SHELL_SORT
  float *sortValue;
  float approx_ln_to_log2, ftmp;
  int logsplit, itmp, psort;
  int i, j, k;

  approx_ln_to_log2 = 1.0 / log( (double) 2.0 ) + 0.000001;

  sortValue = (float *) safe_malloc(sizeof(float) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  logsplit = ( log( (double) size ) * approx_ln_to_log2 );
  i = size;

  for (psort = 1; psort <= logsplit; psort++) {
    i >>= 1;
    for (j = i+1; j <= size; j++) {
      k = j - i;
      ftmp = sortValue[j-1];
      itmp = index[j-1];
      while (k >= 1 && sortValue[k-1] > ftmp ) {
	sortValue[k+i-1] = sortValue[k-1];
	index[k+i-1] = index[k-1];
	k -= i;
      }
      sortValue[k+i-1] = ftmp;
      index[k+i-1] = itmp;
    }
  }

  safe_free(sortValue);
#else
  float *sortValue;
  float ftmp;
  int itmp;
  int i,j;
    
  sortValue = (float *) safe_malloc(sizeof(float) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  for (i=1; i < size; i++) {
    ftmp = sortValue[i];
    itmp = index[i];
    j = i-1;
    while (j >= 0 && sortValue[j] > ftmp) {
      sortValue[j+1] = sortValue[j];
      index[j+1] = index[j];
      j--;
    }
    sortValue[j+1] = ftmp;
    index[j+1] = itmp;
  }
  safe_free(sortValue);
#endif
}





/*
 *  The goal of this routine is to create a new ptrajState (newstate)
 *  based on the old ptrajState (oldstate) deleting atoms that are
 *  set in the mask array.  This will modify the value of newstate.
 */

#undef  ROUTINE
#define ROUTINE "modifyStateByMask()"

   void
modifyStateByMask(ptrajState **newstatep, ptrajState **oldstatep,
		  int *mask, int strip)
{
  int i, ires, isol, imol;
  int j, jres, jsol, jmol;
  int curres, cursol, curmol; 
  int k;
  Name *atomName, *residueName;
  double *charges, *masses;
  int *startsol, *stopsol, *ipres, *moleculeInfo;
  ptrajState *newstate, *oldstate;

  oldstate = *oldstatep;

     /*
      *  allocate space for the new state
      */
  newstate = (ptrajState *) safe_malloc(sizeof(ptrajState));
  INITIALIZE_ptrajState(newstate);

     /*
      *  allocate space for temporary arrays and perform initialization
      */

  atomName     = (Name *)   safe_malloc(sizeof(Name)   * oldstate->atoms);
  charges      = (double *) safe_malloc(sizeof(double) * oldstate->atoms);
  masses       = (double *) safe_malloc(sizeof(double) * oldstate->atoms);
  residueName  = (Name *)   safe_malloc(sizeof(Name)   * oldstate->residues);
  ipres        = (int *)    safe_malloc(sizeof(int)    * (oldstate->residues+1));
  if (oldstate->molecules) 
    moleculeInfo = (int *)    safe_malloc(sizeof(int)    * oldstate->molecules);
  if (oldstate->solventMolecules) {
    startsol = (int *) safe_malloc(sizeof(int) * oldstate->solventMolecules);
    stopsol  = (int *) safe_malloc(sizeof(int) * oldstate->solventMolecules);

    for (i=0; i < oldstate->solventMolecules; i++)
      startsol[i] = -1;
  }
  
  j = 0; 
  jres = -1; jsol = -1; jmol = -1;
  ires = -1; isol = -1; imol = -1;

  /*
   *  loop over all atoms and set up information for the newstate if the atom is 
   *  not to be deleted...
   */
  for (i=0; i < oldstate->atoms; i++) {

    curres = atomToResidue(i+1, oldstate->residues, oldstate->ipres)-1;

    if ( mask[i] == (strip ? 0 : 1) ) {
      /*
       *  this atom is not to be deleted
       */

         /*
          *  copy over atom information
          */
      strcpy(atomName[j], oldstate->atomName[i]);
      charges[j] = oldstate->charges[i];
      masses[j] = oldstate->masses[i];

         /*
          *  check to see if we are in the same residue or not
          *  and copy relevant information
          */
      if (ires == -1 || ires != curres) {
	jres++;
	strcpy(residueName[jres], oldstate->residueName[curres]);
	ipres[jres] = j+1;
	ires = curres;
	
      }

         /*
          *  deal with the molecule information
          */
      if (oldstate->molecules) {
	curmol = atomToMolecule(i+1, oldstate->molecules, oldstate->moleculeInfo)-1;
	if (imol == -1 || imol != curmol) {
	  jmol++;
	  moleculeInfo[jmol] = oldstate->moleculeInfo[curmol];
	  for (k=1; k < oldstate->moleculeInfo[curmol]; k++)
	    if (mask[i+k] == (strip ? 1 : 0))
	      moleculeInfo[jmol]--;
	  imol = curmol;
	}
      }

         /*
          *  deal with the solvent information
          */
      if (oldstate->solventMolecules) {

	cursol = atomToSolventMolecule(i+1, oldstate->solventMolecules, 
				       oldstate->solventMoleculeStart,
				       oldstate->solventMoleculeStop)-1;

	if (cursol >= 0 && (isol == -1 || isol != cursol)) {
	  jsol++;
	  startsol[jsol] = j;
	  stopsol[jsol] = j; 
	  for (k=i; k < oldstate->solventMoleculeStop[cursol]; k++)
	    if (mask[k] == (strip ? 0 : 1)) stopsol[jsol]++;
	  isol = cursol;
	}
      }
         /*
          *  increment the new atom counter
          */
      j++;
    }
  }

     /*
      *  fix up IPRES
      */
  ipres[jres+1] = j+1;

     /*
      *  set up the newstate using the date placed into the temporary arrays;
      *  free up unneeded space as we go...
      */
  newstate->atoms = j;
  newstate->masses   = (double *) safe_malloc(sizeof(double) * newstate->atoms);
  newstate->charges  = (double *) safe_malloc(sizeof(double) * newstate->atoms);
  newstate->atomName = (Name *)   safe_malloc(sizeof(Name)   * newstate->atoms);

  for (i=0; i < newstate->atoms; i++) {
    newstate->masses[i]  = masses[i];
    newstate->charges[i] = charges[i];
    strcpy(newstate->atomName[i], atomName[i]);
  }
  safe_free(masses);
  safe_free(charges);
  safe_free(atomName);

  newstate->residues = jres+1;
  newstate->ipres = (int *) safe_malloc(sizeof(int) * (newstate->residues+1));
  newstate->residueName = (Name *) safe_malloc(sizeof(Name) * newstate->residues);
  for (i=0; i < newstate->residues; i++) {
    newstate->ipres[i] = ipres[i];
    strcpy(newstate->residueName[i], residueName[i]);
  }
  newstate->ipres[i] = ipres[i];

  safe_free(ipres);
  safe_free(residueName);

  newstate->IFBOX = oldstate->IFBOX;
  if (jsol < 0) {
    jsol = 0;
    newstate->solventMolecules = jsol;
    newstate->solventMoleculeStart = NULL;
    newstate->solventMoleculeStop = NULL;
    newstate->solventMask = NULL;
  } else {
    newstate->solventMolecules = jsol+1;
    newstate->solventMoleculeStart = (int *) safe_malloc(sizeof(int) * (jsol+1));
    newstate->solventMoleculeStop  = (int *) safe_malloc(sizeof(int) * (jsol+1));
    for (i=0; i < newstate->solventMolecules; i++) {
      newstate->solventMoleculeStart[i] = startsol[i];
      newstate->solventMoleculeStop[i]  = stopsol[i];
    }
    newstate->solventMask = (int *) safe_malloc(sizeof(int) * newstate->atoms);
    newstate->solventAtoms = 0;

    cursol = 0;
    for (i=0; i < newstate->atoms; i++) {
      newstate->solventMask[i] = 0;

      if (i == newstate->solventMoleculeStart[cursol]) {
	for (j=i; j < newstate->solventMoleculeStop[cursol]; j++) {
	  newstate->solventAtoms++;
	  newstate->solventMask[j] = 1;
	}
	i = newstate->solventMoleculeStop[cursol]-1;
	cursol++;
      }
    }
    safe_free(startsol);
    safe_free(stopsol);
  }

  newstate->maxFrames = oldstate->maxFrames;
  if (oldstate->molecules) {
    newstate->molecules = jmol+1;
    newstate->moleculeInfo = (int *) safe_malloc(sizeof(int) * newstate->molecules);
    for (i=0; i<newstate->molecules; i++)
      newstate->moleculeInfo[i] = moleculeInfo[i];
    safe_free(moleculeInfo);
  }

  for (i=0; i<6; i++)
    newstate->box[i] = oldstate->box[i];

  *newstatep = newstate;
}




   void
printAtomMask(int *mask, ptrajState *state)
{
  printAtomMaskDetailed(mask, 
			state->atoms, 
			state->residues,
			state->atomName,
			state->residueName,
			state->ipres);
}




   void
printHBondMask(int *mask, int *maskH1, int *maskH2, int *maskH3, ptrajState *state)
{
  int i;

  fprintf(stdout, "    Atom#  Residue#  Name   --   Atom#  Residue#  Name\n");
  for (i=0; i < state->atoms; i++) {

    if (mask[i] != 0) {
      
      fprintf(stdout, "  %5i    %4i       %s ",
	      i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
      
      if (maskH1 != NULL) {
	if (maskH1[i] >= 0) {
	  
	  fprintf(stdout, " -- %5i    %4i       %s ",
		  maskH1[i]+1, atomToResidue(maskH1[i]+1, state->residues, state->ipres),
		  state->atomName[maskH1[i]]);
	}
	if (maskH2[i] >= 0) {
	  
	  fprintf(stdout, "\n  %5i    %4i       %s ",
		  i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
	  fprintf(stdout, " -- %5i    %4i       %s ",
		  maskH2[i]+1, atomToResidue(maskH2[i]+1, state->residues, state->ipres),
		  state->atomName[maskH2[i]]);
	}
	
	if (maskH3[i] >= 0) {
	  fprintf(stdout, "\n  %5i    %4i       %s ",
		  i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
	  fprintf(stdout, " -- %5i    %4i      %s ",
		  maskH3[i]+1, atomToResidue(maskH3[i]+1, state->residues, state->ipres),
		  state->atomName[maskH3[i]]);
	}
      }
      fprintf(stdout, "\n");
    }
  }
}


#undef  ROUTINE
#define ROUTINE "printHBondInfo()"

   void
printHBondInfo(int numdonor, int *donor, int numacceptor, 
	       int *acceptor, int *acceptorH, ptrajState *state)
{
  int i;

  if (numdonor > 0 && donor != NULL) {
    fprintf(stdout, "  HBOND DONOR LIST:\n");
    fprintf(stdout, "    Atom#  Residue#  Name\n");
    for (i=0; i < numdonor; i++) {

      fprintf(stdout, "  %5i    %4i       %s\n",
	      donor[i]+1, atomToResidue(donor[i]+1, state->residues, state->ipres), 
	      state->atomName[donor[i]]);
    }
  }

  if (numacceptor > 0 && acceptor != NULL) {
    fprintf(stdout, "  HBOND ACCEPTOR LIST:\n");
    fprintf(stdout, "    Atom#  Residue#  Name   --   Atom#  Residue#  Name\n");
    for (i=0; i < numacceptor; i++) {

      fprintf(stdout, "  %5i    %4i       %s ",
	      acceptor[i]+1, atomToResidue(acceptor[i]+1, state->residues, state->ipres), 
	      state->atomName[acceptor[i]]);

      fprintf(stdout, " -- %5i    %4i       %s\n",
	      acceptorH[i]+1, atomToResidue(acceptorH[i]+1, state->residues, state->ipres),
	      state->atomName[acceptorH[i]]);
    }
  }
}


   int *
processAtomMaskDetailed( char *maskString, int atoms, int residues,
			Name *atomName, Name *residueName, int *ipres)
{
  int actualAtoms;
  int not = 0;
  char *maskp;
  char *tmp;

  int *residueMask;
  int residueMaskActive = 0;
  int *atomMask;
  int atomMaskActive = 0;
  int res1, res2;
  int atom1, atom2;
  int continuation = 0;
  int i, j;
  stackType *residueStack = NULL;
  stackType *atomStack = NULL;
  parseEntry *pp;

  Name name;

  maskp = maskString;
  skipWhitespace(maskp);

  /*
   *  allocate mask strings
   */
  atomMask = (int *) safe_malloc(sizeof(int) * atoms);
  residueMask = (int *) safe_malloc(sizeof(int) * residues);
  memset(atomMask, 0, sizeof(int) * atoms);
  memset(residueMask, 0, sizeof(int) * residues);

  /*
   *  special case, choose ALL atoms
   */
  if ( maskp[0] == (char) 0 || maskp[0] == '*' ) {
    for (i=0; i < atoms; i++)
      atomMask[i] = 1;
    goto clean_return;
  }


  /*
   *  get rid of all NOT characters "~"; only one is significant
   *  and set NOT status if one is found...
   */
  while ( (tmp = strchr(maskString, '~' )) != NULL ) {
    not = 1;
    tmp[0] = ' ';
  }

  /*
   *  check for error
   */
  if (strchr(maskp, ':') == NULL &&
      strchr(maskp, '@') == NULL) {
    fprintf(stdout, "WARNING: Error in mask string, no \"@\" or \":\" present (%s)\n", 
	    maskString);
    safe_free( atomMask );
    safe_free( residueMask );
    return NULL;
  }

  /*
   *  do the main "parsing"
   */
  skipWhitespace(maskp);
  while ( maskp[0] != (char) 0 ) {
    
    if ( maskp[0] == ':' ) {
      maskp++;
      for (;;) {

	skipWhitespace(maskp);
	pp = parseToken(&maskp, PARSE_RESIDUE);
	pushBottomStack(&residueStack, (void *) pp);
	if (maskp[0] == (char) 0 || 
	    maskp[0] == '@' || 
	    maskp[0] == ':') break;
      }
    }

    if ( maskp[0] == '@' ) {
      maskp++;
      for (;;) {

	skipWhitespace(maskp);
	pp = parseToken(&maskp, PARSE_ATOM);
	pushBottomStack(&atomStack, (void *) pp);
	if ( maskp[0] == (char) 0 || maskp[0] == ':' ) break;
	if ( maskp[0] == '@' )
	  maskp++;
      }
    }

    /*
     *  now process the atomStack and residueStack
     */


    if ( not ) {
      for (i=0; i < atoms; i++)
	atomMask[i] = 1;
      for (i=0; i < residues; i++)
	residueMask[i] = 1;
    }

    while ( residueStack != NULL ) {

      if ( continuation ) {
	res1 = res2;
	res2 = -1;
      }

      pp = (parseEntry *) popStack( &residueStack );
      if ( pp->isnumber ) {
	if ( sscanf(pp->token, "%i", &res2) != 1 ) {
	  fprintf(stdout, "WARNING: error parsing atom mask\n");
	  safe_free( atomMask );
	  safe_free( residueMask );
	  return NULL;
	}
	res2--;

	if (continuation) {
	  continuation = 0;
	  if (res1 < 0) res1 = 0;
	  if (res2 >= residues) res2 = residues-1;
	  if (res2 < res1) res2 = res1;

	  for (i = res1; i <= res2; i++) {
	    residueMask[i] = (not ? 0 : 1);
	    residueMaskActive = 1;
	  }
	} else {
	  residueMask[res2] = (not ? 0 : 1);
	  residueMaskActive = 1;
	}

	if ( pp->operator == PARSE_CONTINUATION )
	  continuation = 1;

      }	else {

	continuation = 0;
	strcpy(name, "    ");
	for (i=0; i < strlen(pp->token) && i < NAME_SIZE; i++) {
	  name[i] = pp->token[i];
	}
	name[NAME_SIZE-1] = (char) 0;

	for (i=0; i < residues; i++)
	  if ( isMatch(name, residueName[i]) ) {
	    residueMask[i] = (not ? 0 : 1);
	    residueMaskActive = 1;
	  } 
      }
      safe_free(pp->token);
      pp->token=NULL;
      safe_free(pp);
    }

    while ( atomStack != NULL ) {

      if ( continuation ) {
	atom1 = atom2;
	atom2 = -1;
      }

      pp = (parseEntry *) popStack( &atomStack );
      if ( pp->isnumber ) {
	if ( sscanf(pp->token, "%i", &atom2) != 1 ) {
	  fprintf(stdout, "WARNING: error parsing atom mask\n");
	  safe_free( atomMask );
	  safe_free( residueMask );
	  return NULL;
	}
	atom2--;

	if (continuation) {
	  continuation = 0;
	  if (atom1 < 0) atom1 = 0;
	  if (atom2 > atoms) atom2 = atoms-1;
	  if (atom2 < atom1) atom2 = atom1;

	  for (i = atom1; i <= atom2; i++) {
	    atomMask[i] = (not ? 0 : 1);
	    atomMaskActive = 1;
	  }
	} else {
	  atomMask[atom2] = (not ? 0 : 1);
	  atomMaskActive = 1;
	}

	if ( pp->operator == PARSE_CONTINUATION )
	  continuation = 1;

      }	else {

	continuation = 0;
	strcpy(name, "    ");
	for (i=0; i < strlen(pp->token) && i < NAME_SIZE; i++) {
	  name[i] = pp->token[i];
	}
	name[NAME_SIZE-1] = (char) 0;

	for (i=0; i < atoms; i++)
	  if ( isMatch(name, atomName[i]) ) {
	    atomMask[i] = (not ? 0 : 1);
	    atomMaskActive = 1;
	  }
      }
      safe_free(pp->token);
      pp->token = NULL;
      safe_free(pp);
    }
    
    if ( atomMaskActive && residueMaskActive ) {
      for (i=0; i < residues; i++) {

	if ( residueMask[i] == 0 ) {
	  for (j = ipres[i]-1; 
	       j < ipres[i+1]-1;
	       j++) {
	    atomMask[j] = 0;
	  }
	}
      }
    } else if ( residueMaskActive ) {
      for (i=0; i < residues; i++) {

	for (j = ipres[i]-1;
	     j < ipres[i+1] - 1;
	     j++) {
	  if ( residueMask[i] )
	    atomMask[j] = 1;
	  else if (not)
	    atomMask[j] = 0;
	}
      }
    }
    
    atomMaskActive = 0;
    residueMaskActive = 0;
  }



 clean_return:

  actualAtoms = 0;
  for (i=0; i < atoms; i++)
    if ( atomMask[i] ) actualAtoms++;


  if ( tmp = strchr(maskString, '\n') )
    tmp[0] = (char) 0;

  if (actualAtoms > 0) {
    fprintf(stdout, "Mask %s%s] represents %i atoms\n", 
	    (not ? "[~" : "["), maskString, actualAtoms);
  } else  {
    fprintf(stdout, "Mask %s%s] represents %i atoms ",
            (not ? "[~" : "["), maskString, actualAtoms);
    fprintf(stdout, "!!!NO ATOMS DETECTED!!!\n");
    safe_free(atomMask);
    atomMask = NULL;
  }

  safe_free( residueMask );
  return ( atomMask );

}


   int *
processAtomMaskDetailedVH( char *maskString, int atoms, int residues,
			   Name *atomName, Name *residueName, int *ipres)
{
  char *charmask;
  int *intmask;
  int i, actualAtoms;

  intmask = (int *) safe_malloc( atoms * sizeof(int) );

  charmask = parseMaskString(maskString, atoms, residues, atomName, residueName, ipres);
  
  /*
   *  the new mask parsing routine returns a character array
   *  rather than integer; this makes more sense, however in the meantime
   *  we need to convert between the two.
   */

  actualAtoms = 0;
  for (i = 0; i < atoms; i++) {
    if (charmask[i] == 'T') {
      intmask[i] = 1;
      actualAtoms++;
    } else
      intmask[i] = 0;
  }

  if (actualAtoms > 0) {
    fprintf(stdout, "Mask [%s] represents %i atoms\n", 
	    maskString, actualAtoms);
  } else  {
    fprintf(stdout, "Mask [%s] represents %i atoms ",
            maskString, actualAtoms);
    fprintf(stdout, "!!!NO ATOMS DETECTED!!!\n");
    safe_free(intmask);
    intmask = NULL;
  }

  if (prnlev > 2) {
    fprintf(stdout, "Parsed mask string matches:\n");
    printAtomMaskDetailed(intmask, atoms, residues, atomName, residueName, ipres);
  }

  safe_free(charmask);
  return(intmask);

}


   int *
processAtomMask( char *maskString, ptrajState *state )
{

  if ( maskString && (strchr(maskString, '!') != NULL ||
		      strchr(maskString, '&') != NULL ||
		      strchr(maskString, '|') != NULL)) {

    return (processAtomMaskDetailedVH(maskString, state->atoms, state->residues, 
				      state->atomName, state->residueName, state->ipres));

  } else {

    return (processAtomMaskDetailed(maskString, state->atoms, state->residues,
				    state->atomName, state->residueName, state->ipres));
  }

}


   int *
processAtomMaskWrapper(char *buffer, ptrajState *state, int printit, int returnit)
{
  int *mask;


  if ( buffer && (strchr(buffer, '!') != NULL ||
		  strchr(buffer, '&') != NULL ||
		  strchr(buffer, '|') != NULL)) {

    mask = processAtomMaskDetailedVH(buffer, state->atoms, state->residues,
				     state->atomName, state->residueName, state->ipres);

  } else {

    mask = processAtomMaskDetailed(buffer, state->atoms, state->residues,
				   state->atomName, state->residueName, state->ipres);

  }

  if (printit) {
    printAtomMaskDetailed(mask,
			  state->atoms,
			  state->residues,
			  state->atomName,
			  state->residueName,
			  state->ipres);
  }

  if (returnit) {
    return(mask);
  } else {
    safe_free(mask);
    return(NULL);
  }

}



   int *
processAtomMaskWrapperOLD(char *buffer, int printit, int returnit)
{
  int *mask;
  int i;
  ptrajState *state;

  state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState(state);

  state->atoms = parm->NTOTAT;
  state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
  state->masses = (double *) safe_malloc(sizeof(double) * state->atoms);
  state->charges = NULL;
  for (i=0; i < state->atoms; i++) {
    strcpy(state->atomName[i], parm->atom[i].igraph);
    state->masses[i] = parm->atom[i].amass;
  }

  state->residues = parm->NTOTRS;
  state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
  for (i=0; i < state->residues; i++) {
    strcpy(state->residueName[i], parm->residue[i].labres);
  }
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  for (i=0; i <= state->residues; i++) {
    state->ipres[i] = parm->residue[i].ipres;
  }

  mask = processAtomMaskDetailed(buffer,
				 state->atoms,
				 state->residues,
				 state->atomName,
				 state->residueName,
				 state->ipres);

  if (printit) {
    printAtomMaskDetailed(mask,
			  state->atoms,
			  state->residues,
			  state->atomName,
			  state->residueName,
			  state->ipres);
  }

  if ( state != NULL ) {
    safe_free( state->atomName );
    safe_free( state->residueName );
    safe_free( state->ipres );
    safe_free( state->masses );
    safe_free( state );
  }

  if (returnit) {
    return(mask);
  } else {
    safe_free(mask);
    return(NULL);
  }

}




  scalarInfo *
scalarStackGetName(stackType **scalarStackP, char *name)
{
  stackType *s;
  scalarInfo *info, *match;

  match = NULL;
  for (s = *scalarStackP; s != NULL; s = s->next) {
    info = (scalarInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

  transformHBondInfo *
hbondInfoStackGetName(stackType **scalarStackP, char *name)
{
  stackType *s;
  transformHBondInfo *info, *match;

  match = NULL;
  for (s = *scalarStackP; s != NULL; s = s->next) {
    info = (transformHBondInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}


  transformMatrixInfo *
matrixInfoStackGetName(stackType **matrixStackP, char *name)
{
  stackType *s;
  transformMatrixInfo *info, *match;

  match = NULL;
  for (s = *matrixStackP; s != NULL; s = s->next) {
    info = (transformMatrixInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

  modesInfo *
modesInfoStackGetName(stackType **modesStackP, char *name)
{
  stackType *s;
  modesInfo *info, *match;

  match = NULL;
  for (s = *modesStackP; s != NULL; s = s->next) {
    info = (modesInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

   void
boxToRecip(double box[6], double ucell[9], double recip[9])
{
  double u12x,u12y,u12z;
  double u23x,u23y,u23z;
  double u31x,u31y,u31z;
  double volume;

  ucell[0] = box[0];
  ucell[1] = 0.0;
  ucell[2] = 0.0;
  ucell[3] = box[1]*cos(DEGRAD*box[5]);
  ucell[4] = box[1]*sin(DEGRAD*box[5]);
  ucell[5] = 0.0;
  ucell[6] = box[2]*cos(DEGRAD*box[4]);
  ucell[7] = (box[1]*box[2]*cos(DEGRAD*box[3]) - ucell[6]*ucell[3]) / ucell[4];
  ucell[8] = sqrt(box[2]*box[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);

  u23x = ucell[4]*ucell[8] - ucell[5]*ucell[7];
  u23y = ucell[5]*ucell[6] - ucell[3]*ucell[8];
  u23z = ucell[3]*ucell[7] - ucell[4]*ucell[6];
  u31x = ucell[7]*ucell[2] - ucell[8]*ucell[1];
  u31y = ucell[8]*ucell[0] - ucell[6]*ucell[2];
  u31z = ucell[6]*ucell[1] - ucell[7]*ucell[0];
  u12x = ucell[1]*ucell[5] - ucell[2]*ucell[4];
  u12y = ucell[2]*ucell[3] - ucell[0]*ucell[5];
  u12z = ucell[0]*ucell[4] - ucell[1]*ucell[3];
  volume=ucell[0]*u23x + ucell[1]*u23y + ucell[2]*u23z;

  recip[0] = u23x/volume;
  recip[1] = u23y/volume;
  recip[2] = u23z/volume;
  recip[3] = u31x/volume;
  recip[4] = u31y/volume;
  recip[5] = u31z/volume;
  recip[6] = u12x/volume;
  recip[7] = u12y/volume;
  recip[8] = u12z/volume;

}

/*
 *  Calculate the minimum possible distance between periodic images.
 *  This routine assumes that the ucell and recip information have
 *  previously been set via a call to boxToRecip
 */

   double
calculateMinImagedDistance2(double *box, double *ucell, double *recip,
			    double x1, double y1, double z1,
			    double x2, double y2, double z2,
			    int *rix, int *riy, int *riz)
{
  double X, Y, Z;
  double min, dist;
  double fx, fy, fz;
  int ix, iy, iz, i;

  if (box == NULL)
    return(0.0);

  min = 100.0 * (box[0]*box[0]+box[1]*box[1]+box[2]*box[2]);

  if (prnlev > 6) {
    fprintf(stdout, "ATOM      0  XXX A1      1     %7.3f %7.3f %7.3f\n",
	    x1, y1, z1);
    fprintf(stdout, "ATOM      1  XXX A2      1     %7.3f %7.3f %7.3f\n",
	    x2, y2, z2);

  }

  fx = x1*recip[0] + y1*recip[1] + z1*recip[2];
  fy = x1*recip[3] + y1*recip[4] + z1*recip[5];
  fz = x1*recip[6] + y1*recip[7] + z1*recip[8];

  if (prnlev > 6) {
    i = 2;
    for (ix = -1; ix <= 1; ix++) {
      for (iy = -1; iy <= 1; iy++) {
	for (iz = -1; iz <= 1; iz++) {

	  X = (fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6];
	  Y = (fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7];
	  Z = (fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8];

	  fprintf(stdout, "ATOM    %3i  XXX B%-2i     1     %7.3f %7.3f %7.3f\n",
		  i, i, X, Y, Z);
	  i++;
	}
      }
    }
  }

  i = 1;
  for (ix = -1; ix <= 1; ix++) {
    for (iy = -1; iy <= 1; iy++) {
      for (iz = -1; iz <= 1; iz++) {

	X = ((fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6]) - x2;
	Y = ((fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7]) - y2;
	Z = ((fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8]) - z2;
	dist = X*X + Y*Y + Z*Z;

	if (prnlev > 6) {
	  fprintf(stdout, "  IMAGE FAMILIAR  distance %3i: %6.3f  (%5i %5i %5i)\n",
		  i++, dist, ix, iy, iz);
	}

	if (dist < min) {
	  min = dist;
	  *rix = ix;
	  *riy = iy;
	  *riz = iz;
	}
      }
    }
  }
  if (prnlev > 4) {
    fprintf(stdout, "  IMAGE FAMILIAR, min distance is %6.3f (%5i %5i %5i)\n",
	    min, *rix, *riy, *riz);
  }

  return(min);
}



/*
 *  This routine will calculate a distance squared performing imaging,
 *  such that the minimum imaged distance (squared) is returned.  If
 *  the box is orthorhomic, simply subtract off multiple of the box lengthes.
 *  For non-orthorhomic, this is a little more tricky since this procedure
 *  is no longer applicable.  In this case, the ucell, recip and closest2
 *  information is used.  "closest2" represents a cutoff distance, such that if
 *  the calculate distance**2 is less than this return it (without calculating
 *  all possible image distances).  All possible images in each direction are 
 *  investigated currently.  This should be made "smarter" to only search in the
 *  appropriate octant based on the angle values.
 */
   double
calculateDistance2(int i, int j, double *x, double *y, double *z, 
		   double *box, double *ucell, double *recip, double closest2)
     /*
       i -- index for first atom
       j -- index for second atom
       x -- coordinate arrays
       y
       z
       box -- box coordinates
       ucell
       recip
       closest2 -- minimum distance between a periodic image 
     */
{
  double X, Y, Z, W;
  double fx, fy, fz, min, dist;
  int ix, iy, iz;

  X = x[i] - x[j];
  Y = y[i] - y[j];
  Z = z[i] - z[j];

  if (box == NULL || box[0] == 0.0)
    return(X*X + Y*Y + Z*Z);
    
  if (box[3] == 90.0 && box[4] == 90.0 && box[5] == 90.0) {
    /*
     *  DO ORTHORHOMIBIC IMAGING (this is faster!)
     */

       /*
        *  rid sign information
        */
    if ( X < 0 ) X = -X;
    if ( Y < 0 ) Y = -Y;
    if ( Z < 0 ) Z = -Z;

       /*
        *  rid multiples of the box length and image
        */
    while ( X > box[0] ) X = X - box[0];
    while ( Y > box[1] ) Y = Y - box[1];
    while ( Z > box[2] ) Z = Z - box[2];

       /*
        *  find shortest distance in the periodic reference
        */
    W = box[0] - X;
    if ( W < X ) X = W;

    W = box[1] - Y;
    if ( W < Y ) Y = W;

    W = box[2] - Z;
    if ( W < Z ) Z = W;

    return ( X*X + Y*Y + Z*Z );

  } else {

    /*
     *  NON-ORTHORHOMBIC CASE: find shortest distance in periodic reference
     *  This is a brute force check requiring up to 26 distance evaluations.
     *  It has been adapted to be smarter by returning the first distance that
     *  is shorter than the minimum possible distance between images.
     */

    fx = x[j]*recip[0] + y[j]*recip[1] + z[j]*recip[2];
    fy = x[j]*recip[3] + y[j]*recip[4] + z[j]*recip[5];
    fz = x[j]*recip[6] + y[j]*recip[7] + z[j]*recip[8];

    X = fx*ucell[0] + fy*ucell[3] + fz*ucell[6] - x[i];
    Y = fx*ucell[1] + fy*ucell[4] + fz*ucell[7] - y[i];
    Z = fx*ucell[2] + fy*ucell[5] + fz*ucell[8] - z[i];
    
    min = X*X + Y*Y + Z*Z;

    if (prnlev > 3) {
      if (min > 0.0)
	printf("DISTANCE: NON IMAGED        is %8.3f\n", sqrt(min));
      printf("DISTANCE BOX:               is %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
	     box[0],box[1],box[2],box[3],box[4],box[5]);
    }

    if (closest2 != 0.0 && min < closest2) return (min);

    for (ix = -1; ix <= 1; ix++) {
      for (iy = -1; iy <= 1; iy++) {
	for (iz = -1; iz <= 1; iz++) {

	  if (! (ix == 0 && iy == 0 && iz == 0) ) {
	    X = ((fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6]) - x[i];
	    Y = ((fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7]) - y[i];
	    Z = ((fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8]) - z[i];
	    dist = X*X + Y*Y + Z*Z;

	    if (prnlev > 3)
	      printf("DISTANCE: %3i %3i %3i  is %8.3f\n", ix, iy, iz, sqrt(dist));

	    if (dist < min) {
	      min = dist;
	      if (closest2 != 0.0 && min < closest2) 
		return(min);
	    }
	  }
	}
      }
    }
    return(min);
  }
}



   ptrajState **
ptrajCurrentState()
{
  return(&ptrajCurrentStatePointer);
}

   ptrajState *
ptrajCopyState(ptrajState **stateinp)
{
  ptrajState *state, *statein;
  int i;

  /*
   *  Make a copy of the current state
   */
  statein = *stateinp;

  state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState(state);

  state->atoms = statein->atoms;
  state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
  state->masses = (double *) safe_malloc(sizeof(double) * state->atoms);
  state->charges = (double *) safe_malloc(sizeof(double) * state->atoms);
  for (i=0; i < state->atoms; i++) {
    strcpy(state->atomName[i], statein->atomName[i]);
    state->masses[i] = statein->masses[i];
    state->charges[i] = statein->charges[i];
  }

  state->residues = statein->residues;
  state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
  for (i=0; i < state->residues; i++) {
    strcpy(state->residueName[i], statein->residueName[i]);
  }
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  for (i=0; i <= state->residues; i++) {
    state->ipres[i] = statein->ipres[i];
  }

  state->solventMolecules = statein->solventMolecules;
  state->solventAtoms     = statein->solventAtoms;

  if (statein->solventMolecules) {

    state->solventMask = (int *) safe_malloc(sizeof(int) * statein->atoms);
    for (i=0; i < state->atoms; i++)
      state->solventMask[i] = statein->solventMask[i];

    state->solventMoleculeStart = (int *) 
      safe_malloc(sizeof(int) * statein->solventMolecules);
    state->solventMoleculeStop  = (int *) 
      safe_malloc(sizeof(int) * statein->solventMolecules);
    for (i=0; i < state->solventMolecules; i++) {
      state->solventMoleculeStart[i] = statein->solventMoleculeStart[i];
      state->solventMoleculeStop[i]  = statein->solventMoleculeStop[i];
    }
  } else {
    state->solventMoleculeStart = NULL;
    state->solventMoleculeStop  = NULL;
    state->solventMask = NULL;
  }

  state->IFBOX = statein->IFBOX;
  for (i=0; i < 6; i++)
    state->box[i] = statein->box[i];

  if ( statein->molecules > 0 ) {

    state->molecules = statein->molecules;
    state->moleculeInfo = (int *)
      safe_malloc( sizeof(int) * state->molecules );
    for (i=0; i < state->molecules; i++) {
      state->moleculeInfo[i] = statein->moleculeInfo[i];
    }
  }

  return(state);
}



   void
checkAtomMask(char *buffer)
{
  ptrajState **statep;

  statep = ptrajCurrentState();
  processAtomMaskWrapper(buffer, *statep, 1, 0);
}



   int *
returnAtomMask(char *buffer)
{
  ptrajState **statep;

  statep = ptrajCurrentState();
  return( processAtomMaskWrapper(buffer, *statep, 0, 1) );
}


   void
atomMaskIsActive(int *mask, ptrajState *state, int *activep, int *firstp)
{
  int i, active, first;

  if (mask == NULL) {
    *activep = 0;
    return;
  }

  active = 0;
  for (i=0; i < state->atoms; i++) {
    if (mask[i] != 0) {
      active++;
      if (active == 1) first = i;
    }
  }
  *activep = active;
  *firstp = first;
}



#undef  ROUTINE
#define ROUTINE "ptrajPrintState()"

   void
ptrajPrintState(ptrajState *state)
{
  int i,j;
  int curres;

  printf("Dumping state information...\n");
  printf("  atoms:      %i\n", state->atoms);
  printf("  residues:   %i\n", state->residues);
  printf("  box length: %8.3f %8.3f %8.3f\n", state->box[0], state->box[1], state->box[2]);
  printf("  box angles: %8.3f %8.3f %8.3f\n", state->box[3], state->box[4], state->box[5]);
  printf("  molecules:  %i\n", state->molecules);
  printf("  max frames: %i\n", state->maxFrames);
  if (prnlev > 1) {
    printf("  ATOM   NAME  RESIDUE NAME  CHARGE    MASS\n");
    curres = 0;
    for (i=0; i < state->atoms; i++) {
    
      if (i == state->ipres[curres+1]-1) curres++;

      printf("  %8i %s %8i %s %8.3f %8.3f\n", 
	     i+1, state->atomName[i], 
	     curres+1, state->residueName[curres],
	     state->charges[i], state->masses[i]);
    }
    printf("  Molecule information (atoms in each molecule):\n");
    for (i=0; i < state->molecules; i++) {
      printf("  %5i", state->moleculeInfo[i]);
      if (i!=0 && i%10==0) printf("\n");
    }
    printf("\n");
  }
  if (state->solventMolecules > 0) {
    printf("  solvent molecules: %i (%i atoms)\n", 
	   state->solventMolecules, state->solventAtoms);
    printf("  solvent mask is: ");
    printAtomMask(state->solventMask, state);
    if (prnlev > 1) {
      for (i=0; i < state->solventMolecules; i++) {
	printf("  SOLVENT %4i: ", i+1);
	for (j=state->solventMoleculeStart[i]; j < state->solventMoleculeStop[i]; j++) {
	  printf(" (%s %5i)", state->atomName[j], j+1);
	}
	printf("\n");
      }
    }
  }
}


   void
ptrajClearState(ptrajState **stateinp)
{
  ptrajState *state;

  state = *stateinp;

  if (state != NULL) {
    safe_free( state->atomName );
    safe_free( state->residueName );
    safe_free( state->ipres );
    safe_free( state->masses );
    safe_free( state->charges );
    if (state->solventMolecules) {
      safe_free( state->solventMoleculeStart );
      safe_free( state->solventMoleculeStop );
      safe_free( state->solventMask );
    }
    safe_free( state->moleculeInfo );
    INITIALIZE_ptrajState(state);
    safe_free( state );
  }
}




/*
 *  checkCoordinates(): a routine to determine what type of coordinate
 *  file specified by "filename" is and how many coordinate frames the 
 *  file represents.  The file is opened up, checked, then the file
 *  is closed.  It is the responsibility of later routines to reopen
 *  the file.  [This requirement is so the file limit will not be blown
 *  in the case of a user processing a boatload of files; in this way,
 *  sequential file access is maintained.]
 *  
 *  On success, a "coordinateInfo *" will be returned with the filename,
 *  format and start/stop represented.  
 *
 *  On failure, NULL is returned.
 */

#undef  ROUTINE
#define ROUTINE "checkCoordinates()"

   coordinateInfo *
checkCoordinates(char *filename, int totalAtoms)
{
  double *x = NULL;
  float junk1, junk2, junk3, junk4, junk5, junk6, junk7, junk8, junk9;
  int i, actualAtoms, binposeof, ncid, err;
  int spatial = 0;
  int isBox = 0;
  int isVelocity = 0;
  int lines_per_set;
  int start = 1;
  int stop = 1;
  size_t ist;
  float *binposScratch;
  FILE *fp;
  pdb_record r;
  coordType type = COORD_UNKNOWN;
  char buffer1[BUFFER_SIZE], buffer2[BUFFER_SIZE];
  Restart *restrt;
  coordinateInfo *trajInfo;
  charmmTrajectoryInfo **charmmTrajectoryp = NULL;
  charmmTrajectoryInfo  *charmmTrajectory = NULL;
  netcdfTrajectoryInfo *NCInfo = NULL;

  trajInfo = (coordinateInfo *) safe_malloc(sizeof(coordinateInfo));
  INITIALIZE_coordinateInfo(trajInfo);
  trajInfo->filename = copyString(filename);

  /*
   *  Check to see if the file specified by "filename" is a NetCDF file;
   *  if it is NOT open the file in the regular manner...
   */

#ifdef BINTRAJ
  err = nc_open(filename, NC_NOWRITE, &ncid);
  if (err == NC_NOERR) {
    /*
     *  this appears to be a NetCDF file so initialize necessary data
     */
    type = COORD_AMBER_NETCDF;
    NCInfo = (netcdfTrajectoryInfo *) safe_malloc(sizeof(netcdfTrajectoryInfo));
    INITIALIZE_netcdfTrajectoryInfo( NCInfo );
    NCInfo->ncid = ncid;
  } else
#endif
    if ( openFile(&fp, filename, "r") == 0 ) {

      fprintf(stdout, "trajin %s ignored; could not open file (%s)\n", filename, filename);
      safe_free(trajInfo->filename);
      safe_free(trajInfo);
      return NULL;

  }

  /*
   *  Try to determine what "type" of file the input file is if it is
   *  -not- a NetCDF file.
   *
   *  NOTE: if new coordType's are implemented, this section may
   *  need to be modified.
   *
   *    (a) read two lines of text
   *    (b) check for pdb  [this requires reading both lines of
   *        text since the comment from an AMBER restrt or 
   *        trajectory may have text which "fools" the pdb routines...]
   *    (c) look for the BINPOS magic string
   *    (d) look for CHARMM trajectory magic string
   *    (e) if not a pdb or BINPOS, check the second line to see if it 
   *        has three numbers which implies it is an AMBER trajectory file.
   */

  if (type == COORD_UNKNOWN) {

    if (fgets(buffer1, BUFFER_SIZE, fp) == NULL ||
	fgets(buffer2, BUFFER_SIZE, fp) == NULL ) {
      fprintf(stdout, "%s: EOF encountered prematurely (%s)\n", ROUTINE, filename);
      safe_free(trajInfo->filename);
      safe_free(trajInfo);
      return NULL;
    }

      /*
       *  Is this a NetCDF file (and NetCDF support is not compiled in)
       */
#ifndef BINTRAJ
    if (buffer1[0] == 'C' && 
	buffer1[1] == 'D' && 
	buffer1[2] == 'F') {
      fprintf(stdout, "    *** NetCDF trajectory file found but support for this type of file\n");
      fprintf(stdout, "    is not compiled into this version of ptraj!  Recompile with NetCDF\n");
      fprintf(stdout, "    binary support (see configure and make sure -DBINTRAJ is defined)\n");
      fprintf(stdout, "trajin %s ignored; support for NetCDF file not enabled, file (%s)\n", filename, filename);
      safe_free(trajInfo->filename);
      safe_free(trajInfo);
      return NULL;
    } 
#endif
      /*
       *  Is this file a PDB file?
       */
    r = pdb_read_string(buffer1);


    if (r.record_type != PDB_UNKNOWN) {
      r = pdb_read_string(buffer2);
      printf("Reading string (%s) R = %i\n", buffer2, r);
      if (r.record_type != PDB_UNKNOWN)
	type = COORD_PDB;
    } else if (strncmp(buffer1, "REMARK", 6) == 0)
      /*
       *  this is a hack to allow CHARMM pdb files to be read or other PDB's 
       *  which have malformed remarks
       */
      type = COORD_PDB;


      /*
       *  Does this file have the magic header for the Scripps binary format?
       */
    if (buffer1[0] == 'f' && 
	buffer1[1] == 'x' && 
	buffer1[2] == 'y' &&
	buffer1[3] == 'z') {
      type = COORD_BINPOS;
    }

      /*
       *  Does file file have the magic header for a CHARMM trajectory???
       */
    if (buffer1[4] == 'C' &&
	buffer1[5] == 'O' &&
	buffer1[6] == 'R' &&
	buffer1[7] == 'D') {
      type = COORD_CHARMM_TRAJECTORY;
    }

      /*
       *  Can we scan in three numbers from the first line?  If so this is 
       *  an AMBER trajectory otherwise assume it is an AMBER restrt file
       */
    if ( type == COORD_UNKNOWN ) {
      if ( (i = sscanf(buffer2, "%f%f%f", &junk1, &junk2, &junk3)) == 3 )
	type = COORD_AMBER_TRAJECTORY;
      else
	type = COORD_AMBER_RESTART;
    }
  }

  /*
   *  EXIT if we have not figured out what kind of file this is
   */

  if ( type == COORD_UNKNOWN ) {
    fprintf(stdout, "  WARNING: trajin %s, unsupported file format...\n",
	    filename);
    safe_free(trajInfo->filename);
    safe_free(trajInfo);
    return NULL;
  }

  /*
   *  REOPEN THE FILE (if the file is not a NetCDF file)
   */

  if ( type != COORD_AMBER_NETCDF && (fp = safe_freopen(fp)) == NULL ) {
    safe_free(trajInfo->filename);
    safe_free(trajInfo);
    return NULL;
  }

     /*
      *  Read through the file, depending on type, and make sure
      *  that they have the same number of atoms and casually check
      *  to make sure the file isn't corrupted.  In the case of a 
      *  trajectory, check to see how many frames are present in the
      *  file...
      */

  actualAtoms = 0;

  switch ( type ) {

  case COORD_AMBER_NETCDF:

#ifdef BINTRAJ
    /*
     *  Get global attributes (after initializing a structure to hold the data), 
     */
    
    trajInfo->title =           netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "title");
    trajInfo->application =     netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "application");
    trajInfo->program =         netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "program");
    trajInfo->version =         netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "programVersion");
    NCInfo->Conventions =       netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "Conventions");
    NCInfo->ConventionVersion = netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "ConventionVersion");
    if (strstr(NCInfo->Conventions, "AMBER") == NULL) {
      fprintf(stdout, "WARNING: NetCDF file has Conventions that do not include the string \"AMBER\"\n");
    }
    if (strcmp(NCInfo->ConventionVersion, "1.0") != 0) {
      fprintf(stdout, "WARNING: NetCDF file has ConventionVersion differing from \"1.0\"\n");
    }

    /*
     *  get the NetCDF dimension ID's and sizes for the frames, spatial and atoms
     */
    
    NCInfo->frameDID   =  netcdfGetDimensionInfo(NCInfo->ncid, AMBER_NETCDF_FRAME, &stop);
    NCInfo->spatialDID =  netcdfGetDimensionInfo(NCInfo->ncid, AMBER_NETCDF_SPATIAL, &spatial);
    NCInfo->atomDID    =  netcdfGetDimensionInfo(NCInfo->ncid, AMBER_NETCDF_ATOM, &actualAtoms);

    if (spatial != 3) {
      type = COORD_UNKNOWN;
      fprintf(stdout, "ptraj cannot handle NetCDF files with other than 3 dims\n");
    }

    /*
     *  perform a sanity check on time variable and units
     */
    err = nc_inq_varid(NCInfo->ncid, AMBER_NETCDF_TIME, &NCInfo->timeVID);
    if (err != NC_NOERR) error(ROUTINE, "NetCDF time variable, error: %s", nc_strerror(err));
    err = nc_inq_varnatts(NCInfo->ncid, NCInfo->timeVID, &i);
    if ( err ) {
      error(ROUTINE, "Getting number of time attributes in NetCDF file %s\n", filename);
    } 

    if ( i != 1 ) error(ROUTINE, "Only one time attribute is expected in NetCDF file %s\n", filename);

    err = nc_inq_attlen(NCInfo->ncid, NCInfo->timeVID, "units", &ist);
    if (err == NC_NOERR)
      NCInfo->timeUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
    else
      fprintf(stdout, "NetCDF grabbing time units attribute length error: %s\n", nc_strerror(err));

    err = nc_get_att_text(NCInfo->ncid, NCInfo->timeVID, "units", NCInfo->timeUnits);
    if (err != NC_NOERR) {
      fprintf(stdout, "%s: Could not get time units from NetCDF file %s: %s", ROUTINE,
	      filename, nc_strerror(err));
      safe_free(NCInfo->timeUnits);
    }

    if (strcmp("picosecond",NCInfo->timeUnits) != 0)
      fprintf(stdout, "WARNING: Expecting time units in picoseconds, got -%s-\n", NCInfo->timeUnits);
    


    err = nc_inq_varid(NCInfo->ncid, AMBER_NETCDF_SPATIAL, &NCInfo->spatialVID);
    if (err != NC_NOERR)
      fprintf(stdout, "NetCDF error getting spatialVID in the NetCDF file %s: %s\n", 
	      filename, nc_strerror(err));


    /*
     *  NETCDF: check to see what trajectory data is within the NetCDF file
     *  and perform sanity checks...
     */

       /*
        *  ARE COORDINATES PRESENT?
        */ 
    err = nc_inq_varid(NCInfo->ncid, AMBER_NETCDF_COORDS, &NCInfo->coordinateVID);
    if (err != NC_NOERR) {
      fprintf(stdout, "No coordinates are present in the NetCDF file %s\n", filename);
    } else {

      err = nc_inq_varnatts(NCInfo->ncid, NCInfo->coordinateVID, &i);
      if ( err ) {
	error(ROUTINE, "Getting number of coordinate attributes in NetCDF file %s\n", filename);
      } else {
	if ( i != 1 ) {
	  error(ROUTINE, "Only a single coordinate attribute is expected in NetCDF file %s\n", filename);
	} else {
	  err = nc_inq_attlen(NCInfo->ncid, NCInfo->coordinateVID, "units", &ist);
	  if (err == NC_NOERR) {
	    NCInfo->coordinateUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
	    err = nc_get_att_text(NCInfo->ncid, NCInfo->coordinateVID, "units", NCInfo->coordinateUnits);
	    if (err != NC_NOERR) {
	      fprintf(stdout, "%s: Could not get coordinate units from NetCDF file %s", ROUTINE, filename);
	      safe_free(NCInfo->coordinateUnits);
	    }
	    if (strcmp("angstrom",NCInfo->coordinateUnits) != 0)
	    fprintf(stdout, "WARNING: Expecting coordinate units in angstroms, got %s\n", NCInfo->coordinateUnits);
	  } else
	    fprintf(stdout, "NetCDF grabbing coordinate units attribute error: %s\n", nc_strerror(err));
	}
      }
    }

       /*
        *  ARE CELL_LENGTHS PRESENT?
        */ 
    err = nc_inq_varid(NCInfo->ncid, "cell_angles", &NCInfo->cellAngleVID);
    if (err != NC_NOERR) fprintf(stdout, "NetCDF inquiring about cell angle variable ID: %s", nc_strerror(err));
    err = nc_inq_varid(NCInfo->ncid, "cell_lengths", &NCInfo->cellLengthVID);
    if (err != NC_NOERR) fprintf(stdout, "NetCDF inquiring about cell length variable ID: %s", nc_strerror(err));

    if (err == NC_NOERR) {

      isBox = 1;

      err = nc_inq_varnatts(NCInfo->ncid, NCInfo->cellAngleVID, &i);
      if ( i > 0 ) {
	err = nc_inq_attlen(NCInfo->ncid, NCInfo->cellAngleVID, "units", &ist);
	if (err == NC_NOERR) {
	  NCInfo->cellAngleUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
	  err = nc_get_att_text(NCInfo->ncid, NCInfo->cellAngleVID, "units", NCInfo->cellAngleUnits);
	}
      }
      err = nc_inq_varnatts(NCInfo->ncid, NCInfo->cellLengthVID, &i);
      if ( i > 0 ) {
	err = nc_inq_attlen(NCInfo->ncid, NCInfo->cellLengthVID, "units", &ist);
	if (err == NC_NOERR) {
	  NCInfo->cellLengthUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
	  err = nc_get_att_text(NCInfo->ncid, NCInfo->cellLengthVID, "units", NCInfo->cellLengthUnits);
	}
      }
    }

       /*
        *  ARE VELOCITIES PRESENT?
        */ 
    err = nc_inq_varid(NCInfo->ncid, "velocities", &NCInfo->velocityVID);
    if (err == NC_NOERR) {

      isVelocity = 1;

      err = nc_inq_varnatts(NCInfo->ncid, NCInfo->velocityVID, &i);
      if ( i > 1 ) {

	err = nc_inq_attlen(NCInfo->ncid, NCInfo->velocityVID, "units", &ist);
	if (err == NC_NOERR) {
	  NCInfo->velocityUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
	  err = nc_get_att_text(NCInfo->ncid, NCInfo->velocityVID, "units", NCInfo->velocityUnits);
	}
	err = nc_get_att_double(NCInfo->ncid, NCInfo->velocityVID, "scale_factor", &NCInfo->velocityScale);
      }
    }
#endif
    break;

  case COORD_PDB:

       /*
        *  read PDB records until the EOF is encountered counting
        *  the number of atoms encountered along the way...
        *  NOTE: this code will NOT properly handle pdbs that are
        *  strung together into a single file.
        */

    for (;;) {
      r = pdb_read_record(fp);

      if ( r.record_type == PDB_ATOM || r.record_type == PDB_HETATM )
	actualAtoms++;

      if ( r.record_type == PDB_END ) break;
    }

    break;

  case COORD_AMBER_RESTART:


    restrt = readAmberRestart(totalAtoms, filename, fp);
    if (restrt == NULL) return NULL;
    actualAtoms = restrt->natoms;
    isBox = restrt->isbox;
    isVelocity = restrt->restart;

    break;

  case COORD_AMBER_TRAJECTORY:

       /*
        *  there is no perfect, easy, way to make sure that the 
        *  trajectory file has the "correct" number of atoms 
        *  in a general way.  If the number of atoms * 3 is not
        *  a multiple of 10, it is possible; however we do not check
        *  this currently.
        */

    actualAtoms = totalAtoms;

    if (fgets(buffer1, BUFFER_SIZE, fp) == NULL) {
      fprintf(stdout, 
	      "WARNING in %s: EOF encountered during reading of\n", ROUTINE);
      fprintf(stdout, 
	      "   title from (%s)\n", filename);
      return NULL;
    }


    lines_per_set = (int) (actualAtoms * 3) / 10;
    if ( (actualAtoms * 3) % 10 ) lines_per_set += 1;

    for (i=0; fgets(buffer1, BUFFER_SIZE, fp) != NULL; i++) {
      if ( strchr(buffer1, (int) '*') != NULL ) {
	fprintf(stdout, 
		"WARNING in %s, AMBER trajectory file is corrupted: '*' detected (%s)\n", 
		ROUTINE, filename);
	break;
      }

      if (i == lines_per_set) {
	if (actualAtoms > 2 && 
	    sscanf(buffer1, "%f%f%f%f%f%f%f%f%f",
		   &junk1, &junk2, &junk3, &junk4, &junk5, &junk6,
		   &junk7, &junk8, &junk9) < 9)
	  isBox = 1;
      }
    }

    stop = i / (lines_per_set + isBox);
    break;


  case COORD_BINPOS:

    if (openbinpos(fp) < 0) break;

    binposScratch = safe_malloc(sizeof(float) * totalAtoms * 3);
    binposeof = 0;
    for (i=0; binposeof == 0; i++) {

      if (readbinpos(fp, &actualAtoms, binposScratch, &binposeof) < 0)
        binposeof = 1;
    }
    stop = i-1;
    safe_free(binposScratch);
	break; 

  case COORD_CHARMM_TRAJECTORY:

    /*
     *  pre-process
     */
    stop = -1;
    x = NULL;
    charmmTrajectoryp = (charmmTrajectoryInfo **) safe_malloc(sizeof(charmmTrajectoryInfo *));
    *charmmTrajectoryp = NULL;
    readCharmmTrajectory(fp, charmmTrajectoryp, x, x, x, x, stop);
    charmmTrajectory = *charmmTrajectoryp;

    x = (double *) safe_malloc(sizeof(double) * charmmTrajectory->natrec);

    /*
     *  read in sets until there are no more
     */
    stop = 0;
    while ( readCharmmTrajectory(fp, charmmTrajectoryp, x, x, x, x, stop+1) )
      stop++;
    actualAtoms = charmmTrajectory->natrec;
    isVelocity = 0;
    isBox = charmmTrajectory->icntrl[10];

    if (charmmTrajectory->icntrl[0] != stop) {
      fprintf(stdout, "NOTE: this charmm trajectory contains %i sets, expecting %i\n",
	      stop, charmmTrajectory->icntrl[0]);
    }

    safe_free(x);
    x = NULL;
  }


  if (actualAtoms != totalAtoms) {
    fprintf(stdout, 
	    "  WARNING in %s: The actual number of atoms (%d)\n", ROUTINE, actualAtoms );
    fprintf(stdout, "  does not match the expected number of atoms (%d) in (%s)\n", 
	    totalAtoms, filename);
    fprintf(stdout, "  With this version of the code, this will likely lead to program failure!!!\n");
  }


  trajInfo->start = start;
  trajInfo->stop = stop;
  trajInfo->offset = 1;
  trajInfo->isBox = isBox;
  trajInfo->isVelocity = isVelocity;
  trajInfo->type = type;

  switch (type) {
  case COORD_CHARMM_TRAJECTORY:
    trajInfo->info = (void *) charmmTrajectory;
    break;
#ifdef BINTRAJ
  case COORD_AMBER_NETCDF:
    trajInfo->info = (void *) NCInfo;
    err = nc_close(NCInfo->ncid);
    if (err)
      fprintf(stdout, "Error closing NetCDF file %s (nc_close), error: %s\n",
	      filename, nc_strerror(err));
#endif
    break;
  }

  if (type != COORD_AMBER_NETCDF)
    safe_fclose(fp);
  
  return ( trajInfo );
}

#undef  ROUTINE
#define ROUTINE "checkCoordinatesWrap()"

   void
checkCoordinatesWrap(char *filename)
{
  coordinateInfo *info;
  charmmTrajectoryInfo *charmmTraj;
  byte u;

  info = checkCoordinates(filename, parm->NTOTAT);

  if (info == NULL) {
    fprintf(stdout, "WARNING in %s.  Encountered a problem analyzing\n", ROUTINE);
    fprintf(stdout, "coordinates from file %s\n", filename);
    return;
  }

  switch ( info->type ) {

  case COORD_AMBER_NETCDF:
    fprintf(stdout, "File (%s) is a NetCDF AMBER trajectory%s",
	    info->filename, (info->isBox ? " with box coordinates" : ""));
    if (info->isVelocity > 0)
      fprintf(stdout, " with velocities");
    if (info->stop > 0)
      fprintf(stdout, " representing %i sets\n", info->stop);
    else
      fprintf(stdout, "\n");
    break;

  case COORD_PDB:
    fprintf(stdout, "File (%s) is a PDB file\n", info->filename);
    break;

  case COORD_AMBER_TRAJECTORY:
    fprintf(stdout, "File (%s) is an AMBER trajectory%s",
	    info->filename, (info->isBox ? " with box coordinates" : ""));
    if (info->stop > 0)
      fprintf(stdout, " representing %i sets\n", info->stop);
    else
      fprintf(stdout, "\n");
    break;

  case COORD_CHARMM_TRAJECTORY:

    charmmTraj = (charmmTrajectoryInfo *) info->info;
    fprintf(stdout,
	    "File (%s) is a CHARMM trajectory in %s endian binary format %s",
	    info->filename, (charmmTraj->byteorder ? "little" : "big"),
	    (info->isBox ? "with box coordinates" : ""));
    if (info->stop > 0)
      fprintf(stdout, "representing %i sets\n", info->stop);
    else
      fprintf(stdout, "\n");
    if (prnlev > 2) {
      printf("  NFILE = %i\n", charmmTraj->icntrl[0]); /* number of coordinate sets in file   */
      printf("  ISTEP = %i\n", charmmTraj->icntrl[1]); /* number of previous dynamics steps   */
      printf("  NINTV = %i\n", charmmTraj->icntrl[2]); /* frequency for saving coordinates    */
      printf("  NSAVC = %i\n", charmmTraj->icntrl[3]); /* number of steps for creation runs   */
      printf("  NSAVV = %i\n", charmmTraj->icntrl[4]);
      printf("  NDEGF = %i\n", charmmTraj->icntrl[7]);
      printf("  NFREA = %i\n", charmmTraj->icntrl[8]);
      u.i = charmmTraj->icntrl[9];
      printf("  DELTA = %f\n", u.f);
      printf("  QCRYS = %i\n", charmmTraj->icntrl[10]);
      printf("  QDIM4 = %i\n", charmmTraj->icntrl[11]);
      printf("  VERNU = %i\n", charmmTraj->icntrl[19]);
      printf("  Dumping the title:\n");
      printStack(&charmmTraj->titleStack, printString, NULL);
    }
    break;

  case COORD_AMBER_RESTART:
    fprintf(stdout, "File (%s) is an AMBER restart file ", info->filename);
    if ( info->isBox && info->isVelocity )
      fprintf(stdout, "with box and velocity information\n");
    else if (info->isBox)
      fprintf(stdout, "with box information\n");
    else if (info->isVelocity)
      fprintf(stdout, "with velocity information\n");
    else
      fprintf(stdout, "\n");
    break;
  }
}



#undef  ROUTINE
#define ROUTINE "loadCharmmPSF()"

   ptrajState *
loadCharmmPSF(FILE *fd, int skipheader)
{
  char *bufferp, *buffer, *pbuffer;
  int natom;

  double *charges, *masses;
  float charge, mass;
  int  *tmp_ipres, *tmp_molinfo;
  Name *atomName;
  Name residueNumber;
  Name *tmp_resName;
  Name segid, cursegid, vdwt;
  int atom;
  int residue;
  int fixed;
  int vdwp;
  int i,j,k,len, scannedValues;
  int curres, curresidx, curmol, numcurmol;
  int solventMolecules, *solventMoleculeStart, *solventMoleculeStop;
  int *solventMask, solventAtoms;
  int *tmp_solventMoleculeStart, *tmp_solventMoleculeStop;
  int ntitle;
  ptrajState *state;

  buffer = (char *) safe_malloc(sizeof(char) * BUFFER_SIZE);
  pbuffer = (char *) safe_malloc(sizeof(char) * BUFFER_SIZE);

     /*
      *  Open up the PSF file
      */
  if (fd == NULL) {
    error(ROUTINE, "Couldn't open file psf\n");
  }

     /*
      *  grab header and make sure it matches "PSF"
      */
  if (skipheader != 1) {
    bufferp = fgets(buffer, BUFFER_SIZE, fd);
    if (bufferp == NULL || strncmp(buffer, "PSF", 3) != 0) {
      error(ROUTINE, "File is not a CHARMM PSF file!\n");
    }
  }
  fprintf(stdout, "Reading in CHARMM PSF file\n");

     /*
      *  skip the next line, read the title and the blank line that follows
      */
  bufferp = fgets(buffer, BUFFER_SIZE, fd);
  bufferp = fgets(buffer, BUFFER_SIZE, fd);
  if (strstr(bufferp, "!NTITLE")) {
    if ( sscanf(bufferp, "%i", &ntitle) != 1 ) 
      ntitle = 0;
  }


  fprintf(stdout, "Reading in the title...\n\n");
  if (ntitle) {
    for (i=0; i<=ntitle; i++) {
      bufferp = fgets(buffer, BUFFER_SIZE, fd);
      fprintf(stdout, "%s", buffer);
    }

  } else {

    bufferp = fgets(buffer, BUFFER_SIZE, fd);
    while (buffer[0] == '*')
      bufferp = fgets(buffer, BUFFER_SIZE, fd);
      fprintf(stdout, "%s", buffer);
  }

     /*
      *  get the total number of atoms
      */
  bufferp = fgets(buffer, BUFFER_SIZE, fd);

  if (sscanf(buffer, "%i", &natom) != 1) error(ROUTINE, "Scanning natoms\n");
  fprintf(stdout, "Total number is atoms is %d\n", natom);


  /*
   *  read in the atom segid, residue, name and charge/mass info
   */
  charges     = (double *) safe_malloc(sizeof(double) * natom);
  masses      = (double *) safe_malloc(sizeof(double) * natom);
  tmp_ipres   = (int *)    safe_malloc(sizeof(double) * natom);
  tmp_molinfo = (int *)    safe_malloc(sizeof(double) * natom);

  atomName    = (Name *)   safe_malloc(sizeof(Name) * natom);
  tmp_resName = (Name *)   safe_malloc(sizeof(Name) * natom);

  tmp_solventMoleculeStart = (int *) safe_malloc(sizeof(int) * natom);
  tmp_solventMoleculeStop  = (int *) safe_malloc(sizeof(int) * natom);
  solventMask = (int *) safe_malloc(sizeof(int) * natom);
  memset(solventMask, 0, sizeof(int) * natom);
  memset(tmp_solventMoleculeStart, 0, sizeof(int) * natom);
  memset(tmp_solventMoleculeStop,  0, sizeof(int) * natom);

  curres = 1;
  curresidx = 1;
  tmp_ipres[0] = 1;
  solventMolecules = 0;
  solventAtoms = 0;

  fprintf(stdout, "Reading in the atom information...\n");

  for (i=0; i < natom; i++) {

    bufferp = fgets(buffer, BUFFER_SIZE, fd);

       /*
        * atom  seg res# res atom  vdwp  charge    mass    fixed?
        *  1    WAT  1   WAT  OH    2   -0.820000 15.9994    0
        */
    scannedValues = sscanf(buffer, "%8i %4s %4i %4s %4s %4i %f %f%8i",
			   &atom, (char *) &segid, &residue, (char *) &tmp_resName[i], 
			   (char *) &atomName[i], &vdwp, &charge, &mass, &fixed);

    if (scannedValues < 9) {

      /*
       *  assume this is an XPLOR psf file w/ text for the vdw type...
       *  NOTE: there was a funny problem with reading residues as numbers
       *  since when they get over 9999 they revert to 0001 and the sscanf choked.
       *  Replace the number with a string and then convert and it seems to work.
       */

      scannedValues = sscanf(buffer, "%8i %4s %4s %4s %4s %4s %f %f %i",
			     &atom, (char *) &segid, &residueNumber, (char *) &tmp_resName[i], 
			     (char *) &atomName[i], (char *) &vdwt, &charge, &mass, &fixed);
      if (scannedValues == 9) {	
	residue = atoi(residueNumber);
      } else {

	/*
	 *  assume this is an XPLOR psf file w/o SEGID information...
	 */
	if (sscanf(buffer, "%8i %8i %4s %4s %4s %f %f%8i",
		   &atom, &residue, 
		   (char *) &tmp_resName[i], (char *) &atomName[i], (char *) &vdwt,
		   &charge, &mass, &fixed) < 8) {
	  
	  error(ROUTINE, "Scanning atom information\n");
	}
      }
    }

    atomName[i][NAME_SIZE-1] = (char) 0;
    tmp_resName[i][NAME_SIZE-1] = (char) 0;
    for (j=NAME_SIZE-2; atomName[i][j] == (char) 0; j--) {
      atomName[i][j] = ' ';
    }
    for (j=NAME_SIZE-2; tmp_resName[i][j] == (char) 0; j--) {
      tmp_resName[i][j] = ' ';
    }

    if (i==0) {
      strcpy(cursegid, segid);
      curmol = 0;
      numcurmol = 1;
    } else {

      if (strcmp(cursegid, segid) == 0) {
	numcurmol++;
      } else {
	tmp_molinfo[curmol] = numcurmol;
	curmol++;
	numcurmol = 1;
	strcpy(cursegid, segid);
      }
    }


    charges[i] = charge;
    masses[i] = mass;

    if (prnlev > 1) 
      printf("Atom %4i (%s) seg (%4s) res %4i (%s) charge: %7.3f mass: %7.3f fixed? %s %i\n",
	     atom, atomName[i], segid, residue, tmp_resName[i], charge, mass,
	     (fixed ? "yes" : "no"), vdwp);


    if (residue != curres) {
      tmp_ipres[curresidx] = i+1;
      curresidx++;
      curres = residue;
    }
  }

  tmp_ipres[curresidx] = i+1;
  tmp_molinfo[curmol] = numcurmol;

  for (i=0; i < curresidx; i++) {
    j = tmp_ipres[i]-1;
    if ( strcmp(tmp_resName[j], "WAT ") == 0 || 
	 strcmp(tmp_resName[j], "IP3 ") == 0 ||
	 strcmp(tmp_resName[j], "TIP3") == 0) {
      tmp_solventMoleculeStart[solventMolecules] = j;
      tmp_solventMoleculeStop[solventMolecules] = j+3;
      solventMolecules++;
    }
  }

  solventMoleculeStart = (int *) safe_malloc(sizeof(int) * (solventMolecules+1));
  solventMoleculeStop  = (int *) safe_malloc(sizeof(int) * (solventMolecules+1));
  for (i = 0; i < solventMolecules; i++) {

    solventMoleculeStart[i] = tmp_solventMoleculeStart[i];
    solventMoleculeStop[i]  = tmp_solventMoleculeStop[i];

    for (j=solventMoleculeStart[i]; j < solventMoleculeStop[i]; j++) {
      solventMask[j] = 1;
      solventAtoms++;
    }
  }

  if (prnlev > 1) {
    printf("The total number of molecules is %i\n", curmol+1);
    printf("The total number of residues is %i\n", curresidx);
  }

  fprintf(stdout, "Dumping out residue names:\n");
  j = 0;
  k = 0;
  len = 0;
  pbuffer[0] = (char) 0;
  for (i=0; i < curresidx; i++) {

    k = tmp_ipres[i]-1;
    sprintf(buffer+j, "%5s", tmp_resName[k]);
    j += 5;
    if ( i == curresidx - 1 || ( i != 0 && (i+1) % 10 == 0 ) ) {
      buffer[j] = (char) 0;
      if ( strcmp(pbuffer, buffer) == 0 ) {
	if ( ! len ) {
	  fprintf(stdout, " ...\n");
	  len = 1;
	}
      } else {
	fprintf(stdout, "%s\n", buffer);
	len = 0;
	strcpy(pbuffer, buffer);
      }
      j = 0;
    }
  }

  state = (ptrajState *) safe_malloc(sizeof(ptrajState));
  INITIALIZE_ptrajState(state);

  state->box[0] = 0.0;
  state->box[1] = 0.0;
  state->box[2] = 0.0;
  state->box[3] = 90.0;
  state->box[4] = 90.0;
  state->box[5] = 90.0;
  state->atoms = natom;
  state->residues = curresidx;
  state->charges = charges;
  state->masses = masses;
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  state->residueName = (Name *) safe_malloc(sizeof(Name) * state->residues);
  for (i=0; i <= state->residues; i++)
    state->ipres[i] = tmp_ipres[i];
  for (i=0; i < state->atoms; i++) {
    curresidx = atomToResidue(i+1, state->residues, state->ipres)-1;
    strcpy(state->residueName[curresidx], tmp_resName[i]);
    i = state->ipres[curresidx+1]-2;
  }

  state->IFBOX = 0;
  state->molecules = curmol+1;
  state->moleculeInfo = (int *) safe_malloc(sizeof(int) * state->molecules);
  for (i=0; i < state->molecules; i++) {
    state->moleculeInfo[i] = tmp_molinfo[i];
  }

  state->solventMask = solventMask;
  state->solventAtoms = solventAtoms;
  state->solventMolecules = solventMolecules;
  state->solventMoleculeStart = solventMoleculeStart;
  state->solventMoleculeStop = solventMoleculeStop;

  state->atomName = atomName;

  safe_free(buffer);
  safe_free(tmp_ipres);
  safe_free(tmp_molinfo);
  safe_free(tmp_resName);
  safe_free(tmp_solventMoleculeStart);
  safe_free(tmp_solventMoleculeStop);

  if (prnlev >= 0) 
    ptrajPrintState(state);
  return state;
  
}



#undef  ROUTINE
#define ROUTINE "ptrajInitializeState()"

   void
ptrajInitializeState(ptrajState **statep, char *filename)
{
  FILE *fp;
  char *buffer;
  ptrajState *state;
  int i, j, k, *start, *stop;

  /*
   *  open up the filename (if it exists, otherwise prompt the user)
   *  and determine whether it is a AMBER prmtop or CHARMM PSF file
   */

  if (filename == NULL) {
    filename = promptToOpenFile(&fp, "", "r", 
             "Input the name of an AMBER prmtop or CHARMM PSF: ");
  } else {
    if ( openFile(&fp, filename, "r") == 0 )
      error(ROUTINE, "Attempting to open parameter/topology file %s",
	    filename);
  }
  buffer = (char *) safe_malloc(sizeof(char)* BUFFER_SIZE);
  if (fgets(buffer, BUFFER_SIZE, fp) == NULL)
    error(ROUTINE, "Attempting to read parameter/topology file");

  if ( strcmp(buffer, "PSF \n") == 0 ||
       strcmp(buffer, "PSF\n") == 0 )
    state = loadCharmmPSF(fp, 1);
  else {
    safe_freopen(fp);

    if ( parm != NULL ) 
      clearParm( parm );
    parm = safe_malloc( sizeof( Parm ) );
    initializeParm(parm);
    parm->filename = filename;
    parm->fp = fp;
    readParm();

    state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
    INITIALIZE_ptrajState(state);

    state->atoms = parm->NTOTAT;
    state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
    state->masses = (double *) safe_malloc(sizeof(double) * state->atoms);
    state->charges = (double *) safe_malloc(sizeof(double) * state->atoms);
    for (i=0; i < state->atoms; i++) {
      strcpy(state->atomName[i], parm->atom[i].igraph);
      state->masses[i] = parm->atom[i].amass;
      state->charges[i] = parm->atom[i].chrg / CHARGE_TO_KCALS;
    }

    state->residues = parm->NTOTRS;
    state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
    for (i=0; i < state->residues; i++) {
      strcpy(state->residueName[i], parm->residue[i].labres);
    }
    state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
    for (i=0; i <= state->residues; i++) {
      state->ipres[i] = parm->residue[i].ipres;
    }

    state->box[0] = 0.0;
    state->box[1] = 0.0;
    state->box[2] = 0.0;
    state->box[3] = 90.0;
    state->box[4] = 90.0;
    state->box[5] = 90.0;
    state->maxFrames = 0;

    state->IFBOX = parm->IFBOX;
    if ( state->IFBOX ) {

      state->molecules = parm->box->nspm;
      state->moleculeInfo = (int *)
	safe_malloc( sizeof(int) * state->molecules );
      for (i=0; i < state->molecules; i++) {
	state->moleculeInfo[i] = parm->box->nsp[i];
      }

      state->box[0] = parm->box->box[0];
      state->box[1] = parm->box->box[1];
      state->box[2] = parm->box->box[2];
      if (parm->box->beta != 90.0) {

	if (parm->box->beta > 109.47 && parm->box->beta < 109.48) {

	  state->box[3] = 2.0*acos(1.0/sqrt(3.0))*RADDEG;
	  state->box[4] = state->box[3];
	  state->box[5] = state->box[3];
	  fprintf(stdout, " Setting box to be an exact truncated octahedron, angle is %f\n",
		  state->box[3]);

	} else if (parm->box->beta == 60.0) {

	  fprintf(stdout, " Setting box to be a rhombic dodecahedron, i.e. alpha=gamma=60.0, beta=90.0\n");
	  state->box[3] = 60.0;
	  state->box[4] = 90.0;
	  state->box[5] = 60.0;
	}

      }
    }

    state->solventMolecules = 0;
    state->solventAtoms = 0;
    state->solventMask = NULL;
    state->solventMoleculeStart = NULL;
    state->solventMoleculeStop = NULL;


    start = (int *) safe_malloc(sizeof(int) * state->atoms);
    stop  = (int *) safe_malloc(sizeof(int) * state->atoms);
    state->solventMask = (int *) safe_malloc(sizeof(int) * state->atoms);
    for (i=0; i < state->atoms; i++)
      state->solventMask[i] = 0;

    if ( ! parm->IFBOX) {
      /*
       *  treat all the molecules starting with parm->box->nspsol
       *  as solvent molecules IF there is box information
       */

      j = 0;
      for (i=0; i < state->molecules; i++) {

	if (i+1 >= parm->box->nspsol) {
	  /*
	   *  add this molecule to the solvent list
	   */

	  state->solventAtoms += state->moleculeInfo[i];

	  for (k = j; k < j+state->moleculeInfo[i]; k++)
	    state->solventMask[k] = 1;
	  start[state->solventMolecules] = j;
	  stop[ state->solventMolecules] = j+state->moleculeInfo[i];
	  state->solventMolecules++;
	}

	j += state->moleculeInfo[i];
      }


    } else {
      /*
       *  treat all the residues named "WAT " as solvent molecules
       */

      for (i=0; i < state->residues; i++) {
	if (strcmp("WAT ", state->residueName[i]) == 0) {
	  /*
	   * add this residue to the list of solvent
	   */
	  j = state->ipres[i+1]-state->ipres[i];

	  state->solventAtoms += j;
	  start[state->solventMolecules] = state->ipres[i]-1;
	  stop[ state->solventMolecules] = state->ipres[i+1]-1;
	  state->solventMolecules++;

	  for (k=state->ipres[i]-1; k < state->ipres[i+1]-1; k++)
	    state->solventMask[k] = 1;
	}
      }
    }

    /*
    state->solventMoleculeStart = (int *) 
      safe_malloc(sizeof(int) * state->solventMolecules);
    state->solventMoleculeStop = (int *) 
      safe_malloc(sizeof(int) * state->solventMolecules);
    for (i=0; i < state->solventMolecules; i++) {
      state->solventMoleculeStart[i] = start[i];
      state->solventMoleculeStop[i]  = stop[i];
    }

    safe_free(stop);
    safe_free(start);
    */

    state->solventMoleculeStart = start;
    state->solventMoleculeStop  = stop;

  }

  *statep = state;
}



#undef  ROUTINE
#define ROUTINE printCoordinateInfo

   void
printCoordinateInfo( void *entry )
{
  coordinateInfo *p;
  charmmTrajectoryInfo *charmmTraj;
  byte u;
  int i;

  p = (coordinateInfo *) entry;
  if (p == NULL) {
    fprintf(stdout, "  *** No coordinate file of this type has been defined\n");
    return;
  }

  switch ( p->type ) {

  case COORD_BINPOS:
    fprintf(stdout, "  File (%s) is a BINPOS file ", p->filename);
    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
	fprintf(stdout, " with %i sets (processing only %i)\n", p->stop, i);
      else
	fprintf(stdout, " with %i sets\n", i);
    } else {
      fprintf(stdout, "\n");
    }
    break;

  case COORD_PDB:
    fprintf(stdout, "  File (%s) is a PDB file", p->filename);
    if (p->option1 == 1)
      fprintf(stdout, ": AMBER charges and radii to occupancy and temp factor columns");
    else if (p->option1 == 2)
      fprintf(stdout, ": AMBER charges and PARSE radii to occupancy and temp factor columns");
    fprintf(stdout, "\n");

    break;

  case COORD_AMBER_TRAJECTORY:
    fprintf(stdout, "  File (%s) is an AMBER trajectory%s",
	    p->filename, (p->isBox ? " (with box info)" : ""));
    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
	fprintf(stdout, " with %i sets (processing only %i)\n", p->stop, i);
      else
	fprintf(stdout, " with %i sets\n", i);
    } else {
      fprintf(stdout, "\n");
    }
    break;

  case COORD_AMBER_NETCDF:

    fprintf(stdout, "  File (%s) is a NetCDF AMBER trajectory%s%s",
	    p->filename, 
	    (p->isBox ? " with box info" : ""),
	    (p->isVelocity ? " and velocities" : ""));

    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
	fprintf(stdout, " with %i sets (processing only %i)\n", p->stop, i);
      else
	fprintf(stdout, " with %i sets\n", i);
    } else {
      fprintf(stdout, "\n");
    }
    if (prnlev > 2) {
      if (p->title != NULL)
	fprintf(stdout, "    title:        \"%s\"\n", p->title);
      if (p->application != NULL) 
	fprintf(stdout, "    application:  \"%s\"\n", p->application);
      if (p->program != NULL) 
	fprintf(stdout, "    program:      \"%s\"\n", p->program);
      if (p->version != NULL) 
	fprintf(stdout, "    version:      \"%s\"\n", p->version);
    }

    break;

  case COORD_CHARMM_TRAJECTORY:

    charmmTraj = (charmmTrajectoryInfo *) p->info;
    fprintf(stdout,
	    "  File (%s) is a CHARMM trajectory in %s endian binary format\n",
	    p->filename, (charmmTraj->byteorder ? "little" : "big"));
    fprintf(stdout, "%s", (p->isBox ? "    with box information " : "    "));
    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
	fprintf(stdout, "representing %i sets (processing only %i)\n", p->stop, i);
      else
	fprintf(stdout, "representing %i sets\n", p->stop);
    } else
      fprintf(stdout, "\n");

    /*
     *  NFILE  -- number of coordinate sets in file
     *  ISTEP  -- number of previous dynamics steps
     *  NINTV  -- frequency for saving coordinates
     *  NSAVC  -- number of steps for creation run
     */  
    printf("  NFILE = %8i ISTEP = %8i NINTV = %8i NSAVC = %8i NSAVV = %8i\n", 
	   charmmTraj->icntrl[0], charmmTraj->icntrl[1],
	   charmmTraj->icntrl[2], charmmTraj->icntrl[3],
	   charmmTraj->icntrl[4]);

    u.i = charmmTraj->icntrl[9];
    printf("  NDEGF = %8i NFREA = %8i DELTA = %8.4f QCRYS = %8i QDIM4 = %8i\n", 
	   charmmTraj->icntrl[7], charmmTraj->icntrl[8], u.f,
	   charmmTraj->icntrl[10], charmmTraj->icntrl[11]);
    printf("  VERNU = %8i NATREC= %8i NFREAT = %8i\n", 
	   charmmTraj->icntrl[19], charmmTraj->natrec, charmmTraj->nfreat);
    printf("  Dumping the title:\n");
    printStack(&charmmTraj->titleStack, printString, NULL);

    break;

  case COORD_AMBER_RESTART:
    fprintf(stdout, "  File (%s) is an AMBER restart file ", p->filename);
    if ( p->isBox && p->isVelocity )
      fprintf(stdout, "with box and velocity information\n");
    else if (p->isBox)
      fprintf(stdout, "with box information\n");
    else if (p->isVelocity)
      fprintf(stdout, "with velocity information\n");
    else
      fprintf(stdout, "\n");
    break;
  }

  if (prnlev > 4) {
    if ( p->file == NULL ) {
      fprintf(stdout, "    [the FILE is currently closed]\n");
    } else if (p->file == stdin) {
      fprintf(stdout, "    [the FILE is standard input]\n");
    } else if (p->file == stdout) {
      fprintf(stdout, "    [the FILE is standard output]\n");
    } else if (p->file == stderr) {
      fprintf(stdout, "    [the FILE is standard error]\n");
    }
  }

  if ( p->mask != NULL )
    fprintf(stdout, "    [The file has an active atom mask]\n");

}


#undef  ROUTINE
#define ROUTINE "ptrajCleanup()"

   void
ptrajCleanup()
{
  actionInformation *a;
  analyzeInformation *an;
  coordinateInfo *f;
  int *mask;

  /*
   *  Clean up the INPUT files, transformFileStack
   */

  while (transformFileStack != NULL  &&
	 (f = (coordinateInfo *) popStack(&transformFileStack)) != NULL ) {
    
    safe_free(f->filename);
    safe_free(f->info);
    safe_free(f->mask);
    safe_free(f->x);
    safe_free(f->y);
    safe_free(f->z);
    INITIALIZE_coordinateInfo(f);
    safe_free(f);
    
  }
  transformFileStack = NULL;

  /*
   *  Clean up the ACTIONS
   */

  while (transformActionStack != NULL &&
	 (a = (actionInformation *) popStack(&transformActionStack)) != NULL) {

       /*
        *  Free any associated mask
        */

    if (a->mask) 
      safe_free(a->mask);

       /*
        *  Free any associated state information
        */

    if (a->state != NULL) {
      ptrajClearState(&a->state);
      a->state = NULL;
    }

       /*
        *  Clean up any of the complex arguments as necessary; this
        *  is done by the associated action function in the PTRAJ_CLEANUP
        *  mode
        */

    if (a->type != TRANSFORM_TRANSFORM &&
	a->type != TRANSFORM_NOOP) {

      if (a->fxn != NULL) {
	a->fxn(a, NULL, NULL, NULL, NULL, PTRAJ_CLEANUP);
      }
      INITIALIZE_actionInformation(a);
    }

       /*
        *  Free up the action
        */

    safe_free(a);

  }

  transformActionStack = NULL;

  /*
   *  Clean up the ANALYZEs
   */

  while (transformAnalyzeStack != NULL &&
	 (an = (analyzeInformation *) popStack(&transformAnalyzeStack)) != NULL) {

       /*
        *  Clean up any of the complex arguments as necessary; this
        *  is done by the associated action function in the PTRAJ_CLEANUP
        *  mode
        */

    if (an->fxn != NULL) {
      an->fxn(an, NULL, PTRAJ_CLEANUP);
    }
    INITIALIZE_analyzeInformation(an);

       /*
        *  Free up the action
        */
    safe_free(an);

  }

  transformAnalyzeStack = NULL;

  /*
   *  Free up outInfo
   */

  if (globalOutInfo != NULL) {
    safe_free(globalOutInfo->filename);
    globalOutInfo->filename = NULL;
    safe_free(globalOutInfo->mask);
    globalOutInfo->mask = NULL;
    safe_free(globalOutInfo->info);
    globalOutInfo->info = NULL;
    safe_free(globalOutInfo);
    globalOutInfo = NULL;
  }

  /*
   *  Free up referenceInfo
   */

  while (transformReferenceStack != NULL  &&
	 (f = (coordinateInfo *) popStack(&transformReferenceStack)) != NULL ) {
    
    safe_free(f->filename);
    safe_free(f->info);
    safe_free(f->mask);
    safe_free(f->x);
    safe_free(f->y);
    safe_free(f->z);
    INITIALIZE_coordinateInfo(f);
    safe_free(f);
    
  }
  transformReferenceStack = NULL;
  referenceInfo = NULL;

  /*
   *  Clean up hbond donor/acceptor information
   */

  while (hbondDonorStack != NULL &&
	 (mask = (int *) popStack(&hbondDonorStack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorStack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorStack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorH1Stack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorH1Stack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorH2Stack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorH2Stack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorH3Stack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorH3Stack)) != NULL) {
    safe_free(mask);
  }

}


#undef  ROUTINE
#define ROUTINE "ptrajPreprocessInputCoordinates()"

   int
ptrajPreprocessInputCoordinates(coordinateInfo *currentCoordinateInfo)
{
  charmmTrajectoryInfo **charmmTrajp;
  netcdfTrajectoryInfo *NCInfo;
  char buffer[BUFFER_SIZE];
  int err;


     /*
      *  open up the file
      */
#ifdef BINTRAJ
  if ( currentCoordinateInfo->type == COORD_AMBER_NETCDF ) {

    NCInfo = (netcdfTrajectoryInfo *) currentCoordinateInfo->info;
    err = nc_open(currentCoordinateInfo->filename, NC_NOWRITE, &NCInfo->ncid);
    if (err != NC_NOERR) error(ROUTINE, "%s", nc_strerror(err));

  } else 
#endif
    if ( ! openFile(&currentCoordinateInfo->file,
		  currentCoordinateInfo->filename, "r") ) {

    fprintf(stdout, "WARNING in ptrajPreprocessInputCoordinates(): Error on opening\n");
    fprintf(stdout, "input coordinate file (%s)\n", 
	    currentCoordinateInfo->filename);
    return 1;
  }

     /*
      *  preprocess the input coordinates as necessary 
      *  (i.e. to remove titles, etc.)
      */

  switch(currentCoordinateInfo->type) {
  case COORD_PDB:
    fprintf(stdout, "\nProcessing PDB file %s\n",
	    currentCoordinateInfo->filename);
    break;
  case COORD_BINPOS:
    fprintf(stdout, "\nProcessing BINPOS file %s\n",
	    currentCoordinateInfo->filename);
    break;
  case COORD_AMBER_RESTART:
    fprintf(stdout, "\nProcessing AMBER restart file %s\n",
	    currentCoordinateInfo->filename);
    break;
  case COORD_AMBER_TRAJECTORY:
    fprintf(stdout, "\nProcessing AMBER trajectory file %s\n",
	    currentCoordinateInfo->filename);
    break;
  case COORD_CHARMM_TRAJECTORY:
    fprintf(stdout, "\nProcessing CHARMM trajectory file %s\n",
	    currentCoordinateInfo->filename);
    break;
  }
      
  switch ( currentCoordinateInfo->type ) {

  case COORD_AMBER_TRAJECTORY:

    /*
     *  we need to read the title line to set up for processing
     *  by readAmberTrajectory()
     */

    if ( fgets(buffer, BUFFER_SIZE, currentCoordinateInfo->file) == NULL ) {
      fprintf(stdout, "WARNING: Error on processing the title from the AMBER\n");
      fprintf(stdout, "trajectory (%s)\n", currentCoordinateInfo->filename);
      return 1;
    }
    break;

  case COORD_BINPOS:

	if ( openbinpos( currentCoordinateInfo->file) < 0 ){
      fprintf(stdout, "Error on opening BINPOS file %s\n",
         currentCoordinateInfo->filename);
      return 1;
    }
    break;

  case COORD_CHARMM_TRAJECTORY:

    /*
     *  read in all the header information; this is done by calling with a
     *  negative value for the current set...
     */

    charmmTrajp = (charmmTrajectoryInfo **) safe_malloc(sizeof(charmmTrajectoryInfo *));
    *charmmTrajp = NULL;
    readCharmmTrajectory(currentCoordinateInfo->file, charmmTrajp,
			 NULL, NULL, NULL, NULL, -1);


    break;

  case COORD_UNKNOWN:

    fprintf(stdout, 
	    "WARNING: Attempting to process a coordinate file of unknown type\n");
    fprintf(stdout, "in %s, ignoring...\n", ROUTINE);
    return 1;
  }
  
  return 0;

}    

#undef  ROUTINE
#define ROUTINE "ptrajProcessInputCoordinates()"     

   void
ptrajProcessInputCoordinates(coordinateInfo *currentCoordinateInfo, 
			     ptrajState *state, 
			     double *X, double *Y, double *Z,
			     double *box, int set,
			     int *readCoordinates, int *processCoordinates)
{
  pdb_record *pdb;
  float *binpos;
  float fbox[6];
  int i,j,n_atoms,eof,err;
  charmmTrajectoryInfo **charmmTrajp;
  netcdfTrajectoryInfo *NCInfo;
  size_t start[3], count[3];
  char xyz[3];
  float time;

  /*
   *  NOTE: it is assumed that the box information is set to valid values on entry
   */



  /*
   *  READ IN THE CURRENT FILE, ONE SET AT A TIME
   */
  switch( currentCoordinateInfo->type ) {

  case COORD_PDB:

    i = loadPdb(currentCoordinateInfo->file, &pdb);
    j = getCoordinatesFromPdb(pdb, X, Y, Z);
    if (j != state->atoms) {

      /*
       *  on error, print warning but do not stop processing of this set.  This
       *  allows the specification of multiple reference sets
       */
      fprintf(stdout, 
	      "\nWARNING in %s: Unexpected number of atoms\n", ROUTINE);
      fprintf(stdout, " encountered when reading PDB!\n");
      fprintf(stdout, "  coordinates read for %i atoms, expecting %i\n", j, state->atoms);
      
    }


    /*
     *  We assume that a PDB only contains ONE set of coordinates!!!
     */
    safe_fclose(currentCoordinateInfo->file);
    currentCoordinateInfo->file = NULL;
    *readCoordinates = 0;
    *processCoordinates = 1;
    break;
  
  case COORD_AMBER_RESTART:

    if ( getCoordinatesFromRestart(currentCoordinateInfo->file, 
				   X, Y, Z,
				   &box[0], &box[1], &box[2],
				   &box[3], &box[4], &box[5]) 
	 != state->atoms ) {
      fprintf(stdout, "WARNING in ptrajProcessInputCoordinates(): Unexpected\n");
      fprintf(stdout, "number of atoms encountered reading restrt file\n");
      ptrajCleanup();
      *readCoordinates = 0;
      *processCoordinates = 0;
      return;
    }

    /*
     *  get default box information if necessary
     */
    if ( state->IFBOX && (box[0] == 0.0 && box[1] == 0.0 && box[2] == 0.0) ) {
      box[0] = state->box[0];
      box[1] = state->box[1];
      box[2] = state->box[2];
    }

    if ( state->IFBOX && (box[3] == 0.0 && box[4] == 0.0 && box[5] == 0.0) ) {
      box[3] = state->box[3];
      box[4] = state->box[4];
      box[5] = state->box[5];
    }

    /*
     *  AMBER restrt files only contain a single set of coordinates
     */
    safe_fclose(currentCoordinateInfo->file);
    currentCoordinateInfo->file = NULL;
    *readCoordinates = 0;
    *processCoordinates = 1;
    break;

  case COORD_AMBER_TRAJECTORY:

    /*
     *  check to see if we've already loaded enough of the current
     *  trajectory file
     */

    if ( set > currentCoordinateInfo->stop ) {
      *readCoordinates = 0;
      
    } else {

      /*
       *  read in next frame
       */

      *readCoordinates = readAmberTrajectory(currentCoordinateInfo->file,
					     state->atoms, X, Y, Z, box, set,
					     currentCoordinateInfo->isBox);
      
      if ( state->IFBOX && (box[0] == 0.0 && box[1] == 0.0 && box[2] == 0.0) ) {
	box[0] = state->box[0];
	box[1] = state->box[1];
	box[2] = state->box[2];
      }

      if ( state->IFBOX && (box[3] == 0.0 && box[4] == 0.0 && box[5] == 0.0) ) {
	box[3] = state->box[3];
	box[4] = state->box[4];
	box[5] = state->box[5];
      }

      *processCoordinates = 1;
    }

    if ( *readCoordinates == 0 ) {
      *processCoordinates = 0;
      safe_fclose(currentCoordinateInfo->file);
      currentCoordinateInfo->file = NULL;
    }

    break;

  case COORD_AMBER_NETCDF:

    /*
     *  Prepare to load up a frame of NetCDF data and
     *  assume initially that we have successfully
     *  loaded up a set of coordinates
     */

#ifdef BINTRAJ
    NCInfo = (netcdfTrajectoryInfo *) currentCoordinateInfo->info;
    *readCoordinates = 1;
    *processCoordinates = 1;

      /*
       *  Check to see if we've loaded enough frames
       */

    if ( set > currentCoordinateInfo->stop ) {

      *readCoordinates = 0;
      *processCoordinates = 0;
  
    } else {

        /*
         *  unnecessary (but useful) sanity check: we should never trigger this since 
         *  this routine should not be called if coordinates are not present
         */

      if (NCInfo->coordinateVID == 0) {
	*readCoordinates = 0;
	fprintf(stdout, "Expecting coordinates in NetCDF file %s, none present...\n",
		currentCoordinateInfo->filename);
      }

      if ( *readCoordinates ) {
	if ( X == NULL || Y == NULL || Z == NULL ) 
	  error(ROUTINE, "coordinate arrays are NULL\n");

	/*
	 *  allocate space for coordinates from NetCDF file if necessary
	 */
	if ( NCInfo->R == NULL ) {
	  NCInfo->R = (float *) safe_malloc( sizeof(float) * state->atoms * 3 );
	}
    
	/*
	 *  scan in coordinates for this frame
	 *
	 *  !! NOTE: C arrays are opposite the F90/Fortran so reverse
	 *     the array bounds on start and count.
	 */


	start[0] = 0;
	count[0] = 3;
	xyz[0] = (char) 0;
	err = nc_get_vara_text(NCInfo->ncid, NCInfo->spatialVID, start, count, xyz);
	if (err != NC_NOERR) {
	  fprintf(stdout, "Yikes!  We should see X, Y and Z chars in the spatial VID: %s\n",
		  nc_strerror(err));
	} else if ( xyz[0] != 'x' || xyz[1] != 'y' || xyz[2] != 'z' ) {
	  fprintf(stdout, "Whoa Nellie & Cripes!!!  We should see 'x', 'y' and 'y' characters\n");
	  fprintf(stdout, "in the spatial variables of the AMBER NetCDF trajectory file unless\n");
	  fprintf(stdout, "someone altered the spec.  We see '%c', '%c' and '%c'\n",
		  xyz[0], xyz[1], xyz[2]);
	}

	start[0] = set-1;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = state->atoms;
	count[2] = 3;
      
	err = nc_get_vara_float(NCInfo->ncid, NCInfo->coordinateVID, start, count, NCInfo->R);
	if (err != NC_NOERR) {
	  fprintf(stdout, "Error detected upon reading NetCDF coordinates on set %i: %s", 
		  set, nc_strerror(err));
	  *readCoordinates = 0;
	  *processCoordinates = 0;
	} else {
	
	  for (i = 0, j = 0; i < state->atoms; i++, j += 3) {
	    X[i] = NCInfo->R[ j   ];
	    Y[i] = NCInfo->R[ j+1 ];
	    Z[i] = NCInfo->R[ j+2 ];
	  }
	}
      }

      if ( *readCoordinates && state->IFBOX) {

	start[0] = set-1;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = 3;
	count[2] = 0;

	err = nc_get_vara_double(NCInfo->ncid, NCInfo->cellLengthVID, start, count, box);
	if (err != NC_NOERR) 
	  error(ROUTINE, "Getting cell lengths: %s", nc_strerror(err));
	err = nc_get_vara_double(NCInfo->ncid, NCInfo->cellAngleVID, start, count, &box[3]);
	if (err != NC_NOERR)
	  error(ROUTINE, "Getting cell angles: %s", nc_strerror(err));

	/*
	fprintf(stdout, "\nGOT BOX INFO: %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		box[0], box[1], box[2], box[3], box[4], box[5]);
	*/


      }

      /*
       *  get time
       */
      start[0] = set-1;
      count[0] = 1;
      err = nc_get_vara_float(NCInfo->ncid, NCInfo->timeVID, start, count, &time);
      if (err != NC_NOERR)
	fprintf(stdout, "Warning in %s, getting the time from the NetCDF file: %s\n",
		ROUTINE, nc_strerror(err));

      /*
      fprintf(stdout, "\nREAD IN TIME %10.5f\n", time);
      */


    }
    
    if ( *readCoordinates == 0 ) {
      *processCoordinates = 0;

      err = nc_close(NCInfo->ncid);
      if (err) fprintf(stdout, "Error closing NetCDF file %s (nc_close), error: %s\b",
		       currentCoordinateInfo->filename, nc_strerror(err));

      FREE_netcdfTrajectoryInfo(NCInfo);
      currentCoordinateInfo->info = NULL;
      currentCoordinateInfo->file = NULL;

    }
#endif

    break;

  case COORD_CHARMM_TRAJECTORY:

    charmmTrajp = (charmmTrajectoryInfo ** ) safe_malloc(sizeof(charmmTrajectoryInfo *));
    *charmmTrajp = (charmmTrajectoryInfo *) currentCoordinateInfo->info;

    *readCoordinates = 
      readCharmmTrajectory(currentCoordinateInfo->file, charmmTrajp,
			   X, Y, Z, box, set);
    break;

  case COORD_BINPOS:

    n_atoms = state->atoms;
    binpos = safe_malloc(sizeof(float) * n_atoms * 3);
    *readCoordinates = readbinpos( currentCoordinateInfo->file, 
				   &n_atoms, binpos, &eof);
    j=0;
    for (i=0; i<n_atoms; i++) {
		X[i] = binpos[j];
		Y[i] = binpos[j+1];
		Z[i] = binpos[j+2];
        j += 3;
    }
    safe_free(binpos);

    if ( (eof == 1) || (set > currentCoordinateInfo->stop) ) {
	*readCoordinates = 0;
    }
    if ( *readCoordinates == 0 ) {
      *processCoordinates = 0;
      safe_fclose(currentCoordinateInfo->file);
      currentCoordinateInfo->file = NULL;
    } else {
      *processCoordinates = 1;
    }
    break;

  default:
    *readCoordinates = 0;
    *processCoordinates = 0;
    fprintf(stdout, "WARNING in ptrajProcessInputCoordinates(): Attempting to process\n");
    fprintf(stdout, "a file of unknown type (%s)", currentCoordinateInfo->filename);
    fprintf(stdout, "Ignoring this file...\n");
    
  }
}

   int
lesSize( int atoms )
{
    int natomCL = 0, i=0;
    
    for( i=0; i < atoms; i++ )
    {
        if( parm->lescnum[i] == 0 || parm->lescnum[i] == 1 ) 
        {
            natomCL++;
        }
    }
    
    return natomCL;
}
   int*
lesMask( int atoms )
{
    int i;
    int* mask;

    mask = (int *) safe_malloc( sizeof(int) * atoms );

    for( i=0; i < atoms; i++ )
    {    
        if( parm->lescnum[i] == 0 || parm->lescnum[i] == 1 )
        {
            mask[i] = 1;
        }
        else
        { 
            mask[i] = 0;
        }
    }

    return mask;
}

   void
lesSplit( int atoms, double *X, double *Y, double *Z, 
          int icopy, double* xrep, double* yrep, double* zrep )
{
    int ia=0, i=0;
    
    for( i=0; i < atoms; i++ )
    {
        if( parm->lescnum[i] == 0 || parm->lescnum[i] == icopy + 1 )
        {
            xrep[ia] = X[i];
            yrep[ia] = Y[i];
            zrep[ia] = Z[i];
            ia++;
        }
    }
}

   void
lesAverage( int atoms, int nlescopy, double *X, double *Y, double *Z, 
	    double* xrep, double* yrep, double* zrep )
{
    int ia=0, i=0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    
    for( i=0; i < atoms; i++ )
    {
        if( parm->lescnum[i] == 0 )
        {
            xrep[ia] = X[i];
            yrep[ia] = Y[i];
            zrep[ia] = Z[i];
            ia++;
        }
        else
        {
            if( parm->lescnum[i] == 1 )
            {
                sum_x = 0.0;
                sum_y = 0.0;
                sum_z = 0.0;
            }
            
            sum_x += X[i];
            sum_y += Y[i];
            sum_z += Z[i];

            if( parm->lescnum[i] == nlescopy )
            {
                xrep[ia] = sum_x / nlescopy;
                yrep[ia] = sum_y / nlescopy;
                zrep[ia] = sum_z / nlescopy;
                ia++;
            }
        }
    }
}



#undef  ROUTINE
#define ROUTINE "ptrajOutputCoordinates()"

   void
ptrajOutputCoordinates(coordinateInfo *outInfo, ptrajState *state,
		       int set, int append, int first, int last, int atoms,
		       double *X, double *Y, double *Z, double *box)
{
  /*
   *  outInfo:  a pointer to the coordinateInfo file structure for the file to be output
   *  state:    the current ptrajState information
   *  set:      the current set or frame being processed
   *  append:   a flag (if non-zero) specifying that the outInfo file should be appended
   *  first:    a flag (if non-zero) specifying that this is the first call for this file
   *  last:     a flag (if non-zero) specifying that this is the last call for this file
   *  atoms:    the number of atoms
   *  X, Y, Z:  the coordinates
   *  box:      the box information (if not NULL)
   */
 
  char buffer[BUFFER_SIZE];
  pdb_record *pdb;
  netcdfTrajectoryInfo *NCInfo;
  int i, j, err, status;
#ifdef BINTRAJ
  int dimensionID[NC_MAX_VAR_DIMS];
  size_t start[3], count[3];
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a', 
		   'b', 'e', 't', 'a', ' ', 
		   'g', 'a', 'm', 'm', 'a' };
#endif
  float time;

  /*
   *  Four basic tasks are performed by this routine that enable output of the
   *  processed coordinates into different trajectory formats as specified by the
   *  (coordType) outInfo->type values which are defined in trajectory.h.
   *
   *  If you add a new coordType, implement code for each of these functions!
   *
   *    If (last) {
   *        (1) POST-PROCESS AND CLOSE FILE
   *        return;
   *    }
   *
   *    If (first) {
   *        (2) OPEN FILE AND PRE-PROCESS AND/OR SET UP TO APPEND [AND DO NOT return]
   *    }
   *    (3) DUMP CURRENT FRAME/SET
   */

 

    /*
     *  Perform final clean-up (closing) of trajectory file
     */
  if (last) {

    switch(outInfo->type) {

    case COORD_CHARMM_TRAJECTORY:
    case COORD_AMBER_TRAJECTORY:
    case COORD_BINPOS:

      safe_fclose(outInfo->file);
      outInfo->file = NULL;
      return;

#ifdef BINTRAJ
    case COORD_AMBER_NETCDF:
      
      NCInfo = (netcdfTrajectoryInfo *) outInfo->info;
      err = nc_close(NCInfo->ncid);
      if (err) fprintf(stdout, "Error closing NetCDF file %s (nc_close): %s\n",
		       outInfo->filename, nc_strerror(err));
      return;
#endif

    default:
      return;
    }

  }

    /*
     *  Do recursive calls as necessary for outputing LES trajectories
     */
  if ( outInfo->les_action != LES_NONE && outInfo->les_status == LES_READY ) {
    outInfo->les_status = LES_DONE;
    int natomCL = lesSize( atoms );
  
    double* xrep = safe_malloc( sizeof( double ) * natomCL );
    double* yrep = safe_malloc( sizeof( double ) * natomCL );
    double* zrep = safe_malloc( sizeof( double ) * natomCL );

    if ( outInfo->les_action == LES_SPLIT ) {
      int icopy=0;
      for ( icopy=0; icopy < outInfo->nlescopy; icopy++ ) {
	if ( icopy > 0 && first != 0 ) first = 0;
	lesSplit( atoms, X, Y, Z, icopy, xrep, yrep, zrep );

	ptrajOutputCoordinates( outInfo, state, (set-1) * outInfo->nlescopy + icopy + 1, 
				append, first, last, natomCL, xrep, yrep, zrep, box );
      }
    } else {
      assert( outInfo->les_action == LES_AVERAGE );
          
      lesAverage( atoms, outInfo->nlescopy, X, Y, Z, xrep, yrep, zrep );

      ptrajOutputCoordinates( outInfo, state, set, append, first, last, natomCL, xrep, yrep, zrep, box );
    }

    safe_free( xrep );
    safe_free( yrep );
    safe_free( zrep );

    outInfo->les_status = LES_READY;

    return;
  }

  if (first) {

    switch( outInfo->type ) {
	    
    case COORD_PDB:
      /*
       *  file handling is done per frame
       */
      break;

    case COORD_AMBER_NETCDF:

       /*
        *  On first call, preprocess to create the NetCDF file and to dump header information,
        */

      if (append)
	error(ROUTINE, "Appending to NetCDF files is not yet implemented\n");

#ifdef BINTRAJ
      NCInfo = (netcdfTrajectoryInfo *) safe_malloc(sizeof(netcdfTrajectoryInfo));
      INITIALIZE_netcdfTrajectoryInfo(NCInfo);

        /*
         *  create the NetCDF file
         */
      err = nc_create(outInfo->filename, NC_64BIT_OFFSET, &NCInfo->ncid);
      if (err != NC_NOERR) {
	fprintf(stdout, "WARNING in ptrajOutputCoordinates(): Error opening\n");
	fprintf(stdout, "NetCDF output coordinate file (%s): %s\n", 
		outInfo->filename, nc_strerror(err));
	outInfo->type = COORD_UNKNOWN;
	break;
      }

        /*
         *  define the global dimensions of the NetCDF file
         */
      err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_FRAME, NC_UNLIMITED, &NCInfo->frameDID);
      if (err != NC_NOERR)
	fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF frame dimension ID: %s\n",
		nc_strerror(err));

      err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_SPATIAL, 3, &NCInfo->spatialDID);
      if (err != NC_NOERR)
	fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF spatial dimension ID: %s\n",
		nc_strerror(err));
      
      err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_ATOM, atoms, &NCInfo->atomDID);
      if (err != NC_NOERR)
	fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF atom dimension ID: %s\n",
		nc_strerror(err));

      err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_LABEL, AMBER_NETCDF_LABELLEN, &NCInfo->labelDID);
      if (err != NC_NOERR)
	fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF label dimension ID: %s\n",
		nc_strerror(err));

      err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, 3, &NCInfo->cell_spatialDID);
      if (err != NC_NOERR)
	fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF cell_spatial dimension ID: %s\n",
		nc_strerror(err));

      err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, 3, &NCInfo->cell_angularDID);
      if (err != NC_NOERR)
	fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF cell_angular dimension ID: %s\n",
		nc_strerror(err));

        /*
         *  put global attributes
         */

      netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "title", outInfo->title);
      netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "application", outInfo->application);
      netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "program", outInfo->program);
      netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "programVersion", outInfo->version);
      netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "Conventions", "AMBER");
      netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "ConventionVersion", "1.0");

        /*
         *  handle definition of non-optional variables
         */

      dimensionID[0] = NCInfo->spatialDID;
      err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->spatialVID);
      if (err != NC_NOERR) fprintf(stdout, "NetCDF defining spatial VID: %s\n", nc_strerror(err));

      dimensionID[0] = NCInfo->frameDID;
      err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_TIME, NC_FLOAT, 1, dimensionID, &NCInfo->timeVID);
      if (err != NC_NOERR) fprintf(stdout, "NetCDF defining time VID: %s\n", nc_strerror(err));
      netcdfPutAttributeText(NCInfo->ncid, NCInfo->timeVID, "units", "picosecond");
      if (err != NC_NOERR) fprintf(stdout, "NetCDF putting time VID units: %s\n", nc_strerror(err));

      dimensionID[0] = NCInfo->frameDID;
      dimensionID[1] = NCInfo->atomDID;
      dimensionID[2] = NCInfo->spatialDID;
      err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_COORDS, NC_FLOAT, 3, dimensionID, &NCInfo->coordinateVID);
      if (err != NC_NOERR) fprintf(stdout, "NetCDF defining coordinate VID: %s\n", nc_strerror(err));
      netcdfPutAttributeText(NCInfo->ncid, NCInfo->coordinateVID, "units", "angstrom");
      if (err != NC_NOERR) fprintf(stdout, "NetCDF putting coordinate VID units: %s\n", nc_strerror(err));

      dimensionID[0] = NCInfo->cell_spatialDID;
      err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->cell_spatialVID);
      if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell spatial VID: %s\n", nc_strerror(err));

      dimensionID[0] = NCInfo->cell_angularDID;
      dimensionID[1] = NCInfo->labelDID;
      err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, NC_CHAR, 2, dimensionID, &NCInfo->cell_angularVID);
      if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell spatial VID: %s\n", nc_strerror(err));


      if (outInfo->isBox) {

	dimensionID[0] = NCInfo->frameDID;
	dimensionID[1] = NCInfo->cell_spatialDID;
	err = nc_def_var(NCInfo->ncid, "cell_lengths", NC_DOUBLE, 2, dimensionID, &NCInfo->cellLengthVID);
	if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell length VID: %s\n", nc_strerror(err));
	netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellLengthVID, "units", "angstrom");
	if (err != NC_NOERR) fprintf(stdout, "NetCDF putting cell length VID units: %s\n", nc_strerror(err));

	dimensionID[1] = NCInfo->cell_angularDID;
	err = nc_def_var(NCInfo->ncid, "cell_angles", NC_DOUBLE, 2, dimensionID, &NCInfo->cellAngleVID);
	if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell angle VID: %s\n", nc_strerror(err));
	netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellAngleVID, "units", "degree");

	if (err != NC_NOERR) fprintf(stdout, "NetCDF putting cell angle VID units: %s\n", nc_strerror(err));
    }

        /*
         *  set fill mode
         */
      err = nc_set_fill(NCInfo->ncid, NC_NOFILL, &i);
      if (err != NC_NOERR) fprintf(stdout, "NetCDF setting fill value: %s\n", nc_strerror(err));

        /*
         *  end of NetCDF definitions
         */
      err = nc_enddef(NCInfo->ncid);
      if (err != NC_NOERR) fprintf(stdout, "NetCDF error on ending definitions: %s\n", nc_strerror(err));
      
        /*
         *  specify spatial dimension labels
         */
      start[0] = 0;
      count[0] = 3;

      xyz[0] = 'x';
      xyz[1] = 'y';
      xyz[2] = 'z';
      err = nc_put_vara_text(NCInfo->ncid, NCInfo->spatialVID, start, count, xyz);
      if (err != NC_NOERR) 
	fprintf(stdout, "Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n", 
		nc_strerror(err));

      xyz[0] = 'a';
      xyz[1] = 'b';
      xyz[2] = 'c';
      err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_spatialVID, start, count, xyz);
      if (err != NC_NOERR) 
	fprintf(stdout, "Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n", 
		nc_strerror(err));


      start[0] = 0;
      start[1] = 0;
      count[0] = 3;
      count[1] = 5;
      err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_angularVID, start, count, abc);
      if (err != NC_NOERR) 
	fprintf(stdout, "Error on NetCDF output of angular VID 'alpha', 'beta ' and 'gamma': %s\n", 
		nc_strerror(err));

      outInfo->info = (void *) NCInfo;
      /*
       *  end of pre-processing upon first visit
       */
#endif
      break;

    case COORD_AMBER_TRAJECTORY:

       /*
        *  preprocessing to open file, if necessary
        */
      if (append) {

	if ( ! openFile(&outInfo->file, outInfo->filename, "a") ) {
	  fprintf(stdout, "WARNING in ptrajOutputCoordinates(): Error opening\n");
	  fprintf(stdout, "output coordinate file in append mode (%s)\n", outInfo->filename);
	  safe_fclose(outInfo->file);
	  outInfo->file = NULL;
	}

      } else {

	if ( ! openFile(&outInfo->file, outInfo->filename, "w") ) {
	  fprintf(stdout, "WARNING in ptrajOutputCoordinates(): Error opening\n");
	  fprintf(stdout, "output coordinate file (%s)\n", outInfo->filename);
	  safe_fclose(outInfo->file);
	  outInfo->file = NULL;
	}

	if (outInfo->title == NULL) 
	  fprintf(outInfo->file, "ptraj generated trajectory\n");
	else
	  fprintf(outInfo->file, "%s\n", outInfo->title);
      }
    
      break;

    case COORD_CHARMM_TRAJECTORY:

       /*
        *  preprocessing to open file and write header
        */
      if (append)
	error(ROUTINE, "Appending CHARMM trajectories is not implemented\n");

      if ( ! openFile(&outInfo->file, outInfo->filename, "w") ) {
	fprintf(stdout, "WARNING in ptrajOutputCoordinates(): Error opening\n");
	fprintf(stdout, "output coordinate file (%s)\n", outInfo->filename);
      }
      dumpCharmmTrajectory(outInfo->file, (charmmTrajectoryInfo *) outInfo->info,
			   atoms, NULL, NULL, NULL, NULL, -1);
      break;

    case COORD_AMBER_RESTART:

      if (append) {
	sprintf(buffer, "%s", outInfo->filename);
	status = openFile(&outInfo->file, buffer, "a");
      } else {
	sprintf(buffer, "%s.%i", outInfo->filename, set);
	status = openFile(&outInfo->file, buffer, "w");
      }

      if (status == 0) {
	fprintf(stdout, "WARNING in ptrajOutputCoordinates(): Could not open\n");
	fprintf(stdout, "output coordinate file %s for %s.  Not dumping to output file.\n",
		buffer, (append ? "appending" : "writing")); 
	safe_fclose(outInfo->file);
	outInfo->file = NULL;
      }
  
      break;
    
    case COORD_BINPOS:

       /*
        *  preprocessing; open input file and write magic header
        */

      if ( ! openFile(&outInfo->file, outInfo->filename, (append ? "a" : "w")) ) {
	fprintf(stdout, "WARNING in ptrajOutputCoordinates(): Error opening\n");
	fprintf(stdout, "output coordinate file (%s) for %s\n", outInfo->filename,
		(append ? "appending" : "writing"));
      }
      fwrite( "fxyz", 4, 1, outInfo->file );

      break;
    
    default:

      break;
    }
    /*
     *  end of pre-processing on first call
     */
  }
  

  switch( outInfo->type ) {
	    
  case COORD_PDB:

    if (append) {
      sprintf(buffer, "%s", outInfo->filename);
      status = openFile(&outInfo->file, buffer, "a");
    } else {
      sprintf(buffer, "%s.%i", outInfo->filename, set);
      status = openFile(&outInfo->file, buffer, "w");
    }
    if (status == 0) {
      fprintf(stdout, "ptrajOutputCoordinates(): Could not %s PDB file %s\n", 
	      (append ? "append" : "write" ), buffer);
    } else {

      if (outInfo->title != NULL)
	fprintf(outInfo->file, "REMARK %3d %47s (set %5d)\n", 1, outInfo->title, set);
      /*
	fprintf(outInfo->file, "TITLE     %s (set %i)\n", outInfo->title, set);
      */
      if( outInfo->les_action != LES_NONE )
      {
          int* mask = lesMask( state->atoms );
          pdb = ptrajStateToPdb(state, mask, outInfo->option2);
          free( mask );
      }
      else
      {   
          pdb = ptrajStateToPdb(state, NULL, outInfo->option2);
      }

      putCoordinatesInPdb(pdb, atoms, X, Y, Z);

      if (outInfo->option1 > 0) {
	   /*
            *  include charges and radii and dump out in higher precision
            */
	putQandRInPdb(pdb, outInfo->option1, outInfo->option2, 0);
	savePdbHigh(outInfo->file, pdb);
      } else {
	savePdb(outInfo->file, pdb);
      }
      if (append) fprintf(outInfo->file, "TER\n");
      safe_fclose(outInfo->file);
      outInfo->file = NULL;
      safe_free( (void *) pdb );

    }
    break;

  case COORD_AMBER_NETCDF:

#ifdef BINTRAJ
    NCInfo = (netcdfTrajectoryInfo *) outInfo->info;
    if (NCInfo == NULL) return;

    if (NCInfo->R == NULL)
      NCInfo->R = safe_malloc(sizeof(double) * atoms * 3);

    for (i = 0, j = 0; i < atoms; i++, j += 3) {
      NCInfo->R[ j   ] = X[i];
      NCInfo->R[ j+1 ] = Y[i];
      NCInfo->R[ j+2 ] = Z[i];
    }

    start[0] = NCInfo->currentFrame;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = atoms;
    count[2] = 3;
    err = nc_put_vara_float(NCInfo->ncid, NCInfo->coordinateVID, start, count, NCInfo->R);
    if (err != NC_NOERR)
      fprintf(stdout, "NetCDF error on output of coordinates in set %i: %s\n", set, nc_strerror(err));

    if (outInfo->isBox) {

      start[0] = NCInfo->currentFrame;
      start[1] = 0;
      start[2] = 0;
      count[0] = 1;
      count[1] = 3;
      count[2] = 0;

      if (prnlev > 4) {
	fprintf(stdout, "\nPUT BOX INFO: %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		box[0], box[1], box[2], box[3], box[4], box[5]);
      }

      err = nc_put_vara_double(NCInfo->ncid, NCInfo->cellLengthVID, start, count, box);
      if (err != NC_NOERR)
	fprintf(stdout, "NetCDF error on output of cell lengths in set %i: %s\n", set, nc_strerror(err));

      err = nc_put_vara_double(NCInfo->ncid, NCInfo->cellAngleVID, start, count, &box[3]);
      if (err != NC_NOERR)
	fprintf(stdout, "NetCDF error on output of cell angles in set %i: %s\n", set, nc_strerror(err));

    }

    start[0] = NCInfo->currentFrame;
    count[0] = 1;
    /*
     *  TODO: Figure out time
     */
    time = 0.0;
    err = nc_put_vara_float(NCInfo->ncid, NCInfo->timeVID, start, count, &time);

    if (err != NC_NOERR)
	fprintf(stdout, "NetCDF error on output of time in set %i: %s\n", set, nc_strerror(err));
    NCInfo->currentFrame++;
    nc_sync(NCInfo->ncid);
#endif
    break;


  case COORD_AMBER_TRAJECTORY:

    if (outInfo->isBox) {
      dumpAmberTrajectory(outInfo, outInfo->file, atoms, X, Y, Z, box);
    } else {
      dumpAmberTrajectory(outInfo, outInfo->file, atoms, X, Y, Z, NULL);
    }
	    
    break;

  case COORD_CHARMM_TRAJECTORY:

    dumpCharmmTrajectory(outInfo->file, (charmmTrajectoryInfo *) outInfo->info,
			 atoms, X, Y, Z, box, set);
    break;

  case COORD_AMBER_RESTART:

    if (append) {
      sprintf(buffer, "%s", outInfo->filename);
      status = openFile(&outInfo->file, buffer, "a");
    } else {
      sprintf(buffer, "%s.%i", outInfo->filename, set); 
      status = openFile(&outInfo->file, buffer, "w");
    }
    
    if (status == 0) {
      fprintf(stdout, "ptrajOutputCoordinates(): Could not %s coordinate file %s, not dumping...\n",
	      (append ? "append" : "write"), outInfo->filename);
      break;
    }

    dumpAmberRestart(outInfo->file, atoms, X, Y, Z, NULL, NULL, NULL, 
		     (outInfo->isBox ? box : NULL), outInfo->title);
    safe_fclose(outInfo->file);
    outInfo->file = NULL;
    break;
    
  case COORD_BINPOS:

       /*
        *  preprocessing; open input file and write magic header
        */

    writebinpos(outInfo->file, atoms, X, Y, Z);

    break;
    
  default:
	      
    return;
  }
}



#undef  ROUTINE
#define ROUTINE "parseHBondDonor()"

   void
parseHBondDonor(stackType **argumentStackPointer, ptrajState *state, int *hbondDonor)
{
  char *buffer, *buffer1;
  int *mask, i, j;
  Name atom, res;

  buffer = argumentStackKeyToString( argumentStackPointer, "mask", NULL );
     /*
      *  mask
      */
  if (buffer != NULL) {
    mask = processAtomMask(buffer, state);
    atomMaskIsActive(mask, state, &i, &j);
    if (i==0) {
      fprintf(stdout, "WARNING in ptraj, donor: No atoms selected (%s), ignoring...\n",
	      buffer);
      safe_free(buffer);
      safe_free(mask);
      return;
    }

    for (i=0; i < state->atoms; i++) {
      if (mask[i] != 0) {
	hbondDonor[i] = 1;
	if (prnlev > 0)
	  fprintf(stdout, "  DONOR: adding atom %5i, residue %4i, atom name %s\n", 
		  i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
      }
    }
    safe_free(buffer);
    safe_free(mask);
    return;
  }
    
     /*
      *  resname atomname
      */
  buffer  = getArgumentString( argumentStackPointer, NULL );
  buffer1 = getArgumentString( argumentStackPointer, NULL );
  if (buffer != NULL && buffer1 != NULL) {

       /*
        *  copy in atom and residue names and pad with spaces as appropriate
        */
    for (i=0; i < 4; i++) {
      res[i] = (char) 0;
      atom[i] = (char) 0;
    }
    strncpy(res,  buffer,  5);
    strncpy(atom, buffer1, 5);
    res[4]  = (char) 0;
    atom[4] = (char) 0;
    i = 3;
    while (res[i] == (char) 0) {
      res[i] = ' ';
      i--;
    }
    i = 3;
    while (atom[i] == (char) 0) {
      atom[i] = ' ';
      i--;
    }
    
    for (i=0; i < state->atoms; i++) {
      j = atomToResidue(i+1, state->residues, state->ipres)-1;
      if (strcmp(atom, state->atomName[i]) == 0 &&
	  strcmp(res, state->residueName[j]) == 0) {
	hbondDonor[i] = 1;
	if (prnlev > 0)
	  fprintf(stdout, "  DONOR: adding atom %5i, residue %4i, atom name %s\n", 
		  i+1, j+1, state->atomName[i]);
      }
    }
    safe_free(buffer);
    safe_free(buffer1);
  } else {
    fprintf(stdout, 
	    "WARNING in ptraj, donor: either the residue or atom name specification is blank!\n");
  }
  return;
}


#undef  ROUTINE
#define ROUTINE "parseHBondAcceptor()"

   void
parseHBondAcceptor(stackType **argumentStackPointer, ptrajState *state, int *hbondAcceptor,
		   int *hbondAcceptorH1, int *hbondAcceptorH2, int *hbondAcceptorH3)
{
  char *buffer, *buffer1;
  int *mask, *mask1, i, j, i1, j1, k, stop;
  Name atom, atom1, res;

  buffer1 = NULL;
  buffer  = argumentStackKeyToString( argumentStackPointer, "mask", NULL );
  if (buffer != NULL) {
    buffer1 = getArgumentString( argumentStackPointer, NULL );
  }
     /*
      *  mask
      */
  if (buffer != NULL && buffer1 != NULL) {
    mask = processAtomMask(buffer, state);
    atomMaskIsActive(mask, state, &i, &i1);
    mask1 = processAtomMask(buffer1, state);
    atomMaskIsActive(mask1, state, &j, &j1);

    stop = 0;
    if (i==0) {
      fprintf(stdout, 
	      "WARNING in ptraj, acceptor: No heavy atom was selected (%s), ignoring...\n",
	      buffer);
      stop = 1;
    }
    if (j==0) {
      fprintf(stdout,
	      "WARNING in ptraj, acceptor: No hydrogen atom was selected (%s), ignoring...\n",
	      buffer1);
      stop = 1;
    }

    if (i != j) {
      fprintf(stdout,
	      "WARNING in ptraj, acceptor: There is not a 1-1 correspondence between the\n");
      fprintf(stdout, "atom selection in the two masks %s and %s which contain %i and %i\n",
	      buffer, buffer1, i, j);
      fprintf(stdout, "atoms respectively.  Ignoring...\n");
      stop = 1;
    }

    if (stop) {
      safe_free(buffer);
      safe_free(buffer1);
      safe_free(mask);
      safe_free(mask1);
      return;
    }


    stop = i;
    for (i=0; i < stop; i++) {
    
      hbondAcceptor[i1] = 1;
      if (hbondAcceptorH1[i1] >= 0) {
	if (hbondAcceptorH2[i1] >= 0) {
	  if (hbondAcceptorH3[i1] >= 0) {
	    fprintf(stdout,
		    "WARNING in ptraj, acceptor: More than three hydrogens have been selected\n");
	    fprintf(stdout, "as acceptors for the atom %5i.  Ignoring the latest...\n",
		    i1+1);
	  } else
	    hbondAcceptorH3[i1] = j1;
	} else
	  hbondAcceptorH2[i1] = j1;
      } else
	hbondAcceptorH1[i1] = j1;
    
      if (prnlev > 0) {
	fprintf(stdout, "  ACCEPTOR: adding atom %5i, res %4i, atom name %s -- ",
		i1+1, atomToResidue(i1+1, state->residues, state->ipres), state->atomName[i1]);
	fprintf(stdout, "atom %5i, res %4i, atom name %s\n",
		j1+1, atomToResidue(j1+1, state->residues, state->ipres), state->atomName[j1]);
      }
      mask[i1] = 0;
      mask1[j1] = 0;
      while (mask[i1] == 0 && i1 < state->atoms) i1++;
      while (mask1[j1] == 0 && j1 < state->atoms) j1++;
      if (i1 == state->atoms || j1 == state->atoms) return;
    }
    safe_free(mask);
    safe_free(mask1);
    safe_free(buffer);
    safe_free(buffer1);
    return;
  } else {
    if (buffer1 != NULL || buffer != NULL) {
      
      fprintf(stdout, 
	      "WARNING in ptraj, acceptor: Error in mask specification.  Ignoring...\n");
      safe_free(buffer);
      safe_free(buffer1);
      return;
    }
  }
    
     /*
      *  resname atomname atomname
      */
  stop = 0;
  buffer  = getArgumentString( argumentStackPointer, NULL );
  if (stop || buffer == NULL)
    stop = 1;
  else {
    for (i=0; i < 4; i++) {
      res[i] = (char) 0;
    }
    strncpy(res,  buffer,  5);
    res[4]  = (char) 0;
    i = 3;
    while (res[i] == (char) 0) {
      res[i] = ' ';
      i--;
    }
    safe_free(buffer);
  }
  
  buffer  = getArgumentString( argumentStackPointer, NULL );
  if (stop || buffer == NULL)
    stop = 1;
  else {
    for (i=0; i < 4; i++) {
      atom[i] = (char) 0;
    }
    strncpy(atom,  buffer,  5);
    atom[4]  = (char) 0;
    i = 3;
    while (atom[i] == (char) 0) {
      atom[i] = ' ';
      i--;
    }
    safe_free(buffer);
  }

  buffer  = getArgumentString( argumentStackPointer, NULL );
  if (stop || buffer == NULL)
    stop = 1;
  else {
    for (i=0; i < 4; i++) {
      atom1[i] = (char) 0;
    }
    strncpy(atom1,  buffer,  5);
    atom1[4]  = (char) 0;
    i = 3;
    while (atom1[i] == (char) 0) {
      atom1[i] = ' ';
      i--;
    }
    safe_free(buffer);
  }

  if (stop == 0) {
    
    for (i=0; i < state->atoms; i++) {
      j = atomToResidue(i+1, state->residues, state->ipres)-1;
      if (strcmp(atom, state->atomName[i]) == 0 &&
	  strcmp(res, state->residueName[j]) == 0) {
	
	for (k=state->ipres[j]-1; k < state->ipres[j+1]-1; k++) {
	  if (strcmp(atom1, state->atomName[k]) == 0) {
	    hbondAcceptor[i] = 1;
	    if (hbondAcceptorH1[i] >= 0)
	      if (hbondAcceptorH2[i] >= 0)
		if (hbondAcceptorH3[i] >= 0) {
		  fprintf(stdout, "WARNING in ptraj, acceptor: More than three hydrogens ");
		  fprintf(stdout, "have been selected\nas acceptors for atom %i. ", i+1);
		  fprintf(stdout, "Ignoring the latest...\n");
		} else
		  hbondAcceptorH3[i] = k;
	      else
		hbondAcceptorH2[i] = k;
	    else
	      hbondAcceptorH1[i] = k;

	    if (prnlev > 0) {
	      fprintf(stdout, "  ACCEPTOR: adding atom %5i, res %4i, atom name %s -- ",
		      i+1, atomToResidue(i+1, state->residues, state->ipres), 
		      state->atomName[i]);
	      fprintf(stdout, "atom %5i, res %4i, atom name %s\n",
		      k+1, atomToResidue(k+1, state->residues, state->ipres), 
		      state->atomName[k]);
	    }
	  }
	}
      }
    }
  } else 
    fprintf(stdout, 
	    "WARNING in ptraj, acceptor: error in specification of res/atom selection\n");
}




/*
 *  ptrajSetupIO(): This is called upon receipt of the trigger strings
 *  that specify input and output files, reference structures and
 *  global variables that may be used by various actions (such as solvent
 *  information or hydrogen bonding donor/acceptor atoms).
 *  For I/O, the transformFileStack list of input files is set up as well
 *  as the output trajectory file.
 */

#undef  ROUTINE
#define ROUTINE "ptrajSetupIO()"

   void
ptrajSetupIO(stackType *argumentStack, actionType type)
{
  char *filename;
  char *buffer;
  char *title, *titlenew;
  int i, j, k, start, stop, offset;
  coordinateInfo *info;
  ptrajState *state, **statep;
  int *mask, byres, bytype, curres, *startsol, *stopsol;
  Name res;
  charmmTrajectoryInfo *charmmTraj, *ctrj;
  stackType *sp;
  coordinateInfo *infiles;
  double boxtmp;

  statep = ptrajCurrentState();
  state = *statep;

  /*
   *  Input/Output setup for ptraj.  The following command/arguments are parsed:
   *
   *  TRANSFORM_BOX
   *
   *    box [x value] [y value] [z value] [alpha value] [beta value] [gamma value] 
   *        [fixx] [fixy] [fixz] [fixalpha] [fixbeta] [fixgamma]
   *
   *  TRANSFORM_PRNLEV
   *
   *    prnlev <value>
   *
   *  TRANSFORM_REFERENCE
   *
   *    reference filename
   *
   *  TRANSFORM_SOLVENT
   *
   *    solvent mask1 [mask2] ... [maskN] [byres | bytype | byname]
   *
   *  TRANSFORM_TRAJIN
   *
   *    trajin filename [start] [stop] [delta]
   *
   *  TRANSFORM_TRAJOUT
   *
   *    trajout filename [nobox] [ PDB | RESTART | BINPOS | NETCDF ] [append]
   *      [title <title>] [application <application>] [program <program>]
   *
   *  TRANSFORM_DONOR
   *
   *    donor {print | clear | <resname> <atomname> | mask <mask>}
   *
   *  TRANSFORM_ACCEPTOR
   *
   *    acceptor print | clear | 
   *       { <resname> <atomname1> <atomname2> } | 
   *       { mask <mask1> <mask2> }
   *
   */


  switch ( type ) {

  case TRANSFORM_PRNLEV:

    prnlev = getArgumentInteger(&argumentStack, 0.0);
    fprintf(stdout, "  PRNLEV: value is now %i\n", prnlev);
    return;

  case TRANSFORM_BOX:

    boxtmp = argumentStackKeyToDouble(&argumentStack, "x", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[0] = boxtmp;
      fprintf(stdout, "  BOX X = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "y", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[1] = boxtmp;
      fprintf(stdout, "  BOX Y = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "z", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[2] = boxtmp;
      fprintf(stdout, "  BOX Z = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "alpha", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[3] = boxtmp;
      fprintf(stdout, "  BOX ALPHA = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "beta", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[4] = boxtmp;
      fprintf(stdout, "  BOX BETA  = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "gamma", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[5] = boxtmp;
      fprintf(stdout, "  BOX GAMMA = %f\n", boxtmp);
    }

    if (argumentStackContains( &argumentStack, "fixx" ))
      state->boxfixed[0] = 1;
    if (argumentStackContains( &argumentStack, "fixy" ))
      state->boxfixed[1] = 1;
    if (argumentStackContains( &argumentStack, "fixz" ))
      state->boxfixed[2] = 1;
    if (argumentStackContains( &argumentStack, "fixalpha" ))
      state->boxfixed[3] = 1;
    if (argumentStackContains( &argumentStack, "fixbeta" ))
      state->boxfixed[4] = 1;
    if (argumentStackContains( &argumentStack, "fixgamma" ))
      state->boxfixed[5] = 1;

    return;

  case TRANSFORM_DONOR:

       /*
        *  print
        */
    if (argumentStackContains( &argumentStack, "print" )) {
      if (hbondDonor != NULL) {
	printHBondMask(hbondDonor, NULL, NULL, NULL, state);
      }
      return;
    }

       /*
        *  clear
        */
    if (argumentStackContains( &argumentStack, "clear" )) {
      if (hbondDonor != NULL) {
	pushBottomStack(&hbondDonorStack, (void *) hbondDonor);
	hbondDonor = NULL;
      }
      return;
    }

       /*
        *  allocate space for the hbondDonor if necessary...
        */
    if (hbondDonor == NULL) {
      hbondDonor = (int *) safe_malloc(sizeof(int) * state->atoms);
      for (k=0; k < state->atoms; k++)
	hbondDonor[k] = 0;
    }

    parseHBondDonor(&argumentStack, state, hbondDonor);

       /*
        *  make sure at least one donor has been chosen or free the memory
        */
    atomMaskIsActive(hbondDonor, state, &i, &j);
    if (i == 0) {
      safe_free(hbondDonor);
      hbondDonor = NULL;
    }
    return;


  case TRANSFORM_ACCEPTOR:

       /*
        *  print
        */
    if (argumentStackContains( &argumentStack, "print" )) {
      if (hbondAcceptor != NULL) {
	printHBondMask(hbondAcceptor, hbondAcceptorH1, hbondAcceptorH2, hbondAcceptorH3, state);
      }
      return;
    }

       /*
        *  clear
        */
    if (argumentStackContains( &argumentStack, "clear" )) {
      if (hbondAcceptor != NULL) {
	pushBottomStack(&hbondAcceptorStack, (void *) hbondAcceptor);
	pushBottomStack(&hbondAcceptorH1Stack, (void *) hbondAcceptorH1);
	pushBottomStack(&hbondAcceptorH2Stack, (void *) hbondAcceptorH2);
	pushBottomStack(&hbondAcceptorH3Stack, (void *) hbondAcceptorH3);
	hbondAcceptor = NULL;
	hbondAcceptorH1 = NULL;
	hbondAcceptorH2 = NULL;
	hbondAcceptorH3 = NULL;
      }
      return;
    }

       /*
        *  allocate space for the hbondAcceptor lists if necessary...
        */
    if (hbondAcceptor == NULL) {
      hbondAcceptor = (int *) safe_malloc(sizeof(int) * state->atoms);
      hbondAcceptorH1 = (int *) safe_malloc(sizeof(int) * state->atoms);
      hbondAcceptorH2 = (int *) safe_malloc(sizeof(int) * state->atoms);
      hbondAcceptorH3 = (int *) safe_malloc(sizeof(int) * state->atoms);
      for (k=0; k < state->atoms; k++) {
	hbondAcceptor[k] = 0;
	hbondAcceptorH1[k] = -1;
	hbondAcceptorH2[k] = -1;
	hbondAcceptorH3[k] = -1;
      }
    }

    parseHBondAcceptor(&argumentStack, state, hbondAcceptor, 
		       hbondAcceptorH1, hbondAcceptorH2, hbondAcceptorH3);

       /*
        *  make sure at least one acceptor has been chosen or free the memory
        */
    atomMaskIsActive(hbondAcceptor, state, &i, &j);
    if (i == 0) {
      safe_free(hbondAcceptor);
      safe_free(hbondAcceptorH1);
      safe_free(hbondAcceptorH2);
      safe_free(hbondAcceptorH3);
      hbondAcceptor = NULL;
      hbondAcceptorH1 = NULL;
      hbondAcceptorH2 = NULL;
      hbondAcceptorH3 = NULL;
    }
    return;

  case TRANSFORM_REFERENCE:

    filename = getArgumentString( &argumentStack, NULL );
    if (filename == NULL) {
      fprintf(stdout, "WARNING in ptrajSetupIO(): reference command lacks a filename!\n");
      return;
    }

    info = checkCoordinates(filename, state->atoms);
    safe_free(filename);

    if ( info == NULL )
      return;

    info->mask = NULL;
    info->offset = 1;
    info->start = 1;
    info->stop = -1;
    info->option1 = 0;
    info->option2 = 0;

    if ( ptrajPreprocessInputCoordinates(info) ) {
      safe_free(info->filename);
      safe_free(info->info);
      return;
    }

    info->x = (double *) safe_malloc(sizeof(double) * state->atoms);
    info->y = (double *) safe_malloc(sizeof(double) * state->atoms);
    info->z = (double *) safe_malloc(sizeof(double) * state->atoms);

    /*
     *  note, we initialize the memory here since the checks to disable processing
     *  of truncated sets have been removed, i.e. it is possible to read in a truncated
     *  PDB reference structure noting that the truncated part will be replaced by zeros.
     *  If the set is truncated, the user will be warned...
     */
    for (i=0; i < state->atoms; i++) {
      info->x[i] = 0.0;
      info->y[i] = 0.0;
      info->z[i] = 0.0;
    }
    ptrajProcessInputCoordinates(info, state, info->x, info->y, info->z, state->box,
				 1, &i, &j);
    if (j == 0) {
      safe_free(info->x);
      safe_free(info->y);
      safe_free(info->z);
      safe_free(info->filename);
      safe_free(info->info);
      return;
    }

    referenceInfo = info;
    pushBottomStack( &transformReferenceStack, (void *) info );
    
    return;

  case TRANSFORM_SOLVENT:

    /*
     *  solvent [byres | bytype | byname] mask1 [mask2] [mask3] ...
     */

    byres  = argumentStackContains( &argumentStack, "byres");
    bytype = argumentStackContains( &argumentStack, "bytype");
    if ( argumentStackContains( &argumentStack, "byname") == 1) {
      byres  = 0;
      bytype = 0;
    }

    if (bytype == 1) {
      fprintf(stderr, "WARNING in ptrajSetupIO(): solvent \"bytype\" option has not yet\n");
      fprintf(stderr, "been implemented.  Defaulting to byname...\n");
      bytype = 0;
    }

    buffer = (char *) argumentStack->entry;

       /*
        *  if we've exhausted the argument stack, assume this is a call to
	*  solvent with no arguments; in this case we do not alter the solvent
	*  information, we just print a summary of it...
        */
    if (buffer[0] == (char) 0) {
      ptrajPrintState(state);
      return;
    }

    /*
     *  re-initialize the solvent information
     */
    state->solventMolecules = 0;
    state->solventAtoms = 0;
    if (state->solventMask != NULL) safe_free(state->solventMask);
    state->solventMask = NULL;
    if (state->solventMoleculeStart != NULL) safe_free(state->solventMoleculeStart);
    state->solventMoleculeStart = NULL;
    if (state->solventMoleculeStop != NULL) safe_free(state->solventMoleculeStop);
    state->solventMoleculeStop = NULL;

    state->solventMask = (int *) safe_malloc(sizeof(int) * state->atoms);
    for (i=0; i < state->atoms; i++)
      state->solventMask[i] = 0;
    startsol = (int *) safe_malloc(sizeof(int) * state->atoms);
    stopsol  = (int *) safe_malloc(sizeof(int) * state->atoms);

    buffer = NULL;
    while ( (buffer = getArgumentString(&argumentStack, NULL)) != NULL ) {

      if (byres) {
	/*
	 *  load up solvent by residue using the residues specified in
	 *  the input masks...
	 */

	mask = processAtomMask(buffer, state);
	fprintf(stdout, "       Searching for solvent by mask ");
	printAtomMask(mask, state);
	
	for (i=0; i < state->atoms; i++) {

	  if (mask[i]) {
	    curres = atomToResidue(i+1, state->residues, state->ipres)-1;
	    j = isActiveResidue(curres, mask, state);

	    if (j > 0) {
	      /*
	       *  all the atoms in "curres" are active, hence add this
	       *  solvent molecule to the list of solvent.  Note that
	       *  isActiveResidue returns the number of active atoms
	       *  found (if not zero)...
	       */

	      state->solventAtoms += j;
	      startsol[state->solventMolecules] = i;
	      stopsol[ state->solventMolecules] = i+j;
	      state->solventMolecules++;
	      for (k=i; k < i+j; k++)
		state->solventMask[k] = 1;
	      i += j-1;
	    }
	  }
	}
    
      } else if (bytype) {

	/*
	 *  search for solvent using the mask to define
	 *  what a representative solvent molecule is...
	 *
	 *  CURRENTLY NOT IMPLEMENTED
	 */
	mask = processAtomMask(buffer, state);

      } else {
	/*
	 *  search for solvent by residue name
	 */

           /*
            *  copy the residue name from the buffer into "res" and pad
            *  with spaces to conform to standard atom/residue naming
            */
	for (i=0; i < 4; i++)
	  res[i] = (char) 0;
	strncpy(res, buffer, 5);
	res[4] = (char) 0;
	i = 3;
	while (res[i] == (char) 0) {
	  res[i] = ' ';
	  i--;
	}

	fprintf(stdout, "       Searching for solvent by residue name %s\n", res);

           /*
            *  loop over all residues and check for matches
            */
	for (i=0; i < state->residues; i++) {
	  if (strcmp(res, state->residueName[i]) == 0) {
	    /*
	     * add this residue to the list of solvent
	     */
	    j = state->ipres[i+1]-state->ipres[i];

	    state->solventAtoms += j;
	    startsol[state->solventMolecules] = state->ipres[i]-1;
	    stopsol[ state->solventMolecules] = state->ipres[i+1]-1;
	    state->solventMolecules++;
	    for (k=state->ipres[i]-1; k < state->ipres[i+1]-1; k++)
	      state->solventMask[k] = 1;
	  }
	}

      }

      safe_free(mask);
      safe_free(buffer);
    }

       /*
        *  update the solvent information and free any unnecessary memory
        */
    state->solventMoleculeStart = (int *) 
      safe_malloc(sizeof(int) * state->solventMolecules);
    state->solventMoleculeStop = (int *) 
      safe_malloc(sizeof(int) * state->solventMolecules);
    for (i=0; i < state->solventMolecules; i++) {
      state->solventMoleculeStart[i] = startsol[i];
      state->solventMoleculeStop[i]  = stopsol[i];
    }

    safe_free(buffer);
    safe_free(stopsol);
    safe_free(startsol);
    if (prnlev > 0)
      ptrajPrintState(state);
    return;

  case TRANSFORM_TRAJIN:

    filename = getArgumentString( &argumentStack, NULL );
    if (filename == NULL) {
      fprintf(stdout, "WARNING in ptrajSetupIO: trajin command lacks a filename!\n");
      return;
    }

    /*
     *  check the coordinates and push the information to the bottom
     *  of the transformFileStack
     */
    fprintf(stderr, "  Checking coordinates: %s\n", filename);
    info = checkCoordinates(filename, state->atoms);

    if (info == (coordinateInfo *) NULL) {
      safe_free(filename);
      return;
    }
    safe_free(filename);

    /*
     *  set start, stop and offset if they were specified and relevent
     */

    start  = getArgumentInteger( &argumentStack,  1 );
    stop   = getArgumentInteger( &argumentStack, -1 );
    offset = getArgumentInteger( &argumentStack,  1 );

    if (stop > 0) {
      if (stop > info->stop) {
	fprintf(stdout, 
		"FYI %s: trajin stop value (%i) is greater than the number of sets read in.\n", 
		ROUTINE, stop);
	fprintf(stdout, "Setting stop to the maximum value (%i)\n",
		info->stop);
      } else {
	info->stop = stop;
      }
    }

    if (start > 1) {
      if (info->stop > 0 && start > info->stop) {
	fprintf(stdout, "WARNING in %s: trajectory start is > stop; no\n", ROUTINE);
	fprintf(stdout, "configurations will be processed\n");
      }
      info->start = start;
    }

    info->offset = offset;
    if (info->offset != 1 && info->offset > info->stop - info->start) {
      fprintf(stdout, "WARNING in %s: Offset is so large that only 1 set\n", ROUTINE);
      fprintf(stdout, "  will be processed...\n");
    }

    pushBottomStack( &transformFileStack, (void *) info );
    return;

  case TRANSFORM_TRAJOUT:

    /*
     *  Set up the coordinateInfo structure for the file to be
     *  output.  Note that we do not actually open the filename
     *  yet as this is done by ptrajOutputCoordinates()
     */

    filename = getArgumentString( &argumentStack, NULL );
    if (filename == NULL) {
      fprintf(stdout, "WARNING in %s: trajout command lacks a filename!\n", ROUTINE);
      safe_free(filename);
      return;
    }

    info = (coordinateInfo *) safe_malloc(sizeof(coordinateInfo));
    INITIALIZE_coordinateInfo(info);
    info->filename = copyString(filename);
    if (argumentStackContains( &argumentStack, "nobox" ))
      info->isBox = 0;
    else
      info->isBox = state->IFBOX;
    info->append =      argumentStackContains( &argumentStack, "append" );

    /*
     *  TODO: Integrate LES options into processing aka strip so that actions can act on
     *  LES average or subset trajectories
     */


    /* 
     *  check if there are les options
     */

    info->les_action = LES_NONE;
    if (argumentStackContains( &argumentStack, "les" ) )
    {
        info->nlescopy = ( parm->nlestyp == 1 ) ? (int)( parm->lesfac[0] + 0.1 ) : (int)( parm->lesfac[3] + 0.1 );

        if( argumentStackContains( &argumentStack, "split" ) )
        {
            info->les_action = LES_SPLIT;
	    info->les_status = LES_READY;
        }
        else if( argumentStackContains( &argumentStack, "average" ) )
        {
            info->les_action = LES_AVERAGE;
	    info->les_status = LES_READY;
        }
        else
        {
            error( "setup_les_output", "unknown les action" );
        } 

    }
    
       /*
        *  check to see if a format other than amber trajectory is wanted
        */
    if (argumentStackContains( &argumentStack, "pdb" ))
      info->type = COORD_PDB;
    else if (argumentStackContains( &argumentStack, "restart" ))
      info->type = COORD_AMBER_RESTART;
    else if (argumentStackContains( &argumentStack, "restrt" ))
      info->type = COORD_AMBER_RESTART;
    else if (argumentStackContains( &argumentStack, "rest" ))
      info->type = COORD_AMBER_RESTART;
    else if (argumentStackContains( &argumentStack, "binpos" ))
      info->type = COORD_BINPOS;
    else if (argumentStackContains( &argumentStack, "charmm" ))
      info->type = COORD_CHARMM_TRAJECTORY;
#ifdef BINTRAJ
    else if (argumentStackContains( &argumentStack, "netcdf" ))
      info->type = COORD_AMBER_NETCDF;
    else if (argumentStackContains( &argumentStack, "cdf" ))
      info->type = COORD_AMBER_NETCDF;
#else
    else if (argumentStackContains( &argumentStack, "netcdf" )) {
      info->type = COORD_AMBER_TRAJECTORY;
      fprintf(stdout, "trajout: NetCDF support is not compiled into this version (Add -DBINTRAJ)\n");
      fprintf(stdout, "         defaulting to AMBER trajectory format...\n");
    } else if (argumentStackContains( &argumentStack, "cdf" )) {
      info->type = COORD_AMBER_TRAJECTORY;
      fprintf(stdout, "trajout: NetCDF support is not compiled into this version (Add -DBINTRAJ)\n");
      fprintf(stdout, "         defaulting to AMBER trajectory format...\n");
    }
#endif
    else
      info->type = COORD_AMBER_TRAJECTORY;

      /*
       *  set up program/version information
       */
    switch (info->type) {
    case COORD_PDB:
      info->title = argumentStackKeyToString( &argumentStack, "title", "PDB file generated by ptraj" );
      break;
    case COORD_AMBER_NETCDF:
      /*
      info->title = argumentStackKeyToString( &argumentStack, "title", "NetCDF trajectory generated by ptraj" );
      */
      info->title = argumentStackKeyToString( &argumentStack, "title", "" );
      break;

    default:
      info->title = argumentStackKeyToString( &argumentStack, "title", "trajectory generated by ptraj" );
    }
    info->application = argumentStackKeyToString( &argumentStack, "application", "AMBER" );
    /*
    info->program =     argumentStackKeyToString( &argumentStack, "program", "ptraj" );
    */
    info->program =     argumentStackKeyToString( &argumentStack, "program", "sander" );
    /*
    info->version =     argumentStackKeyToString( &argumentStack, "version", PTRAJ_VERSION_STRING );
    */
    info->version =     argumentStackKeyToString( &argumentStack, "version", "9.0" );

       /*
        *  check to see if we want charges/radii dumped to pdb
        */
    if (info->type == COORD_PDB) {

      if (argumentStackContains( &argumentStack, "dumpq" ))
	info->option1 = 1;
      if (argumentStackContains( &argumentStack, "parse" ))
	info->option1 = 2;
      if (argumentStackContains( &argumentStack, "nowrap" ))
	info->option2 = 1;
    }

       /*
        *  check to see if other CHARMM related information is present
        */
    if (info->type == COORD_CHARMM_TRAJECTORY) {

      /*
       *  search through the list of input files to setup the CHARMM information
       *  structure defaults; if none is present, make it up and/or modify according
       *  to what the user specifies...
       */
      charmmTraj = (charmmTrajectoryInfo *) safe_malloc(sizeof(charmmTrajectoryInfo));
      INITIALIZE_charmmTrajectoryInfo(charmmTraj);

      ctrj = NULL;
      for (sp = transformFileStack; ctrj == NULL && sp != NULL; sp = sp->next) {
	infiles = (coordinateInfo *) sp->entry;
	if (infiles->type == COORD_CHARMM_TRAJECTORY) 
	  ctrj = (charmmTrajectoryInfo *) infiles->info;
      }

      if (ctrj != NULL) {
	charmmTraj->byteorder = ctrj->byteorder;
	charmmTraj->magic = ctrj->magic;
	for (i=0; i < 20; i++)
	  charmmTraj->icntrl[i] = ctrj->icntrl[i];
	charmmTraj->ntitle = ctrj->ntitle;
	for (sp = ctrj->titleStack; sp != NULL; sp = sp->next) {
	  title = (char *) sp->entry;
	  titlenew = (char *) safe_malloc(sizeof(char) * (strlen(title) + 1));
	  strcpy(titlenew, title);
	  pushBottomStack(&charmmTraj->titleStack, (void *) title);
	}
	charmmTraj->natrec = ctrj->natrec;
	charmmTraj->nfreat = ctrj->nfreat;
	if (ctrj->nfreat != ctrj->natrec) {
	  charmmTraj->freeat = (int *) safe_malloc(sizeof(int) * ctrj->nfreat);
	  for (i=0; i < ctrj->nfreat; i++) 
	    charmmTraj->freeat[i] = ctrj->freeat[i];
	}
	for (i=0; i<6; i++)
	  charmmTraj->xtlabc[i] = ctrj->xtlabc[i];

      } else {
	charmmTraj->byteorder = 0;
	charmmTraj->magic.c[0] = 'C';
	charmmTraj->magic.c[1] = 'O';
	charmmTraj->magic.c[2] = 'R';
	charmmTraj->magic.c[3] = 'D';
	for (i=0; i<20; i++)
	  charmmTraj->icntrl[i] = 0;
	charmmTraj->icntrl[19] = 26;

	if (state->IFBOX) {
	  charmmTraj->icntrl[10] = 1;  /* QCRYS */
	}
      }

      if (info->isBox == 0) charmmTraj->icntrl[10] = 0;


      if (argumentStackContains( &argumentStack, "big" ))
	charmmTraj->byteorder = 0;
      else if (argumentStackContains( &argumentStack, "little" ))
	charmmTraj->byteorder = 1;
	
      info->info = (void *) charmmTraj;
    }

    safe_free(filename);
    globalOutInfo = info;
    return;
  }
}


/*
 *  ptrajSetup(): This routine is called for every trigger that is related
 *  to coordinate processing (i.e. not those commands that are I/O 
 *  related, such as trajin, trajout or reference and not those that
 *  are involved with postprocessing any acculated data).  This creates
 *  the transformActionStack stack of "actions" to be performed and 
 *  is called by dispatchToken() upon receipt of the appropriate trigger.
 *  Most of the actual setup of the action function is performed by the
 *  actual action function itself in the PTRAJ_SETUP mode.  See the
 *  detailed comments in actions.c for more information.
 */

#undef  ROUTINE
#define ROUTINE "ptrajSetup()"

   void
ptrajSetup(stackType *argumentStack, actionType type)
{
  actionInformation *action;
  int ierr;

     /*
      *  Allocate and initialize a actionInformation structure
      */

  action = (actionInformation *)
    safe_malloc(sizeof(actionInformation));
  INITIALIZE_actionInformation(action);

     /*
      *  Place a copy of the current state into the action
      */

  action->state = ptrajCopyState(ptrajCurrentState());

     /*
      *  Set the action type
      */

  action->type = type;

     /*
      *  Set the action function
      */

  switch ( type ) {

  case TRANSFORM_ANGLE:

    action->type = TRANSFORM_ANGLE;
    action->fxn  = (actionFunction) transformAngle;
    break;

  case TRANSFORM_ATOMICFLUCT:

    action->type = TRANSFORM_ATOMICFLUCT;
    action->fxn  = (actionFunction) transformAtomicFluct;
    break;

  case TRANSFORM_AVERAGE:

    action->type = TRANSFORM_AVERAGE;
    action->fxn  = (actionFunction) transformAverage;
    break;

  case TRANSFORM_CENTER:

    action->type = TRANSFORM_CENTER;
    action->fxn  = (actionFunction) transformCenter;
    break;

  case TRANSFORM_CHECKOVERLAP:

    action->type = TRANSFORM_CHECKOVERLAP;
    action->fxn  = (actionFunction) transformCheckOverlap;
    break;

  case TRANSFORM_CLOSESTWATERS:

    action->type = TRANSFORM_CLOSESTWATERS;
    action->fxn  = (actionFunction) transformClosestWaters;
    break;

  case TRANSFORM_CORRELATION:

    action->type = TRANSFORM_CORRELATION;
    action->fxn  = (actionFunction) transformCorr;
    break;

  case TRANSFORM_CONTACTS:

    action->type = TRANSFORM_CONTACTS;
    action->fxn  = (actionFunction) transformContacts;
    break;

  case TRANSFORM_DIHEDRAL:

    action->type = TRANSFORM_DIHEDRAL;
    action->fxn  = (actionFunction) transformDihedral; 
    break;

  case TRANSFORM_DIFFUSION:

    action->type = TRANSFORM_DIFFUSION;
    action->fxn  = (actionFunction) transformDiffusion; 
    break;

  case TRANSFORM_DIPOLE:

    action->type = TRANSFORM_DIPOLE;
    action->fxn  = (actionFunction) transformDipole;
    break;

  case TRANSFORM_DISTANCE:

    action->type = TRANSFORM_DISTANCE;
    action->fxn  = (actionFunction) transformDistance;
    break;

  case TRANSFORM_DNAIONTRACKER:

    action->type = TRANSFORM_DNAIONTRACKER;
    action->fxn  = (actionFunction) transformDNAiontracker;
    break;

  case TRANSFORM_ECHO:

    action->type = TRANSFORM_ECHO;
    action->fxn = (actionFunction) transformEcho;
    break;

  case TRANSFORM_ENERGY:

    action->type = TRANSFORM_ENERGY;
    action->fxn  = (actionFunction) transformEnergy;
    break;

  case TRANSFORM_GRID:

    action->type = TRANSFORM_GRID;
    action->fxn  = (actionFunction) transformGrid;
    break;

  case TRANSFORM_HBOND:

    action->type = TRANSFORM_HBOND;
    action->fxn  = (actionFunction) transformHBond;
    break;

  case TRANSFORM_IMAGE:

    action->type = TRANSFORM_IMAGE;
    action->fxn  = (actionFunction) transformImage;
    break;

  case TRANSFORM_MATRIX:

    action->type = TRANSFORM_MATRIX;
    action->fxn  = (actionFunction) transformMatrix;
    break;

  case TRANSFORM_PRINCIPAL:

    action->type = TRANSFORM_PRINCIPAL;
    action->fxn  = (actionFunction) transformPrincipal;
    break;

  case TRANSFORM_PROJECTION:

    action->type = TRANSFORM_PROJECTION;
    action->fxn  = (actionFunction) transformProjection;
    break;

  case TRANSFORM_PUCKER:

    action->type = TRANSFORM_PUCKER;
    action->fxn  = (actionFunction) transformPucker;
    break;

  case TRANSFORM_RADIAL:

    action->type = TRANSFORM_RADIAL;
    action->fxn  = (actionFunction) transformRadial;
    break;

  case TRANSFORM_RADIUSOFGYRATION:

    action->type = TRANSFORM_RADIUSOFGYRATION;
    action->fxn  = (actionFunction) transformRadiusOfGyration;
    break;

  case TRANSFORM_RANDOMIZEIONS:

    action->type = TRANSFORM_RANDOMIZEIONS;
    action->fxn  = (actionFunction) transformRandomizeIons;
    break;

  case TRANSFORM_RMS:
      
    action->type = TRANSFORM_RMS;
    action->fxn  = (actionFunction) transformRMS;
    break;

  case TRANSFORM_RUNNINGAVERAGE:
      
    action->type = TRANSFORM_RUNNINGAVERAGE;
    action->fxn  = (actionFunction) transformRunningAverage;
    break;

  case TRANSFORM_SCALE:
      
    action->type = TRANSFORM_SCALE;
    action->fxn  = (actionFunction) transformScale;
    break;

  case TRANSFORM_SECONDARYSTRUCT:

    action->type = TRANSFORM_SECONDARYSTRUCT;
    action->fxn  = (actionFunction) transformSecondaryStruct;
    break;

  case TRANSFORM_STRIP:
      
    action->type = TRANSFORM_STRIP;
    action->fxn  = (actionFunction) transformStrip;
    break;

  case TRANSFORM_TRANSLATE:

    action->type = TRANSFORM_TRANSLATE;
    action->fxn  = (actionFunction) transformTranslate;
    break;

  case TRANSFORM_TRUNCOCT:

    action->type = TRANSFORM_TRUNCOCT;
    action->fxn  = (actionFunction) transformTruncOct;
    break;

  case TRANSFORM_VECTOR:

    action->type = TRANSFORM_VECTOR;
    action->fxn  = (actionFunction) transformVector;
    break;

  case TRANSFORM_WATERSHELL:

    action->type = TRANSFORM_WATERSHELL;
    action->fxn  = (actionFunction) transformWatershell;
    break;

  case TRANSFORM_2DRMS:

    action->type = TRANSFORM_2DRMS;
    action->fxn  = (actionFunction) transform2dRMS;
    break;

  case TRANSFORM_TRANSFORM:
  case TRANSFORM_NOOP:

    action->type = type;
    action->fxn = NULL;
    break;

  default:

    fprintf(stdout, "%s: Attempting to setup an unknown action type %i\n", ROUTINE, type);
    error(ROUTINE, "There is no way you should be here!  Terminating...\n");

  }

     /*
      *  Parse the arguments.  This is done by the action->fxn in the
      *  PTRAJ_SETUP mode.  The argumentStack is placed into the
      *  complex argument 1 slot for this.
      */

  ierr = 0;
  if (action->type != TRANSFORM_TRANSFORM && 
      action->type != TRANSFORM_NOOP) {

    action->carg1 = (void *) &argumentStack;
    ierr = action->fxn(action, NULL, NULL, NULL, NULL, PTRAJ_SETUP);

  }

     /*
      *  If the setup fails, -1 is returned and therefore this action
      *  should not be placed on the action stack and the associated
      *  memory should be freed
      */

  if (ierr < 0) {
    safe_free(action->state);
    action->state= NULL;
    safe_free(action);
    action = NULL;
  } else {

     /*
      *  Place the now setup action structure onto the transformActionStack
      */

    pushBottomStack( &transformActionStack, (void *) action );
  
  }


}


#undef  ROUTINE
#define ROUTINE "ptrajSetupAnalyze()"

   void
ptrajSetupAnalyze(stackType *argumentStack, actionType type)
{
  analyzeInformation *analyze;
  char *buffer;
  int ierr;

     /*
      *  Make sure that this is indeed an "analyze" action (which really isn't
      *  necessary)
      */
  if (type != TRANSFORM_ANALYZE) {
    fprintf(stdout, "Error in ptrajSetupAnalyze(): Called with the wrong type!\n");
    fprintf(stdout, "Ignoring this command...\n");
    return;
  }

     /*
      *  Grab the first argument off the argument stack.  This is the "trigger" for
      *  the analyze function
      */
  buffer = getArgumentStringLower(&argumentStack, NULL);
  if (buffer == NULL) {
    fprintf(stdout, "ptrajSetupAnalyze(): No command passed to analyze, ignoring...\n");
    return;
  }

     /*
      *  Allocate and initialize a analyzeInformation structure and set the type
      */
  analyze = (analyzeInformation *)
    safe_malloc(sizeof(analyzeInformation));
  INITIALIZE_analyzeInformation(analyze);

     /*
      *  search for a match to the trigger (stored in buffer)
      */

  if (strncmp(buffer, "correlationcoe", 14) == 0) {
    analyze->type = ANALYZE_CORRELATIONCOEFFICIENT;
    analyze->fxn  = (analyzeFunction) analyzeCorrelationCoefficient;

  } else if (strncmp(buffer, "hbond", 5) == 0) {
    analyze->type = ANALYZE_HBOND;
    analyze->fxn  = (analyzeFunction) analyzeHBond;

  } else if (strncmp(buffer, "matrix", 6) == 0) {
    analyze->type = ANALYZE_MATRIX;
    analyze->fxn  = (analyzeFunction) analyzeMatrix;

  } else if (strncmp(buffer, "modes", 5) == 0) {
    analyze->type = ANALYZE_MODES;
    analyze->fxn  = (analyzeFunction) analyzeModes;

  } else if (strncmp(buffer, "set", 3) == 0) {
    analyze->type = ANALYZE_SET;
    analyze->fxn  = (analyzeFunction) analyzeSet;

  } else if (strncmp(buffer, "stat", 4) == 0) {
    analyze->type = ANALYZE_STATISTICS;
    analyze->fxn  = (analyzeFunction) analyzeStatistics;

  } else if (strncmp(buffer, "timecorr", 8) == 0) {
    analyze->type = ANALYZE_TIMECORR;
    analyze->fxn  = (analyzeFunction) analyzeTimecorr;

  } else if (strncmp(buffer, "test", 4) == 0) {
    analyze->type = ANALYZE_TEST;
    analyze->fxn  = (analyzeFunction) analyzeTest;

  } else {

    fprintf(stdout, "WARNING in ptrajSetupAnalyze(): unknown analyze type %i\n", type);
    safe_free(analyze);
    return;

  }

     /*
      *  Parse the arguments.  This is done by the analyze->fxn itself in the
      *  PTRAJ_SETUP mode.  The argumentStack is placed into the
      *  complex argument 1 slot for this.
      */

  ierr = 0;
  if (analyze->type != ANALYZE_NOOP) {

    analyze->carg1 = (void *) &argumentStack;
    ierr = analyze->fxn(analyze, scalarStack, PTRAJ_SETUP);

  }

     /*
      *  If the setup fails, -1 is returned and therefore this action
      *  should not be placed on the action stack and the associated
      *  memory should be freed
      */

  if (ierr < 0) {
    safe_free(analyze);
    analyze = NULL;
  } else {

     /*
      *  Place the now setup analyze structure onto the transformAnalyzeStack
      */

    pushBottomStack( &transformAnalyzeStack, (void *) analyze );
  
  }


}


#undef  ROUTINE
#define ROUTINE "ptraj()"

   void
ptraj(char *filenamep)
{
  FILE *infile;
  char buffer[BUFFER_SIZE];
  char *bufferp;
  coordinateInfo *currentCoordinateInfo;
  double *X, *Y, *Z;
  double box[6], boxnew[6];
  int boxfixed[6];
  actionInformation *action;
  analyzeInformation *analyze;
  int set;
  int local_set;
  int processed;
  int i;
  int suppressProcessing;
  stackType *actionStackTemp = NULL;
  int outputTrajectory = 0;
  stackType *sp, *argumentStack;
  ptrajState *startingState, *currentState, **statep;
  char *continuation;

  int readCoordinates;
  int processCoordinates;
  int firstOutput;


  /*
   *  --------------- INPUT FILE PROCESSING --------------------
   *
   *  if an input file was specified, open it up, else use
   *  standard input...
   */
  if ( filenamep == NULL || 
       strcmp(filenamep, "") == 0 || 
       strcmp(filenamep, "stdin") == 0 ) {

    fprintf(stdout, "\nPTRAJ: Processing input from \"STDIN\" ...\n");
    infile = stdin;
  } else if ( openFile(&infile, filenamep, "r") == 0 ) {
    fprintf(stdout, 
	    "WARNING in %s: Could not open input file (%s), exiting\n", ROUTINE,
	    filenamep);
    ptrajCleanup();
    return;
  } else {
    fprintf(stdout, "\nPTRAJ: Processing input from file %s\n", filenamep);
  }


  /*
   *  ----------------- SETUP INITIAL STATE --------------------
   *
   *  this gives a snapshot of the current "state" based on the 
   *  appropriate parameter/topology information (i.e. AMBER prmtop)
   *  upon entry.  Note this assumes that the GLOBAL ptrajState
   *  information was previously set in the main routine (main.c)
   *  via a call to ptrajInitializeState().
   */

  statep = ptrajCurrentState();
  startingState = *statep;
  currentState = startingState;

  /*
   *  --------------- PROCESS THE INPUT FILE -------------------
   *
   *  Read the input file, line by line, using the "dispatchToken()"
   *  routine (dispatch.c) to find a match for each command, ignoring
   *  comments in the input file.  Each "command" typed by the user 
   *  has an associated routine that is called for setup.  Currently,
   *  this is ptrajSetup() for "actions" and ptrajSetupIO for input/
   *  output functions.  Note that ptrajSetup() will actually call the
   *  "action" routine to perform the setup and parse arguments.  Note
   *  also that the global state pointer may be altered!!!
   *
   *  Lines of text are processed until EOF is encountered.
   */

  argumentStack = NULL;
  while ( (bufferp = fgets(buffer, BUFFER_SIZE, infile)) != NULL &&
	  strcmp(bufferp, "go\n") != 0 ) {

    continuation = bufferp;
    while ( continuation != NULL ) {
      if (strlen(buffer) >= BUFFER_SIZE)
	continuation = NULL;
      else {
	continuation = strrchr(buffer, '\\');
   
	if (continuation)
	  continuation = fgets(continuation, (BUFFER_SIZE - strlen(buffer) - 1), infile);
      }
    }

    skipWhitespace(bufferp);

    /*
     *  skip blank lines and/or comments denoted by "#" or "%"
     */
    if (bufferp[0] != (char) 0 && 
	bufferp[0] != '#' &&
	bufferp[0] != '%') {  

      fprintf(stdout, "\nPTRAJ: %s", buffer);
      dispatchToken((Token *) &ptrajTokenlist, argumentStack, (char *) buffer);
    }
  }
  if (infile != stdin)
    safe_fclose(infile);

    /*
     *  ----------------- ERROR CHECKING ---------------------
     */

  if (transformFileStack == NULL) {
    fprintf(stdout, "WARNING in %s: No input trajectories specified (trajin), aborting...\n",ROUTINE);
    ptrajCleanup();
    return;
  }

  if (globalOutInfo == NULL) {
    fprintf(stdout, "[No output trajectory specified (trajout)]\n");
  } else {
    outputTrajectory = 1;
  }

  /*
   *  set default box information.  This will allow, at setup time, specification
   *  of the current box sizes (if none were specified such as can happen when
   *  reading CHARMM PSF files) and also the option to FIX box sizes.
   */
  if (currentState->IFBOX) {
    for (i=0; i<6; i++) {
      box[i] = currentState->box[i];
      boxfixed[i] = currentState->boxfixed[i];
    }
  } else {
    box[0] = 0.0;
    box[1] = 0.0;
    box[2] = 0.0;
    box[3] = 90.0;
    box[4] = 90.0;
    box[5] = 90.0;
    boxfixed[0] = 0; boxfixed[1] = 0; boxfixed[2] = 0;
    boxfixed[3] = 0; boxfixed[4] = 0; boxfixed[5] = 0;
  }

  /*
   * --------- CHECK HOW MANY FRAMES WILL BE PROCESSED -------
   */

  /*
   *  TO DO: FIX THIS AS IT IS BROKEN
   */

  startingState->maxFrames = 0;
  for (sp = transformFileStack; sp != NULL; sp = sp->next) {
    currentCoordinateInfo = (coordinateInfo *) sp->entry;

    startingState->maxFrames += (currentCoordinateInfo->stop - 
				 currentCoordinateInfo->start) /
                                 currentCoordinateInfo->offset + 1;
  }

  fprintf(stdout, "\nPTRAJ: Successfully read the input file.\n");
  fprintf(stdout, "       Coordinate processing will occur on %i frames.\n", 
	  startingState->maxFrames);
  fprintf(stdout, "       Summary of I/O and actions follows:\n\n");


  /*
   * ----- PRINT A SUMMARY OF THE FILES IN/OUT AND ACTIONS ------
   */

  fprintf(stdout, "INPUT COORDINATE FILES\n");
  printStack( &transformFileStack, printCoordinateInfo, NULL );

  if ( globalOutInfo == NULL ) {
    fprintf(stdout, "\nNO OUTPUT COORDINATE FILE WAS SPECIFIED\n");
  } else {
    fprintf(stdout, "\nOUTPUT COORDINATE FILE\n");
    printCoordinateInfo( (void *) globalOutInfo );
  }

  if (transformReferenceStack != NULL) {
    fprintf(stdout, "\nREFERENCE FILE\n");
    printStack( &transformReferenceStack, printCoordinateInfo, NULL );
  }

  if ( transformActionStack == NULL ) {

    fprintf(stdout, "\nNO ACTIONS WERE SPECIFIED\n");

  } else {

    fprintf(stdout, "\nACTIONS\n");
    i = 1;
    for (actionStackTemp = transformActionStack;
	 actionStackTemp != NULL;
	 actionStackTemp = actionStackTemp->next) {

      action = (actionInformation *) actionStackTemp->entry;
         /*
          *  With each action's state variable local, it is necessary
          *  to update each with the current states maxFrames value!!!
          */
      action->state->maxFrames = startingState->maxFrames;
         /*
          *  call the action function with the mode PTRAJ_STATUS
          */
      if (action->type != TRANSFORM_NOOP  &&
	  action->type != TRANSFORM_TRANSFORM) {
	fprintf(stdout, "%3i>", i);
	for (i=0; i < 6; i++)
	  boxnew[i] = box[i]; /* protect the current box info */
	action->fxn(action, X, Y, Z, boxnew, PTRAJ_STATUS);
	i++;
      }
    }
    fprintf(stdout, "\n");
  }

  if (transformAnalyzeStack != NULL) {
    fprintf(stdout, "\nANALYZE\n");
    i = 1;
    for (actionStackTemp = transformAnalyzeStack;
	 actionStackTemp != NULL;
	 actionStackTemp = actionStackTemp->next) {

      analyze = (analyzeInformation *) actionStackTemp->entry;
      if (analyze->type != ANALYZE_NOOP) {
	fprintf(stdout, "%3i>", i);
	analyze->fxn(analyze, scalarStack, PTRAJ_STATUS);
	i++;
      }
    }
    fprintf(stdout, "\n");
  }

  /*
   *  ---------- ALLOCATE SPACE FOR COORDINATES ------------
   *
   *  perform and initial setup necessary prior to reading in
   *  the coordinates
   */

  X = (double *) safe_malloc(sizeof(double) * startingState->atoms);
  Y = (double *) safe_malloc(sizeof(double) * startingState->atoms);
  Z = (double *) safe_malloc(sizeof(double) * startingState->atoms);

  /*
y   *  --------- MAIN LOOP FOR COORDINATE PROCESSING -------------
   *
   *  loop over each of the files representing the coordinates.
   *
   *  "set"        -- the global counter over all sets, all files
   *  "local_set"  -- the counter over sets in each individual file
   */

  processed = 0;
  firstOutput = 1;
  set = 1;
  for (sp = transformFileStack; sp != NULL; sp = sp->next) {

    currentCoordinateInfo = (coordinateInfo *) sp->entry;

    /*
     *  ------------- PREPROCESS COORDINATE FILES-----------------
     */

       /*
        *  open up the file and preprocess the coordinates
        */

    if ( ptrajPreprocessInputCoordinates(currentCoordinateInfo) ) {
      ptrajCleanup();
      return;
    }

    local_set = 1;

    /*
     *  -------- READ IN, PROCESS, AND OUTPUT COORDINATES --------
     *
     *  The following two variables control whether to continue
     *  reading configurations from this file or not and also
     *  whether the coordinates should be processed after this
     *  read...
     *
     *  "readCoordinates"    > 0 if there are more coordinate sets 
     *                           to read in this file
     *
     *  "processCoordinates" > 0 if this coordinate set should be
     *                           processed
     */

    readCoordinates = 1;
    processCoordinates = 1;

    while ( readCoordinates ) {

      /*
       *  ADD STOP/START CHECK!!!
       */


         /*
          *  read in the current file of coordinates, a single set at a time.
          */
      for (i=0; i<6; i++)
	boxnew[i] = box[i];

      ptrajProcessInputCoordinates(currentCoordinateInfo, startingState, X, Y, Z, boxnew,
				   local_set, &readCoordinates, &processCoordinates);

      for (i=0; i<6; i++)
	if (boxfixed[i] == 0)
	  box[i] = boxnew[i];

      
         /*
          *  process coordinates if necessary
          */

      if ( processCoordinates ) {

           /*
            *  check to see if this snapshot is within bounds/offset
            */

	if ((local_set >= currentCoordinateInfo->start &&
	     local_set <= currentCoordinateInfo->stop) &&
	    ((currentCoordinateInfo->offset == 1) ||
	     ((local_set - currentCoordinateInfo->start) %
	     currentCoordinateInfo->offset == 0)) ) {

             /*
	      *  TRAVERSE THE ACTION STACK to perform each
              *  action on each coordinate set.  Note that a 
	      *  particular action can suppress processing of
	      *  further actions on the stack and prevent output
	      *  by setting the suppressProcessing flag in the
	      *  actionInformation structure.  This is useful 
	      *  when various sets are to be accumulated (for
	      *  example when calculating running averages) prior
	      *  to output and further processing...
              */
	  suppressProcessing = 0;
	  for (actionStackTemp = transformActionStack;
	       actionStackTemp != NULL;
	       actionStackTemp = actionStackTemp->next) {
	    action = (actionInformation *)
	      actionStackTemp->entry;
	    if (action->type != TRANSFORM_NOOP &&
		action->type != TRANSFORM_TRANSFORM &&
		suppressProcessing == 0) {

	      for (i=0; i<6; i++)
		boxnew[i] = box[i]; /* protect box coordinates */

	      action->fxn(action, X, Y, Z, boxnew, PTRAJ_ACTION);

              for (i=0; i<6; i++)
                if (boxfixed[i] == 0)
                  box[i] = boxnew[i];


                 /*
                  *  update the current state
                  */
	      currentState = action->state;
                 /*
                  *  check if any of the actions have suppressed output
                  */
	      if (suppressProcessing == 0)
		suppressProcessing = action->suppressProcessing;
	    }
	  }

             /*
              *  perform output as necessary
              */

	  if ( outputTrajectory && suppressProcessing == 0) {

	    ptrajOutputCoordinates(globalOutInfo, currentState, set, globalOutInfo->append, firstOutput, 0,
				   currentState->atoms, X, Y, Z, box);
	    firstOutput = 0;

	  }

             /*
              *  print out a dot for each snapshot processed to show progress
              */
	  processed++;
	  if ( local_set % 50 == 0 || local_set == 1 )
	    fprintf(stdout, "\nSet %6i ", local_set);
	  if ( local_set % 50 != 0 ) 
	    fprintf(stdout, "."); 

	} else {

             /*
              *  for coordinates that were not processed, print a space
              */
	  if ( local_set % 50 == 0 || local_set == 1 )
	    fprintf(stdout, "\nSet %6i ", local_set);
	  if ( local_set % 50 != 0 ) 
	    fprintf(stdout, " "); 

	} 
	fflush(stdout);
	set += 1;
	local_set += 1;

      }   /* IF (processCoordinates) */
    }     /* WHILE (readCoordinates) */
  }     /* FOR over transformFileStack */


  /*
   *  -------------- FINAL POSTPROCESSING -----------------
   */

  if (outputTrajectory) {
    ptrajOutputCoordinates(globalOutInfo, currentState, set, globalOutInfo->append, 0, 1, 
			   currentState->atoms, NULL, NULL, NULL, NULL);
  }

  /*
   *  -------------- DUMP ACCUMULATED DATA -------------------
   */

  fprintf(stdout, "\n\nPTRAJ: Successfully read in %i sets and processed %i sets.\n",
	  set-1, processed);
  fprintf(stdout, "       Dumping accumulated results (if any)\n\n");
  for (actionStackTemp = transformActionStack;
       actionStackTemp != NULL;
       actionStackTemp = actionStackTemp->next) {

    action = (actionInformation *) actionStackTemp->entry;
    if (action->type != TRANSFORM_NOOP &&
	action->type != TRANSFORM_TRANSFORM) {

      for (i=0; i<6; i++)
	boxnew[i] = box[i];
      action->fxn(action, X, Y, Z, boxnew, PTRAJ_PRINT);
    }
  }

  /*
   *  ----------- PERFORM ANY REQUESTED ANALYSIS -------------
   */

  if (transformAnalyzeStack != NULL) {

    fprintf(stdout, "\nPTRAJ: Analyzing accumulated data\n");
    for (actionStackTemp = transformAnalyzeStack;
	 actionStackTemp != NULL;
	 actionStackTemp = actionStackTemp->next) {

      analyze = (analyzeInformation *) actionStackTemp->entry;
      if (analyze->type)
	analyze->fxn(analyze, scalarStack, PTRAJ_ACTION);
    }
  }

  ptrajCleanup();
  ptrajClearState(&startingState);
  safe_free(X);
  safe_free(Y);
  safe_free(Z);
  X = NULL;
  Y = NULL;
  Z = NULL;

}




