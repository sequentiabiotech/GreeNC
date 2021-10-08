/**
 * Source code making the folding of a sequence and all its subsequences
 * containing a putative miRNA sequence.
 * This source code is a modified version of RNAfold developped by Ivo Hofacker,
 * Chrisoph Flamm (original implementation by Walter Fontana). 
 *
 *
 * The source code
 * has been modified by Anthony Mathelier.
*/

/* Original header : */
/* Last changed Time-stamp: <2005-07-24 10:54:43 ivo> */
/*                
		  minimum free energy
		  RNA secondary structure prediction

		  c Ivo Hofacker, Chrisoph Flamm
		  original implementation by
		  Walter Fontana
		  
		  Vienna RNA package
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"

/*@unused@*/
/* static char rcsid[] UNUSED = "$Id: fold.c,v 1.33 2005/08/07 11:44:43 ivo Exp $";
 */

#define PAREN

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1		/* stacking energies for bulges of size 1 */
#define NEW_NINIO     1		/* new asymetry penalty */

PUBLIC void fold (const char *string, char *structure);
PUBLIC float energy_of_struct (const char *string, const char *structure);
PUBLIC int energy_of_struct_pt (const char *string, short *ptable,
				short *s, short *s1);
PUBLIC void free_arrays (void);
PUBLIC void initialize_fold (int length);
PUBLIC void update_fold_params (void);

PUBLIC int logML = 0;		/* if nonzero use logarithmic ML energy in
				   energy_of_struct */
PUBLIC int uniq_ML = 0;		/* do ML decomposition uniquely (for subopt) */
/*@unused@*/
/* PRIVATE void  letter_structure(char *structure, int length) UNUSED;
 */
PRIVATE void parenthesis_structure (char *structure, int length);
PRIVATE void get_arrays (unsigned int size);
/* PRIVATE void  scale_parameters(void); */
PRIVATE int stack_energy (int i, const char *string);
PRIVATE int ML_Energy (int i, int is_extloop);
PRIVATE void make_ptypes (const short *S, const char *structure);
PRIVATE void encode_seq (const char *sequence);
PRIVATE int backtrack (const char *sequence, int s, int begin, int length);
PRIVATE void fill_arrays (const char *sequence);
/*@unused@*/
inline PRIVATE int oldLoopEnergy (int i, int j, int p, int q, int type,
				  int type_2);
inline int LoopEnergy (int n1, int n2, int type, int type_2, int si1,
		       int sj1, int sp1, int sq1);
inline int HairpinE (int size, int type, int si1, int sj1,
		     const char *string);

#define MAXSECTORS      500	/* dimension for a backtrack array */
#define LOCALITY        0.	/* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))

PRIVATE paramT *P = NULL;

PRIVATE int *indx;		/* index for moving in the triangle matrices c[] and fMl[] */

PRIVATE int *c;			/* energy array, given that i-j pair */
PRIVATE int *cc;		/* linear array for calculating canonical structures */
PRIVATE int *cc1;		/*   "     "        */
/**
 *
 * 	Modifications made by A.Mathelier
 *
**/

PRIVATE int **f;		/* min energy tables of (i,j) */

/**
 *
 * 	End of modifications
 *
**/
PRIVATE int *f5;		/* energy of 5' end */
PRIVATE int *fML;		/* multi-loop auxiliary energy array */
PRIVATE int *fM1;		/* second ML array, only for subopt */
PRIVATE int *Fmi;		/* holds row i of fML (avoids jumps in memory) */
PRIVATE int *DMLi;		/* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int *DMLi1;		/*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int *DMLi2;		/*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE char *ptype;		/* precomputed array of pair types */
PRIVATE short *S, *S1;
PRIVATE int init_length = -1;

/* needed by cofold/eval */
PRIVATE int cut_in_loop (int i);
PUBLIC int cut_point = -1;	/* set to first pos of second seq for cofolding */
PUBLIC int eos_debug = 0;	/* verbose info from energy_of_struct */

/*--------------------------------------------------------------------------*/

void
initialize_fold (int length)
{
  unsigned int n;
  if (length < 1)
    nrerror ("initialize_fold: argument must be greater 0");
  if (init_length > 0)
    free_arrays ();
  get_arrays ((unsigned) length);
  init_length = length;

  for (n = 1; n <= (unsigned) length; n++)
    indx[n] = (n * (n - 1)) >> 1;	/* n(n-1)/2 */

  update_fold_params ();
}

/*--------------------------------------------------------------------------*/

PRIVATE void
get_arrays (unsigned int size)
{
  indx = (int *) space (sizeof (int) * (size + 1));
  c = (int *) space (sizeof (int) * ((size * (size + 1)) / 2 + 2));
  fML = (int *) space (sizeof (int) * ((size * (size + 1)) / 2 + 2));
  if (uniq_ML)
    fM1 = (int *) space (sizeof (int) * ((size * (size + 1)) / 2 + 2));

  ptype = (char *) space (sizeof (char) * ((size * (size + 1)) / 2 + 2));

  /**
   *
   * 	Modifications made by A.Mathelier
   *
  **/

  f = (int **) space (sizeof (int *) * (size + 1));
  int a;
  for (a = 0; a < size + 1; a++)
    f[a] = (int *) space (sizeof (int) * (size + 1));

  /**
   *
   * 	End of modifications
   *
  **/

  f5 = (int *) space (sizeof (int) * (size + 2));
  cc = (int *) space (sizeof (int) * (size + 2));
  cc1 = (int *) space (sizeof (int) * (size + 2));
  Fmi = (int *) space (sizeof (int) * (size + 1));
  DMLi = (int *) space (sizeof (int) * (size + 1));
  DMLi1 = (int *) space (sizeof (int) * (size + 1));
  DMLi2 = (int *) space (sizeof (int) * (size + 1));
  base_pair = (struct bond *) space (sizeof (struct bond) * (1 + size / 2));
}

/*--------------------------------------------------------------------------*/

void
free_arrays (void)
{
  free (indx);
  free (c);
  free (fML);
  free (f5);
  free (cc);
  free (cc1);
  free (ptype);
  if (fM1)
    {
      free (fM1);
      fM1 = NULL;
    }

  free (base_pair);
  base_pair = NULL;
  free (Fmi);
  free (DMLi);
  free (DMLi1);
  free (DMLi2);
  init_length = 0;
}

/*--------------------------------------------------------------------------*/

void
export_fold_arrays (int **f5_p, int **c_p, int **fML_p, int **fM1_p,
		    int **indx_p, char **ptype_p)
{
  /* make the DP arrays available to routines such as subopt() */
  *f5_p = f5;
  *c_p = c;
  *fML_p = fML;
  *fM1_p = fM1;
  *indx_p = indx;
  *ptype_p = ptype;
}

/*--------------------------------------------------------------------------*/

PRIVATE int *BP;		/* contains the structure constrainsts: BP[i]
				   -1: | = base must be paired
				   -2: < = base must be paired with j<i
				   -3: > = base must be paired with j>i
				   -4: x = base must not pair
				   positive int: base is paired with int      */

void
fold (const char *string, char *structure)
{
  int length;

  length = (int) strlen (string);
  if (length > init_length)
    initialize_fold (length);
  if (fabs (P->temperature - temperature) > 1e-6)
    update_fold_params ();

  encode_seq (string);

  BP = (int *) space (sizeof (int) * (length + 2));
  make_ptypes (S, structure);

  fill_arrays (string);
}

PRIVATE void
fill_arrays (const char *string)
{
  /* fill "c", "fML" and "f5" arrays and return  optimal energy */

  int i, j, k, length, energy;
  int decomp, new_fML, max_separation;
  int no_close, type, type_2, tt;
  int bonus = 0;

  length = (int) strlen (string);

  max_separation = (int) ((1. - LOCALITY) * (double) (length - 2));	/* not in use */

  for (j = 1; j <= length; j++)
    {
      Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;
    }

  for (j = 1; j <= length; j++)
    for (i = (j > TURN ? (j - TURN) : 1); i < j; i++)
      {
	c[indx[j] + i] = fML[indx[j] + i] = INF;
	if (uniq_ML)
	  fM1[indx[j] + i] = INF;
      }

  for (i = length - TURN - 1; i >= 1; i--)
    {				/* i,j in [1..length] */

      for (j = i + TURN + 1; j <= length; j++)
	{
	  int p, q, ij;
	  ij = indx[j] + i;
	  bonus = 0;
	  type = ptype[ij];

	  /* enforcing structure constraints */
	  if ((BP[i] == j) || (BP[i] == -1) || (BP[i] == -2))
	    bonus -= BONUS;
	  if ((BP[j] == -1) || (BP[j] == -3))
	    bonus -= BONUS;
	  if ((BP[i] == -4) || (BP[j] == -4))
	    type = 0;

	  no_close = (((type == 3) || (type == 4)) && no_closingGU
		      && (bonus == 0));

	  if (j - i - 1 > max_separation)
	    type = 0;		/* forces locality degree */

	  if (type)
	    {			/* we have a pair */
	      int new_c = 0, stackEnergy = INF;
	      /* hairpin ---------------------------------------------- */

	      if (no_close)
		new_c = FORBIDDEN;
	      else
		new_c =
		  HairpinE (j - i - 1, type, S1[i + 1], S1[j - 1],
			    string + i - 1);

	/*--------------------------------------------------------
	  check for elementary structures involving more than one
	  closing pair.
	  --------------------------------------------------------*/

	      for (p = i + 1; p <= MIN2 (j - 2 - TURN, i + MAXLOOP + 1); p++)
		{
		  int minq = j - i + p - MAXLOOP - 2;
		  if (minq < p + 1 + TURN)
		    minq = p + 1 + TURN;
		  for (q = minq; q < j; q++)
		    {
		      type_2 = ptype[indx[q] + p];

		      if (type_2 == 0)
			continue;
		      type_2 = rtype[type_2];

		      if (no_closingGU)
			if (no_close || (type_2 == 3) || (type_2 == 4))
			  if ((p > i + 1) || (q < j - 1))
			    continue;	/* continue unless stack */

#if 1
		      energy =
			LoopEnergy (p - i - 1, j - q - 1, type, type_2,
				    S1[i + 1], S1[j - 1], S1[p - 1],
				    S1[q + 1]);
#else
		      /* duplicated code is faster than function call */

#endif
		      new_c = MIN2 (energy + c[indx[q] + p], new_c);
		      if ((p == i + 1) && (j == q + 1))
			stackEnergy = energy;	/* remember stack energy */

		    }		/* end q-loop */
		}		/* end p-loop */

	      /* multi-loop decomposition ------------------------ */


	      if (!no_close)
		{
		  int MLenergy;
		  decomp = DMLi1[j - 1];
		  if (dangles)
		    {
		      int d3 = 0, d5 = 0;
		      tt = rtype[type];
		      d3 = P->dangle3[tt][S1[i + 1]];
		      d5 = P->dangle5[tt][S1[j - 1]];
		      if (dangles == 2)	/* double dangles */
			decomp += d5 + d3;
		      else
			{	/* normal dangles */
			  decomp =
			    MIN2 (DMLi2[j - 1] + d3 + P->MLbase, decomp);
			  decomp =
			    MIN2 (DMLi1[j - 2] + d5 + P->MLbase, decomp);
			  decomp =
			    MIN2 (DMLi2[j - 2] + d5 + d3 +
				  2 * P->MLbase, decomp);
			}
		    }

		  MLenergy = P->MLclosing + P->MLintern[type] + decomp;

		  new_c = MLenergy < new_c ? MLenergy : new_c;
		}

	      /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

	      if (dangles == 3)
		{
		  decomp = INF;
		  for (k = i + 2 + TURN; k < j - 2 - TURN; k++)
		    {
		      type_2 = ptype[indx[k] + i + 1];
		      type_2 = rtype[type_2];
		      if (type_2)
			decomp =
			  MIN2 (decomp,
				c[indx[k] + i + 1] +
				P->stack[type][type_2] +
				fML[indx[j - 1] + k + 1]);
		      type_2 = ptype[indx[j - 1] + k + 1];
		      type_2 = rtype[type_2];
		      if (type_2)
			decomp =
			  MIN2 (decomp,
				c[indx[j - 1] + k + 1] +
				P->stack[type][type_2] + fML[indx[k] +
							     i + 1]);
		    }
		  /* no TermAU penalty if coax stack */
		  decomp += 2 * P->MLintern[1] + P->MLclosing;
		  new_c = MIN2 (new_c, decomp);
		}

	      new_c = MIN2 (new_c, cc1[j - 1] + stackEnergy);
	      cc[j] = new_c + bonus;
	      if (noLonelyPairs)
		c[ij] = cc1[j - 1] + stackEnergy + bonus;
	      else
		c[ij] = cc[j];

	    }
	  /* end >> if (pair) << */
	  else
	    c[ij] = INF;


	  /* done with c[i,j], now compute fML[i,j] */
	  /* free ends ? ----------------------------------------- */

	  new_fML = fML[ij + 1] + P->MLbase;
	  new_fML = MIN2 (fML[indx[j - 1] + i] + P->MLbase, new_fML);
	  energy = c[ij] + P->MLintern[type];
	  if (dangles == 2)
	    {			/* double dangles */
	      energy += (i == 1) ?	/* works also for circfold */
		P->dangle5[type][S1[length]] : P->dangle5[type][S1[i - 1]];
	      /* if (j<length) */ energy += P->dangle3[type][S1[j + 1]];
	    }
	  new_fML = MIN2 (energy, new_fML);
	  if (uniq_ML)
	    fM1[ij] = MIN2 (fM1[indx[j - 1] + i] + P->MLbase, energy);

	  if (dangles % 2 == 1)
	    {			/* normal dangles */
	      tt = ptype[ij + 1];	/* i+1,j */
	      new_fML = MIN2 (c[ij + 1] + P->dangle5[tt][S1[i]]
			      + P->MLintern[tt] + P->MLbase, new_fML);
	      tt = ptype[indx[j - 1] + i];
	      new_fML = MIN2 (c[indx[j - 1] + i] + P->dangle3[tt][S1[j]]
			      + P->MLintern[tt] + P->MLbase, new_fML);
	      tt = ptype[indx[j - 1] + i + 1];
	      new_fML =
		MIN2 (c[indx[j - 1] + i + 1] + P->dangle5[tt][S1[i]] +
		      P->dangle3[tt][S1[j]] + P->MLintern[tt] +
		      2 * P->MLbase, new_fML);
	    }

	  /* modular decomposition ------------------------------- */

	  for (decomp = INF, k = i + 1 + TURN; k <= j - 2 - TURN; k++)
	    decomp = MIN2 (decomp, Fmi[k] + fML[indx[j] + k + 1]);

	  DMLi[j] = decomp;	/* store for use in ML decompositon */
	  new_fML = MIN2 (new_fML, decomp);

	  /* coaxial stacking */
	  if (dangles == 3)
	    {
	      /* additional ML decomposition as two coaxially stacked helices */
	      for (decomp = INF, k = i + 1 + TURN; k <= j - 2 - TURN; k++)
		{
		  type = ptype[indx[k] + i];
		  type = rtype[type];
		  type_2 = ptype[indx[j] + k + 1];
		  type_2 = rtype[type_2];
		  if (type && type_2)
		    decomp = MIN2 (decomp,
				   c[indx[k] + i] + c[indx[j] + k + 1] +
				   P->stack[type][type_2]);
		}

	      decomp += 2 * P->MLintern[1];	/* no TermAU penalty if coax stack */
#if 0
	      /* This is needed for Y shaped ML loops with coax stacking of
	         interior pairts, but backtracking will fail if activated */
	      DMLi[j] = MIN2 (DMLi[j], decomp);
	      DMLi[j] = MIN2 (DMLi[j], DMLi[j - 1] + P->MLbase);
	      DMLi[j] = MIN2 (DMLi[j], DMLi1[j] + P->MLbase);
	      new_fML = MIN2 (new_fML, DMLi[j]);
#endif
	      new_fML = MIN2 (new_fML, decomp);
	    }

	  fML[ij] = Fmi[j] = new_fML;	/* substring energy */

	}

      {
	int *FF;		/* rotate the auxilliary arrays */
	FF = DMLi2;
	DMLi2 = DMLi1;
	DMLi1 = DMLi;
	DMLi = FF;
	FF = cc1;
	cc1 = cc;
	cc = FF;
	for (j = 1; j <= length; j++)
	  {
	    cc[j] = Fmi[j] = DMLi[j] = INF;
	  }
      }
    }

  /* calculate energies of 5' and 3' fragments */

      /**
       *
       * 	Modifications made by A.Mathelier
       *
      **/
  f5[TURN + 1] = 0;
  for (i = 0; i <= length; i++)
    f[i][TURN + 1] = 0;

  for (j = TURN + 2; j <= length; j++)
    {
      for (i = j - TURN - 1; i >= 1; i--)
	{
	  f5[j] = f5[j - 1];
	  f[i][j] = f[i][j - 1];
	  type = ptype[indx[j] + i];
	  if (type)
	    {
	      energy = c[indx[j] + i];
	      if (type > 2)
		energy += P->TerminalAU;
	      if ((dangles == 2) && (j < length))	/* double dangles */
		energy += P->dangle3[type][S1[j + 1]];
	      f5[j] = MIN2 (f5[j], energy);
	      f[i][j] = MIN2 (f[i][j], energy);
	    }
	  type = ptype[indx[j - 1] + i];
	  if ((type) && (dangles % 2 == 1))
	    {
	      energy = c[indx[j - 1] + i] + P->dangle3[type][S1[j]];
	      if (type > 2)
		energy += P->TerminalAU;
	      f5[j] = MIN2 (f5[j], energy);
	      f[i][j] = MIN2 (f[i][j], energy);
	    }
	  int k, energy2;
	  for (k = j - TURN - 1; k > i; k--)
	    {
	      type = ptype[indx[j] + k];
	      if (type)
		{
		  energy = f5[k - 1] + c[indx[j] + k];
		  energy2 = f[i][k - 1] + c[indx[j] + k];
		  if (type > 2)
		    {
		      energy += P->TerminalAU;
		      energy2 += P->TerminalAU;
		    }
		  if (dangles == 2)
		    {
		      energy += P->dangle5[type][S1[k - 1]];
		      energy2 += P->dangle5[type][S1[k - 1]];
		      if (j < length)
			{
			  energy += P->dangle3[type][S1[j + 1]];
			  energy2 += P->dangle3[type][S1[j + 1]];
			}
		    }
		  f5[j] = MIN2 (f5[j], energy);
		  f[i][j] = MIN2 (f[i][j], energy2);
		  if (dangles % 2 == 1)
		    {
		      energy =
			f5[k - 2] + c[indx[j] + k] +
			P->dangle5[type][S1[k - 1]];
		      energy2 =
			f[i][k - 2] + c[indx[j] + k] +
			P->dangle5[type][S1[k - 1]];
		      if (type > 2)
			{
			  energy += P->TerminalAU;
			  energy2 += P->TerminalAU;
			}
		      f5[j] = MIN2 (f5[j], energy);
		      f[i][j] = MIN2 (f[i][j], energy2);
		    }
		}
	      type = ptype[indx[j - 1] + k];
	      if ((type) && (dangles % 2 == 1))
		{
		  energy = c[indx[j - 1] + k] + P->dangle3[type][S1[j]];
		  energy2 = c[indx[j - 1] + k] + P->dangle3[type][S1[j]];
		  if (type > 2)
		    {
		      energy += P->TerminalAU;
		      energy2 += P->TerminalAU;
		    }
		  f5[j] = MIN2 (f5[j], f5[k - 1] + energy);
		  f[i][j] = MIN2 (f[i][j], f[i][k - 1] + energy2);
		  f5[j] =
		    MIN2 (f5[j],
			  f5[k - 2] + energy + P->dangle5[type][S1[k - 1]]);
		  f[i][j] =
		    MIN2 (f[i][j],
			  f[i][k - 2] + energy2 +
			  P->dangle5[type][S1[k - 1]]);
		}
	    }
	}
    }
/** base version**/

  /*f5[TURN+1]=0;
     for (j=TURN+2; j<=length; j++) {
     f5[j] = f5[j-1];
     type=ptype[indx[j]+1];
     if (type) {
     energy = c[indx[j]+1];
     if (type>2) energy += P->TerminalAU;
     if ((dangles==2)&&(j<length))  * double dangles *
     energy += P->dangle3[type][S1[j+1]];
     f5[j] = MIN2(f5[j], energy);
     }
     type=ptype[indx[j-1]+1];
     if ((type)&&(dangles%2==1)) {
     energy = c[indx[j-1]+1]+P->dangle3[type][S1[j]];
     if (type>2) energy += P->TerminalAU;
     f5[j] = MIN2(f5[j], energy);
     }
     for (i=j-TURN-1; i>1; i--) {
     type = ptype[indx[j]+i];
     if (type) {
     energy = f5[i-1]+c[indx[j]+i];
     if (type>2) energy += P->TerminalAU;
     if (dangles==2) {
     energy += P->dangle5[type][S1[i-1]];
     if (j<length) energy += P->dangle3[type][S1[j+1]];
     }
     f5[j] = MIN2(f5[j], energy);
     if (dangles%2==1) {
     energy = f5[i-2]+c[indx[j]+i]+P->dangle5[type][S1[i-1]];
     if (type>2) energy += P->TerminalAU;
     f5[j] = MIN2(f5[j], energy);
     }
     }
     type = ptype[indx[j-1]+i];
     if ((type)&&(dangles%2==1)) {
     energy = c[indx[j-1]+i]+P->dangle3[type][S1[j]];
     if (type>2) energy += P->TerminalAU;
     f5[j] = MIN2(f5[j], f5[i-1]+energy);
     f5[j] = MIN2(f5[j], f5[i-2]+energy+P->dangle5[type][S1[i-1]]);
     }
     return f5[length]; */

      /**
       *
       * 	End of modifications
       *
      **/

}

struct sect
{
  int i;
  int j;
  int ml;
} static sector[MAXSECTORS];	/* stack of partial structures for backtracking */

PRIVATE int
backtrack (const char *string, int s, int begin, int end)
{

  /*------------------------------------------------------------------
    trace back through the "c", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/

  /* normally s=0. 
     If s>0 then s items have been already pushed onto the sector stack */
  int i, j, k, energy, new;
  int no_close, type, type_2, tt;
  int bonus;
  int b = 0;

  if (s == 0)
    {
	  /**
	   *
	   * 	Modifications made by A.Mathelier
	   *
	  **/

      sector[++s].i = begin;

      /*base version */
/* 	sector[++s].i = 1; 
 */

      sector[s].j = end;
  /**
     *
     * 	End of modifications
     *
    **/

      sector[s].ml =
	(backtrack_type == 'M') ? 1 : ((backtrack_type == 'C') ? 2 : 0);
    }
  while (s > 0)
    {
      int ml, fij, fi, cij, traced, i1, j1, d3, d5, mm, p, q, jj = 0;
      int canonical = 1;	/* (i,j) closes a canonical structure */
      i = sector[s].i;
      j = sector[s].j;
      ml = sector[s--].ml;	/* ml is a flag indicating if backtracking is to 
				   occur in the fML- (1) or in the f-array (0) */
      if (ml == 2)
	{
	  base_pair[++b].i = i;
	  base_pair[b].j = j;
	  goto repeat1;
	}

      if (j < i + TURN + 1)
	continue;		/* no more pairs in this interval */

    /**
     *
     * 	Modifications made by A.Mathelier
     *
    **/

      fij = (ml) ? fML[indx[j] + i] : f[begin][j];
      fi = (ml) ? (fML[indx[j - 1] + i] + P->MLbase) : f[begin][j - 1];

      if (fij == fi)
	{			/* 3' end is unpaired */
	  sector[++s].i = i;
	  sector[s].j = j - 1;
	  sector[s].ml = ml;
	  continue;
	}
    /**
     *
     * 	End of modifications
     *
    **/

      if (ml == 0)
	{			/* backtrack in f5 */
	  /* j or j-1 is paired. Find pairing partner */
	  for (k = j - TURN - 1, traced = 0; k >= begin; k--)
	    {
	      int cc, en;
	      jj = k - 1;
	      type = ptype[indx[j - 1] + k];
	      if ((type) && (dangles % 2 == 1))
		{
		  cc = c[indx[j - 1] + k] + P->dangle3[type][S1[j]];
		  if (type > 2)
		    cc += P->TerminalAU;
		  if (fij == cc + f[begin][k - 1])
		    traced = j - 1;
		  if (k > i)
		    if (fij ==
			f[begin][k - 2] + cc + P->dangle5[type][S1[k - 1]])
		      {
			traced = j - 1;
			jj = k - 2;
		      }
		}
	      type = ptype[indx[j] + k];
	      if (type)
		{
		  cc = c[indx[j] + k];
		  if (type > 2)
		    cc += P->TerminalAU;
		  en = cc + f[begin][k - 1];
		  if (dangles == 2)
		    {
		      if (k > begin)
			en += P->dangle5[type][S1[k - 1]];
		      if (j < end)
			en += P->dangle3[type][S1[j + 1]];
		    }
		  if (fij == en)
		    traced = j;
		  if ((dangles % 2 == 1) && (k > begin))
		    if (fij ==
			f[begin][k - 2] + cc + P->dangle5[type][S1[k - 1]])
		      {
			traced = j;
			jj = k - 2;
		      }
		}
	      if (traced)
		break;
	    }

	  if (!traced)
	    nrerror ("backtrack failed in f5");
      /**
       *
       * 	Modifications made by A.Mathelier
       *
      **/

	  /*base version */
	  sector[++s].i = begin;

      /**
       *
       * 	End of modifications
       *
      **/
	  sector[s].j = jj;
	  sector[s].ml = ml;

	  i = k;
	  j = traced;
	  base_pair[++b].i = i;
	  base_pair[b].j = j;
	  goto repeat1;
	}

/** base version**/

/*     fij = (ml)? fML[indx[j]+i] : f5[j];
 *     fi  = (ml)?(fML[indx[j-1]+i]+P->MLbase):f5[j-1];
 * 
 *     if (fij == fi) {  // 3' end is unpaired 
 *       sector[++s].i = i;
 *       sector[s].j   = j-1;
 *       sector[s].ml  = ml;
 *       continue;
 *     }
 *      
 *     if (ml == 0) { // backtrack in f5 
 *       // j or j-1 is paired. Find pairing partner 
 * 
 *       int cc, en;
 *       jj = k-1; 
 *       type = ptype[indx[j-1]+k];
 *       if((type)&&(dangles%2==1)) {
 *         cc = c[indx[j-1]+k]+P->dangle3[type][S1[j]];
 *         if (type>2) cc += P->TerminalAU;
 *         if (fij == cc + f5[k-1]) 
 *           traced=j-1;
 *         if (k>i) 
 *           if (fij == f5[k-2] + cc + P->dangle5[type][S1[k-1]]) {
 *             traced=j-1; jj=k-2;
 *           }
 *       }
 *       type = ptype[indx[j]+k];
 *       if (type) {
 *         cc = c[indx[j]+k];
 *         if (type>2) cc += P->TerminalAU; 
 *         en = cc + f5[k-1];
 *         if (dangles==2) {
 * 
 *           if (fij == f5[k-2]+cc+P->dangle5[type][S1[k-1]]) {
 *             traced=j; jj=k-2;
 *           }
 *       }
 *       if (traced) break;
 *       }
 * 
 *       if (!traced) nrerror("backtrack failed in f5");
 *       sector[s].j   = jj;
 *       sector[s].ml  = ml;
 *        
 *       i=k; j=traced;
 *       base_pair[++b].i = i;
 *       base_pair[b].j   = j;
 *       goto repeat1;
 *     }
 */

/*     End of base version
 */

      else
	{			/* trace back in fML array */
	  int cij1 = INF, ci1j = INF, ci1j1 = INF;
	  if (fML[indx[j] + i + 1] + P->MLbase == fij)
	    {			/* 5' end is unpaired */
	      sector[++s].i = i + 1;
	      sector[s].j = j;
	      sector[s].ml = ml;
	      continue;
	    }

	  tt = ptype[indx[j] + i];
	  cij = c[indx[j] + i] + P->MLintern[tt];
	  if (dangles == 2)
	    {			/* double dangles, works also for circfold */
	      cij += (i == 1) ?
		P->dangle5[tt][S1[end]] : P->dangle5[tt][S1[i - 1]];
	      /* if (j<length) */ cij += P->dangle3[tt][S1[j + 1]];
	    }
	  else if (dangles % 2 == 1)
	    {			/* normal dangles */
	      tt = ptype[indx[j] + i + 1];
	      ci1j =
		c[indx[j] + i + 1] + P->dangle5[tt][S1[i]] +
		P->MLintern[tt] + P->MLbase;
	      tt = ptype[indx[j - 1] + i];
	      cij1 =
		c[indx[j - 1] + i] + P->dangle3[tt][S1[j]] +
		P->MLintern[tt] + P->MLbase;
	      tt = ptype[indx[j - 1] + i + 1];
	      ci1j1 =
		c[indx[j - 1] + i + 1] + P->dangle5[tt][S1[i]] +
		P->dangle3[tt][S1[j]] + P->MLintern[tt] + 2 * P->MLbase;
	    }

	  if ((fij == cij) || (fij == ci1j) || (fij == cij1)
	      || (fij == ci1j1))
	    {
	      /* found a pair */
	      if (fij == ci1j)
		i++;
	      else if (fij == cij1)
		j--;
	      else if (fij == ci1j1)
		{
		  i++;
		  j--;
		}
	      base_pair[++b].i = i;
	      base_pair[b].j = j;
	      goto repeat1;
	    }

	  for (k = i + 1 + TURN; k <= j - 2 - TURN; k++)
	    if (fij == (fML[indx[k] + i] + fML[indx[j] + k + 1]))
	      break;

	  if ((dangles == 3) && (k > j - 2 - TURN))
	    {			/* must be coax stack */
	      ml = 2;
	      for (k = i + 1 + TURN; k <= j - 2 - TURN; k++)
		{
		  type = ptype[indx[k] + i];
		  type = rtype[type];
		  type_2 = ptype[indx[j] + k + 1];
		  type_2 = rtype[type_2];
		  if (type && type_2)
		    if (fij ==
			c[indx[k] + i] + c[indx[j] + k + 1] +
			P->stack[type][type_2] + 2 * P->MLintern[1])
		      break;
		}
	    }

	  sector[++s].i = i;
	  sector[s].j = k;
	  sector[s].ml = ml;
	  sector[++s].i = k + 1;
	  sector[s].j = j;
	  sector[s].ml = ml;

	  if (k > j - 2 - TURN)
	    nrerror ("backtrack failed in fML");
	  continue;
	}

    repeat1:

    /*----- begin of "repeat:" -----*/
      if (canonical)
	cij = c[indx[j] + i];

      type = ptype[indx[j] + i];

      bonus = 0;

      if ((BP[i] == j) || (BP[i] == -1) || (BP[i] == -2))
	bonus -= BONUS;
      if ((BP[j] == -1) || (BP[j] == -3))
	bonus -= BONUS;

      if (noLonelyPairs)
	if (cij == c[indx[j] + i])
	  {
	    /* (i.j) closes canonical structures, thus
	       (i+1.j-1) must be a pair                */
	    type_2 = ptype[indx[j - 1] + i + 1];
	    type_2 = rtype[type_2];
	    cij -= P->stack[type][type_2] + bonus;
	    base_pair[++b].i = i + 1;
	    base_pair[b].j = j - 1;
	    i++;
	    j--;
	    canonical = 0;
	    goto repeat1;
	  }
      canonical = 1;


      no_close = (((type == 3) || (type == 4)) && no_closingGU
		  && (bonus == 0));
      if (no_close)
	{
	  if (cij == FORBIDDEN)
	    continue;
	}
      else
	if (cij ==
	    HairpinE (j - i - 1, type, S1[i + 1], S1[j - 1],
		      string + i - 1) + bonus)
	continue;

      for (p = i + 1; p <= MIN2 (j - 2 - TURN, i + MAXLOOP + 1); p++)
	{
	  int minq;
	  minq = j - i + p - MAXLOOP - 2;
	  if (minq < p + 1 + TURN)
	    minq = p + 1 + TURN;
	  for (q = j - 1; q >= minq; q--)
	    {

	      type_2 = ptype[indx[q] + p];
	      if (type_2 == 0)
		continue;
	      type_2 = rtype[type_2];
	      if (no_closingGU)
		if (no_close || (type_2 == 3) || (type_2 == 4))
		  if ((p > i + 1) || (q < j - 1))
		    continue;	/* continue unless stack */

	      /* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
	      energy = LoopEnergy (p - i - 1, j - q - 1, type, type_2,
				   S1[i + 1], S1[j - 1], S1[p - 1],
				   S1[q + 1]);

	      new = energy + c[indx[q] + p] + bonus;
	      traced = (cij == new);
	      if (traced)
		{
		  base_pair[++b].i = p;
		  base_pair[b].j = q;
		  i = p, j = q;
		  goto repeat1;
		}
	    }
	}

      /* end of repeat: -------------------------------------------------- */

      /* (i.j) must close a multi-loop */
      tt = rtype[type];
      mm = bonus + P->MLclosing + P->MLintern[tt];
      d5 = P->dangle5[tt][S1[j - 1]];
      d3 = P->dangle3[tt][S1[i + 1]];
      i1 = i + 1;
      j1 = j - 1;
      sector[s + 1].ml = sector[s + 2].ml = 1;

      for (k = i + 2 + TURN; k < j - 2 - TURN; k++)
	{
	  int en;
	  en = fML[indx[k] + i + 1] + fML[indx[j - 1] + k + 1] + mm;
	  if (dangles == 2)	/* double dangles */
	    en += d5 + d3;
	  if (cij == en)
	    break;
	  if (dangles % 2 == 1)
	    {			/* normal dangles */
	      if (cij ==
		  (fML[indx[k] + i + 2] + fML[indx[j - 1] + k + 1] + mm +
		   d3 + P->MLbase))
		{
		  i1 = i + 2;
		  break;
		}
	      if (cij ==
		  (fML[indx[k] + i + 1] + fML[indx[j - 2] + k + 1] + mm +
		   d5 + P->MLbase))
		{
		  j1 = j - 2;
		  break;
		}
	      if (cij ==
		  (fML[indx[k] + i + 2] + fML[indx[j - 2] + k + 1] + mm +
		   d3 + d5 + P->MLbase + P->MLbase))
		{
		  i1 = i + 2;
		  j1 = j - 2;
		  break;
		}
	    }
	  /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
	  /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
	  if (dangles == 3)
	    {
	      type_2 = ptype[indx[k] + i + 1];
	      type_2 = rtype[type_2];
	      if (type_2)
		{
		  en = c[indx[k] + i + 1] + P->stack[type][type_2] +
		    fML[indx[j - 1] + k + 1];
		  if (cij == en + 2 * P->MLintern[1] + P->MLclosing)
		    {
		      ml = 2;
		      sector[s + 1].ml = 2;
		      break;
		    }
		}
	      type_2 = ptype[indx[j - 1] + k + 1];
	      type_2 = rtype[type_2];
	      if (type_2)
		{
		  en = c[indx[j - 1] + k + 1] + P->stack[type][type_2] +
		    fML[indx[k] + i + 1];
		  if (cij == en + 2 * P->MLintern[1] + P->MLclosing)
		    {
		      sector[s + 2].ml = 2;
		      break;
		    }
		}
	    }

	}
      if (k <= j - 3 - TURN)
	{			/* found the decomposition */
	  sector[++s].i = i1;
	  sector[s].j = k;
	  sector[++s].i = k + 1;
	  sector[s].j = j1;
	}
      else
	{
#if 0
	  /* Y shaped ML loops fon't work yet */
	  if (dangles == 3)
	    {
	      /* (i,j) must close a Y shaped ML loop with coax stacking */
	      if (cij ==
		  fML[indx[j - 2] + i + 2] + mm + d3 + d5 + P->MLbase +
		  P->MLbase)
		{
		  i1 = i + 2;
		  j1 = j - 2;
		}
	      else if (cij == fML[indx[j - 2] + i + 1] + mm + d5 + P->MLbase)
		j1 = j - 2;
	      else if (cij == fML[indx[j - 1] + i + 2] + mm + d3 + P->MLbase)
		i1 = i + 2;
	      else /* last chance */ if (cij !=
					 fML[indx[j - 1] + i + 1] + mm +
					 P->MLbase)
		fprintf (stderr, "backtracking failed in repeat");
	      /* if we arrive here we can express cij via fML[i1,j1]+dangles */
	      sector[++s].i = i1;
	      sector[s].j = j1;
	    }
	  else
#endif
	    nrerror ("backtracking failed in repeat");
	}

    }

  base_pair[0].i = b;		/* save the total number of base pairs */

  return f[begin][end];
}

/*---------------------------------------------------------------------------*/

inline int
HairpinE (int size, int type, int si1, int sj1, const char *string)
{
  int energy;
  energy = (size <= 30) ? P->hairpin[size] :
    P->hairpin[30] + (int) (P->lxc * log ((size) / 30.));
  if (tetra_loop)
    if (size == 4)
      {				/* check for tetraloop bonus */
	char tl[7] = { 0 }, *ts;
	strncpy (tl, string, 6);
	if ((ts = strstr (P->Tetraloops, tl)))
	  energy += P->TETRA_ENERGY[(ts - P->Tetraloops) / 7];
      }
  if (size == 3)
    {
      char tl[6] = { 0, 0, 0, 0, 0, 0 }, *ts;
      strncpy (tl, string, 5);
      if ((ts = strstr (P->Triloops, tl)))
	energy += P->Triloop_E[(ts - P->Triloops) / 6];

      if (type > 2)		/* neither CG nor GC */
	energy += P->TerminalAU;	/* penalty for closing AU GU pair */
    }
  else				/* no mismatches for tri-loops */
    energy += P->mismatchH[type][si1][sj1];

  return energy;
}

/*---------------------------------------------------------------------------*/

inline PRIVATE int
oldLoopEnergy (int i, int j, int p, int q, int type, int type_2)
{
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int n1, n2, m, energy;
  n1 = p - i - 1;
  n2 = j - q - 1;

  if (n1 > n2)
    {
      m = n1;
      n1 = n2;
      n2 = m;
    }
  /* so that n2>=n1 */
  if (n2 == 0)
    energy = P->stack[type][type_2];	/* stack */

  else if (n1 == 0)
    {				/* bulge */
      energy = (n2 <= MAXLOOP) ? P->bulge[n2] :
	(P->bulge[30] + (int) (P->lxc * log (n2 / 30.)));

#if STACK_BULGE1
      if (n2 == 1)
	energy += P->stack[type][type_2];
#endif
    }
  else
    {				/* interior loop */

      if ((n1 + n2 == 2) && (james_rule))
	/* special case for loop size 2 */
	energy = P->int11[type][type_2][S1[i + 1]][S1[j - 1]];
      else
	{
	  energy = (n1 + n2 <= MAXLOOP) ? (P->internal_loop[n1 + n2]) :
	    (P->internal_loop[30] + (int) (P->lxc * log ((n1 + n2) / 30.)));

#if NEW_NINIO
	  energy += MIN2 (MAX_NINIO, (n2 - n1) * P->F_ninio[2]);
#else
	  m = MIN2 (4, n1);
	  energy += MIN2 (MAX_NINIO, ((n2 - n1) * P->F_ninio[m]));
#endif
	  energy += P->mismatchI[type][S1[i + 1]][S1[j - 1]] +
	    P->mismatchI[type_2][S1[q + 1]][S1[p - 1]];
	}
    }
  return energy;
}

/*--------------------------------------------------------------------------*/

inline int
LoopEnergy (int n1, int n2, int type, int type_2,
	    int si1, int sj1, int sp1, int sq1)
{
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, energy;

  if (n1 > n2)
    {
      nl = n1;
      ns = n2;
    }
  else
    {
      nl = n2;
      ns = n1;
    }

  if (nl == 0)
    return P->stack[type][type_2];	/* stack */

  if (ns == 0)
    {				/* bulge */
      energy = (nl <= MAXLOOP) ? P->bulge[nl] :
	(P->bulge[30] + (int) (P->lxc * log (nl / 30.)));
      if (nl == 1)
	energy += P->stack[type][type_2];
      else
	{
	  if (type > 2)
	    energy += P->TerminalAU;
	  if (type_2 > 2)
	    energy += P->TerminalAU;
	}
      return energy;
    }
  else
    {				/* interior loop */
      if (ns == 1)
	{
	  if (nl == 1)		/* 1x1 loop */
	    return P->int11[type][type_2][si1][sj1];
	  if (nl == 2)
	    {			/* 2x1 loop */
	      if (n1 == 1)
		energy = P->int21[type][type_2][si1][sq1][sj1];
	      else
		energy = P->int21[type_2][type][sq1][si1][sp1];
	      return energy;
	    }
	}
      else if (n1 == 2 && n2 == 2)	/* 2x2 loop */
	return P->int22[type][type_2][si1][sp1][sq1][sj1];
      {				/* generic interior loop (no else here!) */
	energy = (n1 + n2 <= MAXLOOP) ? (P->internal_loop[n1 + n2]) :
	  (P->internal_loop[30] + (int) (P->lxc * log ((n1 + n2) / 30.)));

	energy += MIN2 (MAX_NINIO, (nl - ns) * P->F_ninio[2]);

	energy += P->mismatchI[type][si1][sj1] +
	  P->mismatchI[type_2][sq1][sp1];
      }
    }
  return energy;
}


/*---------------------------------------------------------------------------*/

PRIVATE void
encode_seq (const char *sequence)
{
  unsigned int i, l;

  l = strlen (sequence);
  S = (short *) space (sizeof (short) * (l + 2));
  S1 = (short *) space (sizeof (short) * (l + 2));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = S1[0] = (short) l;

  for (i = 1; i <= l; i++)
    {				/* make numerical encoding of sequence */
      S[i] = (short) encode_char (toupper (sequence[i - 1]));
      S1[i] = alias[S[i]];	/* for mismatches of nostandard bases */
    }
  /* for circular folding add first base at position n+1 */
  S[l + 1] = S[1];
  S1[l + 1] = S1[1];
}

/*---------------------------------------------------------------------------*/

PRIVATE void
parenthesis_structure (char *structure, int length)
{
  int n, k;

  for (n = 0; n <= length - 1; structure[n++] = '.');
  structure[length] = '\0';

  for (k = 1; k <= base_pair[0].i; k++)
    {
      structure[base_pair[k].i - 1] = '(';
      structure[base_pair[k].j - 1] = ')';
    }
}

/*---------------------------------------------------------------------------*/

PUBLIC void
update_fold_params (void)
{
  P = scale_parameters ();
  make_pair_matrix ();
  if (init_length < 0)
    init_length = 0;
}

/*---------------------------------------------------------------------------*/
PRIVATE short *pair_table;

float
energy_of_struct (const char *string, const char *structure)
{
  int energy;
  short *ss, *ss1;

  if ((init_length < 0) || (P == NULL))
    update_fold_params ();
  if (fabs (P->temperature - temperature) > 1e-6)
    update_fold_params ();

  if (strlen (structure) != strlen (string))
    nrerror ("energy_of_struct: string and structure have unequal length");

  /* save the S and S1 pointers in case they were already in use */
  ss = S;
  ss1 = S1;
  encode_seq (string);

  pair_table = make_pair_table (structure);

  energy = energy_of_struct_pt (string, pair_table, S, S1);

  free (pair_table);
  free (S);
  free (S1);
  S = ss;
  S1 = ss1;
  return (float) energy / 100.;
}

int
energy_of_struct_pt (const char *string, short *ptable, short *s, short *s1)
{
  /* auxiliary function for kinfold,
     for most purposes call energy_of_struct instead */

  int i, length, energy;

  pair_table = ptable;
  S = s;
  S1 = s1;

  length = S[0];
  energy = backtrack_type == 'M' ? ML_Energy (0, 0) : ML_Energy (0, 1);
  if (eos_debug > 0)
    printf ("External loop                           : %5d\n", energy);
  for (i = 1; i <= length; i++)
    {
      if (pair_table[i] == 0)
	continue;
      energy += stack_energy (i, string);
      i = pair_table[i];
    }
  for (i = 1; !SAME_STRAND (i, length); i++)
    {
      if (!SAME_STRAND (i, pair_table[i]))
	{
	  energy += P->DuplexInit;
	  break;
	}
    }
  return energy;
}

/*---------------------------------------------------------------------------*/
PRIVATE int
stack_energy (int i, const char *string)
{
  /* calculate energy of substructure enclosed by (i,j) */
  int ee, energy = 0;
  int j, p, q, type;

  j = pair_table[i];
  type = pair[S[i]][S[j]];
  if (type == 0)
    {
      type = 7;
      if (eos_debug >= 0)
	fprintf (stderr,
		 "WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
		 string[i - 1], string[j - 1]);
    }

  p = i;
  q = j;
  while (p < q)
    {				/* process all stacks and interior loops */
      int type_2;
      while (pair_table[++p] == 0);
      while (pair_table[--q] == 0);
      if ((pair_table[q] != (short) p) || (p > q))
	break;
      type_2 = pair[S[q]][S[p]];
      if (type_2 == 0)
	{
	  type_2 = 7;
	  if (eos_debug >= 0)
	    fprintf (stderr,
		     "WARNING: bases %d and %d (%c%c) can't pair!\n", p,
		     q, string[p - 1], string[q - 1]);
	}
      /* energy += LoopEnergy(i, j, p, q, type, type_2); */
      if (SAME_STRAND (i, p) && SAME_STRAND (q, j))
	ee = LoopEnergy (p - i - 1, j - q - 1, type, type_2,
			 S1[i + 1], S1[j - 1], S1[p - 1], S1[q + 1]);
      else
	ee = ML_Energy (cut_in_loop (i), 1);
      if (eos_debug > 0)
	printf ("Interior loop (%3d,%3d) %c%c; (%3d,%3d) %c%c: %5d\n",
		i, j, string[i - 1], string[j - 1], p, q, string[p - 1],
		string[q - 1], ee);
      energy += ee;
      i = p;
      j = q;
      type = rtype[type_2];
    }				/* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p > q)
    {				/* hair pin */
      if (SAME_STRAND (i, j))
	ee = HairpinE (j - i - 1, type, S1[i + 1], S1[j - 1], string + i - 1);
      else
	ee = ML_Energy (cut_in_loop (i), 1);
      energy += ee;
      if (eos_debug > 0)
	printf ("Hairpin  loop (%3d,%3d) %c%c              : %5d\n",
		i, j, string[i - 1], string[j - 1], ee);

      return energy;
    }

  /* (i,j) is exterior pair of multiloop */
  while (p < j)
    {
      /* add up the contributions of the substructures of the ML */
      energy += stack_energy (p, string);
      p = pair_table[p];
      /* search for next base pair in multiloop */
      while (pair_table[++p] == 0);
    }
  {
    int ii;
    ii = cut_in_loop (i);
    ee = (ii == 0) ? ML_Energy (i, 0) : ML_Energy (ii, 1);
  }
  energy += ee;
  if (eos_debug > 0)
    printf ("Multi    loop (%3d,%3d) %c%c              : %5d\n",
	    i, j, string[i - 1], string[j - 1], ee);

  return energy;
}

/*---------------------------------------------------------------------------*/

PRIVATE int
ML_Energy (int i, int is_extloop)
{
  /* i is the 5'-base of the closing pair (or 0 for exterior loop)
     loop is scored as ML if extloop==0 else as exterior loop

     since each helix can coaxially stack with at most one of its
     neighbors we need an auxiliarry variable  cx_energy
     which contains the best energy given that the last two pairs stack.
     energy  holds the best energy given the previous two pairs do not
     stack (i.e. the two current helices may stack)
     We don't allow the last helix to stack with the first, thus we have to
     walk around the Loop twice with two starting points and take the minimum
   */

  int energy, cx_energy, best_energy = INF;
  int i1, j, p, q, u, x, type, count;
  int mlintern[NBPAIRS + 1], mlclosing, mlbase;

  if (is_extloop)
    {
      for (x = 0; x <= NBPAIRS; x++)
	mlintern[x] = P->MLintern[x] - P->MLintern[1];	/* 0 or TerminalAU */
      mlclosing = mlbase = 0;
    }
  else
    {
      for (x = 0; x <= NBPAIRS; x++)
	mlintern[x] = P->MLintern[x];
      mlclosing = P->MLclosing;
      mlbase = P->MLbase;
    }

  for (count = 0; count < 2; count++)
    {				/* do it twice */
      int ld5 = 0;		/* 5' dangle energy on prev pair (type) */
      if (i == 0)
	{
	  j = pair_table[0] + 1;
	  type = 0;		/* no pair */
	}
      else
	{
	  j = pair_table[i];
	  type = pair[S[j]][S[i]];
	  if (type == 0)
	    type = 7;

	  if (dangles == 3)
	    {			/* prime the ld5 variable */
	      if (SAME_STRAND (j - 1, j))
		{
		  ld5 = P->dangle5[type][S1[j - 1]];
		  if ((p = pair_table[j - 2]) && SAME_STRAND (j - 2, j - 1))
		    if (P->dangle3[pair[S[p]][S[j - 2]]][S1[j - 1]] < ld5)
		      ld5 = 0;
		}
	    }
	}
      i1 = i;
      p = i + 1;
      u = 0;
      energy = 0;
      cx_energy = INF;
      do
	{			/* walk around the multi-loop */
	  int tt, new_cx = INF;

	  /* hope over unpaired positions */
	  while (p <= pair_table[0] && pair_table[p] == 0)
	    p++;

	  /* memorize number of unpaired positions */
	  u += p - i1 - 1;
	  /* get position of pairing partner */
	  if (p == pair_table[0] + 1)
	    q = tt = 0;		/* virtual root pair */
	  else
	    {
	      q = pair_table[p];
	      /* get type of base pair P->q */
	      tt = pair[S[p]][S[q]];
	      if (tt == 0)
		tt = 7;
	    }

	  energy += mlintern[tt];
	  cx_energy += mlintern[tt];

	  if (dangles)
	    {
	      int dang5 = 0, dang3 = 0, dang;
	      if ((SAME_STRAND (p - 1, p)) && (p > 1))
		dang5 = P->dangle5[tt][S1[p - 1]];	/* 5'dangle of pq pair */
	      if ((SAME_STRAND (i1, i1 + 1)) && (i1 < S[0]))
		dang3 = P->dangle3[type][S1[i1 + 1]];	/* 3'dangle of previous pair */

	      switch (p - i1 - 1)
		{
		case 0:	/* adjacent helices */
		  if (dangles == 2)
		    energy += dang3 + dang5;
		  else if (dangles == 3 && i1 != 0)
		    {
		      if (SAME_STRAND (i1, p))
			{
			  new_cx = energy + P->stack[rtype[type]][rtype[tt]];
			  /* subtract 5'dangle and TerminalAU penalty */
			  new_cx +=
			    -ld5 - mlintern[tt] - mlintern[type] +
			    2 * mlintern[1];
			}
		      ld5 = 0;
		      energy = MIN2 (energy, cx_energy);
		    }
		  break;
		case 1:	/* 1 unpaired base between helices */
		  dang =
		    (dangles == 2) ? (dang3 + dang5) : MIN2 (dang3, dang5);
		  if (dangles == 3)
		    {
		      energy = energy + dang;
		      ld5 = dang - dang3;
		      /* may be problem here: Suppose
		         cx_energy>energy, cx_energy+dang5<energy
		         and the following helices are also stacked (i.e.
		         we'll subtract the dang5 again */
		      if (cx_energy + dang5 < energy)
			{
			  energy = cx_energy + dang5;
			  ld5 = dang5;
			}
		      new_cx = INF;	/* no coax stacking with mismatch for now */
		    }
		  else
		    energy += dang;
		  break;
		default:	/* many unpaired base between helices */
		  energy += dang5 + dang3;
		  if (dangles == 3)
		    {
		      energy = MIN2 (energy, cx_energy + dang5);
		      new_cx = INF;	/* no coax stacking possible */
		      ld5 = dang5;
		    }
		}
	      type = tt;
	    }
	  if (dangles == 3)
	    cx_energy = new_cx;
	  i1 = q;
	  p = q + 1;
	}
      while (q != i);
      best_energy = MIN2 (energy, best_energy);	/* don't use cx_energy here */
      /* fprintf(stderr, "%6.2d\t", energy); */
      if (dangles != 3 || is_extloop)
	break;			/* may break cofold with co-ax */
      /* skip a helix and start again */
      while (pair_table[p] == 0)
	p++;
      if (i == pair_table[p])
	break;
      i = pair_table[p];
    }
  energy = best_energy;
  energy += mlclosing;
  /* logarithmic ML loop energy if logML */
  if ((!is_extloop) && logML && (u > 6))
    energy += 6 * mlbase + (int) (P->lxc * log ((double) u / 6.));
  else
    energy += mlbase * u;
  /* fprintf(stderr, "\n"); */
  return energy;
}

/*---------------------------------------------------------------------------*/

PRIVATE int
cut_in_loop (int i)
{
  /* walk around the loop;  return j pos of pair after cut if
     cut_point in loop else 0 */
  int p, j;
  p = j = pair_table[i];
  do
    {
      i = pair_table[p];
      p = i + 1;
      while (pair_table[p] == 0)
	p++;
    }
  while (p != j && SAME_STRAND (i, p));
  return SAME_STRAND (i, p) ? 0 : pair_table[p];
}

/*---------------------------------------------------------------------------*/

PRIVATE void
make_ptypes (const short *S, const char *structure)
{
  int n, i, j, k, l;

  n = S[0];
  for (k = 1; k < n - TURN; k++)
    for (l = 1; l <= 2; l++)
      {
	int type, ntype = 0, otype = 0;
	i = k;
	j = i + TURN + l;
	if (j > n)
	  continue;
	type = pair[S[i]][S[j]];
	while ((i >= 1) && (j <= n))
	  {
	    if ((i > 1) && (j < n))
	      ntype = pair[S[i - 1]][S[j + 1]];
	    if (noLonelyPairs && (!otype) && (!ntype))
	      type = 0;		/* i.j can only form isolated pairs */
	    ptype[indx[j] + i] = (char) type;
	    otype = type;
	    type = ntype;
	    i--;
	    j++;
	  }
      }

  if (fold_constrained && (structure != NULL))
    {
      int hx, *stack;
      char type;
      stack = (int *) space (sizeof (int) * (n + 1));

      for (hx = 0, j = 1; j <= n; j++)
	{
	  switch (structure[j - 1])
	    {
	    case '|':
	      BP[j] = -1;
	      break;
	    case 'x':		/* can't pair */
	      for (l = 1; l < j - TURN; l++)
		ptype[indx[j] + l] = 0;
	      for (l = j + TURN + 1; l <= n; l++)
		ptype[indx[l] + j] = 0;
	      break;
	    case '(':
	      stack[hx++] = j;
	      /* fallthrough */
	    case '<':		/* pairs upstream */
	      for (l = 1; l < j - TURN; l++)
		ptype[indx[j] + l] = 0;
	      break;
	    case ')':
	      if (hx <= 0)
		{
		  fprintf (stderr, "%s\n", structure);
		  nrerror ("unbalanced brackets in constraints");
		}
	      i = stack[--hx];
	      type = ptype[indx[j] + i];
	      for (k = i + 1; k <= n; k++)
		ptype[indx[k] + i] = 0;
	      /* don't allow pairs i<k<j<l */
	      for (l = j; l <= n; l++)
		for (k = i + 1; k <= j; k++)
		  ptype[indx[l] + k] = 0;
	      /* don't allow pairs k<i<l<j */
	      for (l = i; l <= j; l++)
		for (k = 1; k <= i; k++)
		  ptype[indx[l] + k] = 0;
	      for (k = 1; k < j; k++)
		ptype[indx[j] + k] = 0;
	      ptype[indx[j] + i] = (type == 0) ? 7 : type;
	      /* fallthrough */
	    case '>':		/* pairs downstream */
	      for (l = j + TURN + 1; l <= n; l++)
		ptype[indx[l] + j] = 0;
	      break;
	    }
	}
      if (hx != 0)
	{
	  fprintf (stderr, "%s\n", structure);
	  nrerror ("unbalanced brackets in constraint string");
	}
      free (stack);
    }
}

/* Last changed Time-stamp: <2005-09-22 10:10:55 ivo> */
/*                
		  Ineractive Access to folding Routines

		  c Ivo L Hofacker
		  Vienna RNA package
*/

#include <stdlib.h>
#include "part_func.h"
#include "PS_dot.h"
extern void read_parameter_file (const char fname[]);
/*@unused@*/
/* static char rcsid[] = "$Id: RNAfold.c,v 1.18 2005/10/30 20:07:12 ivo Exp $";
 */

#define PRIVATE static

struct Hairpin
{
  char *name;			/* nom de l'hairpin */
  char *seq;			/* sequence de l'hairpin */
  char *shape;			/* forme de l'hairpin en point, parentheses */
  float mfe;			/* mfe de l'hairpin */
  float mfei;			/* mfei de l'hairpin */
  int i;			/* nb de nucleotide avant le mirna putatif */
  int j;			/* nb de nucleotide apres le mirna putatif */
};

/*--------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  char *string, *name;
  char *structure = NULL;
  char fname[40];
  int i, j, length, l, stems, minlen, nt_before, nt_after;
  double min_en;
  int noconv = 0;
  FILE *output;
  double max_amfe, max_mfei, max_unpaired, min_rapp, max_rapp;

  if (argc != 13)
    {
      printf
	("\nUse %s <sequence> <nb_nt_before_mirna> <nb_nt_after_mirna> <max amfe> <max mfei> <max unpaired> <ratiomin> <ratiomax> <output_file> <stems>  <min pre-Mi len> <name>\n\n",
	 argv[0]);
      return EXIT_FAILURE;
    }
  string = argv[1];
  nt_before = atoi (argv[2]);
  nt_after = atoi (argv[3]);
  max_amfe = atof (argv[4]);
  max_mfei = atof (argv[5]);
  max_unpaired = atof (argv[6]);
  min_rapp = atof(argv[7]);
  max_rapp = atof(argv[8]);
  (void) sscanf (argv[9], "%39s", fname);
  stems = atoi(argv[10]);
  minlen = atoi(argv[11]);
  name = argv[12];
  output = fopen (fname, "a");

  if (output == NULL)
    {
      printf ("Cannot open %s, exit\n", fname);
      return EXIT_FAILURE;
    }


  do_backtrack = 1;
  length = (int) strlen (string);

  structure = (char *) space ((unsigned) length + 1);
  for (l = 0; l < length; l++)
    {
      string[l] = toupper (string[l]);
      if (!noconv && string[l] == 'T')
	string[l] = 'U';
    }

  fold (string, structure);
  int end_mirna = length - nt_after;	/* la numerotation commence a 1 */
  struct Hairpin hairpin;
  hairpin.mfei = 0.;
  hairpin.name = NULL;
  hairpin.seq = NULL;
  hairpin.shape = NULL;
  hairpin.mfe = 0.;
  hairpin.i = 0;
  hairpin.j = 0;
  int boolean = 0;

  for (i = nt_before + 1; i >= 1; i--)
    {
      for (j = end_mirna; j <= length; j++)
	{
	  if (j - i + 1 >= minlen)
	    {
	      min_en = (backtrack (string, 0, i, j)) / 100.;
	      parenthesis_structure (structure, length);
	      char *substring = (char *) malloc (sizeof (char) * (j - i + 2));
	      substring = strncpy (substring, string + i - 1, j - i + 1);
	      substring[j - i + 1] = '\0';
	      char *substruct = (char *) malloc (sizeof (char) * (j - i + 2));
	      substruct = strncpy (substruct, structure + i - 1, j - i + 1);
	      substruct[j - i + 1] = '\0';
	      int z, nbGC = 0;
	      for (z = 0; z < j - i + 1; z++)
		{
		  if (substring[z] == 'g' || substring[z] == 'G'
		      || substring[z] == 'c' || substring[z] == 'C')
		    nbGC++;
		}
	      double propGC = (double) nbGC / (j - i + 1);
	      propGC *= 100;

	      char *mirna_ss =
		(char *) malloc (sizeof (char) *
				 (length - nt_before - nt_after + 1));
	      mirna_ss =
		strncpy (mirna_ss, structure + nt_before,
			 length - nt_before - nt_after);
	      mirna_ss[length - nt_before - nt_after] = '\0';
	      int nb_point = 0;
	      for (z = 0; z < length - nt_before - nt_after; z++)
		{
		  if (mirna_ss[z] == '.')
		    nb_point++;
		}

	      int b = base_pair[0].i;
	      int k, min = length, max = 0;
	      for (k = 1; k <= b; k++)
		{
		  int bpi = base_pair[k].i;
		  int bpj = base_pair[k].j;
		  if (bpi > nt_before && bpi <= end_mirna)
		    {
		      if (bpj < min)
			min = bpj;
		      if (bpj > max)
			max = bpj;
		    }
		  else if (bpj > nt_before && bpj <= end_mirna)
		    {
		      if (bpi < min)
			min = bpi;
		      if (bpj > max)
			max = bpi;
		    }
		}
              
              int nbStem = 0, closing = 0;
              for(k=0; k < j - i + 1; k++){
                if(substruct[k] == ')' && !closing){
                  closing = 1;
                  nbStem++;
                }else{
                  if(substruct[k] == '(')
                    closing = 0;
                }
              }
              
	      double pl =
		(double) nb_point / (double) (length - nt_before - nt_after);
	      double rapp = (end_mirna - nt_before) / (abs (max - min) + 1.0);
	      double mfei = ((min_en / (j - i + 1)) * 100) / propGC;
	      if ((nbStem <= stems || stems == 0) && ((min_en / (j - i + 1)) * 100 < max_amfe) &&
		  (pl <= max_unpaired)
		  && (!(strchr (mirna_ss, '(') != NULL
			&& strchr (mirna_ss, ')') != NULL))
		  && mfei < max_mfei && rapp >= min_rapp && rapp <= max_rapp)
		{		/* une nouvelle hairpin valide les conditions */
		  if (!boolean || (hairpin.mfei > mfei))
		    {		/* premiere hairpin ou hairpin avec un meilleur mfei */
		      hairpin.name = name;
		      hairpin.seq =
			(char *) malloc ((strlen (substring) + 1) *
					 (sizeof (char)));
		      strcpy (hairpin.seq, substring);
		      hairpin.shape =
			(char *) malloc ((strlen (substruct) + 1) *
					 (sizeof (char)));
		      strcpy (hairpin.shape, substruct);
		      hairpin.mfe = min_en;
		      hairpin.mfei = mfei;
		      hairpin.i = i;
		      hairpin.j = j;
		      boolean = 1;	/* il y a maintenant au moins une hairpin valide */
		    }
		}
              free(substring);
              free(substruct);
              free(mirna_ss);
	    }
	}
    }				/* fin de la double boucle parcourrant toutes les sous-sequences */

  if (boolean)
    {				/* il y a au moins une hairpin validant les criteres */
/*                     printf ("i= %d, j= %d, end_mir= %d\n", hairpin.i, hairpin.j, end_mirna);
 */
      fprintf (output, ">%s_%d_%d\n", hairpin.name, nt_before - hairpin.i + 1,
	       hairpin.j - end_mirna);
      fprintf (output, "%s\n%s (%6.2f)\n", hairpin.seq,
	       hairpin.shape, hairpin.mfe);
      free(hairpin.seq);
      free(hairpin.shape);
/* 		    fprintf (output, "%s\n%s\n", mirna, mirna_ss);
 * 		    if (fname[0] != '\0')
 * 		    {
 * 		    strcpy (ffname, fname);
 * 		    strcat (ffname, "_ss.ps");
 * 
 * nopostscript
 * 		    int begin = nt_before - hairpin.i + 1;
 * 		    int after = hairpin.j - end_mirna;
 * 		    sprintf(ffname, "%s_%d_%d_ss.ps", hairpin.name, begin, after);
 * 
 * 		    strcpy (gfname, fname);
 * 		    strcat (gfname, "_ss.g");
 * 		    }
 * 
 * nopostscript
 *                     if (length < 2000) {
 * 			begin += 1;
 * 			int end = end_mirna - hairpin.i + 1;
 * 			//printf ("ffname= %s\n", ffname);
 * 			(void) PS_rna_plot(hairpin.seq, hairpin.shape, ffname,
 * 					   ffname, begin, end);
 * 		    } else {
 * 			struct bond *bp;
 * 			fprintf(stderr,
 * 				"INFO: structure too long, not doing xy_plot\n");
 * 
 * 			// free mfe arrays but preserve base_pair for PS_dot_plot 
 * 			bp = base_pair;
 * 			base_pair = space(16);
 * 			free_arrays();	// free's base_pair 
 * 			base_pair = bp;
 * 		    }
  */
      return EXIT_SUCCESS;
    }

  if (length >= 2000)
    free (base_pair);
  (void) fflush (stdout);
  free (structure);
  fflush(output);
  fclose(output);
  return EXIT_SUCCESS;
}
