/*
 * Made by G. Myers (1999)
 * Modified by A. Mathelier
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/file.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "parse.i"
#define WORD long
#define BUF_MAX 2048				/** taille maximale du buffer pour la lecture des lettres du texte */

static int W;					/** taille d'un mot memoire */
static unsigned WORD All = -1;			/** tous les bits sont mis a 1 */
static unsigned WORD Ebit;			/** masque pour le bloc representant la derniere lettre du motif */
static unsigned WORD *Pc[5];
static unsigned WORD *first[5];
static int seg;					/** nombre de mot(s) memoire necessaire(s) */
static int rem;					/** nombre de blocs non codant dans le dernier mot memoire */
static unsigned WORD Un = 0L;
static unsigned WORD Ov = 0L, Z, ZOv;

void
setup_search (int bsize, char *pat)
{
  register int a, i;
  W = sizeof (unsigned WORD) * 8;
  int blocw = W / bsize;			/** nombre de bloc maximal dans un mot memoire */
  if (blocw > patlen)
    {
      seg = 1;
      rem = patlen;
    }
  else
    {
      seg = patlen / blocw + 1;
      rem = patlen % blocw;
    }
  if (rem == 0)
    {
      rem = blocw;
      seg -= 1;
    }
  rem = W / bsize - rem;

  for (i = 0; i < 5; i++)
    {
      first[i] = (unsigned WORD *) malloc (seg * sizeof (unsigned WORD));
      Pc[i] = (unsigned WORD *) malloc ((seg+1) * sizeof (unsigned WORD));
      for (a = 0; a < seg; a++)
	{
	  Pc[i][a] = -1;
	  first[i][a] = 0L;
	}
       Pc[i][seg] = -1;
    }
  int boola = 1, boolc = 1, boolg = 1, boolt = 1;
  for (a = 0; a < patlen; a++)
    {
      int mw = a / blocw;
      int decalage = a % blocw;
      register unsigned WORD mask =
	((((1L << bsize) - 1L) << (decalage * bsize)));
      switch (pat[a])
	{
	case 'a':
	case 'A':
	  if (boola)
	    {
	      boola = 0;
	      first[0][mw] = mask;
	    }
	  Pc[1][mw] -= mask;
	  Pc[2][mw] -= mask;
	  Pc[3][mw] -= mask;
	  Pc[4][mw] -= mask;
	  break;
	case 'c':
	case 'C':
	  if (boolc)
	    {
	      boolc = 0;
	      first[1][mw] = mask;
	    }
	  Pc[0][mw] -= mask;
	  Pc[2][mw] -= mask;
	  Pc[3][mw] -= mask;
	  Pc[4][mw] -= mask;
	  break;
	case 'g':
	case 'G':
	  if (boolg)
	    {
	      boolg = 0;
	      first[2][mw] = mask;
	    }
	  Pc[0][mw] -= mask;
	  Pc[1][mw] -= mask;
	  Pc[3][mw] -= mask;
	  Pc[4][mw] -= mask;
	  break;
        case 'u':
        case 'U':
	case 't':
	case 'T':
	  if (boolt)
	    {
	      boolt = 0;
	      first[3][mw] = mask;
	    }
	  Pc[0][mw] -= mask;
	  Pc[1][mw] -= mask;
	  Pc[2][mw] -= mask;
	  Pc[4][mw] -= mask;
	  break;
	}
    }

  Ebit = ((1L << bsize) - 1L) << ((blocw - 1) * bsize);
  for (i = 0; i < blocw; i++)
    {
      Un += 1L << (i * bsize);
      Ov += (1L << (bsize - 1)) << (i * bsize);
    }
  Z = Ov - Un;
  ZOv = Z | Ov;
}

typedef struct
{
  unsigned WORD P;				/** correspond a Pv */
  unsigned WORD M;				/** correspond a Mv */
  unsigned WORD Ph;
  unsigned WORD Mh;
  unsigned WORD L;				/** correspond a la matrice des longueurs */
  int V;					/** correspond au score de la cellule pour la matrice des erreurs */
} Scell;

void
search (ifile, dif, bsize, patlen, ml)
     int ifile, dif, bsize, patlen, ml;
{
  int num, i, base, diw, a, Cscore, blocw;
  int* diff = (int *) malloc (12 * sizeof (int));
  for (i = 0; i < 12; i++)
    diff[i] = 0;
  register Scell *s, *sd;
  unsigned WORD pc, mc;
  register unsigned WORD *e;
  register unsigned WORD P, M, U = 0L, X, Y, Mh, Ph,
           mask, lbloc = 0L, Md, Mf, fm = -1;
  register Scell *S, *SE;
  static char buf[BUF_MAX];

  S = (Scell *) malloc (sizeof (Scell) * seg);
  unsigned WORD bloc = (1L << bsize) - 1L;	/** que des 1 dans le premier bloc */
  blocw = W / bsize;				/** nombre de bloc maximal dans un mot memoire */
  int dec = (blocw * bsize) - bsize;            /* decalage utilisé */
  SE = S + (seg - 1);				/** pointe sur la derniere cellule, i.e. le dernier mot memoire */
  diw = dif + blocw;				/** nombre d'erreurs autorisee + nombre de bloc dans un mot memoire */

  mask = bloc << ((blocw - rem - 1) * bsize);

  for (i = 0; i < rem; i++)
    lbloc += mask << ((i + 1) * bsize);

  sd = S + (dif - 1) / blocw;			/** pointe sur la cellule qui contient le nombre maximal d'erreurs autorisees */
  if (sd == SE - 1)
    sd++;
  for (s = S; s <= sd; s++)
    {
      s->P = All;				/** on place tous les bits a 1 */
      s->M = 0;					/** on place tous les bits a 0 */
      s->Ph = 0;
      s->Mh = 0;
      s->L = 0;					/** toutes les longueurs sont initialisees a 0 */
      s->V = ((s - S) + 1) * blocw;		/** equivalent a b*w */
    }

  Mf = 0L;
  fm -= bloc;
  register unsigned WORD *Mf_old, *Lold, *M1old;
  Mf_old = (unsigned WORD *) malloc (seg * sizeof (unsigned WORD));
  Lold = (unsigned WORD *) malloc (seg * sizeof (unsigned WORD));
  M1old = (unsigned WORD *) malloc (seg * sizeof (unsigned WORD));

  int c;
  for (c = 0; c < seg; c++){
    Mf_old[c] = 0L;
    M1old[c] = 0L;
    Lold[c] = Un;
  }

  for (base = 1 - rem; (num = read (ifile, buf, BUF_MAX)) > 0; base += num)
    {                                           /** on va parcourir toutes les lettres du texte */
      i = 0;
      if (sd == S)
	{					/** le nombre d'erreurs max se trouve sur le premier mot memoire */
	  P = S->P;				/** P = Pv */
	  M = S->M;				/** M = Mv */
	  Cscore = S->V;			/** score de la cellule (a la fin de la cellule */
	  for (; i < num; i++)
	    {					/** on parcourt les dernieres lettres lues dans le texte */
	      a = buf[i];			/** lettre lue */
	      int letter;
	      switch (a)
		{
		case 'a':
		case 'A':
		  letter = 0;
		  break;
		case 'c':
		case 'C':
		  letter = 1;
		  break;
		case 'g':
		case 'G':
		  letter = 2;
		  break;
                case 'u':
                case 'U':
		case 't':
		case 'T':
		  letter = 3;
		  break;
                case '\n':
                  letter = 4;
                  base--;
                  break;
		default:
                  letter = 4;
		}
              
              if (letter == 4)
                continue;

	      U = Pc[letter][0];		/** masque pour la position de la lettre lue dans le premier mot memoire */
	      X = ((U & P) + P) & Un;
	      X = (((X << bsize) - X) ^ P) | U;	/** X = Xh */
	      U |= M;				/** U devient Xv */
	      Y = P;				/** Y = Pv */
	      P = M | ~(X | Y);			/** P devient Ph */
	      Ph = P & ZOv;
	      M = Y & X;			/** M devient Mh */
	      Mh = M & ZOv;
	      Md = (~(P | M)) << bsize;
	      if (P & Ebit)
		Cscore += 1;
	      else if (M & Ebit)
		Cscore -= 1;

	      Y = P << bsize;			/** Y devient Ph */
	      P = (M << bsize) | ~(U | Y);	/** P devient Pv */
	      M = Y & U;			/** M devient Mv */
              Mf = Mf_old[0];
	      register unsigned WORD Ltamp, M1, M2, M3, L;
	      M1 = (P | (Pc[letter][0] & first[letter][0] & Md)) & fm;
	      M3 = (~P & ~Pc[letter][0] & ~Mf) & fm;
	      M2 = (~(M1 | M3)) & fm;
              
/*               Ltamp = (M1old[0] << bsize) & M1;
 *               L = Ltamp | M3;
 *               M3 = ((((L + (Un & M3)) ^ L) & L) | M3) & fm;
 *               L = Ltamp | M2;
 *               M2 = ((((L + (Un & M2)) ^ L) & L) | M2) & fm;
 *               M1old[0] = M1;
 *               M1 = (~(M2 | M3)) & fm;
 */

	      L = ((S->L << bsize) + Un) & M2;
	      L += (S->L + Un) & M3;
	      L += Un & ~fm;
	      Ltamp = L & ~M1;
              diff[0] += 3;

	      while (M1 != 0L)
		{
		  L = (L << bsize) & M1;
		  M1 = M1 & (M1 << bsize);
		  Ltamp += L & ~M1;
		  L = Ltamp;
		}
	      Ltamp -= ((Ltamp & ~Z) >> (bsize - 1)) & Un;
	      if (S == SE)
		{
		  S->L <<= bsize;
		  S->L = (S->L & lbloc) + (Ltamp & ~lbloc);
		  S->L &= ZOv;
		}
	      else
		S->L = Ltamp & ZOv;

	      Mf = ~M;
              Mf_old[0] = Mf;

	      if (Cscore <= dif)
		break;
	    }					/** fin de for (;i < num; i++) */

	  S->P = P;				/** on attribue les masques de la premiere cellule */
	  S->M = M;				/** idem */
	  S->V = Cscore;			/** idem */

	  if (i >= num)
	    continue;				/** le break n'a pas ete execute donc on lit les lettres suivantes s'il en reste */
	  {
	    /* les 1e et dernier mot mem sont egaux donc nous sommes sur un seul mot mem */
             int length = (S->L & Ebit) >> (blocw - 1) * bsize;
             if ((sd == SE) && (length <= ml))
               printf ("  Match at %d with length %d\n", base + i, length);
	  }
	  i += 1;
	}					/** fin du if (sd == S) */

      for (; i < num; i++)
	{					/** le nombre max d'erreurs ne se trouve pas dans le premier mot mem */
	  a = buf[i];
	  int letter;
	  switch (a)
	    {
	    case 'a':
	    case 'A':
	      letter = 0;
	      break;
	    case 'c':
	    case 'C':
	      letter = 1;
	      break;
	    case 'g':
	    case 'G':
	      letter = 2;
	      break;
            case 'u':
            case 'U':
	    case 't':
	    case 'T':
	      letter = 3;
	      break;
            case '\n':
              base--;
              letter = 4;
              break;
	    default:
	      /*fprintf (stderr, "letter %c not in [a,c,g,t]\n", a);*/
              /*exit (EXIT_FAILURE);*/
              letter = 4;
	    }
          if (letter == 4)
            continue;

	  e = Pc[letter];			/** e represente les mots mem decrivant les pos de la lettre lue dans le motif */
	  pc = mc = 0;				/** P et M initialises a 0, mc correspond au Mv de l'algo de base qui est initialise a 0 */
	  s = S;				/** s est place sur le premier mot mem */
	  int k = 0;
	  while (s <= sd)
	    {					/** on avance sur les mots mem jusqu'a celui qui contient le nombre max d'erreurs */
	      U = *e++;				/** U contient les pos de la lettre et on positionne e sur le prochain mot mem, U = Eq */
	      P = s->P;				/** on prend l'ancienne valeur de Pv pour le mot mem */
	      M = s->M;				/** idem avec Mv */
	      Y = U | mc;			/** correspond au between l.7 and 8 of basic algo : Y = (Eq |= 1) if hin < 0 */
	      X = ((Y & P) + P) & Un;
	      X = (((X << bsize) - X) ^ P) | Y;	/** X = Xh */
	      U |= M;				/** U = Xv */
	      Y = P;				/** Y = Pv */
	      P = M | ~(X | Y);			/** P = Ph */
	      Ph = P & ZOv;
	      M = Y & X;			/** M = Mh */
	      Mh = M & ZOv;
	      Y = (P << bsize) | pc;		/** Y = Ph, correspondant a l'ajout between l.16 and 17 */
	      s->P = (M << bsize) | mc | ~(U | Y);/** s->P prend Pv */
	      s->M = Y & U;			/** s->M prend Mv */
	      U = s->V;				/** U prend la valeur du score = hout */
	      pc = mc = 0;
	  /** les valeurs des scores sont mises a jour et
	   *  pc et mc egalement pour les lignes ajoutees au basic algo
	   */
	      if (P & Ebit)
		{
		  pc = bloc;
		  s->V = U + 1;
		}
	      else if (M & Ebit)
		{
		  mc = bloc;
		  s->V = U - 1;
		}

	      Mf = Mf_old[k];
	      if (k == 0)
		Md = (~(P | M)) << bsize;
	      else
                Md =
                  ((~(Ph | Mh)) << bsize) +
                  (((~((s - 1)->Ph | (s - 1)->Mh)) & ZOv) >> dec);
	      register unsigned WORD Ltamp, M1, M2, M3, L;
              int c;
	      M1 = s->P | (Pc[letter][k] & first[letter][k] & Md);
	      M3 = ~(s->P) & ~Pc[letter][k] & ~Mf;
	      M2 = ~(M1 | M3);

	      if (k == 0)
		{
		  M1 &= fm;
/*                   M2 &= fm;
 *                   M3 &= fm;
 *                   Ltamp = (M1old[k] << bsize) & M1;
 *                   L = Ltamp | M3;
 *                   M3 = ((((L + (Un & M3)) ^ L) & L) | M3) & fm;
 *                   L = Ltamp | M2;
 *                   M2 = ((((L + (Un & M2)) ^ L) & L) | M2) & fm;
 *                   M1old[k] = M1;
 *                   M1 = (~(M2 | M3)) & fm;
 */

		  L = ((s->L << bsize) + Un) & M2 & fm;
		  L += (s->L + Un) & M3 & fm;
		  L += Un & ~fm;
                  diff[3] += 1;
		}
	      else
		{
/* 		  Ltamp = (M1old[k] << bsize) & M1;
 *                   L = Ltamp | M3;
 *                   M3 = ((((L + (Un & M3)) ^ L) & L) | M3);
 *                   L = Ltamp | M2;
 *                   M2 = ((((L + (Un & M2)) ^ L) & L) | M2);
 *                   M1old[k] = M1;
 *                   M1 = (~(M2 | M3));
 */
                  register unsigned WORD Ljdecal =
		    (s->L << bsize) + (Lold[k - 1] >> dec);
		  L = ((Ljdecal) + Un) & M2;
		  L += (s->L + Un) & M3;
		}
	      Ltamp = L & ~M1;
              diff[5] += 1;
	      c = 0;

	      while (M1 != 0L)
		{
		  if (c == 0 && s > S)
		    L =
		      ((L << bsize) +
		       ((s - 1)->L >> ((blocw - 1) * bsize))) & M1;
		  else
		    L = (L << bsize) & M1;
		  M1 = M1 & (M1 << bsize);
		  Ltamp += L & ~M1;
		  L = Ltamp;
		  c++;
		}
	      Ltamp -= ((Ltamp & ~Z) >> (bsize - 1)) & Un;
	      Lold[k] = s->L;
	      if (s == SE)
		{
		  s->L <<= bsize;
		  s->L = (s->L & lbloc) + (Ltamp & ~lbloc);
		  s->L &= ZOv;
		}
	      else
		s->L = Ltamp & ZOv;

	      Mf_old[k] = ~(s->M);
	      s->Mh = Mh;
	      s->Ph = Ph;

	      s++;
	      k++;
	    }					/** fin du while (s <= sd) */

	  if ((U <= dif) && ((*e & bloc) | mc) && (s <= SE))
	    {					/** l.10 - block base algo */
	      s->P = All;
	      s->M = 0;
	      if (pc == bloc)
		s->M = bloc;
	      if (mc != bloc)
		s->P <<= bsize;
	      s->V = U = U + blocw - 1;
	      sd = s;
	      Ph = pc & ZOv;
	      Mh = mc & ZOv;
	      Mf = Mf_old[k];

	      if (k == 0)
		Md = (~(pc | mc)) << bsize;
	      else
		Md =
		  ((~(pc | mc)) << bsize) +
		  (((~((s - 1)->Ph | (s - 1)->Mh)) & ZOv) >> dec);

	      register unsigned WORD Ltamp, M1, M2, M3, L;
	      M1 = s->P | (Pc[letter][k] & first[letter][k] & Md);
	      M3 = ~(s->P) & ~Pc[letter][k] & ~Mf;
	      M2 = ~(M1 | M3);
              int c;
	      if (k == 0)
		{
		  M1 &= fm;
/*                   M2 &= fm;
 *                   M3 &= fm;              
 *                   Ltamp = (M1old[k] << bsize) & M1;
 *                   L = Ltamp | M3;
 *                   M3 = ((((L + (Un & M3)) ^ L) & L) | M3) & fm;
 *                   L = Ltamp | M2;
 *                   M2 = ((((L + (Un & M2)) ^ L) & L) | M2) & fm;
 *                   M1old[k] = M1;
 *                   M1 = (~(M2 | M3)) & fm;
 */

		  L = ((s->L << bsize) + Un) & M2 & fm;
		  L += (s->L + Un) & M3 & fm;
		  L += Un & ~fm;
		}
	      else
		{
/*                   Ltamp = (M1old[k] << bsize) & M1;
 *                   L = Ltamp | M3;
 *                   M3 = ((((L + (Un & M3)) ^ L) & L) | M3);
 *                   L = Ltamp | M2;
 *                   M2 = ((((L + (Un & M2)) ^ L) & L) | M2);
 *                   M1old[k] = M1;
 *                   M1 = (~(M2 | M3));
 */

		  register unsigned WORD Ljdecal =
		    (s->L << bsize) + (Lold[k - 1] >> dec);
		  L = ((Ljdecal) + Un) & M2;
		  L += (s->L + Un) & M3;
		}
	      Ltamp = L & ~M1;
	      c = 0;

	      while (M1 != 0L)
		{
		  if (c == 0 && s > S)
		    L =
		      ((L << bsize) +
		       ((s - 1)->L >> ((blocw - 1) * bsize))) & M1;
		  else
		    L = (L << bsize) & M1;
		  M1 = M1 & (M1 << bsize);
		  Ltamp += L & ~M1;
		  L = Ltamp;
		  c++;
		}
	      Ltamp -= ((Ltamp & ~Z) >> (bsize - 1)) & Un;
	      Lold[k] = s->L;
	      if (s == SE)
		{
		  s->L <<= bsize;
		  s->L = (s->L & lbloc) + (Ltamp & ~lbloc);
		  s->L &= ZOv;
		}
	      else
		s->L = Ltamp & ZOv;
	      Mf_old[k] = ~(s->M);
	      s->Mh = Mh;
	      s->Ph = Ph;
	    }
	  else
	    {
	      U = sd->V;
	      while (U > diw)
		{
		  U = (--sd)->V;
		}
	    }
	  {
             int length = (SE->L & Ebit) >> (blocw - 1) * bsize;
             if ((sd == SE) && (U <= dif) && (length <= ml))
 	    printf ("  Match at %d with length %d\n", base + i, length);
	  }
	}					/** fin du for (;i < num; i++) */

      while (sd > S && sd < SE)
	{
	  i = sd->V;
          if (i <= dif)
            break;
	  P = sd->P;
	  M = sd->M;
	  Y = Ebit;
	  for (X = 0; X < blocw; X++)
	    {
	      if (P & Y)
		{
		  i -= 1;
		  if (i <= dif)
		    break;
		}
	      else if (M & Y)
		{
		  i += 1;
		  if (i <= dif)
		    break;			/** ajout */
		}
	      Y >>= bsize;
	    }
	  if (i <= dif)
	    break;
	  sd -= 1;
	}					/** fin du while (sd > S) */
    }						/** fin du for (base = 1 - rem; ...) */

  if (sd == SE)
    {
      P = sd->P;
      M = sd->M;
      U = sd->V;
      SE->L <<= bsize;
      for (i = 0; i < rem; i++)
	{
	  if (P & Ebit)
	    U -= 1;
	  else if (M & Ebit)
	    U += 1;
	  P <<= bsize;
	  M <<= bsize;
	  {
             int length = (SE->L & Ebit) >> ((blocw - 1) * bsize);
             if ((U <= dif) && (length <= ml))
               printf ("  Match at %d with length %d\n", base + i, length);
	  }
	  SE->L <<= bsize;
	}					/** fin du for (i = 0; i < rem; i++) */
    }						/** fin du if (sd == SE) */
  free (Mf_old);
  free (M1old);
  free (Lold);
  free (diff);
  free(S);
  for(i = 0; i < 5; i++){
    free(first[i]);
    free(Pc[i]);
  }
}						/** fin de la procedure search */

#include "main.i"
