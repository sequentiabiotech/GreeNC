/*
 * Made by G. Myers (1999)
 * Modified by A. Mathelier
 */

static char *progname;

void error(char* msg)
{ fprintf(stderr,"%s: ",progname);
  fprintf(stderr,msg, 0);
  fprintf(stderr,"\n");
}

#define SIG  70
#define BUF 256

void   srand48();
double drand48();

int main(argc,argv) int argc; char *argv[];
{ int   ifile;
  int   dif;
  char *pat;
  int   ml;
  int   len, alpha;
  char *msg = malloc (BUF * sizeof(char));
  static char *cset = 
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$%^&*()";

  progname = argv[0];

  if (argc < 5){
    sprintf(msg, "Usage is '%s <pat> <dif> <file> <maxlength>'",progname);
    error(msg);
    free(msg);
    exit(1);
  }

  if (isdigit(argv[1][0])){
    if (sscanf(argv[1],"%d,%d",&len,&alpha) != 2){
      error("pattern spec should be of the form #l,#a");
      free(msg);
      exit(1);
    }
    if (len < 0 || alpha < 0){
      error("only positive arguments in pattern spec.");
      free(msg);
      exit(1);
    }
    if (alpha > SIG){
      sprintf (msg, "alphabet has more than %d symbols",SIG);
      error(msg);
      free(msg);
      exit(1);
    }
    pat = (char *) malloc(len+1);
    srand48(SIG*len + alpha);
    pat[len] = '\0';
    while (len-- > 0)
      pat[len] = cset[(int) (drand48()*alpha)];
  }
  else
    pat = argv[1];

  if (argc >= 3)
    { dif = atoi(argv[2]);
      if (dif < 0){
        sprintf(msg, "Threshold, %d, is negative\n",dif);
        error(msg);
        free(msg);
        exit(1);
      }
    }
  else
    dif = 0;

  if (argc >= 4)
    { ifile = open(argv[3],O_RDONLY);
      if (ifile == -1){
        sprintf(msg, "Can't open file %s",argv[3]);
        error(msg);
        free(msg);
        exit(1);
      }
    }
  else
    ifile = 1;

  if (argc >= 5){
      ml = atoi (argv[4]);
      if (ml < 1){
        sprintf(msg, "maxlength, %d, must be greater than 0\n", ml);
        error(msg);
        free(msg);
        exit(1);
      }
  }
  else
      ml = 1;

  encode_pattern(pat);
  if (dif > 100){
    error("Threshold greater than pattern length");
    free(msg);
    exit(1);
  }
  dif = (int) (patlen * dif / 100);

#ifdef SHOW
  printf("Pat = '%s'(%d) dif=%d\n",pat,patlen,dif);
  show_pat();
#endif

#ifdef STATS
  setup_search();
  search(ifile,dif);

#else
  
  int bsize;
  int u = ml + 1;
  bsize = 0;
  while (u > 0){
      u = u >> 1;
      bsize += 1;
  }
  if (!bsize)
      bsize ++;
  bsize++;

  setup_search(bsize, pat);
  if (argc >= 4){
    search(ifile,dif,bsize,patlen,ml);
    close(ifile);
    ifile = open(argv[3],O_RDONLY);
  }
#endif
return EXIT_SUCCESS;
}
