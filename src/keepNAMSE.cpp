/* keepNAMSE.cpp:
* The codes are built on mse.c from https://urldefense.proofpoint.com/v2/url?u=https-3A__www.physionet.org_physiotools_mse_mse.c&d=BQIGAg&c=KXXihdR8fRNGFkKiMQzstu-8MbOxd1NuZkcSBymGmgo&r=TcM8b82z3tegOJXIRc2Vz7clBmJyyKx9Zj12bbp3ukQ&m=4nLbigf5_1qJiGHdFI8ZiGDhvou2n5icBydzoa9wGqA&s=blB40K2y5Alb2amu_G4HOWJ1LecEzeB5MZJ2Xn_YsTQ&e= .
* The improvement is in handling with missing value. mse.c does not allow the input data
*  to have missing values whereas keepNAMSE.cpp does not only allow the input data to
*  missing value but also contains codes to handle with missing values in the calculation
*  process.
*/
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#include <Rcpp.h>
using namespace Rcpp;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

#define MAXSTR 1250     /* maximum string size */
#define DATA_MAX 600000 /* maximum number of data points */
#define M_MAX 10        /* maximum value of the parameter m */
#define SCALE_MAX 100    /* maximum scale value */
#define R_STEP_MAX 10   /* maximum number of different r values */
#define FILEMAX 100     /* maximum number of data files in the list file*/


/* Global variables */
char *prog, line[MAXSTR];
double SE[FILEMAX][M_MAX][R_STEP_MAX][SCALE_MAX];
double *u, *y, r_min, r_max, r_step;
int c, nlin, m_min, m_max, m_step, scale_max, scale_step, i_min, i_max;

static char file[500];
FILE *fl;
FILE *pin;
FILE *pout;

void KeepNASampEn(char** cmdstrp);
void KeepNASampEnDo(int argc, char* argv[]);
void usage();

// [[Rcpp::export]]
void KeepNASampEn(String cmdstrp)
{
  char cmdstr[200];
  strcpy_s(cmdstr,200,cmdstrp.get_cstring());
  int argc = 0;
  char *argv[30];
  const char *delim = " ";
  char *next_token;
  char *p = strtok_s(cmdstr, delim, &next_token);

  do
  {
    argv[argc++] = p;
    p = strtok_s(NULL, delim, &next_token);
  }while (p != NULL);

  KeepNASampEnDo(argc,argv);
}


void KeepNASampEnDo(int argc, char* argv[])
{

    int i, j, l, nfile, flag;
    double r, sd, tmp;

    /* Function prototypes */
    void ReadData(void);
    double *array(int n);
    double StandardDeviation(void);
	void CoarseGraining (int j);
    void SampleEntropy (int ll, double r, double sd, int j);
    void PrintResults (int nfile);

    /* Initialize variables. */
  prog = argv[0];
  scale_max = 1;
  scale_step = 1;
  m_min = 2;
  m_max = 2;
  m_step = 1;
  i_min = 0;
  i_max = 39999;
  r_min = 0.15;
  r_max = 0.15;
  r_step = 0.05;
  c = 0;
  nfile = 0;
  flag = 0;


    /* Read and interpret the command line. */
    i = 0;
    while (++i < argc && *argv[i] == '-') {
	switch(argv[i][1]) {

	  case 'F':
	    if ((fl = fopen(argv[++i], "r")) == NULL) {
		fprintf(stderr, "%s [-F]: can't open input file %s\n",
			prog, argv[i]);
		exit(1);
	    }
	    flag=1;
	    break;
case 'n':
      if ((scale_max=atoi(argv[++i])) <= 0 || scale_max > SCALE_MAX) {
        fprintf(stderr,
                "%s [-n]: maximum scale must be between 1 and %d (default: 20)\n",
                prog, SCALE_MAX);
        exit(1);
      }
      break;

    case 'a':
      if ((scale_step=atoi(argv[++i])) <= 0 || scale_step >= SCALE_MAX) {
        fprintf(stderr,
                "%s [-a]: scale increment must be between 1 and %d (default: 1)\n",
                prog, SCALE_MAX);
        exit(1);
      }
      break;

    case 'm':
      if ((m_min = atoi(argv[++i])) >= M_MAX || m_min <= 0 ) {
        fprintf(stderr,
                "%s [-m]: minimum m value must be between 1 and %d (default: 2)\n",
                prog, M_MAX);
        exit(1);
      }
      break;

    case 'M':
      if ((m_max = atoi(argv[++i])) >= M_MAX || m_max <= 0) {
        fprintf(stderr,
                "%s [-M]: maximum m value must be between 1 and %d (default: 2)\n",
                prog, M_MAX);
        exit(1);
      }
      break;

    case 'b':
      if ((m_step = atoi(argv[++i])) <= 0 || m_step > M_MAX) {
        fprintf(stderr,
                "%s [-b]: m increment must be between 1 and %d (default: 1)\n",
                prog, M_MAX);
        exit(1);
      }
      break;

    case 'r':
      if ((r_min = atof(argv[++i])) <= 0) {
        fprintf(stderr,
                "%s [-r]: minimum r must be greater than 0 (default: 0.15)\n",
                prog);
        exit(1);
      }
      break;

    case 'R':
      if ((r_max = atof(argv[++i])) <= 0) {
        fprintf(stderr,
                "%s [-R]: maximum r must be greater than 0 (default: 0.15)\n",
                prog);
        exit(1);
      }
      break;

    case 'c':
      if ((r_step=atof(argv[++i]))<=0||r_step<(r_max-r_min)/R_STEP_MAX) {
        fprintf(stderr,
                "%s [-c]: r increment must be greater than %g (default: 0.05)\n",
                prog, (r_max-r_min)/R_STEP_MAX);
      }
      exit(1);
      break;
	  
	  case 'i':
	    if ((i_min = atoi(argv[++i])) < 0) {
		fprintf(stderr,
	           "%s [-i]: minimum i must not be less than 0 (default: 0)\n",
			prog);
		exit(1);
	    }
	    break;

	  case 'I':
	    if ((i_max = atoi(argv[++i])) <= 0 || i_max <= i_min) {
		fprintf(stderr,
			"%s [-I]: maximum i must be greater than %d "
			"(default: number of data points)\n",
			prog, i_min);
		exit(1);
	    }
	    break;
    case 'z':
      if ((pin = fopen(argv[++i], "r")) == NULL) {
        fprintf(stderr, "%s [-z]: can't open input file %s\n",
                prog, argv[i]);
        exit(1);
      }
      break;
    case 'Z':
      if ((pout = fopen(argv[++i], "w")) == NULL) {
        fprintf(stderr, "%s [-z]: can't open output file %s\n",
                prog, argv[i]);
        exit(1);
      }
      break;

	  default:
	    usage();
	  exit(1);
	}
    }
  if (m_max < m_min) {
    tmp = m_max;
    m_max = m_min;
    m_min = tmp;
  }

  if (r_max < r_min){
    tmp = r_max;
    r_max = r_min;
    r_min = tmp;
  }
  
    /* Memory allocation. */
    u = array(DATA_MAX);
    y = array(DATA_MAX);

    /* Process a single data file. */
    if (flag == 0) {
	l = 0;
	nfile = 1;

	/* Read data from stdin. */
	// pin = stdin;
	ReadData();

	/* Calculate standard deviation. */
	sd = StandardDeviation();

        for (j = 1; j <= scale_max; j += scale_step){
      CoarseGraining(j);

      /* Calculate SampEn for each scale and each r value. */
      c = 0;
      for (r = r_min; r <= (r_max*1.0000000001); r += r_step){
        SampleEntropy(l, r, sd, j);
        c++;
      }
    }

	/* Print results. */
	PrintResults(nfile);
	
	fclose(pin);
	fclose(pout);
    }


//    /* Process multiple data files. */
//    if (flag == 1) {

//	/* Read the list of data files. */
//	for (l = 0; fscanf(fl, "%s", file) == 1; l++) {
//	    nfile++;   /*count the number of data files*/

//	    if ((pin = fopen(file, "r")) == NULL) {   /* open each data file */
//		fprintf(stderr, "%s : Cannot open file %s\n", prog, file);
//		exit(1);
//	    }


//	    /* Read the data. */
//	    ReadData();

	    /* Calculate the standard deviation. */
//	    sd = StandardDeviation ();

//       SampleEntropy (l, 0.15, sd);
//	}

	/* Print results. */
//	PrintResults(nfile);

//	fclose(pin);
//	fclose(fl);
//   }

    free(u);
    free(y);
//    exit(0);
}

double *array(int n)
{
    double *a;

    if ((a = (double*)calloc (n, sizeof (double))) == NULL) {
	fprintf(stderr, "%s : insufficient memory\n", prog);
	exit(2);
    }
    return (a);
}

void ReadData(void)
{

    nlin=0;
    int j = -1;

    int nmissing=0;


    while (fgets(line, MAXSTR, pin) != NULL) {
	j++;
	if (j >= i_min && j <= i_max) {
	    if(sscanf(line, "%lf", &u[nlin])==0) //sscanf return 0 if no valid number found in line
        {
            u[nlin] = DBL_MAX;
            nmissing++;
        }
	    nlin++;
	}
    }

}

double StandardDeviation(void)
{
    double sum=0.0, sum2=0.0, sd;
    int j;

    int count=0;

    for (j = 0; j < nlin; j++) {
        if(u[j]==DBL_MAX)
        {
            ;
        }
        else
        {
            sum += u[j];
            sum2 += u[j] * u[j];
            count++;
        }
    }
    sd = sqrt((sum2 - sum*sum/count)/(count - 1));
    return (sd);
}

void CoarseGraining(int j)
{
	int i, k;
	
	for (i = 0; i < nlin/j; i++) {
		y[i] = 0;
		for (k = 0; k < j; k++)
			y[i] += u[i*j+k];
		y[i] /= j;
	}
}

void SampleEntropy(int ll, double r, double sd, int j)
{
    int i, k, l, nlin_j, kk;
    int cont[M_MAX+1];
    double r_new;

    nlin_j = (nlin/j) - m_max;

    // It's meaningless to compare when you have only one point
    if(nlin_j<2)
        return;
    r_new = r*sd;

    for (i = 0; i < M_MAX; i++)
	cont[i]=0;

    for (i = 0; i < nlin_j; ++i) {
	for (l = i+1; l < nlin_j; ++l) { /*self-matches are not counted*/
	    k = 0;
	    while (k <= m_max && u[i+k]!=DBL_MAX && u[l+k]!=DBL_MAX)
            k++;

        if(k==m_max+1)
        {
            kk = 0;
            while(kk <= m_max && fabs(u[i+kk] - u[l+kk]) <= r_new)
                cont[++kk]++;
        }
	}
    }


    for (i = 1; i <= m_max; i++)

        if (cont[i+1] == 0 || cont[i] == 0)
            SE[ll][c][j][i] = -log((double)1/((nlin_j)*(nlin_j-1)));
        else
            SE[ll][c][j][i] = -log((double)cont[i+1]/cont[i]);

}



void PrintResults(int nfile)
{
  int j, m, k, l;

  fprintf (pout, "\n");
  fprintf (pout, "\nmax_line_read = %12d", nlin);
  fprintf (pout, "\n");

  for (m = m_min; m <= m_max; m += m_step)
    for (k = 0; k < c; k++) {
      fprintf (pout, "\nm = %d,  r = %.3f\n\n", m, r_min+k*r_step);
      if (nfile > 1) {
        fseek(fl, 0, SEEK_SET);
        for(l = 0; fscanf(fl, "%s", file) == 1; l++)
          fprintf (pout, "\t%.6s", file);
        fprintf (pout, "\n");
      }
      for (j = 1; j <= scale_max; j += scale_step) {
        fprintf (pout, "%d\t", j);
        for (l=0; l<nfile; l++)
        {
          fprintf (pout, "%.3lf\t", SE[l][k][j][m]);
        }
        fprintf (pout, "\n");
      }
    }
}


void usage()
{
  fprintf(stderr, "usage: %s [options]\n", prog);
  fprintf(stderr, "\nTo calculate MSE for a single data file:\n"
            "    %s <datafile >outputfile\n"
            "To calculate MSE for multiple data files:\n"
            "    %s -F listfile >outputfile\n"
            "(where listfile contains a list of data files).\n\n", prog, prog);
  fprintf(stderr, "Data files should contain a single column of numbers\n");
  fprintf(stderr, "Options may include:\n");
  fprintf(stderr, "  -a N   set scale increment to N [1-%d; default: 1]\n",
          SCALE_MAX);
  fprintf(stderr, "  -b N   set m increment to N [1-%d; default: 1]\n",
          M_MAX);
  fprintf(stderr, "  -c X   set r increment to X [>%g; default: 0.05]\n",
          (r_max-r_min)/R_STEP_MAX);
  fprintf(stderr, "  -i N   set minimum i to N [0-39998; default: 0]\n");
  fprintf(stderr, "  -I N   set maximum i to N [1-39999: default: 39999]\n");
  fprintf(stderr, "  -m N   set minimum m to N [1-%d; default: 2]\n", M_MAX);
  fprintf(stderr, "  -M N   set maximum m to N [1-%d; default: 2]\n", M_MAX);
  fprintf(stderr, "  -n N   set maximum scale to N [1-%d; default: 20]\n",
          SCALE_MAX);
  fprintf(stderr, "  -r X   set minimum r to X [>0; default: 0.15]\n");
  fprintf(stderr, "  -R X   set maximum r to X [>0; default: 0.15]\n");
  fprintf(stderr,
          "Option arguments indicated as N are integers; those shown as X may be given"
          "\nin any floating point format. \n");
}
