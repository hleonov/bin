#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PULL_GROUPS 20
#define SIGN(x) (x < 0 ? -1 : 1)

static float Temperature=298.0;
static real Tolerance;
static real Roundofftol;

typedef struct {
  int nSkip;
  char Reference[256];
  int nPull;
  int nDim;
  rvec Dims;
  char PullName[MAX_PULL_GROUPS][256];
  double UmbPos[MAX_PULL_GROUPS][3];
  double UmbCons[MAX_PULL_GROUPS][3];
} t_UmbrellaInfo;

typedef struct {
  double **Histo;
  double *k;
  double *pos;
  double *z;
  double * N;
} t_UmbrellaWindow;

typedef struct {
  int tag;
  double min;
  double max;
  double *profile;
} t_UmbrellaResult;

void tokenizer(char *in, char *token, char *rest)
{
  char ch=0;
  char buffer[256];
  int count=0, rest_count=0;
 
  while(ch != ' ' && ch != '\t' && ch != '\n' && ch!= '\'')
    {
      buffer[count]=in[count];
      count++;
      ch=in[count];
    }
  buffer[count]='\0';
  strcpy(token,buffer);
 
  while(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\'')
    {
      ch=in[count];
      count++;
    }
  count--;
 
  while(ch != '\0' || count <= 100)
    {
      rest[rest_count]=in[count];
      ch=in[count];
      count++;
      rest_count++;
    }
  rest[rest_count]='\0';
}

void read_umbrella_info(FILE *file,t_UmbrellaInfo *header)
{
  char tester[256];
  char token[256];
  char line[256];
  int i,j;

  /*  line 1 */
  fgets(line,257,file);
  for(i=0;i<2;i++)
    {tokenizer(line,token,line);}
  if(strcmp(token,"UMBRELLA")) 
    gmx_fatal(FARGS,"This does not appear to be a valid pdo file");
  tokenizer(line,token,line);
  if(strcmp(token,"3.0"))
    gmx_fatal(FARGS,"This does not appear to be a version 3.0 pdo file");

  /* line2 */
  fgets(line,257,file);
  for(i=0;i<3;i++)
    {tokenizer(line,token,line);}
#ifdef DOUBLE
  sscanf(line,"%lf%lf%lf",&(header->Dims[0]),&(header->Dims[1]),&(header->Dims[2]));
#else
  sscanf(line,"%f%f%f",&(header->Dims[0]),&(header->Dims[1]),&(header->Dims[2]));
#endif

  header->nDim = header->Dims[0] + header->Dims[1] + header->Dims[2];
  if(header->nDim!=1)
    gmx_fatal(FARGS,"Currently only supports one dimension");

  /* line3 */
  fgets(line,257,file);
  for(i=0;i<2;i++)
    {tokenizer(line,token,line);}
  sscanf(line,"%d",&(header->nSkip));

  /* line4 */
  fgets(line,257,file);
  for(i=0;i<3;i++)
    {tokenizer(line,token,line);}
  sscanf(line,"%s",header->Reference);

  /* line5 */
  fgets(line,257,file);
  for(i=0;i<5;i++)
    {tokenizer(line,token,line);}
  sscanf(line,"%d",&(header->nPull));

  if(header->nPull > MAX_PULL_GROUPS)
    gmx_fatal(FARGS,"Only %d groups are allowed",MAX_PULL_GROUPS);

  /* The rest of the header */
  for(i=0;i<header->nPull;++i) 
    {
      fgets(line,257,file);
      for(j=0;j<4;j++)
	{tokenizer(line,token,line);}
      sscanf(token,"%s",header->PullName[i]);
      tokenizer(line,token,line);
      tokenizer(line,token,line);
      for(j=0;j<header->nDim;++j) 
	{
	  tokenizer(line,token,line);
	  sscanf(token,"%lf",&(header->UmbPos[i][j]));
	}
      tokenizer(line,token,line);
      tokenizer(line,token,line);
      for(j=0;j<header->nDim;++j) 
	{
	  tokenizer(line,token,line);
	  sscanf(token,"%lf",&(header->UmbCons[i][j]));
	}
    }	
  fgets(line,257,file);
  /* fscanf(file,"%s",Buffer3); */
  /* printf("%s\n",Buffer3); */
}

void find_bounds(FILE *file, int nprof, t_UmbrellaInfo header, t_UmbrellaResult *profile)
{
  char token[256];
  char line[2048];
  int i,j,k;
  double temp;
  while (fgets(line,2049,file)!=NULL)
    {
      tokenizer(line,token,line);
      if(strcmp(token,"#") != 0)  /* if the first field is not '#' */
	for(i=0;i<header.nPull;i++)  /* now do the rest of the npull fields */
	  {
	    tokenizer(line,token,line);
	    for(j=0;j<nprof;j++)
	      {
		if(profile[j].tag == i)
		  {
		    sscanf(token,"%lf",&temp);
		    temp += header.UmbPos[i][0];
		    if(temp < profile[j].min)
		      profile[j].min = temp;
		    if(temp > profile[j].max)
		      profile[j].max = temp;
		  }
	      }
	  }
    }
}

void read_umbrella_data(FILE * file, int nprof, t_UmbrellaInfo *header,
			int bins, t_UmbrellaResult *profile,
			int fileno, t_UmbrellaWindow * win)
{
  double temp;
  char token[256];
  char line[2048];
  int i,j,k;

  /* Below we define a new structure array and use pointer 
     arithmetic.  This makes things cleaner.  We could have
     used simply win[fileno-1] for everything if we wanted, but it is
     more typing and looks messy. */

  /* Now... we set up the window structure array element, fileno-1*/
  t_UmbrellaWindow *window;
  window=win+(fileno-1);

  snew(window->Histo,nprof);
  snew(window->z,nprof);
  snew(window->k,nprof);
  snew(window->pos,nprof);
  snew(window->N,nprof);

  for(i=0;i<nprof;++i) 
    {
      window->z[i]=1;
      snew(window->Histo[i],bins);
      window->k[i]=header->UmbCons[profile[i].tag][0];
      window->pos[i]=header->UmbPos[profile[i].tag][0];
      window->N[i]=0;
    }

  /* Done with setup, now read the file */
			
  while (fgets(line,2049,file)!=NULL)
    {
      tokenizer(line,token,line);
      if(strcmp(token,"#") != 0)  /*  if the first field is not '#' */
	for(i=0;i<header->nPull;i++)  /* now do the rest of the npull fields */
	  {
	    tokenizer(line,token,line);
	    for(j=0;j<nprof;j++)
	      {
		if(profile[j].tag == i)
		  {
		    sscanf(token,"%lf",&temp);
		    temp += window->pos[j];
		    temp -= profile[j].min;
		    temp /= (profile[j].max-profile[j].min);
		    temp *= bins;
		    temp = floor(temp);
		    if((int)temp >= 0 && (int)temp < bins)
		      {
			window->Histo[j][(int)temp]++;
			window->N[j]++;
		      }
		  }
	      }
	  }
    }
/*  for(i=0;i<nprof;i++)
    for(j=0;j<bins;j++)
	printf("%f  %f\n",(double)(j+0.5)/bins*(profile[i].max-profile[i].min)+profile[i].min,window->Histo[i][j]);
*/
        
}

void calc_profile(t_UmbrellaResult *profile,t_UmbrellaWindow *window, int nWindows, int profnum, int bins)
{
  int i,j,k,l,m,n;
  double num;
  double denom;
  double U=0,temp=0;
  double TOTAL=0;
	
  for(i=0;i<bins;++i) 
    {	
      num=denom=0;
      for(j=0;j<nWindows;++j) 
	{
	  temp=(double)(i+0.5)*(profile[profnum].max-profile[profnum].min)/bins+profile[profnum].min;
	  U=0.5*window[j].k[profnum]*(window[j].pos[profnum]-temp)*(window[j].pos[profnum]-temp);
	  num+=window[j].Histo[profnum][i];  
	  /*num += window[j].Histo[profnum][i] * exp(-U/(8.31451*Temperature));*/ 
	  denom+=window[j].N[profnum]*exp(- U/(8.31451e-3*Temperature) + window[j].z[profnum]);
	}
      profile[profnum].profile[i]=num/denom;
      TOTAL+=profile[profnum].profile[i];
      /*printf("PROFILE %e\n",profile[profnum].profile[i]); */
    }
  
  /* for(i=0;i<window[0].nBin;++i) {
   *  profile[i]/=TOTAL;
   * }
   */
}

void calc_combined_profile(t_UmbrellaResult *profile,t_UmbrellaWindow *window, int nWindows, int nprof, int bins)
{
  int i,j,k,l,m,n;
  double num;
  double denom;
  double U=0,temp=0;
  double TOTAL=0;
	
  for(i=0;i<bins;++i) 
    {	
      num=denom=0;
      for(k=0;k<nprof;++k)
	for(j=0;j<nWindows;++j) 
	  {
	    temp=(double)(i+0.5)*(profile[0].max-profile[0].min)/bins+profile[0].min;
	    U=0.5*window[j].k[k]*(window[j].pos[k]-temp)*(window[j].pos[k]-temp);
	    num+=window[j].Histo[k][i];  
	    /*num += window[j].Histo[profnum][i] * exp(-U/(8.31451*Temperature));*/ 
	    denom+=window[j].N[k]*exp(- U/(8.31451e-3*Temperature) + window[j].z[k]);
	  }
      profile[0].profile[i]=num/denom;
    }
}

double calc_z(t_UmbrellaResult *profile,t_UmbrellaWindow *window, int nWindows, int profnum, int bins)
{
  int i,j,k;
  double U=0,dist=0;
  double MAX=-1e20;
  double total=0;
  double log_total;
  double temp;
  double crap;
	
  for(i=0;i<nWindows;++i) 
    {
      total=0;
      for(k=0;k<bins;++k) 
	{
	  dist=(double)(k+0.5)*(profile[profnum].max-profile[profnum].min)/bins+profile[profnum].min;
	  U=0.5*window[i].k[profnum]*(window[i].pos[profnum]-dist)*(window[i].pos[profnum]-dist);
	  total+=profile[profnum].profile[k]*exp(-U/(8.31451e-3*Temperature));
	}
      /* printf("tot %e\n",total); */
      /* log_total=-log(total); */
      /* crap=-8.314e-3*Temperature*log(total); */
      /* temp=fabs(log_total+log(window[i].z[j])); */
      /* if(temp>MAX) MAX=temp; */
      /* printf("%e\n",temp); */
      /* window[i].z[j]=total; */
      total = -log(total);
      temp = fabs(total - window[i].z[profnum]);
      if(temp > MAX) 
	MAX=temp;
      window[i].z[profnum] = total;
    }
  printf("NEW Maximum change:%e\n",MAX);
  return (MAX);
}

double calc_combined_z(t_UmbrellaResult *profile,t_UmbrellaWindow *window, int nWindows, int nprof, int bins)
{
  int i,j,k;
  double U=0,dist=0;
  double MAX=-1e20;
  double total=0;
  double log_total;
  double temp;
  double crap;
	
  for(j=0;j<nprof;++j)
    {
      for(i=0;i<nWindows;++i) 
	{
	  total=0;
	  for(k=0;k<bins;++k) 
	    {
	      dist=(double)(k+0.5)*(profile[0].max-profile[0].min)/bins+profile[0].min;
	      U=0.5*window[i].k[j]*(window[i].pos[j]-dist)*(window[i].pos[j]-dist);
	      total+=profile[0].profile[k]*exp(-U/(8.31451e-3*Temperature));
	    }
	  total = -log(total);
	  temp = fabs(total - window[i].z[j]);
	  if(temp > MAX) 
	    MAX=temp;
	  window[i].z[j] = total;
	}
    }
  printf("NEW Maximum change:%e\n",MAX);
  return(MAX);
}

int gmx_wham(int argc,char *argv[])
{
  static char *desc[] = {
    "This is an analysis program that implements the Weighted",
    "Histogram Analysis Method (WHAM).  It is intended to analyze",
    "pdo files generated by mdrun using umbrella sampling to",
    "create a potential of mean force (PMF). The options are[PAR]",
    "  [TT]-o[tt]      name of the PMF output file(s)[BR]",
    "  [TT]-hist[tt]   name of the histogram output file(s)[BR]",
    "  [TT]-min[tt]    minimum coordinate to use (single profiles only)[BR]",
    "  [TT]-max[tt]    maximum coordinate to use (single profiles only)[PAR] ",
    "Note: the program will throw out any data that is outside",
    "of the interval, [min,max]. You can use: [PAR]",
    "  [TT]-noprof[tt]  only calculate min and max[BR]",
    "  [TT]-bins[tt]  number of bins to use in calculation[BR]",
    "  [TT]-auto[tt]  automatic selection of min and max[PAR]",
    "the automatic determination of bounds is an option when only",
    "calculating the PMF for one group, but is the default if more", 
    "PMFs are calculated. The auto option makes use of the roundoff",
    "precision given by -roundtol in order to determine the way the",
    "boundaries are chosen. When using -auto, maximum usage of the",
    "data is ensured.[PAR]",
    "The parameter -tol sets a convergence level for calculating the",
    "probability distribution from the individual histograms (see",
    "the WHAM literature for details). Both -roundtol and -tol,",
    "although unrelated parameters, may be set to 0.001 for most",
    "cases (this is the default). The boolean option, -log, will",
    "automatically take the negative natural logarithm of the",
    "probability distribution function(s).[PAR]",
    "The boolean option -combine allows you to combine all of the",
    "data from the selected groups in order to calculate a single PMF.",
    "Finally, the boolean option -forceread foregoes any checks to",
    "ensure the consistency between pdo files. Use this option only if",
    "you know the details of what the pdo files contain. This option is",
    "useful if you sampled a reaction coordinate with equivalent groups having",
    "different names in each window simulation. Thus if, say, grp1 and grp2",
    "are equivalent groups in two different window simulations, and each of their",
    "data occupy the first column of the pdo file in their own respective window",
    "simulations, g_wham will pay no attention to their difference in naming and",
    "will use the data from column one without complaining. Without the -forceread",
    "option, the program will exit if it encounters inconsistent naming for a given",
    "column in separate pdo files. Also if the number of pulled",
    "groups in any input pdo file differs from the first, the program will likely",
    "give a segmentation error. It is always safer to use the default mode without",
    "forced reading and to have uniform input pdo files. Use the -forceread option",
    "at your own disgression.[PAR]",
    "NOTE: though this program calculates the PMF for multiple groups, each PMF",
    "is one-dimensional. The assumption here is that there is no correlation",
    "among the groups. This applies also to combining data from multiple groups.",
    "All groups sampled must have no correlation in the window simulations."};
  
  static float min=0;
  static float max=0;
  static int bins=100;
  static bool Bnoprof=TRUE;
  static bool Bautobounds=FALSE;
  static bool Blog=FALSE;
  static bool Bcombine=FALSE;	
  static bool Bforceread=FALSE;
  static int n=1;
  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-min", FALSE, etREAL, {&min},
      "Lower bound in profile. Only for calculating a single profile."},
    { "-max", FALSE, etREAL, {&max},
      "Upper bound in profile. Only for calculating a single profile."},
    { "-auto", FALSE, etBOOL, {&Bautobounds},
      "Automatic determination of bounds (min and max)."},
    { "-bins",FALSE, etINT, {&bins},
      "Number of bins in profile"},
    { "-prof", FALSE, etBOOL, {&Bnoprof},
      "Only calculate min and max"},
    { "-temp", FALSE, etREAL, {&Temperature},
      "Temperature"},
    { "-tol", FALSE, etREAL, {&Tolerance},
      "Tolerance"},
    { "-roundtol", FALSE, etREAL, {&Roundofftol},
      "Precision for rounding off histogram boundaries (nm)."},
    { "-combine", FALSE, etBOOL, {&Bcombine},
      "Combine data from all selected groups into first group"},
    { "-forceread", FALSE, etBOOL, {&Bforceread},
      "Force the program to read all pdo files based on 1st input file without checking for consistency"},
    { "-log", FALSE, etBOOL, {&Blog},
      "Calculate the log of each profile before printing"}
  };
  
  t_filenm fnm[] = {
    { efXVG, "-o", "profile", ffWRITE },    /* output file for profile */
    { efXVG, "-hist","histo", ffWRITE}	    /* output file for histograms */
  };
  
  int i,j,k,l,nprof=0;
  real temp=0.0;
  t_UmbrellaInfo header, header1;
  t_UmbrellaWindow * window=NULL;
  t_UmbrellaResult *profile=NULL;
  /*double *profile;*/
  bool flag=FALSE;

    /* output histograms */
  FILE * histout;
  FILE * profout;
  char histfilename[255], proffilename[255], tempstr[255];
  int counter=0;

   /* for reading in headers of pdo files */
  FILE  * file;    
  char Buffer[256];    
  char * Path=NULL;
  
#define NFILE asize(fnm)

  Tolerance=0.001;
  Roundofftol=0.001;
  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_VIEW|PCA_NOEXIT_ON_ARGS,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
		    
  /* Check to see what we get for command line arguments */
  if(argc<2) gmx_fatal(FARGS,"You need to specify a series of pdo files as input");
  for(i=1;i<argc;++i) {
    if(!fexist(argv[i]))
      gmx_fatal(FARGS,"Unable to open file %s b/c it does not exist.",argv[i]);
  }

  /* Now we need to process the files */
  if(Bnoprof)
    snew(window,argc-1);

  /* Loop over all files */
  for(i=1;i<argc;++i) {
    /* read in the headers */
		
    if(!(Path=getenv("GMX_PATH_GZIP")))
      sprintf(Buffer,"/bin/gunzip -c < %s",argv[i]);
    else
      sprintf(Buffer,"%s/gunzip -c < %s",Path,argv[i]);

    if((file=popen(Buffer,"r"))==NULL)
      gmx_fatal(FARGS,"Unable to open file %s.\n",argv[i]);
    printf("Opening file %s.\n",argv[i]);
    read_umbrella_info(file,&header);

    /* Check for pdo file consistency */
    if(i == 1)
      {
	header1 = header;
	if (header1.nPull > 1)
	  {
	    if (!Bcombine)
	       printf("There are %d pulled groups. How many profiles?: ",header1.nPull);
	    else
	       {
	         printf("\n** Combine option is TRUE:\n");
	         printf("   There are %d pulled groups. How many pull groups\n",header1.nPull); 
	         printf("   to be combined into one profile?: "); 
	       }
	    do { 
	      scanf("%d",&nprof);
	      if(nprof > header1.nPull || nprof <= 0)
		printf("You must choose a number less than or equal to %d, but greater than 0.\n",header1.nPull);
	    } while (nprof <= 0 || nprof > header1.nPull); 
	    /* nprof = 1; */
	    snew(profile,nprof);
	    printf("Groups Sampled:\n");
	    for(j=0;j<header1.nPull;j++)
	      printf("%d --> %s\n",j,header1.PullName[j]);
	    printf("Which group(s)?:\n");
	    for(j=0;j<nprof;++j)
	      {
	      scanf("%d",&(profile[j].tag));	    
	      printf("I got %s\n",header1.PullName[profile[j].tag]);
	      profile[j].min = 1e20;
	      profile[j].max = -1e20;
	      }
	    printf("\nReading file(s) ... Please wait.\n");
	  }
	else if(header1.nPull == 1)
	  {
	    nprof = 1;
	    snew(profile,nprof);
	    profile[0].min = 1e20;
	    profile[0].max = -1e20;
	    profile[0].tag = 0;
	  }
	else
	  gmx_fatal(FARGS,"No groups sampled according to input file(s).\n");
      }  
    else if (!Bforceread)
      {
	if (header.nSkip != header1.nSkip)
	  gmx_fatal(FARGS,"pdo file %d does not appear to jive with pdo file 1 -- see nSkip",i);
	if (strcmp(header.Reference,header1.Reference) != 0)
	  gmx_fatal(FARGS,"pdo file %d does not appear to jive with pdo file 1 -- see reference",i);
	if (header.nPull != header1.nPull)
	  gmx_fatal(FARGS,"pdo file %d does not appear to jive with pdo file 1 -- see nPull",i);
	if (header.nDim != header1.nDim)
	  gmx_fatal(FARGS,"pdo file %d does not appear to jive with pdo file 1 -- see pulldims",i);
	for(j=0;j<header1.nDim;++j)
	  if(header.Dims[j] != header1.Dims[j])
	    gmx_fatal(FARGS,"pdo file %d does not appear to jive with pdo file 1 -- see pulldims",i);
	for(j=0;j<header1.nPull;++j)
	  if(strcmp(header.PullName[j],header1.PullName[j]) != 0)
	    gmx_fatal(FARGS,"pdo file %d does not appear to jive with pdo file 1 -- see pull groups",i);
      }
	
    find_bounds(file,nprof,header,profile);
    pclose(file); 
  }

  if(!Bnoprof) 
    {
      printf("\nOutputting boundaries for profiles:\n");
      for(i=0;i<nprof;i++)
	printf("%s : Min = %e : Max = %e\n",header1.PullName[profile[i].tag],profile[i].min,profile[i].max);
      exit(0);
    }
  else if(Bautobounds || nprof > 1)
    {
      printf("\nNo bounds specified.\nUsing automatic determination of profile boundaries, min and max.\n");
      for(i=0;i<nprof;i++)
	{
	  if(profile[i].min > 0)
	    profile[i].min = profile[i].min - fmod(profile[i].min,Roundofftol) + Roundofftol;
	  if(profile[i].min <= 0)
	    profile[i].min = profile[i].min - fmod(profile[i].min,Roundofftol);
	  if(profile[i].max >= 0)
	    profile[i].max=profile[i].max - fmod(profile[i].max,Roundofftol);
	  if(profile[i].max < 0)
	    profile[i].max=profile[i].max - fmod(profile[i].max,Roundofftol) - Roundofftol;
	  printf("%s : Min = %e : Max = %e\n",header1.PullName[profile[i].tag],profile[i].min,profile[i].max);
	}
    }
  else
    {
      printf("\nUsing input min and max values for profile\n");
      profile[0].min = min;
      profile[0].max = max;
    }

  if(Bcombine)
    {
      for(i=0;i<nprof;i++)
	{
	  if (profile[i].min < profile[0].min)
	    profile[0].min = profile[i].min;
	  if (profile[i].max > profile[0].max)
	    profile[0].max = profile[i].max;
	}
      for(i=0;i<nprof;i++)
	{
	  profile[i].min = profile[0].min;
	  profile[i].max = profile[0].max;
	}
      printf("\n*** Combine option is set to TRUE ***\n");
      printf("*** I will combine data from all selected groups into the first group, %s\n",header1.PullName[profile[0].tag]);
      printf("*** New combined bounds for the profile are : Min = %e : Max = %e\n", profile[0].min, profile[0].max);
    }
    
  printf("\nCalculating stuff ... Please wait.\n\n");
  /* Loop over all files AGAIN*/
  for(i=1;i<argc;++i) 
    {
      /* read in the headers */
      if(!(Path=getenv("GMX_PATH_GZIP")))
	sprintf(Buffer,"/bin/gunzip -c < %s",argv[i]);
      else
	sprintf(Buffer,"%s/gunzip -c < %s",Path,argv[i]);

      if((file=popen(Buffer,"r"))==NULL)
	gmx_fatal(FARGS,"Unable to open file %s",argv[i]);
      read_umbrella_info(file,&header);
      read_umbrella_data(file,nprof,&header,bins,profile,i,window);
      pclose(file);
    }



  /* Do output here */
  
  if(Bcombine)
    {    
      printf("\n*** Combine option is TRUE. I am now combining all groups into the first group\n");  
      strcpy(tempstr,header1.PullName[profile[0].tag]);
      strcpy(histfilename,strcat(strcat(tempstr,"_"),opt2fn("-hist",NFILE,fnm)));
      histout=ffopen(histfilename,"w");
      for(j=0;j<bins;++j)
	{
	  for(i=0;i<nprof;++i)
       	    {
	      fprintf(histout,"%e\t",(double)(j+0.5)/bins*(profile[i].max-profile[i].min)+profile[i].min);
	      for(k=1;k<argc;++k) 
		{
		  fprintf(histout,"%e\t",window[k-1].Histo[i][j]);
		}
	    }
	  fprintf(histout,"\n");
	}
      ffclose(histout);
    }
  else
    {
      for(i=0;i<nprof;++i)
	{
	  strcpy(tempstr,header1.PullName[profile[i].tag]);
	  strcpy(histfilename,strcat(strcat(tempstr,"_"),opt2fn("-hist",NFILE,fnm)));
	  histout=ffopen(histfilename,"w");
	  for(j=0;j<bins;++j) 
	    {
	      fprintf(histout,"%e\t",(double)(j+0.5)/bins*(profile[i].max-profile[i].min)+profile[i].min);
	      for(k=1;k<argc;++k) 
		{
		  fprintf(histout,"%e\t",window[k-1].Histo[i][j]);
		}
	      fprintf(histout,"\n");
	    }
	  ffclose(histout);
	}
    }
  
  /* Calculate profile */
  if(Bcombine)
     {
      snew(profile[0].profile,bins);
      printf("Doing Combined Profile of All Groups Into Group 0 --> %s\n",header1.PullName[profile[0].tag]);
      do 
	{
	  calc_combined_profile(profile,window,argc-1,nprof,bins);
	} while(calc_combined_z(profile, window, argc-1, nprof, bins) > Tolerance);
      nprof = 1;
     }
  else
    {
      for(i=0;i<nprof;i++)
	{
	  snew(profile[i].profile,bins);	
	  printf("Doing Profile %d --> %s\n",i,header1.PullName[profile[i].tag]);
	  do 
	    {
	      calc_profile(profile,window,argc-1,i,bins);
	    } while(calc_z(profile, window, argc-1, i, bins) > Tolerance);
	}
    }

  for(i=0;i<nprof;i++)
    {
      strcpy(tempstr,header1.PullName[profile[i].tag]);
      strcpy(proffilename,strcat(tempstr,"_"));
      strcpy(proffilename,strcat(proffilename,opt2fn("-o",NFILE,fnm)));
      printf("profile name is %s\n",proffilename);
      profout=ffopen(proffilename,"w");
      for(j=0;j<bins;j++)
	{
	  if(Blog)
	    {
	      if(profile[i].profile[j] != 0.0)
		temp = -log(profile[i].profile[j]);
	      fprintf(profout,"%e\t%e\n",(double)(j+0.5)/bins*(profile[i].max-profile[i].min)+profile[i].min,temp);
	    }
	  else
	    {
	      fprintf(profout,"%e\t%e\n",(double)(j+0.5)/bins*(profile[i].max-profile[i].min)+profile[i].min,profile[i].profile[j]);
	      /* printf("i = %d : j = %d : %e\t\n",i,j,(double)(j+0.5)/bins*(profile[i].max-profile[i].min)+profile[i].min);*/
	    }
	}
      ffclose(profout);	
    }
  
  thanx(stderr);
  
  return 0;
}
