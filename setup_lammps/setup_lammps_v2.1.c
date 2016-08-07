#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <stddef.h>
#include <ctype.h>

#define MAX 48
#define RECORD 80
#define BUFFER 83
#define BPAD 5.0
#define MAXPARM 10
#define PBUFF 50
typedef struct params
{

  int nparms;
  char parms[MAXPARM][PBUFF];

}PARAMS;



typedef struct sysdat
{
  int nats,nbnds,nangs,ntops,total_ats,total_bnds,total_angs,total_improps,foundatoms,boxinfo;
  int uniq_nats,uniq_nbnds,uniq_nangs,uniq_imps;
  int *param_bnds,*param_angs,*param_imps;
  int ndihes,total_dihes,uniq_ndihes,*param_dihes;
  int nimp;
  double *coordx,*coordy,*coordz;
  double boxx,boxy,boxz;
} SYSDAT;

typedef struct topdat
{
  int *bndndx1,*bndndx2,*bndtype,*angndx1,*angndx2,*angndx3,*angtype;
  int *improp_func,*impropndx1,*impropndx2,*impropndx3,*impropndx4;
  int *index, nat,nbnd,nang, nimprop, nmol,*parm_atomtype;
  int ndihe,*dihendx1,*dihendx2,*dihendx3,*dihendx4,*dihetype;
  int nimp, *imptype;
  double *mass,*charge;
  char (*type)[5],(*segid)[5],(*resname)[5];
}TOPDAT;

typedef struct database
{
  
  int nvdwtype,nbndtype,nangtype,nimptype;
  PARAMS* p_dihe;
  PARAMS* p_imp;
  PARAMS* p_bnd;
  PARAMS* p_ang;
  PARAMS* p_pair;
  int ndihetype,*d_dihe,*n_dihe;
  char (*vdwtype1)[5],(*vdwtype2)[5],(*vdwstyle)[7];
  char (*bndtype1)[5],(*bndtype2)[5];
  char (*angtype1)[5],(*angtype2)[5],(*angtype3)[5];
  char (*dihetype1)[5],(*dihetype2)[5],(*dihetype3)[5],(*dihetype4)[5];
  char (*imptype1)[5],(*imptype2)[5],(*imptype3)[5],(*imptype4)[5];
 //lammps style parameters (pair_style, bond_style, etc.)
  PARAMS pstyle;
  PARAMS bstyle;
  PARAMS astyle;
  PARAMS dstyle;
  PARAMS istyle;
  PARAMS kstyle;
  PARAMS sbonds;			
}DATABASE;

/* CODE STARTS HERE */
/* CODE STARTS HERE */
/* CODE STARTS HERE */
/* CODE STARTS HERE */

//functions to operate on PARAMS struct
set_PARAMnum(PARAMS* prm, int n){
	prm->nparms = n;
	return;
}

push_toPARAM(PARAMS* prm, char* w, int ind){
	//printf("copying string %s to parms index %d \n",w,ind);
	strcpy(prm->parms[ind],w);
//	printf("prm->parms[ind]->w %s \n",prm->parms[ind]);
	return;
}
copySinglePARAM_toPARAMstarAt(PARAMS* prma, PARAMS* prmb, int ind){
	prmb[ind].nparms = prma->nparms;
	//printf("prmb[ind].nparms %d \n",prmb[ind].nparms);
	int i;
	for (i = 0; i < prma->nparms; i += 1)
	{
		//printf("i %d\n",i);
		//printf("%s\n",prma->parms[i]);
		strcpy(prmb[ind].parms[i], prma->parms[i]);
	}
	return;
}
//------------------

guess_box(SYSDAT *sdat)
{
	printf("NOTE: Guessing box dimensions based on provided coordinates.\n");
	double xlo=0.0, xhi=0.0, ylo=0.0, yhi=0.0, zlo=0.0, zhi=0.0;
	double xcur=0.0, ycur=0.0, zcur=0.0;
	int i;
	int nat = sdat->foundatoms;
	//sysdat->foundatoms=0;
	//printf("found atoms %i \n",nat);
	for (i = 0; i < nat; i += 1)
	{
		//printf("i %i", i);
		//get values of coordinates
		
		xcur = sdat->coordx[i];
		ycur = sdat->coordy[i];
		zcur = sdat->coordz[i];
		//printf("xcur %f \n",xcur);
		if (xcur<xlo)
		{
			xlo=xcur;
		}
		else if (xcur>xhi)
		{
			xhi=xcur;
		}
		if (ycur<ylo)
		{
			ylo=ycur;
		}
		else if (ycur>yhi)
		{
			yhi=ycur;
		}
		if (zcur<zlo)
		{
			zlo=zcur;
		}
		else if (zcur>zhi)
		{
			zhi=zcur;
		}
		//sysdat->foundatoms++;
	}
	sdat->boxx=xhi-xlo + 2.0*BPAD;
	sdat->boxy=yhi-ylo + 2.0*BPAD;
	sdat->boxz=zhi-zlo + 2.0*BPAD;
	//printf("boxx %f",sdat->boxx);
	sdat->boxinfo=1;
	return;
}


read_pdb(char *filename,SYSDAT *sysdat)
{
  char string[BUFFER];
  FILE *fpin;
  char tmp[RECORD+3]; /* space for line + cr + lf + NUL */
  char keep, *s;
  char numstr[50]; /* store all fields in one array to save memset calls */
  
  sysdat->boxinfo=0;
  sysdat->foundatoms=0;
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }

  while(fgets(string, RECORD + 2, fpin) != NULL) {
    if (!strncmp(string, "CRYST1", 6))
    {
      printf("FOUND BOXSIZE DATA\n");
      sysdat->boxinfo=1;
      memset(tmp, 0, sizeof(tmp));
      strncpy(tmp, string, RECORD);
      /* advance to the 6 pos and then store and zero */
      s = tmp+6 ; keep = tmp[15]; tmp[15] = 0;
      sysdat->boxx = (float) atof(s);
      /* advance 15 places and replace the first char */
      s = tmp+15; *s = keep; keep = tmp[24]; tmp[24] = 0;
      sysdat->boxy = (float) atof(s);
      
      s = tmp+24; *s = keep; tmp[33] = 0;
      sysdat->boxz = (float) atof(s);      

    }
    else if (!strncmp(string, "ATOM ",  5) || !strncmp(string, "HETATM", 6))
    {
      memset(numstr, 0, sizeof(numstr));
      
      strncpy(numstr, string + 30, 8);
      sysdat->coordx[sysdat->foundatoms] = (float) atof(numstr);
      
      strncpy(numstr+10, string + 38, 8);
      sysdat->coordy[sysdat->foundatoms] = (float) atof(numstr+10);
      
      strncpy(numstr+20, string + 46, 8);
      sysdat->coordz[sysdat->foundatoms] = (float) atof(numstr+20);
      sysdat->foundatoms++;
    }
  }
  if(!sysdat->foundatoms){
    printf("ERROR: DID NOT FIND ANY ATOMS IN THE PDB FILE\n");
    exit(1);
  }
  if(!sysdat->boxinfo){
    printf("WARNING: DID NOT FIND CELL SIZE!!\n");
    printf("BOX SIZE WILL HAVE BE SET BY HAND\n");
  }
  return;
}


getwords(char *line, char *words[], int maxwords)
{
char *p = line;
int nwords = 0;

while(1)
	{
	while(isspace(*p))
		p++;

	if(*p == '\0')
		return nwords;

	words[nwords++] = p;

	while(!isspace(*p) && *p != '\0')
		p++;

	if(*p == '\0')
		return nwords;

	*p++ = '\0';

	if(nwords >= maxwords)
		return nwords;
	}
}

getwordstochar(char *line, char *words[], int maxwords, char c)
{
char *p = line;
int nwords = 0;

while(1)
	{
	while(isspace(*p))
		p++;

	if(*p == '\0')
		return nwords;
	if(*p == c)
		return nwords;

	words[nwords++] = p;

	while( !isspace(*p) && (*p != '\0' && *p != c ) )
		p++;

	if(*p == '\0')
		return nwords;
	if(*p == c)
		return nwords;
	*p++ = '\0';

	if(nwords >= maxwords)
		return nwords;
	}
}


void read_xyz(char *filename,SYSDAT *sysdat)
{
  char string[BUFFER];
  FILE *fpin;
  int boxinfo,foundatoms;
  char numstr[50]; /* store all fields in one array to save memset calls */
  
  sysdat->boxinfo=0;
  sysdat->foundatoms=0;
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  int line = 0;
  while(fgets(string, RECORD + 2, fpin) != NULL) {
	if (line==0)
	{
		printf("Note: standard xyz format does not contain box size data. \n");
		//printf("Box sizes will need to set by hand. \n");
		
	}
    else if (line>1)
    {
    	
        int nwords;
		char *words[5];
		nwords = getwords(string, words, 5);
		if (nwords>4){
			fprintf(stderr,"ERROR: TOO MANY ARGUMENTS IN XYZ FILE \n");
			printf("Check that coordfile is correct or has correct xyz format. \n");
			exit(1);
		}
		else if (nwords<4){
			fprintf(stderr,"ERROR: TOO FEW ARGUMENTS IN XYZ FILE \n");
			printf("Check that coordfile is correct or has correct xyz format. \n");
			if (nwords==0 && line>2)
			{
				printf("Hint: check for trailing (empty) lines in the coordinate file");
			}
			exit(1);
		}
	
      	sysdat->coordx[sysdat->foundatoms] = (float) atof(words[1]);
		//double xc = sysdat->coordx[sysdat->foundatoms];
		//printf("xc in read xyz: %d",xc);
       	sysdat->coordy[sysdat->foundatoms] = (float) atof(words[2]);
       	sysdat->coordz[sysdat->foundatoms] = (float) atof(words[3]);
      	sysdat->foundatoms++;
    }
    line++;
  }
  if(!sysdat->foundatoms){
    printf("ERROR: DID NOT FIND ANY ATOMS IN THE XYZ FILE\n");
    exit(1);
  }
 // if(!sysdat->boxinfo){
  // 	guess_box(sysdat);
  //}
  return;
}

void read_gcd(char *filename,SYSDAT *sysdat)
{
  char string[BUFFER];
  FILE *fpin;
  char numstr[50]; /* store all fields in one array to save memset calls */
  
  sysdat->boxinfo=0;
  sysdat->foundatoms=0;
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  int line = 0;
  while(fgets(string, RECORD + 2, fpin) != NULL) {
	if (line==0)
	{
		printf("Note: unformatted coordinates do not contain box size data. \n");
		printf("Box sizes will need to set by hand. \n");
		
	}
    
    	
        int nwords;
		char *words[5];
		nwords = getwords(string, words, 5);
		//printf("nwords %i\n",nwords);
		if (nwords>3){
			fprintf(stderr,"ERROR: TOO MANY ARGUMENTS IN COORD FILE \n");
			printf("Check that coordfile is correct or has general coordinate format. \n");
			printf("i.e. :\n");
			printf(" x1 y1 z1 \n");
			printf(" x2 y2 z2 \n");
			printf(".........\n");
			printf(" xN yN zN \n");
			exit(1);
		}
		else if(nwords<3){
			fprintf(stderr,"ERROR: TOO FEW ARGUMENTS IN COORD FILE \n");
			printf("Check that coordfile is correct or has general coordinate format. \n");
			printf("i.e. :\n");
			printf(" x1 y1 z1 \n");
			printf(" x2 y2 z2 \n");
			printf(".........\n");
			printf(" xN yN zN \n");
			if (nwords==0 && line>0)
			{
				printf("Hint: check for trailing (empty) lines in the coordinate file");
			}
			exit(1);

		}
	
      	sysdat->coordx[sysdat->foundatoms] = (float) atof(words[0]);
       	sysdat->coordy[sysdat->foundatoms] = (float) atof(words[1]);
       	sysdat->coordz[sysdat->foundatoms] = (float) atof(words[2]);
      	sysdat->foundatoms++;
    
    line++;
  }
  if(!sysdat->foundatoms){
    printf("ERROR: DID NOT FIND ANY ATOMS IN THE GCD FILE\n");
    exit(1);
  }
 /* if(!sysdat->boxinfo){
    printf("WARNING: DID NOT FIND CELL SIZE!!\n");
    printf("BOX SIZE WILL HAVE BE SET BY HAND\n");
  *///}
  return;
}

void read_coords(char *filename,DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED, char *ctype)
{
  FILE *fpout;
  int i,j,k,atindex,molindex,bondindex,angleindex,offset,impropindex;
  int diheindex;
  
/* will call modules depending on file type to be read */
/* right now only pdb capability */

  if(ISCHARGED)
    printf("FOUND CHARGES!!\n");
  else
    printf("DID NOT FIND CHARGES!!\n");
    
  
  if((fpout = fopen("DATA.FILE","w")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open DATA.FILE\n");
    exit(1);
  }

  //printf("Will read coordinates from file type ",ctype,"\n");
 // fprintf(stderr,"Will read coordinates from file type %s\n",ctype);

  if (!strncmp(ctype, "xyz",  3))
  {
  	/* call for reading xyz */
	 read_xyz(filename,sysdat);
     printf("READ XYZ FILE\n");
  }
  else if(!strncmp(ctype, "pdb",  3))
  {	
	/* call for reading pdb */
  	read_pdb(filename,sysdat);
  	printf("READ PDB FILE\n");
  }
  else if(!strncmp(ctype, "gcd",  3))
  {	
	/* call for reading pdb */
  	read_gcd(filename,sysdat);
  	printf("READ GCD FILE\n");
  }
  else{
	fprintf(stderr,"ERROR: COORDINATE FILE TYPE %s\n",ctype);
	printf("--DOES NOT MATCH ALLOWED TYPES \n");
	printf("--VALID TYPES: pbd, xyz \n");
	exit(1);
  }
  if(!sysdat->boxinfo){
    printf("WARNING: DID NOT FIND CELL SIZE!!\n");
   // printf("BOX SIZE WILL HAVE BE SET BY HAND\n");
	guess_box(sysdat);
  }
  if(sysdat->total_ats!=sysdat->foundatoms)
  {
    printf("ERROR: NUMBER OF ATOMS READ FROM TOPOLOGY AND COMMAND LINE\n");
    printf("DOES NOT MATCH THE NUMBER OF ATOMS FOUND IN THE COORD FILE\n");
  }
  
  fprintf(fpout,"LAMMPS description\n");
  fprintf(fpout,"\n");
  fprintf(fpout,"%-7d atoms\n",sysdat->total_ats);

  if(sysdat->nbnds>0)
    fprintf(fpout,"%-7d bonds\n",sysdat->total_bnds);

  if(sysdat->nangs>0)
    fprintf(fpout,"%-7d angles\n",sysdat->total_angs);

  if(sysdat->ndihes>0)
    fprintf(fpout,"%-7d dihedrals\n",sysdat->total_dihes);

  if(sysdat->total_improps>0)
    fprintf(fpout,"%-7d impropers\n",sysdat->total_improps);
  
  fprintf(fpout,"\n");

  fprintf(fpout,"%-7d atom types\n",sysdat->uniq_nats);

  if(sysdat->nbnds>0)
    fprintf(fpout,"%-7d bond types\n",sysdat->uniq_nbnds);

  if(sysdat->nangs>0)
    fprintf(fpout,"%-7d angle types\n",sysdat->uniq_nangs);

  if(sysdat->ndihes>0)
    fprintf(fpout,"%-7d dihedral types\n",sysdat->uniq_ndihes);

  if(sysdat->total_improps>0)
    fprintf(fpout,"%-7d improper types\n",sysdat->uniq_imps);

  fprintf(fpout,"\n");

  if(sysdat->boxinfo){
    fprintf(fpout,"%lf %lf xlo xhi\n",(sysdat->boxx/2.0)*-1.0,sysdat->boxx/2.0);
    fprintf(fpout,"%lf %lf ylo yhi\n",(sysdat->boxy/2.0)*-1.0,sysdat->boxy/2.0);
    fprintf(fpout,"%lf %lf zlo zhi\n",(sysdat->boxz/2.0)*-1.0,sysdat->boxz/2.0);
  } else {
    fprintf(fpout,"0.0 0.0 xlo xhi\n");
    fprintf(fpout,"0.0 0.0 ylo yhi\n");
    fprintf(fpout,"0.0 0.0 zlo zhi\n");
  }

  fprintf(fpout,"\n");
  fprintf(fpout,"Atoms\n");
  fprintf(fpout,"\n");
  atindex=0;
  molindex=0;
    
  for(i=0;i<sysdat->ntops;i++) {
    for(j=0;j<topdat[i].nmol;j++){
      molindex++;
      for(k=0;k<topdat[i].nat;k++){
        atindex++;

        if(ISCHARGED){
          if(sysdat->total_bnds>0){
            fprintf(fpout,"%u %d %d %7.4f %7.4f %7.4f %7.4f # %s\n",atindex,molindex,topdat[i].parm_atomtype[k]+1,topdat[i].charge[k] ,sysdat->coordx[atindex-1],
               sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]);
          }
          else{
            fprintf(fpout,"%u %d %7.4f %7.4f %7.4f %7.4f # %s\n",atindex,topdat[i].parm_atomtype[k]+1,topdat[i].charge[k],sysdat->coordx[atindex-1],
               sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]);
          }
        } else {
          if(sysdat->total_bnds>0){
            fprintf(fpout,"%u %d %d %7.4f %7.4f %7.4f # %s\n",atindex,molindex,topdat[i].parm_atomtype[k]+1 ,sysdat->coordx[atindex-1],
               sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]);
          }
          else{
            fprintf(fpout,"%u %d %7.4f %7.4f %7.4f # %s\n",atindex,topdat[i].parm_atomtype[k]+1,sysdat->coordx[atindex-1],
               sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]);
          }
        }
      }
    }
  }

  if(sysdat->total_bnds>0){
    fprintf(fpout,"\n");
    fprintf(fpout,"Bonds\n");
    fprintf(fpout,"\n");
    bondindex=0;
    offset=0;
    for(i=0;i<sysdat->ntops;i++) {
      for(j=0;j<topdat[i].nmol;j++){
        molindex++;
        for(k=0;k<topdat[i].nbnd;k++){
          bondindex++;
          fprintf(fpout,"%d %d %d %d # %s %s\n",bondindex,topdat[i].bndtype[k]+1,
             topdat[i].bndndx1[k]+(j*topdat[i].nat)+offset,topdat[i].bndndx2[k]+(j*topdat[i].nat)+offset,
             topdat[i].type[topdat[i].bndndx1[k]-1],topdat[i].type[topdat[i].bndndx2[k]-1]);
        }
      }
      offset+=topdat[i].nmol*topdat[i].nat;
    }
  }

  /* if ther are angles....print them out */
  
  if(sysdat->total_angs>0){
    fprintf(fpout,"\n");
    fprintf(fpout,"Angles\n");
    fprintf(fpout,"\n");
    angleindex=0;
    offset=0;
    for(i=0;i<sysdat->ntops;i++) {
      for(j=0;j<topdat[i].nmol;j++){
        molindex++;
        for(k=0;k<topdat[i].nang;k++){
          angleindex++;
          fprintf(fpout,"%d %d %d %d %d # %s %s %s\n",angleindex,topdat[i].angtype[k]+1,
             topdat[i].angndx1[k]+(j*topdat[i].nat)+offset,topdat[i].angndx2[k]+(j*topdat[i].nat)+offset,
             topdat[i].angndx3[k]+(j*topdat[i].nat)+offset,topdat[i].type[topdat[i].angndx1[k]-1],
             topdat[i].type[topdat[i].angndx2[k]-1],topdat[i].type[topdat[i].angndx3[k]-1]);
        }
      }
      offset+=topdat[i].nmol*topdat[i].nat;
    }
  }


  /* if ther are dihedrals....print them out */
  
  if(sysdat->total_dihes>0){
    fprintf(fpout,"\n");
    fprintf(fpout,"Dihedrals\n");
    fprintf(fpout,"\n");
    diheindex=0;
    offset=0;
    for(i=0;i<sysdat->ntops;i++) {
      for(j=0;j<topdat[i].nmol;j++){
        molindex++;
        for(k=0;k<topdat[i].ndihe;k++){
          diheindex++;
          fprintf(fpout,"%d %d %d %d %d %d # %s %s %s %s\n",diheindex,topdat[i].dihetype[k]+1,
             topdat[i].dihendx1[k]+(j*topdat[i].nat)+offset,topdat[i].dihendx2[k]+(j*topdat[i].nat)+offset,
             topdat[i].dihendx3[k]+(j*topdat[i].nat)+offset,topdat[i].dihendx4[k]+(j*topdat[i].nat)+offset,
	     topdat[i].type[topdat[i].dihendx1[k]-1],topdat[i].type[topdat[i].dihendx2[k]-1],
	     topdat[i].type[topdat[i].dihendx3[k]-1],topdat[i].type[topdat[i].dihendx4[k]-1]);
        }
      }
      offset+=topdat[i].nmol*topdat[i].nat;
    }
  }

    /* if ther are impropers....print them out */
  
  if(sysdat->total_improps>0){
    fprintf(fpout,"\n");
    fprintf(fpout,"Impropers\n");
    fprintf(fpout,"\n");
    impropindex=0;
    offset=0;
    for(i=0;i<sysdat->ntops;i++) {
      for(j=0;j<topdat[i].nmol;j++){
        molindex++;
        for(k=0;k<topdat[i].nimprop;k++){
          impropindex++;
          fprintf(fpout,"%d %d %d %d %d %d # %s %s %s %s\n",impropindex,topdat[i].imptype[k]+1,
             topdat[i].impropndx1[k]+(j*topdat[i].nat)+offset,topdat[i].impropndx2[k]+(j*topdat[i].nat)+offset,
             topdat[i].impropndx3[k]+(j*topdat[i].nat)+offset,topdat[i].impropndx4[k]+(j*topdat[i].nat)+offset , topdat[i].type[topdat[i].impropndx1[k]-1],
             topdat[i].type[topdat[i].impropndx2[k]-1],topdat[i].type[topdat[i].impropndx3[k]-1],topdat[i].type[topdat[i].impropndx4[k]-1]);
        }
      }
      offset+=topdat[i].nmol*topdat[i].nat;
    }
  }

return;
}


/* Decide which params will be needed */


void get_unique(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat,int ISCHARGED)
{
  FILE *fpout;
  int i,j,k;
  char (*uniq_atype)[5];
  double *uniq_mass,*uniq_charge;
  int uniq_nats,uniq_bnds,uniq_angs,ikeep,ifound,vdwtmp;
  int *bnd_params,*ang_params,*ang_vdw,datndx;
  int uniq_dihes,*dihe_params;
  int uniq_imps, *imp_params;	

  if((fpout = fopen("PARM.FILE","w")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open PARM.FILE\n");
    exit(1);
  }
  
  uniq_atype=calloc(sysdat->nats,sizeof(*uniq_atype)); 
  uniq_mass=(double *)calloc(sysdat->nats,sizeof(double));
  uniq_charge=(double *)calloc(sysdat->nats,sizeof(double));
  bnd_params=(int *)calloc(database->nbndtype,sizeof(int));
  ang_params=(int *)calloc(database->nangtype,sizeof(int));
  ang_vdw=(int *)calloc(database->nangtype,sizeof(int));
  dihe_params=(int *)calloc(database->ndihetype,sizeof(int));
  imp_params=(int *)calloc(database->nimptype,sizeof(int));

  fprintf(fpout,"# Generated by setup_lammps\n");
  fprintf(fpout,"\n");


  fprintf(fpout,"pair_style    ");
  int o;	
  for (o = 0; o < database->pstyle.nparms; o += 1)
  {
  	fprintf(fpout,"%s ",database->pstyle.parms[o]);
  }
  fprintf(fpout,"\n");
  if(sysdat->nbnds > 0){
    fprintf(fpout,"bond_style    ");
	
	  for (o = 0; o < database->bstyle.nparms; o += 1)
	  {
	  	fprintf(fpout,"%s ",database->bstyle.parms[o]);
	  }
	  fprintf(fpout,"\n");
  }	
  if(sysdat->nangs > 0){
      fprintf(fpout,"angle_style    ");
		
	  for (o = 0; o < database->astyle.nparms; o += 1)
	  {
	  	fprintf(fpout,"%s ",database->astyle.parms[o]);
	  }
	  fprintf(fpout,"\n");
  }	
  if(sysdat->ndihes > 0){
      fprintf(fpout,"dihedral_style    ");
		
	  for (o = 0; o < database->dstyle.nparms; o += 1)
	  {
	  	fprintf(fpout,"%s ",database->dstyle.parms[o]);
	  }
	  fprintf(fpout,"\n");
  }	
  if(sysdat->total_improps>0){
    fprintf(fpout,"improper_style    ");
		
	  for (o = 0; o < database->istyle.nparms; o += 1)
	  {
	  	fprintf(fpout,"%s ",database->istyle.parms[o]);
	  }
	  fprintf(fpout,"\n");
  }	
  if(ISCHARGED){
    fprintf(fpout,"kspace_style    ");
		
	  for (o = 0; o < database->kstyle.nparms; o += 1)
	  {
	  	fprintf(fpout,"%s ",database->kstyle.parms[o]);
	  }
	  fprintf(fpout,"\n");

  }	

	fprintf(fpout,"special_bonds    ");

	for (o = 0; o < database->sbonds.nparms; o += 1)
	{
	fprintf(fpout,"%s ",database->sbonds.parms[o]);
	}
	fprintf(fpout,"\n");

  	
 

  
  fprintf(fpout,"\n");

  /* first gather the unique atom types */
  
  uniq_nats=0;
  for(i=0;i<sysdat->ntops;i++){
    for(j=0;j<topdat[i].nat;j++){
      ikeep=1;
      for(k=0;k<uniq_nats;k++){
        if(!strcmp(topdat[i].type[j],uniq_atype[k]))
        {
          ikeep=0;
          topdat[i].parm_atomtype[j]=k;
          k=uniq_nats;
        }
      }
      if(ikeep==1) {
        strcpy(uniq_atype[uniq_nats],topdat[i].type[j]);
        uniq_mass[uniq_nats]=topdat[i].mass[j];
        uniq_charge[uniq_nats]=topdat[i].charge[j];
        topdat[i].parm_atomtype[j]=uniq_nats;
        uniq_nats++;
      }
    }
  }
  
  printf("FOUND %d UNIQUE ATOMS\n",uniq_nats);
  sysdat->uniq_nats=uniq_nats;
  for(i=0;i<uniq_nats;i++)
    fprintf(fpout,"mass   %-5d  %6.4f # %s\n",i+1, uniq_mass[i],uniq_atype[i]);
  
  fprintf(fpout,"\n");
  
/* get pair interactions */
  
  
  for(i=0;i<uniq_nats;i++){
    for(j=i;j<uniq_nats;j++){
      ifound=0;
      for(k=0;k<database->nvdwtype;k++){
        if(!strcmp(database->vdwtype1[k],uniq_atype[i]) && !strcmp(database->vdwtype2[k],uniq_atype[j])) {
/*           printf("Found database entry\n"); */
          ifound=1;
          vdwtmp=k;
          k=database->nvdwtype;
        }
        else if (!strcmp(database->vdwtype2[k],uniq_atype[i]) && !strcmp(database->vdwtype1[k],uniq_atype[j])) {
/*           printf("Found database entry\n"); */
          ifound=1;
          vdwtmp=k;
          k=database->nvdwtype;
        }
      }
      if(ifound==0){
        printf("*********************\n");
        printf("ERROR:No params for VDW interaction between %s and %s\n",uniq_atype[i],uniq_atype[j]);
        printf("UPDATE DATABASE!!!\n");
        exit(1);
      } else if (ifound==1){
/*         printf("Found database entry for VDW\n"); */
        //fprintf(fpout,"pair_coeff  %-5d %-5d %-6s %5.4f %5.4f # %-4s %-4s\n",i+1,j+1,&database->vdwstyle[vdwtmp], database->eps[vdwtmp],database->sig[vdwtmp], &database->vdwtype1[vdwtmp], &database->vdwtype2[vdwtmp]);
/*	fprintf(fpout,"pair_coeff  %-5d %-5d %5.4f %5.4f # %-4s %-4s\n",i+1,j+1, database->eps[vdwtmp],database->sig[vdwtmp], &database->vdwtype1[vdwtmp], &database->vdwtype2[vdwtmp]);*/
	fprintf(fpout,"pair_coeff  %-5d %-5d",i+1,j+1);
	fprintf(fpout," ");
	int o;
	for (o = 0; o < database->p_pair[vdwtmp].nparms; o += 1)
	{
		fprintf(fpout,"%s ",database->p_pair[vdwtmp].parms[o]);
	//	printf("%s ",database->p_pair[vdwtmp].parms[o]);
	}
	fprintf(fpout,"# %-4s %-4s\n",&database->vdwtype1[vdwtmp], &database->vdwtype2[vdwtmp]);
      }
    }
  }

/* get bond interactions  */
  fprintf(fpout,"\n");

  if(sysdat->nbnds > 0) {
    uniq_bnds=0;
    
    for(i=0;i<sysdat->ntops;i++){
      for(j=0;j<topdat[i].nbnd;j++){
        ikeep=1;
        datndx=-1;
/* now compare to the database */
        for(k=0;k<database->nbndtype;k++){
          if(!strcmp(database->bndtype1[k], topdat[i].type[topdat[i].bndndx1[j]-1])
             && !strcmp(database->bndtype2[k], topdat[i].type[topdat[i].bndndx2[j]-1])){
            datndx=k;
            k=database->nbndtype;
          }
          if(!strcmp(database->bndtype2[k], topdat[i].type[topdat[i].bndndx1[j]-1])
             && !strcmp(database->bndtype1[k], topdat[i].type[topdat[i].bndndx2[j]-1])){
            datndx=k;
            k=database->nbndtype;
          }
        }
/* did not find any params to die */
        if(datndx==-1){
          printf("ERROR: DID NOT FIND BOND PARAMETERS IN DATABASE %d %d %s %s\n",topdat[i].bndndx1[j],
             topdat[i].bndndx2[j],topdat[i].type[topdat[i].bndndx1[j]-1],topdat[i].type[topdat[i].bndndx2[j]-1]);
          printf("%d %s\n",topdat[i].bndndx1[j],topdat[i].type[0]);
          
          exit(1);
        }
        
/* Now make sure we do not already know we have this interaction */
        
        for(k=0;k<uniq_bnds;k++){
          if(datndx==bnd_params[k])
          {
            ikeep=0;
            topdat[i].bndtype[j]=k;
            k=uniq_bnds;
          }
        }
        if(ikeep==1) {
          bnd_params[uniq_bnds]=datndx;
          sysdat->param_bnds[uniq_bnds]=datndx;
          topdat[i].bndtype[j]=uniq_bnds;
          uniq_bnds++;
        }
      }
    }
    for(i=0;i<uniq_bnds;i++){
    //  fprintf(fpout,"bond_coeff  %-6d %8.4f %8.4f # %s %s\n",i+1,database->fbnd[bnd_params[i]],
     //    database->bnde[bnd_params[i]],database->bndtype1[bnd_params[i]],database->bndtype2[bnd_params[i]]);
		fprintf(fpout,"bond_coeff  %-6d",i+1);
		fprintf(fpout," ");
		int o;
		for (o = 0; o < database->p_bnd[bnd_params[i]].nparms; o += 1)
		{
			fprintf(fpout,"%s ",database->p_bnd[bnd_params[i]].parms[o]);
			//printf("%s ",database->p_bnd[i].parms[o]);
		}
		fprintf(fpout,"# %-4s %-4s\n",&database->bndtype1[bnd_params[i]], &database->bndtype2[bnd_params[i]]);

    }
  }
  sysdat->uniq_nbnds=uniq_bnds;
  
/* get angle interactions if needed */
  fprintf(fpout,"\n");
  
  if(sysdat->nangs > 0) {
    
    uniq_angs=0;
    
    for(i=0;i<sysdat->ntops;i++){
      for(j=0;j<topdat[i].nang;j++){
        ikeep=1;
        datndx=-1;
        
/* now compare to the database */
        
        for(k=0;k<database->nangtype;k++){
          if(!strcmp(database->angtype2[k], topdat[i].type[topdat[i].angndx2[j]-1]))
          {
            if(!strcmp(database->angtype1[k], topdat[i].type[topdat[i].angndx1[j]-1])
               && !strcmp(database->angtype3[k], topdat[i].type[topdat[i].angndx3[j]-1])){
              datndx=k;
              k=database->nangtype;
            }
            if(!strcmp(database->angtype3[k], topdat[i].type[topdat[i].angndx1[j]-1])
               && !strcmp(database->angtype1[k], topdat[i].type[topdat[i].angndx3[j]-1])){
              datndx=k;
              k=database->nangtype;
            }
          }
        }
        

/* RHD Get the VDW for the CG angles */
/*        ifound=0;*/
/*        for(k=0;k<database->nvdwtype;k++){*/
/*          if(!strcmp(database->vdwtype1[k],topdat[i].type[topdat[i].angndx1[j]-1]) && !strcmp(database->vdwtype2[k],*/
/*                topdat[i].type[topdat[i].angndx3[j]-1])) {*/
/*            ifound=1;*/
/*            vdwtmp=k;*/
/*            k=database->nvdwtype;*/
/*          }*/
/*          else if (!strcmp(database->vdwtype1[k],topdat[i].type[topdat[i].angndx3[j]-1]) && !strcmp(database->vdwtype2[k],*/
/*                topdat[i].type[topdat[i].angndx1[j]-1])) {*/
/*            ifound=1;*/
/*            vdwtmp=k;*/
/*            k=database->nvdwtype;*/
/*          }*/
/*        }*/
/*        if(ifound==0){*/
/*          printf("*********************\n");*/
/*          printf("ERROR:No params for VDW interaction between %s and %s for angle\n",topdat[i].type[topdat[i].angndx1[j]-1],*/
/*             topdat[i].type[topdat[i].angndx3[j]-1]);*/
/*          printf("UPDATE DATABASE!!!\n");*/
/*          exit(1);*/
/*        }*/
/*        */
/* end VDW for CG angles */
        
/* No params for this interaction in the database */
        
        if(datndx==-1){
          printf("ERROR: DID NOT FIND ANGLE PARAMETERS IN DATABASE %s %s %s\n",topdat[i].type[topdat[i].angndx1[j]-1],
             topdat[i].type[topdat[i].angndx2[j]-1],topdat[i].type[topdat[i].angndx3[j]-1]);
          exit(1);
        }
        
/* Now make sure we do not already have this one */
        
        for(k=0;k<uniq_angs;k++){
          if(datndx==ang_params[k])
          {
            ikeep=0;
            topdat[i].angtype[j]=k;
            k=uniq_angs;
          }
        }
        if(ikeep==1) {
          ang_params[uniq_angs]=datndx;
         // ang_vdw[uniq_angs]=vdwtmp;
          sysdat->param_angs[uniq_angs]=datndx;
          topdat[i].angtype[j]=uniq_angs;
          uniq_angs++;
        }
      }
    }
    for(i=0;i<uniq_angs;i++){
/*     fprintf(fpout,"angle_coeff %-6d %8.4f %8.4f %8.4f %8.4f # %s %s %s\n",i+1,database->fang[ang_params[i]],*/
/*        database->ange[ang_params[i]],database->kr[ang_params[i]],database->ro[ang_params[i]],database->angtype1[ang_params[i]],database->angtype2[ang_params[i]],database->angtype3[ang_params[i]]);*/
		fprintf(fpout,"angle_coeff  %-6d",i+1);
		fprintf(fpout," ");
		int o;
		for (o = 0; o < database->p_ang[ang_params[i]].nparms; o += 1)
		{
			fprintf(fpout,"%s ",database->p_ang[ang_params[i]].parms[o]);
			//printf("%s ",database->p_ang[i].parms[o]);
		}
		fprintf(fpout,"# %-4s %-4s %-4s\n",&database->angtype1[ang_params[i]], &database->angtype2[ang_params[i]],&database->angtype3[ang_params[i]]);

    }
  }
  sysdat->uniq_nangs=uniq_angs;

  /* get dihedral angle interactions if needed */
  fprintf(fpout,"\n");
  
  if(sysdat->ndihes > 0) {
    
    uniq_dihes=0;
    
    for(i=0;i<sysdat->ntops;i++){
      for(j=0;j<topdat[i].ndihe;j++){
        ikeep=1;
        datndx=-1;
        
/* now compare to the database */
        
        for(k=0;k<database->ndihetype;k++){
          if( (!strcmp(database->dihetype2[k], topdat[i].type[topdat[i].dihendx2[j]-1]) && !strcmp(database->dihetype3[k], topdat[i].type[topdat[i].dihendx3[j]-1])) || (!strcmp(database->dihetype3[k], topdat[i].type[topdat[i].dihendx2[j]-1]) && !strcmp(database->dihetype2[k], topdat[i].type[topdat[i].dihendx3[j]-1])) )
          {
            if(!strcmp(database->dihetype1[k], topdat[i].type[topdat[i].dihendx1[j]-1])
               && !strcmp(database->dihetype4[k], topdat[i].type[topdat[i].dihendx4[j]-1])){
              datndx=k;
              k=database->ndihetype;
            }
            if(!strcmp(database->dihetype4[k], topdat[i].type[topdat[i].dihendx1[j]-1])
               && !strcmp(database->dihetype1[k], topdat[i].type[topdat[i].dihendx4[j]-1])){
              datndx=k;
              k=database->ndihetype;
            }
          }
        }
             
/* No params for this interaction in the database */
        
        if(datndx==-1){
          printf("ERROR: DID NOT FIND DIHEDRAL PARAMETERS IN DATABASE %s %s %s %s\n",topdat[i].type[topdat[i].dihendx1[j]-1],
             topdat[i].type[topdat[i].dihendx2[j]-1],topdat[i].type[topdat[i].dihendx3[j]-1],topdat[i].type[topdat[i].dihendx4[j]-1]);
          exit(1);
        }
        
/* Now make sure we do not already have this one */
        
        for(k=0;k<uniq_dihes;k++){
          if(datndx==dihe_params[k])
          {
            ikeep=0;
            topdat[i].dihetype[j]=k;
            k=uniq_dihes;
          }
        }
        if(ikeep==1) {
          dihe_params[uniq_dihes]=datndx;
          sysdat->param_dihes[uniq_dihes]=datndx;
          topdat[i].dihetype[j]=uniq_dihes;
          uniq_dihes++;
        }
      }
    }
    for(i=0;i<uniq_dihes;i++){
/*     fprintf(fpout,"dihedral_coeff %-6d %8.4f %-2d %-2d 0.0 # %s %s %s %s\n",i+1,database->k_dihe[dihe_params[i]],database->d_dihe[dihe_params[i]],database->n_dihe[dihe_params[i]],*/
/*        database->dihetype1[dihe_params[i]],database->dihetype2[dihe_params[i]],database->dihetype3[dihe_params[i]],database->dihetype4[dihe_params[i]]);*/
		fprintf(fpout,"dihedral_coeff  %-6d",i+1);
		fprintf(fpout," ");
		int o;
		for (o = 0; o < database->p_dihe[dihe_params[i]].nparms; o += 1)
		{
			fprintf(fpout,"%s ",database->p_dihe[dihe_params[i]].parms[o]);
			//printf("%s ",database->p_dihe[i].parms[o]);
		}
		fprintf(fpout,"# %-4s %-4s %-4s %-4s\n",&database->dihetype1[dihe_params[i]], &database->dihetype2[dihe_params[i]],&database->dihetype3[dihe_params[i]],&database->dihetype4[dihe_params[i]]);
    }
  }
  sysdat->uniq_ndihes=uniq_dihes;


 /* get improper angle interactions if needed */
  fprintf(fpout,"\n");
  
  if(sysdat->nimp > 0) {
    
    uniq_imps=0;
    
    for(i=0;i<sysdat->ntops;i++){
      for(j=0;j<topdat[i].nimp;j++){
        ikeep=1;
        datndx=-1;
        
/* now compare to the database */
        
        for(k=0;k<database->nimptype;k++){
          if( (!strcmp(database->imptype2[k], topdat[i].type[topdat[i].impropndx2[j]-1]) && !strcmp(database->imptype3[k], topdat[i].type[topdat[i].impropndx3[j]-1])) || (!strcmp(database->imptype3[k], topdat[i].type[topdat[i].impropndx2[j]-1]) && !strcmp(database->imptype2[k], topdat[i].type[topdat[i].impropndx3[j]-1])) )
          {
            if(!strcmp(database->imptype1[k], topdat[i].type[topdat[i].impropndx1[j]-1])
               && !strcmp(database->imptype4[k], topdat[i].type[topdat[i].impropndx4[j]-1])){
              datndx=k;
              k=database->nimptype;
            }
            if(!strcmp(database->imptype4[k], topdat[i].type[topdat[i].impropndx1[j]-1])
               && !strcmp(database->imptype1[k], topdat[i].type[topdat[i].impropndx4[j]-1])){
              datndx=k;
              k=database->nimptype;
            }
          }
        }
             
/* No params for this interaction in the database */
        
        if(datndx==-1){
          printf("ERROR: DID NOT FIND IMPROPER PARAMETERS IN DATABASE %s %s %s %s\n",topdat[i].type[topdat[i].impropndx1[j]-1],
             topdat[i].type[topdat[i].impropndx2[j]-1],topdat[i].type[topdat[i].impropndx3[j]-1],topdat[i].type[topdat[i].impropndx4[j]-1]);
          exit(1);
        }
        
/* Now make sure we do not already have this one */
        
        for(k=0;k<uniq_imps;k++){
          if(datndx==imp_params[k])
          {
            ikeep=0;
            topdat[i].imptype[j]=k;
            k=uniq_imps;
          }
        }
        if(ikeep==1) {
          imp_params[uniq_imps]=datndx;
          sysdat->param_imps[uniq_imps]=datndx;
          topdat[i].imptype[j]=uniq_imps;
          uniq_imps++;
        }
      }
    }
    for(i=0;i<uniq_imps;i++){
/*     fprintf(fpout,"improper_coeff %-6d %8.4f %8.4f # %s %s %s %s\n",i+1,database->k_imp[imp_params[i]],database->d_imp[imp_params[i]],*/
/*        database->imptype1[imp_params[i]],database->imptype2[imp_params[i]],database->imptype3[imp_params[i]],database->imptype4[imp_params[i]]);*/

		fprintf(fpout,"improper_coeff  %-6d",i+1);
		fprintf(fpout," ");
		int o;
		for (o = 0; o < database->p_imp[imp_params[i]].nparms; o += 1)
		{
			fprintf(fpout,"%s ",database->p_imp[imp_params[i]].parms[o]);
			//printf("%s ",database->p_imp[i].parms[o]);
		}
		fprintf(fpout,"# %-4s %-4s %-4s %-4s\n",&database->imptype1[imp_params[i]], &database->imptype2[imp_params[i]],&database->imptype3[imp_params[i]],&database->imptype4[imp_params[i]]);
    }
  }
  sysdat->uniq_imps=uniq_imps;

return;

}


/* Read the database file and store unique params */
/* Warn if you find duplicates */


void read_database(char *filename, DATABASE *database) 
{
  printf("READING DATABASE: %s\n",filename);	
  FILE *fpin;
  int c,nvdw,nbnd,nang,i,ikeep,nimp;
  char col1[10],vdwtype1[10],vdwtype2[10];
  char bndtype1[5],bndtype2[5],angtype1[5],angtype2[5],angtype3[5];
   //number of paramters to read
  int np_pair, np_bnd, np_ang, np_dihe, np_imp;
  //set the default values -- count matches CHARMM
  np_pair=2,np_bnd=2,np_ang=4,np_dihe=2,np_imp=2;
   
  //Temporary containers for params
  PARAMS p_pair, p_bnd, p_ang, p_dihe, p_imp;		  

  char line[350];

  char dihetype1[5],dihetype2[5],dihetype3[5],dihetype4[5];
   char imptype1[5],imptype2[5],imptype3[5],imptype4[5];
 
  int ndihe,d_dihe,n_dihe;
 
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  nvdw=0;
  nbnd=0;
  nang=0;
  ndihe=0;
  nimp=0;

 // printf("reading...\n");	
  while(c=fscanf(fpin,"%s",&col1)!= EOF)
  {
	//check for style commands
   if (!strcmp((char *)&col1, "pstyle"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#')){
		
			 int nwords;
			 char *words[PBUFF];
			 char hash = '#';
			 nwords = getwordstochar(line, words, MAXPARM,hash);
			//copy params to temp param holder
			database->pstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->pstyle, words[i],i);
			}
		  }	
     	 else{
			int nwords;
			 char *words[PBUFF];
			 nwords = getwords(line, words, MAXPARM);
			//copy params to temp param holder
			database->pstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->pstyle, words[i],i);
			}
		}
				
	}

	if (!strcmp((char *)&col1, "bstyle"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#')){
		
			 int nwords;
			 char *words[PBUFF];
			 char hash = '#';
			 nwords = getwordstochar(line, words, MAXPARM,hash);
			//copy params to temp param holder
			database->bstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->bstyle, words[i],i);
			}
		  }	
     	 else{
			int nwords;
			 char *words[PBUFF];
			 nwords = getwords(line, words, MAXPARM);
			//copy params to temp param holder
			database->bstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->bstyle, words[i],i);
			}
		}
				
	}


	if (!strcmp((char *)&col1, "astyle"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#')){
		
			 int nwords;
			 char *words[PBUFF];
			 char hash = '#';
			 nwords = getwordstochar(line, words, MAXPARM,hash);
			//copy params to temp param holder
			database->astyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->astyle, words[i],i);
			}
		  }	
     	 else{
			int nwords;
			 char *words[PBUFF];
			 nwords = getwords(line, words, MAXPARM);
			//copy params to temp param holder
			database->astyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->astyle, words[i],i);
			}
		}
				
	}


	if (!strcmp((char *)&col1, "dstyle"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#')){
		
			 int nwords;
			 char *words[PBUFF];
			 char hash = '#';
			 nwords = getwordstochar(line, words, MAXPARM,hash);
			//copy params to temp param holder
			database->dstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->dstyle, words[i],i);
			}
		  }	
     	 else{
			int nwords;
			 char *words[PBUFF];
			 nwords = getwords(line, words, MAXPARM);
			//copy params to temp param holder
			database->dstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->dstyle, words[i],i);
			}
		}
				
	}

	if (!strcmp((char *)&col1, "istyle"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#')){
		
			 int nwords;
			 char *words[PBUFF];
			 char hash = '#';
			 nwords = getwordstochar(line, words, MAXPARM,hash);
			//copy params to temp param holder
			database->istyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->istyle, words[i],i);
			}
		  }	
     	 else{
			int nwords;
			 char *words[PBUFF];
			 nwords = getwords(line, words, MAXPARM);
			//copy params to temp param holder
			database->istyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->istyle, words[i],i);
			}
		}
				
	}

	if (!strcmp((char *)&col1, "kstyle"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#')){
		
			 int nwords;
			 char *words[PBUFF];
			 char hash = '#';
			 nwords = getwordstochar(line, words, MAXPARM,hash);
			//copy params to temp param holder
			database->kstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->kstyle, words[i],i);
			}
		  }	
     	 else{
			int nwords;
			 char *words[PBUFF];
			 nwords = getwords(line, words, MAXPARM);
			//copy params to temp param holder
			database->kstyle.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->kstyle, words[i],i);
			}
		}
				
	}

	if (!strcmp((char *)&col1, "sbonds"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#')){
		
			 int nwords;
			 char *words[PBUFF];
			 char hash = '#';
			 nwords = getwordstochar(line, words, MAXPARM,hash);
			//copy params to temp param holder
			database->sbonds.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->sbonds, words[i],i);
			}
		  }	
     	 else{
			int nwords;
			 char *words[PBUFF];
			 nwords = getwords(line, words, MAXPARM);
			//copy params to temp param holder
			database->sbonds.nparms = nwords;
			int i;
			for (i = 0; i < nwords;i += 1)
			{
				push_toPARAM(&database->sbonds, words[i],i);
			}
		}
				
	}


	//check for paramter number update commands in database
	if (!strcmp((char *)&col1, "npairargs"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#'))
			sscanf(line,"%d %*[^\n]",&np_pair);	 
		else
			sscanf(line,"%d",&np_pair);
				
	}
	if (!strcmp((char *)&col1, "nbondargs"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#'))
			sscanf(line,"%d %*[^\n]",&np_bnd);	 
		else
			sscanf(line,"%d",&np_bnd); 		
	}
	if (!strcmp((char *)&col1, "nangleargs"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#'))
			sscanf(line,"%d %*[^\n]",&np_ang);	 
		else
			sscanf(line,"%d",&np_ang);	 		
	}
	if (!strcmp((char *)&col1, "ndihedralargs"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#'))
			sscanf(line,"%d %*[^\n]",&np_dihe);	 
		else
			sscanf(line,"%d",&np_dihe);	 		
	}
	if (!strcmp((char *)&col1, "nimproperargs"))
	{
		fscanf(fpin,"%[^\n]",line);
		if(strchr(line,'#'))
			sscanf(line,"%d %*[^\n]",&np_imp);	 
		else
			sscanf(line,"%d",&np_imp);	 		
	}
  //  printf("column 1 %s \n",col1);
    if(!strcmp((char *)&col1, "pair"))
    {

	 // printf("found pair!\n");
		fscanf(fpin,"%[^\n]",line);
	 	
      if(strchr(line,'#')){
		//printf("reading pair line with comment hash\n");	
       // sscanf(line,"%s %s %lf %lf %*[^\n]",&vdwtype1,&vdwtype2,&eps,&sig);
		 int nwords;
		 char *words[PBUFF];
		 char hash = '#';
		 nwords = getwordstochar(line, words, MAXPARM+2,hash);
		  strcpy(vdwtype1,words[0]);	
		   strcpy(vdwtype2,words[1]);	
         if (nwords-2>np_pair){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR PAIR %s %s\n",&vdwtype1,&vdwtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_pair,nwords-2);
            exit(0);
         }
		 else if(nwords-2<np_pair){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR PAIR %s %s\n",&vdwtype1,&vdwtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_pair,nwords-2);
            exit(0);
		}
		//copy params to temp param holder
		p_pair.nparms = np_pair;
		int i;
		for (i = 0; i < np_pair;i += 1)
		{
			push_toPARAM(&p_pair, words[i+2],i);
		}
	  }	
      else{
		// printf("reading pair line\n");	
        //sscanf(line,"%s %s %lf %lf",&vdwtype1,&vdwtype2,&eps,&sig);
	     int nwords;
		 char *words[PBUFF];
		 nwords = getwords(line, words, MAXPARM+2);
		// printf("Number of words is %d \n",nwords);	
		  strcpy(vdwtype1,words[0]);	
		   strcpy(vdwtype2,words[1]);	
		// printf("checking against number of pair parms\n");	
         if (nwords-2>np_pair){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR PAIR %s %s\n",&vdwtype1,&vdwtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_pair,nwords-2);
            exit(0);
         }
		 else if(nwords-2<np_pair){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR PAIR %s %s\n",&vdwtype1,&vdwtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_pair,nwords-2);
            exit(0);
		}
		//printf("passed number parms check\n");
		//copy params to temp param holder
		p_pair.nparms = np_pair;
		//printf("set the p_pair nparms\n");
		//printf("p_pair.nparms %d \n",p_pair.nparms);
		int i;
		for (i = 0; i < np_pair;i += 1)
		{
			push_toPARAM(&p_pair, words[i+2],i);
		}
		//printf("pushed all params to p_pair\n");
	  }
	 // printf("eps sig %f %f",eps,sig);
      ikeep=1;
      for(i=0;i<nvdw;i++){
        if(!strcmp(vdwtype1,database->vdwtype1[i]) && !strcmp(vdwtype2,database->vdwtype2[i]))
        {
          printf("WARNING: FOUND DUP PAIR PARAM %s %s\n",&vdwtype1,&vdwtype2);
          ikeep=0;
        }
        else if (!strcmp(vdwtype1,database->vdwtype2[i]) && !strcmp(vdwtype2,database->vdwtype1[i]))
        {
          printf("WARNING: FOUND DUP PAIR PARAM %s %s\n",&vdwtype1,&vdwtype2);
          ikeep=0;
        }
      }
      if(ikeep==1)
      {
        strcpy(database->vdwtype1[nvdw],vdwtype1);
        strcpy(database->vdwtype2[nvdw],vdwtype2);
       // strcpy(database->vdwstyle[nvdw],vdwstyle);
/*        database->eps[nvdw]=eps;*/
/*        database->sig[nvdw]=sig;*/
		//printf("attempting copy from p_pair to database ppair\n");
		copySinglePARAM_toPARAMstarAt(&p_pair, database->p_pair, nvdw);
		//printf(" copy successful\n");
        nvdw++;
      }
    }

    if(!strcmp((char *)&col1, "bond"))
    {

      fscanf(fpin,"%[^\n]",line);
      
      if(strchr(line,'#')){
 //     sscanf(line,"%s %s %lf %lf %*[^\n]",&bndtype1,&bndtype2,&fbnd,&bnde);
		 int nwords;
		 char *words[PBUFF];
		 char hash = '#';
		 nwords = getwordstochar(line, words, MAXPARM+2,hash);
		  strcpy(bndtype1,words[0]);	
		   strcpy(bndtype2,words[1]);	
         if (nwords-2>np_bnd){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR BOND %s %s\n",&bndtype1,&bndtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_bnd,nwords-2);
            exit(0);
         }
		 else if(nwords-2<np_bnd){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR BOND %s %s\n",&bndtype1,&bndtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_bnd,nwords-2);
            exit(0);
		}
		//copy params to temp param holder
		p_bnd.nparms = np_bnd;
		int i;
		for (i = 0; i < np_bnd;i += 1)
		{
			push_toPARAM(&p_bnd, words[i+2],i);
		}
	   }	
      else{
       // sscanf(line,"%s %s %lf %lf",&bndtype1,&bndtype2,&fbnd,&bnde);
		int nwords;
		 char *words[PBUFF];
		 nwords = getwords(line, words, MAXPARM+2);
		// printf("Number of words is %d \n",nwords);	
		  strcpy(bndtype1,words[0]);	
		   strcpy(bndtype2,words[1]);	
		// printf("checking against number of pair parms\n");	
         if (nwords-2>np_bnd){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR BOND %s %s\n",&bndtype1,&bndtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_bnd,nwords-2);
            exit(0);
         }
		 else if(nwords-2<np_bnd){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR BOND %s %s\n",&bndtype1,&bndtype2);
			printf("--EXPECTED %d BUT GOT %d\n",np_bnd,nwords-2);
            exit(0);
		}
		//printf("passed number parms check\n");
		//copy params to temp param holder
		p_bnd.nparms = np_bnd;
		
		int i;
		for (i = 0; i < np_bnd;i += 1)
		{
			push_toPARAM(&p_bnd, words[i+2],i);
		}
		//printf("pushed all params to p_pair\n");
	  }
      ikeep=1;
      for(i=0;i<nbnd;i++){
        if(!strcmp(bndtype1,database->bndtype1[i]) && !strcmp(bndtype2,database->bndtype2[i]))
        {
          printf("WARNGIN: FOUND DUP BOND PARAM %s %s\n",&bndtype1,&bndtype2);
          ikeep=0;
        }
        else if (!strcmp(bndtype1,database->bndtype2[i]) && !strcmp(bndtype2,database->bndtype1[i]))
        {
          printf("WARNING: FOUND DUP BOND PARAM %s %s\n",&bndtype1,&bndtype2);
          ikeep=0;
        }
      }
      if(ikeep==1)
      {
        strcpy(database->bndtype1[nbnd],bndtype1);
        strcpy(database->bndtype2[nbnd],bndtype2);
/*        database->fbnd[nbnd]=fbnd;*/
/*        database->bnde[nbnd]=bnde;*/
		copySinglePARAM_toPARAMstarAt(&p_bnd, database->p_bnd, nbnd);
        nbnd++;
      }
    }

    if(!strcmp((char *)&col1, "angle"))
    {

      fscanf(fpin,"%[^\n]",line);

      if(strchr(line,'#')){
 //     sscanf(line,"%s %s %lf %lf %*[^\n]",&bndtype1,&bndtype2,&fbnd,&bnde);
		 int nwords;
		 char *words[PBUFF];
		 char hash = '#';
		 nwords = getwordstochar(line, words, MAXPARM+3,hash);
		  strcpy(angtype1,words[0]);	
		   strcpy(angtype2,words[1]);
			strcpy(angtype3,words[2]);	
         if (nwords-3>np_ang){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR ANGLE %s %s %s\n",&angtype1,&angtype2,&angtype3);
			printf("--EXPECTED %d BUT GOT %d\n",np_ang,nwords-3);
            exit(0);
         }
		 else if(nwords-3<np_ang){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR ANGLE %s %s %s\n",&angtype1,&angtype2,&angtype3);
			printf("--EXPECTED %d BUT GOT %d\n",np_ang,nwords-3);
            exit(0);
		 }
		//copy params to temp param holder
		p_ang.nparms = np_ang;
		int i;
		for (i = 0; i < np_ang;i += 1)
		{
			push_toPARAM(&p_ang, words[i+3],i);
		}
	   }	
      else{
       // sscanf(line,"%s %s %lf %lf",&bndtype1,&bndtype2,&fbnd,&bnde);
		int nwords;
		 char *words[PBUFF];
		 nwords = getwords(line, words, MAXPARM+3);
		 strcpy(angtype1,words[0]);	
		   strcpy(angtype2,words[1]);
			strcpy(angtype3,words[2]);	
         if (nwords-3>np_ang){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR ANGLE %s %s %s\n",&angtype1,&angtype2,&angtype3);
			printf("--EXPECTED %d BUT GOT %d\n",np_ang,nwords-3);
            exit(0);
         }
		 else if(nwords-3<np_ang){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR ANGLE %s %s %s\n",&angtype1,&angtype2,&angtype3);
			printf("--EXPECTED %d BUT GOT %d\n",np_ang,nwords-3);
            exit(0);
		 }
		//copy params to temp param holder
		p_ang.nparms = np_ang;
		int i;
		for (i = 0; i < np_ang;i += 1)
		{
			push_toPARAM(&p_ang, words[i+3],i);
		}
	  }

      ikeep=1;
      for(i=0;i<nang;i++){
        if(!strcmp(angtype2,database->angtype2[i])) {
          if(!strcmp(angtype1,database->angtype1[i]) && !strcmp(angtype3,database->angtype3[i]))
          {
            printf("WARNGIN: FOUND DUP ANGLE PARAM %s %s\n",&angtype1,&angtype2,&angtype3);
            ikeep=0;
          }
          else if (!strcmp(angtype3,database->angtype1[i]) && !strcmp(angtype1,database->angtype3[i]))
          {
            printf("WARNING: FOUND DUP ANGLE PARAM %s %s\n",&angtype1,&angtype2,&angtype3);
            ikeep=0;
          }
        }
      }
      if(ikeep==1)
      {
        strcpy(database->angtype1[nang],angtype1);
        strcpy(database->angtype2[nang],angtype2);
        strcpy(database->angtype3[nang],angtype3);
/*        database->fang[nang]=fang;*/
/*        database->ange[nang]=ange;*/
/*		 database->kr[nang]=kr;*/
/*		 database->ro[nang]=ro;*/
		copySinglePARAM_toPARAMstarAt(&p_ang, database->p_ang, nang);
        nang++;
      }
    }

    if(!strcmp((char *)&col1, "dihedral"))
    {

      fscanf(fpin,"%[^\n]",line);

      if(strchr(line,'#')){
		 int nwords;
		 char *words[PBUFF];
		 char hash = '#';
		 nwords = getwordstochar(line, words, MAXPARM+4,hash);
		  strcpy(dihetype1,words[0]);	
		   strcpy(dihetype2,words[1]);
			strcpy(dihetype3,words[2]);
			strcpy(dihetype4,words[3]);	
         if (nwords-4>np_dihe){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR DIHEDRAL %s %s %s %s\n",&dihetype1,&dihetype2,&dihetype3,&dihetype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_dihe,nwords-4);
            exit(0);
         }
		 else if(nwords-4<np_dihe){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR DIHEDRAL %s %s %s %s\n",&dihetype1,&dihetype2,&dihetype3,&dihetype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_dihe,nwords-4);
            exit(0);
		 }
		//copy params to temp param holder
		p_dihe.nparms = np_dihe;
		int i;
		for (i = 0; i < np_dihe;i += 1)
		{
			push_toPARAM(&p_dihe, words[i+4],i);
		}
	   }	
      else{
       // sscanf(line,"%s %s %lf %lf",&bndtype1,&bndtype2,&fbnd,&bnde);
		int nwords;
		 char *words[PBUFF];
		 nwords = getwords(line, words, MAXPARM+3);
		  strcpy(dihetype1,words[0]);	
		   strcpy(dihetype2,words[1]);
			strcpy(dihetype3,words[2]);
			strcpy(dihetype4,words[3]);	
         if (nwords-4>np_dihe){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR DIHEDRAL %s %s %s %s\n",&dihetype1,&dihetype2,&dihetype3,&dihetype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_dihe,nwords-4);
            exit(0);
         }
		 else if(nwords-4<np_dihe){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR DIHEDRAL %s %s %s %s\n",&dihetype1,&dihetype2,&dihetype3,&dihetype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_dihe,nwords-4);
            exit(0);
		 }
		//copy params to temp param holder
		p_dihe.nparms = np_dihe;
		int i;
		for (i = 0; i < np_dihe;i += 1)
		{
			push_toPARAM(&p_dihe, words[i+4],i);
		}
	  }

      ikeep=1;
      for(i=0;i<ndihe;i++){
        if( (!strcmp(dihetype2,database->dihetype2[i]) && !strcmp(dihetype3,database->dihetype3[i])) || (!strcmp(dihetype3,database->dihetype2[i]) && !strcmp(dihetype2,database->dihetype3[i])) ) {
          if(!strcmp(dihetype1,database->dihetype1[i]) && !strcmp(dihetype4,database->dihetype4[i]))
          {
            printf("WARNING: FOUND DUP DIHEDRAL PARAM %s %s %s %s\n",&dihetype1,&dihetype2,&dihetype3,&dihetype4);
            ikeep=0;
          }
          else if (!strcmp(dihetype4,database->dihetype1[i]) && !strcmp(dihetype1,database->dihetype4[i]))
          {
            printf("WARNING: FOUND DUP DIHEDRAL PARAM %s %s %s %s\n",&dihetype1,&dihetype2,&dihetype3,&dihetype4);
            ikeep=0;
          }
        }
      }
      if(ikeep==1)
      {
        strcpy(database->dihetype1[ndihe],dihetype1);
        strcpy(database->dihetype2[ndihe],dihetype2);
        strcpy(database->dihetype3[ndihe],dihetype3);
        strcpy(database->dihetype4[ndihe],dihetype4);
    
		copySinglePARAM_toPARAMstarAt(&p_dihe, database->p_dihe, ndihe);
        ndihe++;
      }
    }
	
	 if(!strcmp((char *)&col1, "improper"))
    {

      fscanf(fpin,"%[^\n]",line);

      if(strchr(line,'#')){
         int nwords;
		 char *words[PBUFF];
		 char hash = '#';
		 nwords = getwordstochar(line, words, MAXPARM+4,hash);
		  strcpy(imptype1,words[0]);	
		   strcpy(imptype2,words[1]);
			strcpy(imptype3,words[2]);
			strcpy(imptype4,words[3]);	
         if (nwords-4>np_imp){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR IMPROPER %s %s %s %s\n",&imptype1,&imptype2,&imptype3,&imptype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_imp,nwords-4);
            exit(0);
         }
		 else if(nwords-4<np_imp){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR IMPROPER %s %s %s %s\n",&imptype1,&imptype2,&imptype3,&imptype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_imp,nwords-4);
            exit(0);
		 }
		//copy params to temp param holder
		p_imp.nparms = np_imp;
		int i;
		for (i = 0; i < np_imp;i += 1)
		{
			push_toPARAM(&p_imp, words[i+4],i);
		}
	   }	
      else{
       // sscanf(line,"%s %s %lf %lf",&bndtype1,&bndtype2,&fbnd,&bnde);
		int nwords;
		 char *words[PBUFF];
		 nwords = getwords(line, words, MAXPARM+3);
		    strcpy(imptype1,words[0]);	
		   strcpy(imptype2,words[1]);
			strcpy(imptype3,words[2]);
			strcpy(imptype4,words[3]);	
         if (nwords-4>np_imp){
         	fprintf(stderr,"ERROR: TOO MANY PARAMS FOR IMPROPER %s %s %s %s\n",&imptype1,&imptype2,&imptype3,&imptype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_imp,nwords-4);
            exit(0);
         }
		 else if(nwords-4<np_imp){
			fprintf(stderr,"ERROR: TOO FEW PARAMS FOR IMPROPER %s %s %s %s\n",&imptype1,&imptype2,&imptype3,&imptype4);
			printf("--EXPECTED %d BUT GOT %d\n",np_imp,nwords-4);
            exit(0);
		 }
		//copy params to temp param holder
		p_imp.nparms = np_imp;
		int i;
		for (i = 0; i < np_imp;i += 1)
		{
			push_toPARAM(&p_imp, words[i+4],i);
		}
	  }

      ikeep=1;
      for(i=0;i<nimp;i++){
        if( (!strcmp(imptype2,database->imptype2[i]) && !strcmp(imptype3,database->imptype3[i])) || (!strcmp(imptype3,database->imptype2[i]) && !strcmp(imptype2,database->imptype3[i])) ) {
          if(!strcmp(imptype1,database->imptype1[i]) && !strcmp(imptype4,database->imptype4[i]))
          {
            printf("WARNGIN: FOUND DUP IMPROPER PARAM %s %s %s %s\n",&imptype1,&imptype2,&imptype3,&imptype4);
            ikeep=0;
          }
          else if (!strcmp(imptype4,database->imptype1[i]) && !strcmp(imptype1,database->imptype4[i]))
          {
            printf("WARNING: FOUND DUP IMPROPER PARAM %s %s %s %s\n",&imptype1,&imptype2,&imptype3,&imptype4);
            ikeep=0;
          }
        }
      }
      if(ikeep==1)
      {
        strcpy(database->imptype1[nimp],imptype1);
        strcpy(database->imptype2[nimp],imptype2);
        strcpy(database->imptype3[nimp],imptype3);
        strcpy(database->imptype4[nimp],imptype4);
        copySinglePARAM_toPARAMstarAt(&p_imp, database->p_imp, nimp);
	   
        nimp++;
      }
    }
  	

  }


     
  fclose(fpin);
  
  database->nvdwtype=nvdw;
  database->nbndtype=nbnd;
  database->nangtype=nang;
  database->ndihetype=ndihe;
  database->nimptype=nimp;

return;
}


/* Count the number of params in the database so we can allocate for storage */


void count_params(char *filename, DATABASE *database)
{
  FILE *fpin;
  int c;
  char col1[10];
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  database->nvdwtype=0;
  database->nbndtype=0;
  database->nangtype=0;
  database->ndihetype=0;
  
  while( (c=fscanf(fpin,"%s",&col1))!= EOF)
  {
    if(!strcmp((char *)&col1, "pair"))
    {

      fscanf(fpin,"%*[^\n]");
      database->nvdwtype++;
    }
    if(!strcmp((char *)&col1, "bond"))
    {

      fscanf(fpin,"%*[^\n]");
      database->nbndtype++;
    }
    if(!strcmp((char *)&col1, "angle"))
    {
      
      fscanf(fpin,"%*[^\n]");
      database->nangtype++;
    }
    if(!strcmp((char *)&col1, "dihedral"))
    {
      
      fscanf(fpin,"%*[^\n]");
      database->ndihetype++;
    }
	 if(!strcmp((char *)&col1, "improper"))
    {
      
      fscanf(fpin,"%*[^\n]");
      database->nimptype++;
    }
  }
  fclose(fpin);

return;
}


/* Read the topology file and store the data */


void read_top(char *filename, TOPDAT *topdat, int ntop, int *ISCHARGED)
{
  FILE *fpin;
  int c,ndx,bndx,andx,indx,dihendx;
  char col1[10];
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  ndx=0;
  bndx=0;
  andx=0;
  indx=0;
  dihendx=0;
  printf("READ %s\n", filename);
  
  while(c=fscanf(fpin,"%s",&col1)!= EOF)
  {
    if(!strcmp((char *)&col1, "atom"))
    {
      fscanf(fpin,"%d %s %*s %s %lf %lf %s",&topdat[ntop].index[ndx],&topdat[ntop].resname[ndx], &topdat[ntop].type[ndx],&topdat[ntop].mass[ndx],
         &topdat[ntop].charge[ndx],&topdat[ntop].segid[ndx]);
      if((topdat[ntop].charge[ndx]*topdat[ntop].charge[ndx]) > 0.00001){
        printf("CHARGE IN TOP FILE: %s %lf\n", filename, topdat[ntop].charge[ndx]);
        *ISCHARGED=1;
      }
      
        ndx++;
    }
    if(!strcmp((char *)&col1, "bond"))
    {
      fscanf(fpin,"%d %d", &topdat[ntop].bndndx1[bndx], &topdat[ntop].bndndx2[bndx]);
      bndx++;
    }
    if(!strcmp((char *)&col1, "angle"))
    {
      fscanf(fpin,"%d %d %d", &topdat[ntop].angndx1[andx], &topdat[ntop].angndx2[andx], &topdat[ntop].angndx3[andx]);
      andx++;
    }
    if(!strcmp((char *)&col1, "dihedral"))
    {
      fscanf(fpin,"%d %d %d %d", &topdat[ntop].dihendx1[dihendx], &topdat[ntop].dihendx2[dihendx], &topdat[ntop].dihendx3[dihendx], &topdat[ntop].dihendx4[dihendx]);
      dihendx++;
    }
    if(!strcmp((char *)&col1, "improper"))
    {
      fscanf(fpin,"%d %d %d %d", &topdat[ntop].impropndx1[indx], &topdat[ntop].impropndx2[indx], &topdat[ntop].impropndx3[indx],&topdat[ntop].impropndx4[indx]);
      indx++;
    }
  }
  fclose(fpin);

return;
}

/* count the number of things in the topology files so we can allocate */

void count_atoms(char *filename, TOPDAT *topdat, int ntop)
{
  FILE *fpin;
  int c;
  char col1[10];
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  topdat[ntop].nat=0;
  topdat[ntop].nbnd=0;
  topdat[ntop].nang=0;
  topdat[ntop].ndihe=0;
  topdat[ntop].nimprop=0;
  
  while( (c=fscanf(fpin,"%s",&col1))!= EOF )
  {
    if(!strcmp((char *)&col1, "atom"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nat++;
    }
    if(!strcmp((char *)&col1, "bond"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nbnd++;
    }
    if(!strcmp((char *)&col1, "angle"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nang++;
    }
    if(!strcmp((char *)&col1, "dihedral"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].ndihe++;
    }
    if(!strcmp((char *)&col1, "improper"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nimprop++;
    }
  }
  topdat[ntop].nimp = topdat[ntop].nimprop;	
  fclose(fpin);
}


/* Main routine. Call and allocate                                       */
/* The idea is to read in the topologies and then check the database for */
/* all of the required interaction params.                               */


main(int argc, char **argv )
{
  int ntops,rdtp,i,ISCHARGED;
  TOPDAT topdat[10];
  DATABASE database;
  SYSDAT sysdat;
  
  if(argc < 6 || argc%2!=0){
    printf("USAGE: setup_lammps_v2 <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <param file> <coordfile type> <coordfile>\n");
 //   printf("Prints out input files for a lammps run.  Takes pdb, xyz, or general coordinate file as the coordfile\n");
	printf("Prints out input files for a lammps run.  Takes pdb, xyz, or general coordinate (gcd) file as the coordfile\n");
    exit(1);
  }

  ISCHARGED=0;
  
//   printf ("argc=%d\n",argc); 
/*   printf("0= %s\n",argv[0]); */
/*   printf("1= %s\n",argv[1]); */
/*   printf("2= %s\n",argv[2]); */
/*   printf("3= %s\n",argv[3]); */
/*   printf("4= %s\n",argv[4]); */
//  printf("4= %s\n",argv[4]);   
//  exit(0);
  ntops=(argc-4)/2;
  printf("WILL READ %d TOPOLOGY FILE(S).\n",ntops);
  rdtp=0;

/* Loop through the topologies and count the number of atoms, bonds and bends */
  
  while(rdtp < ntops){
    topdat[rdtp].nmol=atoi(argv[(2*rdtp)+2]);
    count_atoms(argv[(2*rdtp)+1],topdat,rdtp);
    rdtp++;
  }
  sysdat.ntops=ntops;
  sysdat.nats=0;
  sysdat.nbnds=0;
  sysdat.nangs=0;
  sysdat.ndihes=0;
  sysdat.total_ats=0;
  sysdat.total_bnds=0;
  sysdat.total_angs=0;
  sysdat.total_dihes=0;
  sysdat.total_improps=0;
  
  for(i=0;i<ntops;i++){
    printf("BOOKKEEPING:\n");
    printf("TOPFILE %s\n",argv[(2*i)+1]);
    printf("FOUND: %d atoms\n",topdat[i].nat);
    printf("FOUND: %d bonds\n",topdat[i].nbnd);
    printf("FOUND: %d angles\n",topdat[i].nang);
    printf("FOUND: %d dihedrals\n",topdat[i].ndihe);
    printf("FOUND: %d improperss\n",topdat[i].nimprop);
    
    topdat[i].index   = (int *)calloc(topdat[i].nat,sizeof(int));
    topdat[i].mass    = (double *)calloc(topdat[i].nat,sizeof(double));
    topdat[i].charge  = (double *)calloc(topdat[i].nat,sizeof(double));

    topdat[i].bndndx1 = (int *)calloc(topdat[i].nbnd,sizeof(int));
    topdat[i].bndndx2 = (int *)calloc(topdat[i].nbnd,sizeof(int));
    topdat[i].bndtype = (int *)calloc(topdat[i].nbnd,sizeof(int));

    topdat[i].angndx1 = (int *)calloc(topdat[i].nang,sizeof(int));
    topdat[i].angndx2 = (int *)calloc(topdat[i].nang,sizeof(int));
    topdat[i].angndx3 = (int *)calloc(topdat[i].nang,sizeof(int));
    topdat[i].angtype = (int *)calloc(topdat[i].nang,sizeof(int));

    topdat[i].dihendx1 = (int *)calloc(topdat[i].ndihe,sizeof(int));
    topdat[i].dihendx2 = (int *)calloc(topdat[i].ndihe,sizeof(int));
    topdat[i].dihendx3 = (int *)calloc(topdat[i].ndihe,sizeof(int));
    topdat[i].dihendx4 = (int *)calloc(topdat[i].ndihe,sizeof(int));
    topdat[i].dihetype = (int *)calloc(topdat[i].ndihe,sizeof(int));

    topdat[i].improp_func = (int *)calloc(topdat[i].nimprop,sizeof(int));
    topdat[i].impropndx1 = (int *)calloc(topdat[i].nimprop,sizeof(int));
    topdat[i].impropndx2 = (int *)calloc(topdat[i].nimprop,sizeof(int));
    topdat[i].impropndx3 = (int *)calloc(topdat[i].nimprop,sizeof(int));
    topdat[i].impropndx4 = (int *)calloc(topdat[i].nimprop,sizeof(int));
	topdat[i].imptype = (int *)calloc(topdat[i].nimprop,sizeof(int));

    topdat[i].parm_atomtype = (int *)calloc(topdat[i].nat,sizeof(int));
    topdat[i].type    = calloc(topdat[i].nat,sizeof(*topdat[i].type));
    topdat[i].resname = calloc(topdat[i].nat,sizeof(*topdat[i].resname));
    topdat[i].segid   = calloc(topdat[i].nat,sizeof(*topdat[i].segid));

    sysdat.nats+=topdat[i].nat;
    sysdat.nbnds+=topdat[i].nbnd;
    sysdat.nangs+=topdat[i].nang;
    sysdat.ndihes+=topdat[i].ndihe;
    sysdat.total_ats+=topdat[i].nat*topdat[i].nmol;
    sysdat.total_bnds+=topdat[i].nbnd*topdat[i].nmol;
    sysdat.total_angs+=topdat[i].nang*topdat[i].nmol;
    sysdat.total_dihes+=topdat[i].ndihe*topdat[i].nmol;
    sysdat.total_improps+=topdat[i].nimprop*topdat[i].nmol;
	sysdat.nimp=sysdat.total_improps;
  }
/*     printf("BOOKKEEPING:\n"); */
/*     printf("TOTAL NONUNIQUE\n",i); */
/*     printf("FOUND: %d atoms\n",sysdat.nats); */
/*     printf("FOUND: %d bonds\n",sysdat.nbnds); */
    printf("FOUND: %d dihedrals\n",sysdat.total_dihes);
    printf("FOUND: %d impropers\n",sysdat.total_improps);

  
  rdtp=0;  
  while(rdtp < ntops){
    topdat[rdtp].nmol=atoi(argv[(2*rdtp+2)]);
    read_top(argv[(2*rdtp)+1],topdat,rdtp,&ISCHARGED);
    rdtp++;
  } 
  
  count_params(argv[argc-3],&database);
  
/*  database.fbnd=(double *)calloc(database.nbndtype,sizeof(double));*/
/*  database.bnde=(double *)calloc(database.nbndtype,sizeof(double));*/
/*  database.fang=(double *)calloc(database.nangtype,sizeof(double));*/
/*  database.ange=(double *)calloc(database.nangtype,sizeof(double));*/
/*   database.kr=(double *)calloc(database.nangtype,sizeof(double));*/
/*	 database.ro=(double *)calloc(database.nangtype,sizeof(double));*/
/*  database.eps=(double *)calloc(database.nvdwtype,sizeof(double));*/
/*  database.sig=(double *)calloc(database.nvdwtype,sizeof(double));*/
/*  database.k_dihe=(double *)calloc(database.ndihetype,sizeof(double));*/
/*  database.d_dihe=(int *)calloc(database.ndihetype,sizeof(int));*/
/*  database.n_dihe=(int *)calloc(database.ndihetype,sizeof(int));*/
/*  database.k_imp=(double *)calloc(database.nimptype,sizeof(double));*/
/*  database.d_imp=(double *)calloc(database.nimptype,sizeof(double));*/

  database.p_pair=(PARAMS *)calloc(database.nvdwtype,sizeof(PARAMS));
  database.p_bnd=(PARAMS *)calloc(database.nbndtype,sizeof(PARAMS));
  database.p_ang=(PARAMS *)calloc(database.nangtype,sizeof(PARAMS));
  database.p_dihe=(PARAMS *)calloc(database.ndihetype,sizeof(PARAMS));	
  database.p_imp=(PARAMS *)calloc(database.nimptype,sizeof(PARAMS));

  database.vdwtype1=calloc(database.nvdwtype,sizeof(*database.vdwtype1));
  database.vdwtype2=calloc(database.nvdwtype,sizeof(*database.vdwtype2));
  database.vdwstyle=calloc(database.nvdwtype,sizeof(*database.vdwstyle));
  database.bndtype1=calloc(database.nbndtype,sizeof(*database.bndtype1));
  database.bndtype2=calloc(database.nbndtype,sizeof(*database.bndtype2));
  database.angtype1=calloc(database.nangtype,sizeof(*database.angtype1));
  database.angtype2=calloc(database.nangtype,sizeof(*database.angtype2));
  database.angtype3=calloc(database.nangtype,sizeof(*database.angtype3));
  database.dihetype1=calloc(database.ndihetype,sizeof(*database.dihetype1));
  database.dihetype2=calloc(database.ndihetype,sizeof(*database.dihetype2));
  database.dihetype3=calloc(database.ndihetype,sizeof(*database.dihetype3));
  database.dihetype4=calloc(database.ndihetype,sizeof(*database.dihetype4));

  database.imptype1=calloc(database.nimptype,sizeof(*database.imptype1));
  database.imptype2=calloc(database.nimptype,sizeof(*database.imptype2));
  database.imptype3=calloc(database.nimptype,sizeof(*database.imptype3));
  database.imptype4=calloc(database.nimptype,sizeof(*database.imptype4));

  sysdat.param_bnds=(int *)calloc(database.nbndtype,sizeof(int));
  sysdat.param_angs=(int *)calloc(database.nangtype,sizeof(int));
  sysdat.param_dihes=(int *)calloc(database.ndihetype,sizeof(int));
  sysdat.param_imps=(int *)calloc(database.nimptype,sizeof(int));	
  sysdat.coordx=(double *)calloc(sysdat.total_ats,sizeof(double));
  sysdat.coordy=(double *)calloc(sysdat.total_ats,sizeof(double));
  sysdat.coordz=(double *)calloc(sysdat.total_ats,sizeof(double));
  
  //initialize default styles (CHARMM) in database
//  push_toPARAM(PARAMS* prm, char* w, int ind){
  push_toPARAM(&database.pstyle, "lj/cut/coul/long 10",0);
  set_PARAMnum(&database.pstyle, 1);	

  push_toPARAM(&database.bstyle, "harmonic",0);
  set_PARAMnum(&database.bstyle, 1);

  push_toPARAM(&database.astyle, "charmm",0);
  set_PARAMnum(&database.astyle, 1);

  push_toPARAM(&database.dstyle, "charmm",0);
  set_PARAMnum(&database.dstyle, 1);

  push_toPARAM(&database.istyle, "harmonic",0);
  set_PARAMnum(&database.istyle, 1);

  push_toPARAM(&database.kstyle, "ewald",0);
  push_toPARAM(&database.kstyle, "0.0001",1);
  set_PARAMnum(&database.kstyle, 2);
  if(sysdat.total_improps>0)
     push_toPARAM(&database.sbonds, "lj 0.0 0.0 0.0",0);
  else
     push_toPARAM(&database.sbonds, "lj 0.0 0.0 1.0",0);	
 
  set_PARAMnum(&database.sbonds, 1);	

  read_database(argv[argc-3],&database);

  printf("FOUND %d UNIQUE VDW PAIR PARAMS\n",database.nvdwtype);
  printf("FOUND %d UNIQUE BOND PARAMS\n",database.nbndtype);
  printf("FOUND %d UNIQUE ANGLE PARAMS\n",database.nangtype);
  printf("FOUND %d UNIQUE DIHEDRAL PARAMS\n",database.ndihetype);
  printf("FOUND %d UNIQUE IMPROPER PARAMS\n",database.nimptype);
 get_unique(&database,topdat,&sysdat,ISCHARGED);

  read_coords(argv[argc-1],&database,topdat,&sysdat,ISCHARGED,argv[argc-2]);
  
  return (0);
}
