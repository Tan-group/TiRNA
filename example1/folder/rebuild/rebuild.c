/**********
This program is used to convert CG conf. to All-Atom structure
***********/
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <arpa/inet.h>
#include <stdlib.h>
#include <time.h>
#include <dirent.h>
#define     w   0.02         
#define  RMSD_limit 3.0      
#define cut 1.60
#define nm_frag 5
FILE *fragA[20],*fragG[20],*fragC[20],*fragU[20],*CG,*out_AA,*out_rmsd;
FILE *inpdb, *outpdb;
char atom[20][1000],typen_A[20][100][4],typen_G[20][100][4],typen_C[20][100][4],typen_U[20][100][4],type[500][4],typen_AA[20][500][4];
char type_A[20][100][5],type_G[20][100][5],type_C[20][100][5],type_U[20][100][5],type_AA[20][500][5];
float x_A[20][5][100],x_G[20][5][100],x_C[20][5][100],x_U[20][5][100],x[5][500];
int id[20],idnu[20],duo1,duo2,NA[20],NG[20],NC[20],NU[20],N,N_AA[20],nucl_AA,atom_AA[20][500],nucl_AA1,duange;
char mark1,mark2,mark3,chain[20];
float r,Q,f;
float x_CG[10][5][10],y_CG[5][10],R[5][5],R_T[5][5], xmiddle[20][5],ymiddle[5];
double eigenvalue[5], lambda, Tao[10][10], A[10][10], U[20][10][10], B[10][10], dis, dis_r;
float x_AA[20][10][100];
char atomCG[1000][5],typeACG[1000][5],chainCG[1000][5];
int idCG[1000],idnuCG[1000];
char QatomCG[1000][5],QtypeACG[1000][5],QchainCG[1000][5],Qtype[1000][5];
int QidCG[1000],QidnuCG[1000];
float Qx[5][1000];
int atom_number;
int condition;
/******************************************************************/
FILE *op;
char Oatom[10000],Otype_A[10000][5],Otypen_A[10000][5],Ochain[10000];
int Oid[10000],Oidnu[10000];
float Ox[10000],Oy[10000],Oz[10000],Omark1,Omark2;
int N_S;
int Npdb,Nii;
char infileN[100][50],temp[4],fileN[500][4];
float RMSD[20];
FILE *fragGC[10],*fragGU[10],*fragAU[10],*fragCG[10],*fragUG[10],*fragUA[10];
/*******************************************************************/

/********************************stem**************************************/
char type_GC[10][100][5],type_GU[10][100][5],type_AU[10][100][5],type_CG[10][100][5],type_UA[10][100][5],type_UG[10][100][5];
char typen_GC[10][100][5],typen_GU[10][100][5],typen_AU[10][100][5],typen_CG[10][100][5],typen_UA[10][100][5],typen_UG[10][100][5];
float xx_GC[10][5][100],xx_GU[10][5][100],xx_AU[10][5][100],xx_CG[10][5][100],xx_UG[10][5][100],xx_UA[10][5][100];
int c[500],se[500][500],NGC[10],NGU[10],NAU[10],NCG[10],NUG[10],NUA[10],A_s;
char sse[10000],seq[1000];
/********************************************************************/
FILE *fse;
int main(int argc, char *argv[])
{
  
    void Input(void),optimize(),reconstruction(),secondary();
    int i,readpdb();
    FILE *fpcs;
    Input();
    inpdb=fopen("CG.pdb","r+");   
 //   fse=fopen("ss.dat","r+");  
    outpdb=fopen("All_atom.pdb","w+");
    N=readpdb();
    fclose(inpdb);
    fpcs=fopen("cs.dat","r+");
    int j,ai,aj;
    char ax,ay;
    for(i=1;i<N;i++)
    {
	for(j=i+1;j<=N;j++)
	{
		se[i][j]=0;
		c[i]=0;c[j]=0;
	}
    }
    while(!feof(fpcs))
    {
    	fscanf(fpcs,"%d %c %d %c\n",&ai,&ax,&aj,&ay);
    	se[ai*3][aj*3]=1;
    	c[ai*3]=1;c[aj*3]=1;
    }
    fclose(fpcs);
   /* printf("Sequence:");
    scanf("%s",seq);
    printf("Secondary structure:");
    scanf("%s",sse);*/
  /*  while(!feof(fse))
    {
         fscanf(fse,"%s\n",sse);
    }
    fclose(fse);*/
    if(QtypeACG[1][0]!='P')
    {
         condition=1;
         for(i=1;i<=N;i++)
         {     
            strcpy(atomCG[i+1],QatomCG[i]);
            idCG[i+1]=QidCG[i];
            strcpy(typeACG[i+1],QtypeACG[i]);
            strcpy(type[i+1],Qtype[i]);
            strcpy(chainCG[i+1],QchainCG[i]);
            idnuCG[i+1]=QidnuCG[i];
            x[1][i+1]=Qx[1][i];
            x[2][i+1]=Qx[2][i];
            x[3][i+1]=Qx[3][i];
          }
          strcpy(atomCG[1],QatomCG[1]);
          idCG[1]=QidCG[1];
          strcpy(typeACG[1],"P");
          strcpy(type[1],Qtype[1]);
          strcpy(chainCG[1],QchainCG[1]);
          idnuCG[1]=QidnuCG[1];
          x[1][1]=Qx[1][1]+2.00;
          x[2][1]=Qx[2][1]+2.00;
          x[3][1]=Qx[3][1]+2.00;
          strcpy(atomCG[N+2],QatomCG[N]);
          idCG[N+2]=QidCG[N];
          strcpy(typeACG[N+2],"P");
          strcpy(type[N+2],Qtype[N]);
          strcpy(chainCG[N+2],QchainCG[N]);
          idnuCG[N+1]=QidnuCG[N];
          x[1][N+2]=Qx[1][N]+2.00;
          x[2][N+2]=Qx[2][N]+2.00;
          x[3][N+2]=Qx[3][N]+2.00;
          N=N+2;
      
   }
   else
   { 
       condition=0;
       for(i=1;i<=N;i++)
       {       
          strcpy(atomCG[i],QatomCG[i]);
          idCG[i]=QidCG[i];
          strcpy(typeACG[i],QtypeACG[i]);
          strcpy(type[i],Qtype[i]);
          strcpy(chainCG[i],QchainCG[i]);
          idnuCG[i]=QidnuCG[i];
          x[1][i]=Qx[1][i];
          x[2][i]=Qx[2][i];
          x[3][i]=Qx[3][i];
        }
        strcpy(atomCG[N+1],QatomCG[N]);
        idCG[N+1]=QidCG[N];
        strcpy(typeACG[N+1],"P");
        strcpy(type[N+1],Qtype[N]);
        strcpy(chainCG[N+1],QchainCG[N]);
        idnuCG[N+1]=QidnuCG[N];
        x[1][N+1]=Qx[1][N]+2.00;
        x[2][N+1]=Qx[2][N]+2.00;
        x[3][N+1]=Qx[3][N]+2.00;
        N=N+1;              
   }
   // secondary();
    reconstruction();
    optimize(); 
    fclose(outpdb);
   // printf("Rebuilding completed!\n");
    //printf("File name: All_atom.pdb\n"); 
    return 0;
}


void reconstruction()
{


    void Assemble_l(),Assemble_s();
    atom_number=1;
    A_s=0;
    Assemble_l();
    A_s=1;
    Assemble_s();
    

}
void Input(void)
{
    int i,ii;
     char filename[40];
     for(ii=1;ii<=nm_frag;ii++)
     {
                sprintf(filename,"fragment/A/fragA%d.pdb",ii); fragA[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/G/fragG%d.pdb",ii); fragG[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/C/fragC%d.pdb",ii); fragC[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/U/fragU%d.pdb",ii); fragU[ii]=fopen(filename,"r+");

        i = 0;
    while(!feof(fragA[ii]))
    {
        i++;
        fscanf(fragA[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_A[ii][i],typen_A[ii][i],&chain[ii],&idnu[ii],&x_A[ii][1][i],&x_A[ii][2][i],&x_A[ii][3][i],&mark1,&mark2,&mark3);
    }
    NA[ii]=i;
  
    i = 0;
    while(!feof(fragG[ii]))
    {
        i++;
        fscanf(fragG[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_G[ii][i],typen_G[ii][i],&chain[ii],&idnu[ii],&x_G[ii][1][i],&x_G[ii][2][i],&x_G[ii][3][i],&mark1,&mark2,&mark3); 
    }
    NG[ii]=i;
    i = 0;
    while(!feof(fragC[ii]))
    {
         i++;
         fscanf(fragC[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_C[ii][i],typen_C[ii][i],&chain[ii],&idnu[ii],&x_C[ii][1][i],&x_C[ii][2][i],&x_C[ii][3][i],&mark1,&mark2,&mark3);
    }
    NC[ii]=i;
    i = 0;
    while(!feof(fragU[ii]))
    {
         i++;
         fscanf(fragU[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_U[ii][i],typen_U[ii][i],&chain[ii],&idnu[ii],&x_U[ii][1][i],&x_U[ii][2][i],&x_U[ii][3][i],&mark1,&mark2,&mark3);
    }
    NU[ii]=i;          
   fclose(fragA[ii]);
   fclose(fragG[ii]);
   fclose(fragC[ii]);
   fclose(fragU[ii]);           
     }  

   for(ii=1;ii<=nm_frag;ii++)
     {
                sprintf(filename,"fragment/stem/GC/GC%d.pdb",ii); fragGC[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/stem/GU/GU%d.pdb",ii); fragGU[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/stem/AU/AU%d.pdb",ii); fragAU[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/stem/CG/CG%d.pdb",ii); fragCG[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/stem/UG/UG%d.pdb",ii); fragUG[ii]=fopen(filename,"r+");
                sprintf(filename,"fragment/stem/UA/UA%d.pdb",ii); fragUA[ii]=fopen(filename,"r+");
   i=0;
   while(!feof(fragGC[ii]))
   {
        i++;
        fscanf(fragGC[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_GC[ii][i],typen_GC[ii][i],&chain[ii],&idnu[ii],&xx_GC[ii][1][i],&xx_GC[ii][2][i],&xx_GC[ii][3][i],&mark1,&mark2,&mark3);
    }
    NGC[ii]=i;
   i=0;
   while(!feof(fragGU[ii]))
   {
        i++;
        fscanf(fragGU[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_GU[ii][i],typen_GU[ii][i],&chain[ii],&idnu[ii],&xx_GU[ii][1][i],&xx_GU[ii][2][i],&xx_GU[ii][3][i],&mark1,&mark2,&mark3);
    }
    NGU[ii]=i;
  i=0;
   while(!feof(fragAU[ii]))
   {
        i++;
        fscanf(fragAU[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_AU[ii][i],typen_AU[ii][i],&chain[ii],&idnu[ii],&xx_AU[ii][1][i],&xx_AU[ii][2][i],&xx_AU[ii][3][i],&mark1,&mark2,&mark3);
    }
    NAU[ii]=i;
i=0;
   while(!feof(fragCG[ii]))
   {
        i++;
        fscanf(fragCG[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_CG[ii][i],typen_CG[ii][i],&chain[ii],&idnu[ii],&xx_CG[ii][1][i],&xx_CG[ii][2][i],&xx_CG[ii][3][i],&mark1,&mark2,&mark3);
    }
    NCG[ii]=i;
   i=0;
   while(!feof(fragUG[ii]))
   {
        i++;
        fscanf(fragUG[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_UG[ii][i],typen_UG[ii][i],&chain[ii],&idnu[ii],&xx_UG[ii][1][i],&xx_UG[ii][2][i],&xx_UG[ii][3][i],&mark1,&mark2,&mark3);
    }
    NUG[ii]=i;
   i=0;
   while(!feof(fragUA[ii]))
   {
        i++;
        fscanf(fragUA[ii],"%s %d %s %s %s %d %f %f %f %s %s %s\n",&atom[ii][i],&id[ii],type_UA[ii][i],typen_UA[ii][i],&chain[ii],&idnu[ii],&xx_UA[ii][1][i],&xx_UA[ii][2][i],&xx_UA[ii][3][i],&mark1,&mark2,&mark3);
    }
    NUA[ii]=i;
  } 
     
   
}
void Assemble_l()
{
  int i,j,k,kk=1,ii;
  void Rotate();
  void Replace();
  long double min;
  for(i=1;i<=N;i++)
  {
      if(fmod(i,3)==0)
      {
//////////////G////////////
       if(c[i]==0)
       {        
        for(j=1;j<=3;j++)
         {    
                     y_CG[j][1]=x[j][i-2]; y_CG[j][2]=x[j][i-1]; y_CG[j][3]=x[j][i]; y_CG[j][4]=x[j][i+1];             
          } 
        for(ii=1;ii<=nm_frag;ii++)
        {
          if(strcmp(type[i],"G")==0) 
          {
                 N_AA[ii]=NG[ii];
                 for(j=1; j<=3; j++)
                 {
                    
                      x_CG[ii][j][1]=x_G[ii][j][1]; 
                      x_CG[ii][j][2]=x_G[ii][j][6]; 
                      x_CG[ii][j][3]=x_G[ii][j][13];  
                      x_CG[ii][j][4]=x_G[ii][j][24];                 
 
                  }
                  for(k=1;k<=N_AA[ii];k++)
                  {
                        atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_G[ii][k]); strcpy(typen_AA[ii][k],typen_G[ii][k]); nucl_AA=i/3;
                        for(j=1;j<=3;j++)  
                        {
                                x_AA[ii][j][k]=x_G[ii][j][k];
                        }
                  }
           }
//////////////A////////////
           else if(strcmp(type[i],"A")==0) 
           { 
                  N_AA[ii]=NA[ii];
                  for(j=1; j<=3; j++)
                  { 
                             x_CG[ii][j][1]=x_A[ii][j][1]; 
                             x_CG[ii][j][2]=x_A[ii][j][6]; 
                             x_CG[ii][j][3]=x_A[ii][j][13];  
                             x_CG[ii][j][4]=x_A[ii][j][23];                       
                  }
                  for(k=1;k<=N_AA[ii];k++)
                  {
                           atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_A[ii][k]); strcpy(typen_AA[ii][k],typen_A[ii][k]); nucl_AA=i/3;
                           for(j=1;j<=3;j++)  
                           {
                                  x_AA[ii][j][k]=x_A[ii][j][k];
                           }
                   }
            }
//////////////C////////////
            else if(strcmp(type[i],"C")==0) 
            {
                     N_AA[ii]=NC[ii];
                     for(j=1; j<=3; j++)
                     {

                     
                           x_CG[ii][j][1]=x_C[ii][j][1]; 
                           x_CG[ii][j][2]=x_C[ii][j][6]; 
                           x_CG[ii][j][3]=x_C[ii][j][13];  
                           x_CG[ii][j][4]=x_C[ii][j][21]; 
                            
                     }
                     for(k=1;k<=N_AA[ii];k++)
                     {
                            atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_C[ii][k]); strcpy(typen_AA[ii][k],typen_C[ii][k]); nucl_AA=i/3;
                            for(j=1;j<=3;j++)  
                            {
                                   x_AA[ii][j][k]=x_C[ii][j][k];
                             }
                      }
             }
//////////////U////////////
             else if(strcmp(type[i],"U")==0) 
             {
                      N_AA[ii]=NU[ii];
                      for(j=1; j<=3; j++)
                      {
        
                              x_CG[ii][j][1]=x_U[ii][j][1]; 
                              x_CG[ii][j][2]=x_U[ii][j][6]; 
                              x_CG[ii][j][3]=x_U[ii][j][13]; 
                              x_CG[ii][j][4]=x_U[ii][j][21]; 
                            
                      }
                      for(k=1;k<=N_AA[ii];k++)
                      {
                             atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_U[ii][k]); strcpy(typen_AA[ii][k],typen_U[ii][k]); nucl_AA=i/3;
                             for(j=1;j<=3;j++)  
                             {
                                         x_AA[ii][j][k]=x_U[ii][j][k];
                             }
                       }
               }
               else 
               {
                       // printf("Error: There is a usual nucleotide\n"); break;
               }
                Rotate(x_CG,y_CG,ii,4);    
    }
       min=RMSD[1];
       for(ii=1;ii<=nm_frag;ii++)
       {
             if(RMSD[ii]<=min) { min=RMSD[ii]; kk=ii;}
       }
       Replace(kk,3);
      }   
    }
  } 
}


void Assemble_s()
{
  int i,sk,k,j,ii=1,kk=1;
  void Rotate();
  void Replace1(); 
  float min;
  for(i=1;i<=N;i++)
  {
    if(fmod(i,3)==0)
    {
      for(sk=i+3;sk<=N;sk++)
      {
        if(fmod(sk,3)==0)
        {   
          if(se[i][sk]==1)
          {         
            for(j=1;j<=3;j++)
            {   
               y_CG[j][1]=x[j][i-2];  y_CG[j][2]=x[j][i-1];  y_CG[j][3]=x[j][i]; y_CG[j][4]=x[j][i+1];        
               y_CG[j][5]=x[j][sk-2]; y_CG[j][6]=x[j][sk-1]; y_CG[j][7]=x[j][sk]; y_CG[j][8]=x[j][sk+1];
             }                             
              for(ii=1;ii<=nm_frag;ii++)
              {
                if(strcmp(type[i],"G")==0&&strcmp(type[sk],"C")==0) 
                {                      
                  N_AA[ii]=NGC[ii];
                  for(j=1; j<=3; j++)
                  {                     
                            x_CG[ii][j][1]=xx_GC[ii][j][1]; 
                            x_CG[ii][j][2]=xx_GC[ii][j][6]; 
                            x_CG[ii][j][3]=xx_GC[ii][j][13]; 
                            x_CG[ii][j][4]=xx_GC[ii][j][24]; 
                            x_CG[ii][j][5]=xx_GC[ii][j][27];
                            x_CG[ii][j][6]=xx_GC[ii][j][32];
                            x_CG[ii][j][7]=xx_GC[ii][j][39];
                            x_CG[ii][j][8]=xx_GC[ii][j][47];                               
                   }
                   for(k=1;k<=N_AA[ii];k++)
                   {
                     atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_GC[ii][k]); strcpy(typen_AA[ii][k],typen_GC[ii][k]); nucl_AA=i/3; 
                     duange=24;nucl_AA1=sk/3;
                      for(j=1;j<=3;j++)  { x_AA[ii][j][k]=xx_GC[ii][j][k];}
                    }
                 }
                 if(strcmp(type[i],"G")==0&&strcmp(type[sk],"U")==0) 
                 {
                   N_AA[ii]=NGU[ii];
                   for(j=1; j<=3; j++)
                   {                     
                           x_CG[ii][j][1]=xx_GU[ii][j][1]; 
                           x_CG[ii][j][2]=xx_GU[ii][j][6]; 
                           x_CG[ii][j][3]=xx_GU[ii][j][13];  
                           x_CG[ii][j][4]=xx_GU[ii][j][24];
                           x_CG[ii][j][5]=xx_GU[ii][j][27];
                           x_CG[ii][j][6]=xx_GU[ii][j][32]; 
                           x_CG[ii][j][7]=xx_GU[ii][j][39]; 
                           x_CG[ii][j][8]=xx_GU[ii][j][47];                        
                   }
                   for(k=1;k<=N_AA[ii];k++)
                   {
                      atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_GU[ii][k]); strcpy(typen_AA[ii][k],typen_GU[ii][k]); nucl_AA=i/3;
                      duange=24;nucl_AA1=sk/3;
                      for(j=1;j<=3;j++)  { x_AA[ii][j][k]=xx_GU[ii][j][k]; }
                   }
                  }
                  if(strcmp(type[i],"A")==0&&strcmp(type[sk],"U")==0) 
                  {
                     N_AA[ii]=NAU[ii];
                     for(j=1; j<=3; j++)
                     {                     
                           x_CG[ii][j][1]=xx_AU[ii][j][1]; 
                           x_CG[ii][j][2]=xx_AU[ii][j][6]; 
                           x_CG[ii][j][3]=xx_AU[ii][j][13];  
                           x_CG[ii][j][4]=xx_AU[ii][j][23];
                           x_CG[ii][j][5]=xx_AU[ii][j][26];
                           x_CG[ii][j][6]=xx_AU[ii][j][31]; 
                           x_CG[ii][j][7]=xx_AU[ii][j][38];
                           x_CG[ii][j][8]=xx_AU[ii][j][46];                        
                     }
                     for(k=1;k<=N_AA[ii];k++)
                     {
                        atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_AU[ii][k]); strcpy(typen_AA[ii][k],typen_AU[ii][k]); nucl_AA=i/3;
                        duange=23;nucl_AA1=sk/3;
                        for(j=1;j<=3;j++) { x_AA[ii][j][k]=xx_AU[ii][j][k]; }
                     }
                  }
                  if(strcmp(type[i],"C")==0&&strcmp(type[sk],"G")==0) 
                  {
                     N_AA[ii]=NCG[ii];
                     for(j=1; j<=3; j++)
                     {                     
                           x_CG[ii][j][1]=xx_CG[ii][j][1]; 
                           x_CG[ii][j][2]=xx_CG[ii][j][6]; 
                           x_CG[ii][j][3]=xx_CG[ii][j][13];  
                           x_CG[ii][j][4]=xx_CG[ii][j][21];
                           x_CG[ii][j][5]=xx_CG[ii][j][24];
                           x_CG[ii][j][6]=xx_CG[ii][j][29]; 
                           x_CG[ii][j][7]=xx_CG[ii][j][36];
                           x_CG[ii][j][8]=xx_CG[ii][j][47];                    
                     }
                     for(k=1;k<=N_AA[ii];k++)
                     {
                         atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_CG[ii][k]); strcpy(typen_AA[ii][k],typen_CG[ii][k]); nucl_AA=i/3;
                         duange=21;nucl_AA1=sk/3;
                         for(j=1;j<=3;j++) { x_AA[ii][j][k]=xx_CG[ii][j][k]; }
                     }
                  }
                  if(strcmp(type[i],"U")==0&&strcmp(type[sk],"G")==0) 
                  {
                     N_AA[ii]=NUG[ii];
                     for(j=1; j<=3; j++)
                     {                     
                           x_CG[ii][j][1]=xx_UG[ii][j][1]; 
                           x_CG[ii][j][2]=xx_UG[ii][j][6]; 
                           x_CG[ii][j][3]=xx_UG[ii][j][13];  
                           x_CG[ii][j][4]=xx_UG[ii][j][21];
                           x_CG[ii][j][5]=xx_UG[ii][j][24];
                           x_CG[ii][j][6]=xx_UG[ii][j][29];
                           x_CG[ii][j][7]=xx_UG[ii][j][36];
                           x_CG[ii][j][8]=xx_UG[ii][j][47];                     
                     }
                     for(k=1;k<=N_AA[ii];k++)
                     {
                          atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_UG[ii][k]); strcpy(typen_AA[ii][k],typen_UG[ii][k]); nucl_AA=i/3;
                          duange=21;nucl_AA1=sk/3;
                          for(j=1;j<=3;j++) { x_AA[ii][j][k]=xx_UG[ii][j][k];}
                     }
                  }
                  if(strcmp(type[i],"U")==0&&strcmp(type[sk],"A")==0) 
                  {
                     N_AA[ii]=NUA[ii];
                     for(j=1; j<=3; j++)
                     {                     
                           x_CG[ii][j][1]=xx_UA[ii][j][1]; 
                           x_CG[ii][j][2]=xx_UA[ii][j][6]; 
                           x_CG[ii][j][3]=xx_UA[ii][j][13];  
                           x_CG[ii][j][4]=xx_UA[ii][j][21];
                           x_CG[ii][j][5]=xx_UA[ii][j][24];
                           x_CG[ii][j][6]=xx_UA[ii][j][29]; 
                           x_CG[ii][j][7]=xx_UA[ii][j][36];
                           x_CG[ii][j][8]=xx_UA[ii][j][47];                
                     }
                     for(k=1;k<=N_AA[ii];k++)
                     {
                           atom_AA[ii][k]=i/3+k-1; strcpy(type_AA[ii][k],type_UA[ii][k]); strcpy(typen_AA[ii][k],typen_UA[ii][k]); nucl_AA=i/3;
                           duange=21;nucl_AA1=sk/3;
                           for(j=1;j<=3;j++) { x_AA[ii][j][k]=xx_UA[ii][j][k]; }
                     }
                  }
                  Rotate(x_CG,y_CG,ii,8);
              }
              min=RMSD[1];
              for(ii=1;ii<=nm_frag;ii++)
              {
                  if(RMSD[ii]<=min) { min=RMSD[ii]; kk=ii;}
              }
              Replace1(kk,3);
            }
           }
         }
      } 
    }
}




void Rotate(float xxx_CG[10][5][10],float yyy_CG[5][10],int ii,int nn)
{
  int         i,j,k,n;
  float       max, min, unit, test, xnew_CG[10][10];
  double      sigma, sigma0,xnew_CG_i;
   
   for(i=0;i<=4;i++)       //ininalization of matrix;
   {
          xmiddle[ii][i]=0.0; ymiddle[i]=0.0;
          for(j=0;j<=4;j++)
          {
                  R[i][j]=0.0; R_T[i][j]=0.0; Tao[i][j]=0.0; //U[i][j]=0.0; B[i][j]=0.0;
          }
    }

         RMSD[ii]=0.0;
    
   for(i=1;i<=nn;i++)
   {
         xmiddle[ii][1] += xxx_CG[ii][1][i]; xmiddle[ii][2] += xxx_CG[ii][2][i]; xmiddle[ii][3] += xxx_CG[ii][3][i];    
         ymiddle[1] += yyy_CG[1][i]; ymiddle[2] += yyy_CG[2][i]; ymiddle[3] += yyy_CG[3][i];    
   } 
   xmiddle[ii][1] /= nn; xmiddle[ii][2] /= nn; xmiddle[ii][3] /= nn; 
   ymiddle[1] /= nn; ymiddle[2] /= nn; ymiddle[3] /= nn; 
   for(i=1; i<=nn; i++)
   {   
        for(j=1; j<=3; j++)
        {        
            xxx_CG[ii][j][i] -= xmiddle[ii][j];   
            yyy_CG[j][i] -= ymiddle[j];    
         } 

    }  
    for(i=1; i<=3; i++)
    {
        for(j=1; j<=3; j++)
        {            
            for(k=1; k<=nn; k++)
            {
                R[i][j] += yyy_CG[i][k]*xxx_CG[ii][j][k]*w; 
            }       
         }   
     }
     for(i=1; i<=3; i++)
     {
        for(j=1; j<=3; j++)
        {
            
            R_T[i][j] = R[j][i];
         }
     }    
     for(i=1; i<=3; i++)
     {
        for(j=1; j<=3; j++)
        {
            
            for(k=1; k<=3; k++)
            {
                Tao[i][j] += R_T[i][k]*R[k][j];
            }            
         }
      }
      min = -100.; max = 10000.; unit = 0.0001; sigma0 = 0; n = 0;   dis=0.0; 
      lambda = min; 
      do
      {
        sigma = (Tao[1][1]-lambda)*(Tao[2][2]-lambda)*(Tao[3][3]-lambda) + Tao[1][2]*Tao[2][3]*Tao[3][1] + Tao[1][3]*Tao[2][1]*Tao[3][2] - Tao[1][3]*Tao[3][1]*(Tao[2][2]-lambda) - Tao[1][2]*Tao[2][1]*(Tao[3][3]-lambda) - Tao[2][3]*Tao[3][2]*(Tao[1][1]-lambda);              
        
        if(sigma * sigma0 < 0)      
        {
            n++;
            eigenvalue[n] = lambda - unit/2.;        
        }
        sigma0 = sigma;
        lambda += unit;
        
       }
       while( !(n==3 || lambda > max) );      
      for(i=0; i<=9; i++)
      {
           for(j=0; j<=9; j++)
           {
                A[i][j] = 0.;  B[i][j]=0.;U[ii][i][j]=0.; xnew_CG[i][j]=0.; 
           }
       }
       for(i=1; i<=3; i++)	 
       {
          A[1][i] = 1.;     
          A[3][i] = ( Tao[2][1] + (Tao[2][2] - eigenvalue[i])*(eigenvalue[i]-Tao[1][1])/Tao[1][2] ) / ( Tao[1][3]*(Tao[2][2]-eigenvalue[i])/Tao[1][2] - Tao[2][3]  );        
          A[2][i] = ( eigenvalue[i]-Tao[1][1]-Tao[1][3] * A[3][i] ) / Tao[1][2];            
       }
       for(i=1; i<=3; i++)
       {
        
           dis = sqrt( A[1][i]*A[1][i] + A[2][i]*A[2][i] + A[3][i]*A[3][i]  ); 
           A[1][i] /= dis;      
           A[2][i] /= dis;
           A[3][i] /= dis;       
           test = Tao[3][1] * A[1][i] + Tao[3][2] * A[2][i] + (Tao[3][3] - eigenvalue[i]) * A[3][i];    
           if(test>1.) 
           {
              // printf("Error: the fragment would be inaccurate\n");
           }
        }
    
        for(i=1; i<=3; i++)
        {
             for(j=1; j<=3; j++)
             {
                 for(k=1; k<=3; k++)
                 {
                        B[j][i] += R[j][k]*A[k][i]/ sqrt(eigenvalue[i]);     
                 }               
             }
         }
    
//    printf("\nthe U matrix\n");
         for(i=1; i<=3; i++)
         {
               for(j=1; j<=3; j++)
               {
                     for(k=1; k<=3; k++)
                     {          
                         U[ii][i][j] += B[i][k]*A[j][k];         
                     }
               }
          }
          for(k=1; k<=nn; k++)
          {
                dis_r = 0.;
                for(i=1; i<=3; i++)
                {
                      xnew_CG_i=0.; 
                      for(j=1; j<=3; j++)
                      {
                             xnew_CG_i += U[ii][i][j]*xxx_CG[ii][j][k];  
                      }  
                      xnew_CG[i][k]=xnew_CG_i;
                      dis_r += (xnew_CG[i][k]-yyy_CG[i][k])*(xnew_CG[i][k]-yyy_CG[i][k]);      
                }
                RMSD[ii] += dis_r;    
           }    
           RMSD[ii] = sqrt(RMSD[ii]/nn);    
  
           for(k=1;k<=nn; k++)
           {
               for(j=1; j<=3; j++)
               {
                    xnew_CG[j][k] += ymiddle[j];    
                    yyy_CG[j][k] += ymiddle[j];        
                }   
            }
}

void Replace(int ii,int last)
{
  int     i,j,k; 
  float   xnew[5][100];
      for(k=1;k<=N_AA[ii];k++)
      {
           for(i=1;i<=3;i++)
           {
                x_AA[ii][i][k] -= xmiddle[ii][i];   xnew[i][k]=0.;
           }  
      }

     for(k=1;k<=N_AA[ii];k++)
     {     
           for(i=1;i<=3;i++)
           {
    
                for(j=1;j<=3;j++)
                {
                     xnew[i][k] += U[ii][i][j]*x_AA[ii][j][k]; 
                }   
                xnew[i][k]=xnew[i][k]+ymiddle[i];     
            }   
     } 
     for(k=1;k<=N_AA[ii]-last;k++) 
     {
               
              
               Oid[atom_number]=atom_number;
               strcpy(Otype_A[atom_number],type_AA[ii][k]);
               strcpy(Otypen_A[atom_number],typen_AA[ii][k]);          
               Ox[atom_number]=xnew[1][k];
               Oy[atom_number]=xnew[2][k];
               Oz[atom_number]=xnew[3][k];
               if(A_s==0) {Oidnu[atom_number]=nucl_AA;}
               else
               {
                 if(k<duange) { Oidnu[atom_number]=nucl_AA; }
                 else         { Oidnu[atom_number]=nucl_AA1; }
               }
              
               atom_number++;
         
               
     }
     
}

void Replace1(int ii,int last)
{
  int     i,j,k; 
  float   xnew[5][100];
      for(k=1;k<=N_AA[ii];k++)
      {
           for(i=1;i<=3;i++)
           {
                x_AA[ii][i][k] -= xmiddle[ii][i];   xnew[i][k]=0.;
           }  
      }
     for(k=1;k<=N_AA[ii];k++)
     {     
           for(i=1;i<=3;i++)
           {
    
                for(j=1;j<=3;j++)
                {
                     xnew[i][k] += U[ii][i][j]*x_AA[ii][j][k]; 
                }   
                xnew[i][k]=xnew[i][k]+ymiddle[i];     
            }   
     } 
     for(k=1;k<duange;k++) 
     {
               
              
               Oid[atom_number]=atom_number;
               strcpy(Otype_A[atom_number],type_AA[ii][k]);
               strcpy(Otypen_A[atom_number],typen_AA[ii][k]);          
               Ox[atom_number]=xnew[1][k];
               Oy[atom_number]=xnew[2][k];
               Oz[atom_number]=xnew[3][k];
               Oidnu[atom_number]=nucl_AA;            
               atom_number++;
     }
     for(k=duange+3;k<=N_AA[ii]-last;k++) 
     {
               
              
               Oid[atom_number]=atom_number;
               strcpy(Otype_A[atom_number],type_AA[ii][k]);
               strcpy(Otypen_A[atom_number],typen_AA[ii][k]);          
               Ox[atom_number]=xnew[1][k];
               Oy[atom_number]=xnew[2][k];
               Oz[atom_number]=xnew[3][k];
               Oidnu[atom_number]=nucl_AA1;
               atom_number++;
     }
     
}


void optimize()
{
     int i=1,number,kk,lidnu[10000];
     float centerx,centery,centerz;
     float allx=0.0,ally=0.0,allz=0.0;
      //xx,yy,zz;
     char ltype_A[10000][5],ltypen_A[10000][5];
     float lx[10000],ly[10000],lz[10000];

     number=atom_number-1;
  //   printf("%d\n",number);
     for(i=1;i<=number;i++)
     {
            allx=allx+Ox[i];
            ally=ally+Oy[i];
            allz=allz+Oz[i];
    }
    centerx=allx/number;
    centery=ally/number;
    centerz=allz/number;
    for(i=1;i<=number;i++)
    {
      Ox[i]=Ox[i]-centerx;
      Oy[i]=Oy[i]-centery;
      Oz[i]=Oz[i]-centerz;
     // fprintf(fp1,"%-6s%5d %-4s %2s %s%4d %11.3f%8.3f%8.3f  %.2f  %.2f\n","ATOM",id[i],type_A[i],typen_A[i],"A",idnu[i],x[i],y[i],z[i],1.00,0.00); 
   }   
   int ri=1;       
   for(kk=1;kk<=(N-1)/3;kk++)
   {  
           for(i=1;i<=number;i++)
           {
               if(Oidnu[i]==kk)
               {  
                  strcpy(ltype_A[ri],Otype_A[i]); strcpy(ltypen_A[ri],Otypen_A[i]);lidnu[ri]=Oidnu[i];lx[ri]=Ox[i];ly[ri]=Oy[i];lz[ri]=Oz[i];
                  ri++;
                //  printf("ri   %d\n",ri);
               }
            }
    }
    if(condition==1)
   {       
     for(i=4;i<=number;i++)
     {
         
         fprintf(outpdb,"%-6s%5d  %-4s %2s %s%4d %11.3f%8.3f%8.3f\n","ATOM",Oid[i],ltype_A[i],ltypen_A[i],"A",lidnu[i],lx[i],ly[i],lz[i]); 
     }
   }
   else
   {       
     for(i=1;i<=number;i++)
     {
         
        fprintf(outpdb,"%-6s%5d  %-4s %2s %s%4d %11.3f%8.3f%8.3f\n","ATOM",Oid[i],ltype_A[i],ltypen_A[i],"A",lidnu[i],lx[i],ly[i],lz[i]); 
     }
   }   

}


int readpdb()
{

    char  numl1[10000][10],numl2[10000][10];
    char x1[10000][10], y1[10000][10], z1[10000][10];
    char a[500];
    int i,j;
    i=1;
    memset(x1, 0, sizeof(x1));
    memset(y1, 0, sizeof(y1));
    memset(z1, 0, sizeof(z1));
    while(fgets(a,500,inpdb)!=NULL)
    {
        sprintf(QatomCG[i],"%c%c%c%c",a[0],a[1],a[2],a[3]);
        sprintf(numl1[i],"%c%c%c%c",a[8],a[9],a[10],a[11]);
        QidCG[i]=atof(numl1[i]);
        sprintf(QtypeACG[i],"%c%c%c",a[13],a[14],a[15]);
        sprintf(Qtype[i],"%c",a[19]);//residue_type
        sprintf(QchainCG[i],"%c",a[21]);//chain_type
        sprintf(numl2[i],"%c%c%c%c",a[22],a[23],a[24],a[25]);//residue_number
        QidnuCG[i]=atof(numl2[i]);
        sprintf(x1[i],"%c%c%c%c%c%c%c%c",a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);//x_coordinate
        Qx[1][i]=atof(x1[i]);
        sprintf(y1[i],"%c%c%c%c%c%c%c%c",a[38],a[39],a[40],a[41],a[42],a[43],a[44],a[45]);//y_coordinate
        Qx[2][i]=atof(y1[i]);
        sprintf(z1[i],"%c%c%c%c%c%c%c%c",a[46],a[47],a[48],a[49],a[50],a[51],a[52],a[53]);//z_coordinate
        Qx[3][i]=atof(z1[i]);
        i++;
     }
     memset(a,0,sizeof(a));
     j=i-1;
     return j;
}




void secondary()
{
   char stack[1000],stack1[1000],stack2[1000];
   int xx[1000],kk,kk1,kk2,xx1[1000],xx2[1000];
   int top=-1,top1=-1,top2=-1;
   int index=0;
   int ai,aj,i,sk;
   for(i=0;i<=N;i=i+3) 
   { 
      c[i]=0;
      for(sk=i+9;sk<=N;sk=sk+3)
      {
          se[i][sk]=0;
      }
   }
   while(sse[index]!=0)
   {
          if(sse[index]=='(')
          {
               kk=top;
               stack[++top]=sse[index];
               xx[++kk]=index;
               
          }
          if(sse[index]==')')
          {
                 if(sse[index]==')'&&stack[top]=='(')
                 {
                    //fprintf(fp,"%d %d\n",xx[top]+1,index+1);
                     ai=xx[top]+1;aj=index+1;
                     se[ai*3][aj*3]=1;c[ai*3]=1;c[aj*3]=1;
                   //  printf("%d %d\n",xx[top]+1,index+1);
              
                  }
                  top--;
          }
         if(sse[index]=='[')
          {
               kk1=top1;
               stack1[++top1]=sse[index];
               xx1[++kk1]=index;               
          }
          if(sse[index]==']')
          {
                 if(sse[index]==']'&&stack1[top1]=='[')
                 {
                    //fprintf(fp,"%d %d\n",xx1[top1]+1,index+1);
                    ai=xx1[top1]+1;aj=index+1;
                    se[ai*3][aj*3]=1;c[ai*3]=1;c[aj*3]=1;
                   // printf("%d %d\n",xx1[top1]+1,index+1);
                  }
                  top1--;
           }
           if(sse[index]=='{')
          {
               kk2=top2;
               stack2[++top2]=sse[index];
               xx2[++kk2]=index;
               
          }
          if(sse[index]=='}')
          {
                 if(sse[index]=='}'&&stack2[top2]=='{')
                 {
                    //fprintf(fp,"%d %d\n",xx[top]+1,index+1);
                     ai=xx2[top2]+1;aj=index+1;
                     se[ai*3][aj*3]=1;c[ai*3]=1;c[aj*3]=1;
                   //  printf("%d %d\n",xx[top]+1,index+1);
              
                  }
                  top2--;
          } 
          index++;
       }
}

