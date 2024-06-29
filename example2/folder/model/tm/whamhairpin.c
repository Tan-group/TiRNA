#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define N 15
#define M 200000
#define co 1
#define bins 100
#define KB 0.00198292923700
FILE *fpE[N],*fpBp[N],*fpP,*fp2,*fp10,*fp3,*fp4,*fp5,*fp6;
/****************************************U************************************************************/
int f,k;
float UU[N][M],Ulj[N][M],Ubond[N][M],Ubp[N][M],Ust[N][M],Ujd[N][M],Ucs[N][M];
/*************BP******************************************************/
int bh,bp[N][M],BBP[N][M],bk,Bmax,kk=1,PPmax[N],/*BP[10000000]*/tt,BP[N][M];
/***********************************************************************/

double log_p0[50][M],Ra,pp0[50][M],p0[50][M];
int NNU,NB,MT,ii,NM,KU,Nn[N][50][M],n[50][M];
double TL,max;
float /*SU[1000000]*/ Emin,Emax,gap,EE[M];
double log_cc[N][M];
double dde[M],log_de[N][M];
float CV[10000],BP_tm[10000],max_bp;
double FF[50][50][1000],FFF[100];
float t[20]={25.0,31.0,37.0,43.0,49.0,55.0,63.0,71.0,78.0,86.0,94.0,102.0,110.0,120,130};
int main()
{ 
          
  	void U_N(),BP_N(),Gap(),N_State(),T_bais(),Conf_T_state(),Probability(),Cv_Tm(),Gap(),N_State(),cvv(),BP_TM();
  	float tm;
   	fpP=fopen("Probability.dat","w+");
    	fp2=fopen("thermo.dat","w+");
        //  fp10=fopen("check.dat","w+");
  	fp3=fopen("information.dat","w+");
    	fp4=fopen("cv_tm.dat","w+");
    	fp5=fopen("BP_tm.dat","w+");
   	fp6=fopen("thermal_stability.dat","w+");
  	U_N();           //Read energy and base-pairs
  	BP_N();
   	Gap();
    	N_State();
      	for(MT=1;MT<1400;MT++)
      	{
       		CV[MT]=0.0;
               	TL=273.15+MT*0.1; 
               	T_bais();                   
               	Probability();                    //Record the conformation number of evert state
               	Cv_Tm();
    	}
     	cvv();
       	BP_TM();
       	printf(">>>>>>>>>> Calculating thermal stabilites with WHAM\n");
   	fclose(fpP);
   	fclose(fp2);
    	fclose(fp3);
   	fclose(fp4);
   	fclose(fp5);
    	fclose(fp6);
      	return 0; 
}   

/********************************************************************************************************************************/
void U_N()
{
 	int i,j,nk,nc[20],min;
   	char filename[30];
   	for(i=0;i<N;i++)
    	{
      		sprintf(filename,"fragment/Energy_%d.dat",i); fpE[i]=fopen(filename,"r+");
     	}
   	for(i=0;i<N;i++)
     	{
               	j=1;
               	while(!feof(fpE[i]))
               	{
                     
                     	fscanf(fpE[i],"%d %f %f %f %f %f %f %f\n",&f,&UU[i][j],&Ulj[i][j],&Ubond[i][j],&Ubp[i][j],&Ust[i][j],&Ujd[i][j],&Ucs[i][j]);
               		j++;
               	}
               	NNU=j-1;             
               	nc[i]=NNU;
               	fclose(fpE[i]); 
               	fprintf(fp3,"replica is %d  Energy  the conformation is %d\n",i,NNU);
 	}
 	min=nc[0];
    	for(i=0;i<N;i++)
   	{
            	if(nc[i]<=min)
               	{
                     	min=nc[i];
               	}
   	}
   	NNU=min;
    	nk=1;
}

void BP_N()
{
	int i,j,NB,nk,nc[20],min;
    	char filename[30];
   	for(i=0;i<N;i++)
   	{
       		sprintf(filename,"fragment/Bp_%d.dat",i); fpBp[i]=fopen(filename,"r+");
     	}
   	for(i=0;i<N;i++)
     	{
          	j=1;
          	while(!feof(fpBp[i]))
          	{
            		fscanf(fpBp[i],"%d %d %d\n",&tt,&bp[i][j],&BP[i][j]);
                       	j++;
           	}
       		NB=j-1;
               	nc[i]=NB;
               	fprintf(fp3,"replica: %d  Bp  the conformation: %d\n",i,NB);
               	fclose(fpBp[i]);
	} 
  	min=nc[0];
   	for(i=0;i<N;i++)
   	{
        	if(nc[i]<=min)
               	{
                     	min=nc[i];
               	}
   	}
  	NB=min;
           //printf("the Energy number of min replica is %d\n",NB);
  	if(NB<NNU)
     	{
        	NNU=NB;
    	}
    	else
  	{
      		NNU=NNU;
     	}        
  	nk=1;
   	NM=nk-1;                                
} 

void Gap()
{
	int i,j,NN,N_state,ei,ej;
	float EEmin[10000];
  	Emin=UU[0][1];
   	Emax=UU[0][1];
   	Bmax=0;   
    	for(i=0;i<N;i++)
  	{
      		for(j=1;j<=NNU;j++)
        	{
            		if(UU[i][j]<=Emin)
               		{
                   		Emin=UU[i][j];
                   	}
              		if(UU[i][j]>=Emax)
                	{
                     		Emax=UU[i][j];
                          	ei=i;
                         	ej=j;
                	}
               		if(bp[i][j]>=Bmax)
                	{
                        	Bmax=bp[i][j];
                 	}
      		}
	}
  	if(Bmax<=6) 
    	{
        	for(i=0;i<N;i++)
         	{
              		for(j=1;j<=NNU;j++)
                	{
                    		BBP[i][j]=BP[i][j];
                  	}
        	}
 	}
  	else
  	{ 
        	for(i=0;i<N;i++)
         	{
                	for(j=1;j<=NNU;j++)
                	{
                    		BBP[i][j]=bp[i][j];
                   	}
          	}
   	}	
   	gap=(Emax-Emin)/bins; 
   	for(k=1;k<=bins;k++)
     	{
        	EEmin[k]=Emin+(k-1)*gap;
         	EE[k]=EEmin[k]+gap*0.5;
               //  fprintf(fp10,"%d %f\n",k,EE[k]);
    	}
   	N_state=Bmax*bins;
 	fprintf(fp3,"number of state:  %d   Emin = %f   Emax = %f  gap= %f Bmax = %d %d %d\n",N_state,Emin,Emax,gap,Bmax,ei,ej); 
}
void N_State()
{	
	int m,j,k,i,sum=0;
   	for(i=0;i<N;i++)
  	{
      		for(j=0;j<=Bmax;j++)
             	{
                 	for(k=1;k<=bins;k++)
                 	{
                         	Nn[i][j][k]=0;
                 	}
             	}
         }
  	for(j=0;j<=Bmax;j++)
 	{
         	for(k=1;k<=bins;k++)
           	{
                	n[j][k]=0;
          	}
   	}
  	for(i=0;i<N;i++)
  	{
       		for(j=0;j<=Bmax;j++)
                {
                     	for(m=1;m<=NNU;m++)
                     	{
                        	if(BBP[i][m]==j)
                          	{
                               		if(UU[i][m]>(Emin+gap))
                                	{                                       
                                  		k=((int)((UU[i][m]-Emin)/gap))+1;                                                 
                                       		Nn[i][j][k]++;
                                	}
                                 	else
                              		{
                                     		Nn[i][j][1]++;
                                	}                                              
                    		}                          
                     	}
       		}          
 	}
 	for(j=0;j<=Bmax;j++)
  	{	
  		for(k=1;k<=bins;k++)
                {
                  	for(i=0;i<N;i++)
                         {
                                n[j][k]=n[j][k]+Nn[i][j][k];
                         }
                     //    fprintf(fp10,"the state is j=%d,k=%d    nij = %d\n",j,k,n[j][k]);
              	}
 	}     
}


void T_bais()
{
    	int i,k,j;
    	float T[20];
    	for(i=0;i<N;i++)
    	{           
          	T[i]=t[i]+273.15;  
          	for(k=1;k<=bins;k++)
          	{
                  	if(TL>T[i])
                        {
                             	log_cc[i][k]=-((1.0/(KB*T[i]))-(1.0/(KB*TL)))*EE[k];
                        }
                        else
                        {
                                log_cc[i][k]=-((1.0/(KB*T[i]))-(1.0/(KB*TL)))*EE[k];
                        }
 
           	}          
     	}
}
void Probability()
{  
     	int j,i,nn=0,ii=1,k;
     	double pp,P0(int i),log_fi[20];
     	for(i=0;i<N;i++)
     	{
           	log_fi[i]=0.0;
     	}
     	for(ii=1;ii<=50;ii++)
     	{
            	for(i=0;i<N;i++)
            	{
                	for(k=1;k<=bins;k++)
                       	{
                   		log_de[i][k]=0.0;
                             	dde[k]=0.0;
                             	for(j=0;j<=Bmax;j++)
                             	{
                                    	FF[i][j][k]=0.0;
                             	}
                       	}
                    	FFF[i]=0.0;                      
             	}  
             	for(k=1;k<=bins;k++)
             	{
                     	for(i=0;i<N;i++)
                   	{ 
                      		log_de[i][k]=log(NNU)+log_cc[i][k]+log_fi[i];
                          	dde[k]=dde[k]+exp(log_de[i][k]);
                 	}  
             	} 
              	for(j=0;j<=Bmax;j++)
              	{
                       	for(k=1;k<=bins;k++)
                       	{
                            	log_p0[j][k]=log(n[j][k])-log(dde[k]);
                              	p0[j][k]=exp(log_p0[j][k]);             
                       	}
               	}  
               	for(i=0;i<N;i++)
               	{
                 	for(j=0;j<=Bmax;j++)
               		{
                     		for(k=1;k<=bins;k++)
                       		{
                           		FF[i][j][k]=log_cc[i][k]+log_p0[j][k];
                                        FFF[i]=FFF[i]+exp(FF[i][j][k]);
                       		}
                 	}
           		log_fi[i]=-1.0*log(FFF[i]);
               	}
   	}
    	fprintf(fp3,"Temperature:  %.1f    the P0 iter =  %d  times\n",TL-273.15,ii-1);
    	pp=0.0;
    	for(j=0;j<=Bmax;j++)
     	{
            	for(k=1;k<=bins;k++)
            	{
                  	pp=pp+p0[j][k];
            	}
    	}
    	Ra=1.0/pp;  
    	for(j=0;j<=Bmax;j++)
    	{
           	for(k=1;k<=bins;k++)
           	{
                	pp0[j][k]=p0[j][k]*Ra; 
            	}
    	} 
}


void Cv_Tm()
{
    	int j,k;
       	float xxEE=0.0,ee=0.0,cv=0.0,xxBP=0.0,xxbp=0.0,bp_tm=0.0;
       	for(j=0;j<=Bmax;j++)
       	{
          	for(k=1;k<=bins;k++)
          	{
               		xxEE=xxEE+EE[k]*EE[k]*pp0[j][k];
               		ee=ee+EE[k]*pp0[j][k];
               		xxbp=xxbp+j*pp0[j][k];
          	}
     	}
       	cv=(xxEE-ee*ee)/(KB*TL*TL);
       //   printf("%f %f\n",TL,cv); 
       	if(MT==1)
     	{
           	max_bp=xxbp;
       	} 
       	CV[MT]=cv;
       	BP_tm[MT]=xxbp*1.0/max_bp;
       	float fold,unfold;
       	fold=xxbp/max_bp*1.0;
       	unfold=1-fold;
       	fprintf(fp2,"%f %f %f\n",TL-273.15,cv,fold); fflush(fp2);  
       	float ti,I_state;
       	ti=TL-273.15;
       	I_state=1-unfold-fold;
       	if(ti==25.0||ti==31.0||ti==37.0||ti==43.0||ti==49.0||ti==55.0||ti==63.0||ti==71.0||ti==78.0||ti==86.0||ti==94.0||ti==102.0||ti==110.0||ti==120.0||ti==130.0)
       	{
       		fprintf(fp6,"%f %f %f %f\n",ti,unfold,I_state,fold);
       	}
       	fflush(fp6);   
}
void cvv()
{
      	int i,j,k=1;
      	float tm,max;
      	max=CV[1];
      	for(i=1;i<1199;i++)
      	{
               	if(CV[i]>=max)
                {
                      max=CV[i];
                      tm=i*0.1;
                }   
      	}
      	fprintf(fp4,"%.1f\n",tm);     
}


void BP_TM()
{
      	int i;
      	float bp_min[10000],tm,min;
      	for(i=1;i<=1399;i++)
      	{
             	bp_min[i]=fabs(BP_tm[i]-0.5);
      	}
      	min=bp_min[1];
      	for(i=1;i<=1399;i++)
      	{
           	if(bp_min[i]<=min)
            	{
                    	min=bp_min[i];
                    	tm=i*0.1;
           	}
      	}
      	fprintf(fp3,"******   the RNA--------BP:  tm  =  %.1f   **********\n",tm);
      	fprintf(fp5,"%.1f\n",tm);fflush(fp5);
}
