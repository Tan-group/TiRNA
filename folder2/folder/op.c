/**************** 该程序是用来基于二级结构预测RNA三维结构 ************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#define pi 3.1415
#define total1  500000        //Max steps in Folding annealing (Structure predicton)        
//#define ht 35
#define tran 1            //rotational frequency
#define tconf1 5000       //The frequency of conformational output
#define tenergy 50000     //The frequency of energy calculation
#define tbp 1000          //The frequency of output information for base-pairing
#define fp_bp 100         //The frequency of output the value of BP & PBP & G
#define t_si 5000          //Output frequency of structure information
#define tMove 1           //The frequency of 3-nt fragment moves；
#define bl 5.45
#define IMg 3.0            //2:2 IMg=4.0; 2:1 IMg=3.0;
#define ww 0.5
#define Nthd 16
#define ratio 1.0
#define e 0.26
#define B0 -11.5
float T,Kd,I,Ek,q4,CNa,CMg,fNa,D;//9
int N,N0,N1,total,t,l,salt,Energy,Bulge,Pseudoknot,Move;//14
int s[500][500],ss[500][500],a[10000],c[10000];//6
float rand01,phi,theta,rm;//7
float t0;//24
float x[10000],y[10000],z[10000],xx[10000],yy[10000],zz[10000],xxx[10000],yyy[10000],zzz[10000],Q[1000],R[1000],f[1000];//12
float U,du,uu,u,ulj,uulj,uc,uuc,ub,uub,ue,uue,ud,uud,uN,uuN,us,uus,uco,uuco;//23
char  type[1000];//1
/***************new additional movement****************/
float rx0[1000],ry0[1000],rz0[1000],rR0[1000],rQ0[1000],rf0[1000];
char  rtype0[1000];
int   tdn;
/********The parameters (including bonded & nonbonded) of the CG Model************/
float qq4,frac[1000],fi[1000],frac0[1000],frac1[1000],frac2[1000];
/***************************************parallel********************************************/
int    rt,rrt[Nthd],rNCG[Nthd];
float rx[Nthd][1000],ry[Nthd][1000],rz[Nthd][1000],rR[Nthd][1000],rQ[Nthd][1000],rf[Nthd][1000],RU[Nthd],cca[Nthd];
float rt0[Nthd+20]={25.0,31.0,37.0,45.0,54.0,64.0,74.0,86.0,98.0,110.0};
char  rtype[Nthd][1000];
int rp,cs[Nthd][500][500],as[Nthd][2000],bp_all;
float ucst,uucst;
int xxi1,xxj1,xxi2,xxj2,xxi3,xxj3,xxi4,xxj4,xxi5,xxj5;
int ssi1,ssj1,ssi2,ssj2;
int yi1[Nthd],yj1[Nthd],yi2[Nthd],yj2[Nthd],yi3[Nthd],yj3[Nthd],yi4[Nthd],yj4[Nthd],yi5[Nthd],yj5[Nthd];
int li1[Nthd],lj1[Nthd],li2[Nthd],lj2[Nthd];
int junction[Nthd],condition[Nthd];
float Apcp[Nthd][1000],Acpc[Nthd][1000],ApcN[Nthd][1000],ANcp[Nthd][1000];
float dpcpc[Nthd][1000],dcpcp[Nthd][1000],dcpcN[Nthd][1000],dNcpc[Nthd][1000];
int cycle[Nthd];
float NS1,NS2;
int ccs[1000][1000],aas[1000],ter[Nthd][1000],tertiary[1000];
FILE *fp_conf[Nthd],*fp_sec_stru[Nthd],*fp_energy[Nthd],*Infor,*fps,*fpb,*fpb1,*fpb2;
#pragma omp threadprivate(l,total,t,salt,N0,N1,t0)
#pragma omp threadprivate(T,Kd,I,Ek,q4,D,fNa,CNa,CMg)
#pragma omp threadprivate(phi,theta,rm,rand01)
#pragma omp threadprivate(x,y,z,xx,yy,zz,xxx,yyy,zzz,Q,R,f,type)
#pragma omp threadprivate(a,c,s,ss)
#pragma omp threadprivate(u,uu,ub,ue,ud,uN,us,uc,uub,uue,uud,uuN,uus,uuc,uco,uuco,du,U,ulj,uulj)
#pragma omp threadprivate(Energy,Move)
#pragma omp threadprivate(Bulge,Pseudoknot,tdn)
#pragma omp threadprivate(qq4,frac,fi,frac0,frac1,frac2,ucst,uucst)
#pragma omp threadprivate(xxi1,xxj1,xxi2,xxj2,xxi3,xxj3,xxi4,xxj4,xxi5,xxj5)
#pragma omp threadprivate(ssi1,ssj1,ssi2,ssj2,ccs,aas,tertiary,NS1,NS2)
//FILE *fp4;
int main()
{
	void  openfile(),closefile(),remc(),detailed();
    	openfile();
    	for(rt=1;rt<=1;rt++)
  	{
           	remc();
           	detailed();
    	}
    	closefile();
    	return 0;
}

void remc()
{
    	int i,j;
    	void replica();   //input of the initial conformation        //input of the initial conformation
    	if(rt==1)
    	{            
         	for(i=1;i<=rp;i++)
              	{
                  	for(j=0;j<Nthd;j++)
                  	{
                    		rx[j][i]=rx0[i];ry[j][i]=ry0[i];rz[j][i]=rz0[i];
                    		rR[j][i]=rR0[i];rQ[j][i]=rQ0[i];rf[j][i]=rf0[i];
                    		rtype[j][i]=rtype0[i];rrt[j]=rt;rNCG[j]=rp;
                  	}     
              	}
    	}
    	else 
    	{
             	for(j=0;j<Nthd;j++)
             	{
                      	rrt[j]=rt;
             	}
    	}
    	#pragma omp parallel num_threads(Nthd) 
    	{
          tdn=omp_get_thread_num(); 
          replica();
    	}
}
/*********************************************************************************************************/
void detailed()
{      	int k;
        float rand02;
        void exchange1(),exchange2();
        srand((unsigned)time(NULL)); 
        rand02=rand()/(RAND_MAX+1.);
        k=(int)(rand02*10);
        if(fmod(k,2)==0)
        {
               exchange2();
        }
        else 
        {
               exchange1();
        }
}

void exchange1()
{
      	float rrx[1000],rry[1000],rrz[1000],rand02;
      	int i,ri,rj;
      	float rk,r1t0,rderta;
      	srand((unsigned)time(NULL)); 
      	rand02=rand()/(RAND_MAX+1.);
      	ri=0;
      	while(ri<Nthd-1)
      	{
            	rj=ri+1;
            	rk=RU[rj]-RU[ri];r1t0=1/((rt0[rj]+273.15)*2.0*pow(10,-3))-1/((rt0[ri]+273.15)*2.0*pow(10,-3));rderta=rk*r1t0;
            	rand02=rand()/(RAND_MAX+1.);
            	if(exp(rderta)>=rand02)
            	{
                    	for (i=1;i<=rp;i++)
                    	{  
                                 rrx[i]=rx[ri][i];rx[ri][i]=rx[rj][i];rx[rj][i]=rrx[i];
                                 rry[i]=ry[ri][i];ry[ri][i]=ry[rj][i];ry[rj][i]=rry[i];
                                 rrz[i]=rz[ri][i];rz[ri][i]=rz[rj][i];rz[rj][i]=rrz[i];
                    	}  
                    	//printf("replica:  %d  <----------->  %d,  exchange ratio: %f\n",ri,rj,exp(rderta));                            
             	}         
             	else
             	{
                    	for(i=1;i<=rp;i++)
                    	{
                             	rx[ri][i]=rx[ri][i];rx[rj][i]=rx[rj][i];
                              	ry[ri][i]=ry[ri][i];ry[rj][i]=ry[rj][i];
                            	rz[ri][i]=rz[ri][i];rz[rj][i]=rz[rj][i];
                    	}
                    	//printf("replica:  %d  can't exchange %d\n",ri,rj); 
             	}
             	ri=ri+2; 
 	}
}

void exchange2()
{
      	float rrx[1000],rry[1000],rrz[1000],rand02;
      	int i,ri,rj;
      	float rk,r1t0,rderta;
      	srand((unsigned)time(NULL)); 
      	rand02=rand()/(RAND_MAX+1.);
      	ri=1;
      	while(ri<Nthd-1)
      	{
            	rj=ri+1;
            	rk=RU[rj]-RU[ri];r1t0=1/((rt0[rj]+273.15)*2.0*pow(10,-3))-1/((rt0[ri]+273.15)*2.0*pow(10,-3));rderta=rk*r1t0;
            	rand02=rand()/(RAND_MAX+1.);
            	if(exp(rderta)>=rand02)
            	{
                    	for (i=1;i<=rp;i++)
                    	{  
                          	rrx[i]=rx[ri][i];rx[ri][i]=rx[rj][i];rx[rj][i]=rrx[i];
                             	rry[i]=ry[ri][i];ry[ri][i]=ry[rj][i];ry[rj][i]=rry[i];
                            	rrz[i]=rz[ri][i];rz[ri][i]=rz[rj][i];rz[rj][i]=rrz[i];
                    	}
                    	//printf("replica:  %d  <----------->  %d,  exchange ratio: %f\n",ri,rj,exp(rderta));  
             	}         
             	else
             	{
                    	for(i=1;i<=rp;i++)
                    	{
                          	rx[ri][i]=rx[ri][i];rx[rj][i]=rx[rj][i];
                          	ry[ri][i]=ry[ri][i];ry[rj][i]=ry[rj][i];
                           	rz[ri][i]=rz[ri][i];rz[rj][i]=rz[rj][i];
                    	}
                    	//printf("replica:  %d  can't exchange %d\n",ri,rj); 
             	}
             	ri=ri+2;
        }
}
void replica()
{
 	int i,j;
	void Put_File(),Fixed_Atom(),MC_Annealing();
 	Put_File();         //defined the input parameters and out file names;
 	N0=rNCG[tdn];
 	N=(N0-1)/3;          //N0:Total number of CG beads; N: Num. of nt
 	Fixed_Atom();      //The Centre atom N1 will be fixed;
 //printf("exchange time is %d, the thread is %d,  the RNA length =  %d, CG beads is %d \n",rrt[tdn],tdn,N,N0);
 	xxi1=yi1[tdn];xxj1=yj1[tdn];
 	xxi2=yi2[tdn];xxj2=yj2[tdn];
 	xxi3=yi3[tdn];xxj3=yj3[tdn];
 	xxi4=yi4[tdn];xxj4=yj4[tdn];
 	xxi5=yi5[tdn];xxj5=yj5[tdn];
 	ssi1=li1[tdn];ssj1=lj1[tdn];
 	ssi2=li2[tdn];ssj1=lj2[tdn];
 	//printf("----%d %d %d\n",xxi1,xxj1,tdn);
 	for(i=1;i<N0;i++)
 	{
 		for(j=i+1;j<=N0;j++)
 		{
 				ccs[i][j]=cs[tdn][i][j];
 				aas[i]=as[tdn][i];
 				aas[j]=as[tdn][j];
 		}
 	}
 	for(i=1;i<=N0;i++) 
  	{
  		tertiary[i]=ter[tdn][i];
            	x[i]=rx[tdn][i];
            	y[i]=ry[tdn][i];
            	z[i]=rz[tdn][i];
            	type[i]=rtype[tdn][i];
            	Q[i]=rQ[tdn][i];
            	f[i]=rf[tdn][i];
            	R[i]=rR[tdn][i];
            	if(xxi1==0)
            	{
            		junction[tdn]=1;
            	}
            	else if(xxi1!=0&&xxi2==0)
            	{
            		junction[tdn]=2;
            	}
            	else if(xxi3!=0&&xxi4==0&&xxi5==0)
            	{
            		int loop1,loop2,loop3;
            		int st2,st3;
            		st2=xxj2-xxi2;st3=xxj3-xxi3;
            		loop1=xxi2-xxi1;loop2=xxi3-xxj2;loop3=xxj1-xxj3;
            		if(loop2>=12&&loop2>loop1&&loop2>loop3)	
            		{
            			if((st2>2*st3)||(st3>2*st2))	{condition[tdn]=0;}
            			else				{condition[tdn]=1;}
            		}
            		else					{condition[tdn]=0;}
			junction[tdn]=3;
			if(condition[tdn]!=2)
			{
				if((i>xxi1&&i<xxi2)||(i>xxj2&&i<xxi3)||(i>xxj3&&i<xxj1))
            			{
            				type[i]='X';	
            			}
			}
            	}
            	else if(xxi4!=0&&xxi5==0)
            	{
            		if((i>xxi1&&i<xxi2)||(i>xxj2&&i<xxi3)||(i>xxj3&&i<xxi4)||(i>xxj4&&i<xxj1))
            		{
            			type[i]='X';	
            		}
            		int loop1,loop2;
            		loop1=xxi3-xxj2;loop2=xxi4-xxj3;
            		if(loop1>12||loop2>12)	{condition[tdn]=1;}
            		else			{condition[tdn]=0;}
			junction[tdn]=4;
            	}
            	else
            	{
            		if((i>xxi1&&i<xxi2)||(i>xxj2&&i<xxi3)||(i>xxj3&&i<xxi4)||(i>xxj4&&i<xxi5)||(i>xxj5&&i<xxj1))
            		{
            			type[i]='X';	
            		}
			junction[tdn]=5;
            	}
            	//printf("%c",type[i]); 
  	} 
  	//printf("\n");
  	for(i=1;i<=N0;i++)
      	{
      		Apcp[tdn][i]=1.6;Acpc[tdn][i]=1.6;ApcN[tdn][i]=1.64;ANcp[tdn][i]=1.66;
      		dpcpc[tdn][i]=2.56; dcpcp[tdn][i]=-2.94; dcpcN[tdn][i]=-1.16; dNcpc[tdn][i]=0.88; 
        }
        if(junction[tdn]==2)
        {
       		dpcpc[tdn][xxi1]=2.90; dcpcp[tdn][xxi1]=2.80; 
	    	dpcpc[tdn][xxi1+3]=2.90; dcpcp[tdn][xxi1+3]=2.70; 
	     	dpcpc[tdn][xxi1+6]=-1.80; dcpcp[tdn][xxi1+6]=-3.10; 
        }
        if(junction[tdn]==4&&condition[tdn]==1)
        {
        	dpcpc[tdn][xxi1]=0.73; dcpcp[tdn][xxi1]=-2.87; dcpcN[tdn][xxi1]=-1.41; dNcpc[tdn][xxi1]=-1.12; 
	    	dpcpc[tdn][xxi1+3]=-2.78; dcpcp[tdn][xxi1+3]=2.40; dcpcN[tdn][xxi1+3]=-1.26; dNcpc[tdn][xxi1+3]=1.82; 
	     	dpcpc[tdn][xxi1+6]=0.68; dcpcp[tdn][xxi1+6]=-0.50; dcpcN[tdn][xxi1+6]=-1.80; dNcpc[tdn][xxi1+6]=-1.13; 
        }
       	if(junction[tdn]==5)
        {
             Apcp[tdn][xxi1]=1.56;Acpc[tdn][xxi1]=1.44;ApcN[tdn][xxi1]=1.64;ANcp[tdn][xxi1]=1.80; dpcpc[tdn][xxi1]=-1.12; dcpcp[tdn][xxi1]=-0.02; dcpcN[tdn][xxi1]=-1.19; dNcpc[tdn][xxi1]=-2.77; 
	     Apcp[tdn][xxi1+3]=2.44;Acpc[tdn][xxi1+3]=1.63;ApcN[tdn][xxi1+3]=1.79;ANcp[tdn][xxi1+3]=2.04; dpcpc[tdn][xxi1+3]=-1.97; dcpcp[tdn][xxi1+3]=-2.25; dcpcN[tdn][xxi1+3]=3.01; dNcpc[tdn][xxi1+3]=1.30; 
	     Apcp[tdn][xxi1+6]=2.16;Acpc[tdn][xxi1+6]=1.85;ApcN[tdn][xxi1+6]=1.91;ANcp[tdn][xxi1+6]=1.74; dpcpc[tdn][xxi1+6]=-1.96; dcpcp[tdn][xxi1+6]=0.80; dcpcN[tdn][xxi1+6]=-0.23; dNcpc[tdn][xxi1+6]=2.22; 
	     Apcp[tdn][xxi2]=1.55;Acpc[tdn][xxi2]=1.94;ApcN[tdn][xxi2]=1.49;ANcp[tdn][xxi2]=1.68; dpcpc[tdn][xxi2]=2.30; dcpcp[tdn][xxi2]=-2.90; dcpcN[tdn][xxi2]=2.49; dNcpc[tdn][xxi2]=0.81; 
	     Apcp[tdn][xxj2]=1.72;Acpc[tdn][xxj2]=1.80;ApcN[tdn][xxj2]=1.73;ANcp[tdn][xxj2]=1.57; dpcpc[tdn][xxj2]=2.63; dcpcp[tdn][xxj2]=-3.05; dcpcN[tdn][xxj2]=-1.24; dNcpc[tdn][xxj2]=0.90; 
	     Apcp[tdn][xxi3]=1.64;Acpc[tdn][xxi3]=1.72;ApcN[tdn][xxi3]=1.67;ANcp[tdn][xxi3]=1.65; dpcpc[tdn][xxi3]=2.61; dcpcp[tdn][xxi3]=-2.93; dcpcN[tdn][xxi3]=-1.40; dNcpc[tdn][xxi3]=0.94; 
	     Apcp[tdn][xxj3]=2.30;Acpc[tdn][xxj3]=2.38;ApcN[tdn][xxj3]=1.60;ANcp[tdn][xxj3]=2.02; dpcpc[tdn][xxj3]=1.58; dcpcp[tdn][xxj3]=-3.06; dcpcN[tdn][xxj3]=-1.51; dNcpc[tdn][xxj3]=-0.50; 
	     Apcp[tdn][xxi4]=1.54;Acpc[tdn][xxi4]=1.84;ApcN[tdn][xxi4]=1.52;ANcp[tdn][xxi4]=1.66; dpcpc[tdn][xxi4]=2.47; dcpcp[tdn][xxi4]=-2.87; dcpcN[tdn][xxi4]=-1.40; dNcpc[tdn][xxi4]=0.94; 
	     Apcp[tdn][xxj4]=2.20;Acpc[tdn][xxj4]=1.53;ApcN[tdn][xxj4]=1.49;ANcp[tdn][xxj4]=1.67; dpcpc[tdn][xxj4]=-1.91; dcpcp[tdn][xxj4]=0.79; dcpcN[tdn][xxj4]=-1.36; dNcpc[tdn][xxj4]=2.68; 
	     Apcp[tdn][xxi5]=1.80;Acpc[tdn][xxi5]=1.88;ApcN[tdn][xxi5]=1.48;ANcp[tdn][xxi5]=1.67; dpcpc[tdn][xxi5]=2.38; dcpcp[tdn][xxi5]=-2.92; dcpcN[tdn][xxi5]=2.44; dNcpc[tdn][xxi5]=0.88; 
	     Apcp[tdn][xxj5]=1.89;Acpc[tdn][xxj5]=1.89;ApcN[tdn][xxj5]=1.68;ANcp[tdn][xxj5]=1.60; dpcpc[tdn][xxj5]=2.52; dcpcp[tdn][xxj5]=-2.99; dcpcN[tdn][xxj5]=-1.23; dNcpc[tdn][xxj5]=0.83; 
	     Apcp[tdn][xxj1]=1.60;Acpc[tdn][xxj1]=1.79;ApcN[tdn][xxj1]=1.57;ANcp[tdn][xxj1]=1.66; dpcpc[tdn][xxj1]=2.58; dcpcp[tdn][xxj1]=-2.90; dcpcN[tdn][xxj1]=-1.32; dNcpc[tdn][xxj1]=1.00; 	
       }
       if(junction[tdn]==3)
       {
       		for(i=ssi1;i<=ssj1;i++)
       		{
       			for(j=ssi2;j<=ssj2;j++)
       			{
       				if(ccs[i][j]==1)	{condition[tdn]=1;}
       			}
       		}
       	}
       if(N<=40&&junction[tdn]==1)	{NS1=1.0;NS2=1.0;}
       else				{NS1=0.4;NS2=1.6;}
  	MC_Annealing();
  	for(i=1;i<=N0;i++)
  	{
            	rx[tdn][i]=x[i];
            	ry[tdn][i]=y[i];
            	rz[tdn][i]=z[i];
  	}
  	RU[tdn]=U;
 // rt0[tdn]=t0;
}  



void openfile()
{
     	int i,duo1,duo2,j,in;
     	int xi1,xj1,xi2,xj2,xi3,xj3,xi4,xj4,xi5,xj5;
     	int si1,sj1,si2,sj2;
     	FILE *fp;
     	char filename[20];
     	int ai,aj;
     	char ax,ay;
     	Infor=fopen("Information.dat","w+");
     	fp=fopen("ch.dat","r+");
     	fps=fopen("cs1.dat","r+");
    	fpb=fopen("stem.dat","r+");
     	fpb1=fopen("stem_kissing.dat","r+");
     	fpb2=fopen("config1.dat","r+");
      	i=1;
     	while(!feof(fp)) 
     	{
         	fscanf(fp,"%d %d %s %f %f %f %f %f %f\n",&duo1,&duo2,&rtype0[i],&rx0[i],&ry0[i],&rz0[i],&rR0[i],&rQ0[i],&rf0[i]); 
         	i++;
     	}      
     	fclose(fp); 
     	rp=i-1; 
     	for(i=1;i<rp;i++)
     	{
          	for(j=i+1;j<=rp;j++)
          	{
               		for(in=0;in<Nthd;in++)
               		{
                    		cs[in][i][j]=0;
                    		as[in][i]=0;
                    		as[in][j]=0;
                    		ter[in][i]=0;
                    		ter[in][j]=0;
               		}
          	}
     	}
     	int k=0;
     	int ad;
     	while(!feof(fps))
     	{
            	fscanf(fps,"%d %c %d %c %d\n",&ai,&ax,&aj,&ay,&ad);
            	k++;
            	for(i=0;i<Nthd;i++)
            	{
                	cs[i][ai*3][aj*3]=1;
                	as[i][ai*3]=1;
                	as[i][aj*3]=1;
                	ter[i][ai*3]=ad;
                	ter[i][aj*3]=ad;
                 //printf("%d\n",as[0][3]);
            	}
     	}
     	bp_all=k;                // printf("%d\n",cs[0][ai*3][aj*3])}
     	fclose(fps); 
     	while(!feof(fpb))
     	{
           	fscanf(fpb,"%d %d %d %d %d %d %d %d %d %d\n",&xi1,&xj1,&xi2,&xj2,&xi3,&xj3,&xi4,&xj4,&xi5,&xj5);
     	}
     	fclose(fpb);     	   
     	while(!feof(fpb1))
     	{
           	fscanf(fpb1,"%d %d %d %d\n",&si1,&sj1,&si2,&sj2);
     	}
     	fclose(fpb1);
     	for(i=0;i<Nthd;i++)
     	{
     		yi1[i]=xi1*3;yj1[i]=xj1*3;
     		yi2[i]=xi2*3;yj2[i]=xj2*3;
     		yi3[i]=xi3*3;yj3[i]=xj3*3;
     		yi4[i]=xi4*3;yj4[i]=xj4*3;
     		yi5[i]=xi5*3;yj5[i]=xj5*3;
     		li1[i]=si1*3;lj1[i]=sj1*3;
     		li2[i]=si2*3;lj2[i]=sj2*3;
     	}
     	//printf("----------%d %d\n",xi1[0],xj1[0]);
     	//printf("CG beads number:   %d,  the thread number:   %d\n",rp,Nthd);
     	fprintf(Infor,"CG beads number:   %d,  the thread number:   %d\n",rp,Nthd);fflush(Infor);
     	int ca,cb,cc,cd,ce,ct;
     	while(!feof(fpb2))
     	{
     		fscanf(fpb2,"%d %d %d %d %d %d\n",&ct,&ca,&cb,&cc,&cd,&ce);
     	}
     	fclose(fpb2);
     	for (i=0;i<Nthd;i++) 
     	{
           	sprintf(filename,"conf_%d.dat",i); fp_conf[i]=fopen(filename,"w+");
           	sprintf(filename,"secondary_%d.dat",i); fp_sec_stru[i]=fopen(filename,"w+");
         	sprintf(filename,"energy_%d.dat",i);    fp_energy[i]=fopen(filename,"w+");
      	}         
}


void closefile()
{
         int i;
         for(i=0;i<Nthd;i++)
         {
                  fclose(fp_conf[i]); fclose(fp_sec_stru[i]); fclose(fp_energy[i]);
         }
         fclose(Infor);
}

/* &%$#@!~&%$#@!~&%$#@!~   Some functions or modules for move and calculation    &%$#@!~&%$#@!~&%$#@!~ */
/*********************读入文件，读出文件，输入参数************************/
void Put_File (void) 
{
 
   	CMg=0.0;
   	CNa=0.0;
   	salt=0;
   	t0=rt0[tdn];
   	//printf("thread:  %d,  CNa = %f,  CMg = %f, t0 = %f\n",tdn,CNa,CMg,t0);
   	fprintf(Infor,"thread:  %d,  CNa = %f,  CMg = %f, t0 = %f\n",tdn,CNa,CMg,t0);fflush(Infor);
}
/*********************************************/
void Fixed_Atom(void)
{
   	int N10;
   	N10=floor(N0/2)+1; 
   	if (fmod(N10,3)==0) 		{N1=N10-2;}
   	else if (fmod(N10+1,3)==0) 	{N1=N10-1;}
   	else 				{N1=N10;}                         //N1:Fixed the P in centre of chain
 //printf("N %d N0 %d N1 %d\n",N,N0,N1);  //The number of nucleotides & atoms & the actionless atom
}
//******************************************//
void Parameters_T()
{
    	float tt;
    	total=total1;        //Steps in Folding annealing
    	tt=25.0; 
    	T=273.15+tt*1.0; 
    	D=T*2.0*pow(10,-3);   
    	float qq4,qqq=0.0;
    	fNa=0.001*CNa/(0.001*CNa+(8.1-32.4/(N*0.5))*(5.2-log(0.001*CNa))*0.001*CMg); //The percentage of Na+ in Mixture in TBI_Helix
    	Ek=87.740-0.4008*tt+9.398*1e-4*tt*tt-1.41*1e-6*tt*tt*tt;  //Permittivity
    	I=CNa+IMg*CMg; 
    	Kd=sqrt((0.396*Ek*T)/I);                   //Ionic strength & Debye length
    	qq4=5.998*1e-6*bl*Ek*T*0.5*(fNa+1);
    	if (N*5.5<Kd) 
    	{
            	qqq=log(Kd/bl)/log(N); 
            	if (qqq>1.)
            	{ 
                    	q4=qq4*qqq; 
            	}  
            	else 
            	{
                     	q4=qq4;
            	}
     	}
     	else 
     	{
           	q4=qq4;
    	} // q4=b/lB &稀溶液修正
}
/****************Monte Carlo simulated annealing*********************/
void MC_Annealing()
{
   	void MC_T(),Parameters_T(); 
     	l=0;                 //MC steps independent of T
     	Parameters_T();     //parameters at any t0: steps, T, Debye length, ionic strength, ion fraction etc.
     	MC_T();
}
/*****************Monte Carlo simulation at given Temperature******************************/
void MC_T(void)
{
     	int ii,i;
     	void MC_Each_Step(),StructureInformation(),ENERGY(),disbrute_qq();
     	srand((unsigned)time(NULL)); 
     	for (t=1;t<=total1;t++)
     	{    
           	if(salt==1)
           	{
           		if(fmod(t,50)==0||t==1)
           		{
                   		disbrute_qq();             
             		}
             	} 
           	for(ii=1;ii<=N0;ii++)
           	{
           		cycle[tdn]=ii;
                  	MC_Each_Step();
           	}
           	if (fmod(t,tconf1)==0&&t>(total1*3/4)) 
           	{                   
                	for(i=1;i<=N0;i++) 
                	{
                      		fprintf(fp_conf[tdn],"%d %d %c %f %f %f %f %f %f\n",t,i,rtype[tdn][i],x[i],y[i],z[i],R[i],Q[i],f[i]);
                	}
                	fflush(fp_conf[tdn]);                           
            	}              
            	ENERGY();                //Calculate energy of one conformation        	                          
                 //End of the step t cycle;
      		if(fmod(t,1000)==0)
      		{ 
             		//printf("thread: %d,  step:%6d,     T:%5.1f,    Energy: %f kcal/mol\n",tdn,t,t0,U);
             		fprintf(Infor,"thread: %d,  step:%6d,     T:%5.1f,    Energy: %f kcal/mol\n",tdn,t,t0,U);fflush(Infor); 
      		}   
    	}   
}
/*****************Monte Carlo for each step******************************/
void MC_Each_Step()
{
 	int   i0;
 	float PC(),CP(),CN(),PCP(),CPC(),PCN(),NCP(),PCPC(),CPCP(),PCPCh(),CPCPh(),CPCN(),NCPC(),LJ0(),HB(),St(),QQ(),QQ11(),uys(),HB1();
 	void MoveN(),Rand01(),Translate(),Pivot(),RMSDconf(),FOLD();
 	void BasePairing(),BaseStacking(),ExcludedVolume(),Electrostatic(),Bonded(),CoaxialStacking(),constraint();
 	void Metropolis(); 
   	u=0.0;uu=0.0;                                    //Initialize the energies;
   	ulj=0.0;ub=0.0;ue=0.0;ud=0.0;uN=0.0;us=0.0;uc=0.0;
   	uulj=0.0;uub=0.0;uue=0.0;uud=0.0;uuN=0.0;uus=0.0;uuc=0.0;
   	uco=0.0;uuco=0.0;ucst=0.0;uucst=0.0;
   	Energy=0;                      //if Energy=1, fuction ENERGY() is running 
   	Move=0;     //if Move_3nt=1, 3nt fragment moves;

       /*****************固定中心原子*********************/
   	re:;
   	rand01=rand()/(RAND_MAX+1.);
   	i0=floor(rand01*N0)+1; 
   	if (i0==N1) goto re;      
       /***********如果Folding=0，则执行优化，否则执行折叠***********/
   	FOLD(i0);                 //Folding Process
   	Metropolis();
}
/*********************************************/
void Metropolis(void)
{
   	int i;
   	float p;
   	u=ulj+uN+us+uc+(ub*NS1+ue*NS1+ud*NS2)*0.5963+uco+ucst/*+uTri+umis+udang*/;
   	uu=uulj+uuN+uus+uuc+(uub*NS1+uue*NS1+uud*NS2)*0.5963+uuco+uucst/*+uuTri+uumis+uudang*/;
   	du=uu-u; p=0.0;  //du: Energy changes before & after moves;
   	if(du<=0.0)   
   	{
             	for (i=1;i<=N0;i++) 
              	{
                     	x[i]=xx[i];y[i]=yy[i];z[i]=zz[i];
              	}   
    	}  //To update the conf.
    	else
    	{
              	rand01=rand()/(RAND_MAX+1.);  
              	p=exp((-1)*du/D);
              	if(rand01<=p)   
              	{
                  	for (i=1;i<=N0;i++) 
                   	{
                                        x[i]=xx[i];y[i]=yy[i];z[i]=zz[i];
                        }  
              	}
              	else            
              	{
                     	for (i=1;i<=N0;i++) 
                      	{
                     		x[i]=x[i];y[i]=y[i];z[i]=z[i];
                   	}  
              	}
     	}
}
//********************  Move  *******************//
void Rand01()       //Generate rodom Euler angle;
{      
  	rand01=rand()/(RAND_MAX+1.);phi=rand01*pi;  
   	rand01=rand()/(RAND_MAX+1.);theta=rand01*2.*pi; 
	rand01=rand()/(RAND_MAX+1.);rm=rand01*2.*pi;
    	rand01=rand()/(RAND_MAX+1.);
        
}
void Translate(int i1) //Translation of one atom;
{
       xxx[i1]=x[i1]+0.5*rand01*sin(phi)*cos(theta);
       yyy[i1]=y[i1]+0.5*rand01*sin(phi)*sin(theta);
       zzz[i1]=z[i1]+0.5*rand01*cos(phi);
}

void Pivot(int i1,int i10)  //Pivot moves for one segment;
{
 	xx[i1]=(xxx[i1]-xxx[i10])*(cos(theta)*cos(rm)-cos(phi)*sin(theta)*sin(rm))+(yyy[i1]-yyy[i10])*(sin(theta)*cos(rm)+cos(phi)*cos(theta)*sin(rm))+(zzz[i1]-zzz[i10])*sin(phi)*sin(rm)+xxx[i10];
     	yy[i1]=-(xxx[i1]-xxx[i10])*(cos(theta)*sin(rm)+cos(phi)*sin(theta)*cos(rm))+(yyy[i1]-yyy[i10])*(cos(phi)*cos(theta)*cos(rm)-sin(theta)*sin(rm))+(zzz[i1]-zzz[i10])*sin(phi)*cos(rm)+yyy[i10];
    	zz[i1]=(xxx[i1]-xxx[i10])*sin(phi)*sin(theta)-sin(phi)*cos(theta)*(yyy[i1]-yyy[i10])+(zzz[i1]-zzz[i10])*cos(phi)+zzz[i10];
}
void MoveN(int i1) //Movement of each base;
{
        Rand01(); Translate(i1);
        Rand01(); Pivot(i1,i1-1);
}

/* ~~~~~~~~~~ Details of bonded potential calculation ~~~~~~~~~~~~~~ */
float PC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
	float d,ul;
     	d=sqrt((x1[i1-2]-x1[i1-1])*(x1[i1-2]-x1[i1-1])+(y1[i1-2]-y1[i1-1])*(y1[i1-2]-y1[i1-1])+(z1[i1-2]-z1[i1-1])*(z1[i1-2]-z1[i1-1]));
     	ul=68.0*(d-3.8)*(d-3.8);
     	return ul;
}
float CP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d,ul;
     	d=sqrt((x1[i1+1]-x1[i1-1])*(x1[i1+1]-x1[i1-1])+(y1[i1+1]-y1[i1-1])*(y1[i1+1]-y1[i1-1])+(z1[i1+1]-z1[i1-1])*(z1[i1+1]-z1[i1-1]));
     	ul=33.0*(d-3.8)*(d-3.8);
     	return ul;
}
float CN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d,ul;
     	d=sqrt((x1[i1]-x1[i1-1])*(x1[i1]-x1[i1-1])+(y1[i1]-y1[i1-1])*(y1[i1]-y1[i1-1])+(z1[i1]-z1[i1-1])*(z1[i1]-z1[i1-1]));
     	ul=24.8*(d-3.36)*(d-3.36);
     	return ul;
}
float PCP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
      	float d1,d2,d3,w,a1,ue0;
      	d1=sqrt((x1[i1-2]-x1[i1-1])*(x1[i1-2]-x1[i1-1])+(y1[i1-2]-y1[i1-1])*(y1[i1-2]-y1[i1-1])+(z1[i1-2]-z1[i1-1])*(z1[i1-2]-z1[i1-1]));
      	d2=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
      	d3=sqrt((x1[i1+1]-x1[i1-2])*(x1[i1+1]-x1[i1-2])+(y1[i1+1]-y1[i1-2])*(y1[i1+1]-y1[i1-2])+(z1[i1+1]-z1[i1-2])*(z1[i1+1]-z1[i1-2]));
      	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
      	if (w<=-1.0) 
      	{
          	a1=3.14;
      	}
      	else if (w>=1.0) 
      	{
           	a1=0.;
      	}
      	else  
      	{
           	a1=acos(w);
      	} 
      	ue0=9.3*(a1-Apcp[tdn][i1])*(a1-Apcp[tdn][i1]);
      	return ue0;
}
float CPC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
      	float d1,d2,d3,w,a1,ue0;
        d1=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
        d2=sqrt((x1[i1+1]-x1[i1+2])*(x1[i1+1]-x1[i1+2])+(y1[i1+1]-y1[i1+2])*(y1[i1+1]-y1[i1+2])+(z1[i1+1]-z1[i1+2])*(z1[i1+1]-z1[i1+2]));
        d3=sqrt((x1[i1+2]-x1[i1-1])*(x1[i1+2]-x1[i1-1])+(y1[i1+2]-y1[i1-1])*(y1[i1+2]-y1[i1-1])+(z1[i1+2]-z1[i1-1])*(z1[i1+2]-z1[i1-1]));
        w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
      	if (w<=-1.0) 
       	{
          	a1=3.14;
     	}
      	else if (w>=1.0) 
      	{
            	a1=0.;
      	}
       	else  
      	{
           	a1=acos(w);
      	}
      	ue0=15.3*(a1-Acpc[tdn][i1])*(a1-Acpc[tdn][i1]);
      	return ue0;
}
float PCN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
      	float d1,d2,d3,w,a1,ue0;
      	d1=sqrt((x1[i1-1]-x1[i1-2])*(x1[i1-1]-x1[i1-2])+(y1[i1-1]-y1[i1-2])*(y1[i1-1]-y1[i1-2])+(z1[i1-1]-z1[i1-2])*(z1[i1-1]-z1[i1-2]));
      	d2=sqrt((x1[i1-1]-x1[i1])*(x1[i1-1]-x1[i1])+(y1[i1-1]-y1[i1])*(y1[i1-1]-y1[i1])+(z1[i1-1]-z1[i1])*(z1[i1-1]-z1[i1]));
      	d3=sqrt((x1[i1-2]-x1[i1])*(x1[i1-2]-x1[i1])+(y1[i1-2]-y1[i1])*(y1[i1-2]-y1[i1])+(z1[i1-2]-z1[i1])*(z1[i1-2]-z1[i1]));
      	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
      	if (w<=-1.0) {a1=3.14;}
      	else if (w>=1.0) {a1=0.;}
      	else  {a1=acos(w);}
      	ue0=9.8*(a1-ApcN[tdn][i1])*(a1-ApcN[tdn][i1]);
      	return ue0;
}
float NCP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d1,d2,d3,w,a1,ue0;
     	d1=sqrt((x1[i1-1]-x1[i1])*(x1[i1-1]-x1[i1])+(y1[i1-1]-y1[i1])*(y1[i1-1]-y1[i1])+(z1[i1-1]-z1[i1])*(z1[i1-1]-z1[i1]));
     	d2=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
     	d3=sqrt((x1[i1+1]-x1[i1])*(x1[i1+1]-x1[i1])+(y1[i1+1]-y1[i1])*(y1[i1+1]-y1[i1])+(z1[i1+1]-z1[i1])*(z1[i1+1]-z1[i1]));
     	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
     	if (w<=-1.0) {a1=3.14;}
     	else if (w>=1.0) {a1=0.;}
     	else  {a1=acos(w);}
     	ue0=15.2*(a1-ANcp[tdn][i1])*(a1-ANcp[tdn][i1]);
     	return ue0;
}
float PCPC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
 	c1=((y1[i1-2]-y1[i1-1])*(z1[i1-1]-z1[i1+1])-(z1[i1-2]-z1[i1-1])*(y1[i1-1]-y1[i1+1]));
 	c2=((z1[i1-2]-z1[i1-1])*(x1[i1-1]-x1[i1+1])-(x1[i1-2]-x1[i1-1])*(z1[i1-1]-z1[i1+1]));
 	c3=((x1[i1-2]-x1[i1-1])*(y1[i1-1]-y1[i1+1])-(y1[i1-2]-y1[i1-1])*(x1[i1-1]-x1[i1+1]));
 	p1=((y1[i1-1]-y1[i1+1])*(z1[i1+1]-z1[i1+2])-(z1[i1-1]-z1[i1+1])*(y1[i1+1]-y1[i1+2]));
 	p2=((z1[i1-1]-z1[i1+1])*(x1[i1+1]-x1[i1+2])-(x1[i1-1]-x1[i1+1])*(z1[i1+1]-z1[i1+2]));
 	p3=((x1[i1-1]-x1[i1+1])*(y1[i1+1]-y1[i1+2])-(y1[i1-1]-y1[i1+1])*(x1[i1+1]-x1[i1+2]));
 	e1=sqrt(c1*c1+c2*c2+c3*c3); 
 	f1=sqrt(p1*p1+p2*p2+p3*p3);
 	pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 	g1=(x1[i1-2]-x1[i1+2]); 
 	g2=(y1[i1-2]-y1[i1+2]); 
 	g3=(z1[i1-2]-z1[i1+2]);
 	gg1=sqrt(g1*g1+g2*g2+g3*g3); 
 	hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 	if (pp1<=-1.0) {di=-3.14;}
 	else if (pp1>=1.0) {di=0.;}
 	else if (hh1>=0.) {di=acos(pp1);}
 	else {di=-acos(pp1);}
 	if(aas[i1]==1) 	 	  {ud0=2.7*((1-cos(di-2.51))+0.5*(1-cos(3.*(di-2.51))));}
 	else                      {ud0=1.2*((1-cos(di-dpcpc[tdn][i1]))+0.5*(1-cos(3.*(di-dpcpc[tdn][i1]))));}
 	return ud0;
}
float CPCP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
   	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
     	c1=((y1[i1-1]-y1[i1+1])*(z1[i1+1]-z1[i1+2])-(z1[i1-1]-z1[i1+1])*(y1[i1+1]-y1[i1+2]));
     	c2=((z1[i1-1]-z1[i1+1])*(x1[i1+1]-x1[i1+2])-(x1[i1-1]-x1[i1+1])*(z1[i1+1]-z1[i1+2]));
     	c3=((x1[i1-1]-x1[i1+1])*(y1[i1+1]-y1[i1+2])-(y1[i1-1]-y1[i1+1])*(x1[i1+1]-x1[i1+2]));
     	p1=((y1[i1+1]-y1[i1+2])*(z1[i1+2]-z1[i1+4])-(z1[i1+1]-z1[i1+2])*(y1[i1+2]-y1[i1+4]));
     	p2=((z1[i1+1]-z1[i1+2])*(x1[i1+2]-x1[i1+4])-(x1[i1+1]-x1[i1+2])*(z1[i1+2]-z1[i1+4]));
     	p3=((x1[i1+1]-x1[i1+2])*(y1[i1+2]-y1[i1+4])-(y1[i1+1]-y1[i1+2])*(x1[i1+2]-x1[i1+4]));
     	e1=sqrt(c1*c1+c2*c2+c3*c3); 
     	f1=sqrt(p1*p1+p2*p2+p3*p3);
     	pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
     	g1=(x1[i1-1]-x1[i1+4]); 
     	g2=(y1[i1-1]-y1[i1+4]); 
     	g3=(z1[i1-1]-z1[i1+4]);
     	gg1=sqrt(g1*g1+g2*g2+g3*g3); 
     	hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
     	if (pp1<=-1.0) {di=-3.14;}
     	else if (pp1>=1.0) {di=0.;}
     	else if (hh1>=0.) {di=acos(pp1);}
     	else {di=-acos(pp1);}
     	if(aas[i1]==1)       	 {ud0=10.5*((1-cos(di+2.92))+0.5*(1-cos(3.*(di+2.92))));}
 	else                     {ud0=5.3*((1-cos(di-dcpcp[tdn][i1]))+0.5*(1-cos(3.*(di-dcpcp[tdn][i1]))));}
   	return ud0;
}
float CPCN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
  	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
     	c1=((y1[i1-4]-y1[i1-2])*(z1[i1-2]-z1[i1-1])-(z1[i1-4]-z1[i1-2])*(y1[i1-2]-y1[i1-1]));
     	c2=((z1[i1-4]-z1[i1-2])*(x1[i1-2]-x1[i1-1])-(x1[i1-4]-x1[i1-2])*(z1[i1-2]-z1[i1-1]));
     	c3=((x1[i1-4]-x1[i1-2])*(y1[i1-2]-y1[i1-1])-(y1[i1-4]-y1[i1-2])*(x1[i1-2]-x1[i1-1]));
     	p1=((y1[i1-2]-y1[i1-1])*(z1[i1-1]-z1[i1])-(z1[i1-2]-z1[i1-1])*(y1[i1-1]-y1[i1]));
     	p2=((z1[i1-2]-z1[i1-1])*(x1[i1-1]-x1[i1])-(x1[i1-2]-x1[i1-1])*(z1[i1-1]-z1[i1]));
     	p3=((x1[i1-2]-x1[i1-1])*(y1[i1-1]-y1[i1])-(y1[i1-2]-y1[i1-1])*(x1[i1-1]-x1[i1]));
     	e1=sqrt(c1*c1+c2*c2+c3*c3); 
     	f1=sqrt(p1*p1+p2*p2+p3*p3);
     	pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
     	g1=(x1[i1-4]-x1[i1]); 
     	g2=(y1[i1-4]-y1[i1]); 
     	g3=(z1[i1-4]-z1[i1]);
     	gg1=sqrt(g1*g1+g2*g2+g3*g3); 
     	hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
     	if (pp1<=-1.0) {di=-3.14;}
     	else if (pp1>=1.0) {di=0.;}
     	else if (hh1>=0.) {di=acos(pp1);}
     	else {di=-acos(pp1);}
     	ud0=0.825*((1-cos(di-dcpcN[tdn][i1]))+0.5*(1-cos(3.*(di-dcpcN[tdn][i1]))));
   	return ud0;
}
float NCPC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
   	c1=((y1[i1]-y1[i1-1])*(z1[i1-1]-z1[i1+1])-(z1[i1]-z1[i1-1])*(y1[i1-1]-y1[i1+1]));
   	c2=((z1[i1]-z1[i1-1])*(x1[i1-1]-x1[i1+1])-(x1[i1]-x1[i1-1])*(z1[i1-1]-z1[i1+1]));
   	c3=((x1[i1]-x1[i1-1])*(y1[i1-1]-y1[i1+1])-(y1[i1]-y1[i1-1])*(x1[i1-1]-x1[i1+1]));
   	p1=((y1[i1-1]-y1[i1+1])*(z1[i1+1]-z1[i1+2])-(z1[i1-1]-z1[i1+1])*(y1[i1+1]-y1[i1+2]));
   	p2=((z1[i1-1]-z1[i1+1])*(x1[i1+1]-x1[i1+2])-(x1[i1-1]-x1[i1+1])*(z1[i1+1]-z1[i1+2]));
   	p3=((x1[i1-1]-x1[i1+1])*(y1[i1+1]-y1[i1+2])-(y1[i1-1]-y1[i1+1])*(x1[i1+1]-x1[i1+2]));
   	e1=sqrt(c1*c1+c2*c2+c3*c3);
   	f1=sqrt(p1*p1+p2*p2+p3*p3);
   	pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
   	g1=(x1[i1]-x1[i1+2]); 
   	g2=(y1[i1]-y1[i1+2]); 
   	g3=(z1[i1]-z1[i1+2]);
   	gg1=sqrt(g1*g1+g2*g2+g3*g3); 
   	hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
   	if (pp1<=-1.0) 
   	{
         	di=-3.14;
   	}
   	else if (pp1>=1.0) 
   	{
        	di=0.;
   	}
   	else if (hh1>=0.) 
   	{
        	di=acos(pp1);
   	}
   	else 
   	{
        	di=-acos(pp1);
   	}
   	ud0=0.78*((1-cos(di-dNcpc[tdn][i1]))+0.5*(1-cos(3.*(di-dNcpc[tdn][i1]))));
   	return ud0;
}
/*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
//base pairing between two complementary bases (AU,GC,and GU)
float HB(int i1,int j1,float x1[10000],float y1[10000],float z1[10000])
{
   	float d,hb,UHB,d0,d01,d1,d11;
   	d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); //NN
   	if (d>=8.55&&d<=9.39)   //The pairing formation condition
   	{
        	d0=sqrt((x1[i1]-x1[j1-1])*(x1[i1]-x1[j1-1])+(y1[i1]-y1[j1-1])*(y1[i1]-y1[j1-1])+(z1[i1]-z1[j1-1])*(z1[i1]-z1[j1-1]));  //NiCj
        	d01=sqrt((x1[j1]-x1[i1-1])*(x1[j1]-x1[i1-1])+(y1[j1]-y1[i1-1])*(y1[j1]-y1[i1-1])+(z1[j1]-z1[i1-1])*(z1[j1]-z1[i1-1])); //CiNj
        	d1=sqrt((x1[i1-2]-x1[j1])*(x1[i1-2]-x1[j1])+(y1[i1-2]-y1[j1])*(y1[i1-2]-y1[j1])+(z1[i1-2]-z1[j1])*(z1[i1-2]-z1[j1]));  //PiNj
        	d11=sqrt((x1[j1-2]-x1[i1])*(x1[j1-2]-x1[i1])+(y1[j1-2]-y1[i1])*(y1[j1-2]-y1[i1])+(z1[j1-2]-z1[i1])*(z1[j1-2]-z1[i1])); //NiPj
        	if ((type[i1]=='G'&&type[j1]=='C')||(type[i1]=='C'&&type[j1]=='G'))  {hb=-3.325;}
        	if ((type[i1]=='G'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='G'))  {hb=-2.10;}
        	if ((type[i1]=='A'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='A'))  {hb=-2.10;}
//UHB=hb/(1+3.6*(d-8.94)*(d-8.94)+1.9*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.7*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));
        	UHB=hb/(1+2.66*(d-8.94)*(d-8.94)+1.37*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.464*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));
    	}
    	else UHB=0.0;
    	return UHB;
}
float HB1(int i1,int j1,float x1[10000],float y1[10000],float z1[10000])
{
   	float d,hb,UHB,d0,d01,d1,d11;
        d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); //NN
        d0=sqrt((x1[i1]-x1[j1-1])*(x1[i1]-x1[j1-1])+(y1[i1]-y1[j1-1])*(y1[i1]-y1[j1-1])+(z1[i1]-z1[j1-1])*(z1[i1]-z1[j1-1]));  //NiCj
        d01=sqrt((x1[j1]-x1[i1-1])*(x1[j1]-x1[i1-1])+(y1[j1]-y1[i1-1])*(y1[j1]-y1[i1-1])+(z1[j1]-z1[i1-1])*(z1[j1]-z1[i1-1])); //CiNj
        d1=sqrt((x1[i1-2]-x1[j1])*(x1[i1-2]-x1[j1])+(y1[i1-2]-y1[j1])*(y1[i1-2]-y1[j1])+(z1[i1-2]-z1[j1])*(z1[i1-2]-z1[j1]));  //PiNj
        d11=sqrt((x1[j1-2]-x1[i1])*(x1[j1-2]-x1[i1])+(y1[j1-2]-y1[i1])*(y1[j1-2]-y1[i1])+(z1[j1-2]-z1[i1])*(z1[j1-2]-z1[i1])); //NiPj
        if ((type[i1]=='G'&&type[j1]=='C')||(type[i1]=='C'&&type[j1]=='G'))  {hb=-3.325;}
        if ((type[i1]=='G'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='G'))  {hb=-2.10;}
        if ((type[i1]=='A'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='A'))  {hb=-2.10;}
//UHB=hb/(1+3.6*(d-8.94)*(d-8.94)+1.9*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.7*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));
        UHB=hb/(1+2.66*(d-8.94)*(d-8.94)+1.37*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.464*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));
    	return UHB;
}
void BasePairing()
{
    	int i,j;
    	float uN1=0.0,uuN1=0.0;
    	for (i=1;i<=N0;i++)  
    	{
            	a[i]=0;c[i]=0; 
            	for(j=i+3;j<=N0;j++) 
            	{
            		s[i][j]=0;ss[i][j]=0;
            	}
     	}
     	for (i=1;i<=(N0-12);i++) 
     	{   	    
         	if (fmod(i,3)==0)
         	{
             		for(j=i+12;j<=N0;j++)
             		{          
                 		if (fmod(j,3)==0)
                 		{     
                    			if ((type[i]=='G'&&type[j]=='C')||(type[i]=='C'&&type[j]=='G')
                    			||(type[i]=='A'&&type[j]=='U')||(type[i]=='U'&&type[j]=='A')
                    			||(type[i]=='G'&&type[j]=='U')||(type[i]=='U'&&type[j]=='G'))
                    			{
                        			uN1=0.0;
                        			uuN1=0.0;  
                        			if (c[i]==0&&c[j]==0) 
                        			{
                            				uN1=HB(i,j,x,y,z);
                            				if (uN1!=0) 
                            				{
                                				c[i]=1;c[j]=1; s[i][j]=1;                                               
                               
                             				}
                          			}
                          			uN=uN+uN1;    
                          			if (a[i]==0&&a[j]==0) 
                          			{
                              				uuN1=HB(i,j,xx,yy,zz); 
                              				if (uuN1!=0) 
                              				{
                                  				a[i]=1;
                                  				a[j]=1; 
                                  				ss[i][j]=1;
                              				}
                           			}
                           			uuN=uuN+uuN1; 
                      			} 
                      			else
                      			{
                          			uN1=0.0;
                          			uuN1=0.0;
                          			uN=uN+uN1; 
                          			uuN=uuN+uuN1;
                       			} 
                    		}        
                 	} 
     		}
	} 
}
//base stacking between two adjacent base pairs
// ~~~~~~~~~~~~~~~~Calculation of base stacking~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 float BSt(int i1,int j1,int k1,int k2)  //或许可以用四维数组书写，目前暂且如此吧！
{
   float B=0.0;
             // 5'-i1-k1-
             // 3'-j1/k2-
     if (((type[i1]=='A'&&type[k1]=='A')&&(type[j1]=='U'&&type[k2]=='U'))
     ||((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='A'&&type[k2]=='A'))) {B=-6.82-T*0.001*(-19.0-B0);}  // AA/UU
else if (((type[i1]=='A'&&type[k1]=='C')&&(type[j1]=='U'&&type[k2]=='G')) 
     ||((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='C'&&type[k2]=='A'))) {B=-11.4-T*0.001*(-29.5-B0);}  // AC/UG
else if (((type[i1]=='A'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='C'))
     ||((type[i1]=='C'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='A'))) {B=-10.48-T*0.001*(-27.1-B0);} // AG/UC
else if (((type[i1]=='A'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='U'))
     ||((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='A'))) {B=-3.21-T*0.001*(-8.6-B0);}   // AG/UU
else if ((type[i1]=='A'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='A'))   {B=-9.38-T*0.001*(-26.7-B0);}  // AU/UA
else if (((type[i1]=='A'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='G'))
     ||((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='A'))) {B=-8.81-T*0.001*(-24.0-B0);}  // AU/UG
else if (((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='A'&&type[k2]=='C')) 
     ||((type[i1]=='C'&&type[k1]=='A')&&(type[j1]=='G'&&type[k2]=='U'))) {B=-10.44-T*0.001*(-26.9-B0);} // UG/AC
else if (((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='C'&&type[k2]=='C'))
     ||((type[i1]=='C'&&type[k1]=='C')&&(type[j1]=='G'&&type[k2]=='G'))) {B=-13.39-T*0.001*(-32.7-B0);} // GG/CC
else if ((type[i1]=='C'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='C'))   {B=-10.64-T*0.001*(-26.7-B0);} // CG/GC
else if (((type[i1]=='C'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='U'))
     ||((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='C'))) {B=-5.61-T*0.001*(-13.5-B0);}  // CG/GU
else if (((type[i1]=='C'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='G'))
     ||((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='C'))) {B=-12.11-T*0.001*(-32.2-B0);} // CU/GG
else if (((type[i1]=='U'&&type[k1]=='C')&&(type[j1]=='A'&&type[k2]=='G'))
     ||((type[i1]=='G'&&type[k1]=='A')&&(type[j1]=='C'&&type[k2]=='U'))) {B=-12.44-T*0.001*(-32.5-B0);} // UC/AG
else if ((type[i1]=='G'&&type[k1]=='C')&&(type[j1]=='C'&&type[k2]=='G'))   {B=-14.88-T*0.001*(-36.9-B0);} // GC/CG
else if (((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='C'&&type[k2]=='U'))
     ||((type[i1]=='U'&&type[k1]=='C')&&(type[j1]=='G'&&type[k2]=='G'))) {B=-8.33-T*0.001*(-21.9-B0);}  // GG/CU
else if (((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='C'&&type[k2]=='G'))
     ||((type[i1]=='G'&&type[k1]=='C')&&(type[j1]=='U'&&type[k2]=='G'))) {B=-12.59-T*0.001*(-32.5-B0);} // GU/CG
else if ((type[i1]=='U'&&type[k1]=='A')&&(type[j1]=='A'&&type[k2]=='U'))   {B=-7.69-T*0.001*(-20.5-B0);}  // UA/AU
else if (((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='A'&&type[k2]=='U'))
     ||((type[i1]=='U'&&type[k1]=='A')&&(type[j1]=='G'&&type[k2]=='U'))) {B=-6.99-T*0.001*(-19.3-B0);}  // UG/AU
else if (((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='A'&&type[k2]=='G'))
     ||((type[i1]=='G'&&type[k1]=='A')&&(type[j1]=='U'&&type[k2]=='U'))) {B=-12.83-T*0.001*(-37.3-B0);} // UU/AG
else if (((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='G'))
     ||((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='U'))) {B=-13.47-T*0.001*(-44.9-B0);} // UU/GG
else if ((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='U'))   {B=-9.26-T*0.001*(-30.8-B0);}  // UG/GU
else if ((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='G'))   {
         if ((type[i1-3]=='G'&&type[k1+3]=='C')&&(type[j1+3]=='C'&&type[k2-3]=='G')){B=-16.66-T*0.001*(-50.3-B0);}
         else                                                                       {B=-14.59-T*0.001*(-51.2-B0);}}// GU/UG (CGUC/CUGG)
else    {B=0.0*(-6.82-T*0.001*(-19.0-B0));}   //Mismatched stacking
   return B;
}

float St(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
  	float d,d0,kst,ULJ,d1,d5,d6,d10,d12,d01,d05,d06,d010,d012/*,ulj1,ulj2*/;
  	d=sqrt((x1[i1]-x1[i1+3])*(x1[i1]-x1[i1+3])+(y1[i1]-y1[i1+3])*(y1[i1]-y1[i1+3])+(z1[i1]-z1[i1+3])*(z1[i1]-z1[i1+3])); 
  	d1=4.786/d; d5=d1*d1*d1*d1*d1;d6=d1*d5;d10=d5*d5;d12=d6*d6; 
  	d0=sqrt((x1[j1]-x1[j1-3])*(x1[j1]-x1[j1-3])+(y1[j1]-y1[j1-3])*(y1[j1]-y1[j1-3])+(z1[j1]-z1[j1-3])*(z1[j1]-z1[j1-3]));   
	d01=4.786/d0;d05=d01*d01*d01*d01*d01;d06=d01*d05;d010=d05*d05;d012=d06*d06;
  	kst=BSt(i1,j1,i1+3,j1-3);
  	if (kst>=0) 
  	{
         	ULJ=0.0;
  	}
  	else 
  	{
        	ULJ=-0.5*kst*((5*d12-6*d10)+(5*d012-6*d010));
  	}
  	return ULJ;
} 
void BaseStacking()
{
    	int i,j;
    	float us1,uus1,us2,uus2,us0,uus0;
    	for(i=1;i<=(N0-12);i++)
    	{  
     		us0=0.0;uus0=0.0;
      		if (fmod(i,3)==0)
      		{
          		for(j=i+12;j<=N0;j++)
          		{   
              			us1=0.0;
              			uus1=0.0; 
              			us2=0.0; 
              			uus2=0.0;
              			if (fmod(j,3)==0)
              			{
                   			if (s[i][j]==1&&s[i+3][j-3]==1) 
                   			{
                         			us1=St(i,j,x,y,z);
                   			}
                   			if (Energy==0&&ss[i][j]==1&&ss[i+3][j-3]==1) 
                   			{
                         			uus1=St(i,j,xx,yy,zz);
                   			}
                   			us0=us0+us1+us2; 
                   			uus0=uus0+uus1+uus2;  
                		}
            		}
            		us=us+us0; 
            		uus=uus+uus0;
        	}
     	}
}
// ~~~~~~~~~~~~Calculation of Exculded Volume between any two beads~~~~~~~~~~~~
float LJ0(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d,d1,d6,d12,r1,r2,ULJ;
     	r1=0.0;r2=0.0;
     	d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); 
     	if (fmod(i1,3)==0&&fmod(j1,3)==0) 
     	{
             	r1=2.15;
             	r2=2.15;
    	}
     	else 
     	{
             	r1=R[i1];
             	r2=R[j1];
     	}
     	if (d<=(r1+r2))
     	{
              	d1=(r1+r2)/(1.09*d); 
              	d6=d1*d1*d1*d1*d1*d1; 
              	d12=d6*d6;
              	ULJ=(4.0*e*(d12-d6)+e);
     	}
     	else ULJ=0.0;
     	return ULJ;
} 

void ExcludedVolume(int i0)
{
    	int i,j,i01,i01_end,i02,i02_end;    //i01-i01_end：atoms of move;
    	float ulj0,uulj0,ulj1,uulj1;
    	ulj0=0.0;
    	uulj0=0.0; 
    	ulj1=0;uulj1=0; 
    	if(Move==3) 
    	{                   //Moves of 3nt_fragment
      		i01=i0+1; i01_end=i0+8; 
      		i02=1; i02_end=N0;
      		for(i=i01;i<=i01_end;i++)   
      		{
           		ulj0=0.0;uulj0=0.0;
           		for (j=i02;j<=i02_end;j++)  
           		{        
               			if(j<=i0||j>i0+8) 
               			{  
                    			ulj1=0;uulj1=0; 
                    			ulj1=LJ0(i,j,x,y,z); 
                    			ulj0=ulj0+ulj1;                                
                    			uulj1=LJ0(i,j,xx,yy,zz); 
                    			uulj0=uulj0+uulj1; 
                		} 
            		} 
            		ulj=ulj+ulj0; 
            		uulj=uulj+uulj0;    
        	}
     	}
     	else if(Move==4)
     	{
        	for(i=i0+1;i<=i0+11;i++) 
        	{
          		ulj0=0.0;uulj0=0.0;
          		for (j=1;j<=N0;j++)	  
          		{
            			ulj1=0;uulj1=0;
            			if (j<=i0||j>i0+11)
            			{
              				ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                                
              				uulj1=LJ0(i,j,xx,yy,zz); uulj0=uulj0+uulj1;           
            			}
           		}
           		ulj=ulj+ulj0; uulj=uulj+uulj0;
          	}
      	}
     	else if(Move==5)
     	{
         	for(i=i0+1;i<=i0+14;i++) 
        	{
        		ulj0=0.0;uulj0=0.0;
        		for (j=1;j<=N0;j++)	  
        		{
         			ulj1=0;uulj1=0;
         			if (j<=i0||j>i0+14)
         			{
         				ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                                
         				uulj1=LJ0(i,j,xx,yy,zz); uulj0=uulj0+uulj1;           
          			}
         		}
        		ulj=ulj+ulj0; uulj=uulj+uulj0;
         	}
     	}
     	else if(Move==6)
     	{
         	for(i=i0+1;i<=i0+17;i++) 
        	{
        		ulj0=0.0;uulj0=0.0;
        		for (j=1;j<=N0;j++)	  
        		{
         			ulj1=0;uulj1=0;
         			if (j<=i0||j>i0+17)
         			{
         				ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                                
         				uulj1=LJ0(i,j,xx,yy,zz); uulj0=uulj0+uulj1;           
          			}
         		}
        		ulj=ulj+ulj0; uulj=uulj+uulj0;
         	}
     	}
     	else if(Move==7)
     	{
         	for(i=i0+1;i<=i0+20;i++) 
        	{
        		ulj0=0.0;uulj0=0.0;
        		for (j=1;j<=N0;j++)	  
        		{
         			ulj1=0;uulj1=0;
         			if (j<=i0||j>i0+20)
         			{
         				ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                                
         				uulj1=LJ0(i,j,xx,yy,zz); uulj0=uulj0+uulj1;           
          			}
         		}
        		ulj=ulj+ulj0; uulj=uulj+uulj0;
         	}
     	}
     	else 
     	{                             //Moves Pivot
       		if(fmod(i0,3)!=0) 
       		{             //backbone
          		if(i0<N1) 
          		{
               			i01=1; 
               			i01_end=i0; 
               			i02=i0+1; 
               			i02_end=N0;
          		}
          		else 
          		{
               			i01=i0; 
               			i01_end=N0; 
               			i02=1; 
               			i02_end=i0-1;
          		}  
          		for(i=i01;i<=i01_end;i++)   
          		{
              			ulj0=0.0;
              			uulj0=0.0;
              			for (j=i02;j<=i02_end;j++)  
              			{        
                  			ulj1=0;
                  			uulj1=0; 
                  			ulj1=LJ0(i,j,x,y,z); 
                  			ulj0=ulj0+ulj1;                                
                  			uulj1=LJ0(i,j,xx,yy,zz); 
                  			uulj0=uulj0+uulj1;  
               			} 
               			ulj=ulj+ulj0; 
               			uulj=uulj+uulj0;    
            		}
         	}
         	else 
         	{                          //side-chain
             		i01=i0; 
             		i01_end=i0; 
             		i02=1; 
             		i02_end=N0;
             		for(i=i01;i<=i01_end;i++)   
             		{
               			ulj0=0.0;
               			uulj0=0.0;
               			for (j=i02;j<=i02_end;j++)  
               			{        
                  			if(j!=i)  
                  			{  
                      				ulj1=0;
                      				uulj1=0; 
                      				ulj1=LJ0(i,j,x,y,z); 
                      				ulj0=ulj0+ulj1;                                
                      				uulj1=LJ0(i,j,xx,yy,zz); 
                      				uulj0=uulj0+uulj1;
                    			}  
                		} 
                		ulj=ulj+ulj0; 
                		uulj=uulj+uulj0;    
              		}
            	}      
        }
}
// ~~~~~~~Calculation of elsectrostatic between two bead: Debye combining with CC theroy~~~~~~~~~~~~
float QQ(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
   	float d,UQQ=0.0;
   	d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); 
   	UQQ=(1-frac0[(i1+2)/3])*(1-frac0[(j1+2)/3])*(330.9/(Ek*d))*exp(-d/Kd);
   	return UQQ;
}
 
void Electrostatic(int i0)
{
    	int i,j,i01,i01_end,i02,i02_end;
    	float uc0,uuc0,uc1,uuc1;
    	uc=0.0;uuc=0.0; uc1=0.0;uuc1=0.0;
    	if (salt==0) { uc=0.0;uuc=0.0;}
    	else
  	{
          	if (fmod(i0,3)==0)   { uc=0.0;uuc=0.0;}
          	else if(Move==3) 
          	{
             		i01=i0+1; i01_end=i0+8; 
             		i02=4; i02_end=N0-3;
             		for(i=i01;i<=i01_end;i++) 
             		{
               			if(fmod(i+2,3)==0)   
               			{ 
                    			uc0=0.0;
                    			uuc0=0.0; 
                    			for (j=i02;j<=i02_end;j++)  
                    			{        
                         			if((j<=i0||j>i0+8)&&fmod(j+2,3)==0) 
                         			{  
                               				uc1=0.0;
                               				uuc1=0.0; 
                               				uc1=QQ(i,j,x,y,z);    
                               				uc0=uc0+uc1;
                               				uuc1=QQ(i,j,xx,yy,zz); 
                               				uuc0=uuc0+uuc1; 
                           			} 
                     			} 
                     			uc=uc+uc0; 
                     			uuc=uuc+uuc0;    
                    		}    
                 	}
               	}
           	else 
                {
               		if(i0<N1) 
                	{
                        	i01=4; 
                          	i01_end=i0; 
                        	i02=i0+1; 
                         	i02_end=N0-3;
               		}
              		else 
               		{
                        	i01=i0; 
                          	i01_end=N0-3; 
                          	i02=4; 
                        	i02_end=i0-1;
                  	}
                 	for(i=i01;i<=i01_end;i++) 
             		{
                        	if(fmod(i+2,3)==0)   
                         	{ 
                               		uc0=0.0;
                                        uuc0=0.0; 
                                        for (j=i02;j<=i02_end;j++)  
                                        {        
                                             if(fmod(j+2,3)==0) 
                                             {  
                                                  uc1=0.0;
                                                  uuc1=0.0; 
                                                  uc1=QQ(i,j,x,y,z);    
                                                  uc0=uc0+uc1;
                                                  uuc1=QQ(i,j,xx,yy,zz); 
                                                  uuc0=uuc0+uuc1; 
                                              } 
                                         } 
                                         uc=uc+uc0; 
                                         uuc=uuc+uuc0;    
                            	}    
                   	}
     		}
   	}
}
// Bonded potential before and after conformational change: bond, angle and dihedral

/***********************************************************************************************************/

float dist(int i,int j,float x1[1000],float y1[1000],float z1[1000])
{
	return sqrt((x1[i]-x1[j])*(x1[i]-x1[j])+(y1[i]-y1[j])*(y1[i]-y1[j])+(z1[i]-z1[j])*(z1[i]-z1[j]));
}

void disbrute_qq()
{	
  	float dist(int i,int j,float x1[1000],float y1[1000],float z1[1000]);
  	void qqregeneration(float x1[1000],float y1[1000],float z1[1000],int zeta);	
  	int iii;
  	q4=5.998*1e-6*5.45*Ek*T*0.5*2.0;
  	qqregeneration(x,y,z,1);		//calculate for Na+
  	for(iii=1;iii<N+2;iii++)
  	{
   		frac1[iii]=frac[iii];
  	}
  	q4=5.998*1e-6*5.45*Ek*T*0.5*1.0;
  	qqregeneration(x,y,z,2);		//calculate for Mg++  
  	for(iii=1;iii<N+2;iii++)
  	{
   		frac2[iii]=frac[iii];
  	}
  	for(iii=1;iii<N+2;iii++)
  	{
  		frac0[iii]=fNa*frac1[iii]+(1-fNa)*frac2[iii];		//mix
  	}
}


void qqregeneration(float x1[1000],float y1[1000],float z1[1000],int zeta)
{
	
	float sigma1=0.;
	float sigma2=0.;
	float qianfrac[1000];
        int i,k,j,m;
	for(i=1;i<N+2;i++)
	{
		frac[i]=1-q4;
		fi[i]=0.;
	}
	
	for(m=1;m<15;m++)				//m<20
	{
		for(i=1;i<N+2;i++)			//calculate phi, Eq.8
		{
			for(j=1;j<N+2;j++)
			{
				if(j==i)
					continue;
			sigma2=sigma2+((frac[j]-1.)/dist((i-1)*3+1,(j-1)*3+1,x,y,z))*exp(-dist((i-1)*3+1,(j-1)*3+1,x,y,z)/Kd);
			}
			fi[i]=sigma2;
			
			sigma2=0.;
		}
		
		for(i=1;i<N+2;i++)			//calculate f, Eq.7
		{
			for(k=1;k<N+2;k++)
			{
				sigma1=sigma1+exp(-(1/(D))*zeta*fi[k]);
			}
      qianfrac[i]=frac[i];
			frac[i]=((N+1)*(1-q4)/sigma1)*exp(-(1/(D))*zeta*fi[i]);
      frac[i]=qianfrac[i]+ww*(frac[i]-qianfrac[i]);
			sigma1=0.;
		}

	}
}


/*******************************************************************************************************************/
void Bonded(int i0)
{
	if (Move==3)
  	{
     		if (fmod(i0+2,3)==0) 
     		{
        		ub=PC(i0+2,x,y,z)+CP(i0+8,x,y,z); 
        		uub=PC(i0+2,xx,yy,zz)+CP(i0+8,xx,yy,zz);
        		ue=CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+PCN(i0+2,x,y,z)+PCP(i0+8,x,y,z)+CPC(i0+8,x,y,z)+NCP(i0+8,x,y,z);
        		uue=CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+PCN(i0+2,xx,yy,zz)+PCP(i0+8,xx,yy,zz)+CPC(i0+8,xx,yy,zz)+NCP(i0+8,xx,yy,zz);
        		ud=PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+CPCN(i0+2,x,y,z)+CPCP(i0+5,x,y,z)+PCPC(i0+8,x,y,z)+CPCP(i0+8,x,y,z)+NCPC(i0+8,x,y,z)+CPCN(i0+11,x,y,z);
        		uud=PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+CPCN(i0+2,xx,yy,zz)+CPCP(i0+5,xx,yy,zz)+PCPC(i0+8,xx,yy,zz)+CPCP(i0+8,xx,yy,zz)+NCPC(i0+8,xx,yy,zz)+CPCN(i0+11,xx,yy,zz);
      		}
      		if (fmod(i0+1,3)==0)
      		{
         		ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z)+PC(i0+10,x,y,z); 
         		uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz)+PC(i0+10,xx,yy,zz);
         		ue=PCP(i0+1,x,y,z)+PCN(i0+1,x,y,z)+CPC(i0+1,x,y,z)+NCP(i0+1,x,y,z)+CPC(i0+7,x,y,z)+NCP(i0+7,x,y,z)+PCP(i0+10,x,y,z)+PCN(i0+10,x,y,z);
         		uue=PCP(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz)+CPC(i0+7,xx,yy,zz)+NCP(i0+7,xx,yy,zz)+PCP(i0+10,xx,yy,zz)+PCN(i0+10,xx,yy,zz);
         		ud=CPCP(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+1,x,y,z)+CPCN(i0+4,x,y,z)+CPCP(i0+7,x,y,z)+PCPC(i0+7,x,y,z)+PCPC(i0+10,x,y,z)+NCPC(i0+7,x,y,z)+CPCN(i0+10,x,y,z);
         		uud=CPCP(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+1,xx,yy,zz)+CPCN(i0+4,xx,yy,zz)+CPCP(i0+7,xx,yy,zz)+PCPC(i0+7,xx,yy,zz)+PCPC(i0+10,xx,yy,zz)+NCPC(i0+7,xx,yy,zz)+CPCN(i0+10,xx,yy,zz);
       		}
    	}
  	else if (Move==4)
  	{ 
     		if (fmod(i0+2,3)==0) 
     		{
        		ub=PC(i0+2,x,y,z)+CP(i0+11,x,y,z); 
        		uub=PC(i0+2,xx,yy,zz)+CP(i0+11,xx,yy,zz);
        		ue=CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+PCN(i0+2,x,y,z)+PCP(i0+11,x,y,z)+CPC(i0+11,x,y,z)+NCP(i0+11,x,y,z);
        		uue=CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+PCN(i0+2,xx,yy,zz)+PCP(i0+11,xx,yy,zz)+CPC(i0+11,xx,yy,zz)+NCP(i0+11,xx,yy,zz);
        		ud=PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+CPCN(i0+2,x,y,z)+CPCP(i0+8,x,y,z)+PCPC(i0+11,x,y,z)+CPCP(i0+11,x,y,z)+NCPC(i0+11,x,y,z)+CPCN(i0+14,x,y,z);
        		uud=PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+CPCN(i0+2,xx,yy,zz)+CPCP(i0+8,xx,yy,zz)+PCPC(i0+11,xx,yy,zz)+CPCP(i0+11,xx,yy,zz)+NCPC(i0+11,xx,yy,zz)+CPCN(i0+14,xx,yy,zz);
      		}
      		if (fmod(i0+1,3)==0)
      		{
         		ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z)+PC(i0+13,x,y,z); 
         		uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz)+PC(i0+13,xx,yy,zz);
         		ue=PCP(i0+1,x,y,z)+PCN(i0+1,x,y,z)+CPC(i0+1,x,y,z)+NCP(i0+1,x,y,z)+CPC(i0+10,x,y,z)+NCP(i0+10,x,y,z)+PCP(i0+13,x,y,z)+PCN(i0+13,x,y,z);
         		uue=PCP(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz)+CPC(i0+10,xx,yy,zz)+NCP(i0+10,xx,yy,zz)+PCP(i0+13,xx,yy,zz)+PCN(i0+13,xx,yy,zz);
         		ud=CPCP(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+1,x,y,z)+CPCN(i0+7,x,y,z)+CPCP(i0+10,x,y,z)+PCPC(i0+10,x,y,z)+PCPC(i0+13,x,y,z)+NCPC(i0+10,x,y,z)+CPCN(i0+13,x,y,z);
         		uud=CPCP(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+1,xx,yy,zz)+CPCN(i0+7,xx,yy,zz)+CPCP(i0+10,xx,yy,zz)+PCPC(i0+10,xx,yy,zz)+PCPC(i0+13,xx,yy,zz)+NCPC(i0+10,xx,yy,zz)+CPCN(i0+13,xx,yy,zz);
       		}
    	}
  	else if (Move==5)
  	{
     		if (fmod(i0+2,3)==0) 
     		{
        		ub=PC(i0+2,x,y,z)+CP(i0+14,x,y,z); 
        		uub=PC(i0+2,xx,yy,zz)+CP(i0+14,xx,yy,zz);
        		ue=CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+PCN(i0+2,x,y,z)+PCP(i0+14,x,y,z)+CPC(i0+14,x,y,z)+NCP(i0+14,x,y,z);
        		uue=CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+PCN(i0+2,xx,yy,zz)+PCP(i0+14,xx,yy,zz)+CPC(i0+14,xx,yy,zz)+NCP(i0+14,xx,yy,zz);
        		ud=PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+CPCN(i0+2,x,y,z)+CPCP(i0+11,x,y,z)+PCPC(i0+14,x,y,z)+CPCP(i0+14,x,y,z)+NCPC(i0+14,x,y,z)+CPCN(i0+17,x,y,z);
        		uud=PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+CPCN(i0+2,xx,yy,zz)+CPCP(i0+11,xx,yy,zz)+PCPC(i0+14,xx,yy,zz)+CPCP(i0+14,xx,yy,zz)+NCPC(i0+14,xx,yy,zz)+CPCN(i0+17,xx,yy,zz);
      		}
      		if (fmod(i0+1,3)==0)
      		{
         		ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z)+PC(i0+16,x,y,z); 
         		uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz)+PC(i0+16,xx,yy,zz);
         		ue=PCP(i0+1,x,y,z)+PCN(i0+1,x,y,z)+CPC(i0+1,x,y,z)+NCP(i0+1,x,y,z)+CPC(i0+13,x,y,z)+NCP(i0+13,x,y,z)+PCP(i0+16,x,y,z)+PCN(i0+16,x,y,z);
         		uue=PCP(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz)+CPC(i0+13,xx,yy,zz)+NCP(i0+13,xx,yy,zz)+PCP(i0+16,xx,yy,zz)+PCN(i0+16,xx,yy,zz);
         		ud=CPCP(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+1,x,y,z)+CPCN(i0+10,x,y,z)+CPCP(i0+13,x,y,z)+PCPC(i0+13,x,y,z)+PCPC(i0+16,x,y,z)+NCPC(i0+13,x,y,z)+CPCN(i0+16,x,y,z);
         		uud=CPCP(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+1,xx,yy,zz)+CPCN(i0+10,xx,yy,zz)+CPCP(i0+13,xx,yy,zz)+PCPC(i0+13,xx,yy,zz)+PCPC(i0+16,xx,yy,zz)+NCPC(i0+13,xx,yy,zz)+CPCN(i0+16,xx,yy,zz);
       		}
    	}
  	else if(Move==6)
	{
     		if (fmod(i0+2,3)==0) 
     		{
        		ub=PC(i0+2,x,y,z)+CP(i0+17,x,y,z); 
        		uub=PC(i0+2,xx,yy,zz)+CP(i0+17,xx,yy,zz);
        		ue=CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+PCN(i0+2,x,y,z)+PCP(i0+17,x,y,z)+CPC(i0+17,x,y,z)+NCP(i0+17,x,y,z);
        		uue=CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+PCN(i0+2,xx,yy,zz)+PCP(i0+17,xx,yy,zz)+CPC(i0+17,xx,yy,zz)+NCP(i0+17,xx,yy,zz);
        		ud=PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+CPCN(i0+2,x,y,z)+CPCP(i0+14,x,y,z)+PCPC(i0+17,x,y,z)+CPCP(i0+17,x,y,z)+NCPC(i0+17,x,y,z)+CPCN(i0+20,x,y,z);
        		uud=PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+CPCN(i0+2,xx,yy,zz)+CPCP(i0+14,xx,yy,zz)+PCPC(i0+17,xx,yy,zz)+CPCP(i0+17,xx,yy,zz)+NCPC(i0+17,xx,yy,zz)+CPCN(i0+20,xx,yy,zz);
      		}
      		if (fmod(i0+1,3)==0)
      		{
         		ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z)+PC(i0+19,x,y,z); 
         		uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz)+PC(i0+19,xx,yy,zz);
         		ue=PCP(i0+1,x,y,z)+PCN(i0+1,x,y,z)+CPC(i0+1,x,y,z)+NCP(i0+1,x,y,z)+CPC(i0+16,x,y,z)+NCP(i0+16,x,y,z)+PCP(i0+19,x,y,z)+PCN(i0+19,x,y,z);
         		uue=PCP(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz)+CPC(i0+16,xx,yy,zz)+NCP(i0+16,xx,yy,zz)+PCP(i0+19,xx,yy,zz)+PCN(i0+19,xx,yy,zz);
         		ud=CPCP(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+1,x,y,z)+CPCN(i0+13,x,y,z)+CPCP(i0+16,x,y,z)+PCPC(i0+16,x,y,z)+PCPC(i0+19,x,y,z)+NCPC(i0+16,x,y,z)+CPCN(i0+19,x,y,z);
         		uud=CPCP(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+1,xx,yy,zz)+CPCN(i0+13,xx,yy,zz)+CPCP(i0+16,xx,yy,zz)+PCPC(i0+16,xx,yy,zz)+PCPC(i0+19,xx,yy,zz)+NCPC(i0+16,xx,yy,zz)+CPCN(i0+19,xx,yy,zz);
       		}
    	}
    	else if(Move==7)
  	{
     		if (fmod(i0+2,3)==0) 
     		{
        		ub=PC(i0+2,x,y,z)+CP(i0+20,x,y,z); 
        		uub=PC(i0+2,xx,yy,zz)+CP(i0+20,xx,yy,zz);
        		ue=CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+PCN(i0+2,x,y,z)+PCP(i0+20,x,y,z)+CPC(i0+20,x,y,z)+NCP(i0+20,x,y,z);
        		uue=CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+PCN(i0+2,xx,yy,zz)+PCP(i0+20,xx,yy,zz)+CPC(i0+20,xx,yy,zz)+NCP(i0+20,xx,yy,zz);
        		ud=PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+CPCN(i0+2,x,y,z)+CPCP(i0+17,x,y,z)+PCPC(i0+20,x,y,z)+CPCP(i0+20,x,y,z)+NCPC(i0+20,x,y,z)+CPCN(i0+23,x,y,z);
        		uud=PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+CPCN(i0+2,xx,yy,zz)+CPCP(i0+17,xx,yy,zz)+PCPC(i0+20,xx,yy,zz)+CPCP(i0+20,xx,yy,zz)+NCPC(i0+20,xx,yy,zz)+CPCN(i0+23,xx,yy,zz);
      		}
      		if (fmod(i0+1,3)==0)
      		{
         		ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z)+PC(i0+22,x,y,z); 
         		uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz)+PC(i0+22,xx,yy,zz);
         		ue=PCP(i0+1,x,y,z)+PCN(i0+1,x,y,z)+CPC(i0+1,x,y,z)+NCP(i0+1,x,y,z)+CPC(i0+19,x,y,z)+NCP(i0+19,x,y,z)+PCP(i0+22,x,y,z)+PCN(i0+22,x,y,z);
         		uue=PCP(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz)+CPC(i0+19,xx,yy,zz)+NCP(i0+19,xx,yy,zz)+PCP(i0+22,xx,yy,zz)+PCN(i0+22,xx,yy,zz);
         		ud=CPCP(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+1,x,y,z)+CPCN(i0+16,x,y,z)+CPCP(i0+19,x,y,z)+PCPC(i0+19,x,y,z)+PCPC(i0+22,x,y,z)+NCPC(i0+19,x,y,z)+CPCN(i0+22,x,y,z);
         		uud=CPCP(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+1,xx,yy,zz)+CPCN(i0+16,xx,yy,zz)+CPCP(i0+19,xx,yy,zz)+PCPC(i0+19,xx,yy,zz)+PCPC(i0+22,xx,yy,zz)+NCPC(i0+19,xx,yy,zz)+CPCN(i0+22,xx,yy,zz);
       		}
    	}
     	else
     	{
         	if (fmod((i0+1),3)==0)       // Choose C atoms
         	{
               		if (i0<N1)
               		{
                      		ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z);  
                      		uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz);
                      		ue=NCP(i0+1,x,y,z)+PCP(i0+1,x,y,z)+CPC(i0+1,x,y,z)+PCN(i0+1,x,y,z);
                      		uue=NCP(i0+1,xx,yy,zz)+PCP(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz);
                	}
                	else
                	{
                       		ub=PC(i0+1,x,y,z);           
                       		uub=PC(i0+1,xx,yy,zz);
                       		ue=PCN(i0+1,x,y,z)+PCP(i0+1,x,y,z)+CPC(i0-2,x,y,z);                               
                       		uue=PCN(i0+1,xx,yy,zz)+PCP(i0+1,xx,yy,zz)+CPC(i0-2,xx,yy,zz);    
                	}
                	if (i0==2) 
                	{
       		       		ud=PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+4,x,y,z);
                       		uud=PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+4,xx,yy,zz);
                	}	 
                	else if (i0==(N0-2)) 
                	{
                       		ud=PCPC(i0-2,x,y,z)+CPCP(i0-2,x,y,z)+NCPC(i0-2,x,y,z)+CPCN(i0+1,x,y,z);
                       		uud=PCPC(i0-2,xx,yy,zz)+CPCP(i0-2,xx,yy,zz)+NCPC(i0-2,xx,yy,zz)+CPCN(i0+1,xx,yy,zz);
                 	}   
                 	else 
                 	{
                       		if (i0<N1)
                       		{
                           		ud=CPCP(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+4,x,y,z)+CPCN(i0+1,x,y,z);
                           		uud=CPCP(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+4,xx,yy,zz)+CPCN(i0+1,xx,yy,zz);
                       		}
                       		else
                       		{
                           		ud=PCPC(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0-2,x,y,z)+NCPC(i0-2,x,y,z)+CPCN(i0+1,x,y,z);
                           		uud=PCPC(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0-2,xx,yy,zz)+NCPC(i0-2,xx,yy,zz)+CPCN(i0+1,xx,yy,zz);
                       		}
                	 }   
           	}
           	if (fmod((i0+2),3)==0)       // Choose P atoms
           	{
                	if (i0<N1) 
                  	{
                		ub=PC(i0+2,x,y,z);
                     		uub=PC(i0+2,xx,yy,zz);
                  	}
                	else                 
                   	{ 
                    		ub=CP(i0-1,x,y,z);
                         	uub=CP(i0-1,xx,yy,zz);
                 	}
              		if (i0==1) 
                	{
                 		ue=PCN(i0+2,x,y,z)+PCP(i0+2,x,y,z);    
                         	uue=PCN(i0+2,xx,yy,zz)+PCP(i0+2,xx,yy,zz); 
                          	ud=PCPC(i0+2,x,y,z);                   
                          	uud=PCPC(i0+2,xx,yy,zz);
                 	}
                	else if (i0==N0) 
                   	{
                    		ue=NCP(i0-1,x,y,z)+PCP(i0-1,x,y,z); 
                            	uue=NCP(i0-1,xx,yy,zz)+PCP(i0-1,xx,yy,zz);
                           	ud=CPCP(i0-4,x,y,z);                
                         	uud=CPCP(i0-4,xx,yy,zz);
               		} 
              		else  
                	{
                     		if (i0<N1)
                       		{
                              		ue=PCN(i0+2,x,y,z)+PCP(i0+2,x,y,z)+CPC(i0-1,x,y,z);              
                            		uue=PCN(i0+2,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+CPC(i0-1,xx,yy,zz); 
                             		ud=CPCP(i0-1,x,y,z)+PCPC(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+NCPC(i0-1,x,y,z)+CPCN(i0+2,x,y,z);
                            		uud=CPCP(i0-1,xx,yy,zz)+PCPC(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+CPCN(i0+2,xx,yy,zz);
                       		}
                         	else 
                          	{
                           		ue=NCP(i0-1,x,y,z)+PCP(i0-1,x,y,z)+CPC(i0-1,x,y,z);              
                               		uue=NCP(i0-1,xx,yy,zz)+PCP(i0-1,xx,yy,zz)+CPC(i0-1,xx,yy,zz);
                             		ud=CPCP(i0-4,x,y,z)+PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+CPCN(i0+2,x,y,z);
                                	uud=CPCP(i0-4,xx,yy,zz)+PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+CPCN(i0+2,xx,yy,zz);
                        	}
              		}    
    		}
                if (fmod(i0,3)==0)             // Choose N atoms
                {
      			ub=CN(i0,x,y,z);                 
               		uub=CN(i0,xx,yy,zz);
               		ue=NCP(i0,x,y,z)+PCN(i0,x,y,z); 
             		uue=NCP(i0,xx,yy,zz)+PCN(i0,xx,yy,zz); 
              		if (i0==3) 
             		{
              			ud=NCPC(i0,x,y,z); 
                		uud=NCPC(i0,xx,yy,zz);
              		}
      			else if (i0==N0-1) 
         		{
                        	ud=CPCN(i0,x,y,z); 
                       		uud=CPCN(i0,xx,yy,zz);
                       	}
                     	else       
               		{
                   		ud=NCPC(i0,x,y,z)+CPCN(i0,x,y,z); 
                        	uud=NCPC(i0,xx,yy,zz)+CPCN(i0,xx,yy,zz);
                  	}
        	}
 	}
}



float UCoS(int i1,int j1,int k1,int k2,float x[1000],float y[1000],float z[1000])
{
    	float dco1,Ucos,kco,dco2,Ucos1=0.0,Ucos2=0.0,dco_1,dco_2,kdco_1,kdco_2;
    	kco=BSt(i1,j1,k1,k2);
    	dco1=sqrt((x[i1]-x[k1])*(x[i1]-x[k1])+(y[i1]-y[k1])*(y[i1]-y[k1])+(z[i1]-z[k1])*(z[i1]-z[k1])); 
    	dco2=sqrt((x[k2]-x[j1])*(x[k2]-x[j1])+(y[k2]-y[j1])*(y[k2]-y[j1])+(z[k2]-z[j1])*(z[k2]-z[j1])); 
    	if (Bulge==1&&j1-k2>=12) 
    	{
             	dco_1=5.0; dco_2=2*5.0; 
             	kdco_1=2.5; kdco_2=2*2.5; 
    	}
    	else 
    	{
             	dco_1=5.0; dco_2=5.0; 
             	kdco_1=2.5; kdco_2=2.5; 
    	}   // 针对bulge做修正，如果间隔3nt，则最优距离dco=2dco，kdco=2kdco；
     // if (fmod(t,100000)==0) {printf("%f %f %f %f %f %f\n",dco1,dco_1,kdco_1,dco2,dco_2,kdco_2);}
   	Ucos1=-0.5*(kco-1.0)*((1-exp(-(dco1-dco_1)/kdco_1))*(1-exp(-(dco1-dco_1)/kdco_1))-1);
   	Ucos2=-0.5*(kco-1.0)*((1-exp(-(dco2-dco_2)/kdco_2))*(1-exp(-(dco2-dco_2)/kdco_2))-1);
   	Ucos=Ucos1+Ucos2; 
// if (fmod(t,tprint)==0&&Energy==1) printf("%d %d %d %d %f %f %f %f %f %f\n",i1/3,j1/3,k1/3,k2/3,kco,dco1,dco2,Ucos1,Ucos2,Ucos);
   	return Ucos;
}

void CoaxialStacking()
{
   	int i,j,k1,k2;
   	float uco0,uuco0,uco1,uuco1;
   	for(i=9;i<=N0-16;i++)
   	{
      		if (fmod(i,3)==0)
      		{
          		for(j=i+12;j<=N0-4;j++)
          		{
              			if (fmod(j,3)==0&&((s[i][j]==1&&s[i-3][j+3]==1&&s[i+3][j-3]==0)||(ss[i][j]==1&&ss[i-3][j+3]==1&&ss[i+3][j-3]==0)))
              			{
                    			uco0=0.0;
                    			uuco0=0.0;
                    			for(k1=i+3;k1<j;k1++)
                    			{      
                        			if (fmod(k1,3)==0)
                        			{ 
                            				for(k2=k1+12;k2<=N0-4;k2++)
                            				{
                               					if (fmod(k2,3)==0)
                               					{
                                    					if (j>k2&&(abs(i-k1)!=3||abs(k2-j)!=3)&&((s[k1][k2]==1&&s[k1+3][k2-3]==1&&s[k1-3][k2+3]==0)||(ss[k1][k2]==1&&ss[k1+3][k2-3]==1&&ss[k1-3][k2+3]==0)))
                                    					{
                                         					uco1=0.0; 
                                         					uuco1=0.0;
                                         					if (k1==i+3||j==k2+3)          //bulge loop
                                         					{
                                             						Bulge=1;
                                             						if (k1==i+3) 
                                             						{
                                                 						if (s[i][j]==1&&s[k1][k2]==1&&s[i+3][j-3]==0&&s[k1-3][k2+3]==0) 
                                                 						{
                                                    							uco1=UCoS(i,j,k1,k2,x,y,z);
                                                 						} 
                                                 						if (Energy==0&&ss[i][j]==1&&ss[k1][k2]==1&&ss[i+3][j-3]==0&&ss[k1-3][k2+3]==0) 
                                                 						{
                                                    							uuco1=UCoS(i,j,k1,k2,xx,yy,zz);
                                                 						}
                                              						}
                                              						if (j==k2+3) 
                                              						{
                                                   						if (s[i][j]==1&&s[k1][k2]==1&&s[i+3][j-3]==0&&s[k1-3][k2+3]==0) 
                                                   						{
                                                        						uco1=UCoS(k2,k1,j,i,x,y,z);
                                                   						} 
                                                   						if (Energy==0&&ss[i][j]==1&&ss[k1][k2]==1&&ss[i+3][j-3]==0&&ss[k1-3][k2+3]==0) 
                                                   						{
                                                         						uuco1=UCoS(k2,k1,j,i,xx,yy,zz);
                                                   						}
                                               						}
                                               						uco0=uco0+uco1; 
                                               						uuco0=uuco0+uuco1; 
                                             					}
                                             					else {Bulge=0;}
                                         				}
                                         				if (k2>j+9&&j>k1&&k1>i+9&&j<=k1+6&&((s[k1][k2]==1&&s[k1-3][k2+3]==1&&s[k1+3][k2-3]==0)||(ss[k1][k2]==1&&ss[k1-3][k2+3]==1&&ss[k1+3][k2-3]==0)))   //L_loop1>=1 && L_loop0.5>=0 and <=2      Pseudoknot
                                         				{
                                              					uco1=0.0; 
                                              					uuco1=0.0; 
                                              					Pseudoknot=1;
                                              					if (s[i][j]==1&&s[k1][k2]==1&&s[i+3][j-3]==0&&s[k1+3][k2-3]==0) 
                                              					{
                                                   					uco1=UCoS(k1,k2,j,i,x,y,z);
                                              					}
                                              					if (Energy==0&&ss[i][j]==1&&ss[k1][k2]==1&&ss[i+3][j-3]==0&&ss[k1+3][k2-3]==0) 
                                              					{
                                                   					uuco1=UCoS(k1,k2,j,i,xx,yy,zz);
                                              					}
                                              					uco0=uco0+uco1; 
                                              					uuco0=uuco0+uuco1;        
                                          				} 
                                          				else Pseudoknot=0;
                                  				}  
                            				}
                     				}
               				}
               				uco=uco+uco0; 
               				uuco=uuco+uuco0;
       				}
      			}
   		}
	}  
}


/*float uys(int i1,int k1,float x1[1000],float y1[1000],float z1[1000])
{
       	float d,ud;
       	d=sqrt((x1[i1]-x1[k1])*(x1[i1]-x1[k1])+(y1[i1]-y1[k1])*(y1[i1]-y1[k1])+(z1[i1]-z1[k1])*(z1[i1]-z1[k1])); 
       	if(d>9.39||d<8.55)  	{ud=100*(d-9)*(d-9);}
       	else              	{ud=2*HB1(i1,k1,x1,y1,z1);}
       	return ud;
}*/
float uys(int i1,int k1,float x1[1000],float y1[1000],float z1[1000])
{
       	float d,d0,d01,d1,d11,ux,u0,u01,u1,u11,ud;
       	float kk;
       	if(t<100000)	{kk=0.1;}
       	else		{kk=5.0;}
       	d=sqrt((x1[i1]-x1[k1])*(x1[i1]-x1[k1])+(y1[i1]-y1[k1])*(y1[i1]-y1[k1])+(z1[i1]-z1[k1])*(z1[i1]-z1[k1])); 
 	d0=sqrt((x1[i1]-x1[k1-1])*(x1[i1]-x1[k1-1])+(y1[i1]-y1[k1-1])*(y1[i1]-y1[k1-1])+(z1[i1]-z1[k1-1])*(z1[i1]-z1[k1-1]));  //NiCj
        d01=sqrt((x1[k1]-x1[i1-1])*(x1[k1]-x1[i1-1])+(y1[k1]-y1[i1-1])*(y1[k1]-y1[i1-1])+(z1[k1]-z1[i1-1])*(z1[k1]-z1[i1-1])); //CiNj
        d1=sqrt((x1[i1-2]-x1[k1])*(x1[i1-2]-x1[k1])+(y1[i1-2]-y1[k1])*(y1[i1-2]-y1[k1])+(z1[i1-2]-z1[k1])*(z1[i1-2]-z1[k1]));  //PiNj
        d11=sqrt((x1[k1-2]-x1[i1])*(x1[k1-2]-x1[i1])+(y1[k1-2]-y1[i1])*(y1[k1-2]-y1[i1])+(z1[k1-2]-z1[i1])*(z1[k1-2]-z1[i1])); //NiPj
     	ux=kk*(d-8.94)*(d-8.94);
	u0=kk*(d0-12.15)*(d0-12.15);
	u01=kk*(d01-12.15)*(d01-12.15);
	u1=kk*(d1-13.92)*(d1-13.92);
	u11=kk*(d11-13.92)*(d11-13.92);
	ud=ux+u0+u01+u1+u11;
       	return ud;
}
float kissing_34(float x1[1000],float y1[1000],float z1[1000])
{
	float sx1=0.0,sy1=0.0,sz1=0.0,sx2=0.0,sy2=0.0,sz2=0.0;
	float d,uf=0.0;
	int i;
	for(i=ssi1;i<=ssj1;i++)	{sx1=sx1+x1[i];sy1=sy1+y1[i];sz1=sz1+z1[i];}
	sx1=sx1/(ssj1-ssi1+1);sy1=sy1/(ssj1-ssi1+1);sz1=sz1/(ssj1-ssi1+1);
	for(i=ssi2;i<=ssj2;i++)	{sx2=sx2+x1[i];sy2=sy2+y1[i];sz2=sz2+z1[i];}
	sx2=sx2/(ssj2-ssi2+1);sy2=sy2/(ssj2-ssi2+1);sz2=sz2/(ssj2-ssi2+1);
	d=sqrt((sx1-sx2)*(sx1-sx2)+(sy1-sy2)*(sy1-sy2)+(sz1-sz2)*(sz1-sz2));
	//if(fmod(t,100)==0&&cycle==1)	{printf("%f\n",d);}
	uf=5*(d-15)*(d-15);
	return uf;
}

void constraint()
{
       	float u1=0.0,uu1=0.0;
       	float kissing_34(),uk=0.0,uuk=0.0;
        ucst=0.0;uucst=0.0;
        int i,j;  
      	for(i=3;i<=N0-12;i=i+3)
    	{ 
             	for(j=i+12;j<=N0-1;j=j+3)
             	{
             		if(t<400000)
             		{
                    		if(ccs[i][j]==1&&tertiary[i]==0)
                    		{
                    			u1=uys(i,j,x,y,z);  	
                    			uu1=uys(i,j,xx,yy,zz);
                    		}
                    		else
                    		{	 
                         		u1=0.0;uu1=0.0;
                    		}
                    	}
                    	else
                    	{
                    		if(ccs[i][j]==1)
                    		{
                    			u1=uys(i,j,x,y,z);  	
                    			uu1=uys(i,j,xx,yy,zz);
                         	
                    		}
                    		else
                    		{	 
                         		u1=0.0;uu1=0.0;
                    		}
                    	}
                    	ucst=ucst+u1;
                    	uucst=uucst+uu1;
             	}         
     	}
     	if(junction[tdn]==1)		{uk=0.0;uuk=0.0;}
     	if(junction[tdn]==3)	
     	{
     		if(condition[tdn]!=1)	{uk=0.0;uuk=0.0;}
     		else			{uk=kissing_34(x,y,z);uuk=kissing_34(xx,yy,zz);}
     	}
     	ucst=ucst+uk;uucst=uucst+uuk;
       	if((fmod(t,10000)==0)&&cycle[tdn]==1) 
       	{
       		printf("======================== thread %d  ==========================\n",tdn);
       		printf("Predit 3D structure based secondary structure ......\n");
       		printf("\n");
       	}
}
// Energy of one confromation
void ENERGY()
{
   	int i,j;
   	float u0=0.0,uc1=0.0,ulj1=0.0,ulj0=0.0,uc0=0.0,uN1=0.0;
   	Energy=1;
   	if (fmod(t,tenergy)==0||t==1)
   	{
      		ulj=0.0;ub=0.0;
      		uN=0.0;us=0.0;
      		U=0.0; uco=0.0;
      		for(i=1;i<=N0;i++) 
      		{
            		c[i]=0;
            		for(j=1;j<=N0;j++) 
            		{
                  		s[i][j]=0;
            		}
      		}
      		for(i=1;i<=N0-1;i++) 
      		{
            		ulj0=0.0;
            		for (j=i+1;j<=N0;j++)	  
            		{
                  		ulj1=0.0;
                  		uc1=0.0;
                  		ulj1=LJ0(i,j,x,y,z); 
                  		ulj0=ulj0+ulj1;                                          
            		}
            		ulj=ulj+ulj0; 
       		}
/********************************/
       		if (salt==0) 
       		{
            		uc=0.0;uuc=0.0;
        	}
        	else
        	{
           		for(i=4;i<=N0-4;i++) 
           		{
               			uc0=0.0;
               			if (fmod(i+2,3)==0)
              			{
                    			for (j=i+1;j<=N0-3;j++)	  
                    			{
                        			uc1=0.0;
                        			if (fmod(j+2,3)==0)
                        			{
                              				uc1=QQ(i,j,x,y,z); 
                              				uc0=uc0+uc1;
                         			}
                     			}
                     			uc=uc+uc0; 
                 		}
              		}
          	}
/*******************************/
         	for (i=1;i<=N0;i++)  
         	{
              		c[i]=0; 
              		for(j=1;j<=N0;j++) 
              		{
                              	s[i][j]=0;
               		}
          	}
          	for(i=1;i<=N0;i++)  
          	{
                	if (fmod(i,3)==0) 
                	{
                      		for(j=i+12;j<=N0;j++)  
                      		{
                            		if (fmod(j,3)==0) 
                            		{ 
                                 		if ((type[i]=='G'&&type[j]=='C')||(type[i]=='C'&&type[j]=='G')
                                 		||(type[i]=='A'&&type[j]=='U')||(type[i]=='U'&&type[j]=='A')
                                 		||(type[i]=='G'&&type[j]=='U')||(type[i]=='U'&&type[j]=='G'))
                                 		{
                                      			uN1=0.0;
                                      			if (c[i]==0&&c[j]==0) 
                                      			{
                                             			uN1=HB(i,j,x,y,z);
                                             			if (uN1!=0) 
                                             			{
                                                    			c[i]=1;
                                                    			c[j]=1; 
                                                    			s[i][j]=1; 
                                                    			if(fmod(t,tbp)==0)
                                                    			{
                                                           			fprintf(fp_sec_stru[tdn],"%d %c %d %c\n",i/3,type[i],j/3,type[j]);fflush(fp_sec_stru[tdn]);
                                                     			}
                                               			}
                                         		}
                                         		uN=uN+uN1;     
                                     		} 
                                 	}
                            	}
                     	} 
                }
                for (i=1;i<=N0;i++)
                {
                     	u0=0.0;
                     	if (fmod(i,3)==0)
                     	{
                           	if (i==3) 
                            	{
                                    	u0=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z)+PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z)
                                    	+PCPC(i,x,y,z)+CPCP(i,x,y,z)+NCPC(i,x,y,z);
                            	}
                            	else if (i==(N0-1))
                            	{
                                     	u0=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z)+PCP(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z)+CPCN(i,x,y,z);
                             	} 
                             	else 
                             	{
                                      	u0=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z)+PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z)
                                      	+PCPC(i,x,y,z)+CPCP(i,x,y,z)+CPCN(i,x,y,z)+NCPC(i,x,y,z);
                              	}
                              	ub=ub+u0*0.5963;
                       	}
            	} 
        	if(fmod(t,tbp)==0) 
          	{
                     	fprintf(fp_sec_stru[tdn],"\n");fflush(fp_sec_stru[tdn]);
             	}
              	BaseStacking(); 
          	CoaxialStacking(); //Triplex(); //TerminalMismatch(); DanglingEnd();
           	constraint();
          	U=ulj+ub+uN+us+uc+uco/*+uTri+umis+udang*/; l++;
             	fprintf(fp_energy[tdn],"%d %d %f %f %f %f %f %f %f %f\n",t,l,U,ulj,ub,uN,us,uc,uco,ucst);fflush(fp_energy[tdn]);
   	}
}

void FOLD(int i0)
{
    	int i;
    	float prob;
    	prob=rand()/(RAND_MAX+1.);
     	if ((prob<0.25)&&i0>12&&i0<N0-12&&fmod(i0,3)!=0)    //3nt fragment translation;
      	{
        	Move=3; 
        	for(i=1;i<=N0;i++)  
        	{      //i0+1-->i0+8:translated
          		if (i>i0&&i-i0<=8) 
          		{
              			Translate(i);xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
          		}
          		else  
          		{
             			xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
           		}  
        	}
      	}
      	else if ((prob<0.5&&prob>=0.25)&&i0>15&&i0<N0-15&&fmod(i0,3)!=0)    //4nt fragment translation;
      	{
        	Move=4;
        	for(i=1;i<=N0;i++)  
        	{      //i0+1-->i0+8:translated
          		if (i>i0&&i-i0<=11) 
          		{
              			Translate(i);xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
          		}
          		else  
          		{
             			xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
           		}  
        	}
      	}
      	else if ((prob<0.75&&prob>=0.5)&&i0>18&&i0<N0-18&&fmod(i0,3)!=0)    //5nt fragment translation;
      	{
        	Move=5;
        	for(i=1;i<=N0;i++)  
        	{      //i0+1-->i0+8:translated
          		if (i>i0&&i-i0<=14) 
          		{
              			Translate(i);xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
          		}
          		else  
          		{
             			xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
           		}  
        	}
      	}
      	else if (i0>21&&i0<N0-21&&fmod(i0,3)!=0)    //6nt fragment translation;
      	{
        	Move=6;
        	for(i=1;i<=N0;i++)  
        	{      //i0+1-->i0+8:translated
          		if (i>i0&&i-i0<=17) 
          		{
              			Translate(i);xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
          		}
          		else  
          		{
             			xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
           		}  
        	}
      	}
       	else
       	{ 
            	Move=1; //pivot 
            	if (i0<N1)                    //Move 5' end chain
            	{  
                	if (fmod(i0,3)!=0)            //P,C4'is chosen
                	{
                      		Rand01();
                      		for (i=1;i<=N0;i++)
                      		{
                            		if (i>i0)  
                            		{
                                  		xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
                                  		xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
                            		}
                            		else       
                            		{
                                  		Translate(i);  
                                  //printf("%f %f %f\n",x[15],y[15],z[15]);
                                  		if (i==i0) 
                                  		{
                                      			xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
                                   		} 
                            		}
                        	}  
                        	if (fmod(t,tran)==0) 
                        	{
                            		for(i=1;i<=i0-1;i++) 
                            		{
                                 		xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
                             		}
                        	} //non pivot
                        	else                 
                        	{
                             		Rand01();   
                             		if (i0!=1) 
                             		{
                                    		for(i=1;i<=i0-1;i++) 
                                    		{
                                      			Pivot(i,i0);
                                    		}
                              		}    
                         	} //pivot
                    	}
                    	else                      //N is chosen
                    	{
                         	for(i=1;i<=N0;i++)
                         	{
                              		if (i!=i0) 
                              		{
                                      		xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
                                      		xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
                              		}
                              		else  
                              		{
                                     		MoveN(i0);
                              		}
                          	}            
                      	}
                }
             	if (i0>N1)                   //Move 3' end chain
             	{  
                      	if (fmod(i0,3)!=0)
                      	{
                         	Rand01();
                               	for (i=1;i<=N0;i++)
                               	{
                               		if (i<i0) 
                                      	{
                                       		xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
                                             	xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
                                       	}
                                       	else      
                                       	{
                                             	Translate(i); 
                                             	if (i==i0) 
                                             	{
                                                    	xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
                                              	} 
                                        }
                              	}  
                             	if (fmod(t,tran)==0) 
                         	{
                               		for(i=i0+1;i<=N0;i++) 
                                        {
                                   		xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
                                        } 
                          	}
                       		else                 
                        	{
                                     	Rand01();  
                                	if (i0!=N0)   
                                      	{
                                             	for (i=i0+1;i<=N0;i++) 
                                             	{ 
                                                  	Pivot(i,i0);
                                             	}
                                     	} 
                     		}
               		}
             		else
                 	{
                     		for (i=1;i<=N0;i++)
                          	{
                              		if (i!=i0) 
                                     	{
                                 		xx[i]=x[i];yy[i]=y[i];zz[i]=z[i]; 
                                               	xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
                                  	}
                                	else 
                                 	{ 
                                                MoveN(i0);
                                     	}
                       		}            
          		}
      		}
 	}
  	ExcludedVolume(i0);   
        //    Electrostatic(i0);   
    	Bonded(i0); //Excluded volume  &&  electrostatic  && bonded potential  
          //  BasePairing(); 
          //  BaseStacking(); 
	constraint();
        //    CoaxialStacking(); //Triplex(); //TerminalMismatch(); DanglingEnd();
}
