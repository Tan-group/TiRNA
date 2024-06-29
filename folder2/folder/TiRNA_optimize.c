/************* 这个程序是用来优化RNA结构 *****************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#define step_angle 0.8     //Angle step in Pivot_angle move in Folding
#define pi 3.1415
#define step1 .5           //Translation step length in Pivot move in Folding
#define step1_1 .5         //Translation step length in 3nt fragment move in Folding
#define step2 1.0          //One atom move step in Optimization
#define step3 .06          //All atoms move step in Optimizatio
#define total1  500000      //Max steps in Folding annealing (Structure predicton)        
#define Alpha 0.992        //Annealing rate
#define Anneal 1           //Annealing shecdle
#define A -3.5            //base-pairing strength of GC
//#define B0 -11.5           //conformational entropy in Base-stacking potential
#define B01 -9.6           //conformational entropy of base-stacking in CoaxialStacking potential
//#define ht 35
#define tran 10            //rotational frequency
#define BetaGC 0.95
#define BetaGU 0.60          //AU=Beta*GC
#define BetaAU 0.60
#define Beta1 0.0          //GA(UU)=Beta1*GC
#define Beta2 0.0          //GG=Beta2*GC
#define Beta3 0.0          //??=Beta3*GC  For non-canoncial base-pairing
#define gamma 0.10         //The energy threshold of base-pairing for statistics (碱基配对的能量Ubp0<gamma*Ubp算作配对形成)
#define tconf1 500       //The frequency of conformational output
//#define trmsd 5000       //类似MD，计算相邻trmsd的构象间的差别
#define tenergy 500     //The frequency of energy calculation
#define tbp 10          //The frequency of output information for base-pairing
#define tG 0            //计算平均base-pair数量的起始
#define tG_s 100           //Calculate pbp   经简单测试貌似还是100靠谱，设置为10000并没有节省计算时间也没能更好的呈现平衡
#define fp_bp 500         //The frequency of output the value of BP & PBP & G
#define t_si 5000          //Output frequency of structure information
#define tMove 10           //The frequency of 3-nt fragment moves；
// Electrostatic
//#define Ek 78.0          // debye parameters
#define bl 5.45
//#define RG 1.1           //TBI,the Rg of RNAs comparing with Rg of A-form helix;
#define IMg 3.0            //2:2 IMg=4.0; 2:1 IMg=3.0;
#define ww 0.5
#define N_thread 16
float T,Kd,I,Ek,q4,CNa,CMg,fNa,D;//9
int N,N0,N1,Nbp=15,t,salt,nm=0,Energy,Bulge,Pseudoknot;//14
int s[500][500],ss[500][500],a[10000],c[10000],bp,bp0,BP0,BP;//6
float rand01,phi,theta,rm,step;//7
float t0;//24
float x[10000],y[10000],z[10000],xx[10000],yy[10000],zz[10000],xxx[10000],yyy[10000],zzz[10000],Q[1000],R[1000],f[1000];//12
float U,Umin,Up,Up0,du,uu,u,ulj,uulj,uc,uuc,ub,uub,ue,uue,ud,uud,uN,uuN,us,uus,uco,uuco;//23
char  type[1000];//1
/***************new additional movement****************/
float rx0[1000],ry0[1000],rz0[1000],rR0[1000],rQ0[1000],rf0[1000];
char  rtype0[1000];
int   thread_num;
/********The parameters (including bonded & nonbonded) of the CG Model************/
float e=0.26;    //Strength in Uexc(LJ potential) 
float qq4,frac[1000],fi[1000],frac0[1000],frac1[1000],frac2[1000],BBL,av,ed;
float B0,a0,b0;
int Move;
/***************************************parallel********************************************/
int    rt,rrt[N_thread],rNCG[N_thread],total[N_thread],cca[N_thread],ccb[N_thread],ccc[N_thread],ccd[N_thread];
float rx[N_thread][1000],ry[N_thread][1000],rz[N_thread][1000],rR[N_thread][1000],rQ[N_thread][1000],rf[N_thread][1000],RU[N_thread];
float rt0[15]={25.0,31.0,37.0,43.0,49.0,55.0,63.0,71.0,78.0,86.0,94.0,102.0,110.0,120.0,130.0};
char  rtype[N_thread][1000];
int rp;
int yi1[N_thread],yj1[N_thread],yi2[N_thread],yj2[N_thread],yi3[N_thread],yj3[N_thread],yi4[N_thread],yj4[N_thread],yi5[N_thread],yj5[N_thread];
int fi1[N_thread],fj1[N_thread],fi2[N_thread],fj2[N_thread];
int junction,condition;
float thetaxx,thetaxx1,thetaxx2,thetaxx3,thetaxx4;
int xi1,xj1,xi2,xj2,xi3,xj3,xi4,xj4,xi5,xj5;
int si1,sj1,si2,sj2;
float utrs,uutrs;
int cycle,Nstrength;
float NS1,NS2,Bs[30];
float length1,length2,length3,length4,length5;
FILE *fpconf[N_thread],*fpsec_struc[N_thread],*fpenergy[N_thread],*fp,*Infor,*fpcfig;

#pragma omp threadprivate(t0)
#pragma omp threadprivate(N0,N1)
#pragma omp threadprivate(t,step,salt)
#pragma omp threadprivate(T,Kd,I,Ek,q4,D,fNa,CNa,CMg)
#pragma omp threadprivate(phi,theta,rm,rand01)
#pragma omp threadprivate(U,ulj,uulj)
#pragma omp threadprivate(x,y,z,xx,yy,zz,xxx,yyy,zzz,Q,R,f,type)
#pragma omp threadprivate(a,c,s,ss)
#pragma omp threadprivate(u,uu,ub,ue,ud,uN,us,uc,uub,uue,uud,uuN,uus,uuc,uco,uuco,du)
#pragma omp threadprivate(Energy,bp,bp0,BP,BP0)
#pragma omp threadprivate(Bulge,Pseudoknot,nm)
#pragma omp threadprivate(thread_num,B0,a0,b0)
#pragma omp threadprivate(frac,fi,frac0,frac1,frac2,BBL,av,ed,Move)
#pragma omp threadprivate(junction,condition,thetaxx,thetaxx1,thetaxx2,thetaxx3,thetaxx4)
#pragma omp threadprivate(utrs,uutrs,cycle)
#pragma omp threadprivate(xi1,xj1,xi2,xj2,xi3,xj3,xi4,xj4,xi5,xj5,si1,sj1,si2,sj2,Nstrength,NS1,NS2)
#pragma omp threadprivate(length1,length2,length3,length4,length5,Bs)
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
    	printf(">>>>>>>>>> The optimization process over!\n");
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
                  	for(j=0;j<N_thread;j++)
                  	{
                    		rx[j][i]=rx0[i];ry[j][i]=ry0[i];rz[j][i]=rz0[i];rR[j][i]=rR0[i];rQ[j][i]=rQ0[i];rf[j][i]=rf0[i];rtype[j][i]=rtype0[i];rrt[j]=rt;
                    		rNCG[j]=rp;
                  	}     
              	}
    	}
    	else 
    	{
             	for(j=0;j<N_thread;j++)
             	{
                      	rrt[j]=rt;
             	}
    	}
    	#pragma omp parallel num_threads(N_thread) 
    	{
             	thread_num=omp_get_thread_num(); 
             	replica();

    	}
}
/*********************************************************************************************************/
void detailed()
{       int k;
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
      	while(ri<N_thread-1)
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
             	}         
             	else
             	{
                    	for(i=1;i<=rp;i++)
                    	{
                                 rx[ri][i]=rx[ri][i];rx[rj][i]=rx[rj][i];
                                 ry[ri][i]=ry[ri][i];ry[rj][i]=ry[rj][i];
                                 rz[ri][i]=rz[ri][i];rz[rj][i]=rz[rj][i];
                    	}
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
      	while(ri<N_thread-1)
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
             	}         
             	else
             	{
                    	for(i=1;i<=rp;i++)
                    	{
                                rx[ri][i]=rx[ri][i];rx[rj][i]=rx[rj][i];
                                ry[ri][i]=ry[ri][i];ry[rj][i]=ry[rj][i];
                                rz[ri][i]=rz[ri][i];rz[rj][i]=rz[rj][i];
                    	}
             	}
            	ri=ri+2;             
        }
}
void replica()
{
 	int i;
 	void Put_File(),Fixed_Atom(),MC_Annealing();
 	Put_File();         //defined the input parameters and out file names;
 	N0=rNCG[thread_num];
 	N=(N0-1)/3;          //N0:Total number of CG beads; N: Num. of nt
 	Fixed_Atom();      //The Centre atom N1 will be fixed;
 	for(i=1;i<=N0;i++) 
  	{
            	x[i]=rx[thread_num][i];	y[i]=ry[thread_num][i];	z[i]=rz[thread_num][i];
            	type[i]=rtype[thread_num][i];	Q[i]=rQ[thread_num][i];	f[i]=rf[thread_num][i];
            	R[i]=rR[thread_num][i];
  	}
  	xi1=yi1[thread_num]; xj1=yj1[thread_num]; 
  	xi2=yi2[thread_num]; xj2=yj2[thread_num]; 
  	xi3=yi3[thread_num]; xj3=yj3[thread_num]; 
  	xi4=yi4[thread_num]; xj4=yj4[thread_num];
  	xi5=yi5[thread_num]; xj5=yj5[thread_num];
  	si1=fi1[thread_num]; sj1=fj1[thread_num]; 
  	si2=fj1[thread_num]; sj2=fj2[thread_num];   
  	float angle=10;
  	angle=3.0/N_thread;
  	thetaxx=0.3+thread_num*angle;
  	for(i=1;i<=N0;i++)
 	{
 		if(xi1==0)
            	{
            		junction=1;
            	}
            	else if(xi1!=0&&xi2==0)
            	{
            		junction=2;
            	}
 		else if(xi3!=0&&xi4==0)
 		{
 			junction=3;
 			if((xi2-xi1)>(xj1-xj3))	{thetaxx1=thetaxx;thetaxx3=2.4;thetaxx2=2.4;}
 			else			{thetaxx3=thetaxx;thetaxx1=2.4;thetaxx2=2.4;}			 	
 			length1=xi1+N0-xj1+1;
 			length2=xj2-xi2+1;
 			length3=xj3-xi3+1;
    		}
    		else if(xi4!=0&&xi5==0)
    		{
    			junction=4;
    			if((xi2-xi1)>=(xj1-xj4)){condition=1;}
    			else			{condition=2;}
    			thetaxx1=thetaxx-0.3;    			 	
 			length1=xi1+N0-xj1+1;
 			length2=xj2-xi2+1;
 			length3=xj3-xi3+1;
 			length4=xj4-xi4+1;
    		}
    		else if(xi5!=0)
    		{
    			junction=5;
    			if((xi2-xi1)>(xj1-xj5))	{condition=1;}
    			else			{condition=2;}
    			thetaxx4=0.0;
    			thetaxx1=thetaxx-0.3;    				
			length1=xi1+N0-xj1+1;
			length2=xj2-xi2+1;
			length3=xj3-xi3+1;
			length4=xj4-xi4+1;
			length5=xj5-xi5+1;
    		}
    		else
    		{
    			junction=0;
    		}
 	}
 	if(N<=80&&junction<3) 	{NS1=1.0;NS2=1.0;}
 	else	  		{NS1=0.4;NS2=1.6;}
  	MC_Annealing();
  	for(i=1;i<=N0;i++)
  	{
            	rx[thread_num][i]=x[i];	ry[thread_num][i]=y[i];	rz[thread_num][i]=z[i];
  	}
  	RU[thread_num]=U;
 // rt0[thread_num]=t0;
}  
void openfile()
{

	FILE *fpb,*fpb1;
     	int i,duo1,duo2;
    	char filename[20];
     	Infor=fopen("Information.dat","w+");
     	fp=fopen("ch_0.dat","r+");
     	fpb=fopen("stem.dat","r+");
     	fpb1=fopen("stem_kissing.dat","r+");
     	fpcfig=fopen("config1.dat","r+");
     	int ya1,ya2,ya3,ya4,ya5,ya6,ya7,ya8,ya9,ya10,ya11,ya12,ya13,ya14;
     	i=1;
     	while(!feof(fp)) 
     	{
         	fscanf(fp,"%d %d %s %f %f %f %f %f %f\n",&duo1,&duo2,&rtype0[i],&rx0[i],&ry0[i],&rz0[i],&rR0[i],&rQ0[i],&rf0[i]); 
         	i++;
     	}      
     	fclose(fp); 
     	rp=i-1; 
     	while(!feof(fpb))
 	{
      		fscanf(fpb,"%d %d %d %d %d %d %d %d %d %d\n",&ya1,&ya2,&ya3,&ya4,&ya5,&ya6,&ya7,&ya8,&ya9,&ya10);
 	}
 	fclose(fpb);
 	while(!feof(fpb1))
 	{
      		fscanf(fpb1,"%d %d %d %d\n",&ya11,&ya12,&ya13,&ya14);
 	}
 	fclose(fpb1);
 	for(i=0;i<N_thread;i++)
 	{
 		yi1[i]=ya1*3;yj1[i]=ya2*3;
 		yi2[i]=ya3*3;yj2[i]=ya4*3;
 		yi3[i]=ya5*3;yj3[i]=ya6*3;
 		yi4[i]=ya7*3;yj4[i]=ya8*3;
 		yi5[i]=ya9*3;yj5[i]=ya10*3;
 		fi1[i]=ya11*3;fj1[i]=ya12*3;
 		fi2[i]=ya13*3;fj2[i]=ya14*3;
 	}
     	//printf("CG beads number:   %d,  the thread number:   %d\n",rp,N_thread);
     	fprintf(Infor,"CG beads number:   %d,  the thread number:   %d\n",rp,N_thread);fflush(Infor);
     	for (i=0;i<N_thread;i++) 
     	{
           	sprintf(filename,"conf_%d.dat",i); fpconf[i]=fopen(filename,"w+");
           	sprintf(filename,"second_stru_%d.dat",i); fpsec_struc[i]=fopen(filename,"w+");
           	//sprintf(filename,"Energy_%d.dat",i);    fpenergy[i]=fopen(filename,"w+");
      	}  
      	int ca,cb,cc,cd,ce,ct;  
      	while(!feof(fpcfig))
      	{
      		fscanf(fpcfig,"%d %d %d %d %d %d\n",&ct,&ca,&cb,&cc,&cd,&ce);
      	} 
        for(i=0;i<N_thread;i++)
      	{
      		total[i]=cb;ccc[i]=cc;ccd[i]=cd;
      	}   
      	fclose(fpcfig);      
}


void closefile()
{
         int i;
         for(i=0;i<N_thread;i++)
         {
         	fclose(fpconf[i]); fclose(fpsec_struc[i]);// fclose(fpenergy[i]);
         }         
         fclose(Infor);
        
}

/* &%$#@!~&%$#@!~&%$#@!~   Some functions or modules for move and calculation    &%$#@!~&%$#@!~&%$#@!~ */
/*********************读入文件，读出文件，输入参数************************/
void Put_File (void) 
{
 	CNa=ccc[thread_num];CMg=ccd[thread_num];
 	if(CNa==0.0&&CMg==0.0)	{salt=0;}
 	else			{salt=1;}
  	t0=rt0[thread_num];
   	//printf("thread:  %d,  CNa = %f,  CMg = %f, t0 = %f\n",thread_num,CNa,CMg,t0);
   	fprintf(Infor,"thread:  %d,  CNa = %f,  CMg = %f, t0 = %f\n",thread_num,CNa,CMg,t0);fflush(Infor);
}
/*********************************************/
void Fixed_Atom(void)
{
   	int N10;
   	N10=floor(N0/2)+1; 
  	if 	(fmod(N10,3)==0) 	{N1=N10-2;}
   	else if (fmod(N10+1,3)==0) 	{N1=N10-1;}
   	else 				{N1=N10;}                 //N1:Fixed the P in centre of chain
 //  printf("N %d N0 %d N1 %d b0 %f\n",N,N0,N1,b0);  //The number of nucleotides & atoms & the actionless atom
}
//******************************************//
void Parameters_T()
{
    	float tt,TT;
    	void Bs_stacking();
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
            	if (qqq>1.){ q4=qq4*qqq; }  
            	else       { q4=qq4;     }
     	}
    	else 
        {
           	q4=qq4;
     	} // q4=b/lB &稀溶液修正 
	TT=T-273.15;
     	if(TT>=50&&TT<=90)        {B0=-0.068*TT-6.5;}
     	else if(TT<50)            {B0=-0.068*50-6.5;}
     	else                      {B0=-11.5;}
        Bs_stacking();
}

void Bs_stacking()
{
	Bs[1]=-6.82-T*0.001*(-19.0-B0);  // AA/UU
	Bs[2]=-11.4-T*0.001*(-29.5-B0);  // AC/UG
	Bs[3]=-10.48-T*0.001*(-27.1-B0); // AG/UC
	Bs[4]=-3.21-T*0.001*(-8.6-B0);  // AG/UU
	Bs[5]=-9.38-T*0.001*(-26.7-B0);  // AU/UA
	Bs[6]=-8.81-T*0.001*(-24.0-B0); // AU/UG
	Bs[7]=-10.44-T*0.001*(-26.9-B0);// UG/AC
	Bs[8]=-13.39-T*0.001*(-32.7-B0); // GG/CC
	Bs[9]=-10.64-T*0.001*(-26.7-B0);// CG/GC
	Bs[10]=-5.61-T*0.001*(-13.5-B0); // CG/GU
	Bs[11]=-12.11-T*0.001*(-32.2-B0);// CU/GG
	Bs[12]=-12.44-T*0.001*(-32.5-B0);// UC/AG
	Bs[13]=-14.88-T*0.001*(-36.9-B0);// GC/CG
	Bs[14]=-8.33-T*0.001*(-21.9-B0); // GG/CU
	Bs[15]=-12.59-T*0.001*(-32.5-B0); // GU/CG
	Bs[16]=-7.69-T*0.001*(-20.5-B0);  // UA/AU
	Bs[17]=-6.99-T*0.001*(-19.3-B0); // UG/AU
	Bs[18]=-12.83-T*0.001*(-37.3-B0);// UU/AG
	Bs[19]=-13.47-T*0.001*(-44.9-B0);// UU/GG
	Bs[20]=-9.26-T*0.001*(-30.8-B0); // UG/GU
	Bs[21]=-16.66-T*0.001*(-50.3-B0);
	Bs[22]=-14.59-T*0.001*(-51.2-B0);
}
/****************Monte Carlo simulated annealing*********************/
void MC_Annealing()
{
     	void MC_T(),Parameters_T();               
     	Parameters_T();     //parameters at any t0: steps, T, Debye length, ionic strength, ion fraction etc.
     	MC_T();
}
/*****************Monte Carlo simulation at given Temperature******************************/
void MC_T(void)
{
     	int ii,i;
     	void MC_Each_Step(),ENERGY(),disbrute_qq();
     	srand((unsigned)time(NULL)); 
     	for (t=1;t<=total[thread_num];t++)
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
           		cycle=ii;
                  	MC_Each_Step();
           	}
           	if (fmod(t,tconf1)==0) 
           	{                   
                	for(i=1;i<=N0;i++) 
                	{
                      		fprintf(fpconf[thread_num],"%d %d %c %f %f %f %f %f %f\n",t,i,type[i],x[i],y[i],z[i],R[i],Q[i],f[i]);
                	}
                	fflush(fpconf[thread_num]);                           
            	}              
            	ENERGY();                //Calculate energy of one conformation                                  
      		if(fmod(t,1000)==0)
      		{	 
             		//printf("thread: %d,  step:%6d,     T:%5.1f,    Energy: %f kcal/mol\n",thread_num,t,t0,U);
             		fprintf(Infor,"thread: %d,  step:%6d,     T:%5.1f,    Energy: %f kcal/mol\n",thread_num,t,t0,U);fflush(Infor); 
      		}   
    	}   
}
/*****************Monte Carlo for each step******************************/
void MC_Each_Step()
{
 	int   i0;
 	float PC(),CP(),CN(),PCP(),CPC(),PCN(),NCP(),PCPC(),CPCP(),PCPCh(),CPCPh(),CPCN(),NCPC(),LJ0(),HB(),St(),QQ(),QQ11();
 	void MoveN(),Rand01(),Translate(),Pivot(),RMSDconf(),FOLD();
 	void BasePairing(),BaseStacking(),ExcludedVolume(),Electrostatic(),Bonded(),CoaxialStacking();
 	void Metropolis(); 
   	u=0.0;uu=0.0;                                    //Initialize the energies;
   	ulj=0.0;ub=0.0;ue=0.0;ud=0.0;uN=0.0;us=0.0;uc=0.0;
   	uulj=0.0;uub=0.0;uue=0.0;uud=0.0;uuN=0.0;uus=0.0;uuc=0.0;
   	uco=0.0;uuco=0.0;
   	utrs=0.0;uutrs=0.0;
   	Energy=0;Move=0;                      //if Energy=1, fuction ENERGY() is running 
   	ed=0.0;
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
   	u=ulj+uN*Nstrength+us+uc+(ub*NS1+ue*NS1+ud*NS2)*0.5963+uco+utrs;
   	uu=uulj+uuN*Nstrength+uus+uuc+(uub*NS1+uue*NS1+uud*NS2)*0.5963+uuco+uutrs;
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
       	xxx[i1]=x[i1]+step*rand01*sin(phi)*cos(theta);
       	yyy[i1]=y[i1]+step*rand01*sin(phi)*sin(theta);
       	zzz[i1]=z[i1]+step*rand01*cos(phi);
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
     	if(c[i1]==1&&c[i1-3]==1&&a[i1]==1&&a[i1-3]==1)	{ul=133.4*(d-3.90)*(d-3.90);}
     	else						{ul=63.2*(d-3.92)*(d-3.92);}
     	return ul;
}
float CP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d,ul;
     	d=sqrt((x1[i1+1]-x1[i1-1])*(x1[i1+1]-x1[i1-1])+(y1[i1+1]-y1[i1-1])*(y1[i1+1]-y1[i1-1])+(z1[i1+1]-z1[i1-1])*(z1[i1+1]-z1[i1-1]));
     	if(c[i1]==1&&c[i1+3]==1&&a[i1]==1&&a[i1+3]==1)	{ul=107.8*(d-3.87)*(d-3.87);}
     	else						{ul=32.6*(d-3.90)*(d-3.90);}
     	return ul;
}
float CN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d,ul;
     	d=sqrt((x1[i1]-x1[i1-1])*(x1[i1]-x1[i1-1])+(y1[i1]-y1[i1-1])*(y1[i1]-y1[i1-1])+(z1[i1]-z1[i1-1])*(z1[i1]-z1[i1-1]));
     	if(c[i1]==1&&a[i1]==1)	{ul=40.0*(d-3.36)*(d-3.36);}
     	else			{ul=24.5*(d-3.36)*(d-3.36);}
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
      	if(c[i1]==1&&a[i1]==1)	{ue0=28.3*(a1-1.70)*(a1-1.70);}
      	else			{ue0=9.7*(a1-1.82)*(a1-1.82);}
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
        if(c[i1]==1&&a[i1]==1)	{ue0=63.9*(a1-1.70)*(a1-1.70);}
        else			{ue0=15.3*(a1-1.80)*(a1-1.80);}
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
      	if(c[i1]==1&&a[i1]==1)	{ue0=35.5*(a1-1.63)*(a1-1.63);}
      	else			{ue0=10.2*(a1-1.64)*(a1-1.64);}
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
     	if(c[i1]==1&&a[i1]==1)	{ue0=99.8*(a1-1.66)*(a1-1.66);}
     	else			{ue0=15.4*(a1-1.66)*(a1-1.66);}
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
 	if(c[i1]==1&&a[i1]==1)	{ud0=2.68*((1-cos(di-2.51))+0.5*(1-cos(3.*(di-2.51))));}
 	else			{ud0=0.8*((1-cos(di-2.56))+0.5*(1-cos(3.*(di-2.56))));}
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
     	if(c[i1]==1&&a[i1]==1)	{ud0=10.5*((1-cos(di+2.92))+0.5*(1-cos(3.*(di+2.92))));}
     	else			{ud0=3.6*((1-cos(di+2.92))+0.5*(1-cos(3.*(di+2.92))));}
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
     	if(c[i1]==1&&a[i1]==1)
     	{
     		ud0=3.8*((1-cos(di+1.18))+0.5*(1-cos(3.*(di+1.18))));
     	}
     	else
     	{
     		ud0=0.8*((1-cos(di+1.16))+0.5*(1-cos(3.*(di+1.16))));
     	}
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
   	if(c[i1]==1&&a[i1]==1)	{ud0=4.2*((1-cos(di-0.77))+0.5*(1-cos(3.*(di-0.77))));}
   	else				{ud0=0.8*((1-cos(di-0.68))+0.5*(1-cos(3.*(di-0.68))));}
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
 void BasePairing()
{
	int i,j;
    	float uN1=0.0,uuN1=0.0;
    	bp=0;BP=0;
    	for (i=1;i<=N0;i++)  
    	{
            	a[i]=0;c[i]=0; 
            	for(j=i+3;j<=N0;j++) 
            	{
                      s[i][j]=0;ss[i][j]=0;
            	}
     	}
     	for (i=3;i<=(N0-12);i=i+3) 
     	{   	    
      		for(j=i+12;j<=N0;j=j+3)
             	{          		     
                    	if ((type[i]=='G'&&type[j]=='C')||(type[i]=='C'&&type[j]=='G')||(type[i]=='A'&&type[j]=='U')||(type[i]=='U'&&type[j]=='A')||(type[i]=='G'&&type[j]=='U')||(type[i]=='U'&&type[j]=='G'))
                    	{
                        	uN1=0.0;uuN1=0.0; bp0=0; BP0=0;
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
                                  		a[i]=1;a[j]=1; ss[i][j]=1;
                              		}
                           	}
                           	uuN=uuN+uuN1; 
                      	} 
                      	else
                      	{
                          	uN1=0.0;uuN1=0.0;
                          	uN=uN+uN1; 
                          	uuN=uuN+uuN1;
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
     	if 	(((type[i1]=='A'&&type[k1]=='A')&&(type[j1]=='U'&&type[k2]=='U'))
		||((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='A'&&type[k2]=='A'))) 	{B=Bs[1];}  // AA/UU
	else if (((type[i1]=='A'&&type[k1]=='C')&&(type[j1]=='U'&&type[k2]=='G')) 
     		||((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='C'&&type[k2]=='A'))) 	{B=Bs[2];}  // AC/UG
	else if (((type[i1]=='A'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='C'))
     		||((type[i1]=='C'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='A'))) 	{B=Bs[3];} // AG/UC
	else if (((type[i1]=='A'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='A'))) 	{B=Bs[4];}   // AG/UU
	else if ((type[i1]=='A'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='A'))   	{B=Bs[5];}  // AU/UA
	else if (((type[i1]=='A'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='A'))) 	{B=Bs[6];}  // AU/UG
	else if (((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='A'&&type[k2]=='C')) 
     		||((type[i1]=='C'&&type[k1]=='A')&&(type[j1]=='G'&&type[k2]=='U'))) 	{B=Bs[7];} // UG/AC
	else if (((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='C'&&type[k2]=='C'))
     		||((type[i1]=='C'&&type[k1]=='C')&&(type[j1]=='G'&&type[k2]=='G'))) 	{B=Bs[8];} // GG/CC
	else if ((type[i1]=='C'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='C'))   	{B=Bs[9];} // CG/GC
	else if (((type[i1]=='C'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='C'))) 	{B=Bs[10];}  // CG/GU
	else if (((type[i1]=='C'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='C'))) 	{B=Bs[11];} // CU/GG
	else if (((type[i1]=='U'&&type[k1]=='C')&&(type[j1]=='A'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='A')&&(type[j1]=='C'&&type[k2]=='U'))) 	{B=Bs[12];} // UC/AG
	else if ((type[i1]=='G'&&type[k1]=='C')&&(type[j1]=='C'&&type[k2]=='G'))   	{B=Bs[13];} // GC/CG
	else if (((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='C'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='C')&&(type[j1]=='G'&&type[k2]=='G'))) 	{B=Bs[14];}  // GG/CU
	else if (((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='C'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='C')&&(type[j1]=='U'&&type[k2]=='G'))) 	{B=Bs[15];} // GU/CG
	else if ((type[i1]=='U'&&type[k1]=='A')&&(type[j1]=='A'&&type[k2]=='U'))   	{B=Bs[16];}  // UA/AU
	else if (((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='A'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='A')&&(type[j1]=='G'&&type[k2]=='U'))) 	{B=Bs[17];}  // UG/AU
	else if (((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='A'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='A')&&(type[j1]=='U'&&type[k2]=='U'))) 	{B=Bs[18];} // UU/AG
	else if (((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='U'))) 	{B=Bs[19];} // UU/GG
	else if ((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='U'))   	{B=Bs[20];}  // UG/GU
	else   if ((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='G'))   	
	{
        	if 	((type[i1-3]=='G'&&type[k1+3]=='C')&&(type[j1+3]=='C'&&type[k2-3]=='G'))	{B=Bs[21];}
        	else                                                                       		{B=Bs[22];}
	}// GU/UG (CGUC/CUGG)
	else B=0.0;
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
    	for(i=3;i<=(N0-12);i=i+3)
    	{  
      		us0=0.0;uus0=0.0;
          	for(j=i+12;j<=N0;j=j+3)
          	{   
              		us1=0.0;uus1=0.0; 
              		us2=0.0; uus2=0.0;
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
            	us=us+us0; 
            	uus=uus+uus0;
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
/*   else if (fmod(i1+2,3)==0&&fmod(j1+2,3)==0&&abs(j1-i1)>12) {r1=7.8;r2=7.8;}*/
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


void ExcludedVolume(int i0)     //Exculed Volume
{
     	int i,j,i01,i02;
	float ulj0,uulj0,ulj1,uulj1;
	if (Move==3)           //when 3nt fragment moves
     	{
       		for(i=i0+1;i<=i0+8;i++) 
      		{
        		ulj0=0.0;uulj0=0.0;
        		for (j=1;j<=N0;j++)	  
        		{
         			ulj1=0;uulj1=0;
         			if (j<=i0||j>i0+8)
         			{
         				ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                                
         				uulj1=LJ0(i,j,xx,yy,zz); uulj0=uulj0+uulj1;           
          			}
         		}
        		ulj=ulj+ulj0; uulj=uulj+uulj0;
         	}
      	}
	if(Move==4)
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
      	if(Move==5)
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
      	if (Move==1)       //when Pivot moves
     	{
      		if (fmod(i0,3)!=0)
      		{
       			if (i0<N1) {i01=i0+1; i02=i0;}
       			else       {i01=i0; i02=i0-1;} 
       			for(i=i01;i<=N0;i++) 
        		{
        			ulj0=0.0;uulj0=0.0;
        			for (j=1;j<=i02;j++)	  
         			{
         				ulj1=0;uulj1=0;
         				ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                                
         				uulj1=LJ0(i,j,xx,yy,zz); uulj0=uulj0+uulj1;           
          			}
        			ulj=ulj+ulj0; uulj=uulj+uulj0;
         		}
        	}
     	 	else //这里的计算量有些浪费，期待更改:直接去掉base的单独运动即可～
       		{	  
       			i=i0;  ulj0=0.0;uulj0=0.0;
       			for (j=1;j<=N0;j++)	  
        		{ 
        			if (j!=i)
         			{  
        				ulj1=0;uulj1=0;
        				ulj1=LJ0(i,j,x,y,z); ulj=ulj+ulj1;                          
        				uulj1=LJ0(i,j,xx,yy,zz); uulj=uulj+uulj1;        
          			} 
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
    	int i,j,i01,i02;
	float uc0=0.0,uuc0=0.0,uc1=0.0,uuc1=0.0;
    	if (salt==0) {uc=0.0;uuc=0.0;}
    	else
    	{
      		if (Move==3)
     		{
       			for(i=i0;i<=i0+8;i++) //N0-3:the last P is uncharged
       			{
        			uc0=0.0;uuc0=0.0;
        			if (fmod(i+2,3)==0)
        			{
        				for (j=4;j<=N0-3;j++)  //4:the first P is still uncharged  
         				{
         					uc1=0.0;uuc1=0.0;
         					if ((j<i0||j>i0+8)&&fmod(j+2,3)==0)
          					{
         						uc1=QQ(i,j,x,y,z); uc0=uc0+uc1;
         						uuc1=QQ(i,j,xx,yy,zz); uuc0=uuc0+uuc1;            
           					}
          				}
        				uc=uc+uc0; uuc=uuc+uuc0;
         			}
        		}
      		}
		if (Move==4)
     		{
       			for(i=i0;i<=i0+11;i++) //N0-3:the last P is uncharged
       			{
        			uc0=0.0;uuc0=0.0;
        			if (fmod(i+2,3)==0)
        			{
        				for (j=4;j<=N0-3;j++)  //4:the first P is still uncharged  
         				{
         					uc1=0.0;uuc1=0.0;
         					if ((j<i0||j>i0+11)&&fmod(j+2,3)==0)
          					{
         						uc1=QQ(i,j,x,y,z); uc0=uc0+uc1;
         						uuc1=QQ(i,j,xx,yy,zz); uuc0=uuc0+uuc1;            
           					}
          				}
        				uc=uc+uc0; uuc=uuc+uuc0;
         			}
        		}
      		}
		if (Move==5)
     		{
       			for(i=i0;i<=i0+14;i++) //N0-3:the last P is uncharged
       			{
        			uc0=0.0;uuc0=0.0;
        			if (fmod(i+2,3)==0)
        			{
        				for (j=4;j<=N0-3;j++)  //4:the first P is still uncharged  
         				{
         					uc1=0.0;uuc1=0.0;
         					if ((j<i0||j>i0+14)&&fmod(j+2,3)==0)
          					{
         						uc1=QQ(i,j,x,y,z); uc0=uc0+uc1;
         						uuc1=QQ(i,j,xx,yy,zz); uuc0=uuc0+uuc1;            
           					}
          				}
        				uc=uc+uc0; uuc=uuc+uuc0;
         			}
        		}
      		}
      		if (Move==1)
     		{
      			if (fmod(i0,3)==0) {uc=0.0;uuc=0.0;}
      			else 
      			{
       				if (i0<N1) {i01=i0+1; i02=i0;}
       				else       {i01=i0; i02=i0-1;} 
       				for(i=i01;i<=N0-3;i++) //N0-3:the last P is uncharged
       				{
        				uc0=0.0;uuc0=0.0;
        				if (fmod(i+2,3)==0)
        				{
        					for (j=4;j<=i02;j++)  //4:the first P is still uncharged  
         					{
         						uc1=0.0;uuc1=0.0;
         						if (fmod(j+2,3)==0)
          						{
         							uc1=QQ(i,j,x,y,z); uc0=uc0+uc1;
         							uuc1=QQ(i,j,xx,yy,zz); uuc0=uuc0+uuc1;            
           						}
          					}
        					uc=uc+uc0; uuc=uuc+uuc0;
         				}
        			}
       			}
      		}
     	}
}
// Bonded potential before and after conformational change: bond, angle and dihedral

/***********************************************************************************************************/

float QQ11(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
 	float d,ddao;
  	d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); 
  	ddao=1.0/d;
  	return ddao;
}

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
    			ub=PC(i0+2,x,y,z)+CP(i0+8,x,y,z); uub=PC(i0+2,xx,yy,zz)+CP(i0+8,xx,yy,zz);
    			ue=CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+PCN(i0+2,x,y,z)+PCP(i0+8,x,y,z)+CPC(i0+8,x,y,z)+NCP(i0+8,x,y,z);
    			uue=CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+PCN(i0+2,xx,yy,zz)+PCP(i0+8,xx,yy,zz)+CPC(i0+8,xx,yy,zz)+NCP(i0+8,xx,yy,zz);
    			ud=PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+CPCN(i0+2,x,y,z)+CPCP(i0+5,x,y,z)+PCPC(i0+8,x,y,z)+CPCP(i0+8,x,y,z)+NCPC(i0+8,x,y,z)+CPCN(i0+11,x,y,z);
    			uud=PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+CPCN(i0+2,xx,yy,zz)+CPCP(i0+5,xx,yy,zz)+PCPC(i0+8,xx,yy,zz)+CPCP(i0+8,xx,yy,zz)+NCPC(i0+8,xx,yy,zz)+CPCN(i0+11,xx,yy,zz);
    		}
  		if (fmod(i0+1,3)==0)
  		{
    			ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z)+PC(i0+10,x,y,z); uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz)+PC(i0+10,xx,yy,zz);
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
	else
 	{
   		if (fmod((i0+1),3)==0)       // Choose C atoms
  		{
       			if (i0<N1)
       			{
       				ub=CP(i0+1,x,y,z)+CN(i0+1,x,y,z);  uub=CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz);
       				ue=NCP(i0+1,x,y,z)+PCP(i0+1,x,y,z)+CPC(i0+1,x,y,z)+PCN(i0+1,x,y,z);
       				uue=NCP(i0+1,xx,yy,zz)+PCP(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz);
        		}
       			else
       			{
       				ub=PC(i0+1,x,y,z);           uub=PC(i0+1,xx,yy,zz);
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
        		if (i0<N1) {ub=PC(i0+2,x,y,z);uub=PC(i0+2,xx,yy,zz);}
        		else       {ub=CP(i0-1,x,y,z);uub=CP(i0-1,xx,yy,zz);}
        		if (i0==1) 
         		{
        			ue=PCN(i0+2,x,y,z)+PCP(i0+2,x,y,z);    uue=PCN(i0+2,xx,yy,zz)+PCP(i0+2,xx,yy,zz); 
        			ud=PCPC(i0+2,x,y,z);                   uud=PCPC(i0+2,xx,yy,zz);
         	 	}
        		else if (i0==N0) 
         		{
        			ue=NCP(i0-1,x,y,z)+PCP(i0-1,x,y,z); uue=NCP(i0-1,xx,yy,zz)+PCP(i0-1,xx,yy,zz);
        			ud=CPCP(i0-4,x,y,z);                uud=CPCP(i0-4,xx,yy,zz);
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
        		ub=CN(i0,x,y,z);                 uub=CN(i0,xx,yy,zz);
        		ue=NCP(i0,x,y,z)+PCN(i0,x,y,z);  uue=NCP(i0,xx,yy,zz)+PCN(i0,xx,yy,zz); 

        		if (i0==3) {ud=NCPC(i0,x,y,z); uud=NCPC(i0,xx,yy,zz);}
        		else if (i0==N0-1) {ud=CPCN(i0,x,y,z); uud=CPCN(i0,xx,yy,zz);}
        		else       {ud=NCPC(i0,x,y,z)+CPCN(i0,x,y,z); uud=NCPC(i0,xx,yy,zz)+CPCN(i0,xx,yy,zz);}
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
             	dco_1=5.0;   dco_2=2*5.0; 
             	kdco_1=2.5; kdco_2=2*2.5; 
    	}
    	else 
    	{
             	dco_1=5.0;  dco_2=5.0; 
             	kdco_1=2.5; kdco_2=2.5; 
    	}   // 针对bulge做修正，如果间隔3nt，则最优距离dco=2dco，kdco=2kdco；
     // if (fmod(t,100000)==0) {printf("%f %f %f %f %f %f\n",dco1,dco_1,kdco_1,dco2,dco_2,kdco_2);}
   	Ucos1=-0.5*(kco-3.0)*((1-exp(-(dco1-dco_1)/kdco_1))*(1-exp(-(dco1-dco_1)/kdco_1))-1);
  	Ucos2=-0.5*(kco-3.0)*((1-exp(-(dco2-dco_2)/kdco_2))*(1-exp(-(dco2-dco_2)/kdco_2))-1);
   	Ucos=Ucos1+Ucos2; 
// if (fmod(t,tprint)==0&&Energy==1) printf("%d %d %d %d %f %f %f %f %f %f\n",i1/3,j1/3,k1/3,k2/3,kco,dco1,dco2,Ucos1,Ucos2,Ucos);
   	return Ucos;
}

void CoaxialStacking()
{
   	int i,j,k1,k2;
   	float uco0,uuco0,uco1,uuco1;
   	for(i=9;i<=N0-16;i=i+3)
   	{
          	for(j=i+12;j<=N0-4;j=j+3)
          	{
              		if ((s[i][j]==1&&s[i-3][j+3]==1&&s[i-6][j+6]==1&&s[i+3][j-3]==0)||(ss[i][j]==1&&ss[i-3][j+3]==1&&ss[i-6][j+6]==1&&ss[i+3][j-3]==0))
              		{
                    		uco0=0.0; uuco0=0.0;
                    		for(k1=i+3;k1<j;k1=k1+3)
                    		{      
                            		for(k2=k1+12;k2<=N0-4;k2=k2+3)
                            		{
                                    		if (j>k2&&(abs(i-k1)!=3||abs(k2-j)!=3)&&((s[k1][k2]==1&&s[k1+3][k2-3]==1&&s[k1+6][k2-6]==1&&s[k1-3][k2+3]==0)||(ss[k1][k2]==1&&ss[k1+3][k2-3]==1&&ss[k1+6][k2-6]==1&&ss[k1-3][k2+3]==0)))
                                    		{
                                         		uco1=0.0; uuco1=0.0;
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
              	 		uco=uco+uco0; 
               			uuco=uuco+uuco0;
         		}
      		}
 	}  
}


// Energy of one confromation
void ENERGY()
{
   	int i,j;
   	float u0=0.0,uc1=0.0,ulj1=0.0,ulj0=0.0,uc0=0.0,uN1=0.0;
   	Energy=1;
   	if (fmod(t,tenergy)==0)
   	{
      		ulj=0.0;ub=0.0;uN=0.0;
      		us=0.0;U=0.0; uco=0.0;
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
            		uc=0.0;
            		uuc=0.0;
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
                                                           			fprintf(fpsec_struc[thread_num],"%d %c %d %c\n",i/3,type[i],j/3,type[j]);fflush(fpsec_struc[thread_num]);
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
                                    	u0=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z)+PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z)+PCPC(i,x,y,z)+CPCP(i,x,y,z)+NCPC(i,x,y,z);
                            	}
                            	else if (i==(N0-1))
                            	{
                                     	u0=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z)+PCP(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z)+CPCN(i,x,y,z);
                             	} 
                             	else 
                             	{
                                      	u0=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z)+PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z)+PCPC(i,x,y,z)+CPCP(i,x,y,z)+CPCN(i,x,y,z)+NCPC(i,x,y,z);
                              	}
                              	ub=ub+u0*0.5963;
                       	}
             	} 
            	if(fmod(t,tbp)==0) 
           	{
                    	fprintf(fpsec_struc[thread_num],"\n");fflush(fpsec_struc[thread_num]);
            	}
             	BaseStacking(); 
          	CoaxialStacking(); //Triplex(); //TerminalMismatch(); DanglingEnd();
           	U=ulj+ub+uN+us+uc+uco/*+uTri+umis+udang*/;
        	//fprintf(fpenergy[thread_num],"%d %f %f %f %f %f %f %f %f %f\n",t,U,Up,Umin,ulj,ub,uN,us,uc,uco);fflush(fpenergy[thread_num]);
     	}
}
//Out put some structural parameters and input value
void FOLD(int i0)
{
	int i;
	float prob=0.0;
	void constraint();
        step=step1; 
	prob=rand()/(RAND_MAX+1.0);
	if (prob<0.2&&i0>12&&i0<N0-12&&fmod(i0,3)!=0)    //3nt fragment translation;
    	{
      		Move=3;step=step1_1; 
      		for(i=1;i<=N0;i++)  
		{      //i0+1-->i0+8:translated
          		if (i>i0&&i-i0<=8) 
			{
				Translate(i);
				xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
			}
          		else  
			{
				xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
			}  
		}      
      	}
	else if (prob<0.4&&prob>=0.2&&i0>15&&i0<N0-15&&fmod(i0,3)!=0)    //4nt fragment translation;
    	{
      		Move=4;step=step1_1; 
      		for(i=1;i<=N0;i++)  
		{      //i0+1-->i0+8:translated
          		if (i>i0&&i-i0<=11) 
			{
				Translate(i);
				xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
			}
          		else  
			{
				xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
			}  
		}      
      	}
	else if (prob<0.7&&prob>=0.5&&i0>18&&i0<N0-18&&fmod(i0,3)!=0)    //5nt fragment translation;
    	{
      		Move=5;step=step1_1; 
      		for(i=1;i<=N0;i++)  
		{      //i0+1-->i0+8:translated
          		if (i>i0&&i-i0<=14) 
			{
				Translate(i);
				xx[i]=xxx[i];yy[i]=yyy[i];zz[i]=zzz[i];
			}
          		else  
			{
				xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
			}  
		}      
      	}
	else
 	{
		Move=1; step=step1; 
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
 	Electrostatic(i0);   	
	Bonded(i0);
  	BasePairing(); 
  	BaseStacking(); 
  	constraint();
     	CoaxialStacking();
}

void constraint()
{
    	float force_3(),force_4(),force_5();
	if(t<10000)	{Nstrength=100;}
	else		{Nstrength=1.0;}
    	utrs=0.0;uutrs=0.0;
    	if(junction==0)	{utrs=0.0;uutrs=0.0;}
    	if(junction==3)	{utrs=force_3(x,y,z);uutrs=force_3(xx,yy,zz);}
    	if(junction==4)	{utrs=force_4(x,y,z);uutrs=force_4(xx,yy,zz);}
    	if(junction==5)	{utrs=force_5(x,y,z);uutrs=force_5(xx,yy,zz);}	
    	if(fmod(t,5000)==0&&cycle==N0)	
    	{
    		printf("====================================================\n");
       		printf("thread: %d || Optimize 3D stucture ---> step: %d\n",thread_num,t);
       		printf("====================================================\n");
       		printf("\n");
    	}

}

float kissing(float x1[1000],float y1[1000],float z1[1000])
{
 
  	float sx1=0.0,sy1=0.0,sz1=0.0,sx2=0.0,sy2=0.0,sz2=0.0,d=0.0,ukiss=0.0;
  	int i;
  	for(i=si1;i<=sj1;i++) {sx1=sx1+x1[i];sy1=sy1+y1[i];sz1=sz1+z1[i];}
  	sx1=sx1/(sj1-si1+1); sy1=sy1/(sj1-si1+1); sz1=sz1/(sj1-si1+1);
  	for(i=si2;i<=sj2;i++) {sx2=sx2+x1[i];sy2=sy2+y1[i];sz2=sz2+z1[i];}
  	sx2=sx2/(sj2-si2+1); sy2=sy2/(sj2-si2+1); sz2=sz2/(sj2-si2+1);
  	d=sqrt((sx1-sx2)*(sx1-sx2)+(sy1-sy2)*(sy1-sy2)+(sz1-sz2)*(sz1-sz2));
  	ukiss=5*(d-15)*(d-15);
  	return ukiss;
}
float force_3(float x1[1000],float y1[1000],float z1[1000])
{
 	float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
 	float vt1x,vt1y,vt1z,vt2x,vt2y,vt2z,vt3x,vt3y,vt3z;
 	float r1,r2,r3,theta1,theta2,theta3;
 	float u_angle=0.0,u_distance=0.0,uf=0;
 	float kissing();
 	int i;

 	c1x=0.0;c1y=0.0;c1z=0.0;
 	for(i=1;i<=xi1;i++) { c1x=c1x+x1[i];c1y=c1y+y1[i];c1z=c1z+z1[i];}
 	for(i=xj1;i<=N0;i++) { c1x=c1x+x1[i];c1y=c1y+y1[i];c1z=c1z+z1[i];}
 	c1x=c1x/length1;  c1y=c1y/length1;  c1z=c1z/length1;    

 	c2x=0.0;c2y=0.0;c2z=0.0;
 	for(i=xi2;i<=xj2;i++) { c2x=c2x+x1[i];c2y=c2y+y1[i];c2z=c2z+z1[i];}
 	c2x=c2x/length2;  c2y=c2y/length2;  c2z=c2z/length2;
 	
 	c3x=0.0;c3y=0.0;c3z=0.0;
 	for(i=xi3;i<=xj3;i++) { c3x=c3x+x1[i];c3y=c3y+y1[i];c3z=c3z+z1[i];}
 	c3x=c3x/length3;  c3y=c3y/length3;c3z=c3z/length3;

 	float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;
 	jx1=(x1[xi1]+x1[xj1]+x1[xi2]+x1[xj2])/4; jy1=(y1[xi1]+y1[xj1]+y1[xi2]+y1[xj2])/4; jz1=(z1[xi1]+z1[xj1]+z1[xi2]+z1[xj2])/4;
 	jx2=(x1[xi2]+x1[xj2]+x1[xi3]+x1[xj3])/4; jy2=(y1[xi2]+y1[xj2]+y1[xi3]+y1[xj3])/4; jz2=(z1[xi2]+z1[xj2]+z1[xi3]+z1[xj3])/4;
 	jx3=(x1[xi1]+x1[xj1]+x1[xi3]+x1[xj3])/4; jy3=(y1[xi1]+y1[xj1]+y1[xi3]+y1[xj3])/4; jz3=(z1[xi1]+z1[xj1]+z1[xi3]+z1[xj3])/4;
/***********************************************************************************/
 	vt1x=c1x-jx1; vt1y=c1y-jy1; vt1z=c1z-jz1;
 	r1=sqrt(vt1x*vt1x+vt1y*vt1y+vt1z*vt1z);

 	vt2x=c2x-jx2; vt2y=c2y-jy2; vt2z=c2z-jz2;
 	r2=sqrt(vt2x*vt2x+vt2y*vt2y+vt2z*vt2z);

 	vt3x=c3x-jx3; vt3y=c3y-jy3; vt3z=c3z-jz3;
 	r3=sqrt(vt3x*vt3x+vt3y*vt3y+vt3z*vt3z);
          
 	theta1=acos((vt1x*vt2x+vt1y*vt2y+vt1z*vt2z)/(r1*r2));  
 	theta2=acos((vt2x*vt3x+vt2y*vt3y+vt2z*vt3z)/(r3*r2));  
 	theta3=acos((vt1x*vt3x+vt1y*vt3y+vt1z*vt3z)/(r1*r3));

     	u_angle=10.0*(theta1-thetaxx1)*(theta1-thetaxx1)+(theta2-thetaxx2)*(theta2-thetaxx2)+(theta3-thetaxx3)*(theta3-thetaxx3);
     	if(thetaxx<2.4)    {u_distance=0.0;}
     	else               {Nstrength=10;u_distance=kissing(x1,y1,z1);}
     	uf=u_angle+u_distance;
 	return uf;
}

float force_4(float x1[1000],float y1[1000],float z1[1000])
{

 	float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,c4x,c4y,c4z;
 	float vt1x,vt1y,vt1z,vt2x,vt2y,vt2z,vt3x,vt3y,vt3z,vt4x,vt4y,vt4z;
 	float r1,r2,r3,r4,theta1,theta2,thetax1,thetax2;
 	float u_angle=0.0,u_distance=0.0,uf=0;
 	float kissing();
 	int i;
/*******************************************************************************/
 	c1x=0.0;c1y=0.0;c1z=0.0;
 	for(i=1;i<=xi1;i++) { c1x=c1x+x1[i];c1y=c1y+y1[i];c1z=c1z+z1[i];}
 	for(i=xj1;i<=N0;i++) { c1x=c1x+x1[i];c1y=c1y+y1[i];c1z=c1z+z1[i];}
 	c1x=c1x/length1;  c1y=c1y/length1;  c1z=c1z/length1;    
 	
 	c2x=0.0;c2y=0.0;c2z=0.0;
 	for(i=xi2;i<=xj2;i++) { c2x=c2x+x1[i];c2y=c2y+y1[i];c2z=c2z+z1[i];}
 	c2x=c2x/length2;  c2y=c2y/length2;  c2z=c2z/length2;

 	c3x=0.0;c3y=0.0;c3z=0.0;
 	for(i=xi3;i<=xj3;i++) { c3x=c3x+x1[i];c3y=c3y+y1[i];c3z=c3z+z1[i];}
 	c3x=c3x/length3;  c3y=c3y/length3;c3z=c3z/length3;

 	c4x=0.0;c4y=0.0;c4z=0.0;
 	for(i=xi4;i<=xj4;i++) { c4x=c4x+x1[i];c4y=c4y+y1[i];c4z=c4z+z1[i];}
 	c4x=c4x/length4;  c4y=c4y/length4;c4z=c4z/length4;

 	float jx1,jy1,jz1,jx2,jy2,jz2;
 	float vt12x,vt12y,vt12z,vt43x,vt43y,vt43z,vt14x,vt14y,vt14z,vt23x,vt23y,vt23z,r12,r43,r14,r23;
/*********************condition 1*****************************************************/
	if(condition==1)
	{
 		jx1=(x1[xi1]+x1[xj1]+x1[xi4]+x1[xj4])/4; jy1=(y1[xi1]+y1[xj1]+y1[xi4]+y1[xj4])/4; jz1=(z1[xi1]+z1[xj1]+z1[xi4]+z1[xj4])/4;
 		jx2=(x1[xi3]+x1[xj3]+x1[xi2]+x1[xj2])/4; jy2=(y1[xi3]+y1[xj3]+y1[xi2]+y1[xj2])/4; jz2=(z1[xi3]+z1[xj3]+z1[xi2]+z1[xj2])/4;
 	
 		vt1x=c1x-jx1; vt1y=c1y-jy1; vt1z=c1z-jz1; r1=sqrt(vt1x*vt1x+vt1y*vt1y+vt1z*vt1z);
 		vt2x=c2x-jx2; vt2y=c2y-jy2; vt2z=c2z-jz2; r2=sqrt(vt2x*vt2x+vt2y*vt2y+vt2z*vt2z);
 		vt3x=c3x-jx2; vt3y=c3y-jy2; vt3z=c3z-jz2; r3=sqrt(vt3x*vt3x+vt3y*vt3y+vt3z*vt3z);
          	vt4x=c4x-jx1; vt4y=c4y-jy1; vt4z=c4z-jz1; r4=sqrt(vt4x*vt4x+vt4y*vt4y+vt4z*vt4z);

 		vt12x=c2x-c1x; vt12y=c2y-c1y; vt12z=c2z-c1z; r12=sqrt(vt12x*vt12x+vt12y*vt12y+vt12z*vt12z);
 		vt43x=c3x-c4x; vt43y=c3y-c4y; vt43z=c3z-c4z; r43=sqrt(vt43x*vt43x+vt43y*vt43y+vt43z*vt43z);
		vt14x=c4x-c1x; vt14y=c4y-c1y; vt14z=c4z-c1z; r14=sqrt(vt14x*vt14x+vt14y*vt14y+vt14z*vt14z);
 		vt23x=c3x-c2x; vt23y=c3y-c2y; vt23z=c3z-c2z; r23=sqrt(vt23x*vt23x+vt23y*vt23y+vt23z*vt23z);
		theta1=acos((vt1x*vt4x+vt1y*vt4y+vt1z*vt4z)/(r1*r4));
 		theta2=acos((vt2x*vt3x+vt2y*vt3y+vt2z*vt3z)/(r3*r2)); 
 	 
 		thetax1=acos((vt12x*vt43x+vt12y*vt43y+vt12z*vt43z)/(r12*r43));
 		thetax2=acos((vt14x*vt23x+vt14y*vt23y+vt14z*vt23z)/(r14*r23));
 	}
 /*****************************condition 2***********************************************************************************/
 	else
 	{
 		jx1=(x1[xi1]+x1[xj1]+x1[xi2]+x1[xj2])/4; jy1=(y1[xi1]+y1[xj1]+y1[xi2]+y1[xj2])/4; jz1=(z1[xi1]+z1[xj1]+z1[xi2]+z1[xj2])/4;
 		jx2=(x1[xi3]+x1[xj3]+x1[xi4]+x1[xj4])/4; jy2=(y1[xi3]+y1[xj3]+y1[xi4]+y1[xj4])/4; jz2=(z1[xi3]+z1[xj3]+z1[xi4]+z1[xj4])/4;
 	
 		vt1x=c1x-jx1; vt1y=c1y-jy1; vt1z=c1z-jz1; r1=sqrt(vt1x*vt1x+vt1y*vt1y+vt1z*vt1z);
 		vt2x=c2x-jx1; vt2y=c2y-jy1; vt2z=c2z-jz1; r2=sqrt(vt2x*vt2x+vt2y*vt2y+vt2z*vt2z);
 		vt3x=c3x-jx2; vt3y=c3y-jy2; vt3z=c3z-jz2; r3=sqrt(vt3x*vt3x+vt3y*vt3y+vt3z*vt3z);
          	vt4x=c4x-jx2; vt4y=c4y-jy2; vt4z=c4z-jz2; r4=sqrt(vt4x*vt4x+vt4y*vt4y+vt4z*vt4z);

 		vt12x=c2x-c1x; vt12y=c2y-c1y; vt12z=c2z-c1z; r12=sqrt(vt12x*vt12x+vt12y*vt12y+vt12z*vt12z);
 		vt43x=c3x-c4x; vt43y=c3y-c4y; vt43z=c3z-c4z; r43=sqrt(vt43x*vt43x+vt43y*vt43y+vt43z*vt43z);
		vt14x=c4x-c1x; vt14y=c4y-c1y; vt14z=c4z-c1z; r14=sqrt(vt14x*vt14x+vt14y*vt14y+vt14z*vt14z);
 		vt23x=c3x-c2x; vt23y=c3y-c2y; vt23z=c3z-c2z; r23=sqrt(vt23x*vt23x+vt23y*vt23y+vt23z*vt23z);

 		theta1=acos((vt1x*vt2x+vt1y*vt2y+vt1z*vt2z)/(r1*r2));
 		theta2=acos((vt3x*vt4x+vt3y*vt4y+vt3z*vt4z)/(r3*r4));  
 	
 		thetax1=acos((vt14x*vt23x+vt14y*vt23y+vt14z*vt23z)/(r14*r23));
 		thetax2=acos((vt12x*vt43x+vt12y*vt43y+vt12z*vt43z)/(r12*r43));
 	}
/*******************************************************************************************************************/	
     	u_angle=10*((theta1-2.4)*(theta1-2.4)+(theta2-2.4)*(theta2-2.4)+(thetax1-thetaxx1)*(thetax1-thetaxx1)+(thetax2-thetaxx)*(thetax2-thetaxx));
     	if(thetaxx>=0.8&&thetaxx<=2.1&&t>50000)
     	{
          	u_distance=kissing(x1,y1,z1);		
     	}
     	else    
     	{
     		u_distance=0.0;
     	}
     	uf=u_angle+u_distance;
 	return uf;
}


float force_5(float x1[1000],float y1[1000],float z1[1000])
{
 	float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,c4x,c4y,c4z,c5x,c5y,c5z;
 	float vt1x,vt1y,vt1z,vt2x,vt2y,vt2z,vt3x,vt3y,vt3z,vt4x,vt4y,vt4z,vt5x,vt5y,vt5z,vt6x,vt6y,vt6z;
 	float r1,r2,r3,r4,r5,r6,theta1,theta2,thetax1,thetax2;
 	float u_angle=0.0,u_distance=0.0,u_angle_4=0.0,uf=0;
 	float kissing();
	int i;
/*******************************************************************************/
	c1x=0.0;c1y=0.0;c1z=0.0;
 	for(i=1;i<=xi1;i++) { c1x=c1x+x1[i];c1y=c1y+y1[i];c1z=c1z+z1[i];}
 	for(i=xj1;i<=N0;i++) { c1x=c1x+x1[i];c1y=c1y+y1[i];c1z=c1z+z1[i];}
 	c1x=c1x/length1;  c1y=c1y/length1;  c1z=c1z/length1;    

 	c2x=0.0;c2y=0.0;c2z=0.0;
 	for(i=xi2;i<=xj2;i++) { c2x=c2x+x1[i];c2y=c2y+y1[i];c2z=c2z+z1[i];}
 	c2x=c2x/length2;  c2y=c2y/length2;  c2z=c2z/length2;

 	c3x=0.0;c3y=0.0;c3z=0.0;
 	for(i=xi3;i<=xj3;i++) { c3x=c3x+x1[i];c3y=c3y+y1[i];c3z=c3z+z1[i];}
 	c3x=c3x/length3;  c3y=c3y/length3;c3z=c3z/length3;

 	c4x=0.0;c4y=0.0;c4z=0.0;
 	for(i=xi4;i<=xj4;i++) { c4x=c4x+x1[i];c4y=c4y+y1[i];c4z=c4z+z1[i];}
 	c4x=c4x/length4;  c4y=c4y/length4;c4z=c4z/length4;

 	c5x=0.0;c5y=0.0;c5z=0.0;
 	for(i=xi5;i<=xj5;i++) { c5x=c5x+x1[i];c5y=c5y+y1[i];c5z=c5z+z1[i];}
 	c5x=c5x/length5;  c5y=c5y/length5;c5z=c5z/length5;

 	float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;
 	
 	if(condition==1)
 	{
 		float vt12x,vt12y,vt12z,vt53x,vt53y,vt53z,vt15x,vt15y,vt15z,vt23x,vt23y,vt23z,r12,r53,r15,r23;
 		jx1=(x1[xi1]+x1[xj1]+x1[xi5]+x1[xj5])/4; jy1=(y1[xi1]+y1[xj1]+y1[xi5]+y1[xj5])/4; jz1=(z1[xi1]+z1[xj1]+z1[xi5]+z1[xj5])/4;
 		jx2=(x1[xi2]+x1[xj2]+x1[xi3]+x1[xj3])/4; jy2=(y1[xi2]+y1[xj2]+y1[xi3]+y1[xj3])/4; jz2=(z1[xi2]+z1[xj2]+z1[xi3]+z1[xj3])/4;
 		jx3=(x1[xi4]+x1[xj4]+x1[xi5]+x1[xj5])/4; jy3=(y1[xi4]+y1[xj4]+y1[xi5]+y1[xj5])/4; jz3=(z1[xi4]+z1[xj4]+z1[xi5]+z1[xj5])/4;
 
 		vt1x=c1x-jx1; vt1y=c1y-jy1; vt1z=c1z-jz1;r1=sqrt(vt1x*vt1x+vt1y*vt1y+vt1z*vt1z);
		vt4x=c5x-jx1; vt4y=c5y-jy1; vt4z=c5z-jz1;r4=sqrt(vt4x*vt4x+vt4y*vt4y+vt4z*vt4z);
 		vt2x=c2x-jx2; vt2y=c2y-jy2; vt2z=c2z-jz2;r2=sqrt(vt2x*vt2x+vt2y*vt2y+vt2z*vt2z);
		vt3x=c3x-jx2; vt3y=c3y-jy2; vt3z=c3z-jz2;r3=sqrt(vt3x*vt3x+vt3y*vt3y+vt3z*vt3z);
    
 		vt5x=c5x-jx3; vt5y=c5y-jy3; vt5z=c5z-jz3;r5=sqrt(vt5x*vt5x+vt5y*vt5y+vt5z*vt5z);
		vt6x=c4x-jx3; vt6y=c4y-jy3; vt6z=c4z-jz3;r6=sqrt(vt6x*vt6x+vt6y*vt6y+vt6z*vt6z);

 		vt12x=c2x-c1x; vt12y=c2y-c1y; vt12z=c2z-c1z;r12=sqrt(vt12x*vt12x+vt12y*vt12y+vt12z*vt12z);
 		vt53x=c3x-c5x; vt53y=c3y-c5y; vt53z=c3z-c5z;r53=sqrt(vt53x*vt53x+vt53y*vt53y+vt53z*vt53z);
		vt15x=c5x-c1x; vt15y=c5y-c1y; vt15z=c5z-c1z;r15=sqrt(vt15x*vt15x+vt15y*vt15y+vt15z*vt15z);
 		vt23x=c3x-c2x; vt23y=c3y-c2y; vt23z=c3z-c2z;r23=sqrt(vt23x*vt23x+vt23y*vt23y+vt23z*vt23z);
 		
 		theta1=acos((vt2x*vt3x+vt2y*vt3y+vt2z*vt3z)/(r3*r2));  
 		theta2=acos((vt1x*vt4x+vt1y*vt4y+vt1z*vt4z)/(r1*r4));
 		thetax1=acos((vt12x*vt53x+vt12y*vt53y+vt12z*vt53z)/(r12*r53));
 		thetax2=acos((vt15x*vt23x+vt15y*vt23y+vt15z*vt23z)/(r15*r23));
 	}
 	else
	{
		float vt12x,vt12y,vt12z,vt54x,vt54y,vt54z,vt15x,vt15y,vt15z,vt24x,vt24y,vt24z,r12,r54,r15,r24;
 		jx1=(x1[xi1]+x1[xj1]+x1[xi2]+x1[xj2])/4; jy1=(y1[xi1]+y1[xj1]+y1[xi2]+y1[xj2])/4; jz1=(z1[xi1]+z1[xj1]+z1[xi2]+z1[xj2])/4;
 		jx2=(x1[xi5]+x1[xj5]+x1[xi4]+x1[xj4])/4; jy2=(y1[xi5]+y1[xj5]+y1[xi4]+y1[xj4])/4; jz2=(z1[xi5]+z1[xj5]+z1[xi4]+z1[xj4])/4;
 		jx3=(x1[xi4]+x1[xj4]+x1[xi3]+x1[xj3])/4; jy3=(y1[xi4]+y1[xj4]+y1[xi3]+y1[xj3])/4; jz3=(z1[xi4]+z1[xj4]+z1[xi3]+z1[xj3])/4;
 
 		vt1x=c1x-jx1; vt1y=c1y-jy1; vt1z=c1z-jz1;r1=sqrt(vt1x*vt1x+vt1y*vt1y+vt1z*vt1z);
		vt2x=c2x-jx1; vt2y=c2y-jy1; vt2z=c2z-jz1;r2=sqrt(vt2x*vt2x+vt2y*vt2y+vt2z*vt2z);
 		vt3x=c3x-jx2; vt3y=c3y-jy2; vt3z=c3z-jz2;r3=sqrt(vt3x*vt3x+vt3y*vt3y+vt3z*vt3z);
		vt4x=c4x-jx2; vt4y=c4y-jy2; vt4z=c4z-jz2;r4=sqrt(vt4x*vt4x+vt4y*vt4y+vt4z*vt4z);
    
 		vt5x=c3x-jx3; vt5y=c3y-jy3; vt5z=c3z-jz3;r5=sqrt(vt5x*vt5x+vt5y*vt5y+vt5z*vt5z);
		vt6x=c4x-jx3; vt6y=c4y-jy3; vt6z=c4z-jz3;r6=sqrt(vt6x*vt6x+vt6y*vt6y+vt6z*vt6z);

 		vt12x=c2x-c1x; vt12y=c2y-c1y; vt12z=c2z-c1z;r12=sqrt(vt12x*vt12x+vt12y*vt12y+vt12z*vt12z);
 		vt54x=c4x-c5x; vt54y=c4y-c5y; vt54z=c4z-c5z;r54=sqrt(vt54x*vt54x+vt54y*vt54y+vt54z*vt54z);
		vt15x=c5x-c1x; vt15y=c5y-c1y; vt15z=c5z-c1z;r15=sqrt(vt15x*vt15x+vt15y*vt15y+vt15z*vt15z);
 		vt24x=c4x-c2x; vt24y=c4y-c2y; vt24z=c4z-c2z;r24=sqrt(vt24x*vt24x+vt24y*vt24y+vt24z*vt24z);
 		
 		theta1=acos((vt1x*vt2x+vt1y*vt2y+vt1z*vt2z)/(r1*r2));  
 		theta2=acos((vt3x*vt4x+vt3y*vt4y+vt3z*vt4z)/(r3*r4));
 		thetax1=acos((vt12x*vt54x+vt12y*vt54y+vt12z*vt54z)/(r12*r54));
 		thetax2=acos((vt15x*vt24x+vt15y*vt24y+vt15z*vt24z)/(r15*r24));
 	}
 	
 	float thetax5;
 	thetax5=acos((vt5x*vt6x+vt5y*vt6y+vt5z*vt6z)/(r5*r6));
 	
     	if(fmod(t,100000)==0&&cycle==1)
     	{
     		thetaxx4=thetaxx4+0.15;
     		if(thetaxx4>thetaxx)
     		{
     			thetaxx4=0.0;
     		}
     	}
     	u_angle=5*((theta1-2.4)*(theta1-2.4)+(theta2-2.4)*(theta2-2.4)+(thetax1-thetaxx1)*(thetax1-thetaxx1)+(thetax2-thetaxx)*(thetax2-thetaxx));
     	u_angle_4=5*(thetax5-thetaxx4)*(thetax5-thetaxx4);
     	if(thetaxx>=0.6&&thetaxx<=2.2) 
     	{
	 	u_distance=kissing(x1,y1,z1);			
     	}
     	else	
     	{
         	u_distance=0.0;
     	}
     	uf=u_angle+u_distance+u_angle_4;
 	return uf;
}
