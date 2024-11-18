/*********** 本程序用来预测RNA三维结构以及热稳定性 ，使用副本交换蒙特卡洛退火算法 *********************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define pi 3.1415
#define step1 .5           //Translation step length in Pivot move in Folding
#define step1_1 .2         //Translation step length in 3nt fragment move in Folding
#define step2 1.0          //One atom move step in Optimization
#define step3 .06          //All atoms move step in Optimization
#define OPTIMIZE_t 1000       //The frequency of all atoms move in the 1+2 Optimization way
#define total0 100000   //Optimization setps
#define total2 5000000   //Max steps in Folding at constant T (melting)         
#define Alpha 0.992        //Annealing rate
#define Anneal 1           //Annealing shecdle
#define A -3.5            //base-pairing strength of GC
#define B01 -9.6           //conformational entropy of base-stacking in CoaxialStacking potential
#define lt 25             //The lowest temperature in annealing folding process
//#define ht 35
#define jt 10              //Cooling with constant T of jt (e.g., jt=5C) 
#define tran 10            //rotational frequency
#define BetaGC 0.93
#define Beta 0.6          //AU=Beta*GC
#define Beta1 0.0          //GA(UU)=Beta1*GC
#define Beta2 0.0          //GG=Beta2*GC
#define Beta3 0.0          //??=Beta3*GC  For non-canoncial base-pairing
#define gamma 0.1         //The energy threshold of base-pairing for statistics (碱基配对的能量Ubp0<gamma*Ubp算作配对形成)
#define tconf1 500       //The frequency of conformational output
#define tconf2 10000        //The frequency of conformational output at 25C in structure prediction
//#define trmsd 5000       //类似MD，计算相邻trmsd的构象间的差别
#define tenergy 500     //The frequency of energy calculation
#define tbp 100          //The frequency of output information for base-pairing
#define tG 0            //计算平均base-pair数量的起始
#define tG_s 100           //Calculate pbp   经简单测试貌似还是100靠谱，设置为10000并没有节省计算时间也没能更好的呈现平衡
#define fp_bp 500         //The frequency of output the value of BP & PBP & G
#define tprint 10000       //The frequency of screen output
#define t_OPT1 100        //The frequency of conformation output in Optimization
#define t_OPT2 1000
#define t_si 1000          //Output frequency of structure information
#define tMove 1           //The frequency of 3-nt fragment moves；
// Electrostatic
//#define Ek 78.0          // debye parameters
#define bl 5.45
//#define RG 1.1           //TBI,the Rg of RNAs comparing with Rg of A-form helix;
#define IMg 3.0            //2:2 IMg=4.0; 2:1 IMg=3.0;
#define thread 15
#define ww 0.5

float T,Kd,I,Ek,q4,CNa,CMg,fNa,D,B0;
int N,N0,N1,total,Nbp,t,l,ll,lll,llll,salt,nm,Folding=1,Energy,Dangling,Bulge,Internal,Pseudoknot,Move=0;
int s[1000][1000],ss[1000][1000],a[10000],c[10000],bp,bp0,BP0,BP,am[1000],cm[1000],sm[1000][1000],ssm[1000][1000];
int OPTIMIZE_W,OPTIMIZE,M_Y,Stem[1000],Internal_1[1000],Internal_b[1000],Internal_a[1000],Bulge_1[1000],Pseudoknot_1[1000],t_OPT;
float rand01,phi,theta,rm,d1,d2,step;
float dis0,Xc,Yc,Zc,Rg,end,pl2,b,pl3,pl4,rgp,rep,pl3p,pl4p,pl0,pl0p,pl2p,pl2,PBP,pbp,G,GG,m,M;
float x[10000],y[10000],z[10000],xx[10000],yy[10000],zz[10000],xxx[10000],yyy[10000],zzz[10000],Q[10000],R[10000],f[10000];
float U,Umin,Up,Up0,du,uu,u,ulj,uulj,uc,uuc,ub,uub,ue,uue,ud,uud,uN,uuN,us,uus,uco,uuco,uTri,uuTri,umis,uumis,udang,uudang;
float uc0,uuc0,uuc1,uc1,ulj0,uulj0,ulj1,uulj1,uN1,uuN1,utrs,uutrs;
char type[10000];
/********The parameters (including bonded & nonbonded) of the CG Model************/
float e=0.26;    //Strength in Uexc(LJ potential) 
float kpc=60,lpc=3.92,kcp=33,lcp=3.9,kcN=24.8,lcN=3.346;                                      //Bond
float kpcp=6.3,Apcp=1.82,kcpc=15.3,Acpc=1.8,kpcN=9.75,ApcN=1.64,kNcp=15.23,ANcp=1.66;        //Angle
float kpcpc=0.8,dpcpc=2.56,kcpcp=3.6,dcpcp=-2.94,kcpcN=0.825,dcpcN=-1.164,kNcpc=0.76,dNcpc=0.88; //Dihedral
float dNN1=8.55,dNN2=9.39,sigma=4.786,dco=5.0,dkco=1.0,kdco=2.5; //base-pairing & base-stacking; dkco:共轴堆积能量与碱基对堆积能量的差异；
float qq4,frac[1000],fi[1000],frac0[1000],frac1[1000],frac2[1000];
float tt0[15]={25.0,31.0,37.0,43.0,49.0,55.0,63.0,71.0,78.0,86.0,94.0,102.0,110.0,120.0,130.0},t0;
int tem_j;
float Bs[24];
//dHmis/dSmis:avg. value of terminal mismatch parameters;Kmis:First Mismatch base pairs (UU/GA, or GG); GUend:GU or AU end penalty
FILE *fp,*fp_conf[thread],*fp_sec_stru[thread],*fp_Energy[thread],*fp8,*fp9,*fp_Bp[thread],*fpcfig;
//FILE *fp4;
int main()
{
	int i,duo1,duo2,result;
	float x0[10000],y0[10000],z0[10000]/*,x1[10000],y1[10000],z1[10000]*/;
	void Put_File(),Fixed_Atom(),MC_Annealing(),OutputPara();

 	Put_File();         //defined the input parameters and out file names;
 	OutputPara();       //the output of the parametes;

 	i=1;
 	while(!feof(fp)) 
 	{
          	result=fscanf(fp,"%d %d %s %f %f %f %f %f %f\n",&duo1,&duo2,&type[i],&x0[i],&y0[i],&z0[i],&R[i],&Q[i],&f[i]); 
          	i++;
 	}
 	fclose(fp);        //input of the initial conformation 
 	N0=i-1;          //N0:Total number of CG beads; N: Num. of nt
 	N=(N0-1)/3;
 	Fixed_Atom();      //The Centre atom N1 will be fixed;
 	MC_Annealing(x0,y0,z0);

	for(i=0;i<thread;i++)
	{
		fclose(fp_conf[i]);fclose(fp_sec_stru[i]); fclose(fp_Energy[i]);fclose(fp_Bp[i]);
	}
	(void)result;
	fclose(fp9);
 	// fclose(fp4);fclose(fp8);fclose(fp_conf0);
 	return 0;
}  


/* &%$#@!~&%$#@!~&%$#@!~   Some functions or modules for move and calculation    &%$#@!~&%$#@!~&%$#@!~ */
/*********************读入文件，读出文件，输入参数************************/
void Put_File (void) 
{
	int i,result;
	char filename[30];
 	fp=fopen("ch_0.dat","r+");                //initial conformation
 	for (i=0;i<thread;i++) 
     	{
           	sprintf(filename,"conf_%d.dat",i); fp_conf[i]=fopen(filename,"w+");
           	sprintf(filename,"sec_stru_%d.dat",i);fp_sec_stru[i]=fopen(filename,"w+");         //Details of base-pairs 
 		sprintf(filename,"Energy_%d.dat",i);fp_Energy[i]=fopen(filename,"w+"); 
 		sprintf(filename,"Bp_%d.dat",i);fp_Bp[i]=fopen(filename,"w+");    
 	
           	//sprintf(filename,"Energy_%d.dat",i);    fpenergy[i]=fopen(filename,"w+");
      	}                //output of conformations
 	         //The Energy of the conformations 
 	fp9=fopen("para.dat","w+");             //Just like the log file for outputing the important paramaters. 
 	fpcfig=fopen("config1.dat","r+");
 	int ca,cb,cc,cd,ce,ct;
 	while(!feof(fpcfig))
 	{
 		result=fscanf(fpcfig,"%d %d %d %d %d %d\n",&ct,&ca,&cb,&cc,&cd,&ce);
 	}
 	fclose(fpcfig);
 	(void)result;
 	total=ca;
 	CNa=cc;
 	CMg=cd;
 	if(CNa==0&&CMg==0) 	{salt=0;}
 	else			{salt=1;}
 //fp4=fopen("energy.dat","w+");           //用于控制能量的输出以调试程序；
 //fp8=fopen("rmsd-t.dat","w+"); //fp_conf0=fopen("Debye-value.dat","w+"); 
}
/*********************************************/
void Fixed_Atom(void)
{
 	int N10;
 	N10=floor(N0/2)+1; 
 	if (fmod(N10,3)==0) {N1=N10-2;}
 	else if (fmod(N10+1,3)==0) {N1=N10-1;}
 	else {N1=N10;}                         //N1:Fixed the P in centre of chain
 	//printf("N %d N0 %d N1 %d\n",N,N0,N1);  //The number of nucleotides & atoms & the actionless atom
}
//******************************************//
void Parameters_T(float t0)
{
   	float tt,TT;
   	void Bs_stacking();
   	//if (Folding==0)      {total=total0;}        //Steps in Optimization
   	//else                 {total=total2;}        //Steps in Folding at constant temperature
   	tt=0.0+t0; 
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
    	//printf("Temp: %f\nI: %ffNa: %f\n",tt,I,fNa); 
    	if(N<=13) 		
     	{
     		if(t0<=55)	{B0=-9.3;}
     		else 		{B0=-11.3;}
     	}
     	else if(N>13&&N<20)     { B0=-10.6;}
     	else 			{ B0=-12.0;}    
     	Bs_stacking();      
    	fprintf(fp9,"Temp %f T %f Kd %f I %f Ek %f q4 %f\n",tt,T,Kd,I,Ek,q4); 
    	fflush(fp9);  //The output of important para.
}
/***********************************/
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
void MC_Annealing(float x0[1000],float y0[1000],float z0[1000])
{
    	int k,i;
    	void MC_T();
    	l=0;ll=0;lll=0;llll=0;                 //MC steps independent of T
    	for (k=thread-1;k>=0;k--)                  //Temperature cycle mechanism
   	{
   		t0=tt0[k];
   		tem_j=k;
            	Parameters_T(t0);     //parameters at any t0: steps, T, Debye length, ionic strength, ion fraction etc.
            	for (i=1;i<=N0;i++) 
            	{
                   	x[i]=x0[i];y[i]=y0[i];z[i]=z0[i];
            	}  //The initconf.=last one at previous T 
            	MC_T();
            	for(i=1;i<=N0;i++) {x0[i]=x[i];y0[i]=y[i];z0[i]=z[i];}
     	}     //End of the T cycle;
}
/*****************Monte Carlo simulation at given Temperature******************************/
void MC_T(void)
{
 	int i,ii;
 	void MC_Each_Step(),StructureInformation(),ENERGY();
 	srand((unsigned)time(NULL));           //Random changes over time
 	void disbrute_qq();
 	for (t=1;t<=total;t++)
 	{    
 		if(salt==1)
		{
           		if(fmod(t,50)==0||t==1)
           		{
                        	disbrute_qq();         
             		}
		} 
             	if (fmod((t-1),100)==0)    
             	{
                      	nm=0;
             	}  //nm: the No. of acceptance each 100 steps;
 
              	for(ii=1;ii<=100;ii++)
              	{    
                    	MC_Each_Step();
              	}                               //End of nucleotide cycle within each t step

            	if (fmod(t,tconf1)==0) 
           	{
              		for (i=1;i<=N0;i++) 
              		{
                         	fprintf(fp_conf[tem_j],"%d %d %c %f %f %f %f %f %f\n",t,i,type[i],x[i],y[i],z[i],R[i],Q[i],f[i]);
                   	}
          	}
          	if(fmod(t,1000)==0)
          	{
          		//printf("MCSA || Folding RNA 3D structure from the sequence\n");
                   	printf("Tem: %.1f || step: %d\n",t0,t);
                   	printf("\n");
          	}
		fflush(fp_conf[tem_j]);
		ENERGY();  
		if (fmod(t,fp_bp)==0) 
      		{
      			fprintf(fp_Bp[tem_j],"%d %d %d\n",t,bp,BP); 
     		}   
      		fflush(fp_Bp[tem_j]);   
 	}           //Calculate energy of one conformation                                  
}
/*****************Monte Carlo for each step******************************/
void MC_Each_Step()
{
 	int i0;
 	float PC(),CP(),CN(),PCP(),CPC(),PCN(),NCP(),PCPC(),CPCP(),PCPCh(),CPCPh(),CPCN(),NCPC(),LJ0(),HB(),St(),QQ();
 	void MoveN(),Rand01(),Translate(),Pivot(),RMSDconf(),FOLD(),Optimize();
 	void BasePairing(),BaseStacking(),ExcludedVolume(),Electrostatic(),Bonded(),CoaxialStacking(),Triplex(),TerminalMismatch(),DanglingEnd(),constraint();
 	void Metropolis(); 

   	u=0.0;uu=0.0;                                    //Initialize the energies;
   	ulj=0.0;ub=0.0;ue=0.0;ud=0.0;uN=0.0;us=0.0;uc=0.0;
   	uulj=0.0;uub=0.0;uue=0.0;uud=0.0;uuN=0.0;uus=0.0;uuc=0.0;
   	uco=0.0;uuco=0.0;uTri=0.0;uuTri=0.0;umis=0.0;uumis=0.0;udang=0.0;uudang=0.0;utrs=0.0;uutrs=0.0;
   	Energy=0;                      //if Energy=1, fuction ENERGY() is running 
   	Move=0;      //if Move_3nt=1, 3nt fragment moves;

       /*****************固定中心原子*********************/
   	re:;
   	rand01=rand()/(RAND_MAX+1.); i0=floor(rand01*N0)+1; if (i0==N1) goto re;       
       /***********如果Folding=0，则执行优化，否则执行折叠***********/
   	if (Folding==0)                     //structure optimization
   	{
         	if (OPTIMIZE_W==0||OPTIMIZE_W==1||OPTIMIZE_W==2) 
         	{
                    	OPTIMIZE=OPTIMIZE_W;
         	}
         	else 
         	{
                  	if (fmod(t,OPTIMIZE_t)==0) 
                  	{
                        	OPTIMIZE=1;
                  	} 
                  	else 
                  	{
                        	OPTIMIZE=2;
                  	}
         	}
         	Optimize(i0);                      //The ways of conformational changes should be improved--2018.01.01;
    	}
    	else    
    	{
             	FOLD(i0);
    	}                 //Folding Process
    	Metropolis();
 }
/*********************************************/
void Metropolis(void)
{
   	int i;
   	float p;
   	u=ulj+uN+us+uc+(ub+ue+ud)*0.5963+uco/*+uTri+umis+udang*/;
   	uu=uulj+uuN+uus+uuc+(uub+uue+uud)*0.5963+uuco/*+uuTri+uumis+uudang*/;
   	du=uu-u; 
   	p=0.0;  //du: Energy changes before & after moves;

   	if(du<=0.0)   
   	{
          	for (i=1;i<=N0;i++) 
          	{
                   	x[i]=xx[i];y[i]=yy[i];z[i]=zz[i];
          	}   
          	nm=nm+1; 
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
                    	nm=nm+1;  
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
 	if (fmod(t,2)==0)  //Euler ratation
 	{
 		xx[i1]=(xxx[i1]-xxx[i10])*(cos(theta)*cos(rm)-cos(phi)*sin(theta)*sin(rm))+(yyy[i1]-yyy[i10])*(sin(theta)*cos(rm)+cos(phi)*cos(theta)*sin(rm))+(zzz[i1]-zzz[i10])*sin(phi)*sin(rm)+xxx[i10];
 		yy[i1]=-(xxx[i1]-xxx[i10])*(cos(theta)*sin(rm)+cos(phi)*sin(theta)*cos(rm))+(yyy[i1]-yyy[i10])*(cos(phi)*cos(theta)*cos(rm)-sin(theta)*sin(rm))+(zzz[i1]-zzz[i10])*sin(phi)*cos(rm)+yyy[i10];
 		zz[i1]=(xxx[i1]-xxx[i10])*sin(phi)*sin(theta)-sin(phi)*cos(theta)*(yyy[i1]-yyy[i10])+(zzz[i1]-zzz[i10])*cos(phi)+zzz[i10];
 	}
 	else              //Z axis rotation 
 	{
 		xx[i1]=(xxx[i1]-xxx[i10])*cos(rm)-(yyy[i1]-yyy[i10])*sin(rm)+xxx[i10];
 		yy[i1]=(xxx[i1]-xxx[i10])*sin(rm)+(yyy[i1]-yyy[i10])*cos(rm)+yyy[i10];
 		zz[i1]=zzz[i1];
  	}
}
void MoveN(int i1) //Movement of each base;
{
  	Rand01();
  	xxx[i1]=x[i1]+step*rand01*sin(phi)*cos(theta); yyy[i1]=y[i1]+step*rand01*sin(phi)*sin(theta);  zzz[i1]=z[i1]+step*rand01*cos(phi);
  	Rand01();
  	Pivot(i1,i1-1);
}

/* ~~~~~~~~~~ Details of bonded potential calculation ~~~~~~~~~~~~~~ */
float PC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float d,ul;
 	d=sqrt((x1[i1-2]-x1[i1-1])*(x1[i1-2]-x1[i1-1])+(y1[i1-2]-y1[i1-1])*(y1[i1-2]-y1[i1-1])+(z1[i1-2]-z1[i1-1])*(z1[i1-2]-z1[i1-1]));
	ul=kpc*(d-lpc)*(d-lpc);
 	return ul;
}
float CP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float d,ul;
 	d=sqrt((x1[i1+1]-x1[i1-1])*(x1[i1+1]-x1[i1-1])+(y1[i1+1]-y1[i1-1])*(y1[i1+1]-y1[i1-1])+(z1[i1+1]-z1[i1-1])*(z1[i1+1]-z1[i1-1]));
	ul=kcp*(d-lcp)*(d-lcp);
 	return ul;
}
float CN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float d,ul;
 	d=sqrt((x1[i1]-x1[i1-1])*(x1[i1]-x1[i1-1])+(y1[i1]-y1[i1-1])*(y1[i1]-y1[i1-1])+(z1[i1]-z1[i1-1])*(z1[i1]-z1[i1-1]));
	ul=kcN*(d-lcN)*(d-lcN);
 	return ul;
}
float PCP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float d1,d2,d3,w,a1,ue0;
 	d1=sqrt((x1[i1-2]-x1[i1-1])*(x1[i1-2]-x1[i1-1])+(y1[i1-2]-y1[i1-1])*(y1[i1-2]-y1[i1-1])+(z1[i1-2]-z1[i1-1])*(z1[i1-2]-z1[i1-1]));
 	d2=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
 	d3=sqrt((x1[i1+1]-x1[i1-2])*(x1[i1+1]-x1[i1-2])+(y1[i1+1]-y1[i1-2])*(y1[i1+1]-y1[i1-2])+(z1[i1+1]-z1[i1-2])*(z1[i1+1]-z1[i1-2]));
 	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
 	if (w<=-1.0) {a1=3.14;}
 	else if (w>=1.0) {a1=0.;}
 	else  {a1=acos(w);}
 	ue0=kpcp*(a1-Apcp)*(a1-Apcp);   //由nab产生的标准AformRNA中的峰值为1.63，以下的cpc，pcN,Ncp分别为：1.8，1.63，1.68
 	return ue0;
}
float CPC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float d1,d2,d3,w,a1,ue0;
 	d1=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
 	d2=sqrt((x1[i1+1]-x1[i1+2])*(x1[i1+1]-x1[i1+2])+(y1[i1+1]-y1[i1+2])*(y1[i1+1]-y1[i1+2])+(z1[i1+1]-z1[i1+2])*(z1[i1+1]-z1[i1+2]));
 	d3=sqrt((x1[i1+2]-x1[i1-1])*(x1[i1+2]-x1[i1-1])+(y1[i1+2]-y1[i1-1])*(y1[i1+2]-y1[i1-1])+(z1[i1+2]-z1[i1-1])*(z1[i1+2]-z1[i1-1]));
 	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
 	if (w<=-1.0) {a1=3.14;}
 	else if (w>=1.0) {a1=0.;}
 	else  {a1=acos(w);}
	ue0=kcpc*(a1-Acpc)*(a1-Acpc);
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
	ue0=kpcN*(a1-ApcN)*(a1-ApcN);
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
	ue0=kNcp*(a1-ANcp)*(a1-ANcp);
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
 	e1=sqrt(c1*c1+c2*c2+c3*c3); f1=sqrt(p1*p1+p2*p2+p3*p3);
 	pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 	g1=(x1[i1-2]-x1[i1+2]); g2=(y1[i1-2]-y1[i1+2]); g3=(z1[i1-2]-z1[i1+2]);
 	gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 	if (pp1<=-1.0) {di=-3.14;}
 	else if (pp1>=1.0) {di=0.;}
 	else if (hh1>=0.) {di=acos(pp1);}
 	else {di=-acos(pp1);}
//ud0=1.88*(1-cos(di-2.6))+0.5*1.88*(1-cos(3.*(di-2.6))); //由nab产生的AformRNA中的峰值为2.6，cpcp，cpcN,Ncpc分别为：-2.97，-1.286，0.97
	ud0=0.81*(1+cos(di+0.48))-0.10*(1+cos(2*di+2.71))+0.26*(1+cos(3*di+1.38))+0.17*(1+cos(4*di+0.56));
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
 	e1=sqrt(c1*c1+c2*c2+c3*c3); f1=sqrt(p1*p1+p2*p2+p3*p3);
 	pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 	g1=(x1[i1-1]-x1[i1+4]); g2=(y1[i1-1]-y1[i1+4]); g3=(z1[i1-1]-z1[i1+4]);
 	gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 	if (pp1<=-1.0) {di=-3.14;}
 	else if (pp1>=1.0) {di=0.;}
 	else if (hh1>=0.) {di=acos(pp1);}
 	else {di=-acos(pp1);}
	ud0=(0.45*(1+cos(di-0.40))+0.40*(1+cos(2*di+2.68))+0.36*(1+cos(3*di-0.84))+0.26*(1+cos(4*di+2.5)))*1.2;
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
 	e1=sqrt(c1*c1+c2*c2+c3*c3); f1=sqrt(p1*p1+p2*p2+p3*p3);
 	pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 	g1=(x1[i1-4]-x1[i1]); g2=(y1[i1-4]-y1[i1]); g3=(z1[i1-4]-z1[i1]);
 	gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 	if (pp1<=-1.0) {di=-3.14;}
 	else if (pp1>=1.0) {di=0.;}
 	else if (hh1>=0.) {di=acos(pp1);}
 	else {di=-acos(pp1);}
	  ud0=(0.47*(1+cos(di-2.04))+0.33*(1+cos(2*di-0.96))+0.25*(1+cos(3*di+0.07))+0.19*(1+cos(4*di+1.98)));
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
 	g1=(x1[i1]-x1[i1+2]); g2=(y1[i1]-y1[i1+2]); g3=(z1[i1]-z1[i1+2]);
 	gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 	if (pp1<=-1.0) {di=-3.14;}
 	else if (pp1>=1.0) {di=0.;}
 	else if (hh1>=0.) {di=acos(pp1);}
 	else {di=-acos(pp1);}
	ud0=0.65*(1+cos(di+2.27))+0.20*(1+cos(2*di+1.57))+0.20*(1+cos(3*di+1.15))+0.06*(1+cos(4*di+0.18));
 	return ud0;
}
/*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
//base pairing between two complementary bases (AU,GC,and GU)
float HB(int i1,int j1,float x1[10000],float y1[10000],float z1[10000])
{
 	float d,hb,UHB,d0,d01,d1,d11;
 	d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); //NN
 	if (d>=dNN1&&d<=dNN2)   //The pairing formation condition
 	{
 		d0=sqrt((x1[i1]-x1[j1-1])*(x1[i1]-x1[j1-1])+(y1[i1]-y1[j1-1])*(y1[i1]-y1[j1-1])+(z1[i1]-z1[j1-1])*(z1[i1]-z1[j1-1]));  //NiCj
 		d01=sqrt((x1[j1]-x1[i1-1])*(x1[j1]-x1[i1-1])+(y1[j1]-y1[i1-1])*(y1[j1]-y1[i1-1])+(z1[j1]-z1[i1-1])*(z1[j1]-z1[i1-1])); //CiNj
 		d1=sqrt((x1[i1-2]-x1[j1])*(x1[i1-2]-x1[j1])+(y1[i1-2]-y1[j1])*(y1[i1-2]-y1[j1])+(z1[i1-2]-z1[j1])*(z1[i1-2]-z1[j1]));  //PiNj
 		d11=sqrt((x1[j1-2]-x1[i1])*(x1[j1-2]-x1[i1])+(y1[j1-2]-y1[i1])*(y1[j1-2]-y1[i1])+(z1[j1-2]-z1[i1])*(z1[j1-2]-z1[i1])); //NiPj
     		if ((type[i1]=='G'&&type[j1]=='C')||(type[i1]=='C'&&type[j1]=='G'))  {hb=BetaGC*A;}
		else if ((type[i1]=='G'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='G'))  {hb=Beta*A;}
		else if ((type[i1]=='A'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='A'))  {hb=Beta*A;} 
		else if ((type[i1]=='A'&&type[j1]=='G')||(type[i1]=='G'&&type[j1]=='A')||(type[i1]=='U'&&type[j1]=='U'))  {hb=Beta1*A;}  //mismatches
		else if (type[i1]=='G'&&type[j1]=='G')  {hb=Beta2*A;}
		else {hb=Beta3*A;}   //Beta1,2,3=0,不考虑mismatches
//UHB=hb/(1+3.6*(d-8.94)*(d-8.94)+1.9*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.7*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));  
  		UHB=hb/(1+2.66*(d-8.94)*(d-8.94)+1.37*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.464*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));  
  	}
 	else UHB=0.0;
 	return UHB;
}
void BasePairing()
{
    	int i,j;
   	bp=0;BP=0;
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
                                		uN1=0.0;uuN1=0.0; bp0=0;BP0=0;
                                		if (c[i]==0&&c[j]==0) 
                                		{
                                    			uN1=HB(i,j,x,y,z);
                                    			if (uN1!=0) 
                                    			{
                                         			c[i]=1;c[j]=1; s[i][j]=1; bp0=1; 
                                         			if (uN1<=-0.66) 
                                		 		{
                                      			 		BP0=1;
                                		 		}               
                                     			}
                                 		}
                                 		uN=uN+uN1;  bp=bp+bp0;  BP=BP+BP0;    
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
      // if (cm[i]==0&&cm[j]==0) {uN1=HB(i,j,x,y,z); if (uN1!=0) {cm[i]=1;cm[j]=1;sm[i][j]=1;}}
       //if (am[i]==0&&am[j]==0) {uuN1=HB(i,j,xx,yy,zz); if (uuN1!=0) {am[i]=1;am[j]=1;ssm[i][j]=1;}}  //Mismatched stacking 不考虑
                                  		uN=uN+uN1; uuN=uuN+uuN1;
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
	else if ((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='G'))   	
	{
        	if 	((type[i1-3]=='G'&&type[k1+3]=='C')&&(type[j1+3]=='C'&&type[k2-3]=='G'))	{B=Bs[21];}
        	else                                                                       		{B=Bs[22];}
	}// GU/UG (CGUC/CUGG) stacking
   	return B;
}

float St(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
  	float d,d0,kst,ULJ,d1,d5,d6,d10,d12,d01,d05,d06,d010,d012/*,ulj1,ulj2*/;
  	d=sqrt((x1[i1]-x1[i1+3])*(x1[i1]-x1[i1+3])+(y1[i1]-y1[i1+3])*(y1[i1]-y1[i1+3])+(z1[i1]-z1[i1+3])*(z1[i1]-z1[i1+3])); 
  	d1=sigma/d; d5=d1*d1*d1*d1*d1;d6=d1*d5;d10=d5*d5;d12=d6*d6; 
  	d0=sqrt((x1[j1]-x1[j1-3])*(x1[j1]-x1[j1-3])+(y1[j1]-y1[j1-3])*(y1[j1]-y1[j1-3])+(z1[j1]-z1[j1-3])*(z1[j1]-z1[j1-3]));   
  	d01=sigma/d0;d05=d01*d01*d01*d01*d01;d06=d01*d05;d010=d05*d05;d012=d06*d06;
   	kst=BSt(i1,j1,i1+3,j1-3);
   	if (kst>=0) {ULJ=0.0;}
   	else {ULJ=-0.5*kst*((5*d12-6*d10)+(5*d012-6*d010));}
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
				us1=0.0;uus1=0.0; us2=0.0; uus2=0.0;
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
         				us0=us0+us1+us2; uus0=uus0+uus1+uus2;  
          			}
         		}
       			us=us+us0; uus=uus+uus0;
        	}
	}
}
// ~~~~~~~~~~~~Calculation of Exculded Volume between any two beads~~~~~~~~~~~~
float LJ0(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
  	float d,d1,d6,d12,r1,r2,ULJ;
  	r1=0.0;r2=0.0;
  	d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); 
  	if (fmod(i1,3)==0&&fmod(j1,3)==0) {r1=2.15;r2=2.15;}
/*   else if (fmod(i1+2,3)==0&&fmod(j1+2,3)==0&&abs(j1-i1)>12) {r1=7.8;r2=7.8;}*/
  	else {r1=R[i1];r2=R[j1];}
  	if (d<=(r1+r2))
  	{
  		d1=(r1+r2)/(1.09*d); d6=d1*d1*d1*d1*d1*d1; d12=d6*d6;
  /*  ULJ=1.04*(d12-d6);*/
  		ULJ=(4.0*e*(d12-d6)+e);
   	}
  	else 
	{
		ULJ=0.0;
	}
  	return ULJ;
} 
void ExcludedVolume(int i0)     //Exculed Volume
{
      	int i,j,i01,i02;
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
  	if (Move==1)
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

float UCoS(int i1,int j1,int k1,int k2,float x[1000],float y[1000],float z[1000])
{
 	float dco1,Ucos,kco,dco2,Ucos1=0.0,Ucos2=0.0,dco_1,dco_2,kdco_1,kdco_2;
 	kco=BSt(i1,j1,k1,k2);
 	dco1=sqrt((x[i1]-x[k1])*(x[i1]-x[k1])+(y[i1]-y[k1])*(y[i1]-y[k1])+(z[i1]-z[k1])*(z[i1]-z[k1])); 
 	dco2=sqrt((x[k2]-x[j1])*(x[k2]-x[j1])+(y[k2]-y[j1])*(y[k2]-y[j1])+(z[k2]-z[j1])*(z[k2]-z[j1])); 
 	if (Bulge==1&&j1-k2>=12) {dco_1=dco; dco_2=2*dco; kdco_1=kdco; kdco_2=2*kdco; }
 	else {dco_1=dco; dco_2=dco; kdco_1=kdco; kdco_2=kdco; }   // 针对bulge做修正，如果间隔3nt，则最优距离dco=2dco，kdco=2kdco；
     // if (fmod(t,100000)==0) {printf("%f %f %f %f %f %f\n",dco1,dco_1,kdco_1,dco2,dco_2,kdco_2);}
 	Ucos1=-0.5*(kco-dkco)*((1-exp(-(dco1-dco_1)/kdco_1))*(1-exp(-(dco1-dco_1)/kdco_1))-1);
 	Ucos2=-0.5*(kco-dkco)*((1-exp(-(dco2-dco_2)/kdco_2))*(1-exp(-(dco2-dco_2)/kdco_2))-1);
 	Ucos=Ucos1+Ucos2; 
 	//if (fmod(t,tprint)==0&&Energy==1) printf("%d %d %d %d %f %f %f %f %f %f\n",i1/3,j1/3,k1/3,k2/3,kco,dco1,dco2,Ucos1,Ucos2,Ucos);
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
    				if (fmod(j,3)==0&&((s[i][j]==1&&s[i-3][j+3]==1&&s[i-6][j+6]==1&&s[i+3][j-3]==0)||(ss[i][j]==1&&ss[i-3][j+3]==1&&ss[i-6][j+6]==1&&ss[i+3][j-3]==0)))
    				{
					uco0=0.0;uuco0=0.0;
     					for(k1=i+3;k1<j;k1++)
     					{      
      						if (fmod(k1,3)==0)
      						{ 
       							for(k2=k1+12;k2<=N0-4;k2++)
       							{
        							if (fmod(k2,3)==0)
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
            										uco0=uco0+uco1; uuco0=uuco0+uuco1; 
           									}
          									else 
										{
											Bulge=0;
										}
          								}
        								if (k2>j+9&&j>k1&&k1>i+9&&(j<=k1+6)&&((s[k1][k2]==1&&s[k1-3][k2+3]==1&&s[k1+3][k2-3]==0)||(ss[k1][k2]==1&&ss[k1-3][k2+3]==1&&ss[k1+3][k2-3]==0)))   //L_loop1>=1 && L_loop0.5>=0 and <=2      Pseudoknot
         								{
         									uco1=0.0; uuco1=0.0; Pseudoknot=1;
         									if (s[i][j]==1&&s[k1][k2]==1&&s[i+3][j-3]==0&&s[k1+3][k2-3]==0) 
										{
											uco1=UCoS(k1,k2,j,i,x,y,z);
										}
         									if (Energy==0&&ss[i][j]==1&&ss[k1][k2]==1&&ss[i+3][j-3]==0&&ss[k1+3][k2-3]==0) 
										{
											uuco1=UCoS(k1,k2,j,i,xx,yy,zz);
										}
         									uco0=uco0+uco1; uuco0=uuco0+uuco1;        
          								} 
        								else 
									{
										Pseudoknot=0;
									}
         							}  
        						}
       						}
      					}
        				uco=uco+uco0; uuco=uuco+uuco0;
    				}
   			}
  		}
 	}  
}
void ENERGY()
{
  	int i,j;
  	float u0=0.0;
  	Energy=1;
  	if (fmod(t,tenergy)==0)
  	{
          	ulj=0.0;ub=0.0;uN=0.0;us=0.0;
          	U=0.0; uco=0.0;utrs=0;
          	for(i=1;i<=N0;i++) 
          	{
                    	c[i]=0;for(j=1;j<=N0;j++) 
                    	{
                          	s[i][j]=0;
                    	}
          	}
          	for(i=1;i<=N0-1;i++) 
          	{
                   	ulj0=0.0;
                   	for (j=i+1;j<=N0;j++)	  
                   	{
                          	ulj1=0.0;uc1=0.0;
                          	ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                                          
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
                 	c[i]=0; cm[i]=0; 
                   	for(j=1;j<=N0;j++) 
                    	{
                    		s[i][j]=0; sm[i][j]=0;
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
                                                     			c[i]=1;c[j]=1; s[i][j]=1; cm[i]=1;cm[j]=1;sm[i][j]=1;                                               										if(fmod(t,tbp)==0||t==1)
                                                     			{
                                                           			fprintf(fp_sec_stru[tem_j],"%d %c %d %c\n",i/3,type[i],j/3,type[j]);fflush(fp_sec_stru[tem_j]);
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
   		if(fmod(t,tbp)==0||t==1) 
   		{
               		fprintf(fp_sec_stru[tem_j],"\n");fflush(fp_sec_stru[tem_j]);
   		}
   		BaseStacking(); 
   		CoaxialStacking(); //Triplex(); //TerminalMismatch(); DanglingEnd();
   		U=ulj+ub+uN+us+uc+uco/*+uTri+umis+udang*/; l++;
   		if (U<=Umin) {Umin=U;}
   		else {Umin=Umin;}
  		Up0=Up0+U; 
   		Up=Up0/(t/tenergy+1);
   		fprintf(fp_Energy[tem_j],"%d %f %f %f %f %f %f %f\n",t,U,ulj,ub,uN,us,uc,uco); fflush(fp_Energy[tem_j]);
   /* printf("%d %f %f %f %f\n",t,U,ub,uN,us);*/
   	}
}

//Out put some structural parameters and input value
void OutputPara()
{
 	time_t timep; time(&timep);fprintf(fp9,"This program starts at :%s\n",ctime(&timep)); //Output time;
 	fprintf(fp9,"Init Temperature %f Alpha %f Steps is %d\nAnnealing T:",t0,Alpha,total);
// fprintf(fp9,"Arithmetic cooling with 5.0 C \n");
 	fprintf(fp9,"tconf1 %d tconf2 %d tenergy %d tbp %d tG %d\n",tconf1,tconf2,tenergy,tbp,tG);
 	fprintf(fp9,"Input parameters: N %d Nbp %d CNa %f CMg %f initial T %f\n",N,Nbp,CNa,CMg,t0);
 	fprintf(fp9,"Folding= %d  OPTIMIZE_W=  %d  salt= %d\n",Folding,OPTIMIZE_W,salt);  fflush(fp9);
}
void FOLD(int i0)
{
    	int i;
    	float prob;
    	prob=rand()/(RAND_MAX+1.); 
      	if ((prob<0.5&&prob>0)&&i0>12&&i0<N0-12&&fmod(i0,3)!=0)    //3nt fragment translation;
      	{
        	Move=3;step=step1_1; 
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
      	else if ((prob<0.75&&prob>=0.5)&&i0>15&&i0<N0-15&&fmod(i0,3)!=0)    //4nt fragment translation;
      	{
        	Move=4;step=step1_1; 
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
      else if ((prob<1.0&&prob>=0.75)&&i0>18&&i0<N0-18&&fmod(i0,3)!=0)    //5nt fragment translation;
      {
        	Move=5;step=step1_1; 
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
                                         	xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
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
                                           	Pivot(i,i0);
                                    	}    
                               	} //pivot
                       	}
                    	else                      //N is chosen
                 	{
                        	for(i=1;i<=N0;i++)
                               	{
                                     	if (i!=i0) 
                                     	{
                                           	xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
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
                                    		xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
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
                                            	Pivot(i,i0);
                                    	} 
                         	}
                  	}
                	else
                     	{
                          	for (i=1;i<=N0;i++)
                              	{
                                  	if (i!=i0) 
                                     	{
                                    		xx[i]=x[i];yy[i]=y[i];zz[i]=z[i]; xxx[i]=x[i];yyy[i]=y[i];zzz[i]=z[i];
                                      	}
                                  	else 
                                     	{
                                         	MoveN(i0);
                                       	}
                         	}            
              		}
      		}	
	}
        if (Folding!=0)
        {
        	ExcludedVolume(i0);    
          	Electrostatic(i0);    
         	Bonded(i0); //Excluded volume  &&  electrostatic  && bonded potential  */
        	BasePairing(); 
          	BaseStacking();  
          	if(t0<30)
          	{
            		CoaxialStacking(); //Triplex(); //TerminalMismatch(); DanglingEnd();
            	}
   	}        //No in Optimization
}
void Optimize(int i0)
{
   	int i,j;
   	float ub1,uub1,ue1,uue1,ud1,uud1;
   	if (OPTIMIZE==0)
   	{
      		step=step2; t_OPT=t_OPT1;
      		for(i=1;i<=N0;i++)
      		{
         		if (i!=i0) 
         		{
             			xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
         		}
         		else
         		{     
             			rand01=rand()/(RAND_MAX+1.);
             			phi=rand01*pi;
             			rand01=rand()/(RAND_MAX+1.);
             			theta=rand01*2.*pi;  
             			rand01=rand()/(RAND_MAX+1.);
             			xx[i]=x[i]+step*rand01*sin(phi)*cos(theta);
             			yy[i]=y[i]+step*rand01*sin(phi)*sin(theta);
             			zz[i]=z[i]+step*rand01*cos(phi);
         		}   
       		}    
     	}
     	if (OPTIMIZE==1)
     	{
           	step=step3;  t_OPT=t_OPT2;
           	for(i=1;i<=N0;i++)
           	{
                	if (i==N1) 
                	{
                     		xx[i]=x[i];yy[i]=y[i];zz[i]=z[i];
                	}
                	else
                	{
                    		rand01=rand()/(RAND_MAX+1.);
                    		phi=rand01*pi;
                    		rand01=rand()/(RAND_MAX+1.);
                    		theta=rand01*2.*pi;  
                    		rand01=rand()/(RAND_MAX+1.);
                    		xx[i]=x[i]+step*rand01*sin(phi)*cos(theta);
                    		yy[i]=y[i]+step*rand01*sin(phi)*sin(theta);
                    		zz[i]=z[i]+step*rand01*cos(phi);
                	}
           	}    
     	}
     	if (OPTIMIZE==2) 
     	{
             	FOLD(i0);
             	t_OPT=t_OPT2;
     	} 
   	BasePairing();    // Base Pairing

   	if (OPTIMIZE==0)
   	{
   //Excluded volume potential 
      		i=i0; ulj0=0.0; uulj0=0.0;
      		for (j=1;j<=N0;j++)	  
      		{ 
        		if (j!=i)
        		{  
            			ulj1=0;uulj1=0;
            			ulj1=LJ0(i,j,x,y,z); ulj=ulj+ulj1;                          
            			uulj1=LJ0(i,j,xx,yy,zz); uulj=uulj+uulj1;   
        		} 
      		}
   //Electrostatic
     		if (salt==0) 
     		{
          		uc=0.0;uuc=0.0;
     		}
     		else
     		{
         		if (fmod(i0+2,3)!=0||i0==4||i0==N0) 
         		{
              			uc=0.0;uuc=0.0;
         		}
         		else 
         		{
             			i=i0; uc0=0.0;uuc0=0.0;
            	 		for (j=4;j<=N0-3;j++)  //4:the first P is still uncharged  
             			{
               				uc1=0.0;uuc1=0.0;
               				if (fmod(j+2,3)==0&&j!=i)
               				{
                   				uc1=QQ(i,j,x,y,z); uc0=uc0+uc1;
                   				uuc1=QQ(i,j,xx,yy,zz); uuc0=uuc0+uuc1;            
               				}
             			}
             			uc=uc+uc0; uuc=uuc+uuc0;
         		}
       		}

     // bonded potential
      		if (fmod((i0+1),3)==0)
     	 	{
           		ub=PC(i0+1,x,y,z)+CP(i0+1,x,y,z)+CN(i0+1,x,y,z);
           		uub=PC(i0+1,xx,yy,zz)+CP(i0+1,xx,yy,zz)+CN(i0+1,xx,yy,zz); 
           		if (i0==2) 
           		{
                 		ue=PCP(i0+1,x,y,z)+CPC(i0+1,x,y,z)+PCN(i0+1,x,y,z)+NCP(i0+1,x,y,z);
                 		uue=PCP(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz);
                 		ud=PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+4,x,y,z);
                 		uud=PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+4,xx,yy,zz);
           		}	 
           		else if (i0==N0-2) 
           		{
                		ue=CPC(i0-2,x,y,z)+PCP(i0+1,x,y,z)+PCN(i0+1,x,y,z)+NCP(i0+1,x,y,z); 
                		uue=CPC(i0-2,xx,yy,zz)+PCP(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz);
                		ud=PCPC(i0-2,x,y,z)+CPCP(i0-2,x,y,z)+CPCN(i0+1,x,y,z)+NCPC(i0-2,x,y,z);
                		uud=PCPC(i0-2,xx,yy,zz)+CPCP(i0-2,xx,yy,zz)+CPCN(i0+1,xx,yy,zz)+NCPC(i0-2,xx,yy,zz);
           		}
           		else
           		{
                 		ue=CPC(i0-2,x,y,z)+PCP(i0+1,x,y,z)+CPC(i0+1,x,y,z)+PCN(i0+1,x,y,z)+NCP(i0+1,x,y,z); 
                 		uue=CPC(i0-2,xx,yy,zz)+PCP(i0+1,xx,yy,zz)+CPC(i0+1,xx,yy,zz)+PCN(i0+1,xx,yy,zz)+NCP(i0+1,xx,yy,zz);
                 		ud=PCPC(i0-2,x,y,z)+CPCP(i0-2,x,y,z)+PCPC(i0+1,x,y,z)+CPCP(i0+1,x,y,z)+NCPC(i0+1,x,y,z)+CPCN(i0+1,x,y,z)+NCPC(i0-2,x,y,z)+CPCN(i0+4,x,y,z);
                 		uud=PCPC(i0-2,xx,yy,zz)+CPCP(i0-2,xx,yy,zz)+PCPC(i0+1,xx,yy,zz)+CPCP(i0+1,xx,yy,zz)+NCPC(i0+1,xx,yy,zz)+CPCN(i0+1,xx,yy,zz)+NCPC(i0-2,xx,yy,zz)+CPCN(i0+4,xx,yy,zz); 
           		}  
       		}
       		if (fmod((i0+2),3)==0)
       		{    
            		if (i0==1) 
            		{
                		ub=PC(i0+2,x,y,z);uub=PC(i0+2,xx,yy,zz); 
                		ue=PCP(i0+2,x,y,z)+PCN(i0+2,x,y,z); uue=PCP(i0+2,xx,yy,zz)+PCN(i0+2,xx,yy,zz); 
                		ud=PCPC(i0+2,x,y,z);uud=PCPC(i0+2,xx,yy,zz);
            		}
           	 	else if (i0==4)  
            		{
                  		ub=PC(i0+2,x,y,z)+CP(i0-1,x,y,z);
                  		uub=PC(i0+2,xx,yy,zz)+CP(i0-1,xx,yy,zz);
                  		ue=PCP(i0-1,x,y,z)+CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+NCP(i0-1,x,y,z)+PCN(i0+2,x,y,z); 
                  		uue=PCP(i0-1,xx,yy,zz)+CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+NCP(i0-1,xx,yy,zz)+PCN(i0+2,xx,yy,zz); 
                  		ud=PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+NCPC(i0-1,x,y,z)+CPCN(i0+2,x,y,z);
                  		uud=PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+CPCN(i0+2,xx,yy,zz);
            		}  
            		else if (i0==N0-3)
            		{
                  		ub=PC(i0+2,x,y,z)+CP(i0-1,x,y,z);uub=PC(i0+2,xx,yy,zz)+CP(i0-1,xx,yy,zz); 
                  		ue=PCP(i0-1,x,y,z)+CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+NCP(i0-1,x,y,z)+PCN(i0+2,x,y,z); 
                  		uue=PCP(i0-1,xx,yy,zz)+CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+NCP(i0-1,xx,yy,zz)+PCN(i0+2,xx,yy,zz); 
                  		ud=CPCP(i0-4,x,y,z)+PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+NCPC(i0-1,x,y,z)+CPCN(i0+2,x,y,z);
                  		uud=CPCP(i0-4,xx,yy,zz)+PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+CPCN(i0+2,xx,yy,zz);
            		}	 
            		else if (i0==N0)
            		{
                  		ub=CP(i0-1,x,y,z);uub=CP(i0-1,xx,yy,zz); 
                  		ue=CPC(i0-1,x,y,z)+NCP(i0-1,x,y,z); uue=CPC(i0-1,xx,yy,zz)+NCP(i0-1,xx,yy,zz); 
                  		ud=CPCP(i0-4,x,y,z);uud=CPCP(i0-4,xx,yy,zz);
            		}
            		else 
            		{
                  		ub=PC(i0+2,x,y,z)+CP(i0-1,x,y,z);
                  		uub=PC(i0+2,xx,yy,zz)+CP(i0-1,xx,yy,zz); 
                  		ue=PCP(i0-1,x,y,z)+CPC(i0-1,x,y,z)+PCP(i0+2,x,y,z)+NCP(i0-1,x,y,z)+PCN(i0+2,x,y,z); 
                  		uue=PCP(i0-1,xx,yy,zz)+CPC(i0-1,xx,yy,zz)+PCP(i0+2,xx,yy,zz)+NCP(i0-1,xx,yy,zz)+PCN(i0+2,xx,yy,zz); 
                  		ud=CPCP(i0-4,x,y,z)+PCPC(i0-1,x,y,z)+CPCP(i0-1,x,y,z)+PCPC(i0+2,x,y,z)+NCPC(i0-1,x,y,z)+CPCN(i0+2,x,y,z);
                  		uud=CPCP(i0-4,xx,yy,zz)+PCPC(i0-1,xx,yy,zz)+CPCP(i0-1,xx,yy,zz)+PCPC(i0+2,xx,yy,zz)+NCPC(i0-1,xx,yy,zz)+CPCN(i0+2,xx,yy,zz);
            		}
       		}
       		if (fmod(i0,3)==0)
       		{
             		ub=CN(i0,x,y,z);uub=CN(i0,xx,yy,zz);
             		ue=NCP(i0,x,y,z)+PCN(i0,x,y,z); uue=NCP(i0,xx,yy,zz)+PCN(i0,xx,yy,zz); 
             		if (i0==3) 
             		{ 
                      		ud=NCPC(i0,x,y,z); uud=NCPC(i0,xx,yy,zz);
             		}
             		else 
             		{
                     		ud=NCPC(i0,x,y,z)+CPCN(i0,x,y,z); uud=NCPC(i0,xx,yy,zz)+CPCN(i0,xx,yy,zz);
             		}
        	}
    	}
    	if (OPTIMIZE==1)
    	{
   //Excluded volume potential 
      		for (i=1;i<=N0-1;i++)
      		{
          		ulj0=0.0; uulj0=0.0;
          		for (j=i+1;j<=N0;j++)	  
          		{	 
              			ulj1=0;uulj1=0;
              			ulj1=LJ0(i,j,x,y,z); ulj0=ulj0+ulj1;                          
              			uulj1=LJ0(i,j,xx,yy,zz); uulj0=uulj0+uulj1;    
          		}
          		ulj=ulj+ulj0; uulj=uulj+uulj0; 
      		}  
   //Electrostatic
      		if (salt==0) {uc=0.0;uuc=0.0;}
      		else
      		{
            		for(i=4;i<=N0-6;i++)       //4:the first P is  uncharged  
            		{
                   		uc0=0.0;uuc0=0.0;
                   		if (fmod(i+2,3)==0)
                   		{
                         		for (j=i+3;j<=N0-3;j++)  //N0-3:the last P is still uncharged  
                         		{
                               			uc1=0.0;uuc1=0.0;
                               			if (fmod(j+2,3)==0)  
                               			{
                                       			uc1=QQ(i,j,x,y,z); uc0=uc0+uc1;   uuc1=QQ(i,j,xx,yy,zz); uuc0=uuc0+uuc1;
                               			}
                         		}
                         		uc=uc+uc0; uuc=uuc+uuc0;
                   		}
             		}
      		}
   //Bonded Energy
      		for(i=1;i<=N0;i++)
      		{
            		ub1=0.0;uub1=0.0;ue1=0.0;uue1=0.0;ud1=0.0;uud1=0.0;
            		if (fmod(i,3)==0)
            		{
                   		ub1=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z);  
                   		uub1=PC(i,xx,yy,zz)+CP(i,xx,yy,zz)+CN(i,xx,yy,zz);
                   		if (i==3) 
                   		{
                           		ue1=PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z);  
                           		uue1=PCP(i,xx,yy,zz)+CPC(i,xx,yy,zz)+PCN(i,xx,yy,zz)+NCP(i,xx,yy,zz); 
                           		ud1=PCPC(i,x,y,z)+CPCP(i,x,y,z)+NCPC(i,x,y,z);  
                           		uud1=PCPC(i,xx,yy,zz)+CPCP(i,xx,yy,zz)+NCPC(i,xx,yy,zz);
                   		}
                   		else if (i==N0-1) 
                   		{
                           		ue1=PCP(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z);  
                           		uue1=PCP(i,xx,yy,zz)+PCN(i,xx,yy,zz)+NCP(i,xx,yy,zz); 
                           		ud1=CPCN(i,x,y,z);  
                           		uud1=CPCN(i,xx,yy,zz);
                   		}	 
                   		else 
                   		{
                           		ue1=PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z);  
                           		uue1=PCP(i,xx,yy,zz)+CPC(i,xx,yy,zz)+PCN(i,xx,yy,zz)+NCP(i,xx,yy,zz); 
                           		ud1=PCPC(i,x,y,z)+CPCP(i,x,y,z)+CPCN(i,x,y,z)+NCPC(i,x,y,z);  
                           		uud1=PCPC(i,xx,yy,zz)+CPCP(i,xx,yy,zz)+CPCN(i,xx,yy,zz)+NCPC(i,xx,yy,zz);
                   		}
                   		ub=ub+ub1; 
                   		uub=uub+uub1; 
                   		ue=ue+ue1; 
                   		uue=uue+uue1; 
                   		ud=ud+ud1; 
                   		uud=uud+uud1; //printf("%d %d %f %f %f\n",t,i,ub1,ue1,ud1);
               		}  
       		}
   	}
    	if (OPTIMIZE==2) 
    	{
             	ExcludedVolume(i0);    
             	Electrostatic(i0);   
             	Bonded(i0);
    	}
    	BaseStacking(); 
        CoaxialStacking(); //Triplex();
}
