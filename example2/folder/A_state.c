/****************** 该程序是用来判断结构态，用于wham当中 ***************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define e 0.26
#define Beta 0.60          //AU=Beta*GC
#define BetaGC 0.93
#define gamma 0.1          //The energy threshold of base-pairing for statistics (碱基配对的能量Ubp0<gamma*Ubp算作配对形成)
// hydrogen bond
#define dNN1 8.55
#define dNN2 9.39
#define A -3.5             //base-pairing strength of GC
#define N_thread 15
FILE *input[N_thread],*fp_bp[N_thread],*fpstate,*fc;
int N,rf,bp,s[1000][1000],c[10000],NF=0,NS1=0,NS2=0,NU=0,N_other=0;
char type[10000];
float x[5000],y[5000],z[5000],Q[5000],f[5000],R[5000];
float uN,temperature[20]={25.0,31.0,37.0,43.0,49.0,55.0,63.0,71.0,78.0,86.0,94.0,102.0,110.0,120.0,130.0};
int Nstep,N0;
int xi1,xj1,xi2,xj2,xi3,xj3,xi4,xj4,xi5,xj5;
int main()
{
   	int id1,t1,i1,j1;
   	float HB();
   	void BasePairing(),Stem1_2();
   	char filename[20];
   	fc=fopen("ch_0.dat","r+");
   	int i=0;
   	while(!feof(fc))
   	{
          	fscanf(fc,"%d %d %s %f %f %f %f %f %f\n",&t1,&id1,&type[i],&x[i],&y[i],&z[i],&R[i],&Q[i],&f[i]);
           	i++;
   	}
   	fclose(fc);
   	N=i;
   	for(rf=0;rf<N_thread;rf++)
   	{
            	sprintf(filename,"conf_%d.dat",rf);   input[rf]=fopen(filename,"r+");     //input: readin the conformational file
            	sprintf(filename,"bp_%d.dat",rf);     fp_bp[rf]=fopen(filename,"w+"); 
   	} 
   	fpstate=fopen("state.dat","r+");
   	while(!feof(fpstate))
   	{
   		fscanf(fpstate,"%d %d %d %d %d %d %d %d %d %d\n",&xi1,&xj1,&xi2,&xj2,&xi3,&xj3,&xi4,&xj4,&xi5,&xj5);
   	}
   	fclose(fpstate);
   	xi1=xi1*3;xj1=xj1*3;
   	xi2=xi2*3;xj2=xj2*3;
   	xi3=xi3*3;xj3=xj3*3;
   	xi4=xi4*3;xj4=xj4*3;
   	xi5=xi5*3;xj5=xj5*3;
   	for(rf=0;rf<N_thread;rf++)
   	{
         	i1=0;j1=0;NF=0;NS1=0;NS2=0;NU=0;N_other=0;
         	while(!feof(input[rf]))
         	{
                	for(i1=j1*N+1;i1<=j1*N+N;i1++)
                	{
                       		fscanf(input[rf],"%d %d %s %f %f %f %f %f %f\n",&t1,&id1,&type[i1-j1*N],&x[i1-j1*N],&y[i1-j1*N],&z[i1-j1*N],&R[i1-j1*N],&Q[i1-j1*N],&f[i1-j1*N]);
                 	} // 读入构象；i1: No. of nt; j1: No. of conf.;
                 	Nstep=j1+1;
                 	BasePairing();
                 	Stem1_2(); 
                 	j1++;     
        	}
    	}
    	for(rf=0;rf<N_thread;rf++)
    	{
         	fclose(input[rf]);fclose(fp_bp[rf]);
    	}
      	return 0;
}

float HB(int i1,int j1,float x1[10000],float y1[10000],float z1[10000])  //The energy calculation of hydrogen bonds;
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
      		if ((type[i1]=='G'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='G'))  {hb=Beta*A;}
      		if ((type[i1]=='A'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='A'))  {hb=Beta*A;} 
      		UHB=hb/(1+2.66*(d-8.94)*(d-8.94)+1.37*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.464*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));  
    	}
    	else 
    	{
         	UHB=0.0;
    	}
    	return UHB;
}

void BasePairing()     //The formation of basepairing;
{
    	int i,j,bp0=0;
     	float uN1;
     	bp=0;
     	for (i=1;i<=N;i++)  {c[i]=0; for(j=1;j<=N;j++) {s[i][j]=0;}}
     	for(i=1;i<=N;i++)  
    	{
       		if (fmod(i,3)==0) 
       		{	
       			for(j=i+12;j<=N;j++)  
     			{
           			if (fmod(j,3)==0) 
          			{ 
              				if ((type[i]=='G'&&type[j]=='C')||(type[i]=='C'&&type[j]=='G')
                			||(type[i]=='A'&&type[j]=='U')||(type[i]=='U'&&type[j]=='A')
                			||(type[i]=='G'&&type[j]=='U')||(type[i]=='U'&&type[j]=='G'))
             				{
               					uN1=0.0;bp0=0;
                 				if (c[i]==0&&c[j]==0) 
                				{
                    					uN1=HB(i,j,x,y,z);
                    					if (uN1!=0) 
                   					{
                       						if (((type[i]=='G'&&type[j]=='C')||(type[i]=='C'&&type[j]=='G'))) 
                      						{
                        							c[i]=1;c[j]=1; s[i][j]=1; bp0=1;
                       						}
                       						if (((type[i]=='A'&&type[j]=='U')||(type[i]=='U'&&type[j]=='A')
                          					||(type[i]=='G'&&type[j]=='U')||(type[i]=='U'&&type[j]=='G'))) 
                      						{
                        						c[i]=1;c[j]=1; s[i][j]=1; bp0=1;
                       						}      
                    					}
                  				}
             			 		bp=bp+bp0; uN=uN+uN1;
              				} 
            			}
          		}
        	} 
     	}
}

/////////////////
void Stem1_2()
{
	int Stem1=0,Stem2=0, Stem3=0,Stem4=0,Stem5=0,State=0;
	if(xi1!=0&&xi2==0)
   	{
   		if ((s[xi1][xj1]==1&&s[xi1+3][xj1-3]==1)||(s[xi1][xj1]==1&&s[xi1-3][xj1+3]==1)) 	{Stem1=1;}
   		if (Stem1==1) 										{State=2; NF=NF+1;}
   		else if ((Stem1!=1)&&(bp<=1)) 								{State=0; NU=NU+1;}
   		else                         								{State=1; N_other++;} 
   	}
   	if(xi2!=0&&xi3==0)
   	{
   		if ((s[xi1][xj1]==1&&s[xi1+3][xj1-3]==1)||(s[xi1][xj1]==1&&s[xi1-3][xj1+3]==1)) 	{ Stem1=1;}
   		if ((s[xi2][xj2]==1&&s[xi2+3][xj2-3]==1)||(s[xi2][xj2]==1&&s[xi2-3][xj2+3]==1)) 	{ Stem2=1;}
   		if (Stem1==1&&Stem2==1) 								{State=2; NF=NF+1;}
   		else if ((Stem1!=1&&Stem2!=1)&&(bp<=2)) 						{State=0; NU=NU+1;}
   		else                         								{State=1; N_other++;} 
   	}
   	if(xi3!=0&&xi4==0)
   	{
   		if ((s[xi1][xj1]==1&&s[xi1+3][xj1-3]==1)||(s[xi1][xj1]==1&&s[xi1-3][xj1+3]==1)) 	{ Stem1=1;}
   		if ((s[xi2][xj2]==1&&s[xi2+3][xj2-3]==1)||(s[xi2][xj2]==1&&s[xi2-3][xj2+3]==1)) 	{ Stem2=1;}
   		if ((s[xi3][xj3]==1&&s[xi3+3][xj3-3]==1)||(s[xi3][xj3]==1&&s[xi3-3][xj3+3]==1)) 	{ Stem3=1;}
   		if (Stem1==1&&Stem2==1&&Stem3==1) 							{State=2; NF=NF+1;}
   		else if ((Stem1!=1&&Stem2!=1&&Stem3!=1)&&(bp<=3)) 					{State=0; NU=NU+1;}
   		else                         								{State=1; N_other++;} 
   	}
   	if(xi4!=0&&xi5==0)
   	{
   		if ((s[xi1][xj1]==1&&s[xi1+3][xj1-3]==1)||(s[xi1][xj1]==1&&s[xi1-3][xj1+3]==1)) 	{ Stem1=1;}
   		if ((s[xi2][xj2]==1&&s[xi2+3][xj2-3]==1)||(s[xi2][xj2]==1&&s[xi2-3][xj2+3]==1)) 	{ Stem2=1;}
   		if ((s[xi3][xj3]==1&&s[xi3+3][xj3-3]==1)||(s[xi3][xj3]==1&&s[xi3-3][xj3+3]==1)) 	{ Stem3=1;}
   		if ((s[xi4][xj4]==1&&s[xi4+3][xj4-3]==1)||(s[xi4][xj4]==1&&s[xi4-3][xj4+3]==1)) 	{ Stem4=1;}
   		if (Stem1==1&&Stem2==1&&Stem3==1&&Stem4==1) 						{State=2; NF=NF+1;}
   		else if ((Stem1!=1&&Stem2!=1&&Stem3!=1&&Stem4!=1)&&(bp<=4)) 				{State=0; NU=NU+1;}
   		else                         								{State=1; N_other++;} 
   	}
   	if(xi5!=0)
   	{
   		if ((s[xi1][xj1]==1&&s[xi1+3][xj1-3]==1)||(s[xi1][xj1]==1&&s[xi1-3][xj1+3]==1)) 	{ Stem1=1;}
   		if ((s[xi2][xj2]==1&&s[xi2+3][xj2-3]==1)||(s[xi2][xj2]==1&&s[xi2-3][xj2+3]==1)) 	{ Stem2=1;}
   		if ((s[xi3][xj3]==1&&s[xi3+3][xj3-3]==1)||(s[xi3][xj3]==1&&s[xi3-3][xj3+3]==1)) 	{ Stem3=1;}
   		if ((s[xi4][xj4]==1&&s[xi4+3][xj4-3]==1)||(s[xi4][xj4]==1&&s[xi4-3][xj4+3]==1)) 	{ Stem4=1;}
   		if ((s[xi5][xj5]==1&&s[xi5+3][xj5-3]==1)||(s[xi5][xj5]==1&&s[xi5-3][xj5+3]==1)) 	{ Stem5=1;}
   		if (Stem1==1&&Stem2==1&&Stem3==1&&Stem4==1&&Stem5==1) 					{State=2; NF=NF+1;}
   		else if ((Stem1!=1&&Stem2!=1&&Stem3!=1&&Stem4!=1&&Stem5!=1)&&(bp<=5)) 			{State=0; NU=NU+1;}
   		else                         								{State=1; N_other++;} 
   	}
   	fprintf(fp_bp[rf],"%d %d %d %f %f %d %f %f\n",Nstep,1,State,0.0,0.0,1,0.0,0.0); fflush(fp_bp[rf]);
}
