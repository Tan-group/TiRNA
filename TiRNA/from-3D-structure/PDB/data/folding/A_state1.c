/*
This program is designed by Ya-Zhou Shi at 2018/06/01 to obtain the structure information of RNA (e.g., HIV-Tar) based on the conformation file;
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define e 0.26
#define Beta 0.60          //AU=Beta*GC
#define BetaGC 0.93
#define gamma 0.1          //The energy threshold of base-pairing for statistics (碱基配对的能量Ubp0<gamma*Ubp算作配对形成)
#define dNN1 8.55
#define dNN2 9.39
#define A -3.5             //base-pairing strength of GC
#define sigma 4.786
#define GUend 0.25         //GU or AU end penalty
#define dco 5.0            //Distance between NN CoaxialStacking
#define dkco 0.0           //0 or 1
#define kdco 2.5
#define bl 5.45
#define IMg 3.0       //2:2 IMg=4.0; 2:1 IMg=3.0;
#define N_thread 15
FILE *input[N_thread],*fp_se[N_thread],*fc;
int N,rf,bp,s[1000][1000],c[10000];
char type[10000];
float x[10000],y[10000],z[10000],Q[10000],f[10000],R[10000];
float uN,temperature[20]={25.0,31.0,37.0,43.0,49.0,55.0,63.0,71.0,78.0,86.0,94.0,102.0,110.0,120.0,130.0};
int Nstep,N0;
int xi1,xj1,xi2,xj2,xi3,xj3,xi4,xj4,xi5,xj5;
int stem_length[1000],stem_i[1000][1000],stem_j[1000][1000],ter[1000][1000],ds[1000][1000];
float u1,u2;
int main()
{
   	int id1,t1,i1,j1;
   	float HB();
   	void BasePairing(),secondary();
   	char filename[30];
   //float P_NF,P_NU,P_I;
   	fc=fopen("ch_0.dat","r+");
   	int i=0;
   	while(!feof(fc))
   	{
           	fscanf(fc,"%d %d %s %f %f %f %f %f %f\n",&t1,&id1,&type[i],&x[i],&y[i],&z[i],&R[i],&Q[i],&f[i]);
           	i++;
   	}
   	fclose(fc);
   	N=i;N0=i;
   	//printf("N0 %d\n",N0);
   	for(rf=0;rf<N_thread;rf++)
   	{
           	sprintf(filename,"conf_%d.dat",rf);     input[rf]=fopen(filename,"r+");     //input: readin the conformational file
           	sprintf(filename,"secondary_stru_%d.dat",rf); fp_se[rf]=fopen(filename,"w+"); 
            //sprintf(filename,"basepairxx_%d.dat",rf); fp1[rf]=fopen(filename,"w+"); 
   	} 
   	j1=0;
   	for(rf=0;rf<N_thread;rf++)
   	{
         	while(!feof(input[rf]))
         	{
                	for(i1=j1*N+1;i1<=j1*N+N;i1++)
               		{
                       		fscanf(input[rf],"%d %d %s %f %f %f %f %f %f\n",&t1,&id1,&type[i1-j1*N],&x[i1-j1*N],&y[i1-j1*N],&z[i1-j1*N],&R[i1-j1*N],&Q[i1-j1*N],&f[i1-j1*N]);
                 	} // 读入构象；i1: No. of nt; j1: No. of conf.;
                 	Nstep=j1+1;
                // fprintf(fp1[rf],"the %dth conformation\n",Nstep);fflush(fp1[rf]);
                 	BasePairing();
		 	secondary();
                 //fprintf(fp1[rf],"\n");fflush(fp1[rf]);   
                 	j1++;     
        	}
    	}
    	for(rf=0;rf<N_thread;rf++)
    	{
         	fclose(input[rf]);fclose(fp_se[rf]);//fclose(fp1[rf]);
    	}
    	//printf("***************** %d\n",j1);
      	return 0;
}

float HB(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])  //The energy calculation of hydrogen bonds;
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
    	else {	UHB=0.0;}
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
void secondary()
{
   	int i,j,order=1,k1,k2,k3,N_stem,kk,j_[20],N_stem_all,kN;
   	int pseudoknot=0;
   	N_stem_all=0;
   	for (i=0;i<20;i++)  
   	{
      		stem_length[i]=0;        //l_stem[i]:length of i-th stem
      		j_[i]=0;
      		for (j=0;j<20;j++) 
      		{
      			stem_i[i][j]=0; stem_j[i][j]=0;
      		} //The j-th base-pair (stem_i-stem_j) in stem i;
    	}      
    	for(i=3;i<=N0-12;i=i+3) 
    	{
    		for(j=i+6;j<=N0;j=j+3)
    		{
    			ter[i][j]=0;
    		}
    	}  
    	for (i=3;i<=(N0-12);i++)  
    	{   	    
      		if (fmod(i,3)==0&&c[i]==1)  
      		{
         		for (j=i+12;j<=N0;j++)  
         		{
           			if (fmod(j,3)==0&&c[j]==1)  
           			{
             				if (s[i][j]==1)  
             				{                            
						if(i>6&&c[3]==1)
						{
							if((s[i+3][j-3]==1&&s[i-3][j+3]==1)||(s[i-3][j+3]==1&&s[i+3][j-3]==0))
                             				{
                             					order=order; 
                             					stem_length[order]=stem_length[order]+1;
                             					k3=stem_length[order];
                             					stem_i[order][k3]=i;  
                             					stem_j[order][k3]=j;
                             							//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
                             				}                        
                             				else
                             				{
                             					kk=0;
                             					if(N0>150) 	{kN=24;}
                             					else		{kN=9;}
                             					for(k1=3;k1<=6;k1=k1+3)
                    						{
                      							for(k2=3;k2<=kN;k2=k2+3)
                      							{
                      								if(i-k1>0&&j+k2<N0)
                      								{
                        								if(s[i-k1][j+k2]==1)
                        								{	 
                        									kk=1; 
                        								}
                        							}
                      							}                                    
                    						}
                    						for(k1=3;k1<=kN;k1=k1+3)
                    						{
                      							for(k2=3;k2<=6;k2=k2+3)
                      							{
                        							if(i-k1>0&&j+k2<N0)
                      								{
                        								if(s[i-k1][j+k2]==1)
                        								{	 
                        									kk=1; 
                        								}
                        							}
                      							}                                    
                    						}
                    						for(k1=3;k1<=kN-6;k1=k1+3)
                    						{
                      							for(k2=3;k2<=kN-6;k2=k2+3)
                      							{
                        							if(i-k1>0&&j+k2<N0)
                      								{
                        								if(s[i-k1][j+k2]==1)
                        								{	 
                        									kk=1; 
                        								}
                        							}
                      							}                                    
                    						}
                    						if(kk==1)	
                    						{
                    							order=order;
                    							stem_length[order]=stem_length[order]+1;
                             						k3=stem_length[order];
                             						stem_i[order][k3]=i;  
                             						stem_j[order][k3]=j;
                             						//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
                    						}
                    						else		
                    						{
                    							if(stem_length[order]==1)
                    							{
                    								order=order;
                    							}
                    							else
                    							{
                    								order++;
                    								stem_length[order]=stem_length[order]+1;
                             							k3=stem_length[order];
                             							stem_i[order][k3]=i;  
                             							stem_j[order][k3]=j;
                             									//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
                    							}
                    						}
                             				}     
						}
						else
						{
                             				if(s[i+3][j-3]==1)
                             				{
                             					order=order; 
                             					stem_length[order]=stem_length[order]+1;
                             					k3=stem_length[order];
                             					stem_i[order][k3]=i;  
                             					stem_j[order][k3]=j;
                             							//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
                             				}  
                             				else
                             				{
                             					kk=0;
                             					if(N0>150) 	{kN=24;}
                             					else		{kN=9;}
                             					for(k1=3;k1<=6;k1=k1+3)
                    						{
                      							for(k2=3;k2<=kN;k2=k2+3)
                      							{
                        							if(i-k1>0&&j+k2<N0)
                      								{
                        								if(s[i-k1][j+k2]==1)
                        								{	 
                        									kk=1; 
                        								}
                        							}
                      							}                                    
                    						}
                    						for(k1=3;k1<=kN;k1=k1+3)
                    						{
                      							for(k2=3;k2<=6;k2=k2+3)
                      							{
                        							if(i-k1>0&&j+k2<N0)
                      								{
                        								if(s[i-k1][j+k2]==1)
                        								{	 
                        									kk=1; 
                        								}
                        							}
                      							}                                    
                    						}
                    						for(k1=3;k1<=kN-6;k1=k1+3)
                    						{
                      							for(k2=3;k2<=kN-6;k2=k2+3)
                      							{
                        							if(i-k1>0&&j+k2<N0)
                      								{
                        								if(s[i-k1][j+k2]==1)
                        								{	 
                        									kk=1; 
                        								}
                        							}
                      							}                                    
                    						}
                    						if(kk==1)	
                    						{
                    							order=order;
                    							stem_length[order]=stem_length[order]+1;
                             						k3=stem_length[order];
                             						stem_i[order][k3]=i;  
                             						stem_j[order][k3]=j;
                             						//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
                    						}
                    						else		{order++;}
                             				} 
                             			}             
                 				break;                              
                			}
              			}
           		}
         	}
       	}
       	N_stem=order;
       	//printf("****************************%d\n",N_stem);
    	if(N_stem<=5)
    	{
      		for(i=1;i<=N_stem;i++)
       		{	
          		N_stem_all=stem_length[i]+N_stem_all;	
       		}      	
      		void pseu (int kkk);
      		pseudoknot=0;
      		if(N_stem>1&&N_stem<4)
       		{	
       			if(N_stem==2)
       			{
       				j_[1]=(stem_j[1][1]-stem_j[2][stem_length[2]])/3;
       				j_[2]=(stem_j[1][1]-stem_i[2][stem_length[2]])/3;
       				if(j_[1]<0&&j_[2]>0)	
       				{
					pseudoknot=1;
       				//printf("RNA is H-type pseudoknot\n");
       					pseu(1);       				       				
       				}
       				else
       				{
       					pseudoknot=0;
       				}               	
       			}
       			if(N_stem==3)
       			{
       				j_[1]=(stem_j[1][1]-stem_j[2][stem_length[2]])/3;
       				j_[2]=(stem_j[1][1]-stem_i[2][stem_length[2]])/3;
       				if(j_[1]<0&&j_[2]>0)	
       				{
					pseudoknot=1;
                	  		pseu(1);
       				//printf("RNA is H-type pseudoknot\n");
       				}
       				//printf("RNA is H-type pseudoknot\n");}
       				j_[3]=(stem_j[1][1]-stem_j[3][stem_length[3]])/3;
       				j_[4]=(stem_j[1][1]-stem_i[3][stem_length[3]])/3;
       				if(j_[3]<0&&j_[4]>0&&pseudoknot==0)	
       				{
					pseudoknot=1;
                	  		pseu(2);
       				//printf("RNA is H-type pseudoknot\n");
       				}
       				j_[5]=(stem_j[2][1]-stem_j[3][stem_length[3]])/3;
       				j_[6]=(stem_j[2][1]-stem_i[3][stem_length[3]])/3;
       				if(j_[5]<0&&j_[6]>0&&pseudoknot==0)	
       				{
					pseudoknot=1;
       				//printf("RNA is H-type pseudoknot\n");
                	  		pseu(3);
       				}      		
       			}       		
       		} 
   	}
   	char dot[1000];
   	for(i=3;i<=N0;i=i+3)
   	{
   		dot[i]='.';
   	}
   	for(i=3;i<N0;i=i+3)
   	{
   		for(j=i+3;j<=N0;j=j+3)
   		{
   			if(s[i][j]==1)
   			{
   				if(ter[i][j]==1)
   				{
   					dot[i]='[';dot[j]=']';
   				}
   				else
   				{
   					dot[i]='(';dot[j]=')';
   				}
   			}
   		}
   	}
   	for(i=3;i<=N0;i=i+3)
   	{
   		fprintf(fp_se[rf],"%c ",dot[i]);
   	}
   	fprintf(fp_se[rf],"\n");			                 
}

void pseu(int kkk)
{
	float BaseStacking();
	int i,j;
	//printf("-------------%d\n",kkk);
	for(i=3;i<N0;i=i+3)
	{
		for(j=i+3;j<=N0;j=j+3)
		{
			ds[i][j]=0;
		}
	}
/**********************H-type pseudoknot or 3'-end hairpin+H-type pseudoknot********************/
	if(kkk==1)
	{
		for(i=3;i<=N0;i=i+3)
       		{
       			for(j=i+3;j<=N0;j=j+3)
       			{
       				if((i>=stem_i[1][1]&&i<=stem_i[1][stem_length[1]])&&(j>=stem_j[1][stem_length[1]]&&j<=stem_j[1][1]))
       				{
       					if(s[i][j]==1)	{ds[i][j]=1;}
       				}
       			}
       		}      	
       		u1=BaseStacking();
       		for(i=3;i<N0;i=i+3)
		{
			for(j=i+3;j<=N0;j=j+3)
			{
				ds[i][j]=0;
			}
		}
       		for(i=3;i<=N0;i=i+3)
       		{
       			for(j=i+3;j<=N0;j=j+3)
       			{
       				if((i>=stem_i[2][1]&&i<=stem_i[2][stem_length[2]])&&(j>=stem_j[2][stem_length[2]]&&j<=stem_j[2][1]))
       				{
       					if(s[i][j]==1)	{ds[i][j]=1;}
       				}
       			}
       		}
       		u2=BaseStacking();
       		if(u1<u2)
       		{
       			for(i=3;i<=N0;i=i+3)
       			{
       			
       				for(j=i+3;j<=N0;j=j+3)
       				{
       					if((i>=stem_i[2][1]&&i<=stem_i[2][stem_length[2]])||(i>=stem_j[2][stem_length[2]]&&i<=stem_j[2][1]))
       					{       					
       						if(s[i][j]==1)	{ter[i][j]=1;}
       						else		{ter[i][j]=0;}
       					}
       				}
       			}		
       		}
       		else
       		{
       			for(i=3;i<=N0;i=i+3)
       			{
       			
       				for(j=i+3;j<=N0;j=j+3)
       				{
       					if((i>=stem_i[1][1]&&i<=stem_i[1][stem_length[1]])||(i>=stem_j[1][stem_length[1]]&&i<=stem_j[1][1]))
       					{       					
       						if(s[i][j]==1)	{ter[i][j]=1;}
       						else		{ter[i][j]=0;}
       					}
       				}
       			}			
       		}
       	}  
 /*************************************** pseudokont + interloop_helix 1-3 stem****************************/
       	if(kkk==2)
	{
		for(i=3;i<=N0;i=i+3)
       		{
       			for(j=i+3;j<=N0;j=j+3)
       			{
       				if((i>=stem_i[1][1]&&i<=stem_i[1][stem_length[1]])&&(j>=stem_j[1][stem_length[1]]&&j<=stem_j[1][1]))
       				{
       					if(s[i][j]==1)	{ds[i][j]=1;}
       				}
       			}
       		}      	
       		u1=BaseStacking();
       		for(i=3;i<N0;i=i+3)
		{
			for(j=i+3;j<=N0;j=j+3)
			{
				ds[i][j]=0;
			}
		}
       		for(i=3;i<=N0;i=i+3)
       		{
       			for(j=i+3;j<=N0;j=j+3)
       			{
       				if((i>=stem_i[3][1]&&i<=stem_i[3][stem_length[3]])&&(j>=stem_j[3][stem_length[3]]&&j<=stem_j[3][1]))
       				{
       					if(s[i][j]==1)	{ds[i][j]=1;}
       				}
       			}
       		}
       		u2=BaseStacking();
       		if(u1<u2)
       		{
       			for(i=3;i<=N0;i=i+3)
       			{       			
       				for(j=i+3;j<=N0;j=j+3)
       				{
       					if((i>=stem_i[3][1]&&i<=stem_i[3][stem_length[3]])||(i>=stem_j[3][stem_length[3]]&&i<=stem_j[3][1]))
       					{       					
       						if(s[i][j]==1)	{ter[i][j]=1;}
       						else		{ter[i][j]=0;}
       					}
       				}
       			}		
       		}
       		else
       		{
       			for(i=3;i<=N0;i=i+3)
       			{       			
       				for(j=i+3;j<=N0;j=j+3)
       				{
       					if((i>=stem_i[1][1]&&i<=stem_i[1][stem_length[1]])||(i>=stem_j[1][stem_length[1]]&&i<=stem_j[1][1]))
       					{       					
       						if(s[i][j]==1)	{ter[i][j]=1;}
       						else		{ter[i][j]=0;}
       					}
       				}
       			}			
       		}
       	}  
   /************************ pseudokont +  5-end hairpin *******************************/
       	if(kkk==3)
	{
		for(i=3;i<=N0;i=i+3)
       		{
       			for(j=i+3;j<=N0;j=j+3)
       			{
       				if((i>=stem_i[2][1]&&i<=stem_i[2][stem_length[2]])&&(j>=stem_j[2][stem_length[2]]&&j<=stem_j[2][1]))
       				{
       					if(s[i][j]==1)	{ds[i][j]=1;}
       				}
       			}
       		}      	
       		u1=BaseStacking();
       		for(i=3;i<N0;i=i+3)
		{
			for(j=i+3;j<=N0;j=j+3)
			{
				ds[i][j]=0;
			}
		}
       		for(i=3;i<=N0;i=i+3)
       		{
       			for(j=i+3;j<=N0;j=j+3)
       			{
       				if((i>=stem_i[3][1]&&i<=stem_i[3][stem_length[3]])&&(j>=stem_j[3][stem_length[3]]&&j<=stem_j[3][1]))
       				{
       					if(s[i][j]==1)	{ds[i][j]=1;}
       				}
       			}
       		}
       		u2=BaseStacking();
       		if(u1<u2)
       		{
       			for(i=3;i<=N0;i=i+3)
       			{       			
       				for(j=i+3;j<=N0;j=j+3)
       				{
       					if((i>=stem_i[3][1]&&i<=stem_i[3][stem_length[3]])||(i>=stem_j[3][stem_length[3]]&&i<=stem_j[3][1]))
       					{       					
       						if(s[i][j]==1)	{ter[i][j]=1;}
       						else		{ter[i][j]=0;}
       					}
       				}
       			}		
       		}
       		else
       		{
       			for(i=3;i<=N0;i=i+3)
       			{
       			
       				for(j=i+3;j<=N0;j=j+3)
       				{
       					if((i>=stem_i[2][1]&&i<=stem_i[2][stem_length[2]])||(i>=stem_j[2][stem_length[2]]&&i<=stem_j[2][1]))
       					{       					
       						if(s[i][j]==1)	{ter[i][j]=1;}
       						else		{ter[i][j]=0;}
       					}
       				}
       			}		
       		}
       	}    		    		    				
}
float BSt(int i1,int j1,int k1,int k2)  //或许可以用四维数组书写，目前暂且如此吧！
{
	float B=0.0,B0=0.0,T;
	T=300.0;
             // 5'-i1-k1-
             // 3'-j1/k2-
     	if 	(((type[i1]=='A'&&type[k1]=='A')&&(type[j1]=='U'&&type[k2]=='U'))
		||((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='A'&&type[k2]=='A'))) 	{B=-6.82-T*0.001*(-19.0-B0);}  // AA/UU
	else if (((type[i1]=='A'&&type[k1]=='C')&&(type[j1]=='U'&&type[k2]=='G')) 
     		||((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='C'&&type[k2]=='A'))) 	{B=-11.4-T*0.001*(-29.5-B0);}  // AC/UG
	else if (((type[i1]=='A'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='C'))
     		||((type[i1]=='C'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='A'))) 	{B=-10.48-T*0.001*(-27.1-B0);} // AG/UC
	else if (((type[i1]=='A'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='A'))) 	{B=-3.21-T*0.001*(-8.6-B0);}   // AG/UU
	else if ((type[i1]=='A'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='A'))   	{B=-9.38-T*0.001*(-26.7-B0);}  // AU/UA
	else if (((type[i1]=='A'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='A'))) 	{B=-8.81-T*0.001*(-24.0-B0);}  // AU/UG
	else if (((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='A'&&type[k2]=='C')) 
     		||((type[i1]=='C'&&type[k1]=='A')&&(type[j1]=='G'&&type[k2]=='U'))) 	{B=-10.44-T*0.001*(-26.9-B0);} // UG/AC
	else if (((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='C'&&type[k2]=='C'))
     		||((type[i1]=='C'&&type[k1]=='C')&&(type[j1]=='G'&&type[k2]=='G'))) 	{B=-13.39-T*0.001*(-32.7-B0);} // GG/CC
	else if ((type[i1]=='C'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='C'))   	{B=-10.64-T*0.001*(-26.7-B0);} // CG/GC
	else if (((type[i1]=='C'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='C'))) 	{B=-5.61-T*0.001*(-13.5-B0);}  // CG/GU
	else if (((type[i1]=='C'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='C'))) 	{B=-12.11-T*0.001*(-32.2-B0);} // CU/GG
	else if (((type[i1]=='U'&&type[k1]=='C')&&(type[j1]=='A'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='A')&&(type[j1]=='C'&&type[k2]=='U'))) 	{B=-12.44-T*0.001*(-32.5-B0);} // UC/AG
	else if ((type[i1]=='G'&&type[k1]=='C')&&(type[j1]=='C'&&type[k2]=='G'))   	{B=-14.88-T*0.001*(-36.9-B0);} // GC/CG
	else if (((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='C'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='C')&&(type[j1]=='G'&&type[k2]=='G'))) 	{B=-8.33-T*0.001*(-21.9-B0);}  // GG/CU
	else if (((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='C'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='C')&&(type[j1]=='U'&&type[k2]=='G'))) 	{B=-12.59-T*0.001*(-32.5-B0);} // GU/CG
	else if ((type[i1]=='U'&&type[k1]=='A')&&(type[j1]=='A'&&type[k2]=='U'))   	{B=-7.69-T*0.001*(-20.5-B0);}  // UA/AU
	else if (((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='A'&&type[k2]=='U'))
     		||((type[i1]=='U'&&type[k1]=='A')&&(type[j1]=='G'&&type[k2]=='U'))) 	{B=-6.99-T*0.001*(-19.3-B0);}  // UG/AU
	else if (((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='A'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='A')&&(type[j1]=='U'&&type[k2]=='U'))) 	{B=-12.83-T*0.001*(-37.3-B0);} // UU/AG
	else if (((type[i1]=='U'&&type[k1]=='U')&&(type[j1]=='G'&&type[k2]=='G'))
     		||((type[i1]=='G'&&type[k1]=='G')&&(type[j1]=='U'&&type[k2]=='U'))) 	{B=-13.47-T*0.001*(-44.9-B0);} // UU/GG
	else if ((type[i1]=='U'&&type[k1]=='G')&&(type[j1]=='G'&&type[k2]=='U'))   	{B=-9.26-T*0.001*(-30.8-B0);}  // UG/GU
	else if ((type[i1]=='G'&&type[k1]=='U')&&(type[j1]=='U'&&type[k2]=='G'))   	
	{
        	if 	((type[i1-3]=='G'&&type[k1+3]=='C')&&(type[j1+3]=='C'&&type[k2-3]=='G'))	{B=-16.66-T*0.001*(-50.3-B0);}
        	else                                                                       		{B=-14.59-T*0.001*(-51.2-B0);}
	}// GU/UG (CGUC/CUGG)
	else    {B=0.0*(-6.82-T*0.001*(-19.0-B0));}   //Mismatched stacking
   	return B;
}

float St(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
	float d,d0,kst,ULJ,d1,d5,d6,d10,d12,d01,d05,d06,d010,d012/*,ulj1,ulj2*/;
	float BSt();
  	d=sqrt((x1[i1]-x1[i1+3])*(x1[i1]-x1[i1+3])+(y1[i1]-y1[i1+3])*(y1[i1]-y1[i1+3])+(z1[i1]-z1[i1+3])*(z1[i1]-z1[i1+3])); 
  	d1=4.79/d; d5=d1*d1*d1*d1*d1;d6=d1*d5;d10=d5*d5;d12=d6*d6; 
  	d0=sqrt((x1[j1]-x1[j1-3])*(x1[j1]-x1[j1-3])+(y1[j1]-y1[j1-3])*(y1[j1]-y1[j1-3])+(z1[j1]-z1[j1-3])*(z1[j1]-z1[j1-3]));   
  	d01=4.79/d0;d05=d01*d01*d01*d01*d01;d06=d01*d05;d010=d05*d05;d012=d06*d06;
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
float BaseStacking()
{
    	int i,j;
    	float St();
    	float us1,us0,us;
    	us0=0.0;us=0.0;
    	for(i=3;i<=(N0-12);i=i+3)
    	{  
          	for(j=i+12;j<=N0;j=j+3)
          	{   
              		us1=0.0; 
                   	if (ds[i][j]==1&&ds[i+3][j-3]==1) 
                   	{
                         	us1=St(i,j,x,y,z);
                   	}
                   	us0=us0+us1;
            	}
            	us=us+us0; 
     	}
     	return us;
}
