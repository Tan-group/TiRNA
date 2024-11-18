#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define BetaGC 0.9
#define BetaAU 0.6
#define BetaGU 0.6
#define A -3.5
int duo1,duo2;
char type[1000];
char atom[1000][5],typeA[1000][5],chain[1000][5];
int id[1000],idnu[1000];
float x[1000],y[1000],z[1000];
int N0,s[1000][1000],c[1000];
int ds[1000][1000],ter[1000][1000];
int stem_i[20][20],stem_j[20][20];
int stem_length[20];
float u1,u2;
FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
int main()
{
    	int i;
    	void secondary(),BasePairing();
    	fp=fopen("CG.pdb","r+");
    	fp1=fopen("stem.dat","w+");
   	fp2=fopen("stem_kissing.dat","w+");
    	fp3=fopen("cs.dat","w+");
    	fp4=fopen("RNA_type","w+");
    	fp5=fopen("state.dat","w+");
    	fp6=fopen("sec_struc.dat","w+");
    	i=1;
    	while(!feof(fp)) 
    	{
         	fscanf(fp,"%s %d %s %s %s %d %f %f %f\n",atom[i],&id[i],typeA[i],&type[i],chain[i],&idnu[i],&x[i],&y[i],&z[i]);
         	i++;
     	}
     	N0=i-1;
     	BasePairing();
     	secondary();  
     	return 0;    
     	fclose(fp);fclose(fp1);fclose(fp2);fclose(fp3);fclose(fp5);
} 
float HB(int i1,int j1,float x1[1000],float y1[1000],float z1[1000])
{
 	float d,hb,UHB,d0,d01,d1,d11;
 	d=sqrt((x1[i1]-x1[j1])*(x1[i1]-x1[j1])+(y1[i1]-y1[j1])*(y1[i1]-y1[j1])+(z1[i1]-z1[j1])*(z1[i1]-z1[j1])); //NN
 	if (d>=8.55&&d<=9.39)   //The pairing formation condition
 	{
 		d0=sqrt((x1[i1]-x1[j1-1])*(x1[i1]-x1[j1-1])+(y1[i1]-y1[j1-1])*(y1[i1]-y1[j1-1])+(z1[i1]-z1[j1-1])*(z1[i1]-z1[j1-1]));  //NiCj
 		d01=sqrt((x1[j1]-x1[i1-1])*(x1[j1]-x1[i1-1])+(y1[j1]-y1[i1-1])*(y1[j1]-y1[i1-1])+(z1[j1]-z1[i1-1])*(z1[j1]-z1[i1-1])); //CiNj
 		d1=sqrt((x1[i1-2]-x1[j1])*(x1[i1-2]-x1[j1])+(y1[i1-2]-y1[j1])*(y1[i1-2]-y1[j1])+(z1[i1-2]-z1[j1])*(z1[i1-2]-z1[j1]));  //PiNj
 		d11=sqrt((x1[j1-2]-x1[i1])*(x1[j1-2]-x1[i1])+(y1[j1-2]-y1[i1])*(y1[j1-2]-y1[i1])+(z1[j1-2]-z1[i1])*(z1[j1-2]-z1[i1])); //NiPj
     		if ((type[i1]=='G'&&type[j1]=='C')||(type[i1]=='C'&&type[j1]=='G'))  		{hb=BetaGC*A;}
		if ((type[i1]=='G'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='G'))  		{hb=BetaGU*A;}
		if ((type[i1]=='A'&&type[j1]=='U')||(type[i1]=='U'&&type[j1]=='A'))  		{hb=BetaAU*A;} 
  		UHB=hb/(1+2.66*(d-8.94)*(d-8.94)+1.37*((d0-12.15)*(d0-12.15)+(d01-12.15)*(d01-12.15))+0.464*((d1-13.92)*(d1-13.92)+(d11-13.92)*(d11-13.92)));  
  	}
 	else {UHB=0.0;}
 	return UHB;
}

void BasePairing()
{
    	int i,j;
    	float HB();
    	float uN1=0.0;
    	for (i=1;i<=N0;i++)  
    	{
		c[i]=0; 
            	for(j=i+3;j<=N0;j++) 
            	{
                	s[i][j]=0;
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
                        			if (c[i]==0&&c[j]==0) 
                        			{
                            				uN1=HB(i,j,x,y,z);
                            				if (uN1!=0) 
                            				{
                                				c[i]=1;c[j]=1; 
                                				s[i][j]=1;                                                
                             				}
                          			}  
                      			} 
                    		}        
                 	} 
             	}
        }        	
}

void secondary()
{
   	int i,j,k1,k2,order=1,k3,N_stem,kk,j_[20],k_l=0,N_stem_all,j_l_all=0,kN;
   	int i1,j1,i2,j2,i3,j3,i4,j4,i5,j5;
   	int si1,si2,sj1,sj2;
   	int junction,pseudoknot=0;
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
    	for (i=0;i<20;i++)  
   	{
      		stem_length[i]=0;        //l_stem[i]:length of i-th stem
      		j_[i]=0;
      		for (j=0;j<20;j++) 
      		{
      			stem_i[i][j]=0; stem_j[i][j]=0;
      		} //The j-th base-pair (stem_i-stem_j) in stem i;
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
                             							fprintf(fp3,"%d %c %d %c\n",stem_i[order][k3]/3,type[i],stem_j[order][k3]/3,type[j]);
                             						}                        
                             						else
                             						{
                             							kk=0;
                             							if(N0>150) 	{kN=24;}
                             							else		{kN=9;}
                             							for(k1=3;k1<=kN;k1=k1+3)
                    								{
                      									for(k2=3;k2<=kN;k2=k2+3)
                      									{
                      										if(i-k1>0&&j+k2<N0)
                      										{
                        										if(s[i-k1][j+k2]==1) { kk=1; }
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
                             								fprintf(fp3,"%d %c %d %c\n",stem_i[order][k3]/3,type[i],stem_j[order][k3]/3,type[j]);
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
                             									fprintf(fp3,"%d %c %d %c\n",stem_i[order][k3]/3,type[i],stem_j[order][k3]/3,type[j]);
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
                             							fprintf(fp3,"%d %c %d %c\n",stem_i[order][k3]/3,type[i],stem_j[order][k3]/3,type[j]);
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
                        										if(s[i-k1][j+k2]==1) { kk=1; }
                        									}
                      									}                                    
                    								}
                    								for(k1=3;k1<=kN;k1=k1+3)
                    								{
                      									for(k2=3;k2<=6;k2=k2+3)
                      									{
                        									if(i-k1>0&&j+k2<N0)
                      										{
                        										if(s[i-k1][j+k2]==1) { kk=1; }
                        									}
                      									}                                    
                    								}
                    								for(k1=3;k1<=kN-6;k1=k1+3)
                    								{
                      									for(k2=3;k2<=kN-6;k2=k2+3)
                      									{
                        									if(i-k1>0&&j+k2<N0)
                      										{
                        										if(s[i-k1][j+k2]==1) { kk=1; }
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
                             								fprintf(fp3,"%d %c %d %c\n",stem_i[order][k3]/3,type[i],stem_j[order][k3]/3,type[j]);
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
	int structure=0;
	for(i=0;i>20;i++)
	{
		for(j=0;j<20;j++)
		{
			if(stem_i[i][j]==1||stem_j[i][j]==1)
			{
				structure=1;
			}
		}
	}
	if(structure==0)
	{
		fprintf(fp3,"%d %c %d %c\n",0,'A',0,'A');
	}
			
    if(N_stem<=5)
    {
      	for(i=1;i<=N_stem;i++)
       	{
          	N_stem_all=stem_length[i]+N_stem_all;
          	
       	}       	
       	if(N_stem==1)	
       	{
       		junction=1;
                fprintf(fp4,"%d\n",1);
       	}
       	else
       	{
       		fprintf(fp4,"%d\n",2);
       	}
       	fclose(fp4);
      	int pse=0;
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
       				pse=stem_j[1][1];pseudoknot=1;
       				pseu(1);
       				
       				
       			}
       			else
       			{
       				pseudoknot=0;pse=0;
       			}
       			fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",pse/3,0,0,0,0,0,0,0,0,0);
                	fprintf(fp2,"%d %d %d %d\n",0,0,0,0);
                	
       		}
       		if(N_stem==3)
       		{
       			j_[1]=(stem_j[1][1]-stem_j[2][stem_length[2]])/3;
       			j_[2]=(stem_j[1][1]-stem_i[2][stem_length[2]])/3;
       			if(j_[1]<0&&j_[2]>0)	
       			{
       				pse=stem_j[1][1];pseudoknot=1;
       				fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",pse/3,0,0,0,0,0,0,0,0,0);
                	  	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 
                	  	pseu(1);
       			}
       			j_[3]=(stem_j[1][1]-stem_j[3][stem_length[3]])/3;
       			j_[4]=(stem_j[1][1]-stem_i[3][stem_length[3]])/3;
       			if(j_[3]<0&&j_[4]>0&&pseudoknot==0)	
       			{
       				pse=stem_j[1][1];pseudoknot=1;
       				fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",pse/3,0,0,0,0,0,0,0,0,0);
                	  	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 
                	  	pseu(2);
       			}
       			j_[5]=(stem_j[2][1]-stem_j[3][stem_length[3]])/3;
       			j_[6]=(stem_j[2][1]-stem_i[3][stem_length[3]])/3;
       			if(j_[5]<0&&j_[6]>0&&pseudoknot==0)	
       			{
       				pse=stem_j[2][1];pseudoknot=1;
       				fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",pse/3,0,0,0,0,0,0,0,0,0);
                	  	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 
                	  	pseu(3);
       			}
      		
       		}
       		int k1,k2;
                k1=(int)(stem_length[1]/2);
                k2=(int)(stem_length[2]/2);
                fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,0,0,0,0,0,0);
       		
       	}
       	if(pseudoknot==0)
       	{
       		if(N_stem==2)
       		{
       			for(i=1;i<=N_stem;i++)
         		{ 
             			if(i==1)           
             			{   
             				j_[i]=stem_i[2][1]-stem_i[1][stem_length[i]]; 
             			}
             			else
             			{   	
             				j_[i]=stem_j[1][stem_length[1]]-stem_j[i][1]; 
             			}
             			if(j_[i]>=0&&stem_length[i]>2) 
             			{
             				k_l++;
             			}
             			else 
             			{
             				k_l=k_l;
             			}   
             			j_l_all=j_l_all+j_[i];   
         		}
       		}
       		if(N_stem>2)
       		{
         		for(i=1;i<=N_stem;i++)
         		{ 
             			if(i==1)           
             			{   
             				j_[i]=stem_i[2][1]-stem_i[1][stem_length[i]]; 
             			}
             			else if(i==N_stem) 
             			{   	
             				j_[i]=stem_j[1][stem_length[1]]-stem_j[i][1]; 
             			}
             			else               
             			{   
             				j_[i]=stem_i[i+1][1]-stem_j[i][1];            
             			}
             			if(j_[i]>=0&&stem_length[i]>=2) 
             			{
             				k_l++;
             			}
             			else 
             			{
             				k_l=k_l;
             			}   
             			j_l_all=j_l_all+j_[i];   
         		}
   		}
     		if(k_l==N_stem&&N_stem_all>10) 
      		{
      			if(k_l==2)	{junction=2;}
    			if(k_l==3)	{junction=3;}
         		if(k_l==4)	{junction=4;}
         		if(k_l==5)	{junction=5;}
   		}
        	else                    {junction=0;} 
        	if(junction==1)
        	{
        		fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",0,0,0,0,0,0,0,0,0,0);
                	fprintf(fp2,"%d %d %d %d\n",0,0,0,0);
                	int k1;
                	k1=(int)(stem_length[1]/2);
                	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,0,0,0,0,0,0,0,0);
        	}
       		if(junction==2)           
      		{ 
            		i1=stem_i[1][stem_length[1]];j1=stem_j[1][stem_length[1]];
               		i2=stem_i[2][1];j2=stem_j[2][1];
               		i3=0;j3=0;i4=0;j4=0;i5=0;j5=0;
              		si1=stem_i[2][stem_length[2]];sj1=stem_j[2][stem_length[2]];
             		si2=stem_i[3][stem_length[3]];sj2=stem_j[3][stem_length[3]];
           		{
                 		fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",i1/3,j1/3,i2/3,j2/3,i3/3,j3/3,i4/3,j4/3,i5/3,j5/3);
                      		fprintf(fp2,"%d %d %d %d\n",si1/3,sj1/3,si2/3,sj2/3);
            		}
            		int k1,k2;
                	k1=(int)(stem_length[1]/2);
                	k2=(int)(stem_length[2]/2);
                	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,0,0,0,0,0,0);	
 		}
     		else if(junction==3)           
      		{	 
            		i1=stem_i[1][stem_length[1]];j1=stem_j[1][stem_length[1]];
               		i2=stem_i[2][1];j2=stem_j[2][1];
               		i3=stem_i[3][1];j3=stem_j[3][1];
             		i4=0;j4=0;i5=0;j5=0;
              		si1=stem_i[2][stem_length[2]];sj1=stem_j[2][stem_length[2]];
             		si2=stem_i[3][stem_length[3]];sj2=stem_j[3][stem_length[3]];
           		{
                 		fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",i1/3,j1/3,i2/3,j2/3,i3/3,j3/3,i4/3,j4/3,i5/3,j5/3);
                      		fprintf(fp2,"%d %d %d %d\n",si1/3,sj1/3,si2/3,sj2/3);
            		}
            		int k1,k2,k3;
                	k1=(int)(stem_length[1]/2);
                	k2=(int)(stem_length[2]/2);
                	k3=(int)(stem_length[3]/2);
                	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,stem_i[3][k3]/3,stem_j[3][k3]/3,0,0,0,0);
 		}
     		else if(junction==4)           
     		{ 
          		i1=stem_i[1][stem_length[1]];j1=stem_j[1][stem_length[1]];
            		i2=stem_i[2][1];j2=stem_j[2][1];
             		i3=stem_i[3][1];j3=stem_j[3][1];
            		i4=stem_i[4][1];j4=stem_j[4][1];
            		i5=0;j5=0;
           		si1=stem_i[2][stem_length[2]];sj1=stem_j[2][stem_length[2]];
             		si2=stem_i[4][stem_length[4]];sj2=stem_j[4][stem_length[4]];
          		{
                 		fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",i1/3,j1/3,i2/3,j2/3,i3/3,j3/3,i4/3,j4/3,i5/3,j5/3);
                   		fprintf(fp2,"%d %d %d %d\n",si1/3,sj1/3,si2/3,sj2/3);
          		}
          		int k1,k2,k3,k4;
                	k1=(int)(stem_length[1]/2);
                	k2=(int)(stem_length[2]/2);
                	k3=(int)(stem_length[3]/2);
                	k4=(int)(stem_length[4]/2);
                	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,stem_i[3][k3]/3,stem_j[3][k3]/3,stem_i[4][k4]/3,stem_j[4][k4]/3,0,0);
    		}
    		else if(junction==5)           
   		{ 
     			i1=stem_i[1][stem_length[1]];j1=stem_j[1][stem_length[1]];
     	 		i2=stem_i[2][1];j2=stem_j[2][1];
            		i3=stem_i[3][1];j3=stem_j[3][1];
          		i4=stem_i[4][1];j4=stem_j[4][1];
           		i5=stem_i[5][1];j5=stem_j[5][1];
            		si1=stem_i[2][stem_length[2]];sj1=stem_j[2][stem_length[2]];
           		si2=stem_i[5][stem_length[5]];sj2=stem_j[5][stem_length[5]];
            		{
                    		fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",i1/3,j1/3,i2/3,j2/3,i3/3,j3/3,i4/3,j4/3,i5/3,j5/3);
                   		fprintf(fp2,"%d %d %d %d\n",si1/3,sj1/3,si2/3,sj2/3);
           		}
           		int k1,k2,k3,k4,k5;
                	k1=(int)(stem_length[1]/2);
                	k2=(int)(stem_length[2]/2);
                	k3=(int)(stem_length[3]/2);
                	k4=(int)(stem_length[4]/2);
                	k5=(int)(stem_length[5]/2);
                	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,stem_i[3][k3]/3,stem_j[3][k3]/3,stem_i[4][k4]/3,stem_j[4][k4]/3,stem_i[5][k5]/3,stem_j[5][k5]/3);
    		} 
    		else
    		{
    			fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",0,0,0,0,0,0,0,0,0,0);
                	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 		
    		}
    		
    	}  
   }
   else
   {
   	fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",0,0,0,0,0,0,0,0,0,0);
       	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 	
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
   	fprintf(fp6,"%c ",dot[i]);
   }
   fclose(fp6);	
   			                 
}

void pseu(int kkk)
{
	float BaseStacking();
	int i,j;
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
	T=300;
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
	}
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
