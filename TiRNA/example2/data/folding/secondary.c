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
FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5;
int main()
{
    	int i;
    	void secondary(),BasePairing();
    	fp=fopen("top1.pdb","r+");
    	fp1=fopen("stem.dat","w+");
   	fp2=fopen("stem_kissing.dat","w+");
    	fp3=fopen("cs.dat","w+");
    	fp4=fopen("RNA_type","w+");
    	fp5=fopen("state.dat","w+");
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
       /* for(i=1;i<N0;i++)
        {
        	for(j=i+1;j<=N0;j++)
        	{
        		if(s[i][j]==1)
        		{
        			fprintf(fp3,"%d %c %d %c\n",i/3,type[i],j/3,type[j]);
        		}
        	}
        }*/
        	
}

void secondary()
{
   	int i,j,k1,k2,stem_i[20][20],stem_j[20][20],order=1,k3,N_stem,kk,j_[20],k_l=0,N_stem_all,stem_length[20],j_l_all=0;
   	int i1,j1,i2,j2,i3,j3,i4,j4,i5,j5;
   	int si1,si2,sj1,sj2;
   	int junction,pseudoknot=0,kN;
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
                             							//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             								//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             									//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             							//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             								//printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
    if(N_stem<=5)
    {
      	for(i=1;i<=N_stem;i++)
       	{
          	N_stem_all=stem_length[i]+N_stem_all;
          	
       	}       	
       	if(N_stem==1)	
       	{
       		junction=1;
               // printf("RNA is hairpin\n");
                fprintf(fp4,"%d\n",1);
                int k1;
                k1=(int)(stem_length[1]/2);
                fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,0,0,0,0,0,0,0,0);
       	}
       	else
       	{
       		fprintf(fp4,"%d\n",2);
       	}
       	fclose(fp4);
      	int pse=0;
      	pseudoknot=0;
      	if(N_stem>1&&N_stem<4)
       	{
       		if(N_stem==2)
       		{
       			j_[1]=(stem_j[1][1]-stem_j[2][stem_length[2]])/3;
       			j_[2]=(stem_j[1][1]-stem_i[2][stem_length[2]])/3;
       			//printf("------- %d %d\n",stem_i[2][1]/3,stem_j[1][1]/3);
       			if(j_[1]<0&&j_[2]>0)	
       			{
       				pse=stem_j[1][1];pseudoknot=1;
       				//printf("RNA is H-type pseudoknot\n");
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
       				//printf("RNA is H-type pseudoknot\n");
       			}
       				//printf("RNA is H-type pseudoknot\n");}
       			j_[3]=(stem_j[1][1]-stem_j[3][stem_length[3]])/3;
       			j_[4]=(stem_j[1][1]-stem_i[3][stem_length[3]])/3;
       			if(j_[3]<0&&j_[4]>0&&pseudoknot==0)	
       			{
       				pse=stem_j[1][1];pseudoknot=1;
       				fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",pse/3,0,0,0,0,0,0,0,0,0);
                	  	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 
       				//printf("RNA is H-type pseudoknot\n");
       			}
       			j_[5]=(stem_j[2][1]-stem_j[3][stem_length[3]])/3;
       			j_[6]=(stem_j[2][1]-stem_i[3][stem_length[3]])/3;
       			if(j_[5]<0&&j_[6]>0&&pseudoknot==0)	
       			{
       				pse=stem_j[2][1];pseudoknot=1;
       				//printf("RNA is H-type pseudoknot\n");
       				fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",pse/3,0,0,0,0,0,0,0,0,0);
                	  	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 
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
        	}
       		else if(junction==2)           
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
    	//printf("RNA is %d-way juncion\n",junction);  
   }
   else
   {
   	fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",0,0,0,0,0,0,0,0,0,0);
       	fprintf(fp2,"%d %d %d %d\n",0,0,0,0); 	
   }                  
}
