#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#define MAX_LENGTH 10000
char rna_sequence[10000];
int s[1000][1000],c[1000],N0;
int tertiary[1000];
FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5;
char type[1000];
int main()
{
    	int i,ai,aj,N,j;
    	char ax,ay;
    	void secondary();
    	fp=fopen("seq.dat","r+");
    	fp1=fopen("stem.dat","w+");
   	fp2=fopen("stem_kissing.dat","w+");
    	fp3=fopen("cs.dat","r+");
    	fp4=fopen("cs1.dat","w+");
    	fp5=fopen("state.dat","w+");
    	if (fgets(rna_sequence, MAX_LENGTH, fp) == NULL) 
    	{
        	perror("Error reading RNA sequence");
        	fclose(fp);
        	return EXIT_FAILURE;
  	}
	N=strlen(rna_sequence)-1;
	N0=N*3+1;
	//printf("%d %d\n",N0,N);
	for(i=1;i<N0;i++)
	{
		for(j=i+1;j<=N0;j++)
		{
			s[i][j]=0;
			c[i]=0;c[j]=0;
			tertiary[i]=0;
			tertiary[j]=0;
		}
		
	}
	int ad;
    	while(!feof(fp3))
    	{
    		fscanf(fp3,"%d %c %d %c %d\n",&ai,&ax,&aj,&ay,&ad);
    		s[ai*3][aj*3]=1;c[ai*3]=1;c[aj*3]=1;tertiary[ai*3]=ad;tertiary[aj*3]=ad;
    		type[ai*3]=ax;type[aj*3]=ay;
    	}
     	secondary();      
     	fclose(fp1);fclose(fp2);fclose(fp3);fclose(fp4);fclose(fp5);
} 


void secondary()
{
   	int i,j,k1,k2,stem_i[20][20],stem_j[20][20],order=1,k3,N_stem,kk,j_[20],k_l=0,N_stem_all,stem_length[20],j_l_all=0,kN;
   	int i1,j1,i2,j2,i3,j3,i4,j4,i5,j5;
   	int si1,si2,sj1,sj2;
   	int junction,pseudoknot=0,pse_base[1000];
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
/***************初次判断二级结构，是为了区分3-way以上且含有三级相互作用，目的是为了更好的判断junction的way数*************************/    	    
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
                    						}
                    						else	{order++;}
                             				} 
                             			}             
                 			
                 				break;                              
                			}
              			}
           		}
		}
       	}
      	N_stem=order;  
      	//printf("----------- %d\n",N_stem);
/**************************如果stem小于4那么这个结构就肯定不是3-way且带有三级相互作用**************************/  
      	if(N_stem<4)
      	{
      		for(i=1;i<N0;i++)
      		{
      			for(j=i+1;j<=N0;j++)
      			{
      				if(s[i][j]==1)
      				{
      					fprintf(fp4,"%d %c %d %c %d\n",i/3,type[i],j/3,type[j],0);
      				}
      			}
      		}
      	}
      	else
      	{
      	
      		for(i=1;i<N0;i++)
      		{
      			for(j=i+1;j<=N0;j++)
      			{
      				if(s[i][j]==1)
      				{
      					fprintf(fp4,"%d %c %d %c %d\n",i/3,type[i],j/3,type[j],tertiary[i]);
      				}
      			}
      		}
      	}
       	order=1;j_l_all=0; N_stem_all=0;
 /****************************再次判断二级结构，目的是去除三级相互作用，这样就是纯的junction判断方便****************************/
       	if(N_stem>=4)
       	{
       	
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
      			if (fmod(i,3)==0&&c[i]==1&&tertiary[i]==0)  
      			{
         			for (j=i+12;j<=N0;j++)  
         			{
           				if (fmod(j,3)==0&&c[j]==1&&tertiary[j]==0)  
           				{
             					if (s[i][j]==1)  
             					{                            
								if(i>6)
								{
									if((s[i+3][j-3]==1&&s[i-3][j+3]==1)||(s[i-3][j+3]==1&&s[i+3][j-3]==0))
                             						{
                             							order=order; 
                             							stem_length[order]=stem_length[order]+1;
                             							k3=stem_length[order];
                             							stem_i[order][k3]=i;  
                             							stem_j[order][k3]=j;
                             							printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             								printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             									printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             							printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
                             								printf("%d %d %d %d\n",order,k3,stem_i[order][k3]/3,stem_j[order][k3]/3);
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
       		   	
       	}
    if(N_stem<=5)
    {   	
       	
      	//printf("*****%d\n",N_stem); 
      	for(i=1;i<=N_stem;i++)
       	{
          	N_stem_all=stem_length[i]+N_stem_all;
          	
       	}       	
       	if(N_stem==1)	
       	{
       		junction=1;
              //  printf("RNA is hairpin\n");
       	}
      	int pse=0,k;
      	pseudoknot=0;
/*******************目前就只能处理H-type假结*************************************/
      	if(N_stem>1&&N_stem<4)
       	{
       		for(i=1;i<=N0;i++)
       		{
       			pse_base[i]=0;
       		}
       		if(N_stem==2)
       		{
       			
       			j_[1]=(stem_j[1][1]-stem_j[2][stem_length[2]])/3;
       			j_[2]=(stem_j[1][1]-stem_i[2][stem_length[2]])/3;
       			//printf("------- %d %d\n",stem_i[2][1]/3,stem_j[1][1]/3);
       			if(j_[1]<0&&j_[2]>0)	
       			{
       				pse=stem_j[1][1];pseudoknot=1;
       				fprintf(fp1,"%d %d %d %d %d %d %d %d %d %d\n",pse/3,0,0,0,0,0,0,0,0,0);
                		fprintf(fp2,"%d %d %d %d\n",0,0,0,0);
       				//printf("RNA is H-type pseudoknot\n");
       			}
       			else
       			{
       				pseudoknot=0;pse=0;
       			}

                	
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
    if(N_stem==1)
   {
   	int k1;
       	k1=(int)(stem_length[1]/2);
       	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,0,0,0,0,0,0,0,0);
   } 
   if(N_stem==2)
   {
   	int k1,k2;
        k1=(int)(stem_length[1]/2);
       	k2=(int)(stem_length[2]/2);
    	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,0,0,0,0,0,0);	
   }
   if(N_stem==3)
   {
   	int k1,k2,k3;
     	k1=(int)(stem_length[1]/2);
     	k2=(int)(stem_length[2]/2);
       	k3=(int)(stem_length[3]/2);
     	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,stem_i[3][k3]/3,stem_j[3][k3]/3,0,0,0,0);
   }
   if(N_stem==4)
   {
   	int k1,k2,k3,k4;
      	k1=(int)(stem_length[1]/2);
     	k2=(int)(stem_length[2]/2);
      	k3=(int)(stem_length[3]/2);
      	k4=(int)(stem_length[4]/2);
      	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,stem_i[3][k3]/3,stem_j[3][k3]/3,stem_i[4][k4]/3,stem_j[4][k4]/3,0,0);
   }
   if(N_stem==5)
   {
   	int k1,k2,k3,k4,k5;
       	k1=(int)(stem_length[1]/2);
    	k2=(int)(stem_length[2]/2);
     	k3=(int)(stem_length[3]/2);
      	k4=(int)(stem_length[4]/2);
       	k5=(int)(stem_length[5]/2);
     	fprintf(fp5,"%d %d %d %d %d %d %d %d %d %d\n",stem_i[1][k1]/3,stem_j[1][k1]/3,stem_i[2][k2]/3,stem_j[2][k2]/3,stem_i[3][k3]/3,stem_j[3][k3]/3,stem_i[4][k4]/3,stem_j[4][k4]/3,stem_i[5][k5]/3,stem_j[5][k5]/3);
   }                                
}
