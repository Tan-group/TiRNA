#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<ctype.h>
#define MAX_LENGTH 10000
char rna_sequence[MAX_LENGTH];
char dot_bracket[MAX_LENGTH];
int N,se[1000][1000],c[1000];
FILE *fp,*fp1,*fp2,*fp3;
int main()
{
	void secondary();
	fp=fopen("seq.dat","r+");
	fp1=fopen("sequence.dat","w+");
	fp2=fopen("cs.dat","w+");
	fp3=fopen("cs_rebuild.dat","w+");

    	// Read RNA sequence
    	if (fgets(rna_sequence, MAX_LENGTH, fp) == NULL) 
    	{
        	perror("Error reading RNA sequence");
        	fclose(fp);
        	return EXIT_FAILURE;
    	}
    	for (int i = 0; rna_sequence[i] != '\0'; i++) 
    	{
        	rna_sequence[i] = toupper(rna_sequence[i]);
    	}
    	//printf("%s",rna_sequence);
	fprintf(fp1,"%s",rna_sequence);
	N=strlen(rna_sequence)-1;
	//printf("rna_length %d\n",N);
    	// Read dot bracket notation
    	if (fgets(dot_bracket, MAX_LENGTH, fp) == NULL) 
    	{
        	perror("Error reading dot bracket notation");
        	fclose(fp);
        	return EXIT_FAILURE;
    	}

	secondary();
	fclose(fp1);fclose(fp2);fclose(fp3);
	return 0;
}


void secondary()
{
 	char stack[1000],stack1[1000],stack2[1000];
   	int xx[1000],kk,kk1,kk2,xx1[1000],xx2[1000];
   	int top=-1,top1=-1,top2=-1;
   	int index=0;
   	int ai,aj,i,j,sk;
   	int tertiary[1000];
   	for(i=0;i<=N;i=i+3) 
   	{ 
      		c[i]=0;
      		tertiary[i]=0;
      		for(sk=i+9;sk<=N;sk=sk+3)
      		{
          		se[i][sk]=0;
      		}
   	}
   	while(dot_bracket[index]!=0)
   	{
          	if(dot_bracket[index]=='(')
          	{
               		kk=top;
               		stack[++top]=dot_bracket[index];
               		xx[++kk]=index; 
          	}
          	if(dot_bracket[index]==')')
          	{
                 	if(dot_bracket[index]==')'&&stack[top]=='(')
                 	{
                    		//fprintf(fp2,"%d %d\n",xx[top]+1,index+1);
                     		ai=xx[top]+1;aj=index+1;
                     		se[ai][aj]=1;c[ai]=1;c[aj]=1;
				tertiary[ai]=0;tertiary[aj]=0;
                   //  printf("%d %d\n",xx[top]+1,index+1);
                  	}
                  	top--;
          	}
         	if(dot_bracket[index]=='[')
          	{
               		kk1=top1;
               		stack1[++top1]=dot_bracket[index];
               		xx1[++kk1]=index;               
          	}
          	if(dot_bracket[index]==']')
          	{
                 	if(dot_bracket[index]==']'&&stack1[top1]=='[')
                 	{
                    		//fprintf(fp2,"%d %d\n",xx1[top1]+1,index+1);
                    		ai=xx1[top1]+1;aj=index+1;
                    		se[ai][aj]=1;c[ai]=1;c[aj]=1;tertiary[ai]=1;tertiary[aj]=1;
                   		// printf("%d %d\n",xx1[top1]+1,index+1);
                  	}
                  	top1--;
           	}
           	if(dot_bracket[index]=='{')
          	{
               		kk2=top2;
               		stack2[++top2]=dot_bracket[index];
               		xx2[++kk2]=index; 
          	}
          	if(dot_bracket[index]=='}')
          	{
                 	if(dot_bracket[index]=='}'&&stack2[top2]=='{')
                 	{
                    		//fprintf(fp2,"%d %d\n",xx[top]+1,index+1);
                     		ai=xx2[top2]+1;aj=index+1;
                     		se[ai][aj]=1;c[ai]=1;c[aj]=1;tertiary[ai]=1;tertiary[aj]=1;
                   		//  printf("%d %d\n",xx[top]+1,index+1); 
                  	}
                  	top2--;
          	} 
          	index++;
       	}
       	for(i=1;i<N;i++)
       	{
       		for(j=i+1;j<=N;j++)
       		{
       			if(se[i][j]==1)
       			{
       				fprintf(fp2,"%d %c %d %c %d\n",i,rna_sequence[i-1],j,rna_sequence[j-1],tertiary[i]);
       			}
       		}
       	}
       	for(i=1;i<N;i++)
       	{
       		for(j=i+1;j<=N;j++)
       		{
       			if(se[i][j]==1)
       			{
       				fprintf(fp3,"%d %c %d %c\n",i,rna_sequence[i-1],j,rna_sequence[j-1]);
       			}
       		}
       	}
}

