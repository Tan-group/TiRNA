#include<stdio.h>
#include<stdlib.h>
#include<math.h>
FILE *input,*output,*fc,*fpt;
int t1[1000],id1[1000],N_conf;
char type[1000];
float x[1000],y[1000],c[1000],z[1000],R[1000],Q[1000],f[1000];
int N;
int main()
{

      	int i1=0,j1=0;
      	fc=fopen("ch.dat","r+");
      	input=fopen("conf_all.dat","r+");
      	output=fopen("ch_0.dat","w+"); 
      	fpt=fopen("top1.dat","r+");
      	int i=0;
      	while(!feof(fc))
      	{
           	fscanf(fc,"%d %d %s %f %f %f %f %f %f\n",&t1[i],&id1[i],&type[i],&x[i],&y[i],&z[i],&R[i],&Q[i],&f[i]);
           	i++;
      	}
      	fclose(fc);
      	N=i;
      	while(!feof(input))
      	{
         	for(i1=j1*N+1;i1<=j1*N+N;i1++)
         	{
               		fscanf(input,"%d %d %s %f %f %f %f %f %f\n",&t1[i1-j1*N],&id1[i1-j1*N],&type[i1-j1*N],&x[i1-j1*N],&y[i1-j1*N],&z[i1-j1*N],&R[i1-j1*N],&Q[i1-j1*N],&f[i1-j1*N]);
         	} // 读入构象；i1: No. of nt; j1: No. of conf.;
         	j1++;
  	}
	fclose(input);
	while(!feof(fpt))
	{
		fscanf(fpt,"%d\n",&N_conf);
	}
	fclose(fpt);
  	input=fopen("conf_all.dat","r+");
  	j1=0;i1=0;
  	while(!feof(input))
      	{
         	for(i1=j1*N+1;i1<=j1*N+N;i1++)
         	{
               		fscanf(input,"%d %d %s %f %f %f %f %f %f\n",&t1[i1-j1*N],&id1[i1-j1*N],&type[i1-j1*N],&x[i1-j1*N],&y[i1-j1*N],&z[i1-j1*N],&R[i1-j1*N],&Q[i1-j1*N],&f[i1-j1*N]);
         	} // 读入构象；i1: No. of nt; j1: No. of conf.;
       		if(j1==N_conf)
      		{
       			for(i1=j1*N+1;i1<=j1*N+N;i1++)
      			{
       		 		fprintf(output,"%d %d %c %f %f %f %f %f %f\n",t1[i1-j1*N],id1[i1-j1*N],type[i1-j1*N],x[i1-j1*N],y[i1-j1*N],z[i1-j1*N],R[i1-j1*N],Q[i1-j1*N],f[i1-j1*N]);
      			}
      		}
      		j1++;
  	}
      	fclose(input);
      	fclose(output);
      	return 0;     
}
