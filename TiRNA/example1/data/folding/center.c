#include<stdio.h>
#include<stdlib.h>
#include<math.h>
FILE *input,*output,*fc;
int t1[1000],id1[1000];
char type[1000];
float x[1000],y[1000],c[1000],z[1000],R[1000],Q[1000],f[1000];
int N;
int main()
{

      	int i1=0,j1=0;
      	void move();
      	fc=fopen("ch_0.dat","r+");
      	input=fopen("conf_all.dat","r+");
      	output=fopen("cconf_0.dat","w+"); 
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
         	move();
         	j1++;
      	}
      	fclose(input);
      	fclose(output);
      	return 0;     
}


void move()
{
        int i,j;
        float cx=0,cy=0,cz=0,acx=0,acy=0,acz=0;
        for(i=1;i<=N;i++)
        {
           cx=cx+x[i];
           cy=cy+y[i];
           cz=cz+z[i];
        }
        acx=cx/N;
        acy=cy/N;
        acz=cz/N;
        for(i=1;i<=N;i++)
        {
            fprintf(output,"%d %d %c %f %f %f %f %f %f\n",t1[i],id1[i],type[i],x[i]-acx,y[i]-acy,z[i]-acz,R[i],Q[i],f[i]);
        }
}
    
