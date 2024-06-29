#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define N
FILE *fp,*fp1,*fpcfig;
char type[50000][15];
float a[100000];
int t[100000],c;
int main()
{
       int i=1,j,n;
       float b;
       fp=fopen("cgRNASP.txt","r+");
       fp1=fopen("tiao.sh","w+");
       fpcfig=fopen("config1.dat","r+");
       int ca,cb,cc,cd,ce,Ncout,ct;
       while(!feof(fpcfig))
       {
       		fscanf(fpcfig,"%d %d %d %d %d %d\n",&ct,&ca,&cb,&cc,&cd,&ce);
       }
       Ncout=ce;
       while(!feof(fp))
       {
           fscanf(fp,"%s %f\n",type[i],&a[i]);
           i++;
       }
       fclose(fp);
     //  printf("%f\n",a[1]);
       n=i-1;
       for(i=1;i<=n;i++)
       {
             t[i]=i;
       }
       for(i=1;i<n;i++)
       {
          for(j=i+1;j<=n;j++)
          {
                 if(a[j]<=a[i])
                 {
                     b=a[i];
                     a[i]=a[j];
                     a[j]=b;
                     c=t[i];
                     t[i]=t[j];
                     t[j]=c;
                 }
           }
        }
        for(i=1;i<=Ncout;i++)
        {
        	fprintf(fp1,"cp %s ../\n",type[t[i]]);
        	fprintf(fp1,"cd ../\n");
        	fprintf(fp1,"mv cf*.pdb CG_top%d.pdb\n",i);
        	fprintf(fp1,"cd xx\n");
        }
        fclose(fp1);
        
         //   printf("%s\n",type[1]);
      return 0;
}
