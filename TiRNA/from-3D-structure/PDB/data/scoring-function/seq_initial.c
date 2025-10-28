#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//int n1 =0, n2 = 0;
using namespace std;
int main()
{
	FILE *fp1, *fp2, *fp3;
	fp1=fopen("seq.dat","r+");
	fp3=fopen("ch.dat","w+");
	
	char seqs[5][200];
	int chn=0;
	 
	while(!feof(fp1))
	{
		fscanf(fp1,"%s\n",seqs[chn]);
		//printf("%s\n",seqs[chn]);
		chn++;
	}
	printf("%d\n",chn);
	int a;
	int n1[1000], n2[1000],j;
	char type;
	float x, y, z, R, Q, f;
	for(int i=0;i<chn;i++)
	{
		fp2=fopen("initial.dat","r+");
		for(j=0;j<((int)strlen(seqs[i])*3+1);j++)
		{
		       
			fscanf(fp2,"%d %d %s %f %f %f %f %f %f\n",&n1[j],&n2[j],&type,&x,&y,&z,&R,&Q,&f);
			printf("%d %d %c %f %f %f %f %f %f\n",n1[j],n2[j],type,x,y,z,R,Q,f);
			
			if(fmod(j+1,3)==0)
			{
				fprintf(fp3,"%d %d %c %f %f %f %f %f %f\n",n1[j],n2[j],seqs[i][j/3],x,y,z,R,Q,f);
			}
			else
			{
				fprintf(fp3,"%d %d %c %f %f %f %f %f %f\n",n1[j],n2[j],type,x,y,z,R,Q,f);
			}
		}
                printf("j %d\n",j-1);
		fclose(fp2);
	}
	fclose(fp1);fclose(fp3);	
	return 0;
}
