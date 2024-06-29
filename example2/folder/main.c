#include<stdio.h>
#include<stdlib.h>
FILE *fp;
int main()
{
	fp=fopen("config1.dat","r+");
	int ca,cb,cc,cd,ce,cf;
	while(!feof(fp))
	{
		fscanf(fp,"%d %d %d %d %d %d\n",&ca,&cb,&cc,&cd,&ce,&cf);
	}
	fclose(fp);
	if(ca==1)
	{
		system("bash run_remc.sh");
	}
	else
	{
		system("bash run_sa.sh");
	}
	return 0;
}
		
