#include<stdio.h>
#include<stdlib.h>
FILE *fp;
int main()
{
	int a1;
	fp=fopen("RNA_type","r+");
	while(!feof(fp))
	{
		fscanf(fp,"%d\n",&a1);
	}
	fclose(fp);
	if(a1==1)
	{
		system("g++ sdCG_ch.c");
		system("./a.out");		
	}
	return 0;
}
