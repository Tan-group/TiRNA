/************* 该程序是用来调用不同的wham，因为hairpin和其他结构对反映坐标的定义不同，所以需要调用不同的程序 ******************/
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
		system("bash tmhairpin.sh");
	}
	else
	{
		system("bash tm.sh");
	}
	return 0;
}
