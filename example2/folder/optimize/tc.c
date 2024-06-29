/*for changing conformation data to standard PDB files*/

#include <stdio.h>
#include <math.h>
#include <string.h>

int main()
{
	FILE *infile, *outfile, *tempfile, *tempfile2;
	FILE *inna, *incl, *tempna, *tempcl;
	float x,y,z,r,q,d1;
	int cn, n, cn0, cnna, cncl, cnna0, cncl0;
	char cc,dd[4];
	int N=1;
	infile=fopen("cconf_0.dat","r+");
	tempfile=fopen("cconf_0.dat","r+");
	tempfile2=fopen("cconf_0.dat","r+");
	outfile=fopen("cf.pdb","w+");
	char btn,btn0;
	
	inna=fopen("confi.dat","r+");
	incl=fopen("confcl.dat","r+");
	tempna=fopen("confi.dat","r+");
	tempcl=fopen("confcl.dat","r+");
	int ionflag;
	ionflag=1;
	int chn;
	chn=0;
	char chnn;
	chnn='A';
	int ptstat;
	ptstat=0;
	int ishead,istail;
	ishead=0;istail=1;
	int gethead,gettail;
	/*printf("Include head P? 1.yes 0.no\n");
	scanf("%d",&gethead);
	printf("Include tail P? 1.yes 0.no\n");
	scanf("%d",&gettail);*/
	gethead=1;gettail=1;
	
	if(inna==NULL)
	{
		//printf("Without IONS.\n");
		ionflag=0;
	}	
	fscanf(tempfile2,"%d %d %c %f %f %f %f %f %f\n",&cn,&n,&btn0,&x,&y,&z,&r,&q,&d1);
	fscanf(tempfile,"%d %d %c %f %f %f %f %f %f\n",&cn0,&n,&cc,&x,&y,&z,&r,&q,&d1);
	if(ionflag)
	{
		fscanf(tempna,"%d %d %f %f %f %f %f %f\n",&cnna0,&n,&x,&y,&z,&r,&q,&d1);
		fscanf(tempcl,"%d %d %f %f %f %f %f %f\n",&cncl0,&n,&x,&y,&z,&r,&q,&d1);
	}
	
	int i=0, j=0;
	for(int k=0;k<N;k++)
	{
		j=0;
		if(i==0)
			fprintf(outfile,"CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n");
		while(1)
		{
			if(feof(infile))
			{
				fprintf(outfile,"TER\n");
				fprintf(outfile,"END");
				break;
			}
			//tempfile=infile;
			
			while(!feof(tempfile2)&&cc=='P')
			{
				fscanf(tempfile2,"%d %d %c %f %f %f %f %f %f\n",&cn,&n,&btn0,&x,&y,&z,&r,&q,&d1);
				
				if(btn0=='P')
				{
					//chnn++;
					
					if(cn!=cn0)
					{
						ptstat=2;
						ishead=-1;istail=0;
					}
						
					else
					{
						ptstat=1;
						ishead=-1;istail=0;
					}
						
					chn++;
										
					break;
				}
					
				if(btn0=='T'||btn0=='U'||btn0=='C'||btn0=='A'||btn0=='G')
				{
					fscanf(tempfile2,"%d %d %c %f %f %f %f %f %f\n",&cn,&n,&btn,&x,&y,&z,&r,&q,&d1);
					btn=btn0;
					ptstat=0;
					break;
				}
			}
			
			fscanf(infile,"%d %d %c %f %f %f %f %f %f\n",&cn,&n,&cc,&x,&y,&z,&r,&q,&d1);
			if(cc=='P')
			{sprintf(dd,"P");}
			if(cc=='S')
			{sprintf(dd,"C4'");}
			if(cc=='T'||cc=='U'||cc=='C')
			{sprintf(dd,"N1");}
			if(cc=='A'||cc=='G')
			{sprintf(dd,"N9");}
			j++;//printf("%-6s%5d  %-3s %3c %c%4d %11.3f%8.3f%8.3f\n","ATOM",n,dd,btn,"A",(i+2)/3,x,y,z);
			if(feof(infile))
				istail=0;
			if(gethead)
				ishead=1;
			if(gettail)
				istail=1;
			if(ptstat==2&&!feof(infile)&&ishead&&istail)
			fprintf(outfile,"%-6s%5d  %-3s %3c %c%4d %11.3f%8.3f%8.3f\n","ATOM",n,dd,btn,chnn+chn-ptstat+1,(n+2-chn)/3,x,y,z);
			else if(!feof(infile)&&ishead&&istail)
			fprintf(outfile,"%-6s%5d  %-3s %3c %c%4d %11.3f%8.3f%8.3f\n","ATOM",n,dd,btn,chnn+chn-ptstat,(n+2-chn)/3,x,y,z);
			if(feof(infile)&&ishead&&istail)
			fprintf(outfile,"%-6s%5d  %-3s %3c %c%4d %11.3f%8.3f%8.3f\n","ATOM",n,dd,btn,chnn+chn,(n+1-chn)/3,x,y,z);
			ishead++;istail++;
			fscanf(tempfile,"%d %d %c %f %f %f %f %f %f\n",&cn,&n,&cc,&x,&y,&z,&r,&q,&d1);
			if(cn!=cn0)
			{	
				
				fprintf(outfile,"TER\n");
				while(ionflag)
				{
					if(feof(inna))
					{
						break;
					}
				
					fscanf(inna,"%d %d %f %f %f %f %f %f\n",&cnna,&n,&x,&y,&z,&r,&q,&d1);
					fprintf(outfile,"%-6s%5d  %-3s %3s %s%4d %11.3f%8.3f%8.3f\n","ATOM",n,"A","A","B",n,x,y,z);
				
					fscanf(tempna,"%d %d %f %f %f %f %f %f\n",&cnna,&n,&x,&y,&z,&r,&q,&d1);
					if(cnna!=cnna0)
					{
						fprintf(outfile,"TER\n");
						while(1)
						{
							if(feof(incl))
							{
								break;
							}
							fscanf(incl,"%d %d %f %f %f %f %f %f\n",&cncl,&n,&x,&y,&z,&r,&q,&d1);
							fprintf(outfile,"%-6s%5d  %-3s %3s %s%4d %11.3f%8.3f%8.3f\n","ATOM",n,"B","A","C",n,x,y,z);
						
							fscanf(tempcl,"%d %d %f %f %f %f %f %f\n",&cncl,&n,&x,&y,&z,&r,&q,&d1);
							if(cncl!=cncl0)
							{
								fprintf(outfile,"TER\n");
								cncl0=cncl;
								break;
							}
					
						}
						cnna0=cnna;
						break;
					}
					
				}
				fprintf(outfile,"END\n");
				fprintf(outfile,"CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n");
				chn=0;chnn='A';
				cn0=cn;
				i++;
				j=0;
				N++;
				break;
			}
		}
	}
	
	return 0;
}
