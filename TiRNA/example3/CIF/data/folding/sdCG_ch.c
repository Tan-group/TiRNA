#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
FILE *inpdb, *outpdb;
#define Rp 2.11
#define Rc 1.8
#define Rn 1.53
#define Qp -1.0
#define Qc 0.0
#define Qn 0.0
#define ss 0.0

int main()
{
	char atom[20],typen[100],typeb,chainC;
	int id0,idn;
	float x,y,z;
	
	int i;
	
	char str1[10],str2[10],str3[10],str4[10],str5[10],str6[10],str7[10],str8[10],str9[10],str10[10],str11[10],str12[10];
 	sprintf(str1,"ATOM");sprintf(str2,"C4'");sprintf(str3,"P");sprintf(str4,"N1");sprintf(str5,"N9");sprintf(str6,"A");sprintf(str7,"T");sprintf(str8,"U");sprintf(str9,"C");sprintf(str10,"G");sprintf(str11,"END");sprintf(str12,"TER");
 	char fileN[500][50],infileN[100][50],temp[5];
	inpdb=fopen("top1.pdb","r+");
	outpdb=fopen("ch.dat","w+");
 		
 		i=1;
 		while(strcmp(atom,str11)&&(!feof(inpdb)))
 		{
 			fscanf(inpdb,"%s %d %s %c %c %d %f %f %f\n",atom,&id0,typen,&typeb,&chainC,&idn,&x,&y,&z);
 			if(atom[0]=='T')
 			{
 				i=-1;
 				continue;
 			}
 			if(!strcmp(typen,str2))
 				{fprintf(outpdb,"1 %d S %f %f %f %f %f %f\n",i,x,y,z,Rc,Qc,ss);}
 			if(!strcmp(typen,str4)||!strcmp(typen,str5))
 				{fprintf(outpdb,"1 %d %c %f %f %f %f %f %f\n",i,typeb,x,y,z,Rn,Qc,ss);}
 			if(!strcmp(typen,str3))
 				{fprintf(outpdb,"1 %d P %f %f %f %f %f %f\n",i,x,y,z,Rp,Qp,ss);}
 			i++;
 		}
 		fclose(outpdb);
		
		fclose(inpdb);
	
}




