#include<stdio.h>
#include<stdlib.h>
#include<math.h>
FILE *input,*output,*fc;
int t1[1000],id1[1000];
char type[1000];
float x[1000],y[1000],c[1000],z[1000],R[1000],Q[1000],f[1000];
int N0;
float PC(),CP(),CN(),PCP(),CPC(),PCN(),NCP(),PCPC(),CPCP(),CPCN(),NCPC();
int main()
{

      	int i1=0,j1=0;
      	float Bond_energy();
      	float u_conf[10000];
      	fc=fopen("ch.dat","r+");
      	input=fopen("conf_all.dat","r+");
      	output=fopen("top1.dat","w+"); 
      	int i=0;
      	while(!feof(fc))
      	{
           	fscanf(fc,"%d %d %s %f %f %f %f %f %f\n",&t1[i],&id1[i],&type[i],&x[i],&y[i],&z[i],&R[i],&Q[i],&f[i]);
           	i++;
      	}
      	fclose(fc);
      	N0=i;
      	while(!feof(input))
      	{
         	for(i1=j1*N0+1;i1<=j1*N0+N0;i1++)
         	{
               		fscanf(input,"%d %d %s %f %f %f %f %f %f\n",&t1[i1-j1*N0],&id1[i1-j1*N0],&type[i1-j1*N0],&x[i1-j1*N0],&y[i1-j1*N0],&z[i1-j1*N0],&R[i1-j1*N0],&Q[i1-j1*N0],&f[i1-j1*N0]);
         	} // 读入构象；i1: No. of nt; j1: No. of conf.;
         	u_conf[j1]=Bond_energy();
         	//printf("%d %f\n",j1,u_conf[j1]);
         	j1++;
      }
      float umin;
      int nmin=0;
      umin=u_conf[0];
      for(i=1;i<j1-1;i++)
      {
      		if(u_conf[i]<=umin)
      		{
      			umin=u_conf[i];
      			nmin=i;
      		}
      }
      fprintf(output,"%d\n",nmin);
      fclose(input);
      fclose(output);
      return 0;     
}


float Bond_energy()
{
        int i;
        float udihedral=0.0,u=0.0,uangle=0,ubond=0.0;
        for (i=1;i<=N0;i++)
        {
             	if (fmod(i,3)==0)
             	{
                     	if (i==3) 
                     	{
                             	ubond=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z);
                             	uangle=PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z);
                              	udihedral=PCPC(i,x,y,z)+CPCP(i,x,y,z)+NCPC(i,x,y,z);
                      	}
                     	else if (i==(N0-1))
                    	{
                         	ubond=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z);
                         	uangle=PCP(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z);
                         	udihedral=CPCN(i,x,y,z);
                      	} 
                     	else 
                     	{
                             	ubond=PC(i,x,y,z)+CP(i,x,y,z)+CN(i,x,y,z);
                             	uangle=PCP(i,x,y,z)+CPC(i,x,y,z)+PCN(i,x,y,z)+NCP(i,x,y,z);
                               	udihedral=PCPC(i,x,y,z)+CPCP(i,x,y,z)+CPCN(i,x,y,z)+NCPC(i,x,y,z);
                      	}
                      	u=u+udihedral+uangle+ubond;
            	}
    	} 
    	return u;
}
  
float PC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
	float d,ul;
     	d=sqrt((x1[i1-2]-x1[i1-1])*(x1[i1-2]-x1[i1-1])+(y1[i1-2]-y1[i1-1])*(y1[i1-2]-y1[i1-1])+(z1[i1-2]-z1[i1-1])*(z1[i1-2]-z1[i1-1]));
     	ul=68*(d-3.8)*(d-3.8);
     	return ul;
}
float CP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d,ul;
     	d=sqrt((x1[i1+1]-x1[i1-1])*(x1[i1+1]-x1[i1-1])+(y1[i1+1]-y1[i1-1])*(y1[i1+1]-y1[i1-1])+(z1[i1+1]-z1[i1-1])*(z1[i1+1]-z1[i1-1]));
     	ul=33*(d-3.8)*(d-3.8);
     	return ul;
}
float CN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d,ul;
     	d=sqrt((x1[i1]-x1[i1-1])*(x1[i1]-x1[i1-1])+(y1[i1]-y1[i1-1])*(y1[i1]-y1[i1-1])+(z1[i1]-z1[i1-1])*(z1[i1]-z1[i1-1]));
     	ul=24.8*(d-3.35)*(d-3.35);
     	return ul;
}
float PCP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
      	float d1,d2,d3,w,a1,ue0;
      	d1=sqrt((x1[i1-2]-x1[i1-1])*(x1[i1-2]-x1[i1-1])+(y1[i1-2]-y1[i1-1])*(y1[i1-2]-y1[i1-1])+(z1[i1-2]-z1[i1-1])*(z1[i1-2]-z1[i1-1]));
      	d2=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
      	d3=sqrt((x1[i1+1]-x1[i1-2])*(x1[i1+1]-x1[i1-2])+(y1[i1+1]-y1[i1-2])*(y1[i1+1]-y1[i1-2])+(z1[i1+1]-z1[i1-2])*(z1[i1+1]-z1[i1-2]));
      	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
      	if (w<=-1.0) 
      	{
          	a1=3.14;
      	}
      	else if (w>=1.0) 
      	{
           	a1=0.;
      	}
      	else  
      	{
           	a1=acos(w);
      	} 
      	ue0=9.3*(a1-1.6)*(a1-1.6);
      	return ue0;
}
float CPC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
      	float d1,d2,d3,w,a1,ue0;
      	if(i1+2>N0)  
      	{
            	ue0=0.0;
      	}
      	else
      	{
           	d1=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
           	d2=sqrt((x1[i1+1]-x1[i1+2])*(x1[i1+1]-x1[i1+2])+(y1[i1+1]-y1[i1+2])*(y1[i1+1]-y1[i1+2])+(z1[i1+1]-z1[i1+2])*(z1[i1+1]-z1[i1+2]));
           	d3=sqrt((x1[i1+2]-x1[i1-1])*(x1[i1+2]-x1[i1-1])+(y1[i1+2]-y1[i1-1])*(y1[i1+2]-y1[i1-1])+(z1[i1+2]-z1[i1-1])*(z1[i1+2]-z1[i1-1]));
           	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
           	if (w<=-1.0) 
           	{
              		a1=3.14;
           	}
           	else if (w>=1.0) 
           	{
              		a1=0.;
           	}
           	else  
           	{
             		a1=acos(w);
           	}
           	ue0=15.3*(a1-1.6)*(a1-1.6);
      	}
      	return ue0;
}
float PCN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
      	float d1,d2,d3,w,a1,ue0;
      	d1=sqrt((x1[i1-1]-x1[i1-2])*(x1[i1-1]-x1[i1-2])+(y1[i1-1]-y1[i1-2])*(y1[i1-1]-y1[i1-2])+(z1[i1-1]-z1[i1-2])*(z1[i1-1]-z1[i1-2]));
      	d2=sqrt((x1[i1-1]-x1[i1])*(x1[i1-1]-x1[i1])+(y1[i1-1]-y1[i1])*(y1[i1-1]-y1[i1])+(z1[i1-1]-z1[i1])*(z1[i1-1]-z1[i1]));
      	d3=sqrt((x1[i1-2]-x1[i1])*(x1[i1-2]-x1[i1])+(y1[i1-2]-y1[i1])*(y1[i1-2]-y1[i1])+(z1[i1-2]-z1[i1])*(z1[i1-2]-z1[i1]));
      	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
      	if (w<=-1.0) {a1=3.14;}
      	else if (w>=1.0) {a1=0.;}
      	else  {a1=acos(w);}
      	ue0=9.75*(a1-1.64)*(a1-1.64);
      	return ue0;
}
float NCP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
     	float d1,d2,d3,w,a1,ue0;
     	d1=sqrt((x1[i1-1]-x1[i1])*(x1[i1-1]-x1[i1])+(y1[i1-1]-y1[i1])*(y1[i1-1]-y1[i1])+(z1[i1-1]-z1[i1])*(z1[i1-1]-z1[i1]));
     	d2=sqrt((x1[i1-1]-x1[i1+1])*(x1[i1-1]-x1[i1+1])+(y1[i1-1]-y1[i1+1])*(y1[i1-1]-y1[i1+1])+(z1[i1-1]-z1[i1+1])*(z1[i1-1]-z1[i1+1]));
     	d3=sqrt((x1[i1+1]-x1[i1])*(x1[i1+1]-x1[i1])+(y1[i1+1]-y1[i1])*(y1[i1+1]-y1[i1])+(z1[i1+1]-z1[i1])*(z1[i1+1]-z1[i1]));
     	w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
     	if (w<=-1.0) {a1=3.14;}
     	else if (w>=1.0) {a1=0.;}
     	else  {a1=acos(w);}
     	ue0=15.23*(a1-1.66)*(a1-1.66);
     	return ue0;
}
float PCPC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
 	if(i1+2>N0) {ud0=0.0;}
 	else
 	{
 		c1=((y1[i1-2]-y1[i1-1])*(z1[i1-1]-z1[i1+1])-(z1[i1-2]-z1[i1-1])*(y1[i1-1]-y1[i1+1]));
 		c2=((z1[i1-2]-z1[i1-1])*(x1[i1-1]-x1[i1+1])-(x1[i1-2]-x1[i1-1])*(z1[i1-1]-z1[i1+1]));
 		c3=((x1[i1-2]-x1[i1-1])*(y1[i1-1]-y1[i1+1])-(y1[i1-2]-y1[i1-1])*(x1[i1-1]-x1[i1+1]));
 		p1=((y1[i1-1]-y1[i1+1])*(z1[i1+1]-z1[i1+2])-(z1[i1-1]-z1[i1+1])*(y1[i1+1]-y1[i1+2]));
 		p2=((z1[i1-1]-z1[i1+1])*(x1[i1+1]-x1[i1+2])-(x1[i1-1]-x1[i1+1])*(z1[i1+1]-z1[i1+2]));
 		p3=((x1[i1-1]-x1[i1+1])*(y1[i1+1]-y1[i1+2])-(y1[i1-1]-y1[i1+1])*(x1[i1+1]-x1[i1+2]));
 		e1=sqrt(c1*c1+c2*c2+c3*c3); 
 		f1=sqrt(p1*p1+p2*p2+p3*p3);
 		pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 		g1=(x1[i1-2]-x1[i1+2]); 
 		g2=(y1[i1-2]-y1[i1+2]); 
 		g3=(z1[i1-2]-z1[i1+2]);
 		gg1=sqrt(g1*g1+g2*g2+g3*g3); 
 		hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 		if (pp1<=-1.0) {di=-3.14;}
 		else if (pp1>=1.0) {di=0.;}
 		else if (hh1>=0.) {di=acos(pp1);}
 		else {di=-acos(pp1);}
 		ud0=4.0*((1-cos(di-2.51))+0.5*(1-cos(3.*(di-2.51))));
 	}
 	return ud0;
}
float CPCP(int i1,float x1[1000],float y1[1000],float z1[1000])
{
   	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
   	if(i1+2>N0) 
   	{
           	ud0=0.0;
   	}
   	else
   	{
     		c1=((y1[i1-1]-y1[i1+1])*(z1[i1+1]-z1[i1+2])-(z1[i1-1]-z1[i1+1])*(y1[i1+1]-y1[i1+2]));
     		c2=((z1[i1-1]-z1[i1+1])*(x1[i1+1]-x1[i1+2])-(x1[i1-1]-x1[i1+1])*(z1[i1+1]-z1[i1+2]));
     		c3=((x1[i1-1]-x1[i1+1])*(y1[i1+1]-y1[i1+2])-(y1[i1-1]-y1[i1+1])*(x1[i1+1]-x1[i1+2]));
     		p1=((y1[i1+1]-y1[i1+2])*(z1[i1+2]-z1[i1+4])-(z1[i1+1]-z1[i1+2])*(y1[i1+2]-y1[i1+4]));
     		p2=((z1[i1+1]-z1[i1+2])*(x1[i1+2]-x1[i1+4])-(x1[i1+1]-x1[i1+2])*(z1[i1+2]-z1[i1+4]));
     		p3=((x1[i1+1]-x1[i1+2])*(y1[i1+2]-y1[i1+4])-(y1[i1+1]-y1[i1+2])*(x1[i1+2]-x1[i1+4]));
     		e1=sqrt(c1*c1+c2*c2+c3*c3); 
     		f1=sqrt(p1*p1+p2*p2+p3*p3);
     		pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
     		g1=(x1[i1-1]-x1[i1+4]); 
     		g2=(y1[i1-1]-y1[i1+4]); 
     		g3=(z1[i1-1]-z1[i1+4]);
     		gg1=sqrt(g1*g1+g2*g2+g3*g3); 
     		hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
     		if (pp1<=-1.0) {di=-3.14;}
     		else if (pp1>=1.0) {di=0.;}
     		else if (hh1>=0.) {di=acos(pp1);}
     		else {di=-acos(pp1);}
     		ud0=15.0*((1-cos(di+2.92))+0.5*(1-cos(3.*(di+2.92))));
   	} 
   	return ud0;
}
float CPCN(int i1,float x1[1000],float y1[1000],float z1[1000])
{
  	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
  	if(i1-4<0) 
  	{
       		ud0=0.0;
  	}
  	else
  	{
     		c1=((y1[i1-4]-y1[i1-2])*(z1[i1-2]-z1[i1-1])-(z1[i1-4]-z1[i1-2])*(y1[i1-2]-y1[i1-1]));
     		c2=((z1[i1-4]-z1[i1-2])*(x1[i1-2]-x1[i1-1])-(x1[i1-4]-x1[i1-2])*(z1[i1-2]-z1[i1-1]));
     		c3=((x1[i1-4]-x1[i1-2])*(y1[i1-2]-y1[i1-1])-(y1[i1-4]-y1[i1-2])*(x1[i1-2]-x1[i1-1]));
     		p1=((y1[i1-2]-y1[i1-1])*(z1[i1-1]-z1[i1])-(z1[i1-2]-z1[i1-1])*(y1[i1-1]-y1[i1]));
     		p2=((z1[i1-2]-z1[i1-1])*(x1[i1-1]-x1[i1])-(x1[i1-2]-x1[i1-1])*(z1[i1-1]-z1[i1]));
     		p3=((x1[i1-2]-x1[i1-1])*(y1[i1-1]-y1[i1])-(y1[i1-2]-y1[i1-1])*(x1[i1-1]-x1[i1]));
     		e1=sqrt(c1*c1+c2*c2+c3*c3); 
     		f1=sqrt(p1*p1+p2*p2+p3*p3);
     		pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
     		g1=(x1[i1-4]-x1[i1]); 
     		g2=(y1[i1-4]-y1[i1]); 
     		g3=(z1[i1-4]-z1[i1]);
     		gg1=sqrt(g1*g1+g2*g2+g3*g3); 
     		hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
     		if (pp1<=-1.0) {di=-3.14;}
     		else if (pp1>=1.0) {di=0.;}
     		else if (hh1>=0.) {di=acos(pp1);}
     		else {di=-acos(pp1);}
     		ud0=0.825*((1-cos(di+1.164))+0.5*(1-cos(3.*(di+1.164))));
   	}
   	return ud0;
}
float NCPC(int i1,float x1[1000],float y1[1000],float z1[1000])
{
 	float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
 	if(i1+2>N0) 
 	{
     		ud0=0.0;
 	}
 	else
 	{
   		c1=((y1[i1]-y1[i1-1])*(z1[i1-1]-z1[i1+1])-(z1[i1]-z1[i1-1])*(y1[i1-1]-y1[i1+1]));
   		c2=((z1[i1]-z1[i1-1])*(x1[i1-1]-x1[i1+1])-(x1[i1]-x1[i1-1])*(z1[i1-1]-z1[i1+1]));
   		c3=((x1[i1]-x1[i1-1])*(y1[i1-1]-y1[i1+1])-(y1[i1]-y1[i1-1])*(x1[i1-1]-x1[i1+1]));
   		p1=((y1[i1-1]-y1[i1+1])*(z1[i1+1]-z1[i1+2])-(z1[i1-1]-z1[i1+1])*(y1[i1+1]-y1[i1+2]));
   		p2=((z1[i1-1]-z1[i1+1])*(x1[i1+1]-x1[i1+2])-(x1[i1-1]-x1[i1+1])*(z1[i1+1]-z1[i1+2]));
   		p3=((x1[i1-1]-x1[i1+1])*(y1[i1+1]-y1[i1+2])-(y1[i1-1]-y1[i1+1])*(x1[i1+1]-x1[i1+2]));
   		e1=sqrt(c1*c1+c2*c2+c3*c3);
   		f1=sqrt(p1*p1+p2*p2+p3*p3);
   		pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
   		g1=(x1[i1]-x1[i1+2]); 
   		g2=(y1[i1]-y1[i1+2]); 
   		g3=(z1[i1]-z1[i1+2]);
   		gg1=sqrt(g1*g1+g2*g2+g3*g3); 
   		hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
   		if (pp1<=-1.0) 
   		{
         		di=-3.14;
   		}
   		else if (pp1>=1.0) 
   		{
        		di=0.;
   		}
   		else if (hh1>=0.) 
   		{
        		di=acos(pp1);
   		}
   		else 
   		{
        		di=-acos(pp1);
   		}
   		ud0=0.76*((1-cos(di-0.88))+0.5*(1-cos(3.*(di-0.88))));
  	}
   	return ud0;
}  
