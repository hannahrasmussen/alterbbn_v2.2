#include<stdio.h>
#include<string.h>

void create_csv(char *filename,double a[][6],int m,int n)
{
	printf("Creating %s.csv file\n",filename);
	FILE *fp;
	int i,j;
	filename=strcat(filename,".csv");
	fp=fopen(filename,"w+");
	fprintf(fp,"bound, Yp, H2/H, He3/H, Li7/H, Li6/H, He7/H");
	
	for(i=0;i<m;i++)
	{
		if (i==0) fprintf(fp,"\n%s","low");
		else if (i==1) fprintf(fp, "\n%s", "cent");
		else fprintf(fp, "\n%s", "high");
    		for(j=0;j<n;j++)
		{
        		fprintf(fp,",%e ",a[i][j]);
    	
		}
	}
 
	fclose(fp);
}
 
int main()
{
    	double a[3][6]={{2.474e-01, 2.526e-05, 1.025e-05, 5.028e-10, 1.689e-15, 4.745e-10},
		{2.473e-01, 2.463e-05, 1.034e-05, 5.376e-10, 1.085e-14, 5.087e-10},
		{2.473e-01, 2.404e-05, 1.044e-05, 5.746e-10, 3.522e-14, 5.454e-10}};
	char str[100] = "trythis";
	create_csv(str,a,3,6);
	return 0;
}
