#include "src/include.h"

/* ----------------------------------------------------- */
/* Calculation of the abundance of the elements from BBN */
/* ----------------------------------------------------- */

void read_csv(int row, int col, char *filename, double **data)
{
	FILE *file;
	file = fopen(filename, "r"); //what's the r?

	int i = 0;
	char line[4098];
	while (fgets(line, 4098, file) && (i<row))
	{
		//double row[ssParams->nreal+1]; //what??
		char *tmp = strdup(line);

		int j = 0;
		const char* tok;
		for (tok = strtok(line, ","); tok && *tok; j++, tok=strtok(NULL, ","))
		{
			data[i][j] = atof(tok);
		}
		free(tmp);
		i++;
	}
}

int main(int argc, char** argv)
{
	struct relicparam paramrelic;
	double ratioH[NNUC+1],cov_ratioH[NNUC+1][NNUC+1];
	double H2_H,He3_H, Yp, Li7_H, Li6_H, Be7_H;
	double sigma_H2_H, sigma_He3_H, sigma_Yp,sigma_Li7_H, sigma_Li6_H, sigma_Be7_H;
	double eta, Tinit;
	int row;
	int col = 3;
	char fname_Tcm[256];
	char fname_a[256];
	char fname_b[256];
	char fname_c[256];
	char fname_d[256];
	char folder[256] = "../";
	double **data_Tcm;
	double **data_a;
	double **data_b;
	double **data_c;
	double **data_d;
	char mass[256];
	char mix[256];

	int failsafe;

	if (argc<4)
	{
		printf(" This program needs 6 parameters:\n"
		"   Tinit       initial temperature in MeV (2.3267 MeV by default)\n"
		"   eta0        initial value of the baryon-to-photon ratio (default: 6.1e-10)\n"
		"   failsafe    0=fast, 1=precise, 6=robust but slow. See stand_cosmo.c for more options.\n"
		"   row         number of rows in the rho_nu file, including the header row.\n"
		"   mass        The mass of the sterile neutrino model of interest.\n"
	  "   mix         The mixing angle of the sterile neutrino model of interset.\n");
			exit(1);
	}
	else
	{
		sscanf(argv[1],"%lf",&Tinit); //default is Tinit=27.*K_to_eV*1.e3
		sscanf(argv[2],"%lf",&eta); //default is eta=6.1e-10;
		sscanf(argv[3],"%d",&failsafe); //default is failsafe=1;
		sscanf(argv[4],"%d",&row); //number of rows in rho_nu and drho_nu
		sscanf(argv[5],"%s",mass);
		sscanf(argv[6],"%s",mix);
	}

	strcat(folder,mass); strcat(folder,"-"); strcat(folder,mix); strcat(folder,"-FullTestNew/"); strcat(folder,mass); strcat(folder,"-"); strcat(folder,mix); strcat(folder,"-FullTestNew/"); strcat(folder, "mass_"); strcat(folder,mass); strcat(folder, "_mix_"); strcat(folder, mix);

	strcpy(fname_Tcm, folder); strcat(fname_Tcm, "_Tcm.csv");
	strcpy(fname_a, folder); strcat(fname_a, "_a.csv");
	strcpy(fname_b, folder); strcat(fname_b, "_b.csv");
	strcpy(fname_c, folder); strcat(fname_c, "_c.csv");
	strcpy(fname_d, folder); strcat(fname_d, "_d.csv");
	//printf("%s \n %s \n %s \n %s \n %s \n",fname_Tcm,fname_a,fname_b,fname_c,fname_d);

	data_Tcm = (double **)malloc(row * sizeof(double *));
	data_a = (double **)malloc(row * sizeof(double *));
	data_b = (double **)malloc(row * sizeof(double *));
	data_c = (double **)malloc(row * sizeof(double *));
	data_d = (double **)malloc(row * sizeof(double *));
	double dataout[row-1][col-1];

	for (int i = 0; i<row; i++)
	{
		data_Tcm[i] = (double *)malloc(col * sizeof(double));
		data_a[i] = (double *)malloc(col * sizeof(double));
		data_b[i] = (double *)malloc(col * sizeof(double));
		data_c[i] = (double *)malloc(col * sizeof(double));
		data_d[i] = (double *)malloc(col * sizeof(double));
	}

	read_csv(row, col, fname_Tcm, data_Tcm);
	read_csv(row, col, fname_a, data_a);
	read_csv(row, col, fname_b, data_b);
	read_csv(row, col, fname_c, data_c);
	read_csv(row, col, fname_d, data_d);

	for (int i=1; i<row; i++)
	{
		dataout[i-1][0] = data_Tcm[i][1];
		dataout[i-1][1] = data_a[i][1];
		dataout[i-1][2] = data_b[i][1];
		dataout[i-1][3] = data_c[i][1];
		dataout[i-1][4] = data_d[i][1];
	}

	for (int i=0; i<row-1; i++)
	{
		paramrelic.Tcm_arr[i] = dataout[i][0];
		paramrelic.arho[i] = dataout[i][1];
		paramrelic.brho[i] = dataout[i][2];
		paramrelic.crho[i] = dataout[i][3];
		paramrelic.drho[i] = dataout[i][4];
		//printf("%e, %e, %e, %e, %e \n", paramrelic.Tcm_arr[i], paramrelic.arho[i], paramrelic.brho[i], paramrelic.crho[i], paramrelic.drho[i]);
	}

	Init_cosmomodel(&paramrelic); //this sets the values of paramrelic that *every* model will use, other init functions set specific values that only that model uses

	paramrelic.failsafe = failsafe;
	paramrelic.row = row;

	Init_cosmomodel_param(eta,paramrelic.Nnu,0.,paramrelic.life_neutron,paramrelic.life_neutron_error,0.,0.,0.,&paramrelic);
	// why are Nnu, life_neutron, and life_neutron_error defined above by Init_cosmomodel_param when alter_eta supposedly just... alters eta?
	paramrelic.Tinit=Tinit*1.e-3/K_to_eV; //setting Tinit, why isn't this done in Init_cosmomodel_param?

	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
	paramrelic.err=2;
	nucl(&paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	printf("  low:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

	paramrelic.err=0;
	nucl(&paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	printf(" cent:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

	paramrelic.err=1;
	nucl(&paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	printf(" high:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

	paramrelic.err=3;
	if(nucl_err(&paramrelic,ratioH,cov_ratioH)==1)
	{
		printf("--------------------\n\n");
		printf("With uncertainties:\n");
		H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");

		printf("value:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

		double corr_ratioH[NNUC+1][NNUC+1];
		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
		printf("Correlation matrix:\n");
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		printf("Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
		printf("H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
		printf("He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
		printf("Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
		printf("Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
		printf("Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
	}
	else printf("Uncertainty calculation failed\n\n");

	/*paramrelic.err=4;
	if(nucl_err(&paramrelic,ratioH,cov_ratioH))
	{
		printf("--------------------\n\n");
		printf("With MC uncertainties:\n");
				H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");

		printf("mean:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

		double corr_ratioH[NNUC+1][NNUC+1];
		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
		printf("Correlation matrix:\n");
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		printf("Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
		printf("H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
		printf("He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
		printf("Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
		printf("Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
		printf("Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
	}
	else printf("Uncertainty calculation failed\n\n");*/

	paramrelic.err=0;
	int compat=bbn_excluded(&paramrelic);

	if(compat==1) printf("Excluded by BBN constraints (chi2 without correlations)\n");
	else if(compat==0) printf("Compatible with BBN constraints (chi2 without correlations)\n");
	else printf("Computation failed (chi2 without correlations)\n");

	paramrelic.err=3;
	compat=bbn_excluded(&paramrelic);

	if(compat==1) printf("Excluded by BBN constraints (chi2 including correlations)\n");
	else if(compat==0) printf("Compatible with BBN constraints (chi2 including correlations)\n");
	else printf("Computation failed (chi2 including correlations)\n");

	return 1;
}
