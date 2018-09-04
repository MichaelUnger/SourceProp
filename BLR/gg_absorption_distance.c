#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"

#define qu(x) ((x)*(x))
#define pc 3.08567758e18 //cm
#define c 2.99792458e10 //cm.s-1
#define nu2eps 8.093299728e-21
#define nu2GeV 4.135666991e-24
#define sigma_T 6.65e-25
#define EOL(fp) fscanf(fp, "%*[^\n]")
#define pi 3.14159265359





double sigma(e, e1, mu)
double e, e1, mu;
{
double w, beta, x, y;

if ((x=e*e1*(1. - mu)) <= 2.) return (0.);

if (x > 1.e3)
     beta = 1. - 1./x;
else beta = sqrt(1. - 2./x);
if (x > 1.e3)
     y = 2.*x;
else y = (1. + beta)/(1. - beta);

w = (1. - qu(beta))*((3. - pow(beta, 4.))*log(y) - 2.*beta*(2. - qu(beta)));

return (w);
}



int w_in(count)
int count;
{
  return(3 + pow(-1,count+1));
}

/**cc_abs(BLR inner radius, BLR outer radius, Injection height [cm], External eps, External nph, nFn_#.dat file, enFn_#.dat, number of lines in nFn_#.dat);**/
/* void cc_abs(double L_D, double R1, double R2, double Z_0, double Arr_e[], double Arr_ne[], char datpool[60], char datext_ab[40], int num_ext) */
void main()
{
  double nu, F;
  double Arr_e[200], Arr_ne[200];
  double R1, R2, L_D=1.e45;
  double eps0=1.e4, eps1=1.e7, de=1.5, epsilon[100];
  double R_BLR, u_BLR, L_BLR, d_range, R3, L_0;
  double eps, eps_th, mu;
  double tau;
  double em_d, d, N_ph, tL, em_end, N_ph_0;
  double d_p, dR, DR;
  double scale;
  double uinput=0.;
  int i, k, l, n, j, p, y, u;
  int q, r, t, num_ext;
  
  int mark;
  
  double beta1, beta2, beta3, beta4;
  double gamma1, gamma2, gamma3, gamma4;
  double L, alpha, beta, gamma, delta;
  double L1, L2, L3, L4;
  double dln;
  double th, dth, th_c_i, th_c_o, sin_th, sin_alpha;
  int dl, dtheta, dlf;
  int dis_flag, f;
  
  double con;
  double sum_ep, sum_mu, sum_l_1, sum_l_2, sum_max, sum_thr, sum_mu_p;
  
/* char specfile[sizeof "dir_@BLZ/spec_100.dat"]; */
//   char evtfile[sizeof "dir_3C279/evt_100.dat"];
  char tvdfile[100], tvd0[100];
  char sourcename[50]="PKS0736+017", logfile[50]="PKS0736+017.log";
  char a[50];
  FILE *fp_tvd, *fp_r_nuFnu, *fp_out, *fp_tau_v_E, *fp_log, *fp_scale, *fpext;
  
  int par1, par2;
  setbuf(stdout, NULL);//Turns off buffered output for debugging purposes. 



  //3C279
  // u_BLR = 1e-2;
  //
  //PKS1510
  //u_BLR = 4.5e-3;
  //
  
  // PKS 0736+017
  // u_BLR = 9.45e-4;
  
  double par3=0.1;

/* L_BLR = par3*L_D;//choose 3 values (PKS1510)  */

/*  L_BLR = 2.e44;  // 3C279  */

//  L_BLR = 1.5e44;  // PKS 0736+017

// CTA 102 //

  printf ("\n Source name = [%s]  ", sourcename);
  gets(a);
  if (strlen(a)) strcpy(sourcename, a);

  L_BLR = 2e44;
  R_BLR = 2.3e17;//2.3e17;
  
/* u_BLR = 2.45e-2; 
   R_BLR = sqrt(L_BLR/(4*pi*u_BLR*c));  */

  printf ("\n L_BLR = [%e] erg/s  ", L_BLR);
  gets(a);
  if (strlen(a)) L_BLR = atof(a);
  
  printf ("\n R_BLR = [%e] cm  ", R_BLR);
  gets(a);
  if (strlen(a)) R_BLR = atof(a);
  
  
  u_BLR = L_BLR/(4.*pi*qu(R_BLR)*c);

  strcpy(logfile, sourcename);
  l = strlen(sourcename);
  logfile[l] = '.';
  logfile[l+1] = 'l';
  logfile[l+2] = 'o';
  logfile[l+3] = 'g';
  logfile[l+4] = '\0';

  strcpy(tvd0, sourcename);
  tvd0[l] = '/';
  tvd0[l+1] = 't';
  tvd0[l+2] = 'v';
  tvd0[l+3] = 'd';
  tvd0[l+4] = '\0';
  
  printf ("\n Log file = %s; tvd0 = %s \n", logfile, tvd0);

  fp_log = fopen(logfile,"w");
  fprintf(fp_log,"Log File for Source %s\n", sourcename);  
  
  
  L_0 = L_BLR;
  y=100;

  double Arr_d[y];


  printf(" R_BLR = %e cm;  L_BLR = %e erg/s;  u_BLR = %e erg/cm^3 \n", R_BLR, L_BLR, u_BLR);
  fprintf(fp_log, "\n R_BLR = %e cm;  L_BLR = %e erg/s;  u_BLR = %e erg/cm^3 \n", R_BLR, L_BLR, u_BLR);

  R1 = 0.9*R_BLR;
  R2 = 1.1*R_BLR;
  R3 = 10.*R_BLR;
  em_d = 0.2*R1;
  em_end = 5.*R2;
  d_range = (em_end - em_d)/y;
  
  fprintf(fp_log,"\n R1 = %e cm \n R2 = %e cm \n R_start = %e cm \n R_end = %e cm \n\n", R1, R2, em_d, em_end);
  printf("\n R1 = %e cm \n R2 = %e cm \n R_start = %e cm \n R_end = %e cm \n",R1,R2,em_d,em_end);
  sum_thr = 1e-8;
  
/**/		
/*  printf("\n\tReading in spectrum file...");
  if (!(fp_r_nuFnu = fopen(datpool, "r"))) file_error(datpool);
  EOL(fp_r_nuFnu);
  r = ftell(fp_r_nuFnu) +1;
  fseek(fp_r_nuFnu, 0, SEEK_END);
  t = ftell(fp_r_nuFnu);
  fseek(fp_r_nuFnu, 0, SEEK_SET);
  q = t/r; //length of the datpool file (datpool must be in block format)
  */

  k=0;
  /* fp_scale = fopen("CTA102/L_scale.dat","w");  */
  
  printf ("\n Reading external photon file... ");
  
  if (!(fpext=fopen("BLR_1e-2b.dat", "r"))) 
    { printf ("\n Error opening photon spectrum file! \n\n");
      exit(0);
      }
  while (((int)(fgetc(fpext))) != EOF)
  {
    fseek(fpext,-1,1);
    fscanf(fpext,"%lf",&nu);
    fscanf(fpext,"%lf",&F);
    Arr_e[k] = nu;
    Arr_ne[k] = F;
    if (k) uinput+=(Arr_ne[k]*Arr_e[k]*8.176e-7*(Arr_e[k] - Arr_e[k-1]));
    k++;
  }
  fclose(fpext);
  printf(" \n u_input = %e erg/cm^3 \n Done! \n", uinput);
  num_ext = k;  

for (k=0; k<num_ext; ++k) Arr_ne[k]*=(u_BLR/uinput);
  
printf("\n\tCalculating normalisation for em_d within BLR ...");
tL = 0;
n = 150;
f = 400;
dth = pi/((double)(n));
for (p=0; p < n+1; p++)
   { th = p*dth;
     sin_th = sin(th);
     L = sqrt(-pow(em_d*sin_th,2)+pow(R2,2))-sqrt(-pow(em_d*sin_th,2)+pow(R1,2));
     if (p == 0 || p == n) tL += L*sin_th;
     else tL += L*sin_th*w_in(p);
     }
tL *= dth/3.;
N_ph = 2*c/tL;
con = (3.*N_ph*sigma_T)/(16.*4.*pi*c);

printf ("\n N_ph = %e, con = %e\n", N_ph, con);

j = 0;
eps = eps0;
do 
  { epsilon[j] = eps;
    eps*=de;
    j++;
    } while (eps < eps1);
q = j;

printf ("\n epsilon array set up. \n");

  double Arr_cat[y][q];  //Arr[d][eps]=tau
/**/
  mark=0;
  
/*   Loop over blob distance   */

  for (u=0; u < y; u++)
   {
    DR = 1.075;
    n = 150;
    f = 400;
    dth = pi/n;
    if (em_d < R1)
      { dln = (R3-em_d)/f;
        }
    else
      { dln = (R3-R1)/f;
        }
//     sprintf(specfile,"dir_@BLZ/spec_%02d.dat",u);
//     sprintf(evtfile,"dir_PKS1510/%d.dat",u);
//     fp_out = fopen(specfile,"w");
//     fp_tau_v_E = fopen(evtfile,"w");
    
/**/

    L_BLR = L_0;
 /* fprintf(fp_scale,"%.8e %.8e\n", em_d, L_BLR); */
       
    //norm end
    printf("\n\tN_ph = %.8e\n\td=%.8e \n\tDone!",N_ph,em_d/R_BLR);
/**/
    con = (3.*N_ph*sigma_T)/(16.*4.*pi*c);
    printf("\n\tCalculating absorption due to BLR...\n");
    
    
/*  Loop over photon energies   */
 
    for (j=0; j < q; j++)
    {
      L=0;
      eps = epsilon[j];
      if (mark == 0)
	  { mark=j;
	    /* fprintf(fp_log,"%d mark\n%d sam\n\n",mark,q-1); */
	    }	
	dl 	= 0;
	sum_l_1	= 0;
	sum_l_2	= 0;
	sum_max	= 0;
	dis_flag= 0;
	tau=0;
	while (dis_flag == 0) //Infinite boundary on the length integral
	{
	  if ((dl > f) && (d > R2))
	  {
	    dlf = dl-f-1;
	    dR = dln*pow(DR,dlf);
	    d += dR;
	  }else{
	    d = em_d + dln*dl;
	  }
	  sum_mu = 0.;
	  for (i=1; i < n+1; i++)
	  {
	    th = pi - dth*i;
	    mu = cos(th);
	    sin_th = sin(th);
	    eps_th = 2./(eps*(1. + mu));
	    if (d < R1) 			//Before BLR, interaction from the front
	    {
	      L = sqrt(-pow(d*sin_th,2)+pow(R2,2))-sqrt(-pow(d*sin_th,2)+pow(R1,2));
	    }else if (d > R2){			//Beyond BLR, interaction from the back
	      th_c_i = pi - asin(R1/d);
	      th_c_o = pi - asin(R2/d);
	      alpha = pi - th;
	      sin_alpha = sin(alpha);
	      if (th > th_c_i) 
	      {
		gamma1 = pi - asin(d*sin_alpha/R2);
		gamma2 = pi - asin(d*sin_alpha/R1);
		gamma3 = asin(d*sin_alpha/R1);
		gamma4 = asin(d*sin_alpha/R2);
		beta1 = pi - gamma1 - alpha;
		beta2 = pi - gamma2 - alpha;
		beta3 = pi - gamma3 - alpha;
		beta4 = pi - gamma4 - alpha;
		L1 = d*sin(beta1)/sin(gamma1);
		L2 = d*sin(beta2)/sin(gamma2);
		L3 = d*sin(beta3)/sin(gamma3);
		L4 = d*sin(beta4)/sin(gamma4);
		L = L4 - L3 + L2 - L1;
	      }else if ((th > th_c_o) && (th < th_c_i)){
		gamma1 = pi - asin(d*sin_alpha/R2);
		gamma4 = asin(d*sin_alpha/R2);
		beta1 = pi - gamma1 - alpha;
		beta4 = pi - gamma4 - alpha;
		L1 = d*sin(beta1)/sin(gamma1);
		L4 = d*sin(beta4)/sin(gamma4);
		L = L4 - L1;
	      }else{
		L = 0;
		break;
	      }
            } else {				//Enters the BLR
	      th_c_i = pi - asin(R1/d);
	      alpha = pi - th;
	      sin_alpha = sin(alpha);
	      if (th > th_c_i)
	      {
		gamma1 = pi - asin(d*sin_alpha/R1);
		gamma2 = asin(d*sin_alpha/R1);
		gamma3 = asin(d*sin_alpha/R2);
		beta1 = pi - gamma1 - alpha;
		beta2 = pi - gamma2 - alpha;
		beta3 = pi - gamma3 - alpha;
		L1 = R1*sin(beta1)/sin_alpha;
		L2 = R1*sin(beta2)/sin_alpha;
		L3 = R2*sin(beta3)/sin_alpha;
		L = L3 - L2 + L1;
	      }else{
		gamma = asin(d*sin_alpha/R2);
		beta = pi - gamma - alpha;
		L = R2*sin(beta)/sin_alpha;
	      }
	    }
            //            printf("asdf %e %e %e\n", L/R_BLR, d/R_BLR, mu);
	    sum_ep = 0;
	    for (k=1; k < num_ext; k++)
	    {
	      if ( eps_th < Arr_e[k] && eps_th > Arr_e[k-1])
	      {
		de = Arr_e[k] - eps_th;
		sum_ep += de*(Arr_ne[k]*L*sigma(eps,Arr_e[k],-mu));
	      }else
	      if ( eps_th < Arr_e[k] && eps_th < Arr_e[k-1]){
		de = Arr_e[k] - Arr_e[k-1];
		sum_ep += de*(Arr_ne[k]*L*sigma(eps,Arr_e[k],-mu));
	      }
	    }// EPSILON INTEGRAL END
	    if (i == 1 || i == n)
	    {
	      sum_mu += (1+mu)*sin_th*sum_ep*con;
	    }else{
	      sum_mu += (1+mu)*sin_th*sum_ep*w_in(i)*con;
	    }
	  }// MU INTEGRAL END
	  sum_mu *= dth/3.;
	  if ((dl < f+1) && (d < R3))
	  {
	    if (dl == 0 || dl == f)
	    {
	      sum_l_1 += sum_mu;  
	    }else{
	      sum_l_1 += sum_mu*w_in(dl);
	    }
	  }else{
	    /*ESCAPE PARAMETERS*/
	    if (sum_mu == 0){tau=0;break;}
	    if (sum_mu > sum_max){sum_max = sum_mu;}
	    if (sum_mu/sum_max < sum_thr){break;}
	    sum_l_2 += dR*(sum_mu + sum_mu_p)/2.;
	  }
	  sum_mu_p = sum_mu;
	  dl++;
	}// LENGTH INTEGRAL END
	sum_l_1 *= dln/3.;
	tau = sum_l_1 + sum_l_2;
        printf("asdf %e %e %e %e\n", eps*5.11e-4, tau,sum_l_1 , sum_l_2);
//       fprintf(fp_out,"%.6e %.6e %.6e\n",Arr_nu[j],Arr_F[j]*exp(-tau));
//       fprintf(fp_tau_v_E,"%.6e %.6e\n",Arr_nu[j],tau);
      Arr_cat[u][j] = tau;
    }
//     fclose(fp_out);
//     fclose(fp_tau_v_E);
    printf(" Done!\n");
    Arr_d[u] = em_d;
    em_d += d_range;
  }
/**/
  for (par1 = mark; par1 < q ;par1++)//energy
  {
    sprintf(tvdfile,"%s_%d.dat", tvd0, par1);
    printf ("\n tvdfile = %s ", tvdfile);
    fprintf(fp_log,"E[%d] = %.8e GeV\n", par1, epsilon[par1]*5.11e-4);
    fp_tvd = fopen(tvdfile,"w");
    for (par2 = 0; par2 < y; par2++)//distance
    {
      fprintf(fp_tvd,"%.8e %.8e\n", Arr_d[par2], Arr_cat[par2][par1]);
//       fprintf(fp_log,"%.6e %.6e\n", Arr_d[par2], Arr_cat[par2][par1]);
    }
    fprintf(fp_log,"\n");
    fclose(fp_tvd);
  }
fclose(fp_log);
printf ("\n All done! \n\n");
}
