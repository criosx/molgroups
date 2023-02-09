#include <iostream>
#include <cassert>
#include "stdio.h"
#include "refl.h"
#include "reflcalc.h"
#include "molgroups.cc"

#define FWHM 2.354820045   // ga_refl uses FWHM. Divide by FWHM to get sigma units.

#define MODELS 6
#define GAUSSSTART 4
#define DIMENSION 400
#define STEPSIZE 0.5

#define RHOPROTPROT 1.73e-6
#define RHOPROTDEUT 2.99e-6

/* initialising non-standard fitting variables */
double aArea[DIMENSION], anSL[DIMENSION];

double normarea, nF_PS; 
double vf_bilayer,  global_rough, rough_au_cr, mult_tether, l_tether, penetration;
double sigma, nf_tether, nf_PS, nf_chol, ExchProt, scalarmult_off1, scalarmult_off2;
double d_on1, d_on2, d_on3, d_on4, d_on5, d_on6, d_on7;
double vf_on1, vf_on2, vf_on3, vf_on4, vf_on5, vf_on6;
double dl_lipid_on1, dl_lipid_off, protnf;
double l_lipid1, l_lipid2, rho_solv_h2o;

tBLM_HC18_DOPC_DOPS_CHOL  bilayer_neat;  
tBLM_HC18_DOPC_DOPS_CHOL  bilayer_on1;
Discrete                  prot_on1(50,0,"DiscreteTriangle.dat");


/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
  int bPeaked,i,k;
  double dMaxArea, rhoprot;
  fitpars *pars = &fit[0].pars;

  for (i=0; i<pars_count(pars); i++)
  {
    if (pars_peek(pars,i)==pars_max(pars,i))
    {pars_poke(pars,i,pars_min(pars,i)+0.99*(pars_max(pars,i)-pars_min(pars,i)));}
    if (pars_peek(pars,i)==pars_min(pars,i))
    {pars_poke(pars,i,pars_min(pars,i)+0.01*(pars_max(pars,i)-pars_min(pars,i)));}
  }

/* Rescue the free parameters from the model. */
  for (i=0; i < fit[1].pars.n; i++)
    fit[1].pars.value[i] = *(fit[1].pars.address[i]);

/* Go through all layers copying parameters from model 0 to other models */
  tied_parameters(fit);
 
 /* copy the global roughness to all interfaces not SIOX interfaces*/
  for (i=0; i< MODELS; i++) { //if(rough_au_cr>fit[0].m.d[2]) {rough_au_cr=fit[0].m.d[2];}
    for (k=1;k<GAUSSSTART; k++) fit[i].m.rough[k]=FWHM*global_rough;
	  fit[i].m.rough[3]=rough_au_cr;
  }    
  
  /* Restore the free parameters to the model. */
  for (i=0; i < fit[1].pars.n; i++){
    *(fit[1].pars.address[i]) = fit[1].pars.value[i];
  }
    
    fit[2].m.rho[fit[2].m.n-1]=rho_solv_h2o;
    fit[3].m.rho[fit[3].m.n-1]=rho_solv_h2o;
    fit[4].m.rho[fit[4].m.n-1]=rho_solv_h2o;
    
  // neat bilayer
  
    bilayer_neat.sigma=sigma;
    bilayer_neat.global_rough=global_rough;
    bilayer_neat.vf_bilayer=vf_bilayer;
    bilayer_neat.rho_substrate=fit[0].m.rho[3];
    bilayer_neat.nf_tether=nf_tether;
    bilayer_neat.mult_tether=mult_tether;
    bilayer_neat.l_tether=l_tether;
    bilayer_neat.l_lipid1=l_lipid1;
    bilayer_neat.l_lipid2=l_lipid1;
    bilayer_neat.nf_lipid_2=nf_PS;
    bilayer_neat.nf_chol=nf_chol;

    bilayer_neat.fnAdjustParameters();
    mult_tether=bilayer_neat.mult_tether;
    l_tether=bilayer_neat.l_tether;
    dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
    dMaxArea=bilayer_neat.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
    normarea=dMaxArea;

    fnWriteCanvas2Model(aArea,anSL,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,0,2);
  
  // on1
    
    bilayer_on1.sigma=sigma;
    bilayer_on1.global_rough=global_rough;
    bilayer_on1.vf_bilayer=vf_bilayer;
    bilayer_on1.rho_substrate=fit[0].m.rho[3];
    bilayer_on1.nf_tether=nf_tether;
    bilayer_on1.mult_tether=mult_tether;
    bilayer_on1.l_tether=l_tether;
    bilayer_on1.l_lipid1=l_lipid1;
    bilayer_on1.l_lipid2=l_lipid1+dl_lipid_on1;
    bilayer_on1.nf_lipid_2=nf_PS;
    bilayer_on1.nf_chol=nf_chol;
    
    bilayer_on1.fnAdjustParameters();
    dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
    dMaxArea=bilayer_on1.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
    normarea=dMaxArea;
 
    if (scalarmult_off2>scalarmult_off1) {scalarmult_off2=scalarmult_off1;}
    
    prot_on1.fnSetNormarea(normarea);
    prot_on1.dProtExchange=ExchProt;
    prot_on1.dnSLDBulkSolvent=fit[3].m.rho[fit[3].m.n-1];
    prot_on1.nf=protnf;
    prot_on1.fnOverlayProfile(aArea,anSL,DIMENSION,STEPSIZE,dMaxArea);
    fnWriteCanvas2Model(aArea,anSL,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,3,3);

  // off1
    bilayer_on1.l_lipid2=l_lipid1+dl_lipid_off;
    bilayer_on1.fnAdjustParameters();
    
    dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
    dMaxArea=bilayer_on1.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
    normarea=dMaxArea;    
    prot_on1.fnSetNormarea(normarea);
    prot_on1.dProtExchange=ExchProt;
    prot_on1.dnSLDBulkSolvent=fit[4].m.rho[fit[4].m.n-1];
    prot_on1.nf=protnf;//*scalarmult_off1;
    prot_on1.fnOverlayProfile(aArea,anSL,DIMENSION,STEPSIZE,dMaxArea);
    fnWriteCanvas2Model(aArea,anSL,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,4,4);
    
  // off2
    bilayer_on1.l_lipid2=l_lipid1+dl_lipid_off;
    bilayer_on1.fnAdjustParameters();
    
    dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
    dMaxArea=bilayer_on1.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
    normarea=dMaxArea;    
    prot_on1.fnSetNormarea(normarea);
    prot_on1.dProtExchange=ExchProt;
    prot_on1.dnSLDBulkSolvent=fit[5].m.rho[fit[5].m.n-1];
    prot_on1.nf=protnf;//*scalarmult_off2;
    prot_on1.fnOverlayProfile(aArea,anSL,DIMENSION,STEPSIZE,dMaxArea);
    fnWriteCanvas2Model(aArea,anSL,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,5,5);
    
}

void save(fitinfo *fit)
{
        FILE *fp;
        fp=fopen("mol.dat","w");        
        bilayer_on1.fnWriteGroup2File(fp,"bilayer",DIMENSION,STEPSIZE);
        prot_on1.fnWriteGroup2File(fp,"protein",DIMENSION,STEPSIZE);
        fclose(fp);
}


/*============ INITIAL SETUP ============================*/
extern "C"
fitinfo* setup_models(int *models)
{
  static fitinfo fit[MODELS];
  int i,j;
  fitpars *pars = &fit[0].pars;
  fitpars *freepars = &fit[1].pars;
  *models = MODELS;

  for (i=0; i < MODELS; i++) fit_init(&fit[i]);
  
  /* Load the data for each model */
    fit_data(&fit[0],"p0047.refl");
    fit_data(&fit[1],"p0048.refl");
    fit_data(&fit[2],"p0049.refl");
    fit_data(&fit[3],"p0050.refl");
    fit_data(&fit[4],"p0051.refl");
    fit_data(&fit[5],"p0052.refl");
  //fit_data(&fit[8],"sy006.refl"); //0.25 muM
  //fit_data(&fit[9],"sy005.refl");

  /* Initialize instrument parameters for each model.*/
  /* setup for NG7 */
  for (i=0; i < MODELS; i++) {

    const double L = 4.75,dLoL=0.02,d=1549.0;
    double Qlo, Tlo, dTlo,dToT,s1,s2;
    Qlo=0.0175;
    Tlo=Q2T(L,Qlo);
    s1=0.1, s2=s1;
    dTlo=resolution_dT(s1,s2,d);
    dToT=resolution_dToT(s1,s2,d,Tlo);
    data_resolution_fv(&fit[i].dataA,L,dLoL,Qlo,dTlo,dToT);
    fit[i].beam.lambda = L;
    interface_create(&fit[i].rm, "erf", erf_interface, 20);
  }

  /*============= MODEL =====================================*/

  /* Add layers: d, rho, mu, rough */
  for (i=0; i < MODELS; i++) {
      model_layer(&fit[i].m, 0.00, 2.07e-6, 0.0e-8, 7.000);  	    /* 0 Si */
      model_layer(&fit[i].m, 8,    3.55e-6, 0.0e-8, 7.000);			/* 1 oxide */
      model_layer(&fit[i].m, 25,   3.95e-6, 0.0e-8, 7.000);			/* 2 Cr */
      model_layer(&fit[i].m, 130,  4.45e-6, 0.0e-8, 7.000);			/* 3 Au */
      for (j=0; j < DIMENSION; j++) {
	      model_layer(&fit[i].m, STEPSIZE, 0.00e-6, 0.0e-8, 0);		/* Gaussian */
      }
      model_layer(&fit[i].m, 100.000, 6.35e-6, 0.0e-8, 0.000);		/* Solvent */
  }
  
  /*correct solvent layers for different models*/
  /* fit[3].m.d[3] = ... */
    
    d_on1=30;
    d_on2=30;
    d_on3=32;
    d_on4=34;
    d_on5=36;
    d_on6=36;
    d_on7=36;
    vf_on1=0.1;
    vf_on2=0.1;
    vf_on3=0.2;
    vf_on4=0.1;
    vf_on5=0.1;
    vf_on6=0.1;
  
    fit[1].m.rho[fit[1].m.n-1]= 4e-6;
    fit[2].m.rho[fit[2].m.n-1]=-0.5666e-6;
    fit[3].m.rho[fit[3].m.n-1]=-0.5666e-6;
    fit[4].m.rho[fit[4].m.n-1]=-0.5666e-6;
    fit[5].m.rho[fit[5].m.n-1]=6.34-6;
    nf_PS=0.30;
    nf_chol=0.03;
    nf_tether=0.4;
    mult_tether=0.1;
    vf_bilayer=0.95;
    l_tether=15;
    l_lipid1=12;
    l_lipid2=12;
    dl_lipid_on1=0;
    penetration=0;
    global_rough=3;
    sigma=2.5;
 
    rho_solv_h2o=-5.5e-7;

  //headgroup.bWrapping=false;
  
  /*=============== FIT PARAMETERS ===============================*/

  /* Specify which parameters are your fit parameters. Parameters are fitted
   * to be the same in all datasets by default
   */
   
    pars_add(pars, "rho_SiOx", &(fit[0].m.rho[1]), 3.4e-6, 3.8e-6);
    pars_add(pars, "rho_Cr", &(fit[0].m.rho[2]), 3.0e-6, 4.3e-6);
    pars_add(pars, "rho_Au", &(fit[0].m.rho[3]), 4.3e-6, 4.6e-6);
      
    pars_add(pars, "d_oxide", &(fit[0].m.d[1]), 5., 50.);
    pars_add(pars, "d_Cr", &(fit[0].m.d[2]), 30., 50);
    pars_add(pars, "d_Au", &(fit[0].m.d[3]), 100., 120);

    pars_add(pars, "l_tether", &(l_tether), 10, 19);
    pars_add(pars, "l_lipid1", &(l_lipid1), 8, 22);
    //pars_add(pars, "l_lipid2", &(l_lipid2), 8, 22);
    pars_add(pars, "dl_lipid_on1", &(dl_lipid_on1), -4., 4.);
    pars_add(pars, "dl_lipid_off", &(dl_lipid_off), -4., 4.);
    //pars_add(pars, "dl_lipid_off2", &(dl_lipid_off2), -4., 4.);
    pars_add(pars, "penetration", &(penetration), -8., 0.);
    
    pars_add(pars, "nf_tether", &(nf_tether), 0.5, 1.0);
    pars_add(pars, "mult_tether", &(mult_tether), 0.1, 10.);  
    pars_add(pars, "vf_bilayer", &(vf_bilayer), 0.6, 1);  
    
    pars_add(pars, "StartPosProt", &(prot_on1.dStartPosition), 50., 140.);
    pars_add(pars, "NumberFracProt", &(protnf), 0.01, 3);
    pars_add(pars, "ExProtons", &(ExchProt), 0.3, 0.9);
    pars_add(pars, "scalarmult_off1", &(scalarmult_off1), 0.1, 1.0);
    pars_add(pars, "scalarmult_off2", &(scalarmult_off2), 0.1, 1.0);
    
    
    pars_add(pars, "rho_solv_0", &(fit[0].m.rho[fit[0].m.n-1]), 5.6e-6, 6.56e-6);
    pars_add(pars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), 3.6e-6, 4.20e-6);
    //pars_add(pars, "rho_solv_2", &(fit[2].m.rho[fit[2].m.n-1]), -.566e-6, .4e-6);
    pars_add(pars, "rho_solv_h2o", &(rho_solv_h2o), -.566e-6, -0.4e-6);
    pars_add(pars, "rho_solv_5", &(fit[5].m.rho[fit[5].m.n-1]), 3.6e-6, 4.0e-6);
 
    pars_add(pars, "global_rough", &(global_rough), 2.0, 6.);
    pars_add(pars, "sigma",        &(sigma), 2.0, 5.);
    pars_add(pars, "rough_au_cr", &(rough_au_cr), 2.0, 6.);

  //pars_add(pars, "background_0", &(fit[0].beam.background), 1e-9, 5e-7);
    pars_add(pars, "background_1", &(fit[1].beam.background), 1e-9, 1e-6);
    pars_add(pars, "background_2", &(fit[2].beam.background), 1e-9, 5e-6);
    pars_add(pars, "background_3", &(fit[3].beam.background), 1e-9, 5e-6);
    pars_add(pars, "background_4", &(fit[4].beam.background), 1e-9, 5e-6);
    pars_add(pars, "background_5", &(fit[5].beam.background), 1e-9, 5e-6);

  /* Build a list of 'free parameters' in fit[1].pars. These are
   * parameters for which the values are aloowed to differ from those
   * in model 0.  By default all values in all models are the same unless 
   * specified here. The range data is not used here, so set it to [0,1].  
   */

    pars_add(freepars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), 0, 1);
    pars_add(freepars, "rho_solv_2", &(fit[2].m.rho[fit[2].m.n-1]), 0, 1);
    pars_add(freepars, "rho_solv_3", &(fit[3].m.rho[fit[3].m.n-1]), 0, 1);
    pars_add(freepars, "rho_solv_4", &(fit[4].m.rho[fit[4].m.n-1]), 0, 1);
    pars_add(freepars, "rho_solv_4", &(fit[5].m.rho[fit[5].m.n-1]), 0, 1);

  constraints = constr_models;
  output_model = save;
  return fit;
}
