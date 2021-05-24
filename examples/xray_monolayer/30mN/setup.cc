#include <iostream>
#include <cassert>
#include "stdio.h"
#include "refl.h"
#include "reflcalc.h"
#include "molgroups.cc"

#define FWHM 2.354820045   // ga_refl uses FWHM. Divide by FWHM to get sigma units.

#define MODELS 1
#define GAUSSSTART 1
#define DIMENSION 1000
#define STEPSIZE 0.5

/* initialising non-standard fitting variables */
double aArea[DIMENSION], anSL[DIMENSION], aAbsorb[DIMENSION];

double normarea, mu_organic; 
double vf_bilayer;
double sigma, global_rough;
double d_on1, d_on2, d_on3, d_on4, d_on5, d_on6, d_on7;
double d_on8, d_on9, d_on10, d_on11, d_on12, d_on13, d_on14, d_on15;
double vf_on1, vf_on2, vf_on3, vf_on4, vf_on5, vf_on6;
double vf_on7, vf_on8, vf_on9, vf_on10, vf_on11, vf_on12, vf_on13, vf_on14;
double l_lipid, penetration, rhoprot;

Monolayer_DPPS  bilayer_neat;  
Hermite         protein(16,0,0,0);


/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
    int bPeaked,i;
    double dMaxArea;
    fitpars *pars = &fit[0].pars;
    
    for (i=0; i<pars_count(pars); i++)
    {
        if (pars_peek(pars,i)==pars_max(pars,i))
        {pars_poke(pars,i,pars_min(pars,i)+0.99*(pars_max(pars,i)-pars_min(pars,i)));}
        if (pars_peek(pars,i)==pars_min(pars,i))
        {pars_poke(pars,i,pars_min(pars,i)+0.01*(pars_max(pars,i)-pars_min(pars,i)));}
    }
    
    /* Rescue the free parameters from the model. */
    //for (i=0; i < fit[1].pars.n; i++)
    //fit[1].pars.value[i] = *(fit[1].pars.address[i]);
    
    /* Go through all layers copying parameters from model 0 to other models */
    //tied_parameters(fit);
    
    /* copy the global roughness to all interfaces not SIOX interfaces*/
    // for (i=0; i< MODELS; i++) { //if(rough_au_cr>fit[0].m.d[2]) {rough_au_cr=fit[0].m.d[2];}
    //   for (k=1;k<GAUSSSTART; k++) fit[i].m.rough[k]=FWHM*global_rough;
    // }    
    
    /* Restore the free parameters to the model. */
    //for (i=0; i < fit[1].pars.n; i++){
    //  *(fit[1].pars.address[i]) = fit[1].pars.value[i];
    //}
    
    // neat bilayer
    
    bilayer_neat.sigma=sigma;
    bilayer_neat.global_rough=global_rough;
    bilayer_neat.vf_bilayer=vf_bilayer;
    bilayer_neat.rho_substrate=fit[0].m.rho[0];
    bilayer_neat.l_lipid=l_lipid;
    bilayer_neat.absorb=mu_organic;
    
    bilayer_neat.fnAdjustParameters();
    
    dMaxArea=fnClearCanvas(aArea, anSL, aAbsorb, DIMENSION);
    dMaxArea=bilayer_neat.fnWriteProfile(aArea, anSL, aAbsorb, DIMENSION, STEPSIZE, dMaxArea);
    normarea=dMaxArea;
    
    /*
    bPeaked=0;
    if (vf_on2<vf_on1*0.9) {vf_on2=vf_on1*0.9;}
    if (vf_on3<vf_on2*0.9) {
        bPeaked=1;
    }
    if (bPeaked==0) {
        if (vf_on4<vf_on3*0.9) {
            bPeaked=1;
        }
    }
    else {
        if (vf_on4*0.9>vf_on3) {
            vf_on4=vf_on3/0.9;
        }
    }
    if (bPeaked==0) {
        if (vf_on5<vf_on4*0.9) {
            bPeaked=1;
        }
    }
    else {
        if (vf_on5*0.9>vf_on4) {
            vf_on5=vf_on4/0.9;
        }
    }
    if (bPeaked==0) {
        if (vf_on6<vf_on5*0.9) {
            bPeaked=1;
        }
    }
    else {
        if (vf_on6*0.9>vf_on5) {
            vf_on6=vf_on5/0.9;
        }
    }
    */
    
    protein.fnSetNormarea(normarea);
    protein.fnSetnSLD(rhoprot);
    protein.absorb=mu_organic;
    protein.dp[0]=bilayer_neat.headgroup->fnGetZ()+0.5*9.56+penetration;
    protein.dp[1]=protein.dp[0]+d_on1;
    protein.dp[2]=protein.dp[1]+d_on2;
    protein.dp[3]=protein.dp[2]+d_on3;
    protein.dp[4]=protein.dp[3]+d_on4;
    protein.dp[5]=protein.dp[4]+d_on5;
    protein.dp[6]=protein.dp[5]+d_on6;
    protein.dp[7]=protein.dp[6]+d_on7;
    protein.dp[8]=protein.dp[7]+d_on8;
    protein.dp[9]=protein.dp[8]+d_on9;
    protein.dp[10]=protein.dp[9]+d_on10;
    protein.dp[11]=protein.dp[10]+d_on11;
    protein.dp[12]=protein.dp[11]+d_on12;
    protein.dp[13]=protein.dp[12]+d_on13;
    protein.dp[14]=protein.dp[13]+d_on14;
    protein.dp[15]=protein.dp[14]+d_on15;
    protein.vf[0]=0;
    protein.vf[1]=vf_on1;
    protein.vf[2]=vf_on2;
    protein.vf[3]=vf_on3;
    protein.vf[4]=vf_on4;
    protein.vf[5]=vf_on5;
    protein.vf[6]=vf_on6;
    protein.vf[7]=vf_on7;
    protein.vf[8]=vf_on8;
    protein.vf[9]=vf_on9;
    protein.vf[10]=vf_on10;
    protein.vf[11]=vf_on11;
    protein.vf[12]=vf_on12;
    protein.vf[13]=vf_on13;
    protein.vf[14]=vf_on14;
    protein.vf[15]=0;
    //printf("Absorption example %e \n", aAbsorb[152]);
    protein.fnOverlayProfile(aArea,anSL,aAbsorb,DIMENSION,STEPSIZE,dMaxArea);
    //printf("Absorption example %e \n", aAbsorb[152]);
    fnWriteCanvas2Model(aArea,anSL,aAbsorb,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,0,0);
    
    
}

void save(fitinfo *fit)
{
    
    FILE *fp;
    fp=fopen("mol.dat","w");        
    bilayer_neat.fnWritePar2File(fp,"bilayer",DIMENSION,STEPSIZE);
    protein.fnWritePar2File(fp,"protein",DIMENSION,STEPSIZE);
    fclose(fp);
    
}


/*============ INITIAL SETUP ============================*/
extern "C"
fitinfo* setup_models(int *models)
{
    static fitinfo fit[MODELS];
    int i,j;
    fitpars *pars = &fit[0].pars;
    //fitpars *freepars = &fit[1].pars;
    *models = MODELS;
    
    for (i=0; i < MODELS; i++) fit_init(&fit[i]);
    
    /* Load the data for each model */
    fit_data(&fit[0],"exp7dpps30mNpm1mMDTT500nMdp12MLV1stfull001.refl");
    
    /* Initialize instrument parameters for each model.*/
    /* setup for NG7 */
    for (i=0; i < MODELS; i++) {
        
        data_resolution_fixed(&fit[i].dataA,1.54,0.00015,0.,0.,0.0003);
        interface_create(&fit->rm, "erf", erf_interface, 21);
    }
    
    /*============= MODEL =====================================*/
    
    /* Add layers: d, rho, mu, rough */
    for (i=0; i < MODELS; i++) {
        model_layer(&fit[i].m, 0.00, 0.0e-6, 0.0e-8, 7.000);  	    /* 0 Air */
        for (j=0; j < DIMENSION; j++) {
            model_layer(&fit[i].m, STEPSIZE, 0.00e-6, 0.0e-8, 0);		/* Gaussian */
        }
        model_layer(&fit[i].m, 100.000, 6.35e-6, 0.0e-8, 0.000);		/* Solvent */
    }
    
    /*correct solvent layers for different models*/
    /* fit[3].m.d[3] = ... */
    
    bilayer_neat.sigma=2.;
    bilayer_neat.global_rough=3.;
    bilayer_neat.vf_bilayer=1.0;
    bilayer_neat.rho_substrate=fit[0].m.rho[0];
    bilayer_neat.l_lipid=10.;
    
    bilayer_neat.headgroup->nSL=4.54E-3;                            //jury-rigging for x-ray scattering lengths
    bilayer_neat.nslacyllipid=6.81E-3;
    bilayer_neat.nslmethyllipid=5.07E-4;
    bilayer_neat.absorbacyllipid=0;
    bilayer_neat.absorbmethyllipid=0;
    
    bilayer_neat.fnAdjustParameters();
    
    
    
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
    
    vf_bilayer=1;
    l_lipid=12;
    sigma=2.5;
    
    rhoprot=1.2e-5;
    
    //headgroup.bWrapping=false;
    
    /*=============== FIT PARAMETERS ===============================*/
    
    /* Specify which parameters are your fit parameters. Parameters are fitted
     * to be the same in all datasets by default
     */
    
    
    pars_add(pars, "l_lipid",             &(l_lipid), 8, 22);
    //pars_add(pars, "vf_bilayer",          &(vf_bilayer), 0.6, 1);  
    pars_add(pars, "sigma",               &(sigma), 2.0, 5.);
    pars_add(pars, "global_rough",        &(global_rough), 2.0, 5.);
    
    pars_add(pars, "rhoprot", &(rhoprot), 1e-5, 1.4e-5);
    pars_add(pars, "penetration",        &(penetration), -10, 10.);
    pars_add(pars, "d_on1", &(d_on1), 10, 35);
    pars_add(pars, "d_on2", &(d_on2), 10, 35);
    pars_add(pars, "d_on3", &(d_on3), 10, 35);
    pars_add(pars, "d_on4", &(d_on4), 10, 35);
    pars_add(pars, "d_on5", &(d_on5), 10, 35);
    pars_add(pars, "d_on6", &(d_on6), 10, 35);
    pars_add(pars, "d_on7", &(d_on7), 10, 35);
    pars_add(pars, "d_on8", &(d_on8), 10, 35);
    pars_add(pars, "d_on9", &(d_on9), 10, 45);
    pars_add(pars, "d_on10", &(d_on10), 10, 35);
    pars_add(pars, "d_on11", &(d_on11), 10, 35);
    pars_add(pars, "d_on12", &(d_on12), 10, 35);
    pars_add(pars, "d_on13", &(d_on13), 10, 35);
    pars_add(pars, "d_on14", &(d_on14), 10, 35);
    pars_add(pars, "d_on15", &(d_on15), 10, 35);
    pars_add(pars, "vf_on1", &(vf_on1), 0.00, 0.9);
    pars_add(pars, "vf_on2", &(vf_on2), 0.00, 0.9);
    pars_add(pars, "vf_on3", &(vf_on3), 0.00, 0.9);
    pars_add(pars, "vf_on4", &(vf_on4), 0.00, 0.9);
    pars_add(pars, "vf_on5", &(vf_on5), 0.00, 0.9);
    pars_add(pars, "vf_on6", &(vf_on6), 0.00, 0.9);
    pars_add(pars, "vf_on7", &(vf_on7), 0.00, 0.9);
    pars_add(pars, "vf_on8", &(vf_on8), 0.00, 0.9);
    pars_add(pars, "vf_on9", &(vf_on9), 0.00, 0.9);
    pars_add(pars, "vf_on10", &(vf_on10), 0.00, 0.9);
    pars_add(pars, "vf_on11", &(vf_on11), 0.00, 0.9);
    pars_add(pars, "vf_on12", &(vf_on12), 0.00, 0.9);
    pars_add(pars, "vf_on13", &(vf_on13), 0.00, 0.9);
    pars_add(pars, "vf_on14", &(vf_on14), 0.00, 0.9);
    
    pars_add(pars, "rho_solv_0", &(fit[0].m.rho[fit[0].m.n-1]), 1e-6, 20e-6);
    pars_add(pars, "mu_solv_0",  &(fit[0].m.mu[fit[0].m.n-1]), 9e-9, 1e-6);
    pars_add(pars, "mu_organic",  &(mu_organic), 1e-9, 1e-7);
    
    
    pars_add(pars, "beamintensity", &(fit[0].beam.intensity), 0.6, 1.1); 
    pars_add(pars, "background_0", &(fit[0].beam.background), 1e-9, 5e-7);
    
    /* Build a list of 'free parameters' in fit[1].pars. These are
     * parameters for which the values are aloowed to differ from those
     * in model 0.  By default all values in all models are the same unless 
     * specified here. The range data is not used here, so set it to [0,1].  
     */
    
    //pars_add(freepars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), 0, 1);
    
    constraints = constr_models;
    output_model = save;
    return fit;
}
