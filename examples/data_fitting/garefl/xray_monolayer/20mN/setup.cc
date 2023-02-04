#include <iostream>
#include <cassert>
#include "stdio.h"
#include "refl.h"
#include "reflcalc.h"
#include "setup.h"
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
double sigma, global_rough, rho_solv_0;
double d_on1, d_on2, d_on3, d_on4, d_on5, d_on6, d_on7;
double d_on8, d_on9, d_on10, d_on11, d_on12, d_on13, d_on14, d_on15;
double d_on16, d_on17, d_on18, d_on19, d_on20, d_on21, d_on22;
double d_on23, d_on24;
double vf_on1, vf_on2, vf_on3, vf_on4, vf_on5, vf_on6;
double vf_on7, vf_on8, vf_on9, vf_on10, vf_on11, vf_on12, vf_on13, vf_on14;
double vf_on15, vf_on16, vf_on17, vf_on18, vf_on19, vf_on20, vf_on21, vf_on22;
double vf_on23;
double l_lipid, penetration, rhoprot;


Monolayer_DPPS_xray  bilayer_neat;  
Hermite              protein(25,0,0,0);

int fnPeakAdjust(double &value1, double &value2, int &iPeaked, int iMaxPeak, int iLastTrendUp) 
{
    int returnvalue;
    const double tolerance=0.9;
    
    if ((value2>=value1*tolerance) && (value2<=value1/tolerance)) {returnvalue=iLastTrendUp;}
    else if (value1<value2) {returnvalue=1;}
    else {returnvalue=0;}
    
    if (iPeaked<iMaxPeak) {if ((iLastTrendUp==1)&&(value2<value1*tolerance)) {iPeaked+=1;}}
    else {if (value2>value1/tolerance) {value2=value1/tolerance;}};
    
    return returnvalue;
}


/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
    
    
    int iPeaked,iMaxPeak,i, i2;
    double dMaxArea,v;
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
    
    fit[0].m.rho[fit[0].m.n-1]=rho_solv_0;
    
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
    

    iPeaked=0; iMaxPeak=4;
    i2=fnPeakAdjust(vf_on1,vf_on2,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on2,vf_on3,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on3,vf_on4,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on4,vf_on5,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on5,vf_on6,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on6,vf_on7,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on7,vf_on8,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on8,vf_on9,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on9,vf_on10,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on10,vf_on11,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on11,vf_on12,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on12,vf_on13,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on13,vf_on14,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on14,vf_on15,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on15,vf_on16,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on16,vf_on17,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on17,vf_on18,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on18,vf_on19,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on19,vf_on20,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on20,vf_on21,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on21,vf_on22,iPeaked,iMaxPeak,i2);
    i2=fnPeakAdjust(vf_on22,vf_on23,iPeaked,iMaxPeak,i2);
    


    const double Spacing =14.0;
    
    protein.fnSetNormarea(normarea);
    protein.fnSetnSLD(rhoprot);
    protein.absorb=mu_organic;
    protein.dp[0]=0;
    protein.dp[1]=protein.dp[0]+1*Spacing+d_on1;
    protein.dp[2]=protein.dp[0]+2*Spacing+d_on2;
    protein.dp[3]=protein.dp[0]+3*Spacing+d_on3;
    protein.dp[4]=protein.dp[0]+4*Spacing+d_on4;
    protein.dp[5]=protein.dp[0]+5*Spacing+d_on5;
    protein.dp[6]=protein.dp[0]+6*Spacing+d_on6;
    protein.dp[7]=protein.dp[0]+7*Spacing+d_on7;
    protein.dp[8]=protein.dp[0]+8*Spacing+d_on8;
    protein.dp[9]=protein.dp[0]+9*Spacing+d_on9;
    protein.dp[10]=protein.dp[0]+10*Spacing+d_on10;
    protein.dp[11]=protein.dp[0]+11*Spacing+d_on11;
    protein.dp[12]=protein.dp[0]+12*Spacing+d_on12;
    protein.dp[13]=protein.dp[0]+13*Spacing+d_on13;
    protein.dp[14]=protein.dp[0]+14*Spacing+d_on14;
    protein.dp[15]=protein.dp[0]+15*Spacing+d_on15;
    protein.dp[16]=protein.dp[0]+16*Spacing+d_on16;
    protein.dp[17]=protein.dp[0]+17*Spacing+d_on17;
    protein.dp[18]=protein.dp[0]+18*Spacing+d_on18;
    protein.dp[19]=protein.dp[0]+19*Spacing+d_on19;
    protein.dp[20]=protein.dp[0]+20*Spacing+d_on20;
    protein.dp[21]=protein.dp[0]+21*Spacing+d_on21;
    protein.dp[22]=protein.dp[0]+22*Spacing+d_on22;
    protein.dp[23]=protein.dp[0]+23*Spacing+d_on23;
    protein.dp[24]=protein.dp[0]+24*Spacing;
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
    protein.vf[15]=vf_on15;
    protein.vf[16]=vf_on16;
    protein.vf[17]=vf_on17;
    protein.vf[18]=vf_on18;
    protein.vf[19]=vf_on19;
    protein.vf[20]=vf_on20;
    protein.vf[21]=vf_on21;
    protein.vf[22]=vf_on22;
    protein.vf[23]=vf_on23;
    protein.vf[24]=0;
    
    //calculate ratio of protein in hydrocarbon over total hydrocarbon ratio which will be the
    //fraction by which the headgroup density will be reduced
    v=protein.fnGetVolume(bilayer_neat.methyl->z-0.5*bilayer_neat.methyl->l,bilayer_neat.lipid->z+0.5*bilayer_neat.lipid->l)/(bilayer_neat.lipid->vol+bilayer_neat.methyl->vol);
    if (v<0) {v=0;}
    bilayer_neat.hc_substitution=v;
    
    dMaxArea=fnClearCanvas(aArea, anSL, aAbsorb, DIMENSION);
    dMaxArea=bilayer_neat.fnWriteProfile(aArea, anSL, aAbsorb, DIMENSION, STEPSIZE, dMaxArea);
    normarea=dMaxArea;
                                                                                                                              
    protein.fnOverlayProfile(aArea,anSL,aAbsorb,DIMENSION,STEPSIZE,dMaxArea);
    fnWriteCanvas2Model(aArea,anSL,aAbsorb,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,0,0);
    
    
    
}

void save(fitinfo *fit)
{
    
    FILE *fp;
    fp=fopen("mol.dat","w");        
    bilayer_neat.fnWriteGroup2File(fp,"bilayer",DIMENSION,STEPSIZE);
    protein.fnWriteGroup2File(fp,"protein",DIMENSION,STEPSIZE);
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
    fit_data(&fit[0],"exp7dpps20mNpm1mMDTT500nMdp12MLV5thfull.refl");
    
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
    
    
    pars_add(pars, "l_lipid",      &(l_lipid), 8, 22);
    //pars_add(pars, "vf_bilayer",          &(vf_bilayer), 0.6, 1);
    //pars_add(pars, "l_headgroup",  &(bilayer_neat.headgroup->l), 5, 18);
    //pars_add(pars, "vol_headgroup",  &(bilayer_neat.headgroup->vol), 200, 350);
    
    pars_add(pars, "sigma",               &(sigma), 0.5, 10.);
    pars_add(pars, "global_rough",        &(global_rough), 0.5, 10.);
    
    pars_add(pars, "rhoprot", &(rhoprot), 1e-5, 1.4e-5);
    //pars_add(pars, "penetration",        &(penetration), -12, 10.);
    pars_add(pars, "d_on1", &(d_on1), -5, 5);
    pars_add(pars, "d_on2", &(d_on2), -5, 5);
    pars_add(pars, "d_on3", &(d_on3), -5, 7);
    pars_add(pars, "d_on4", &(d_on4), -5, 5);
    pars_add(pars, "d_on5", &(d_on5), -5, 5);
    pars_add(pars, "d_on6", &(d_on6), -5, 5);
    pars_add(pars, "d_on7", &(d_on7), -5, 5);
    pars_add(pars, "d_on8", &(d_on8), -5, 5);
    pars_add(pars, "d_on9", &(d_on9), -5, 5);
    pars_add(pars, "d_on10", &(d_on10), -5, 5);
    pars_add(pars, "d_on11", &(d_on11), -5, 5);
    pars_add(pars, "d_on12", &(d_on12), -5, 5);
    pars_add(pars, "d_on13", &(d_on13), -5, 5);
    pars_add(pars, "d_on14", &(d_on14), -5, 5);
    pars_add(pars, "d_on15", &(d_on15), -5, 5);
    pars_add(pars, "d_on16", &(d_on16), -5, 5);
    pars_add(pars, "d_on17", &(d_on17), -5, 5);
    pars_add(pars, "d_on18", &(d_on18), -5, 5);
    pars_add(pars, "d_on19", &(d_on19), -5, 5);
    pars_add(pars, "d_on20", &(d_on20), -5, 5);
    pars_add(pars, "d_on21", &(d_on21), -5, 5);
    pars_add(pars, "d_on22", &(d_on22), -5, 5);
    pars_add(pars, "d_on23", &(d_on23), -5, 5);
    //pars_add(pars, "d_on24", &(d_on24), -5, 5);
    pars_add(pars, "vf_on1", &(vf_on1), 0.00, 0.9);
    pars_add(pars, "vf_on2", &(vf_on2), 0.00, 0.9);
    pars_add(pars, "vf_on3", &(vf_on3), 0.10, 0.9);
    pars_add(pars, "vf_on4", &(vf_on4), 0.10, 0.9);
    pars_add(pars, "vf_on5", &(vf_on5), 0.10, 0.9);
    pars_add(pars, "vf_on6", &(vf_on6), 0.10, 0.9);
    pars_add(pars, "vf_on7", &(vf_on7), 0.10, 0.9);
    pars_add(pars, "vf_on8", &(vf_on8), 0.10, 0.9);
    pars_add(pars, "vf_on9", &(vf_on9), 0.00, 0.9);
    pars_add(pars, "vf_on10", &(vf_on10), 0.00, 0.9);
    pars_add(pars, "vf_on11", &(vf_on11), 0.00, 0.9);
    pars_add(pars, "vf_on12", &(vf_on12), 0.00, 0.5);
    pars_add(pars, "vf_on13", &(vf_on13), 0.00, 0.5);
    pars_add(pars, "vf_on14", &(vf_on14), 0.00, 0.5);
    pars_add(pars, "vf_on15", &(vf_on15), 0.00, 0.5);
    pars_add(pars, "vf_on16", &(vf_on16), 0.00, 0.5);
    pars_add(pars, "vf_on17", &(vf_on17), 0.00, 0.5);
    pars_add(pars, "vf_on18", &(vf_on18), 0.00, 0.5);
    pars_add(pars, "vf_on19", &(vf_on19), 0.00, 0.5);
    pars_add(pars, "vf_on20", &(vf_on20), 0.00, 0.5);
    pars_add(pars, "vf_on21", &(vf_on21), 0.00, 0.5);
    pars_add(pars, "vf_on22", &(vf_on22), 0.00, 0.5);
    pars_add(pars, "vf_on23", &(vf_on23), 0.00, 0.5);
    
    pars_add(pars, "rho_solv_0", &(rho_solv_0), 1e-6, 20e-6);
    //pars_add(pars, "rho_solv_0", &(fit[0].m.rho[fit[0].m.n-1]), 1e-6, 20e-6);
    pars_add(pars, "mu_solv_0",  &(fit[0].m.mu[fit[0].m.n-1]), 9e-9, 1e-6);
    pars_add(pars, "mu_organic",  &(mu_organic), 1e-9, 1e-7);
    //pars_add(pars, "mu_air",  &(fit[0].m.mu[0]), 1e-12, 1e-9);
    
    
    pars_add(pars, "beamintensity", &(fit[0].beam.intensity), 0.9, 1.1); 
    pars_add(pars, "alignment",     &(fit[0].beam.alignment), -0.000, 0.007);
    pars_add(pars, "background_0", &(fit[0].beam.background), 1e-10, 5e-8);
    
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
