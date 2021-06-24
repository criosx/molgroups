#include <iostream>
#include <cassert>
#include <string>
#include <sstream>
#include "setup.h"
#include "stdio.h"
#include "refl.h"
#include "reflcalc.h"
#include "molgroups.cc"
#include <xmmintrin.h>


using namespace std;


#define FWHM 2.354820045   // ga_refl uses FWHM. Divide by FWHM to get sigma units.

//reflectivity
#define MODELS 2

//canvas for continuous distribution model
#define GAUSSSTART 4
#define DIMENSION 200
#define STEPSIZE 0.5


#define PROTDEUT 2.78e-6
#define PROTNONDEUT 1.67e-6
#define PROTDEUT2 6.7754e-6
#define PROTNONDEUT2 5.3137e-6
#define NSLDH2O -0.5666e-6
#define NSLDD2O 6.36e-6

/* initialising non-standard fitting variables */
double aArea[DIMENSION], anSL[DIMENSION];
double background[MODELS];
char str2[2];

double normarea, protexchratio, fraction_rinse1, fraction_rinse2, fraction_inc2; 
double l_lipid1, l_lipid2, radius_defects, thetaoffset;
double gcase_vf1, gcase_vf2, gcrinse1_vf1, gcrinse1_vf2;
double vf_bilayer,  global_rough, rough_cr_au, mult_tether, l_tether, rho_h2o;
double sigma, nf_tether, nf_lipid_2, penetration;
double vf_bilayer_gcase, dl_lipid_gcase, dl_lipid_inc2;
double vf_bilayer_gcrinse1, dl_lipid_gcrinse1;

tBLM_HC18_POPC_POPS  bilayer;

void fnSetContrastBilayer(double sigma, double global_rough, double rho_substrate, double nf_tether, double mult_tether, double l_tether, double l_lipid1, double l_lipid2, double vf_bilayer, double nf_lipid_2, double nf_lipid3, double nf_chol, double radius_defect, fitinfo *fit, int contraststart, int contrastend){
    
    double dMaxArea;
    int i;
    
    for (i=contraststart; i<contrastend+1; i++) {
        //printf("Enter fnSetContrastBilayer \n");
        bilayer.fnSet(sigma, global_rough, rho_substrate, fit[i].m.rho[fit[i].m.n-1], nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, nf_lipid3, nf_chol, 0, 0, radius_defect);
        dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
        //printf("Midpoint fnSetContrastBilayer \n");
        dMaxArea=bilayer.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
        normarea=dMaxArea;
        //printf("Midpoint 2 fnSetContrastBilayer \n");
        fnWriteCanvas2Model(aArea,anSL,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,contraststart,contrastend);
        //printf("Exit fnSetContrastBilayer \n");
    }
}


/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
    int iPeaked,iMaxPeak, i, k, i2, wait;
    fitpars *pars = &fit[0].pars;

    //printf("Enter Constraints \n");

    
    //for (i=0; i<CONTROLPOINTS; i++) {
    //    printf("point %i pos %e vf %e frac %e \n", i, dp_on[i], vf_on[i], frac2_on[i]);
    //    printf("sigma %e global_rough %e nf_tether % e \n", sigma, global_rough, nf_tether);
    //}
    
    
    /*for (i=0; i<pars_count(pars); i++)
    {
        if (pars_peek(pars,i)==pars_max(pars,i))
        {pars_poke(pars,i,pars_min(pars,i)+0.9999*(pars_max(pars,i)-pars_min(pars,i)));}
        if (pars_peek(pars,i)==pars_min(pars,i))
        {pars_poke(pars,i,pars_min(pars,i)+0.0001*(pars_max(pars,i)-pars_min(pars,i)));}
    }
    */

    /* Rescue the free parameters from the model. */
    for (i=0; i < fit[1].pars.n; i++)
        fit[1].pars.value[i] = *(fit[1].pars.address[i]);
    
    /* Go through all layers copying parameters from model 0 to other models */
    tied_parameters(fit);
    
    /* copy the global roughness to all interfaces not SIOX interfaces*/
    for (i=0; i< MODELS; i++) { //if(rough_au_cr>fit[0].m.d[2]) {rough_au_cr=fit[0].m.d[2];}
        for (k=3;k<GAUSSSTART; k++) fit[i].m.rough[k]=FWHM*global_rough;
        fit[i].m.rough[3]=rough_cr_au*FWHM;
        fit[i].beam.alignment=thetaoffset;
    }    
    
    /* Restore the free parameters to the model. */
    for (i=0; i < fit[1].pars.n; i++){
        *(fit[1].pars.address[i]) = fit[1].pars.value[i];
    }  
    
    for (i=0; i<MODELS; i++) {
        fit[i].beam.background=pow(10,background[i]);
    }
    
    //iPeaked=0; iMaxPeak=1, i2=1; vf_on[0]=0;
    //for (i=1; i<CONTROLPOINTS; i++) {
    //    i2=fnPeakAdjust(vf_on[i-1],vf_on[i],iPeaked,iMaxPeak,i2);
    //}
    /*
    iPeaked=0; iMaxPeak=1, i2=1; frac2_on[0]=0;
    for (i=1; i<CONTROLPOINTS; i++) {
        i2=fnPeakAdjust(frac2_on[i-1],frac2_on[i],iPeaked,iMaxPeak,i2);
    }
    */
    //printf("Midpoint Constraints 1 \n");
    
    
    //----Neat -------
    //fnSetContrastBilayer(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, 0.0, 0.0, radius_defects, fit, 0, 2);
    fnSetContrastBilayer(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, 0.0, 0.0, radius_defects, fit, 0, 1);
    mult_tether=bilayer.mult_tether;
    l_tether=bilayer.l_tether;
    //printf("Midpoint Constraints 2 \n");
    
    
    
    
}

void save(fitinfo *fit)
{
    fnSetContrastBilayer(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, 0.0, 0.0, radius_defects, fit, 0, 1);
   
    FILE *fp;
    fp=fopen("mol.dat","w");        
    bilayer.fnWritePar2File(fp,"bilayer",DIMENSION,STEPSIZE);
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
    fit_data(&fit[0],"os060.refl"); /*neat */
    fit_data(&fit[1],"os061.refl");
    
    /* Initialize instrument parameters for each model.*/
    /* setup for NG7 */
    for (i=0; i < MODELS; i++) {
        
        const double L = 5.00,dLoL=0.015,d=1800.0;
        double Qlo, Tlo, dTlo,dToT,s1,s2;
        Qlo=0.008;
        Tlo=Q2T(L,Qlo);
        s1=0.05, s2=s1;
        dTlo=resolution_dT(s1,s2,d);
        dToT=resolution_dToT(s1,s2,d,Tlo);
        data_resolution_fv(&fit[i].dataA,L,dLoL,Qlo,dTlo,dToT);
        fit[i].beam.lambda = L;
        interface_create(&fit[i].rm, "erf", erf_interface, 21);
    }
    
    /*============= MODEL =====================================*/
    
    /* Add layers: d, rho, mu, rough */
    for (i=0; i < MODELS; i++) {
        model_layer(&fit[i].m, 0.00000, 2.07e-6, 0.0e-8, 3.000);  	/* 0 Si */
        model_layer(&fit[i].m, 20.8, 3.55e-6, 0.0e-8, 3.000);			/* 1 oxide */
        model_layer(&fit[i].m, 20.8, 3.22e-6, 0.0e-8, 3.000);			/* 2 Cr */
        model_layer(&fit[i].m, 20.8, 4.45e-6, 0.0e-8, 7.000);			/* 3 Au */
        for (j=0; j < DIMENSION; j++) {
            model_layer(&fit[i].m, STEPSIZE, 0.00e-6, 0.0e-8, 0);		/* Gaussian */
        }
        model_layer(&fit[i].m, 100.000, 6.35e-6, 0.0e-8, 0.000);		/* Solvent */
    }
    
    /*correct solvent layers for different models*/
    /* fit[3].m.d[3] = ... */
    
   //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    
    nf_lipid_2=0.0;
    radius_defects=100;
    
    //headgroup.bWrapping=false;
    
    /*=============== FIT PARAMETERS ===============================*/
    
    /* Specify which parameters are your fit parameters. Parameters are fitted
     * to be the same in all datasets by default
     */
    
    
    pars_add(pars, "d_oxide", &(fit[0].m.d[1]), 5., 60);
    pars_add(pars, "d_Cr", &(fit[0].m.d[2]), 10, 150.);
    pars_add(pars, "d_gold", &(fit[0].m.d[3]), 60., 170.);
    
    pars_add(pars, "l_tether", &(l_tether), 6, 18);
    pars_add(pars, "l_lipid1", &(l_lipid1), 8, 21);
    pars_add(pars, "l_lipid2", &(l_lipid2), 9, 16);
    
    pars_add(pars, "nf_tether", &(nf_tether), 0.2, 1.0);
    pars_add(pars, "mult_tether", &(mult_tether), 0.1, 4.);
    //pars_add(pars, "nf_lipid_2", &(nf_lipid_2), 0.4, 0.6);
    
    pars_add(pars, "vf_bilayer", &(vf_bilayer), 0.80, 1);
    
    pars_add(pars, "rho_SiOx", &(fit[0].m.rho[1]), 3.0e-6, 3.8e-6);
    pars_add(pars, "rho_Cr", &(fit[0].m.rho[2]), 2.7e-6, 4.8e-6);
    pars_add(pars, "rho_Au", &(fit[0].m.rho[3]), 4.2e-6, 4.8e-6);
    
    pars_add(pars, "rho_solv_0", &(fit[0].m.rho[fit[0].m.n-1]), 5.7e-6, 6.4e-6);
    pars_add(pars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), -0.56e-6, 0.5e-6);
    
    pars_add(pars, "global_rough", &(global_rough), 2.0, 8.);
    pars_add(pars, "rough_cr_au",  &(rough_cr_au),  2.0, 14.0);
    pars_add(pars, "sigma",        &(sigma),        2.0, 5.);
    pars_add(pars, "thetaoffset0", &(thetaoffset), -0.02, 0.02);
    
    pars_add(pars, "background_0", &(background[0]), -9, -5);
    pars_add(pars, "background_1", &(background[1]), -9, -5);
    
    /* Build a list of 'free parameters' in fit[1].pars. These are
     * parameters for which the values are aloowed to differ from those
     * in model 0.  By default all values in all models are the same unless 
     * specified here. The range data is not used here, so set it to [0,1].  
     */
    
    pars_add(freepars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), 0, 1);
    
    constraints = constr_models;
    output_model = save;
    return fit;
}
