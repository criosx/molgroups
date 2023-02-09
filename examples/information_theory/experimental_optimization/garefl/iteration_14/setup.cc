#include <iostream>
#include <cassert>
#include "setup.h"
#include "stdio.h"
#include "refl.h"
#include "reflcalc.h"
#include "molgroups.cc"

#define FWHM 2.354820045   // ga_refl uses FWHM. Divide by FWHM to get sigma units.

//reflectivity
#define MODELS 1

//canvas for continuous distribution model
#define CANVASSTART 1
#define DIMENSION 200
#define STEPSIZE 0.5

/* initialising non-standard fitting variables */
double aArea[DIMENSION], anSL[DIMENSION];
double background[MODELS];
char str2[2];

double normarea;
double l_lipid1, l_lipid2, vf_bilayer,  global_rough, rho_siox, l_siox, l_submembrane, sigma;

ssBLM_POPC  bilayer;

void fnSetContrastBilayer(double sigma, double global_rough, double rho_substrate, double rho_siox, double l_siox, double l_submembrane, double l_lipid1, double l_lipid2, double vf_bilayer, fitinfo *fit, int contraststart, int contrastend){
    
    double dMaxArea;
    int i;
    
    for (i=contraststart; i<contrastend+1; i++) {
        bilayer.fnSet(sigma, global_rough, rho_substrate, rho_siox, l_siox, l_submembrane, l_lipid1, l_lipid2, vf_bilayer, 0, 0);
        dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
        dMaxArea=bilayer.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
        normarea=dMaxArea;
        fnWriteCanvas2Model(aArea,anSL,fit,CANVASSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,i,i);
    }
}

/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
    int i;
    //int iPeaked,iMaxPeak, k, i2;
    //fitpars *pars = &fit[0].pars;
    
    
    /* Rescue the free parameters from the model. */
    //for (i=0; i < fit[1].pars.n; i++)
    //    fit[1].pars.value[i] = *(fit[1].pars.address[i]);
    
    /* Go through all layers copying parameters from model 0 to other models */
    //tied_parameters(fit);
    
    
    /* Restore the free parameters to the model. */
    //for (i=0; i < fit[1].pars.n; i++){
    //    *(fit[1].pars.address[i]) = fit[1].pars.value[i];
    //}
    
    for (i=0; i<MODELS; i++) {
        fit[i].beam.background=pow(10,background[i]);
    }
    
    
    //---- d31-POPC bilayer ----
    //----Neat -------
    fnSetContrastBilayer(sigma, global_rough, fit[0].m.rho[0], rho_siox, l_siox, l_submembrane, l_lipid1, l_lipid2, vf_bilayer, fit, 0, 0);
    
    
}

void save(fitinfo *fit)
{
    fnSetContrastBilayer(sigma, global_rough, fit[0].m.rho[0], rho_siox, l_siox, l_submembrane, l_lipid1, l_lipid2, vf_bilayer, fit, 0, 0);
    FILE *fp;
    fp=fopen("mol.dat","w");        
    bilayer.fnWriteGroup2File(fp,"bilayer",DIMENSION,STEPSIZE);
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
    fit_data(&fit[0],"sim0.dat"); /*neat */
    
    /* Initialize instrument parameters for each model.*/
    /* setup for NG7 */
    for (i=0; i < MODELS; i++) {
        
        const double L = 5.00,dLoL=0.015,d=1800.0;
        double Qlo, Tlo, dTlo,dToT,s1,s2;
        Qlo=0.008;
        Tlo=Q2T(L,Qlo);
        s1=0.1, s2=s1;
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
        for (j=0; j < DIMENSION; j++) {
            model_layer(&fit[i].m, STEPSIZE, 0.00e-6, 0.0e-8, 0);		/* Canvas */
        }
        model_layer(&fit[i].m, 100.000, 6.35e-6, 0.0e-8, 0.000);		/* Solvent */
    }
    
    /*correct solvent layers for different models*/
    /* fit[3].m.d[3] = ... */

    //headgroup.bWrapping=false;
    
    /*=============== FIT PARAMETERS ===============================*/
    
    /* Specify which parameters are your fit parameters. Parameters are fitted
     * to be the same in all datasets by default
     */
    
    
    pars_add(pars, "l_siox", &(l_siox), 10.0, 30.0);
    pars_add(pars, "rho_siox", &(rho_siox), 3.2e-06, 3.8e-06);
    pars_add(pars, "l_submembrane", &(l_submembrane), 1.0, 10.0);
    pars_add(pars, "l_lipid1", &(l_lipid1), 10.0, 15.0);
    pars_add(pars, "l_lipid2", &(l_lipid2), 10.0, 15.0);
    pars_add(pars, "vf_bilayer", &(vf_bilayer), 0.9, 1.0);
    
    pars_add(pars, "rho_solv_0", &(fit[0].m.rho[fit[0].m.n-1]), 5.96e-06, 6.559999999999999e-06);
    
    pars_add(pars, "global_rough", &(global_rough), 2.0, 5.0);
    pars_add(pars, "sigma",        &(sigma), 2.0, 5.0);
    
    pars_add(pars, "background_0", &(background[0]), -9.0, -5.0);
    
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
