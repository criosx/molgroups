#include <iostream>
#include <cassert>
#include "setup.h"
#include "stdio.h"
#include "refl.h"
#include "reflcalc.h"
#include "molgroups.cc"
#include <xmmintrin.h>

#define FWHM 2.354820045   // ga_refl uses FWHM. Divide by FWHM to get sigma units.

//reflectivity
#define MODELS 5

//canvas for continuous distribution model
#define GAUSSSTART 4
#define DIMENSION 500
#define STEPSIZE 0.5

//spline
#define CONTROLPOINTS 9
#define SPACING 20.0
#define PENETRATION 2

//nSL normalization
#define DEUTERATIONFRACTION 0.61
#define DEUTERATIONYIELD 0.65


#define PROTDEUT 3.02e-6
#define PROTNONDEUT 1.82e-6
#define PROTDEUT2 7.35e-6
#define PROTNONDEUT2 6.16e-6
#define NSLDH2O -0.5666e-6
#define NSLDD2O 6.36e-6

/* initialising non-standard fitting variables */
double aArea[DIMENSION], anSL[DIMENSION];
double background[MODELS];
char str2[2];

double normarea, protexchratio, fraction_inc1, fraction_inc2, fraction_inc3;
double l_lipid1, l_lipid2, radius_defects, thetaoffset;
double gcase_vf1, gcase_vf2, gcrinse1_vf1, gcrinse1_vf2;
double vf_bilayer,  global_rough, rough_cr_au, mult_tether, l_tether, rho_h2o;
double sigma, nf_tether, nf_lipid_2, penetration;
double vf_bilayer_gcase, dl_lipid_gcase, dl_lipid_inc2;
double vf_bilayer_gcrinse1, dl_lipid_gcrinse1;

tBLM_HC18_POPC_POPS  bilayer;
SLDHermite          protein(CONTROLPOINTS,0,0);
double dp_on[CONTROLPOINTS], vf_on[CONTROLPOINTS], sld_on[CONTROLPOINTS], frac2_on[CONTROLPOINTS];

void fnSetContrastBilayer(double sigma, double global_rough, double rho_substrate, double nf_tether, double mult_tether, double l_tether, double l_lipid1, double l_lipid2, double vf_bilayer, double nf_lipid_2, double nf_lipid3, double nf_chol, double radius_defect, fitinfo *fit, int contraststart, int contrastend){
    
    double dMaxArea;
    
    bilayer.fnSet(sigma, global_rough, rho_substrate, nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, nf_lipid3, nf_chol, 0, 0, radius_defect);
    dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
    dMaxArea=bilayer.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
    normarea=dMaxArea;
    fnWriteCanvas2Model(aArea,anSL,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,contraststart,contrastend);
}

void fnSetContrastProtein(double sigma, double global_rough, double rho_substrate, double nf_tether, double mult_tether, double l_tether, double l_lipid1, double l_lipid2, double vf_bilayer, double nf_lipid_2, double nf_lipid_3, double nf_chol, double radius_defect, double nf_protein, double protexchratio, double protnondeut, double protdeut, fitinfo *fit, int contraststart, int contrastend){
    
    double nSLtotaldiff,nSLnormfactor, nSLmaxdiff;
    double v1,v2,dMaxArea, pdeut, pnondeut;
    bool b;
    int i,j;
    
    
    //nSL normalization
    //printf("Setting up the spline ... \n");
    nSLnormfactor=1;
    for (j=0; j<CONTROLPOINTS; j++) {
        pnondeut=frac2_on[j]*PROTNONDEUT2+(1-frac2_on[j])*PROTNONDEUT;
        sld_on[j]=pnondeut;
        //printf("j %i pnondeut %g vf[j] %g frac2_on[j] %g sld_on[j] %g \n", j, pnondeut, vf_on[j], frac2_on[j], sld_on[j]);
    }
    bilayer.fnSet(sigma, global_rough, rho_substrate, nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, nf_lipid_3, nf_chol, 0,0, radius_defect);
    protein.fnSetNormarea(normarea);
    protein.fnSetRelative(SPACING,bilayer.headgroup2->fnGetZ()+0.5*9.56-PENETRATION*SPACING,dp_on,vf_on,sld_on,1);
    nSLtotaldiff=protein.fnGetnSL(protein.dp[0],protein.dp[CONTROLPOINTS-1])-PROTNONDEUT*protein.fnGetVolume(protein.dp[0],protein.dp[CONTROLPOINTS-1]);
    nSLmaxdiff=DEUTERATIONFRACTION*DEUTERATIONYIELD*(PROTNONDEUT2-PROTNONDEUT)*protein.fnGetVolume(protein.dp[0],protein.dp[CONTROLPOINTS-1]);
    if (nSLtotaldiff>0) {        
        nSLnormfactor=nSLmaxdiff/nSLtotaldiff;
    }
    //printf("Final integration results ... \n");
    //printf("nSL %g nSLprot %g nSLtotaldiff %g nSLmaxdiff %g nSLnormfactor %g volume %g \n", protein.fnGetnSL(protein.dp[0],protein.dp[CONTROLPOINTS-1]), PROTNONDEUT*protein.fnGetVolume(protein.dp[0],protein.dp[CONTROLPOINTS-1]), nSLtotaldiff, nSLmaxdiff, nSLnormfactor, protein.fnGetVolume(protein.dp[0],protein.dp[CONTROLPOINTS-1]));
    
    
    for (i=contraststart; i<contrastend+1; i++) {
        bilayer.fnSet(sigma, global_rough, rho_substrate, nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, nf_lipid_3, nf_chol, 0,0, radius_defect);
        protein.fnSetNormarea(normarea);
        for (j=0; j<CONTROLPOINTS; j++) {
            pnondeut=frac2_on[j]*nSLnormfactor*PROTNONDEUT2+(1-frac2_on[j]*nSLnormfactor)*PROTNONDEUT;
            pdeut=frac2_on[j]*nSLnormfactor*PROTDEUT2+(1-frac2_on[j]*nSLnormfactor)*PROTDEUT;
            sld_on[j]=pnondeut+protexchratio*(fit[i].m.rho[fit[i].m.n-1]-NSLDH2O)/(NSLDD2O-NSLDH2O)*(pdeut-pnondeut);
            //printf("j %i pnondeut %g vf[j] %g frac2_on[j] %g sld_on[j] %g \n", j, pnondeut, vf_on[j], frac2_on[j], sld_on[j]);
        }
        protein.fnSetRelative(SPACING,bilayer.headgroup2->fnGetZ()+0.5*9.56-PENETRATION*SPACING,dp_on,vf_on,sld_on,nf_protein);
        v1=protein.fnGetVolume(bilayer.lipid1->z-0.5*bilayer.lipid1->l,bilayer.methyl1->z+0.5*bilayer.methyl1->l)/((bilayer.lipid1->vol+bilayer.methyl1->vol)*bilayer.vf_bilayer);
        if (v1<0) {v1=0;}
        v2=protein.fnGetVolume(bilayer.methyl2->z-0.5*bilayer.methyl2->l,bilayer.lipid2->z+0.5*bilayer.lipid2->l)/((bilayer.lipid2->vol+bilayer.methyl2->vol)*bilayer.vf_bilayer);
        if (v2<0) {v2=0;}
        bilayer.fnSet(sigma, global_rough, rho_substrate, nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, nf_lipid_3, nf_chol, v1, v2, radius_defect);
        dMaxArea=fnClearCanvas(aArea, anSL, DIMENSION);
        dMaxArea=bilayer.fnWriteProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
        normarea=dMaxArea;
        protein.fnOverlayProfile(aArea, anSL, DIMENSION, STEPSIZE, dMaxArea);
        fnWriteCanvas2Model(aArea,anSL,fit,GAUSSSTART,DIMENSION,STEPSIZE,dMaxArea,normarea,i ,i );
    }
}

int fnPeakAdjust(double &value1, double &value2, int &iPeaked, int iMaxPeak, int iLastTrendUp)
{
    int returnvalue;
    const double tolerance=0.98;
    
    if ((fabs(value2)>=fabs(value1)*tolerance) && (fabs(value2)<=fabs(value1)/tolerance)) {returnvalue=iLastTrendUp;}
    else if (fabs(value1)<fabs(value2)) {returnvalue=1;}
    else {returnvalue=0;}
    
    if (iPeaked<iMaxPeak) {
        if ((iLastTrendUp==1)&&(fabs(value2)<fabs(value1)*tolerance)) {iPeaked+=1;}
    }
    else {if (fabs(value2)>fabs(value1)/tolerance) {value2=fabs(value1)/tolerance;}};
    //printf("value 1 %g value 2 %g iPeaked %i \n", fabs(value1), fabs(value2), iPeaked);
    
    return returnvalue;
}

/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
    int iPeaked,iMaxPeak, i, k, i2;
    fitpars *pars = &fit[0].pars;
    
    for (i=0; i<pars_count(pars); i++)
    {
        if (pars_peek(pars,i)==pars_max(pars,i))
        {pars_poke(pars,i,pars_min(pars,i)+0.9999*(pars_max(pars,i)-pars_min(pars,i)));}
        if (pars_peek(pars,i)==pars_min(pars,i))
        {pars_poke(pars,i,pars_min(pars,i)+0.0001*(pars_max(pars,i)-pars_min(pars,i)));}
    }
    
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
    /*iPeaked=0; iMaxPeak=1, i2=1;
    for (i=1; i<CONTROLPOINTS; i++) {
        i2=fnPeakAdjust(frac2_on[i-1],frac2_on[i],iPeaked,iMaxPeak,i2);
    }
    */
    
    //----Neat -------
    fnSetContrastBilayer(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2, 0.0, 0.0, radius_defects, fit, 0, 1);
    mult_tether=bilayer.mult_tether;
    l_tether=bilayer.l_tether;
    
    
    
    
    //------------- first rinse gcrinse1 ----------------------------------------------------------------
    
    fnSetContrastProtein(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1+dl_lipid_gcrinse1, l_lipid2+dl_lipid_gcrinse1, vf_bilayer_gcrinse1, nf_lipid_2, 0.0, 0.0, radius_defects, 1, protexchratio, PROTNONDEUT, PROTDEUT, fit, 2, 2);
    
    //------------- second incr----------------------------------------------------------------
    
    fnSetContrastProtein(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1+dl_lipid_gcrinse1, l_lipid2+dl_lipid_gcrinse1, vf_bilayer_gcrinse1, nf_lipid_2, 0.0, 0.0, radius_defects, fraction_inc2, protexchratio, PROTNONDEUT, PROTDEUT, fit, 3, 3);
    
    
    //--------------third inc-----------------------------------------------------------------------------
    fnSetContrastProtein(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1+dl_lipid_gcrinse1, l_lipid2+dl_lipid_gcrinse1, vf_bilayer_gcrinse1, nf_lipid_2, 0.0, 0.0, radius_defects, fraction_inc3, protexchratio, PROTNONDEUT, PROTDEUT, fit, 4, 4);
    
    
}

void save(fitinfo *fit)
{
    fnSetContrastProtein(sigma, global_rough, fit[0].m.rho[3], nf_tether, mult_tether, l_tether, l_lipid1+dl_lipid_gcrinse1, l_lipid2+dl_lipid_gcrinse1, vf_bilayer_gcrinse1, nf_lipid_2, 0.0, 0.0, radius_defects, 1, protexchratio, PROTNONDEUT, PROTDEUT, fit, 2, 2);
    
    FILE *fp;
    fp=fopen("mol.dat","w");        
    bilayer.fnWritePar2File(fp,"bilayer",DIMENSION,STEPSIZE);
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
    fitpars *freepars = &fit[1].pars;
    *models = MODELS;
    
    for (i=0; i < MODELS; i++) fit_init(&fit[i]);
    
    /* Load the data for each model */
    fit_data(&fit[0],"sy089.refl"); /*neat */
    fit_data(&fit[1],"sy090.refl");
    fit_data(&fit[2],"sy093.refl");
    fit_data(&fit[3],"sy094.refl"); //second incubation
    fit_data(&fit[4],"sy095.refl"); //second incubation
    
    /* Initialize instrument parameters for each model.*/
    /* setup for NG7 */
    for (i=0; i < MODELS; i++) {
        
        const double L = 5,dLoL=0.02,d=1700.0;
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
    
    protexchratio=0.7;
    
    rho_h2o=-0.5666666e-6;
    
    dp_on[0]=0;
    dp_on[CONTROLPOINTS-2]=0;
    dp_on[CONTROLPOINTS-1]=0;
    vf_on[0]=0;
    vf_on[CONTROLPOINTS-2]=0;
    vf_on[CONTROLPOINTS-1]=0;
    frac2_on[0]=1;
    frac2_on[CONTROLPOINTS-1]=1;
    
    nf_lipid_2=0.5;
    radius_defects=100;
    
    //headgroup.bWrapping=false;
    
    /*=============== FIT PARAMETERS ===============================*/
    
    /* Specify which parameters are your fit parameters. Parameters are fitted
     * to be the same in all datasets by default
     */
    
    
    pars_add(pars, "d_oxide", &(fit[0].m.d[1]), 5., 40);
    pars_add(pars, "d_Cr", &(fit[0].m.d[2]), 10, 60.);
    pars_add(pars, "d_gold", &(fit[0].m.d[3]), 80, 150.);
    
    pars_add(pars, "l_tether", &(l_tether), 6, 18);
    pars_add(pars, "l_lipid1", &(l_lipid1), 8, 22);
    pars_add(pars, "l_lipid2", &(l_lipid2), 7, 20);
    pars_add(pars, "dl_lipid_gcrinse1", &(dl_lipid_gcrinse1), -3, 3);
    //pars_add(pars, "protexchratio",      &(protexchratio),       0.5,  1.0);
    
    pars_add(pars, "nf_tether", &(nf_tether), 0.5, 1.0);
    pars_add(pars, "mult_tether", &(mult_tether), 0.1, 4.);
    //pars_add(pars, "nf_lipid_2", &(nf_lipid_2), 0.4, 0.6);
    
    pars_add(pars, "vf_bilayer", &(vf_bilayer), 0.90, 1);
    pars_add(pars, "vf_bilayer_gcrinse1", &(vf_bilayer_gcrinse1), 0.6, 1);
    pars_add(pars, "fraction_inc2", &(fraction_inc2), 0.1, 2.0);
    pars_add(pars, "fraction_inc3", &(fraction_inc3), 0.1, 2.0);
    pars_add(pars, "dp_on0", &(dp_on[0]), -7, 7);
    pars_add(pars, "dp_on1", &(dp_on[1]), -7, 7);
    pars_add(pars, "dp_on2", &(dp_on[2]), -7, 7);
    pars_add(pars, "dp_on3", &(dp_on[3]), -7, 7);
    pars_add(pars, "dp_on4", &(dp_on[4]), -7, 7);
    pars_add(pars, "dp_on5", &(dp_on[5]), -7, 7);
    pars_add(pars, "dp_on6", &(dp_on[6]), -7, 7);
    pars_add(pars, "dp_on7", &(dp_on[7]), -7, 7);
    pars_add(pars, "vf_on1", &(vf_on[1]), -0.04, 0.2);
    pars_add(pars, "vf_on2", &(vf_on[2]), -0.04, 0.4);
    pars_add(pars, "vf_on3", &(vf_on[3]), -0.04, 0.4);
    pars_add(pars, "vf_on4", &(vf_on[4]), -0.04, 0.4);
    pars_add(pars, "vf_on5", &(vf_on[5]), -0.04, 0.4);
    pars_add(pars, "vf_on6", &(vf_on[6]), -0.04, 0.4);
    pars_add(pars, "frac2_on0", &(frac2_on[0]), 0., 1.);
    pars_add(pars, "frac2_on1", &(frac2_on[1]), 0., 1.);
    pars_add(pars, "frac2_on2", &(frac2_on[2]), 0., 1.);
    pars_add(pars, "frac2_on3", &(frac2_on[3]), 0., 1.);
    pars_add(pars, "frac2_on4", &(frac2_on[4]), 0., 1.);
    pars_add(pars, "frac2_on5", &(frac2_on[5]), 0., 1.);
    pars_add(pars, "frac2_on6", &(frac2_on[6]), 0., 1.);
    pars_add(pars, "frac2_on7", &(frac2_on[7]), 0., 1.);
    pars_add(pars, "frac2_on8", &(frac2_on[8]), 0., 1.);
    
    pars_add(pars, "rho_SiOx", &(fit[0].m.rho[1]), 3.4e-6, 3.8e-6);
    pars_add(pars, "rho_Cr", &(fit[0].m.rho[2]), 3.03e-6, 4.15e-6);
    pars_add(pars, "rho_Au", &(fit[0].m.rho[3]), 4.2e-6, 4.6e-6);
    
    pars_add(pars, "rho_solv_0", &(fit[0].m.rho[fit[0].m.n-1]), 5.8e-6, 6.40e-6);
    pars_add(pars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), -0.56e-6, -0.1e-6);
    pars_add(pars, "rho_solv_2", &(fit[2].m.rho[fit[2].m.n-1]), 5.8e-6, 6.40e-6);
    pars_add(pars, "rho_solv_3", &(fit[3].m.rho[fit[3].m.n-1]), 3.5e-6, 4.5e-6);
    pars_add(pars, "rho_solv_4", &(fit[4].m.rho[fit[4].m.n-1]), -0.56e-6, -0.1e-6);
    
    pars_add(pars, "global_rough", &(global_rough), 1.0, 25.0);
    pars_add(pars, "rough_cr_au",  &(rough_cr_au),  1.0, 10.0);
    pars_add(pars, "sigma",        &(sigma),        2.0, 5.);
    pars_add(pars, "thetaoffset0", &(thetaoffset), -0.01, 0.01);
    
    pars_add(pars, "background_0", &(background[0]), -9, -6);
    pars_add(pars, "background_1", &(background[1]), -7, -5);
    pars_add(pars, "background_2", &(background[2]), -9, -6);
    pars_add(pars, "background_3", &(background[3]), -8, -5);
    pars_add(pars, "background_4", &(background[4]), -7, -5);
    
    /* Build a list of 'free parameters' in fit[1].pars. These are
     * parameters for which the values are aloowed to differ from those
     * in model 0.  By default all values in all models are the same unless 
     * specified here. The range data is not used here, so set it to [0,1].  
     */
    
    pars_add(freepars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), 0, 1);
    pars_add(freepars, "rho_solv_2", &(fit[2].m.rho[fit[2].m.n-1]), 0, 1);
    pars_add(freepars, "rho_solv_3", &(fit[3].m.rho[fit[3].m.n-1]), 0, 1);
    pars_add(freepars, "rho_solv_4", &(fit[4].m.rho[fit[4].m.n-1]), 0, 1);
    
    constraints = constr_models;
    output_model = save;
    return fit;
}
