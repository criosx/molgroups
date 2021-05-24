#ifndef MOLGROUPS_H 
#define MOLGROUPS_H
//---------------abstract base class---------------------------------------------------------------------
class nSLDObj
{
public:
    nSLDObj();
    virtual ~nSLDObj();
    virtual double fnGetLowerLimit() = 0;
    virtual double fnGetUpperLimit() = 0;
    virtual double fnGetArea(double z) = 0;
    virtual double fnGetConvolutedArea(double z);
    virtual double fnGetnSLD(double z) = 0;
    virtual void   fnSetConvolution(double sigma, int n);
    virtual void   fnSetSigma(double sigma) = 0;
    virtual void   fnSetZ(double dz) {z=dz;};
    virtual void   fnSetnSL(double d) {nSL=d;};
    virtual void   fnSetnSL(double d1, double d2) {};
    virtual void   fnSetnSL(double d1, double d2, double d3) {};
    virtual void   fnSetnSL(double d1, double d2, double d3, double d4) {};
    virtual void   fnSetnSL(double d1, double d2, double d3, double d4, double d5) {};
    virtual double fnGetZ(){return z;};
    virtual double fnGetAbsorb(double z);
    virtual double fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
    virtual double fnWriteProfile(double aArea[], double anSLD[], double aAbsorb[], int dimension, double stepsize, double dMaxArea);
    virtual void   fnOverlayProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
    virtual void   fnOverlayProfile(double aArea[], double anSLD[], double aAbsorb[], int dimension, double stepsize, double dMaxArea);
    virtual void   fnWritePar2File (FILE *fp, const char *cName, int dimension, double stepsize) = 0;
    virtual void   fnWriteData2File (FILE *fp, const char *cName, int dimension, double stepsize);
    
    int iNumberOfConvPoints;
    bool bWrapping, bConvolution, bProtonExchange;
    double absorb, z, l, nf, nSL, nSL2, vol, dSigmaConvolution;
    
protected:
    virtual double CatmullInterpolate(double t, double pm1, double p0, double p1, double p2);
    virtual double fnTriCubicCatmullInterpolate(double p[4][4][4],double t[3]);
    virtual double fnQuadCubicCatmullInterpolate(double p[4][4][4][4],double t[4]);

};

class Box2Err : public nSLDObj
{
public:
    Box2Err() {};
    Box2Err(double z, double sigma1, double sigma2, double length, double vol, double nSL, double numberfraction);
    virtual ~Box2Err();
    virtual double fnGetLowerLimit();
    virtual double fnGetUpperLimit();
    virtual double fnGetArea(double z);
    virtual double fnGetnSL(double bulknsld);
    virtual double fnGetnSLD(double z);
    virtual double fnGetnSLD(double z, double bulknsld);
    virtual void   fnSetnSL(double d1, double d2);
    virtual void   fnSetSigma(double sigma);
    virtual void   fnSetSigma(double sigma1, double sigma2);
    virtual void   fnSetZ(double dz);
    virtual void   fnWritePar2File (FILE *fp, const char *cName, int dimension, double stepsize);
    
    double sigma1, sigma2, nsldbulk_store;
};

class BLM_quaternary: public nSLDObj
{
protected:
    double normarea;
public:
    BLM_quaternary();
    virtual ~BLM_quaternary();
    virtual void   fnAdjustParameters();
    virtual double fnGetLowerLimit();
    virtual double fnGetUpperLimit();
    virtual double fnGetArea(double z);
    virtual double fnGetCenter();
    virtual double fnGetnSLD(double z);
    virtual void fnInit(double va1, double na1, double vm1, double nm1, double vh1, double nh1, double lh1, double va2, double na2, double vm2, double nm2, double vh2, double nh2, double lh2, double va3, double na3, double vm3, double nm3, double vh3, double nh3, double lh3, double vc, double nc);
    virtual void   fnSet(double sigma, double bulknsld, double startz, double l_lipid1, double l_lipid2, double vf_bilayer, double nf_lipid_2=0, double nf_lipid3=0, double nf_chol=0, double hc_substitution_1=0, double hc_substitution_2=0, double radius_defect=100);
    virtual void fnSetSigma(double sigma);
    virtual void fnWritePar2File (FILE *fp, const char *cName, int dimension, double stepsize);
    virtual double fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
    
    Box2Err   *headgroup1;                                                 //mirrored PC head group
    Box2Err   *lipid1;
    Box2Err   *methyl1;
    Box2Err   *methyl2;
    Box2Err   *lipid2;
    Box2Err   *headgroup2;                                                //PC head group
    Box2Err	  *headgroup1_2;
    Box2Err   *headgroup2_2;
    Box2Err	  *headgroup1_3;
    Box2Err   *headgroup2_3;
    
    Box2Err   *defect_hydrocarbon;
    Box2Err   *defect_headgroup;
    
    //primary fit parameters
    double sigma, l_lipid1, l_lipid2, vf_bilayer, startz;
    double hc_substitution_1, hc_substitution_2, radius_defect, bulknsld;
    double nf_lipid_2, nf_lipid_3, nf_chol;
    
    //other parameters
    double volacyllipid, nslacyllipid, volmethyllipid, nslmethyllipid;
    double volacyllipid_2, nslacyllipid_2, volmethyllipid_2, nslmethyllipid_2;
    double volacyllipid_3, nslacyllipid_3, volmethyllipid_3, nslmethyllipid_3;
    double volchol, nslchol;
};


void fnWriteConstant(FILE *fp, const char *cName, double area, double nSLD, int dimension, double stepsize);

#endif /* MOLGROUPS_H */
