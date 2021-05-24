// version 21/03/11


class nSLDObj
{
    public:
	    nSLDObj();
		virtual ~nSLDObj();
	    virtual double fnGetLowerLimit();
		virtual double fnGetUpperLimit();
		virtual double fnGetArea(double z);
		virtual double fnGetnSLD(double z);
	    virtual double fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
	    virtual void fnOverlayProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
        virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
        virtual void fnWriteData2File (FILE *fp, int dimension, double stepsize);
		
		bool bWrapping;
};

//------------------------------------------------------------------------------------------------------
class BoxErr : public nSLDObj
{
    public:
	    BoxErr() {};
	    BoxErr(double z, double sigma, double length, double vol, double nSL, double numberfraction);
		virtual ~BoxErr();
	    virtual double fnGetLowerLimit();
		virtual double fnGetUpperLimit();
		virtual double fnGetArea(double z);
		virtual double fnGetnSLD(double z);
        virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
		
		double z, sigma, vol, nSL, nf, l;               //variable used internally
};
//------------------------------------------------------------------------------------------------------
class Box2Err : public nSLDObj
{
    public:
	    Box2Err() {};
	    Box2Err(double z, double sigma1, double sigma2, double length, double vol, double nSL, double numberfraction);
		virtual ~Box2Err();
	    virtual double fnGetLowerLimit();
		virtual double fnGetUpperLimit();
		virtual double fnGetArea(double z);
		virtual double fnGetnSLD(double z);
		virtual void fnSetAllSigma(double sigma);
		virtual void fnSetZ(double dz);
        virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
		
		double z, sigma1, sigma2, vol, nSL, nf, l;               //variable used internally
};

//------------------------------------------------------------------------------------------------------

class Gaussian : public nSLDObj
{
    public:
    
	    Gaussian(double C, double sigma, double vol, double nSL, double numberfraction);
		virtual ~Gaussian();
	    virtual double fnGetLowerLimit();
		virtual double fnGetUpperLimit();
		virtual double fnGetArea(double z);
		virtual double fnGetnSLD(double z);
        virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
		
		double z, sigma, vol, nSL, nf;               //variable used internally
};

//------------------------------------------------------------------------------------------------------

class Parabolic: public nSLDObj
{
    public:
	    Parabolic(double dC, double dH, double dn, double dnSLD, double dnumberfraction);
		virtual ~Parabolic();
	    virtual double fnGetLowerLimit();
		virtual double fnGetUpperLimit();
		virtual double fnGetArea(double z);
		virtual double fnGetnSLD(double z);
        virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
		
		double C, H, n, nSLD, nf;               //variable used internally
};
//------------------------------------------------------------------------------------------------------
class StretchGaussian : public nSLDObj
{
    public:
	    StretchGaussian(double z, double sigma, double length, double vol, double nSL, double numberfraction);
		virtual ~StretchGaussian();
	    virtual double fnGetLowerLimit();
		virtual double fnGetUpperLimit();
		virtual double fnGetArea(double z);
		virtual double fnGetnSLD(double z);
        virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
		
		double z, sigma, vol, nSL, nf, l;               //variable used internally
};


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//Composite groups
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------
//PC headgroup
//------------------------------------------------------------------------------------------------------
class PC: public nSLDObj
{
        protected:
                double z;
        public:
                PC();
                virtual ~PC();
                virtual void   fnAdjustParameters();
                virtual double fnGetLowerLimit();
                virtual double fnGetUpperLimit();
                virtual double fnGetArea(double z);
                virtual double fnGetnSLD(double z);
				virtual double fnGetZ() {return z;};
                virtual void fnSetAllSigma(double sigma);
                virtual void fnSetZ(double dz);
                virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
                
                Box2Err *cg;
                Box2Err *phosphate;
                Box2Err *choline;
                double l, nf;
};
//------------------------------------------------------------------------------------------------------
//mirrored PC headgroup
//------------------------------------------------------------------------------------------------------

class PCm: public PC
{
        public:
			PCm();
			virtual ~PCm();
			virtual void fnAdjustParameters();
			virtual double fnGetLowerLimit();
			virtual double fnGetUpperLimit();
			virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
};



//---------------------------------------------------------------------------------------------------------
//Freeform group 10 boxes
//---------------------------------------------------------------------------------------------------------
class FreeBox: public nSLDObj
{

	public:
	
		FreeBox(int n, double dstartposition, double dnSLD, double dnormarea);
		virtual ~FreeBox();
		virtual void fnAdjustParameters();
		virtual double fnGetArea(double dz);
		virtual double fnGetnSLD(double dz);
		virtual double fnGetLowerLimit();
		virtual double fnGetUpperLimit();
		virtual void fnSetStartposition(double dz);
		virtual void fnSetNormarea(double dnormarea);
		virtual void fnSetnSLD(double dnSLD);
		virtual void fnSetAllSigma(double sigma);
		virtual void fnWritePar2File(FILE *fp, char *cName, int dimension, double stepsize);
		
		Box2Err *box1, *box2, *box3, *box4, *box5, *box6, *box7, *box8, *box9, *box10;
		
		int numberofboxes;
		double vf1, vf2, vf3, vf4, vf5, vf6, vf7, vf8, vf9, vf10;
		double normarea, startposition, nSLD;

	
};
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//AminoAcid
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//general implementation
//-----------------------------------------------------------------------------------------------------------
class AminoAcid: public BoxErr
{
	public:
		AminoAcid();
		virtual ~AminoAcid() {};
		virtual double fnGetnSLD(double z);
		
		int nH, nExch, Deuterated;
		double ExchangeRatio;
};
//-----------------------------------------------------------------------------------------------------------
//Specific implementations
//-----------------------------------------------------------------------------------------------------------
class AA_Lys: public AminoAcid
{
	public:
		AA_Lys();
		virtual ~AA_Lys() {};
};
class AA_Arg: public AminoAcid
{
	public:
		AA_Arg();
		virtual ~AA_Arg() {};
};
class AA_His: public AminoAcid
{
	public:
		AA_His();
		virtual ~AA_His() {};
};
class AA_Asn: public AminoAcid
{
	public:
		AA_Asn();
		virtual ~AA_Asn() {};
};
class AA_Asp: public AminoAcid
{
	public:
		AA_Asp();
		virtual ~AA_Asp() {};
};
class AA_Cys: public AminoAcid
{
	public:
		AA_Cys();
		virtual ~AA_Cys() {};
};
class AA_Thr: public AminoAcid
{
	public:
		AA_Thr();
		virtual ~AA_Thr() {};
};
class AA_Ser: public AminoAcid
{
	public:
		AA_Ser();
		virtual ~AA_Ser() {};
};
class AA_Gln: public AminoAcid
{
	public:
		AA_Gln();
		virtual ~AA_Gln() {};
};
class AA_Glu: public AminoAcid
{
	public:
		AA_Glu();
		virtual ~AA_Glu() {};
};
class AA_Pro: public AminoAcid
{
	public:
		AA_Pro();
		virtual ~AA_Pro() {};
};
class AA_Gly: public AminoAcid
{
	public:
		AA_Gly();
		virtual ~AA_Gly() {};
};
class AA_Ala: public AminoAcid
{
	public:
		AA_Ala();
		virtual ~AA_Ala() {};
};
class AA_Val: public AminoAcid
{
	public:
		AA_Val();
		virtual ~AA_Val() {};
};
class AA_Ile: public AminoAcid
{
	public:
		AA_Ile();
		virtual ~AA_Ile() {};
};
class AA_Leu: public AminoAcid
{
	public:
		AA_Leu();
		virtual ~AA_Leu() {};
};
class AA_Met: public AminoAcid
{
	public:
		AA_Met();
		virtual ~AA_Met() {};
};
class AA_Tyr: public AminoAcid
{
	public:
		AA_Tyr();
		virtual ~AA_Tyr() {};
};
class AA_Phe: public AminoAcid
{
	public:
		AA_Phe();
		virtual ~AA_Phe() {};
};
class AA_Trp: public AminoAcid
{
	public:
		AA_Trp();
		virtual ~AA_Trp() {};
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//Lipid bilayers
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//General Implementations
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//solid supported lipid bilayer, single-lipid bilayer
//------------------------------------------------------------------------------------------------------
class ssBLM: public nSLDObj
{
	protected:
		double normarea;
	public:
                ssBLM();
                virtual ~ssBLM();
                virtual void   fnAdjustParameters();
                virtual double fnGetLowerLimit();
                virtual double fnGetUpperLimit();
                virtual double fnGetArea(double z);
                virtual double fnGetnSLD(double z);
                virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
				virtual double fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
                
				Box2Err   *substrate;
				PCm       *headgroup1;                                                 //mirrored PC head group
				Box2Err   *lipid1;
				Box2Err   *methyl1;
				Box2Err   *methyl2;
				Box2Err   *lipid2;
				PC        *headgroup2;                                                //PC head group
				
				//primary fit parameters
                double global_rough, sigma, l_lipid1, l_lipid2, vf_bilayer, rho_substrate, l_submembrane;
				
				//other parameters
				double volacyllipid, nslacyllipid, volmethyllipid, nslmethyllipid;
};
	
//------------------------------------------------------------------------------------------------------
//tethered lipid bilayer, single-lipid bilayer
//------------------------------------------------------------------------------------------------------
class tBLM: public nSLDObj
{
        protected:
		        double normarea;
        public:
                tBLM();
                virtual ~tBLM();
                virtual void   fnAdjustParameters();
                virtual double fnGetLowerLimit();
                virtual double fnGetUpperLimit();
                virtual double fnGetArea(double z);
                virtual double fnGetnSLD(double z);
                virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
				virtual double fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
                
				Box2Err   *substrate;
				Box2Err   *bME;
				Box2Err   *tether;
				Box2Err   *tetherg; 
				PCm       *headgroup1;                                                 //mirrored PC head group
				Box2Err   *lipid1;
				Box2Err   *methyl1;
				Box2Err   *methyl2;
				Box2Err   *lipid2;
				PC        *headgroup2;                                                //PC head group
				
				//primary fit parameters
                double global_rough, sigma, l_lipid1, l_lipid2, vf_bilayer, l_tether, nf_tether, mult_tether, rho_substrate;
				
				//other parameters
				double volacyllipid, nslacyllipid, volmethyllipid, nslmethyllipid, volmethyltether;
				double nslmethyltether, volacyltether, nslacyltether;

};
//------------------------------------------------------------------------------------------------------
//tethered lipid bilayer, binary lipid bilayer
//------------------------------------------------------------------------------------------------------
class tBLM_binary: public tBLM
{
        protected:
		        double normarea;
        public:
                tBLM_binary();
                virtual ~tBLM_binary();
                virtual void   fnAdjustParameters();
                virtual double fnGetLowerLimit();
                virtual double fnGetUpperLimit();
                virtual double fnGetArea(double z);
                virtual double fnGetnSLD(double z);
                virtual void fnWritePar2File (FILE *fp, char *cName, int dimension, double stepsize);
				virtual double fnWriteProfile(double aArea[], double anSLD[], int dimension, double stepsize, double dMaxArea);
                
				Box2Err	  *headgroup1_2;	  
				Box2Err   *headgroup2_2;
				
				//primary fit parameters
                double nf_lipid_2;
				
				//other parameters
				double volacyllipid_2, nslacyllipid_2, volmethyllipid_2, nslmethyllipid_2;
};

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
// Bilayers, specific implementation
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
class ssBLM_POPC: public ssBLM
{
        public:
                ssBLM_POPC();
                virtual ~ssBLM_POPC() {};
};

class tBLM_WC14_DPhyPC: public tBLM
{
        public:
                tBLM_WC14_DPhyPC();
                virtual ~tBLM_WC14_DPhyPC() {};
};
class tBLM_WC14_DOPC: public tBLM
{
        public:
                tBLM_WC14_DOPC();
                virtual ~tBLM_WC14_DOPC() {};
};
class tBLM_HC18_DOPC: public tBLM
{
        public:
                tBLM_HC18_DOPC();
                virtual ~tBLM_HC18_DOPC() {};
};
class tBLM_HC18_DOPC_DOPS: public tBLM_binary
{
        public:
                tBLM_HC18_DOPC_DOPS();
                virtual ~tBLM_HC18_DOPC_DOPS() {};
};

class tBLM_HC18_POPC_POPA: public tBLM_binary
{
        public:
                tBLM_HC18_POPC_POPA();
                virtual ~tBLM_HC18_POPC_POPA() {};
};

class tBLM_HC18_POPC_POPG: public tBLM_binary
{
        public:
                tBLM_HC18_POPC_POPG();
                virtual ~tBLM_HC18_POPC_POPG() {};
};
//-----------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
// Standalone
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void fnWriteConstant(FILE *fp, char *cName, double area, double nSLD, int dimension, double stepsize);
void fnOverlayCanvasOnCanvas(double aArea[], double anSL[], double aArea2[], double anSL2[], int dimension, double dMaxArea);
double fnClearCanvas(double aArea[], double anSL[], int dimension);
void fnWriteCanvas2Model(double aArea[], double anSL[], fitinfo fit[], int gaussstart, int dimension, double stepsize, double dMaxArea, double normarea, int modelstart, int modelend);
