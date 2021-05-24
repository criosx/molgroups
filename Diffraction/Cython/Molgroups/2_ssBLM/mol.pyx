# distutils: language = c++

from molgroups cimport Box2Err, ssBLM_POPC_POPS_PIP_chol

cdef class PyBox2Err:
    cdef Box2Err c_box2err

    def __cinit__(self, double d0, double d1, double d2, double d3, double d4, double d5):
        self.c_box2err = Box2Err(d0, d1, d2, d3, d4, d5)

    def get_lowerlimit(self):
        return self.c_box2err.fnGetLowerLimit()

cdef class Py_ssBLM_POPC_POPS_PIP_chol:
    cdef ssBLM_POPC_POPS_PIP_chol c_ssBLM_POPC_POPS_PIP_chol

   #  def __cinit__(self):
   #     self.c_ssBLM_POPC_POPS_PIP_chol = ssBLM_POPC_POPS_PIP_chol()

    def AdjustParameters(self):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnAdjustParameters()

    def GetLowerLimit(self):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetLowerLimit()

    def GetUpperLimit(self):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetUpperLimit()

    def GetArea(self, double z):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetArea(z)

    def GetnSLD(self, double z):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetnSLD(z)

    def Set(self, double sigma, double global_rough, double rho_substrate, double bulknsld, double rho_siox, double l_siox,
            double l_submembrane, double l_lipid1, double l_lipid2, double vf_bilayer, double nf_lipid_2=0,
            double nf_lipid3=0, double nf_chol=0, double hc_substitution_1=0, double hc_substitution_2=0, double radius_defect=100):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnSet(sigma, global_rough, rho_substrate,
        bulknsld, rho_siox, l_siox, l_submembrane, l_lipid1, l_lipid2, vf_bilayer,
        nf_lipid_2, nf_lipid3, nf_chol, hc_substitution_1, hc_substitution_2, radius_defect)

    def SetSigma(self, double sigma):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnSetSigma(sigma)
