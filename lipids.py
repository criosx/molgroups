from periodictable.fasta import Molecule
from molgroups import CompositenSLDObj
# in formulas, use T for exchangeable hydrogens, then Molecule.D2Osld(D2O_fraction=1) to get D2O SLD.
pc = Molecule(name='pc', formula='C10 H18 O8 N P', cell_volume=331.00)
oleoyl = Molecule(name='oleoyl', formula='C17 H33', cell_volume=474.90/2.0)

class Lipid(object):
    def __init__(self, name=None, hg=pc, tails=[oleoyl, oleoyl]):
        # hg = Molecule object with headgroup information or molgroups object
        # tails = List of molecule objects containing lipid tail information
        # default is DOPC

        # TODO: should Box2Errs be created here, or in the __init__ of the BLM object?

        # hg section
        if isinstance(hg, Molecule):
            self.hg = hg
        elif isinstance(hg, CompositenSLDObj):
            # do other stuff?
            self.hg = hg
        
        # tails section
        total_volume = sum([t.cell_volume for t in tails])
        if isinstance(tails, list):
            total_formula = ' '.join([str(t.formula) for t in tails])
            self.tails = Molecule(name='tails',formula=total_formula, cell_volume=total_volume)
        elif isinstance(tails, Molecule):
            self.tails = tails
        else:
            TypeError('Lipid.tails must be Molecule or list of Molecules')
        
        self.hg = hg

dopc = Lipid()
nSL, nSL2 = dopc.hg.Hsld * dopc.hg.cell_volume, dopc.hg.Dsld * dopc.hg.cell_volume
print(nSL, nSL2)

