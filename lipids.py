from periodictable.fasta import Molecule
from molgroups import CompositenSLDObj
# in formulas, use H[1] for exchangeable hydrogens, then Molecule.sld, Molecule.Dsld
# to get limiting slds, or Molecule.D2Osld(D2O_fraction=??) to get arbitrary D2O fraction.

class Headgroup(Molecule):
    # Subclasses Molecule to automatically store a headgroup length for use later
    # and calculate total scattering lengths
    def __init__(self, length=9.575, **kwargs):
        super().__init__(**kwargs)
        self.length = length
        self.nSL = self.sld * self.cell_volume * 1e-6 # units of Ang
        self.DnSL = self.Dsld * self.cell_volume * 1e-6 # units of Ang

# headgroup pieces
choline = Molecule(name='choline', formula='C5 H13 N', cell_volume=120.)
phosphate = Molecule(name='phosphate', formula='PO4', cell_volume=54.)
carbonyl_glycerol = Molecule(name='carbonyl + glycerol', formula='C5 O4 H5', cell_volume=147.)

# headgroups
pc = Headgroup(name='pc', formula='C10 H18 O8 N P', cell_volume=331.00, length=9.575)
pe = Headgroup(name='pe', formula='C7 H9 H[1]3 O8 N P', cell_volume=262., length=7.7)
pg = Headgroup(name='pg', formula='C8 H10 H[1]2 O10 P', cell_volume=240., length=7.8)
ps = Headgroup(name='ps', formula='C8 H8 H[1]3 N O10 P', cell_volume=280., length=8.1)
pa = Headgroup(name='pa', formula='C5 H5 H[1] O8 P', cell_volume=174., length=5.0)
pi = Headgroup(name='pi', formula='C11 H7 H[1]5 O13 P', cell_volume=370.7, length=10.7)
pip2 = Headgroup(name='pi(4,5)p2', formula='C11 H7 H[1]5 O19 P3', cell_volume=500., length=12.0) # doesn't match molgroups.cc
cardiolipin = Headgroup(name='cardiolipin', formula='C13 H15 H[1] O17 P2', cell_volume=684.4, length=9.56)

oleoyl = Molecule(name='oleoyl', formula='C17 H33', cell_volume=474.90/2.0)
palmitoyl = Molecule(name='palmitoyl', formula='C15 H31', cell_volume=770./2.0)
myristoyl = Molecule(name='myristoyl', formula='C13 H27', cell_volume=770./2.0)
phytanoyl = Molecule(name='phytanoyl', formula='C19 H39', cell_volume=1095./2.0)

# TODO: make a Tether class that can be passed to tethered bilayers. Perhaps make Tether class none for supported bilayers?

class Lipid(object):
    def __init__(self, name=None, hg=pc, tails=[oleoyl, oleoyl], **kwargs):
        # hg = Molecule object with headgroup information or PC molgroups object
        # tails = List of molecule objects containing lipid tail information
        # default is DOPC
       
        # tails section
        tail_volume = sum([t.cell_volume for t in tails])
        if isinstance(tails, list):
            tail_formula = ' '.join([str(t.formula) for t in tails])
            self.tails = Molecule(name='tails',formula=tail_formula, cell_volume=tail_volume)
        elif isinstance(tails, Molecule):
            self.tails = tails
        else:
            TypeError('Lipid.tails must be Molecule or list of Molecules')
        
        self.hg = hg
        assert(isinstance(hg, Molecule) | isinstance(hg, CompositenSLDObj))
        TypeError('Lipid.hg must be Headgroup or CompositenSLDObj')

if __name__ == '__main__':
    dopc = Lipid()
    nSL, nSL2 = dopc.hg.sld * dopc.hg.cell_volume, dopc.hg.Dsld * dopc.hg.cell_volume
    print(nSL, nSL2)

