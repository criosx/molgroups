from . import components as cmp
from . import mol

# Composite Headgroups
PC = [mol.CompositeHeadgroup, {'components': [cmp.carbonyl_glycerol, cmp.phosphate, cmp.choline],
                               'sigma1': [2.53, 2.29, 2.02],
                               'sigma2': [2.29, 2.02, 2.26],
                               'rel_pos': [0., 0.58, 1.0],
                               'length': 9.575}]

# Lipids
DOPC = cmp.Lipid(name='DOPC', headgroup=PC, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
POPC = cmp.Lipid(name='POPC', headgroup=PC, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
DOPS = cmp.Lipid(name='DOPS', headgroup=cmp.ps, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
POPE = cmp.Lipid(name='POPE', headgroup=cmp.pe, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
POPG = cmp.Lipid(name='POPG', headgroup=cmp.pg, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
chol = cmp.Lipid(name='chol', headgroup=None, tails=[cmp.cholesterol], methyls=None)
cardiolipin = cmp.Lipid(name='cardiolipin', headgroup=cmp.cardiolipin, tails=4 * [cmp.oleoyl],methyls=[cmp.methyl])

# Combined tethers
WC14 = cmp.Tether(name='WC14', tether=cmp.SEO6, tetherg=cmp.tetherg_ether, tails=[cmp.myristoyl, cmp.myristoyl],
                  methyls=[cmp.methyl])
HC18 = cmp.Tether(name='HC18', tether=cmp.SEO6, tetherg=cmp.tetherg_ether, tails=[cmp.oleoyl, cmp.oleoyl],
                  methyls=[cmp.methyl])
HC18SAc = cmp.Tether(name='HC18SAc', tether=cmp.SAcEO6, tetherg=cmp.tetherg_ether, tails=[cmp.oleoyl, cmp.oleoyl],
                     methyls=[cmp.methyl])

