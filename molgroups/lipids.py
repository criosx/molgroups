from molgroups import components as cmp
from molgroups import mol

# Composite Headgroups
PC = [mol.CompositeHeadgroup, {'components': [cmp.carbonyl_glycerol, cmp.phosphate, cmp.choline],
                               'sigma1': [2.53, 2.29, 2.02],
                               'sigma2': [2.29, 2.02, 2.26],
                               'rel_pos': [0., 0.58, 1.0],
                               'length': 9.575}]

# Lipids
DLPE = cmp.Lipid(name='DLPE', headgroup=cmp.pe, tails=2 * [cmp.linoleoyl], methyls=[cmp.methyl])
DLPG = cmp.Lipid(name='DLPG', headgroup=cmp.pg, tails=2 * [cmp.lauroyl], methyls=[cmp.methyl])
DMPC = cmp.Lipid(name='DMPC', headgroup=PC, tails=2 * [cmp.myristoyl], methyls=[cmp.dmethyl])
d54DMPC = cmp.Lipid(name='d54DMPC', headgroup=PC, tails=2 * [cmp.d_myristoyl], methyls=[cmp.dmethyl])
DOPC = cmp.Lipid(name='DOPC', headgroup=PC, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
DOPS = cmp.Lipid(name='DOPS', headgroup=cmp.ps, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
DOPE = cmp.Lipid(name='DOPE', headgroup=cmp.pe, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
DOTAP = cmp.Lipid(name='DOTAP', headgroup=cmp.tap, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
DPSM = cmp.Lipid(name='DPSM', headgroup=cmp.sm, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
# approximation
LPS_short = cmp.Lipid(name='LIPID_A', headgroup=cmp.lps_short, tails=2*[cmp.tridecanoyl]+4*[cmp.lauroyl],
                      methyls=[cmp.methyl])
PAPC = cmp.Lipid(name='PAPC', headgroup=PC, tails=[cmp.palmitoyl, cmp.arachidonoyl], methyls=[cmp.methyl])
PAPS = cmp.Lipid(name='PAPS', headgroup=cmp.ps, tails=[cmp.palmitoyl, cmp.arachidonoyl], methyls=[cmp.methyl])
PLPC = cmp.Lipid(name='PLPC', headgroup=PC, tails=[cmp.linoleoyl, cmp.oleoyl], methyls=[cmp.methyl])
POPC = cmp.Lipid(name='POPC', headgroup=PC, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
POPE = cmp.Lipid(name='POPE', headgroup=cmp.pe, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
POPG = cmp.Lipid(name='POPG', headgroup=cmp.pg, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
POPS = cmp.Lipid(name='POPS', headgroup=cmp.ps, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
POSM = cmp.Lipid(name='POSMP', headgroup=cmp.sm, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
POPIP2 = cmp.Lipid(name='POPIP2', headgroup=cmp.pip2, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
SOPC = cmp.Lipid(name='SOPC', headgroup=PC, tails=[cmp.steraoyl, cmp.oleoyl], methyls=[cmp.methyl])
chol = cmp.Lipid(name='chol', headgroup=None, tails=[cmp.cholesterol], methyls=None)
cardiolipin = cmp.Lipid(name='cardiolipin', headgroup=cmp.cardiolipin, tails=4 * [cmp.oleoyl],methyls=[cmp.methyl])

# Combined tethers
WC14 = cmp.Tether(name='WC14', tether=cmp.SEO6, tetherg=cmp.tetherg_ether, tails=[cmp.myristoyl, cmp.myristoyl],
                  methyls=[cmp.methyl])
HC18 = cmp.Tether(name='HC18', tether=cmp.SEO6, tetherg=cmp.tetherg_ether, tails=[cmp.oleoyl, cmp.oleoyl],
                  methyls=[cmp.methyl])
HC18SAc = cmp.Tether(name='HC18SAc', tether=cmp.SAcEO6, tetherg=cmp.tetherg_ether, tails=[cmp.oleoyl, cmp.oleoyl],
                     methyls=[cmp.methyl])

