import components as cmp

DOPC = cmp.Lipid(name='DOPC', headgroup=cmp.PC, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
POPC = cmp.Lipid(name='POPC', headgroup=cmp.PC, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
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
