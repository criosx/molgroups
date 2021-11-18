import components as cmp
from molgroups import PC

DOPC = cmp.Lipid(name='DOPC', headgroup=PC, tails=2 * [cmp.oleoyl], methyls=cmp.methyl)
POPC = cmp.Lipid(name='POPC', headgroup=PC, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=cmp.methyl)
DOPS = cmp.Lipid(name='DOPS', headgroup=cmp.ps, tails=2 * [cmp.oleoyl], methyls=cmp.methyl)
chol = cmp.Lipid(name='chol', headgroup=None, tails=[cmp.cholesterol], methyls=None)