import components as cmp
import molgroups as mol

DOPC = cmp.Lipid(name='DOPC', headgroup=mol.PC, tails=2 * [cmp.oleoyl], methyls=cmp.methyl)
POPC = cmp.Lipid(name='POPC', headgroup=mol.PC, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=cmp.methyl)
DOPS = cmp.Lipid(name='DOPS', headgroup=cmp.ps, tails=2 * [cmp.oleoyl], methyls=cmp.methyl)
POPE = cmp.Lipid(name='POPE', headgroup=cmp.pe, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=cmp.methyl)
POPG = cmp.Lipid(name='POPG', headgroup=cmp.pg, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=cmp.methyl)
chol = cmp.Lipid(name='chol', headgroup=None, tails=[cmp.cholesterol], methyls=None)
cardiolipin = cmp.Lipid(name='cardiolipin', headgroup=cmp.cardiolipin, tails=4 * [cmp.oleoyl],methyls=cmp.methyl)