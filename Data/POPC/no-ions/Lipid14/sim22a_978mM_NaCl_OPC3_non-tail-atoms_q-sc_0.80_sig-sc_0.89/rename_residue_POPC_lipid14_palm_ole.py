#
# simple python script for conversion of 
# Lipid14 POPC residue names with duplicate atom names 
# to a split-residue convention: 
#   POPC will be split into 3 virtual residues with names
#   PALM - POPC - OLE
#   so that duplicate atom names can be distinquished at least
#   by a different residue name
#
# script assumes the particular ordering as above with 
# the head group segment beginning with O12 and ending with O22
# the rest (begin/end) being the tails PALM/OLE.

import MDAnalysis as mda

# conf_init.gro should always be present in the sim dir
u = mda.Universe("conf_init.gro")
# find some POPC residue
for r in u.residues:
    if r.name == "POPC":
        some_popc_res = r
        break
# find the indices of the first/last head group atoms
headgr_begin_index = some_popc_res.names.tolist().index("O12")
headgr_end_index   = some_popc_res.names.tolist().index("O22")

# resname the resnames of PALM and OLE segments
for r in u.residues:
    if r.name == "POPC":
        for a in r.atoms[:headgr_begin_index]:
            a.resname = "PALM"
        for a in r.atoms[headgr_end_index+1:]:
            a.resname = "OLE"

# write the new topology out
u.SYSTEM.write("conf_palm.gro")
