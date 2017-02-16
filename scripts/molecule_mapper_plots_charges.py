# -*- coding: utf-8 -*-
"""
script for reordering atoms in molecules
according to the names of a reference molecule

Originally meant to reorder POPC atoms from CHARMM-GUI to Slipids ordering
"""

import pmx
import pmx.forcefield
import numpy as np
import matplotlib.pyplot as plt


class MoleculesCompared:
    """
    Class for comparing molecules (pmx.Model instances)
    (meant mostly for two molecules, but written for any number)

    Initiated with a list of molecules to be compared.
    Molecules have to be of the same size (no. atoms),
    and they must be the instance of pmx.Molecule.

    Residue name is also checked, but only warning is given if it differs.

    1st molecule in the list is the reference molecule

    especially: meant for sorting CHARMM-GUI POPC lipids
    to match the order of Slipids

    See the code for further details
    """
    def __init__(self, mols):
        self.natoms = None
        self.resname = None

        if not isinstance(mols, (tuple, list)):
            raise RuntimeError, "given object is neither a molecule nor a list"

        for item in mols:
            if not isinstance(item, pmx.molecule.Molecule):
                raise RuntimeError, "at least one item in the list is not the instance of pmx.Molecule"
            else:
                if self.resname == None:
                    # should be done only for the first item
                    self.resname = item.resname
                elif self.resname != item.resname:
                    print "WARNING: molecules don't have the same resname! ", self.resname, "  -vs-  ", item.resname
                if self.natoms == None:
                    # should be done only for the first item
                    self.natoms = len(item.atoms)
                elif self.natoms != len(item.atoms):
                    raise RuntimeError, "the number of atoms in the given molecules don't agree -- probably different types of molecules"

        self.mols = mols

    def get_eq_at_index(self, atom):
        """
        Returns the index of the equivalent atom in the ref-molecule
        (ref-mol: the first one in the list-of-molecules self.mols)
        """
        refmol=self.mols[0]
        # fetch atom instance(s) with the same name
        eq_at = refmol.fetch(atom.name)
        no_match_id = 0
        if len(eq_at) != 1:
            #raise RuntimeError, "Atom not found or multiple atoms with the same name exist!"
            no_match_id -= 1
            return no_match_id
        else:
            return eq_at[0].id

    def get_ref_sorted_atoms(self):
        """
        Returns self.mols (list of mols) with name-sorted atoms
        after the first molecule in the list
        """
        sorted_list_mols = [self.mols[0]]
        for mol in self.mols[1:]:
            # sort atoms in the current molecule mol
            sorted_list_atoms = sorted(mol.atoms, key=self.get_eq_at_index)
            # turn the list of atoms into a pmx.Molecule instance
            sorted_mol = pmx.molecule.Molecule(resname=mol.resname,id=mol.id)
            for atom in sorted_list_atoms:
                sorted_mol.append(atom)
            # append the sorted_mol (pmx.Molecule instance)
            # at the end of the list of molecules
            sorted_list_mols.append(sorted_mol)
        return sorted_list_mols

#%%

if __name__ == '__main__':

    # load-in the files
    mol_a = pmx.forcefield.ITPFile("l14_POPC.itp")
    mol_b = pmx.forcefield.ITPFile("slipids_POPC.itp")

#%%

    # translating A->B, so
    # B is a reference being a dictionary in order  FFname - Mapping_name
    # whereas A is being translated so a dictionary in opposite order is more practical
    # A: Mapping_name -> FFname
    mapping_xa = {}
    with open("mappingPOPClipid14.txt","r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                items = line.split()
                mapping_xa[items[0]] = items[1]

    mapping_bx = {}
    with open("mappingFILE.txt","r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                items = line.split()
                mapping_bx[items[1]] = items[0]


#%%

    mapping_ab = {}
    for key in mapping_bx.keys():
        mapping_ab[mapping_xa[mapping_bx[key]]] = key

    mapping_ba = {}
    for key in mapping_bx.keys():
        mapping_ba[key] = mapping_xa[mapping_bx[key]]

#%%
    sorted_atoms_list_names = []
    sorted_atoms_list_qa = []
    sorted_atoms_list_qb = []
    for atom in mol_b.atoms:
        other_atom_name = mapping_ba[atom.name]
        for otherat in mol_a.atoms:
            if otherat.name == other_atom_name:
                break
        print atom.name, atom.q, otherat.name, otherat.q
        sorted_atoms_list_names.append(atom.name)
        sorted_atoms_list_qb.append(atom.q)
        sorted_atoms_list_qa.append(otherat.q)


#%%

    plt.xticks(range(len(sorted_atoms_list_names)), sorted_atoms_list_names, rotation=30)
    plt.plot(sorted_atoms_list_qa, label="Lipid14")
    plt.plot(sorted_atoms_list_qb, label="Slipids")
    plt.legend()
    plt.ylabel("Q")
    plt.savefig("slipids_lipid14_charges_aligned.png", dpi=300)

