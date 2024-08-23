import sys
import itertools
from copy import deepcopy

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import BondType
from rdkit import rdBase, RDLogger


def num_rad_elecs(molec: rdkit.Chem.Mol) -> int:
    """Count the total number of radical electrons in a molecule.

    @param molec: rdkit Mol object of the target molecule.
    """
    return sum([a.GetNumRadicalElectrons() for a in molec.GetAtoms()])


def gen_bond_combs(n_bonds: int, n_conns: int) -> list:
    """Generates all possible combinations of single double and triple bonds.

    @param n_bonds: Number of shared pairs of electrons: A double bond counts
        as two.
    @param n_conns: Number of connections between atoms: A double bond counts
        as one.

    The only condition to fulfill at this step is that the sum of each of
    the combinations has to be equal to n_bonds. This accepts configurations
    breaking the octet.
    """
    iters = [[1, 2, 3]] * n_conns
    combs = []
    for comb in itertools.product(*iters):
        if sum(comb) == n_bonds:
            combs.append(comb)
    return combs


def filter_valid_structs(mol: rdkit.Chem.Mol, combs: list, hvy_bond_ids: list) -> list:
    """Removes all configurations exceeding the octet rule, radicals are allowed.

    @param mol: rdkit Mol object of the target molecule.
    @param combs: Possible combinations of bond configuration for the molecule.
    @param hvy_bond_ids: List of bonds connecting non-H atoms (C-C, C-O, etc.)
    """
    # Disable logging
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.ERROR)
    rdBase.DisableLog('rdApp.error')

    atomic_valences = {'C': [4], 
                       'N': [3, 4, 5],
                       'O': [2], 
                       'H': [1], 
                       'S': [2, 4, 6],
                       'F': [1], 
                       'Cl': [1, 3, 5, 7], 
                       'Br': [1, 3, 5, 7],
                       'I': [1, 3, 5, 7],
                       'Xe': [0]}
    valid_mols = []
    for comb in combs:
        new_mol = deepcopy(mol)
        for i, bo in enumerate(comb):
            bond = new_mol.GetBondWithIdx(hvy_bond_ids[i])
            bond.SetBondType(BondType(bo))

        # Don't include unphysical structures/
        try:
            new_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(new_mol, Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                             | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                             | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                             | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                             | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                             | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                             catchErrors=True)
        except rdkit.Chem.rdchem.AtomSanitizeException:
            continue
        if any([a.GetExplicitValence() > max(atomic_valences[a.GetSymbol()]) + a.GetFormalCharge() 
                for a in new_mol.GetAtoms()]):
            continue

        # Correct radical centers
        valid_comb = False
        for a in new_mol.GetAtoms():
            for val in sorted(atomic_valences[a.GetSymbol()]):
                if a.GetExplicitValence() > val + a.GetFormalCharge():
                    valid_comb = False
                    continue
                else:
                    real_val = val
                    valid_comb = True
                    break
            if not valid_comb:
                break
            n_rad_elecs = real_val - a.GetExplicitValence() - a.GetFormalCharge()
            n_rad_elecs = max(n_rad_elecs, 0)
            a.SetNumRadicalElectrons(n_rad_elecs)
        if not valid_comb:
            continue
        if not any([Chem.MolToSmiles(new_mol) == Chem.MolToSmiles(vmol)
                    for vmol in valid_mols]):
            valid_mols.append(new_mol)
    return valid_mols


def gen_reso_structs(smi: str, min_rads=True) -> list:  # C(=C\\1/[C]C1)\\[CH2]
    """Generate all possible resonant structures of a given molecule.

    @param smi: SMILES code of the given molecule
    @param min_rads: Whether to minimize the number of radical electrons or not.
    """
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                     | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                     | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                     | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                     | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                     | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                     catchErrors=True)
    if not mol:
        raise RuntimeError(f'Unable to make rdkit mol structure for {smi}.')
    mol = Chem.AddHs(mol)
    hvy_bond_ids = [b.GetIdx() for b in mol.GetBonds()
                    if b.GetBeginAtom().GetSymbol() != 'H'
                    and b.GetEndAtom().GetSymbol() != 'H']
    num_bonds = int(sum([mol.GetBondWithIdx(bid).GetBondTypeAsDouble() for bid in hvy_bond_ids]))
    num_conns = len(hvy_bond_ids)
    radic_elecs = num_rad_elecs(mol)
    max_bonds = num_bonds + radic_elecs // 2
    valid_structs = []
    for new_bonds in range(num_bonds, max_bonds + 1):
        combs = gen_bond_combs(new_bonds, num_conns)
        valid_structs += filter_valid_structs(mol, combs, hvy_bond_ids)

    if not valid_structs:
        raise RuntimeError

    # Remove Structures with larger number of radical electrons
    if min_rads:
        min_num_rads = min([num_rad_elecs(struct) for struct in valid_structs])
        valid_structs = [struct for struct in valid_structs
                         if num_rad_elecs(struct) == min_num_rads]

    # Remove explicit Hydrogen atoms
    for s, struct in enumerate(valid_structs):
        valid_structs[s] = Chem.RemoveHs(struct, sanitize=False)

    return valid_structs


if __name__ == '__main__':
    smiles = sys.argv[1]
    valid_strs = gen_reso_structs(smiles)
    for j, m in enumerate(valid_strs):
        Draw.MolToFile(m, f'aaa{j}.png')
