"""Simple utilities to convert ligands from SDF/SMILES to PDBQT for AutoDock Vina."""

import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, SmilesMolSupplier
from meeko import MoleculePreparation, PDBQTWriterLegacy

DEFAULT_UFF_MAX_ITERS = 2000

def sanitize_name(name):
    """Turn ligand name into a safe file name."""
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in name)

def write_pdbqt(pdbqt_str, out_dir, ligand_name):
    """Write PDBQT text into <ligand_name>.pdbqt inside out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / f"{ligand_name}.pdbqt"
    with open(out_path, "w") as f:
        f.write(pdbqt_str)
    return out_path

def prepare_3d_mol(mol, max_iters=DEFAULT_UFF_MAX_ITERS, label=None):
    """
    Add hydrogens, generate a 3D conformer and optimize it with UFF.

    Returns optimized molecule or None if something failed.
    """
    tag = f" [{label}]" if label is not None else ""

    try:
        mol = Chem.AddHs(mol)
    except Exception as e:
        print(f"[3D]{tag} Adding hydrogens failed: {e}")
        return None

    mol.RemoveAllConformers()

    status_emb = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if status_emb != 0:
        print(f"[3D]{tag} Generating 3D conformer (embedding) failed")
        return None   
    
    print(f"[3D]{tag} Running UFF optimization (max_iters={max_iters})")
    
    status_opt = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    if status_opt == -1:
        print(f"[3D]{tag} UFF optimization failed")
        return None

    return mol

def mol_to_pdbqt_string(mol, prep):
    """Prepare a molecule with Meeko and return a PDBQT string, or None on failure."""
    mol_setups = prep.prepare(mol)
    if not mol_setups:
        return None

    setup = mol_setups[0]
    pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setup)
    if not is_ok:
        return None

    return pdbqt_str

def read_sdf_as_3d(sdf_path, max_iters=DEFAULT_UFF_MAX_ITERS):
    """Read molecules from an SDF file and return list of (name, 3D_mol)."""
    sdf_path = Path(sdf_path)
    sdf_supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)

    result = []

    for i, mol in enumerate(sdf_supplier):
        if mol is None:
            print(f"[SDF] Skipping invalid molecule at index {i}")
            continue

        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"ligand_{i}"

        mol_3d = prepare_3d_mol(mol, max_iters=max_iters, label=name)
        if mol_3d is None:
            print(f"[SDF] 3D preparation failed for {name}")
            continue

        result.append((name, mol_3d))

    print(f"[SDF] Total SDF molecules read (and prepared): {len(result)}")
    return result

def convert_sdf_to_pdbqt(sdf_path, out_dir, prep, max_iters=DEFAULT_UFF_MAX_ITERS):
    """Convert all ligands from SDF to individual PDBQT files."""
    out_dir = Path(out_dir)
    name_mol_list = read_sdf_as_3d(sdf_path, max_iters=max_iters)

    for name, mol_3d in name_mol_list:
        safe_name = sanitize_name(name)

        pdbqt_str = mol_to_pdbqt_string(mol_3d, prep)
        if pdbqt_str is None:
            print(f"[SDF] Error preparing PDBQT for {safe_name}")
            continue

        out_path = write_pdbqt(pdbqt_str, out_dir, safe_name)
        print(f"[SDF] Wrote {out_path}")

def read_smi_as_3d(smi_path, max_iters=DEFAULT_UFF_MAX_ITERS):
    """Read molecules from a SMILES file and return list of (name, 3D_mol)."""
    smi_path = Path(smi_path)

    smi_supplier = SmilesMolSupplier(
        str(smi_path),
        delimiter=" ",
        smilesColumn=0,
        nameColumn=1,
        titleLine=True,
    )

    result = []

    for i, mol in enumerate(smi_supplier):
        if mol is None:
            print(f"[SMI] Skipping invalid SMILES at index {i}")
            continue

        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"ligand_{i}"

        mol_3d = prepare_3d_mol(mol, max_iters=max_iters, label=name)
        if mol_3d is None:
            print(f"[SMI] 3D preparation failed for {name}")
            continue

        result.append((name, mol_3d))

    print(f"[SMI] Total SMILES molecules read (and prepared): {len(result)}")
    return result

def convert_smi_to_pdbqt(smi_path, out_dir, prep, max_iters=DEFAULT_UFF_MAX_ITERS):
    """Convert all ligands from SMILES to individual PDBQT files."""
    out_dir = Path(out_dir)
    name_mol_list = read_smi_as_3d(smi_path, max_iters=max_iters)

    for name, mol_3d in name_mol_list:
        safe_name = sanitize_name(name)

        pdbqt_str = mol_to_pdbqt_string(mol_3d, prep)
        if pdbqt_str is None:
            print(f"[SMI] Error preparing PDBQT for {safe_name}")
            continue

        out_path = write_pdbqt(pdbqt_str, out_dir, safe_name)
        print(f"[SMI] Wrote {out_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert ligands from SDF/SMILES to PDBQT for AutoDock Vina."
    )
    parser.add_argument("--sdf", type=str, help="Path to input SDF file")
    parser.add_argument("--smi", type=str, help="Path to input SMILES file")
    parser.add_argument("--out-sdf", type=str, default="sdf-pdbqt",
                        help="Output directory for SDF-based PDBQT files")
    parser.add_argument("--out-smi", type=str, default="smi-pdbqt",
                        help="Output directory for SMILES-based PDBQT files")
    parser.add_argument("--uff-max-iters", type=int, default=DEFAULT_UFF_MAX_ITERS,
                        help="Max iterations for UFF optimization")

    args = parser.parse_args()

    prep = MoleculePreparation(
        rigid_macrocycles=False,
        min_ring_size=6,
    )

    if args.sdf:
        convert_sdf_to_pdbqt(args.sdf, args.out_sdf, prep, max_iters=args.uff_max_iters)
    if args.smi:
        convert_smi_to_pdbqt(args.smi, args.out_smi, prep, max_iters=args.uff_max_iters)


if __name__ == "__main__":
    main()