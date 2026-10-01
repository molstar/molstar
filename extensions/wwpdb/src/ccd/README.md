# Chemical Component Dictionary Extension

The [Chemical Component Dictionary (CCD)](https://www.wwpdb.org/data/ccd) describes all small molecules and monomers found in PDB entries. The dictionary provides a plethora of additional information not present in wwPDB archive structures such as chemical descriptors (SMILES & InChI) and stereochemical assignments, information on bond order and more. Most notably, the CCD provides 2 sets of coordinates:
- `ideal`: idealized/minimized coordinates, obtained using Molecular Networks' Corina, and if there are issues, OpenEye's OMEGA
- `model`: coordinates extracted from an archive structure

## How to Load a Component from URL
1. "Download Structure" -- switch "Source" to "URL"
2. Enter URL of component, e.g. https://files.rcsb.org/ligands/view/HEM.cif, leave "Format" as is
3. Click "Apply"

This parses the corresponding component into 2 models (1st: `ideal` coordinates, 2nd: `model` coordinates) and applies the default representaiton to the 1st model. `model` coordinates are available as 2nd model. Click the canvas to re-focus if you don't see anything after switching models due to the coordinates being far away.

## How to Visualize Components
There's a dedicated representation preset that faciliates the comparison of `ideal` and `model` coordinates.

1. Load a component as described above
2. Switch structure preset to "Chemical Component" (button in the top-right, in the "Structure" panel)

This creates a dedicated component for `ideal` as well as `model` coordinates and represents them as ball-and-stick. Initially, only `ideal` coordinates are shown. After toggling the visibility of `model` coordinates, they appear superimposed with the `ideal` coordinates. 

## Examples & Test Cases
Ligand | Description | Details
-- | -- | --
https://files.rcsb.org/ligands/view/HEM.cif | metal coordination |
https://files.rcsb.org/ligands/view/FE.cif | +3 oxidation state |
https://files.rcsb.org/ligands/view/FE2.cif | +2 oxidation state |
https://files.rcsb.org/ligands/view/RUC.cif | transition metal | 
https://files.rcsb.org/ligands/view/SF4.cif | Fe-S cluster | doesn't align nicely
https://files.rcsb.org/ligands/view/TBR.cif | coords identical | 
https://files.rcsb.org/ligands/view/OER.cif | coords identical | 
https://files.rcsb.org/ligands/view/FEA.cif | charges |
https://files.rcsb.org/ligands/view/PR2.cif | orientation differs | 
https://files.rcsb.org/ligands/view/03R.cif | some atoms missing |
https://files.rcsb.org/ligands/view/02U.cif | many atoms missing |
https://files.rcsb.org/ligands/view/HC0.cif | no ideal coords | unrelated: O and H atoms clashing
https://files.rcsb.org/ligands/view/Q6O.cif | no model coords |
https://files.rcsb.org/ligands/view/H0C.cif | big ligand |
https://files.rcsb.org/ligands/view/2NC.cif | dual representation as PRD and CC |
https://files.rcsb.org/birds/view/PRDCC_000001.cif | PRDCC |
https://raw.githubusercontent.com/wwPDB/extended-wwPDB-identifier-examples/main/CCD/BB87Q.cif | extended CCD identifier |
https://raw.githubusercontent.com/wwPDB/extended-wwPDB-identifier-examples/main/CCD/7ZTVU.cif | extended CCD identifier |
https://raw.githubusercontent.com/wwPDB/extended-wwPDB-identifier-examples/main/CCD/9QRZS.cif | extended CCD identifier |
https://raw.githubusercontent.com/wwPDB/extended-wwPDB-identifier-examples/main/CCD/9ABCD.cif | extended CCD identifier |
https://files.rcsb.org/ligands/view/UNK.cif | CCD special: unknown amino acid | unrelated: some model H are placed far away
https://files.rcsb.org/ligands/view/UNX.cif | CCD special: unknown atom/ion | no ideal coordinates
https://files.rcsb.org/ligands/view/UNL.cif | CCD special: unknown ligand | no coordinates whatsoever