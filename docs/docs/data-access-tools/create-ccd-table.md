# Create Table from CCD
The [Chemical Component Dictionary (CCD)](https://www.wwpdb.org/data/ccd) is as an external reference file describing 
all residue and small molecule components found in PDB entries. The 
[Protonation Variants Companion Dictionary (PVCD)](https://www.wwpdb.org/data/ccd) enumerates protonation variants of
canonical amino acids.

This script bundles all `chem_comp_bond` information from the CCD and the PVCD into a single file for later use.
Optionally, it can also generate a second output file that contains all `chem_comp_atom` information.

## Example
```sh
node --max-old-space-size=4096 lib/commonjs/cli/chem-comp-dict/create-table.js build/data/ccb.bcif -b
```

## Usage
| Argument | Description |
| --- | --- |
| `out` | Generated file output path |
| `--forceDownload`, `-f` | Force download of CCD and PVCD |
| `--binary`, `-b` | Output as BinaryCIF |
| `--ccaOut`, `-a` | File output path of optionally generated chem_comp_atom |

```sh
create-table.js [-h] [--forceDownload] [--binary] [--ccaOut CCAOUT] out
```