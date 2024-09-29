# 0.9.12
* add `health-check` endpoint + `healthCheckPath` config prop to report service health

# 0.9.11
# SDF/MOL2 ligand export: fix atom indices when additional atoms are present

# 0.9.10
* /ligand queries: fix atom count reported by SDF/MOL/MOL2 export

# 0.9.9
* /ligand queries: fix behavior for alternate locations
* /ligand queries: handle additional atoms more gracefully
* /ligand queries: better error message for UNL
* /ligand queries: treat deuterium/tritium as hydrogen

# 0.9.8
* fix support for chem_comp_bond and struct_conn categories

# 0.9.7
* add Surrounding Ligands query

# 0.9.6
* optional download parameter

# 0.9.5
* Support molstar_global_model_transform_info category.

# 0.9.4
* bug fix for /ligand queries on metal ions

# 0.9.3
* optional transform parameter

# 0.9.2
* assemblyName in /residueInteraction
* /ligand query
* additional export encoding formats

# 0.9.1
* query-many
* Config overhaul

# 0.9.0
* REST API support.
* Swagger UI support.
* Response schemas.
* Bug fixes.
* Refactored config which can now be provided as a seprate JSON file.

# 0.8.0
* Let's call this an initial version.