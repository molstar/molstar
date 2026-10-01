[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](./LICENSE)
[![npm version](https://badge.fury.io/js/molstar.svg)](https://www.npmjs.com/package/molstar)
[![Build](https://github.com/molstar/molstar/actions/workflows/node.yml/badge.svg)](https://github.com/molstar/molstar/actions/workflows/node.yml)
[![Gitter](https://badges.gitter.im/molstar/Lobby.svg)](https://gitter.im/molstar/Lobby)

# Mol*

The goal of **Mol\*** (*/'mol-star/*) is to provide a technology stack that serves as a basis for the next-generation data delivery and analysis tools for (not only) macromolecular structure data. Mol* development was jointly initiated by PDBe and RCSB PDB to combine and build on the strengths of [LiteMol](https://litemol.org) (developed by PDBe) and [NGL](https://nglviewer.org) (developed by RCSB PDB) viewers.

When using Mol*, please cite:

David Sehnal, Sebastian Bittrich, Mandar Deshpande, Radka Svobodová, Karel Berka, Václav Bazgier, Sameer Velankar, Stephen K Burley, Jaroslav Koča, Alexander S Rose: [Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures](https://doi.org/10.1093/nar/gkab314), *Nucleic Acids Research*, 2021; https://doi.org/10.1093/nar/gkab314.

### Protein Data Bank Integrations

- The [pdbe-molstar](https://github.com/molstar/pdbe-molstar) library is the Mol* implementation used by EMBL-EBI data resources such as [PDBe](https://pdbe.org/), [PDBe-KB](https://pdbe-kb.org/) and [AlphaFold DB](https://alphafold.ebi.ac.uk/). This implementation can be used as a JS plugin and a Web component and supports property/attribute-based easy customisation. It provides helper methods to facilitate programmatic interactions between the web application and the 3D viewer. It also provides a superposition view for overlaying all the observed ligand molecules on representative protein conformations.

- [rcsb-molstar](https://github.com/molstar/rcsb-molstar) is the Mol* plugin used by [RCSB PDB](https://www.rcsb.org). The project provides additional presets for the visualization of structure alignments and structure motifs such as ligand binding sites. Furthermore, [rcsb-molstar](https://github.com/molstar/rcsb-molstar) allows to interactively add or hide of (parts of) chains, as seen in the [3D Protein Feature View](https://www.rcsb.org/3d-sequence/4hhb).


## Project Structure Overview

The v6 prototype separates code into pnpm workspace packages:

- `packages/{core,io,model,graphics}` own the shared library layers.
- `packages/plugin/{core,ui,headless}` own plugin runtime, React UI and Node capture.
- `packages/mvs/{builder,runtime}` separate standalone MVS construction from plugin loading.
- `extensions/`, `apps/`, `examples/`, `servers/` and `cli/` own their dependencies and builds.
- `distributions/molstar/` assembles the classic and browser ESM distribution.
- `smoke/` checks isolated package consumers and browser rendering.

See the [workspace guide](.v6/plans/workspace-usage.md) for package APIs, ESM consumption,
versioning, commands and the [implementation plan](.v6/plans/workspace-prototype.md).

## Previous Work
This project builds on experience from previous solutions:
- [LiteMol Suite](https://www.litemol.org)
- [WebChemistry](https://webchem.ncbr.muni.cz)
- [NGL Viewer](http://nglviewer.org)
- [MMTF](http://mmtf.rcsb.org)
- [MolQL](http://molql.org)
- [PDB Component Library](https://www.ebi.ac.uk/pdbe/pdb-component-library/)
- And many others (list will be continuously expanded).

## Building & Running

Use Node 22+ and the pnpm version specified in `package.json`.

```sh
pnpm install
pnpm build
pnpm dev:viewer
```

Run `pnpm check:workspace`, `pnpm test` and `pnpm smoke` to check package boundaries,
existing behavior and packed ESM consumers. The smoke browser requires Chromium;
optional native capture has its own `pnpm smoke:headless` check.

Serve `distributions/molstar/` to access `build/viewer/` and `build/mvs-stories/`.
Browser ESM entry points are under `build/esm/`. Use `node scripts/clean.js --all`
for a clean rebuild. Detailed commands and code ownership are in the
[workspace guide](.v6/plans/workspace-usage.md).

Code generators are compiled under their owning `cli/<name>/lib/` directories.
For example:

```sh
node cli/cifschema/lib/index.js -mip @molstar/core/data -o packages/io/src/reader/cif/schema/mmcif.ts -p mmCIF
node cli/lipid-params/lib/index.js -o packages/model/src/model/structure/model/types/lipids.ts
node cli/syminfo/lib/index.js
node cli/cif2bcif/lib/index.js input.cif output.bcif
```

## Development

### Editor

To get syntax highlighting for shader files add the following to Visual Code's settings files and make sure relevant extensions are installed in the editor.

    "files.associations": {
        "*.glsl.ts": "glsl",
        "*.frag.ts": "glsl",
        "*.vert.ts": "glsl"
    },

## Publish

### Prerelease
    npm version prerelease # assumes the current version ends with '-dev.X'
    npm publish --tag next

### Release
    npm version 0.X.0 # provide valid semver string
    npm publish

## Deploy
To prepare apps and demos for https://molstar.org deploy, run:

    npm run test
    npm run deploy:local

To commit these changes remotely to the `molstar/molstar.github.io` repo:

    npm run deploy:remote

## Contributing
Just open an issue or make a pull request. All contributions are welcome.

## Funding
Funding sources include but are not limited to:
* [RCSB PDB](https://www.rcsb.org) funding by a grant [DBI-1338415; PI: SK Burley] from the NSF, the NIH, and the US DoE
* [PDBe, EMBL-EBI](https://pdbe.org)
* [CEITEC](https://www.ceitec.eu/)
* [EntosAI](https://www.entos.ai)
