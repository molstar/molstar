[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](./LICENSE)
[![Build Status](https://travis-ci.org/molstar/molstar-proto.svg?branch=master)](https://travis-ci.org/molstar/molstar-proto)
[![Gitter](https://badges.gitter.im/molstar/Lobby.svg)](https://gitter.im/molstar/Lobby)

# Mol*

The goal of **Mol\*** (*/'mol-star/*) is to provide a technology stack that will serve as basis for the next-generation data delivery and analysis tools for macromolecular structure data. This is a collaboration between PDBe and RCSB PDB teams and the development will be open source and available to anyone who wants to use it for developing visualisation tools for macromolecular structure data available from [PDB](https://www.wwpdb.org/) and other institutions.

This particular project is a prototype implementation of this technology (still under development).

## Project Overview

The core of Mol* currently consists of these modules:

- `mol-task` Computation abstraction with progress tracking and cancellation support.
- `mol-data` Collections (integer based sets, interface to columns/tables, etc.)
- `mol-math` Math related (loosely) algorithms and data structures.
- `mol-io` Parsing library. Each format is parsed into an interface that corresponds to the data stored by it. Support for common coordinate, experimental/map, and annotation data formats.
- `mol-model` Data structures and algorithms (such as querying) for representing molecular data (including coordinate, experimental/map, and annotation data).
- `mol-model-formats` Data format parsers for `mol-model`.
- `mol-model-props` Common "custom properties".
- `mol-script` A scriting language for creating representations/scenes and querying (includes the [MolQL query language](https://molql.github.io)).
- `mol-geo` Creating (molecular) geometries.
- `mol-theme` Molecular representation themeing.
- `mol-repr` Molecular representations.
- `mol-gl` A wrapper around WebGL.
- `mol-canvas3d` A low level 3d view component. Uses `mol-geo` to generate geometries.
- `mol-state` State representation tree with state saving and automatic updates.
- `mol-app` Components for builduing UIs.
- `mol-plugin` Allow to define modular Mol* plugin instances utilizing `mol-state` and `mol-canvas3d`.
- `mol-util` Useful things that do not fit elsewhere.

Moreover, the project contains the imlementation of `servers`, including

- `servers/model` A tool for accessing coordinate and annotation data of molecular structures.
- `servers/volume` A tool for accessing volumetric experimental data related to molecular structures.

The project also contains performance tests (`perf-tests`), `examples`, and basic proof of concept `apps` (CIF to BinaryCIF converter and JSON domain annotation to CIF converter).

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

### Build:
    npm install
    npm run build

### Build automatically on file save:
    npm run watch

### Build with debug mode enabled:
    DEBUG=molstar npm run watch

### Build for production:
    npm run build
    NODE_ENV=production npm run build-webpack

**Run**

If not installed previously:

    npm install -g http-server

...or a similar solution.

From the root of the project:

    http-server -p PORT-NUMBER

and navigate to `build/viewer`


### Code generation
**CIF schemas**

    export NODE_PATH="build/src"; node build/src/apps/schema-generator/schema-from-cif-dic.js -ts -o src/mol-io/reader/cif/schema/mmcif.ts --fieldNamesPath data/mmcif-field-names.csv --name mmCIF

    export NODE_PATH="build/src"; node build/src/apps/schema-generator/schema-from-cif-dic.js -ts -o src/mol-io/reader/cif/schema/ccd.ts --fieldNamesPath data/ccd-field-names.csv --name CCD

    export NODE_PATH="build/src"; node build/src/apps/schema-generator/schema-from-cif-dic.js -ts -o src/mol-io/reader/cif/schema/bird.ts --fieldNamesPath data/bird-field-names.csv --name BIRD

**GraphQL schemas**

    node data/rcsb-graphql/codegen.js

### Other scripts
**Create chem comp bond table**

    export NODE_PATH="build/src"; node --max-old-space-size=8192 build/src/apps/chem-comp-bond/create-table.js build/data/ccb.bcif -b

**Test model server**

    export NODE_PATH="build/src"; node build/src/servers/model/test.js

**State Transformer Docs**

    export NODE_PATH="build/src"; node build/state-docs

**Convert any CIF to BinaryCIF**

    node build/model-server/preprocess -i file.cif -ob file.bcif

To see all available commands, use ``node build/model-server/preprocess -h``.

## Development

### Intallation

If node complains about a missine acorn peer dependency, run the following commands

    npm update acorn --depth 20
    npm dedupe

If the `gl` package does not compile on node 12 (there are currently no pre-built binaries) revert back to node 10.

### Editor

To get syntax highlighting for the shader files add the following to Visual Code's settings files

    "files.associations": {
        "*.glsl.ts": "glsl",
        "*.frag.ts": "glsl",
        "*.vert.ts": "glsl"
    },

## Contributing
Just open an issue or make a pull request. All contributions are welcome.

## Roadmap
Continually develop this prototype project. As individual modules become stable, make them into standalone libraries.

## Funding
Funding sources include but are not limited to:
* [RCSB PDB](https://www.rcsb.org) funding by a grant [DBI-1338415; PI: SK Burley] from the NSF, the NIH, and the US DoE
* [PDBe, EMBL-EBI](https://pdbe.org)
* [CEITEC](https://www.ceitec.eu/)