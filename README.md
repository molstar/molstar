
[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/arose/molio/blob/master/LICENSE)

- general, non-opinionated library for reading and writing molecular structure related file formats
- extending on the ideas of the CIFTools.js library


## Module Overview

- `mol-comp` Computation abstraction with progress tracking and cancellation support.
- `mol-data` Collections (integer based sets, inteface to columns/tables, etc.)
- `mol-math` Math related (loosely) algorithms and data structures.
- `mol-io` Parsing library. Each format is parsed into an interface that corresponds to the data stored by it.
- `mol-model` Data structures and algorithms (such as querying) for representing molecular data.
- `mol-ql` Mapping of `mol-model` to the MolQL query language spec.
- `mol-util` Useful things that do not fit elsewhere.

## Building & Running

### Build:

    npm install
    npm run build

### Build automatically on file save:

    npm run watch

### Bundle with rollup (UMD and ES6)

    npm run bundle

### Make distribution files

    npm run dist

### Build everything above

    npm run-script build && npm run-script bundle && npm run-script dist


## Example script

### Build

    npm run script

### Run

    node ./build/js/script.js


TODO
----

- write about unittest (AR)