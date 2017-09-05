
[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/arose/molio/blob/master/LICENSE)

- general, non-opinionated library for reading and writing molecular structure related file formats
- extending on the ideas of the CIFTools.js library


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
- gro reader
  - read more than one block
  - read velocities
  - detect number of decimal places
