# MolViewSpec (TypeScript/JavaScript)


MolViewSpec (*.msvj) is a JSON-based file format that is used to describe visual scenes or views used in molecular visualization. It adopts declarative data-driven approach to describe, load, render, and visually deliver molecular structures, along with 3D representations, coloring schemes, and associated structural, biological, or functional annotations.

This TypeScript/JavaScript toolkit allows for describing the information required for representing a molecular view state as data in a nested tree format that can be consumed by visualization software tools such as [Mol*](https://github.com/molstar/molstar/tree/master/src/extensions/mvs).

When using MolViewSpec, please cite:

- Adam Midlik, Sebastian Bittrich, Jennifer R Fleming, Sreenath Nair, Sameer Velankar, Stephen K Burley, Jasmine Y Young, Brinda Vallat, David Sehnal: MolViewSpec: a Mol* extension for describing and sharing molecular visualizations, Nucleic Acids Research, 2025; https://doi.org/10.1093/nar/gkaf370.
- Sebastian Bittrich, Adam Midlik, Mihaly Varadi, Sameer Velankar, Stephen K. Burley, Jasmine Y. Young, David Sehnal, Brinda Vallat: Describing and Sharing Molecular Visualizations Using the MolViewSpec Toolkit, Current Protocols, 2024; https://doi.org/10.1002/cpz1.1099.

## The Idea behind MolViewSpec

In the long run, MolViewSpec aims to re-imagine how users define molecular scenes by detaching this process from any concrete 3D viewer.

MolViewSpec's workflow is:
1. `define scene using MolViewSpec`
2. `generic state description as .msvj or .mvsx file`
3. `open in any MolViewSpec-compatible 3D viewer`

Opposed to the traditional workflow that locks users into using a specific 3D viewer, such as:
1. `define scene in Mol*`
2. `Mol*-specific state format`
3. `open only in Mol*`

## Installation

### npm
```bash
npm install molviewspec
```

## Usage

### ES6 Modules

```javascript
import { createMVSBuilder, MVSData } from 'molviewspec';

// Create a new MVS builder
const builder = createMVSBuilder();

// Build a simple molecular visualization
const root = builder.getRoot();
const download = root.download({ url: 'https://files.wwpdb.org/download/1cbs.cif' });
const parse = download.parse({ format: 'mmcif' });
const structure = parse.structure({ type: 'model' });
const component = structure.component({ selector: 'all' });
component.representation({ type: 'cartoon' });

// Get the MVS data
const mvsData = builder.getState();

// Convert to MVSJ string
const mvsjString = MVSData.toMVSJ(mvsData, 2);
console.log(mvsjString);
```


### Basic Example

```javascript
import { createMVSBuilder, MVSData } from 'molviewspec';

// Create a builder
const builder = createMVSBuilder();

// Build the visualization tree
const root = builder.getRoot();
root.download({ url: 'https://files.wwpdb.org/download/1cbs.cif' })
    .parse({ format: 'mmcif' })
    .structure({ type: 'model' })
    .component({ selector: 'all' })
    .representation({ type: 'cartoon' });

// Get the state and convert to JSON
const state = builder.getState();
const json = MVSData.toMVSJ(state, 2);

// Validate the data
const isValid = MVSData.isValid(state);
console.log('Valid MVS:', isValid);
```

### Loading and Parsing MVSJ

```javascript
import { MVSData } from 'molviewspec';

// Load from MVSJ string
const mvsjString = '{"root": {...}, "metadata": {...}}';
const mvsData = MVSData.fromMVSJ(mvsjString);

// Validate
if (MVSData.isValid(mvsData)) {
    console.log('Valid MolViewSpec data');
} else {
    console.log('Invalid MolViewSpec data');
}
```

## API Reference

### Core Classes

- **`createMVSBuilder()`** - Creates a new MVS builder instance
- **`MVSData`** - Main data structure and utilities for MolViewSpec
- **`Root`** - Root node of the MVS tree
- **`MVSTree`** - Tree structure definitions

### Utilities

- **JSON utilities** - For handling JSON serialization
- **Color utilities** - For color manipulation and conversion
- **Object utilities** - For object manipulation helpers
- **Tree utilities** - For tree traversal and manipulation

### Available Exports

The library provides the following main exports:

```javascript
// Core API
import {
    MVSData,
    createMVSBuilder,
    Root,
    GlobalMetadata,
    SnapshotMetadata,
    Snapshot
} from 'molviewspec';

// Tree schema system
import {
    TreeSchema,
    validateTree,
    ParamsSchema,
    FieldSchema
} from 'molviewspec/tree/generic/tree-schema';

// MVS-specific trees
import {
    MVSTree,
    MVSBuilder
} from 'molviewspec/tree/mvs/mvs-tree';

// Utilities
import { HexColor, NamedColor, decodeColor } from 'molviewspec/util/color';
import { safePromise, isDefined } from 'molviewspec/util/helpers';
```


## MolViewSpec Tree Structure

The MolViewSpec format describes molecular visualizations as a tree structure. Here's an example:

```json
{
  "root": {
    "kind": "root",
    "children": [
      {
        "kind": "download",
        "params": {
          "url": "https://files.wwpdb.org/download/1cbs.cif"
        },
        "children": [
          {
            "kind": "parse",
            "params": {
              "format": "mmcif"
            },
            "children": [
              {
                "kind": "structure",
                "params": {
                  "type": "model"
                },
                "children": [
                  {
                    "kind": "component",
                    "params": {
                      "selector": "all"
                    },
                    "children": [
                      {
                        "kind": "representation",
                        "params": {
                          "type": "cartoon"
                        }
                      }
                    ]
                  }
                ]
              }
            ]
          }
        ]
      }
    ]
  },
  "metadata": {
    "version": "1.0.0",
    "timestamp": "2023-11-16T11:41:07.421220"
  }
}
```

## Related Projects

- [MolViewSpec Python](https://github.com/molstar/mol-view-spec) - Python implementation
- [Mol*](https://github.com/molstar/molstar) - Reference implementation for viewing MolViewSpec files
- [Mol* MVS Extension](https://github.com/molstar/molstar/tree/master/src/extensions/mvs) - Official Mol* extension for MolViewSpec
