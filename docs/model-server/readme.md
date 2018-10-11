Model Server
============

Model Server is a tool for preprocessing and querying macromolecular structure data.

Installing and Running
=====================

Getting the code (use node 8+):
```
git clone https://github.com/molstar/molstar-proto
npm install
```

Customize configuration at ``src/server/model/config.ts`` to point to your data and which custom properties to include (see the [Custom Properties](#custom-properties) section). Alternatively, the config can be edited in the compiled version in ``build/node_modules/servers/model/config.js``.

Afterwards, build the project:

```
npm run build
```

(or run watch mode for automatic rebuilds: ``npm run watch``)

Running the server locally for testing:
```
npm run model-server
```
or
```
node build/node_modules/servers/model/server
```

In production it is a good idea to use a service that will keep the server running, such as [forever.js](https://github.com/foreverjs/forever).


## Memory issues

Sometimes nodejs might run into problems with memory. This is usually resolved by adding the ``--max-old-space-size=8192`` parameter.

Preprocessor
============

The preprocessor application allows to add custom data to CIF files and/or convert CIF to BinaryCIF. See the [Custom Properties](#custom-properties) section for providing custom properties.

## Usage

The app works in two modes: single files and folders.

Single files:

```
node build\node_modules\servers\model\preprocess -i input.cif [-oc output.cif] [-ob output.bcif] [--cfg config.json]
```

Folder: 
```
node build\node_modules\servers\model\preprocess -fin input_folder [-foc output_cif_folder] [-fob output_bcif_folder] [--cfg config.json]
```

## Config

The config speficies the maximum number of processes to use (in case of folder processing) and defines sources and parameters for custom properties.

Example:
```json
{
    "numProcesses": 4,
    "customProperties": {
        "sources": [
            "./properties/pdbe"
        ],
        "params": {
            "PDBe": {
                "UseFileSource": false,
                "API": {
                    "residuewise_outlier_summary": "https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry",
                    "preferred_assembly": "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary",
                    "struct_ref_domain": "https://www.ebi.ac.uk/pdbe/api/mappings/sequence_domains"
                }
            }
        }
    }
}
```

Custom Properties
=================

It is possible to provide property descriptors that transform data to internal representation and define how it should be exported into one or mode CIF categories. Examples of this are located in the ``mol-model-props`` module and are linked to the server in the config and ``servers/model/properties``.

Local Mode
==========

The server can be run in local/file based mode:

```
node build/node_modules/servers/model/server jobs.json
```

where ``jobs.json`` is an array of 

```ts
type LocalInput = {
    input: string,
    output: string,
    query: QueryName,
    modelNums?: number[],
    params?: any,
    binary?: boolean
}[]
```

For example

```json
[
  {
    "input": "c:/test/quick/1tqn.cif",
    "output": "c:/test/quick/localapi/1tqn_full.cif",
    "query": "full"
  },
  {
    "input": "c:/test/quick/1tqn.cif",
    "output": "c:/test/quick/localapi/1tqn_full.bcif",
    "query": "full",
    "params": {}
  },
  {
    "input": "c:/test/quick/1cbs_updated.cif",
    "output": "c:/test/quick/localapi/1cbs_ligint.cif",
    "query": "residueInteraction",
    "params": {
      "atom_site": { "label_comp_id": "REA" }
    }
  }
]
```