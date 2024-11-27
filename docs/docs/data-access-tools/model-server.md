# Model Server

Provides access to molecular 1D, 2D, and 3D (sub-)structure models of molecules. Substructures are described by the 
mol-script (MolQL) language. It has the ability to include additional data to mmCIF “on the fly”, e.g. integrate 
primary PDB archival data from [Chemical Component Dictionary (CCD)](https://www.wwpdb.org/data/ccd), 
[Protonation Variants Companion Dictionary (PVCD)](https://www.wwpdb.org/data/ccd) and 
[Biologically Interesting moleculeReference Dictionary (BIRD)](https://www.wwpdb.org/data/bird). 

## Example
```sh
node lib/commonjs/servers/model/server --sourceMap pdb-bcif '/opt/data/bcif/${id}.bcif'
```

## Usage
| Argument | Description |
| --- | --- |
| `--version`, `-v` | Show program's version number and exit. |
| `--cfg` | JSON config file path. If a property is not specified, cmd line param/OS variable/default value are used. |
| `--printCfg` | Print current config for validation and exit. |
| `--cfgTemplate` | Prints default JSON config template to be modified and exit. |
| `--apiPrefix` | Specify the prefix of the API, i.e. &lt;host&gt;/&lt;apiPrefix&gt;/&lt;API queries&gt; |
| `--defaultPort` | Specify the port the server is running on |
| `--cacheMaxSizeInBytes` | Read structures are cached, this specifies the cache size, 0 for off. |
| `--cacheEntryTimeoutMs` | Specify in ms how long to keep entries in cache. |
| `--requestTimeoutMs` | The maximum number of ms the server spends on a request. |
| `--queryTimeoutMs` | The maximum time the server dedicates to executing a query in ms. Does not include the time it takes to read and export the data. |
| `--shutdownTimeoutMinutes` | Server will shut down after this amount of minutes, 0 for off. |
| `--shutdownTimeoutVarianceMinutes` | Modifies the shutdown timer by +/- `timeoutVarianceMinutes` (to avoid multiple instances shutting at the same time) |
| `--maxQueryManyQueries` | Maximum number of queries allowed by the query-many at a time |
| `--defaultSource` | modifies which 'sourceMap' source to use by default |
| `--sourceMap` | Map `id`s for a `source` to a file path. Example: `pdb-bcif '../../data/bcif/${id}.bcif'` - JS expressions can be used inside `${}`, e.g. `${id.substr(1, 2)}/${id}.mdb` Can be specified multiple times. The `SOURCE` variable (e.g. `pdb-bcif`) is arbitrary and depends on how you plan to use the server. Supported formats: cif, bcif, cif.gz, bcif.gz |
| `--sourceMapUrl` | Same as `--sourceMap` but for URL. `--sourceMapUrl src url format` Example: `pdb-cif 'https://www.ebi.ac.uk/pdbe/entry-files/download/${id}_updated.cif' cif` Supported formats: cif, bcif, cif.gz, bcif.gz. Supported protocols: http://, https://, gs:// |

```sh
node lib/commonjs/servers/model/server [-h] [-v]
        [--cfg CFG]
        [--printCfg]
        [--cfgTemplate]
        [--apiPrefix PREFIX]
        [--defaultPort PORT]
        [--cacheMaxSizeInBytes CACHE_SIZE]
        [--cacheEntryTimeoutMs CACHE_TIMEOUT]
        [--requestTimeoutMs REQUEST_TIMEOUT]
        [--queryTimeoutMs QUERY_TIMEOUT]
        [--shutdownTimeoutMinutes TIME]
        [--shutdownTimeoutVarianceMinutes VARIANCE]
        [--maxQueryManyQueries QUERY_MANY_LIMIT]
        [--defaultSource DEFAULT_SOURCE]
        [--sourceMap SOURCE PATH]
        [--sourceMapUrl SOURCE PATH SOURCE_MAP_FORMAT]
```

### Production Use
In production, it is required to use a service that will keep the server running, such as [forever.js](https://github.com/foreverjs/forever).

### Memory Issues
Sometimes nodejs might run into problems with memory. This is usually resolved by adding the ``--max-old-space-size=8192`` parameter.

### Preprocessor Example
The preprocessor application allows addiing custom data to CIF files and/or 
[convert CIF to BinaryCIF](./convert-to-bcif.md).
```sh
node lib/commonjs/servers/model/preprocess
```

### Preprocessor Usage
| Argument | Description |
| --- | --- |
| `--input`, `-i` | Input filename |
| `--outCIF`, `-oc` | Output CIF filename |
| `--outBCIF`, `-ob` | Output BinaryCIF filename |
| `--cfg`, `-c` | Config file path |
| `--folderIn`, `-fin` | Convert folder |
| `--folderOutCIF`, `-foc` | Convert folder text output |
| `--folderOutBCIF`, `-fob` | Convert folder binary output |
| `--folderNumProcesses`, `-fp` | Convert folder number processes |

Example cfg.json:
```ts                       
{ 
    "numProcesses": 1, 
    "customProperties": { 
        "sources": [ "wwpdb" ], 
        "params": { 
            "wwPDB": { 
                "chemCompBondTablePath": "./build/data/ccb.bcif"
            }
        }
    }
}
```

### Local Mode
The server can be run in local/file based mode using
```sh
node lib/commonjs/servers/model/query
```

### Custom Properties
This feature is still in development.

It is possible to provide property descriptors that transform data to internal representation and define how it should 
be exported into one or mode CIF categories. Examples of this are located in the ``mol-model-props`` module and are 
linked to the server in the config and ``servers/model/properties``.

## From NPM

```
npm install --production molstar
cd ./model-server 
```

(or ``node node_modules\.bin\model-server`` in Windows).

The NPM package contains all the tools mentioned in the previous sections as "binaries":

- ``model-server``
- ``model-server-query``
- ``model-server-preprocess``
