# VolumeServer

## What is VolumeServer

Provides near-instantaneous access to volumetric data including density maps (for instance, from X-ray crystallography 
or cryo-electron microscopy experiments), spatial distribution data, output from electrostatic calculations. It works by
utilizing adaptive downsampling (similar to how Google Earth works). 

It uses the text based CIF and BinaryCIF formats to deliver the data to the client. 

For quick info about the benefits of using the server, check out the [examples](examples.md).

## Installing and Running

Requires nodejs 8+.

### From GitHub

```
git clone https://github.com/molstar/molstar
npm install
```

Afterwards, build the project source:

```
npm run build:lib
```

and run the server by 

```
node lib/commonjs/servers/volume/server
```

### From NPM

```
npm install --production molstar
./volume-server 
```

(or ``node node_modules\.bin\volume-server`` in Windows).

The NPM package contains all the tools mentioned here as "binaries":

- ``volume-server``
- ``volume-server-pack``
- ``volume-server-query``


#### Production use

In production it is required to use a service that will keep the server running, such as [forever.js](https://github.com/foreverjs/forever).


#### Memory issues

Sometimes nodejs might run into problems with memory. This is usually resolved by adding the ``--max-old-space-size=8192`` parameter.


### Preparing the Data

For the server to work, CCP4/MAP (models 0, 1, 2 are supported) input data need to be converted into a custom block format. 
To achieve this, use the ``pack`` application (``node lib/commonjs/servers/volume/pack`` or ``volume-server-pack`` binary from the NPM package).

### Local Mode

The program  ``lib/commonjs/servers/volume/query`` (``volume-server-query`` in NPM package) can be used to query the data without running a http server.

### Navigating the Source Code

The source code is split into 2 mains parts: ``pack`` and ``server``:

- The ``pack`` part provides the means of converting CCP4 files into the internal block format.
- The ``server`` includes
  - ``query``: the main part of the server that handles a query. ``execute.ts`` is the "entry point".
  - ``algebra``: linear, "coordinate", and "box" algebra provides the means for calculations necessary to concent a user query into a menaningful response.
  - API wrapper that handles the requests.

## Consuming the Data 


The data can be consumed in any (modern) browser using the [ciftools library](https://github.com/molstar/ciftools) (or any other piece of code that can read text or binary CIF).

The [Data Format](./response-data-format.md) document gives a detailed description of the server response format.

As a reference/example of the server usage is available in Mol* ``mol-plugin`` module.

## Hosting the server

### Example

```sh
node lib/commonjs/servers/volume/server --idMap x-ray '/opt/data/xray/${id}.mdb'
```

### Usage
| Argument= | Description |
| --- | --- |
| `--cfg` | JSON config file path. If a property is not specified, cmd line param/OS variable/default value are used. |
| `--printCfg` | Print current config for validation and exit. |
| `--cfgTemplate` | Prints default JSON config template to be modified and exit. |
| `--apiPrefix` | Specify the prefix of the API, i.e. &lt;host&gt;/&lt;apiPrefix&gt;/&lt;API queries&gt; |
| `--defaultPort` | Specify the port the server is running on |
| `--shutdownTimeoutMinutes` | Server will shut down after this amount of minutes, 0 for off. |
| `--shutdownTimeoutVarianceMinutes` | Modifies the shutdown timer by +/- `timeoutVarianceMinutes` (to avoid multiple instances shutting at the same time) |
| `--idMap` | Map `id`s for a `type` to a file path. Example: `x-ray '../../data/mdb/xray/${id}-ccp4.mdb'` - JS expressions can be used inside `${}`, e.g. `${id.substr(1, 2)}/${id}.mdb` - Can be specified multiple times. - The `TYPE` variable (e.g. `x-ray`) is arbitrary and depends on how you plan to use the server. By default, Mol* Viewer uses `x-ray` and `em`, but any particular use case may vary. - If using URL, it can be http://, https://, gs:// or file:// protocol.|
| `--maxRequestBlockCount` | Maximum number of blocks that could be read in 1 query. This is somewhat tied to the ``maxOutputSizeInVoxelCountByPrecisionLevel`` in that the `&lt;maximum number of voxel&gt; = maxRequestBlockCount * &lt;block size&gt;^3`. The default block size is 96 which corresponds to 28,311,552 voxels with 32 max blocks. |
| `--maxFractionalBoxVolume` | The maximum fractional volume of the query box (to prevent queries that are too big). |
| `--maxOutputSizeInVoxelCountByPrecisionLevel` | What is the (approximate) maximum desired size in voxel count by precision level - Rule of thumb: `&lt;response gzipped size&gt;` in `[&lt;voxel count&gt; / 8, &lt;voxel count&gt; / 4]`. The maximum number of voxels is tied to maxRequestBlockCount. |

```sh
node lib/commonjs/servers/volume/server [-h] [-v]
        [--cfg CFG]
        [--printCfg] 
        [--cfgTemplate]
        [--apiPrefix PREFIX]
        [--defaultPort PORT]
        [--shutdownTimeoutMinutes TIME]
        [--shutdownTimeoutVarianceMinutes VARIANCE] 
        [--idMap TYPE PATH]
        [--maxRequestBlockCount COUNT] 
        [--maxFractionalBoxVolume VOLUME]
        [--maxOutputSizeInVoxelCountByPrecisionLevel LEVEL [LEVEL ...]]
```