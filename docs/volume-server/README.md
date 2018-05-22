What is VolumeServer
=====================

VolumeServer is a service for accessing subsets of volumetric density data. It automatically downsamples the data depending on the volume of the requested region to reduce the bandwidth requirements and provide near-instant access to even the largest data sets.

It uses the text based CIF and BinaryCIF formats to deliver the data to the client. 

For quick info about the benefits of using the server, check out the [examples](examples.md).

Installing the Server 
=====================

- Install [Node.js](https://nodejs.org/en/) (tested on Node 6.* and 7.*; x64 version is strongly preferred).
- Get the code.
- Prepare the data.
- Run the server.

Preparing the Data
------------------

For the server to work, CCP4/MAP (models 0, 1, 2 are supported) input data need to be converted into a custom block format. 
To achieve this, use the ``pack`` application.

- To prepare data from x-ray based methods, use: 

    ```
    node pack -xray main.ccp4 diff.ccp4 out.mdb
    ```

- For EM data, use:

    ```
    node pack -em em.map out.mdb
    ```

Running the Server
------------------

- Install production dependencies:

   ```
   npm install --only=production
   ```

- Update ``server-config.js`` to link to your data and optionally tweak the other parameters.

- Run it:

    ```
    node server
    ```

    In production it is a good idea to use a service that will keep the server running, such as [forever.js](https://github.com/foreverjs/forever).

### Local Mode

The program ``local`` in the build folder can be used to query the data without running a http server.

- ``node local`` prints the program usage.
- ``node local jobs.json`` takes a list of jobs to execute in JSON format. A job entry is defined by this interface:

    ```TypeScript
    interface JobEntry {
        source: {
            filename: string,    
            name: string,
            id: string
        },
        query: {
            kind: 'box' | 'cell',
            space?: 'fractional' | 'cartesian',
            bottomLeft?: number[],
            topRight?: number[],
        }
        params: {
            /** Determines the detail level as specified in server-config */
            detail?: number,
            /** 
             * Determines the sampling level:
             * 1: Original data
             * 2: Downsampled by factor 1/2
             * ...
             * N: downsampled 1/2^(N-1)
             */
            forcedSamplingLevel?: number,
            asBinary: boolean,
        },
        outputFolder: string
    }
    ```

    Example ``jobs.json`` file content:

    ```TypeScript
    [{
        source: {
            filename: `g:/test/mdb/emd-8116.mdb`,
            name: 'em',
            id: '8116',
        },
        query: {
            kind: 'cell'
        },
        params: {
            detail: 4,
            asBinary: true
        },
        outputFolder: 'g:/test/local-test'
    }]
    ```

## Navigating the Source Code

The source code is split into 2 mains parts: ``pack`` and ``server``:

- The ``pack`` part provides the means of converting CCP4 files into the internal block format.
- The ``server`` includes
  - ``query``: the main part of the server that handles a query. ``execute.ts`` is the "entry point".
  - ``algebra``: linear, "coordinate", and "box" algebra provides the means for calculations necessary to concent a user query into a menaningful response.
  - API wrapper that handles the requests.

Consuming the Data 
==================

The data can be consumed in any (modern) browser using the [CIFTools.js library](https://github.com/dsehnal/CIFTools.js) (or any other piece of code that can read text or binary CIF).

The [Data Format](DataFormat.md) document gives a detailed description of the server response format.

As a reference/example of the server usage, please see the implementation in [LiteMol](https://github.com/dsehnal/LiteMol) ([CIF.ts + Data.ts](https://github.com/dsehnal/LiteMol/tree/master/src/lib/Core/Formats/Density), [UI](https://github.com/dsehnal/LiteMol/tree/master/src/Viewer/Extensions/DensityStreaming)) or in Mol*.