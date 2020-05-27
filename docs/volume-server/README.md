What is VolumeServer
=====================

VolumeServer is a service for accessing subsets of volumetric density data. It automatically downsamples the data depending on the volume of the requested region to reduce the bandwidth requirements and provide near-instant access to even the largest data sets.

It uses the text based CIF and BinaryCIF formats to deliver the data to the client. 

For quick info about the benefits of using the server, check out the [examples](examples.md).

Installing and Running
=====================

Requires nodejs 8+.

## From GitHub

```
git clone https://github.com/molstar/molstar
npm install
```

Afterwards, build the project source:

```
npm run build-tsc
```

and run the server by 

```
node lib/commonjs/servers/volume/server
```

## From NPM

```
npm install --production molstar
./volume-server 
```

(or ``node node_modules\.bin\volume-server`` in Windows).

The NPM package contains all the tools mentioned here as "binaries":

- ``volume-server``
- ``volume-server-pack``
- ``volume-server-query``


### Production use

In production it is required to use a service that will keep the server running, such as [forever.js](https://github.com/foreverjs/forever).


### Memory issues

Sometimes nodejs might run into problems with memory. This is usually resolved by adding the ``--max-old-space-size=8192`` parameter.


## Preparing the Data

For the server to work, CCP4/MAP (models 0, 1, 2 are supported) input data need to be converted into a custom block format. 
To achieve this, use the ``pack`` application (``node lib/commonjs/servers/volume/pack`` or ``volume-server-pack`` binary from the NPM package).

## Local Mode

The program  ``lib/commonjs/servers/volume/pack`` (``volume-server-query`` in NPM package) can be used to query the data without running a http server.

## Navigating the Source Code

The source code is split into 2 mains parts: ``pack`` and ``server``:

- The ``pack`` part provides the means of converting CCP4 files into the internal block format.
- The ``server`` includes
  - ``query``: the main part of the server that handles a query. ``execute.ts`` is the "entry point".
  - ``algebra``: linear, "coordinate", and "box" algebra provides the means for calculations necessary to concent a user query into a menaningful response.
  - API wrapper that handles the requests.

Consuming the Data 
==================

The data can be consumed in any (modern) browser using the [ciftools library](https://github.com/molstar/ciftools) (or any other piece of code that can read text or binary CIF).

The [Data Format](DataFormat.md) document gives a detailed description of the server response format.

As a reference/example of the server usage is available in Mol* ``mol-plugin`` module.