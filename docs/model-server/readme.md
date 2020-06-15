Model Server
============

Model Server is a tool for preprocessing and querying macromolecular structure data.

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
node lib/servers/model/server/server
```

## From NPM

```
npm install --production molstar
./model-server 
```

(or ``node node_modules\.bin\model-server`` in Windows).

The NPM package contains all the tools mentioned here as "binaries":

- ``model-server``
- ``model-server-query``
- ``model-server-preprocess``


### Production use

In production it is required to use a service that will keep the server running, such as [forever.js](https://github.com/foreverjs/forever).


### Memory issues

Sometimes nodejs might run into problems with memory. This is usually resolved by adding the ``--max-old-space-size=8192`` parameter.

## Preprocessor

The preprocessor application allows to add custom data to CIF files and/or convert CIF to BinaryCIF. ``node lib/commonjs/servers/model/preprocess`` or ``model-server-preprocess`` binary from the NPM package.


## Local Mode

The server can be run in local/file based mode using ``node lib/commonjs/servers/model/query`` (``model-server-query`` binary from the NPM package).

Custom Properties
=================

This feature is still in development.

It is possible to provide property descriptors that transform data to internal representation and define how it should be exported into one or mode CIF categories. Examples of this are located in the ``mol-model-props`` module and are linked to the server in the config and ``servers/model/properties``.