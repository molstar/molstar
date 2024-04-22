# Installation

## NPM Package

```
yarn add molstar
```

or

```
npm install molstar
```

Mol* code can then be imported from the ``molstar/lib/...`` namespace, e.g.

```ts
import { PluginContext } from 'molstar/lib/mol-plugin/context';
```

## Clone from GitHub


```
git clone https://github.com/molstar/molstar.git
cd molstar
npm install
npm build
```

--------------------

For a watch task to automatically rebuild the source code on changes, run

```
npm run watch
```

or if working just with the Viewer app for better performance

```
npm run watch-viewer
```
