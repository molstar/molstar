/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The examples of Seq* (next): a few datasets, many views. Each view is a builder module over a dataset's fixtures,
// exporting `build()`; `views.ts` lists them. The Node assets adapter (`@molstar/seqx-examples/node-assets`, it imports
// `node:fs`) is an entry of its own, so a browser bundle of the views never pulls it in.

export {
  ComplexViews,
  ExampleViews,
  type ExampleView,
  FeatureViews,
  FilesViews,
  TableFiles,
  GlobinsViews,
  LayerViews,
  P53Views,
  PlasmidViews,
  TrnaScaleView,
  TrnaViews,
} from './views.js';
