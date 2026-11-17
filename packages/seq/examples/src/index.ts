/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The examples of Seq* (next): a few datasets, many views. Each view is a builder module over a dataset's fixtures,
// exporting `build()`; `views.ts` lists them. The Node assets adapter reads the fixtures from the package folder.

export { nodeFixtureAssets, PackageFolder } from './node-assets.js';
export { ExampleViews, type ExampleView, P53Views } from './views.js';
