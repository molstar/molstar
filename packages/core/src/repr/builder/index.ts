/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/core/repr/builder: the typed spec builder's engine (_spec/README.md "Builder"): `createBuilder(language)`
// reads a language's kinds, parents, ids and params from its type and builds the normalized tree. It knows no
// language.

export {
  type Builder,
  type BuilderOptions,
  type ChildKinds,
  createBuilder,
  type FormFn,
  type KindFn,
  type NodeBuilder,
  type ParamsFor,
  type Quantity,
} from './builder.js';
