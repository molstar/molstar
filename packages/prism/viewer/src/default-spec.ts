/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { prismViewerPluginParts } from './sources/base.js';
import type { PrismViewerSpec } from './spec.js';

// The viewer's default spec (viewer-spec.md §5): its base, the default sequence and agent (no assistant), its own
// built-in features as extensions (filled as the features are built, plan phase 3), nothing kept. Mol*'s extensions are
// the app's: they depend on the viewer.

/** A new object each call. */
export function DefaultPrismViewerSpec(): PrismViewerSpec {
  return {
    plugin: prismViewerPluginParts(),
    extensions: [],
    sequence: { behaviors: [] },
    agent: {},
  };
}
