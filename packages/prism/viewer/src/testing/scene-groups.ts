/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SceneNode } from '../model/scene.js';

// A structure's grouped rows (viewer-spec.md §3): its components under Components, its layers under Layers. For tests
// that read the rows under a structure.

/** A structure node's components (the Components group's rows); none: []. */
export const componentsOf = (structure: SceneNode): readonly SceneNode[] =>
  structure.children.value.find((c) => c.kind === 'components')?.children.value ?? [];

/** A structure node's layers (the Layers group's rows: the unit cell, volume streaming); none: []. */
export const layersOf = (structure: SceneNode): readonly SceneNode[] =>
  structure.children.value.find((c) => c.kind === 'layers')?.children.value ?? [];
