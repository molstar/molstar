/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SceneNode } from '../model/scene.js';

// A structure's grouped rows (viewer-spec.md §3): its components under Components, the others under Other. For tests
// that read the rows under a structure.

/** A structure node's components (the Components group's rows); none: []. */
export const componentsOf = (structure: SceneNode): readonly SceneNode[] =>
  structure.children.value.find((c) => c.kind === 'components')?.children.value ?? [];

/** A structure node's other rows (the Other group's: the unit cell, volume streaming); none: []. */
export const otherOf = (structure: SceneNode): readonly SceneNode[] =>
  structure.children.value.find((c) => c.kind === 'other')?.children.value ?? [];
