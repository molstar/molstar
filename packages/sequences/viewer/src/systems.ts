/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Coordinate systems at the edge (seqspec:interaction#focus; seqspec:libraries#the-system-registry): `ctx.systems`,
// core's `SystemRegistry` over the shown set (`graph.assembled`) with a reactive `version`, bumped by every
// registration, so channels holding `Coordinates` light them again. A host's structure link registers with it as with
// any registry; nothing here is the viewer's but the version.

import { SystemRegistry } from '@molstar/seq-core/translate';
import { state } from './state.js';

export class ViewerSystems extends SystemRegistry {
  /** Bumped by every registration: channels holding `Coordinates` light them again. */
  readonly version = state(0);

  constructor() {
    super();
    this.onChange = () => this.version.set(this.version.value + 1);
  }
}
