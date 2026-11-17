/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `hover` behavior (_spec/README.md "Behaviors"): pointer moves to the hover model, through the bindings. A small
// policy: what the hover gesture is bound to on the hit's target decides; a move over nothing, or a leave, clears
// what this pointer owns and nothing else (EC-141). Per move it allocates nothing.

import type { BindingTable, ModifierState } from '../bindings.js';
import { type Hit, type HoverModel, hitTarget } from '../models/hover.js';

export interface HoverBehavior {
  /** This pointer's owner token: a host's other hover sources (a 3D view) use their own. */
  readonly owner: object;
  /** A pointer move over `hit` (undefined: over nothing). */
  move(hit: Hit | undefined, modifiers?: ModifierState): void;
  /** The pointer left: clears the hover if this pointer still owns it. */
  leave(): void;
}

const NoModifiers: ModifierState = {};

export function hoverBehavior(hover: HoverModel, bindings: BindingTable, owner: object = {}): HoverBehavior {
  return {
    owner,
    move(hit, modifiers = NoModifiers) {
      if (hit && bindings.resolve('hover', modifiers, hitTarget(hit)) === 'hover') hover.set(hit, owner);
      else hover.clear(owner);
    },
    leave() {
      hover.clear(owner);
    },
  };
}
