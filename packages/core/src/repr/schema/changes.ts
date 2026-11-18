/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a patch changed (_spec/README.md "Patches", the result): the net difference from an old root to a new one, one
// entry a node whatever sequence of ops reached it. Declared here, beside the validation that reads it
// (`validateChanges`), and re-exported by `repr/patch`, so imports go `tree` ← `schema` ← `patch`.

import type { Address } from '../tree/index.js';

/** Which of a node's records a key is in: its params (the default) or its custom data. */
export type Where = 'params' | 'custom';

export interface ParamChange {
  readonly at: Address;
  readonly keys: readonly string[];
  readonly where: Where;
}

export interface Changes {
  /** Inserted or replaced nodes, at their addresses in the new root; their subtrees are new as a whole. */
  readonly made: readonly Address[];
  /** Keys set or unset on nodes that are not made (ref rewrites included), in the new root. */
  readonly params: readonly ParamChange[];
  /** Parents whose child list changed (ids, order, or which node an id holds), in the new root; never a made node. */
  readonly children: readonly Address[];
  /** Nodes no longer in the tree, at their addresses in the old root (the top of each removed subtree). */
  readonly removed: readonly Address[];
  /** Nodes moved or renamed: `from` in the old root, `to` in the new; a reorder within a parent is not a move. */
  readonly moved: readonly { readonly from: Address; readonly to: Address }[];
}

/** No change at all. */
export const NoChanges: Changes = Object.freeze({
  made: Object.freeze([]),
  params: Object.freeze([]),
  children: Object.freeze([]),
  removed: Object.freeze([]),
  moved: Object.freeze([]),
});
