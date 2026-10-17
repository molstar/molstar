/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { State } from '@molstar/prism-ui/reactive/state';

// The scene's structure, as one change stream (view-models.md §3, collections): a view of the whole tree (the scene
// tree's flattening) listens to this once, not to every node's `children` and `isExpanded`, so what it holds is one
// listener whatever the scene's size. A change names the nodes whose children or expansion changed, by id, once per Mol*
// transaction (the source's sync) or per toggle; the roots are their own state (`SceneModel.roots`).

/** What changed in the scene's structure since the last change told. */
export interface SceneChanges {
  /** Counts up, one a change. */
  readonly version: number;
  /** The nodes (by id) whose children or expansion changed. */
  readonly nodes: ReadonlySet<string>;
}

/** The source's side: nodes touched in a batch (a sync), told once at its end. */
export class SceneChangeLog {
  readonly #state = new State<SceneChanges>({ version: 0, nodes: new Set() });
  #pending = new Set<string>();
  #depth = 0;

  /** The change stream (`SceneModel.changes`). */
  get changes(): ReadonlyState<SceneChanges> {
    return this.#state.asReadonly;
  }

  /** `id`'s children or expansion changed: told at the end of the batch, or now outside one. */
  touch(id: string): void {
    this.#pending.add(id);
    if (this.#depth === 0) this.#tell();
  }

  /** What `fn` touches is told once, at its end (a sync: one Mol* transaction). */
  batch<T>(fn: () => T): T {
    this.#depth++;
    try {
      return fn();
    } finally {
      if (--this.#depth === 0) this.#tell();
    }
  }

  #tell() {
    if (this.#pending.size === 0) return;
    const nodes = this.#pending;
    this.#pending = new Set();
    this.#state.set({ version: this.#state.value.version + 1, nodes });
  }
}

/** The same children (the same nodes, in order): a sync that makes the array again changes nothing. */
export function sameChildren<T>(a: readonly T[], b: readonly T[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}
