/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';

/**
 * What the app keeps a pane (shell §2.1, app side): its stage's model, a plugin, whatever shared chrome acts on. Keyed
 * by pane id (the shell's: `content` or a view's ref); `active` is the active pane's entry, from any state of the
 * active pane's id (the shell's `activePane`). The shell never sees it: the layout says which pane, the app knows what
 * it holds. DOM-free.
 */
export class PaneMap<T> {
  readonly #entries = new State<ReadonlyMap<string, T>>(new Map());
  /** The active pane and its entry, or undefined when the app keeps nothing for it. */
  readonly active: ReadonlyState<{ readonly pane: string; readonly value: T } | undefined>;

  constructor(activePane: ReadonlyState<string>) {
    this.active = derive(
      [activePane, this.#entries],
      (pane, m) => (m.has(pane) ? { pane, value: m.get(pane)! } : undefined),
      (a, b) => a === b || (!!a && !!b && a.pane === b.pane && a.value === b.value),
    );
  }

  /** Keeps `value` for `pane`; returns the undo (the last set wins; the undo takes only its own back). */
  set(pane: string, value: T): () => void {
    const map = new Map(this.#entries.value);
    map.set(pane, value);
    this.#entries.set(map);
    return () => {
      if (this.#entries.value.get(pane) !== value) return;
      const next = new Map(this.#entries.value);
      next.delete(pane);
      this.#entries.set(next);
    };
  }

  get(pane: string): T | undefined {
    return this.#entries.value.get(pane);
  }
}
