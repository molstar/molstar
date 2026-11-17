/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The spec's model (_spec/README.md "SpecModel"): the spec's immutable root, its revision (the root's hash), a
// transient preview over it, and the static model of what is shown. It is the history's port: a checkpoint is the
// root (structural sharing makes it cheap), restoring sets it back and drops a preview.

import { derive, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { createResolver, type ResolveOptions, type StaticModel } from '@molstar/seqx-core/resolve';
import type { HistoryContext, HistoryPort } from '@molstar/ui/history/history-model';

/** A change of the root: the next root, or a function of the current one. */
export type SpecChange = Node | ((root: Node) => Node);

/** Why `apply` refused a change. */
export type SpecRefusal = 'stale-revision';

export interface ApplyOptions {
  /** The revision the change was made against: refused when it is not the current one (EC-81). */
  readonly against?: string;
  /** The history run it is part of (the context `history.run` gives its work); without it, a change outside it. */
  readonly context?: HistoryContext;
}

/** What the spec tells its history: whether a run is in progress, and that the root changed outside it. */
export interface SpecHistory {
  isActive(context: HistoryContext | undefined): boolean;
  invalidate(reason: string): void;
}

/**
 * The history's options for a run that applies a change (`history.run(label, origin, async () => spec.apply(…),
 * SpecRunOptions)`): no step when nothing changed, a failure when refused.
 */
export const SpecRunOptions = {
  changed: (v: boolean | SpecRefusal) => v === true,
  failed: (v: boolean | SpecRefusal) => v === 'stale-revision',
} as const;

export class SpecModel extends ReactiveModel implements HistoryPort<Node> {
  private readonly _root: State<Node>;
  private readonly _preview = new State<Node | undefined>(undefined);
  /** The committed root: what the history checkpoints and a save writes. */
  readonly root: ReadonlyState<Node>;
  /** The root's hash: equal roots, equal revisions. */
  readonly revision: ReadonlyState<string>;
  /** A transient overlay (a drag resizing a track, a theme tried on hover): shown, never saved or in history. */
  readonly preview: ReadonlyState<Node | undefined>;
  /** What is shown: the preview when there is one, else the root. */
  readonly shown: ReadonlyState<Node>;
  /** The static model of what is shown, by a resolver that keeps its entries and paths between specs. */
  readonly model: ReadonlyState<StaticModel>;

  private _history: SpecHistory | undefined;

  constructor(root: Node, options: ResolveOptions = {}) {
    super();
    this._root = new State(root);
    this.root = this._root.asReadonly;
    this.preview = this._preview.asReadonly;
    this.revision = derive(this._root, hashNode);
    this.shown = derive([this._root, this._preview], (r, p) => p ?? r);
    const resolver = createResolver(options);
    this.model = derive(this.shown, (s) => resolver.resolve(s));
  }

  /**
   * Sets the root: one step when run inside the history (with `SpecRunOptions`, passing the run's context). Returns
   * whether it changed (an equal tree is no change), or why it was refused. A preview in progress is dropped. A change
   * outside a run of the attached history makes its earlier steps unsafe: undo never restores over it.
   */
  apply(change: SpecChange, options: ApplyOptions = {}): boolean | SpecRefusal {
    if (options.against !== undefined && options.against !== this.revision.value) return 'stale-revision';
    const next = typeof change === 'function' ? change(this._root.value) : change;
    this._preview.set(undefined);
    if (next === this._root.value || hashNode(next) === hashNode(this._root.value)) return false;
    this._root.set(next);
    if (this._history && !this._history.isActive(options.context))
      this._history.invalidate('The spec was changed outside the history');
    return true;
  }

  /** Tells `history` of changes made outside its runs (`HistoryModel` is one); returns the detach. */
  attachHistory(history: SpecHistory): () => void {
    this._history = history;
    return () => {
      if (this._history === history) this._history = undefined;
    };
  }

  /** Shows a change over the root without committing it; a later call replaces it (it is made from the root). */
  setPreview(change: SpecChange): void {
    this._preview.set(typeof change === 'function' ? change(this._root.value) : change);
  }

  /** Drops the preview (Esc): the root is shown again, nothing recorded. */
  cancelPreview(): void {
    this._preview.set(undefined);
  }

  /** Ends the preview in its one change (run it inside the history for one step); false when there was none. */
  commitPreview(options: Pick<ApplyOptions, 'context'> = {}): boolean {
    const p = this._preview.value;
    if (!p) return false;
    return this.apply(p, options) === true;
  }

  // ---- HistoryPort ----

  checkpoint(): Node {
    return this._root.value;
  }

  async restore(checkpoint: Node): Promise<void> {
    this._preview.set(undefined);
    this._root.set(checkpoint);
  }
}
