/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MarkerAction } from '@molstar/core/util/marker-action';
import { isEveryLoci, type Loci } from '@molstar/model/model/loci';
import { Bond, Structure, StructureElement } from '@molstar/model/model/structure';
import type { Representation } from '../repr/representation.js';

/** the actions that only remove markers: on `EveryLoci` they change only representations that hold some */
function isRemoving(action: MarkerAction) {
  return action === MarkerAction.RemoveHighlight || action === MarkerAction.Deselect || action === MarkerAction.Clear;
}

/** the actions that may add markers (`Toggle` both adds and removes) */
function isAdding(action: MarkerAction) {
  return action === MarkerAction.Highlight || action === MarkerAction.Select || action === MarkerAction.Toggle;
}

/**
 * The representations a loci marked on all of them may change: for structure loci (structure, element, bond) those of
 * its root (`Representation.getStructureRootHashes`) and those that may mark any; for `EveryLoci` with an action that
 * only removes markers (deselect all) those that markers were added to (`note`); for other loci all. A hover then
 * visits the representations of one structure, not those of every structure in the scene. Kept by `set` when a
 * representation is added or updated (its data, so its root, may have changed)
 */
export class MarkIndex {
  private readonly all = new Set<Representation.Any>();
  /** the representations by the root hashes they mark */
  private readonly byRoot = new Map<number, Set<Representation.Any>>();
  /** the representations that may mark loci of any root (no hashes) */
  private readonly any = new Set<Representation.Any>();
  private readonly hashesOf = new Map<Representation.Any, ReadonlyArray<number> | undefined>();
  /** the representations of `all` markers were added to since the last `Clear` on `EveryLoci` (a superset) */
  private readonly marked = new Set<Representation.Any>();
  /** those noted while not in `all` (removed, or not added yet): they keep their markers, `marked` again with `set` */
  private markedAway = new WeakSet<Representation.Any>();

  /** adds `repr` or takes its hashes again */
  set(repr: Representation.Any) {
    this.unindex(repr);
    this.all.add(repr);
    if (this.markedAway.has(repr)) {
      this.markedAway.delete(repr);
      this.marked.add(repr);
    }
    const hashes = repr.getStructureRootHashes?.();
    this.hashesOf.set(repr, hashes);
    if (!hashes) {
      this.any.add(repr);
      return;
    }
    for (const h of hashes) {
      let set = this.byRoot.get(h);
      if (!set) this.byRoot.set(h, (set = new Set()));
      set.add(repr);
    }
  }

  delete(repr: Representation.Any) {
    this.unindex(repr);
    this.all.delete(repr);
    if (this.marked.delete(repr)) this.markedAway.add(repr);
  }

  clear() {
    for (const repr of this.marked) this.markedAway.add(repr);
    this.marked.clear();
    this.all.clear();
    this.byRoot.clear();
    this.any.clear();
    this.hashesOf.clear();
  }

  /** `repr.mark(loci, action)` changed `repr`: noted when the action may have added markers */
  note(repr: Representation.Any, action: MarkerAction) {
    if (!isAdding(action)) return;
    if (this.all.has(repr)) this.marked.add(repr);
    else this.markedAway.add(repr);
  }

  /** the representations whose `mark` may change for `loci` and `action`, each once */
  forEach(loci: Loci, action: MarkerAction, cb: (repr: Representation.Any) => void) {
    if (isEveryLoci(loci) && isRemoving(action)) {
      // the others hold no markers: nothing to remove
      const marked = action === MarkerAction.Clear ? [...this.marked] : this.marked;
      if (action === MarkerAction.Clear) this.marked.clear();
      marked.forEach((r) => cb(r));
      return;
    }
    if (!(Structure.isLoci(loci) || StructureElement.Loci.is(loci) || Bond.isLoci(loci))) {
      this.all.forEach(cb);
      return;
    }
    this.byRoot.get(loci.structure.root.hashCode)?.forEach(cb);
    this.any.forEach(cb);
  }

  private unindex(repr: Representation.Any) {
    if (!this.hashesOf.has(repr)) return;
    const hashes = this.hashesOf.get(repr);
    this.hashesOf.delete(repr);
    if (!hashes) {
      this.any.delete(repr);
      return;
    }
    for (const h of hashes) {
      const set = this.byRoot.get(h);
      if (!set) continue;
      set.delete(repr);
      if (set.size === 0) this.byRoot.delete(h);
    }
  }
}
