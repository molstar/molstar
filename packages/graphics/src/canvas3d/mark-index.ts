/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Loci } from '@molstar/model/model/loci';
import { Bond, Structure, StructureElement } from '@molstar/model/model/structure';
import type { Representation } from '../repr/representation.js';

/**
 * The representations a loci marked on all of them may change: for structure loci (structure, element, bond) those of
 * its root (`Representation.getStructureRootHashes`) and those that may mark any; for other loci all. A hover then
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

  /** adds `repr` or takes its hashes again */
  set(repr: Representation.Any) {
    this.unindex(repr);
    this.all.add(repr);
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
  }

  clear() {
    this.all.clear();
    this.byRoot.clear();
    this.any.clear();
    this.hashesOf.clear();
  }

  /** the representations whose `mark` may change for `loci`, each once */
  forEach(loci: Loci, cb: (repr: Representation.Any) => void) {
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
