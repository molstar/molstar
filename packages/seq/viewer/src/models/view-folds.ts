/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { ResolvedSection } from '@molstar/seqx-core/view';

// The view's folded track groups (_spec/README.md "ViewModel", Groups): each section's folded groups by id, session
// state as a band's clade folds are (the clade folds' model, one model of folding), the groups' `collapsed` the opening
// set; one map replaced on each change, a section's set compared by identity (the history and the frame read it so).

/** Each section's folded groups (by group id), the sections the session set; absent: each group's `collapsed`. */
export type GroupFolds = ReadonlyMap<string, readonly string[]>;

export const NoGroupFolds: GroupFolds = new Map();

const NoIds: readonly string[] = [];

/** A group's id: its address's last segment. */
export const groupIdOf = (address: string): string => address.slice(address.lastIndexOf('/') + 1);

/** The groups a section opens folded (its groups' `collapsed`), by id; one array a resolved section. */
const Opening = new WeakMap<ResolvedSection, readonly string[]>();
export function openingFolds(section: ResolvedSection): readonly string[] {
  let out = Opening.get(section);
  if (out) return out;
  const ids = (section.trackGroups ?? []).filter((g) => g.collapsed).map((g) => groupIdOf(g.address));
  out = ids.length ? ids : NoIds;
  Opening.set(section, out);
  return out;
}

export class ViewGroupFolds {
  private readonly _folds = new State<GroupFolds>(NoGroupFolds);
  readonly folds: ReadonlyState<GroupFolds> = this._folds.asReadonly;

  /** A section's own set (the session's), undefined when it has none (its groups' `collapsed` show). */
  get(section: string): readonly string[] | undefined {
    return this._folds.value.get(section);
  }

  get value(): GroupFolds {
    return this._folds.value;
  }

  /** Sets a section's own set (undefined: back to its groups' `collapsed`); false when it already was. */
  set(section: string, ids: readonly string[] | undefined): boolean {
    const prev = this._folds.value.get(section);
    if (prev === ids) return false;
    const map = new Map(this._folds.value);
    if (ids) map.set(section, ids);
    else map.delete(section);
    this._folds.set(map.size ? map : NoGroupFolds);
    return true;
  }

  /** Replaces every section's set (a step of the history restored). */
  replace(map: GroupFolds): void {
    if (map === this._folds.value) return;
    this._folds.set(map.size ? map : NoGroupFolds);
  }
}

/** Two sets of folds the same: by identity a section (a toggle makes a new array). */
export function foldsEqual(a: GroupFolds | undefined, b: GroupFolds | undefined): boolean {
  const x = a ?? NoGroupFolds;
  const y = b ?? NoGroupFolds;
  if (x === y) return true;
  if (x.size !== y.size) return false;
  for (const [k, v] of x) if (y.get(k) !== v) return false;
  return true;
}

/**
 * The folded groups a reveal opens (_spec/README.md "ViewModel", Groups): the groups of the view's sections whose shown
 * tracks alone draw what the data owners make (an item's key: a layer's data, or a form over it); none when a track
 * not folded draws it (the reveal goes there) or nothing does.
 */
export function foldedHolding(
  sections: readonly ResolvedSection[],
  folded: (section: string) => readonly string[],
  owners: ReadonlySet<string>,
): string[] {
  const out: string[] = [];
  for (const s of sections) {
    const ids = folded(s.address);
    for (const t of s.tracks) {
      if (!t.layers.some((l) => l.data && owners.has(ownerOf(l.data.address)))) continue;
      const g = t.group;
      // an authored summary track stands for its folded group: drawn while folded
      if (!g || t.summaryOf || !ids.includes(groupIdOf(g))) return [];
      if (!out.includes(g)) out.push(g);
    }
  }
  return out;
}

/** A data address's node part (`/data/x/(item …)` → `/data/x`). */
const ownerOf = (address: string): string => {
  const k = address.indexOf('/(');
  return k < 0 ? address : address.slice(0, k);
};
