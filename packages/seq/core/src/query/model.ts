/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the query forms read (_spec/README.md "Queries", The model they read): systems with their numberings,
// annotations with their items, dynamic nodes. An interface, so the forms run over any model that implements it (the
// static model, a test's fixture); a source read by region answers with what it has loaded and what is pending.

import type { Atom, Node } from '@molstar/core/repr/tree';
import type { Location } from './location.js';

/**
 * A numbering of a system's positions: a label per position (undefined for a hole), and the first position of a
 * label. A label is the numbering's own spelling (labels.ts); a repeated label names its first position.
 */
export interface Numbering {
  readonly id: string;
  /** How many positions it numbers: the system's length (a unit's numbering of its molecule: the molecule's). */
  readonly length: number;
  labelAt(index: number): string | undefined;
  /** The first position with this label, or -1. */
  indexOf(label: string): number;
}

/** A coordinate system: a sequence, a unit, an alignment's columns. */
export interface QuerySystem {
  /** Its address, printed canonically. */
  readonly address: string;
  readonly length: number;
  /** A cyclic system's range may wrap past its end (as two ranges). */
  readonly cyclic?: boolean;
  /** The numberings declared on the system itself, by id (never `index` or `ordinal`: those are built in). */
  numbering(id: string): Numbering | undefined;
  /** The ids of its own numberings, for an issue's hint. */
  numberingIds(): readonly string[];
  /**
   * The molecule a unit is of, and where the unit starts on it: the molecule's numberings number the unit's positions
   * (`start` + index) where the unit has none of that id.
   */
  readonly parent?: { readonly system: QuerySystem; readonly start: number };
}

/** One row of an annotation: an item's key, its half-open extent on the target, its fields. */
export interface QueryItem {
  readonly key: string | number;
  readonly start: number;
  readonly end: number;
  readonly fields?: Readonly<Record<string, Atom>>;
}

/** What an annotation answers: the items it holds over what was asked, and what it has not loaded. */
export interface ItemsRead {
  /** Rows sharing a key are one item (a wrapped row is two rows of one item). */
  readonly items: readonly QueryItem[];
  /** The regions not loaded yet, absent or empty when complete. */
  readonly pending?: Location;
}

/** Items placed on one system. */
export interface QueryAnnotation {
  readonly address: string;
  /** The address of the system its items are on. */
  readonly target: string;
  /** Its items overlapping a location on the target (every item when undefined), with what is pending. */
  items(over?: Location): ItemsRead;
}

export interface QueryModel {
  system(address: string): QuerySystem | undefined;
  annotation(address: string): QueryAnnotation | undefined;
  /** A dynamic node held under a node by its key (`(node :key T1)`), or undefined. */
  node?(owner: string, key: Atom): Node | undefined;
}
