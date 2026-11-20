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
  /**
   * Optional, for a numbering that answers without a scan (`index`, `ordinal`, an offset): a range start's position,
   * the first whose label is at or after `label` in the order of labels, or a range end's, the last at or before it;
   * -1 when none is. Without it, the forms scan the positions (bounded by the caller's signal and deadline).
   */
  snap?(label: string, side: 'start' | 'end'): number;
  /** Optional, as `snap`: up to `count` labels nearest `label` in the order of labels, for an issue's hint. */
  nearest?(label: string, count: number): string[];
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
   * Optional: whether it declares a numbering of this id whose labels have not arrived (a source's, read by region):
   * a form naming it waits for them, never reports it unknown.
   */
  numberingPending?(id: string): boolean;
  /**
   * The molecule a unit is of, and where the unit starts on it: the molecule's numberings number the unit's positions
   * (`start` + index) where the unit has none of that id.
   */
  readonly parent?: { readonly system: QuerySystem; readonly start: number };
}

/**
 * One row of an annotation: an item's key, its half-open extent on its system, its fields. A pairs annotation's pair is
 * one item with a row per end, each on its own system (_spec/README.md "Queries", Pairs answer `items`).
 */
export interface QueryItem {
  readonly key: string | number;
  readonly start: number;
  readonly end: number;
  /** The system this row is on: the annotation's `target` when absent (a ranges row's always). */
  readonly system?: string;
  readonly fields?: Readonly<Record<string, Atom>>;
}

/** What an annotation answers: the items it holds over what was asked, and what it has not loaded. */
export interface ItemsRead {
  /** Rows sharing a key are one item (a wrapped row is two rows of one item). */
  readonly items: readonly QueryItem[];
  /** The regions not loaded yet, absent or empty when complete. */
  readonly pending?: Location;
}

/** Items placed on a system (a pairs annotation's on two: each row says its own, `QueryItem.system`). */
export interface QueryAnnotation {
  readonly address: string;
  /** The address of the system its items are on (a pairs annotation's `a`: a row on `b` says so). */
  readonly target: string;
  /** The fields its items have: `(item :<field> v)` names one of them or `id`, else it is an issue. */
  readonly fields: readonly string[];
  /** Its rows overlapping a location on their systems (every row when undefined), with what is pending. */
  items(over?: Location): ItemsRead;
}

export interface QueryModel {
  system(address: string): QuerySystem | undefined;
  annotation(address: string): QueryAnnotation | undefined;
  /** A dynamic node held under a node by its key (`(node :key T1)`), or undefined. */
  node?(owner: string, key: Atom): Node | undefined;
  /**
   * A chain of a structure by its id (`(chain :id A)`), or undefined: a system data makes, addressed by its form
   * (`/data/1tup/(chain :id A)`), so `system` answers that address too.
   */
  chain?(owner: string, id: Atom): Node | undefined;
  /**
   * An alignment's row by its key (its label) as a position of its rows system (`rows`, `<alignment>/(rows)`), or
   * undefined: what `(row :key k)` reads (rows are positions, EC-35).
   */
  row?(rows: string, key: Atom): number | undefined;
  /** Whether an alignment's rows are not held yet (read from a source): a row form then waits, never unknown. */
  rowsPending?(rows: string): boolean;
  /** A system's components (a sequence's, a molecule's residues, a unit's): what `(component …)` reads. */
  components?(system: string): QueryComponents | undefined;
  /** A component's class by its code (its `compound` node's, else the host's names entry), or undefined. */
  compoundClass?(code: string): string | undefined;
  /** The annotations whose items are on a system (a ranges annotation's target): `:annotations` reads them. */
  annotationsOn?(system: string): readonly string[];
  /** A complex's units in order, or undefined when the address is no complex: `chain` and `entity` read them. */
  units?(complex: string): readonly QueryUnit[] | undefined;
  /** A molecule by its address (a spec node or a dynamic one), or undefined. */
  molecule?(address: string): QueryMolecule | undefined;
  /** The identifiers a host's resolver gives as the same compound as one (`ccd:ATP`), never including it. */
  equivalents?(identifier: string): readonly string[];
  /** The atoms system of a molecule or a unit (its own address for an atoms system), or undefined. */
  atoms?(system: string): QueryAtoms | undefined;
}

/**
 * A system's chemical components (_spec/README.md "Queries", Forms: `component`): a sequence's sparse non-standard
 * rows beside its letters' standard components, a molecule's residues a code each.
 */
export interface QueryComponents {
  readonly length: number;
  /** Whether its positions have standard components (a sequence or a unit of one): what `:standard` reads. */
  readonly polymer: boolean;
  /** Whether its letters are not held (a source's): a form reading standard components waits over the system. */
  readonly pending?: boolean;
  /** A sequence's non-standard rows, sorted: position, component, the positions it covers. */
  readonly modified?: {
    readonly position: Int32Array;
    readonly component: readonly string[];
    readonly length: Int32Array;
  };
  /** The components at a position: its non-standard ones (alternatives in order), else its standard one or code. */
  at(position: number): readonly string[];
}

/** A unit of a complex, as the node forms read it. */
export interface QueryUnit {
  readonly address: string;
  readonly node: Node;
  readonly id: string;
  /** The sequence or molecule it instantiates. */
  readonly of: string;
  readonly chain?: string;
  readonly authChain?: string;
  readonly operator?: string;
}

/** A molecule, as `(molecule …)` reads it. */
export interface QueryMolecule {
  readonly address: string;
  readonly node: Node;
  /** `ligand` when its node names none. */
  readonly class: string;
  /** `namespace:value` each. */
  readonly identifiers: readonly string[];
}

/** An atoms system: its address and each atom's name. */
export interface QueryAtoms {
  readonly address: string;
  readonly length: number;
  readonly names: { get(i: number): string | null };
}
