/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Expression } from '@molstar/query-language/language/expression';
import type { QueryLanguages } from './languages.js';

// Targets (_spec/targets.md §1). Ported from Stories Untold's `any-molql`. Kept structural (the plain data of an MVS
// selector), so a target is data anywhere it goes; the one Mol* type here is MolQL's `Expression`, which is JSON.

/** The keywords that are fixed MolScript expressions (`keywords.ts`). */
export type Keyword = 'all' | 'polymer' | 'protein' | 'nucleic' | 'branched' | 'ligand' | 'ion' | 'water' | 'coarse';

/** The keywords that are the person's: Mol*'s current selection and focus, resolved by the viewer (targets.md §1). */
export type PersonKeyword = 'selection' | 'focus';

/**
 * One component expression (all given fields must match; `beg_`/`end_` are ranges): MVS's, with the fields of the
 * model's schema items (`StructureElement.Schema`) that MVS does not list, such as insertion codes of a range.
 */
export interface ComponentExpression {
  label_entity_id?: string;
  label_asym_id?: string;
  auth_asym_id?: string;
  label_seq_id?: number;
  auth_seq_id?: number;
  pdbx_PDB_ins_code?: string;
  beg_label_seq_id?: number;
  end_label_seq_id?: number;
  beg_auth_seq_id?: number;
  end_auth_seq_id?: number;
  /** Insertion code of the first residue of the range, which makes it a range in the order of the file. */
  beg_auth_ins_code?: string;
  /** Insertion code of the last residue of the range. */
  end_auth_ins_code?: string;
  label_comp_id?: string;
  auth_comp_id?: string;
  label_atom_id?: string;
  auth_atom_id?: string;
  type_symbol?: string;
  atom_id?: number;
  atom_index?: number;
  residue_index?: number;
  instance_id?: string;
}

/** A MolScript expression as JSON: the one form in which a target is kept (targets.md §2). */
export interface MolQLSelector {
  molql: unknown;
}

/** What MVS calls a selector: a keyword, one component expression, a list of them, or `{ molql }`. */
export type SelectorT = Keyword | ComponentExpression | readonly ComponentExpression[] | MolQLSelector;

/** Which atoms, residues, chains or structures an operation acts on (targets.md §1). */
export type Target = string | SelectorT;

/** A reference to a structure: a state ref, or a label ("1TQN"). A string for now. */
export type StructureRef = string;

export interface Scoped {
  target: Target;
  /** Absent: the structure the person works with (targets.md §1). */
  structure?: StructureRef;
}

/** Text languages MolQL can be written in. `molscript` is MolScript text (MolQL's own S-expressions and its aliases). */
export type QueryLanguage = 'pymol' | 'vmd' | 'jmol' | 'molscript';

/** A failure is a value, not a throw: the part that did not read, and what to write instead. */
export interface Failure {
  ok: false;
  message: string;
  hint: string;
}

/**
 * A string target as it was read: what the text was, before it is made MolQL. The viewer reads `kind` to resolve
 * `person` and to say what a result was ("selection", "ligand", "A:120-130", ...).
 */
export type ParsedTarget =
  | { kind: 'keyword'; keyword: Keyword }
  | { kind: 'person'; which: PersonKeyword }
  | { kind: 'addresses'; items: ComponentExpression[] }
  | { kind: 'query'; language: QueryLanguage; text: string }
  | { kind: 'selector'; selector: Exclude<SelectorT, Keyword> };

export type ParseResult = { ok: true; parsed: ParsedTarget } | Failure;

/** The viewer's answer for `selection` and `focus`: their MolQL, or nothing when there is none. */
export type PersonResolver = (which: PersonKeyword) => Expression | undefined;

export interface CompileOptions {
  /** The language of query text (the first the viewer reads unless said: PyMOL in the default spec). */
  language?: QueryLanguage;
  /** The languages the viewer reads (its own set; the page-wide one when none is given). */
  languages?: QueryLanguages;
  /** Resolves `selection` and `focus`; without one they fail with a message. */
  resolvePerson?: PersonResolver;
}

export type CompiledTarget = { ok: true; molql: Expression; parsed: ParsedTarget } | Failure;

/** What a selection selected in one author chain. */
export interface SelectedChain {
  /** Author chain id. */
  chain: string;
  residues: number;
  atoms: number;
  /** Residues by author number as compact ranges: `['120-130', '135', '140A']` (an insertion code stays single). */
  ranges: string[];
  /** Ranges left out past the limit. */
  more: number;
}

/**
 * "What it selected" (targets.md §2): comes back with every resolved target, so results, the history and the agent's
 * account need no second evaluation.
 */
export interface Selected {
  atoms: number;
  /** Selected residues (a residue in two instances of an assembly counts twice). */
  residues: number;
  /** Author chain ids, in structure order. */
  chains: string[];
  /** Per author chain, in structure order. */
  byChain: SelectedChain[];
  /** The ranges of all chains as addresses: `['A:120-130', 'A:135', 'B:7']`. */
  addresses: string[];
  /** The instance ids of the assembly or crystal copies the selection is in, when there are several. */
  instances: string[];
  /** One line for a person or a model: "42 residues (350 atoms) in chains A and B". */
  line: string;
}
