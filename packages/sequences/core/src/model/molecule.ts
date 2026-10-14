/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Molecules: sequences and ligands (seqspec:model#molecules, the one owner of the kind rules).

import {
  CustomData,
  Issue,
  type MoleculeId,
  type NameId,
  Provenance,
  type ProvenanceInput,
  type StringTable,
} from "./common.js";
import { Numbering, type NumberingInput, type NumberingSet } from "./numbering.js";

export const MoleculeKind = { polymer: 0, ligand: 1 } as const;
export type Alphabet = "protein" | "dna" | "rna" | "generic";

export interface Molecule {
  readonly id: MoleculeId;
  /** Its ref ("p53", "1brs/entity:1"). */
  readonly key: NameId;
  /** MoleculeKind. */
  readonly kind: number;
  /** A polymer's; null for a ligand. */
  readonly alphabet: Alphabet | null;
  /** Residues, or heavy atoms. */
  readonly length: number;
  /** Upper-case ASCII letters; a ligand: `X` per atom; lazy: empty. */
  readonly code: Uint8Array;
  /** Letters served by a streamed source (seqspec:sources#streamed-sources): a length and an alphabet, no letters. */
  readonly lazy: boolean;
  /** Its coordinate system is cyclic (seqspec:model#cyclic-systems): positions modulo `length`; a polymer's only. */
  readonly cyclic: boolean;
  /** A ligand: the atomic number per atom. */
  readonly element: Uint8Array | null;
  /** NameId per position: a chemical component code (MSE; a ligand's own). */
  readonly chemCompId: Int32Array | null;
  readonly label: NameId;
  readonly identifiers: readonly { namespace: NameId; value: NameId; version: NameId }[];
  /** Molecule-level (seqspec:model#numberings), `ordinal` first. */
  readonly numberings: NumberingSet;
  readonly provenance: Provenance | null;
  readonly custom: CustomData;
}

export interface MoleculeInput {
  key: string;
  /** Polymer by default. */
  kind?: "polymer" | "ligand";
  alphabet?: Alphabet;
  residues?: string;
  /** Letters streamed: no `residues`, an alphabet required. */
  lazy?: { length: number };
  /** A polymer's coordinate system is cyclic (a ligand's never is). */
  cyclic?: boolean;
  elements?: string[];
  chemCompIds?: (string | null)[];
  label?: string;
  numberings?: NumberingInput[];
  ownUnit?: boolean;
  identifiers?: { namespace: string; value: string; version?: string }[];
  provenance?: ProvenanceInput;
  custom?: object;
  /**
   * Runtime properties (seqspec:model#custom-data), on the built molecule's `custom.props` (a lazy one's streamed
   * source).
   */
  props?: Record<string, unknown>;
}

const AMBIGUITY = "RYSWKMBDHV";
const LETTERS: Record<Alphabet, string> = {
  protein: "ACDEFGHIKLMNPQRSTVWYBZJUOX*",
  dna: `ACGTN${AMBIGUITY}`,
  rna: `ACGUN${AMBIGUITY}`,
  generic: "",
};

// 128-entry tables, built once: residues are ASCII, so a character is a table index, never a string.
/** Upper case by ASCII rules only: a to z, nothing else. */
const UPPER = Uint8Array.from({ length: 128 }, (_, c) => (c >= 97 && c <= 122 ? c - 32 : c));
const flags = (has: (c: number) => boolean) => Uint8Array.from({ length: 128 }, (_, c) => (has(c) ? 1 : 0));
const inLetters = (letters: string) => flags((c) => letters.includes(String.fromCharCode(c)));
/** Generic: printable ASCII but space and the gap characters `-` and `.`. */
const VALID: Record<Alphabet, Uint8Array> = {
  protein: inLetters(LETTERS.protein),
  dna: inLetters(LETTERS.dna),
  rna: inLetters(LETTERS.rna),
  generic: flags((c) => c > 32 && c < 127 && c !== 45 && c !== 46),
};
/** Detection's classes: 1 for `ACGTUN`, 2 for an ambiguity code (a nucleotide code, but not `ACGTUN`), 0 else. */
const CLASS = Uint8Array.from({ length: 128 }, (_, c) => {
  const ch = String.fromCharCode(c);
  return "ACGTUN".includes(ch) ? 1 : AMBIGUITY.includes(ch) ? 2 : 0;
});
/**
 * How many times each code occurs in the residues being checked: one scratch, filled by the pass that reads them
 * (`build`'s, or `count`) and read at once, so detection and validation cost a 256-entry loop, not a pass each.
 * The residue pass is the only one a valid sequence takes; the index of a bad letter is looked up on failure only.
 * Integers, so that a count is a small integer before the loop is optimized too (a Float64Array boxes each one).
 */
const COUNTS = new Uint32Array(256);

function count(code: Uint8Array): void {
  COUNTS.fill(0);
  for (let i = 0; i < code.length; i++) COUNTS[code[i]]++;
}

/** The first index whose code is not flagged in `table` (a 128-entry table: codes above 127 never are), or -1. */
function firstNotIn(code: Uint8Array, table: Uint8Array): number {
  for (let i = 0; i < code.length; i++) if (code[i] > 127 || table[code[i]] !== 1) return i;
  return -1;
}
const NUCLEOTIDE = Uint8Array.from(CLASS, (k) => (k > 0 ? 1 : 0));

/** Where detection found the letter that rules a nucleotide alphabet out, or -1. */
export interface Detection {
  at: number;
}

/** `detectAlphabet` over COUNTS, already filled for `code`. */
function detectCounted(code: Uint8Array, out: Detection | undefined): Alphabet {
  if (out) out.at = -1;
  let plain = 0;
  let other = 0;
  for (let c = 0; c < 256; c++) {
    const k = c < 128 ? CLASS[c] : 0;
    if (k === 1) plain += COUNTS[c];
    else if (k === 0) other += COUNTS[c];
  }
  const n = code.length;
  if (n === 0 || plain * 10 < n * 9) return "protein";
  if (other > 0) {
    if (out) out.at = firstNotIn(code, NUCLEOTIDE);
    return "protein";
  }
  return COUNTS[85] > 0 && COUNTS[84] === 0 ? "rna" : "dna";
}

/**
 * Nucleotide when every letter is a nucleotide code and at least 90% of all letters are `ACGTUN` (an ambiguity code
 * counts as a letter that is not); RNA when `U` is there and `T` is not. Else protein, and when at least 90% are
 * `ACGTUN` but a letter is no nucleotide code (`AAAAAAAAAGE`), `out.at` is that letter's index, for the warning.
 */
function detectAlphabet(code: Uint8Array, out?: Detection): Alphabet {
  count(code);
  return detectCounted(code, out);
}

/** The first letter of `code` outside the alphabet, or -1, from COUNTS (already filled for `code`). */
function invalidCounted(alphabet: Alphabet, code: Uint8Array): number {
  const valid = VALID[alphabet];
  for (let c = 0; c < 256; c++) if (COUNTS[c] > 0 && (c > 127 || valid[c] !== 1)) return firstNotIn(code, valid);
  return -1;
}

/** `at`: detection's warning (-1 when the alphabet is declared); `invalid`: the first letter outside it, or -1. */
export interface ResidueCheck extends Detection {
  invalid: number;
}

/**
 * The alphabet of upper-case residue codes, `declared` or detected, and in `out` where detection warns and the first
 * letter outside it. One counting pass over `code`: an alignment's residues take this.
 */
function checkResidues(code: Uint8Array, declared: Alphabet | undefined, out: ResidueCheck): Alphabet {
  count(code);
  return checkCounted(code, declared, out);
}

function checkCounted(code: Uint8Array, declared: Alphabet | undefined, out: ResidueCheck): Alphabet {
  out.at = -1;
  const alphabet = declared ?? detectCounted(code, out);
  out.invalid = invalidCounted(alphabet, code);
  return alphabet;
}

function isLetter(alphabet: Alphabet, c: number): boolean {
  return c >= 0 && c < 128 && VALID[alphabet][c] === 1;
}

// biome-ignore format: the periodic table, one line per row
const SYMBOLS = (
  "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr " +
  "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu " +
  "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr " +
  "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split(" ");
const ATOMIC_NUMBER = new Map(SYMBOLS.map((s, i) => [s.toUpperCase(), i + 1]));

/** The atomic number of an element symbol, case-insensitive; 0 when unknown. */
function atomicNumber(symbol: string): number {
  return ATOMIC_NUMBER.get(symbol.toUpperCase()) ?? 0;
}
/** The symbol of an atomic number ("Fe"); "" when unknown. */
function elementSymbol(z: number): string {
  return SYMBOLS[z - 1] ?? "";
}

/** A molecule from its input, or null with errors. */
function build(input: MoleculeInput, id: number, strings: StringTable, issues: Issue[], path: string): Molecule | null {
  const err = (code: string, field: string, message: string) => {
    issues.push(Issue.error(code, `${path}/${field}`, `Molecule "${input.key}": ${message}`));
    return null;
  };
  const ligand = input.kind === "ligand";
  let code: Uint8Array;
  let element: Uint8Array | null = null;
  let alphabet: Alphabet | null = null;
  if (ligand) {
    if (!input.elements || input.residues != null)
      return err("model.molecule.content", "elements", "a ligand has elements only.");
    element = new Uint8Array(input.elements.length);
    for (let i = 0; i < element.length; i++) {
      element[i] = atomicNumber(input.elements[i]);
      if (element[i] <= 1) {
        const what = element[i] ? "hydrogens are left out of a ligand" : "an unknown element symbol";
        return err("model.ligand.element", `elements/${i}`, `"${input.elements[i]}" is ${what}.`);
      }
    }
    code = new Uint8Array(element.length).fill(88);
  } else if (input.lazy) {
    if (input.residues != null || input.elements || !input.alphabet || !(input.lazy.length > 0))
      return err("model.molecule.content", "lazy", "a lazy sequence has a positive length and an alphabet only.");
    code = new Uint8Array(0);
    alphabet = input.alphabet;
  } else {
    if (input.residues == null || input.elements)
      return err("model.molecule.content", "residues", "a sequence has residues only.");
    const text = input.residues;
    if (!text.length) return err("model.molecule.content", "residues", "a sequence has residues, and this has none.");
    code = new Uint8Array(text.length);
    // The one pass a valid sequence takes: ASCII, upper case, and the counts detection and validation read.
    COUNTS.fill(0);
    for (let i = 0; i < text.length; i++) {
      const c = text.charCodeAt(i);
      if (c > 127) {
        // A surrogate pair is named by its code point (U+1F600), at the index of its first half.
        const point = (text.codePointAt(i) ?? c).toString(16).toUpperCase();
        return err("model.residues.invalid", "residues", `U+${point} at ${i} is not ASCII.`);
      }
      const u = UPPER[c];
      code[i] = u;
      COUNTS[u]++;
    }
    const check: ResidueCheck = { at: -1, invalid: -1 };
    alphabet = checkCounted(code, input.alphabet, check);
    if (check.at >= 0) {
      const m = `"${String.fromCharCode(code[check.at])}" at ${check.at} is no nucleotide code, so it is read as protein.`;
      issues.push(Issue.warning("model.alphabet.detected", `${path}/residues`, `Molecule "${input.key}": ${m}`));
    }
    const bad = check.invalid;
    if (bad >= 0) {
      const m = `"${String.fromCharCode(code[bad])}" at ${bad} is not a ${alphabet} letter.`;
      return err("model.residues.invalid", "residues", m);
    }
  }
  const length = input.lazy && !ligand ? input.lazy.length : code.length;
  let chemCompId: Int32Array | null = null;
  const ccd = input.chemCompIds;
  if (ccd) {
    if (ccd.length !== length && !(ligand && ccd.length === 1)) {
      return err("model.molecule.length", "chemCompIds", `${ccd.length} chemical components for ${length} positions.`);
    }
    chemCompId = Int32Array.from({ length }, (_, i) => {
      const c = ccd.length === 1 ? ccd[0] : ccd[i];
      return c == null ? -1 : strings.intern(c);
    });
  }
  return {
    id: id as MoleculeId,
    key: strings.intern(input.key),
    kind: ligand ? MoleculeKind.ligand : MoleculeKind.polymer,
    alphabet,
    length,
    code,
    lazy: !!input.lazy && !ligand,
    cyclic: input.cyclic === true && !ligand,
    element,
    chemCompId,
    label: strings.intern(input.label ?? input.key),
    identifiers: (input.identifiers ?? []).map((x) => ({
      namespace: strings.intern(x.namespace),
      value: strings.intern(x.value),
      version: (x.version == null ? -1 : strings.intern(x.version)) as NameId,
    })),
    numberings: Numbering.buildSet(input.numberings, length, strings, issues, `${path}/numberings`, true),
    provenance: input.provenance ? Provenance.create(input.provenance, strings) : null,
    custom: CustomData.create(input.custom, input.props),
  };
}

export const Molecules = { build, detectAlphabet, checkResidues, isLetter, atomicNumber, elementSymbol };
