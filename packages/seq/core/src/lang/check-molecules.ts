/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The language's checks of molecules, complexes, units and links (_spec/README.md, Language, Molecules and
// complexes): what the description cannot say. Each takes what it reads as arguments, so the spec's checks and a
// dynamic node's (the overlay) share them.

import { CategoryColumn, isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import type { Code } from './check.js';
import { OpenWord } from './molecule-language.js';

/** Takes one issue of the node checked, at a step path under its params. */
export type Say = (code: Code, message: string, ...tail: (string | number)[]) => void;

const param = (n: Node, key: string): Value | undefined => n.params?.[key];
const ints = (v: Value | undefined): ArrayLike<number> | undefined =>
  isColumn(v) && !(v instanceof CategoryColumn) ? (v as ArrayLike<number>) : undefined;
const Word = new RegExp(OpenWord);
/** Hydrogen and its isotopes: refused in a molecule's atoms (EC-64). */
const Hydrogens = new Set(['H', 'D', 'T']);

/** The positions of a sequence (its inline letters) or a molecule (its inline residues): undefined when not inline. */
export function inlineLength(n: Node | undefined): number | undefined {
  if (!n) return undefined;
  if (n.kind === 'sequence') {
    const l = param(n, 'letters');
    return typeof l === 'string' ? l.length : undefined;
  }
  if (n.kind === 'molecule') {
    const r = param(n, 'residues');
    return r instanceof CategoryColumn ? r.length : undefined;
  }
  return undefined;
}

/**
 * A unit's range on what it instantiates, `[start, end]`, or undefined when it has none (all of it). Its own positions
 * are `[0, end - start)`: its position i is the molecule's `start + i` (_spec/README.md, Queries, Numberings and
 * labels).
 */
export function unitRange(unit: Node): readonly [number, number] | undefined {
  const r = param(unit, 'range');
  return Array.isArray(r) && r.length === 2 ? [r[0] as number, r[1] as number] : undefined;
}

/**
 * A sequence's `modified` table: positions in order (rows at one position are alternatives), each within the inline
 * letters with its `length`, a length of at least 1 that stops before the next row's position.
 */
export function checkModified(sequence: Node, say: Say): void {
  const m = param(sequence, 'modified') as Record<string, Value> | undefined;
  const pos = ints(m?.position);
  if (!pos) return;
  const len = ints(m?.length);
  const total = inlineLength(sequence);
  for (let i = 1; i < pos.length; i++) {
    if (pos[i] < pos[i - 1]) {
      say('seq.lang.modified.order', 'The modified positions are not in order', 'modified', 'position', i);
      return;
    }
  }
  for (let i = 0; i < pos.length; i++) {
    const l = len ? len[i] : 1;
    if (pos[i] < 0 || (total !== undefined && pos[i] + Math.max(l, 1) > total)) {
      say('seq.lang.modified.bounds', 'A modified position is past the sequence', 'modified', 'position', i);
      return;
    }
    let next = i + 1;
    while (next < pos.length && pos[next] === pos[i]) next++;
    if (l < 1 || (next < pos.length && pos[i] + l > pos[next])) {
      say('seq.lang.modified.length', 'A component’s length runs into the next row', 'modified', 'length', i);
      return;
    }
  }
}

/** Identifiers' namespaces are open lowercase words (an identifier is written `namespace:value` in queries). */
function checkIdentifiers(n: Node, say: Say): void {
  const ns = (param(n, 'identifiers') as Record<string, Value> | undefined)?.namespace;
  if (!(ns instanceof CategoryColumn)) return;
  for (let i = 0; i < ns.length; i++) {
    const v = ns.get(i);
    if (v === null || !Word.test(v)) {
      say('seq.lang.molecule.identifier', 'A namespace is no lowercase word', 'identifiers', 'namespace', i);
      return;
    }
  }
}

/**
 * A molecule's own checks: one of `residues` and `source` (atoms beside inline residues), atoms whose residues are
 * its positions, non-decreasing, no hydrogen (EC-64), identifiers' namespaces open words.
 */
export function checkMolecule(n: Node, say: Say): void {
  const residues = param(n, 'residues');
  const atoms = param(n, 'atoms') as Record<string, Value> | undefined;
  if ((residues === undefined) === (param(n, 'source') === undefined) || (atoms && residues === undefined))
    say('seq.lang.molecule.content', 'A molecule has one of `residues` and `source`, its atoms beside residues');
  checkIdentifiers(n, say);
  if (!atoms) return;
  const count = inlineLength(n);
  const res = ints(atoms.residue);
  if (res) {
    for (let i = 0; i < res.length; i++) {
      if (res[i] < 0 || (count !== undefined && res[i] >= count) || (i > 0 && res[i] < res[i - 1])) {
        say(
          'seq.lang.molecule.atoms',
          'An atom’s residue is out of order or not of the molecule',
          'atoms',
          'residue',
          i,
        );
        break;
      }
    }
  }
  const el = atoms.element;
  if (el instanceof CategoryColumn) {
    for (let i = 0; i < el.length; i++) {
      if (Hydrogens.has((el.get(i) ?? '').toUpperCase())) {
        say('seq.lang.molecule.hydrogen', 'A hydrogen in a molecule’s atoms', 'atoms', 'element', i);
        break;
      }
    }
  }
}

/** A compound's identifiers' namespaces are open words. */
export function checkCompound(n: Node, say: Say): void {
  checkIdentifiers(n, say);
}

/**
 * A unit's checks, given the length of what it instantiates (undefined when not inline): its range never wraps
 * (EC-13) and lies within it; its observed runs, on its own positions, in order within it.
 */
export function checkUnit(unit: Node, length: number | undefined, say: Say): void {
  const range = unitRange(unit);
  if (range && (range[1] < range[0] || (length !== undefined && range[1] > length)))
    say('seq.lang.unit.range', 'The unit’s range wraps or runs past its molecule', 'range');
  const obs = param(unit, 'observed') as Record<string, Value> | undefined;
  const s = ints(obs?.start);
  const e = ints(obs?.end);
  if (!s || !e) return;
  const hi = unitLength(unit, length);
  for (let i = 0; i < Math.min(s.length, e.length); i++) {
    if (e[i] < s[i] || s[i] < 0 || (hi !== undefined && e[i] > hi) || (i > 0 && s[i] < e[i - 1])) {
      say('seq.lang.unit.observed', 'An observed run is out of order or outside the unit', 'observed', 'start', i);
      return;
    }
  }
}

/** What a link's checks read of its parent: a complex's unit by id (how many positions it has), or a molecule. */
export interface LinkParent {
  readonly kind: 'data' | 'complex' | 'molecule' | string;
  /** Under a complex: the unit of an id, its `length` undefined when not known; undefined when there is no such unit. */
  unit?(id: string): { readonly length: number | undefined } | undefined;
  /** Under a molecule: its residues' count, when inline. */
  readonly length?: number;
}

/**
 * The checks of a `pairs` node as links: under `data` it names `a`; under a complex or a molecule neither `a` nor `b`
 * (the parent is both ends); unit columns only under a complex, each end naming a unit of it and lying in its range;
 * under a molecule each end within its residues; occupancy 0 to 1.
 */
export function checkLinks(pairs: Node, parent: LinkParent, say: Say): void {
  const under = parent.kind === 'complex' || parent.kind === 'molecule';
  if (under ? param(pairs, 'a') !== undefined || param(pairs, 'b') !== undefined : param(pairs, 'a') === undefined)
    say('seq.lang.link.target', under ? 'Links name no `a` or `b`: their parent is both' : 'The pairs name no `a`');
  const rows = param(pairs, 'rows') as Record<string, Value> | undefined;
  if (!rows) return;
  const occ = rows.occupancy;
  if (isColumn(occ)) {
    const o = occ as ArrayLike<number>;
    for (let i = 0; i < o.length; i++) {
      if (o[i] < 0 || o[i] > 1) {
        say('seq.lang.link.occupancy', 'An occupancy is not between 0 and 1', 'rows', 'occupancy', i);
        break;
      }
    }
  }
  for (const side of ['a', 'b'] as const) {
    const units = rows[`${side}_unit`];
    if (parent.kind !== 'complex') {
      if (units !== undefined) say('seq.lang.link.unit', 'A unit column outside a complex', 'rows', `${side}_unit`);
    }
    const s = ints(rows[`${side}_start`]);
    const e = ints(rows[`${side}_end`]);
    if (!s || !e) continue;
    const n = Math.min(s.length, e.length);
    for (let i = 0; i < n; i++) {
      let hi = parent.kind === 'molecule' ? parent.length : undefined;
      if (parent.kind === 'complex') {
        const id = units instanceof CategoryColumn ? units.get(i) : null;
        const u = id === null ? undefined : parent.unit?.(id);
        if (!u) {
          say('seq.lang.link.unit', 'A link’s end names no unit of the complex', 'rows', `${side}_unit`, i);
          break;
        }
        hi = u.length;
      }
      if (hi !== undefined && e[i] > hi) {
        say('seq.lang.link.bounds', 'A link’s end lies outside its unit or molecule', 'rows', `${side}_start`, i);
        break;
      }
    }
  }
}

/** How many positions a unit has: its range's, else what it instantiates has (`length`, when known). */
export function unitLength(unit: Node, length: number | undefined): number | undefined {
  const r = unitRange(unit);
  return r ? Math.max(0, r[1] - r[0]) : length;
}

/**
 * A block's `form: atoms` needs atoms: its system (a molecule, or the molecule a unit instantiates, `system`) has
 * them; a sequence never does. A molecule from a source is not judged (its atoms are its source's).
 */
export function checkBlockForm(block: Node, system: Node | undefined, say: Say): void {
  if (param(block, 'form') !== 'atoms' || !system) return;
  const none =
    system.kind === 'sequence' ||
    (system.kind === 'molecule' && param(system, 'source') === undefined && param(system, 'atoms') === undefined);
  if (none) say('seq.lang.molecule.form', 'The block is drawn as atoms, and its system has none', 'form');
}
