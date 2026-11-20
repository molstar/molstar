/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IdPattern, type Issue, issueOf } from '@molstar/core/repr/tree';
import { StructureIssues } from './issues.js';
import type { EntityType, StructureDescription } from './structure-description.js';

// A structure description checked (_spec/README.md "Seq io", Checked and bounded): every table and every cell read as
// its type, bounded before anything is sized by it, and what is left keyed by the data's identifiers and ordered by
// them, so the nodes made from it depend on no row order. A host is code Seq* does not own: what does not read is said
// (a coded issue at the source) and left out, never thrown.

/** The bounds a description is read within (Configuration). */
export interface StructureLimits {
  /** Positions of one entity: 1,000,000. */
  readonly maxLength: number;
  /** Units of the complex: 10,000. */
  readonly maxUnits: number;
  /** Rows of one table: 10,000,000. */
  readonly maxRows: number;
}

export const DefaultStructureLimits: StructureLimits = { maxLength: 1_000_000, maxUnits: 10_000, maxRows: 10_000_000 };

export interface ReadEntity {
  readonly id: string;
  readonly type: EntityType;
  readonly polymerType?: string;
  readonly name?: string;
  readonly letters?: string;
  /** Its positions: its letters', else its components' (a molecule's run of positions from 0). */
  readonly length: number;
  /** By position; a position's alternatives in their rows' order (the first gives the letter). */
  readonly components: readonly { readonly position: number; readonly component: string }[];
  /** By namespace, then value. */
  readonly identifiers: readonly { readonly namespace: string; readonly value: string }[];
  /** By residue, then name, then element. */
  readonly atoms: readonly { readonly residue: number; readonly name: string; readonly element: string }[];
}

export interface ReadUnit {
  readonly id: string;
  readonly chain: string;
  readonly auth?: string;
  readonly operator?: string;
  readonly entity: ReadEntity;
  /** Its listed positions (one row each: a position listed twice is left out): observed, and its `auth` label. */
  readonly residues: ReadonlyMap<
    number,
    { readonly observed: boolean; readonly label: string | null; readonly row: number }
  >;
}

export interface ReadLink {
  readonly a: string;
  readonly aPosition: number;
  readonly b: string;
  readonly bPosition: number;
  readonly type: string | null;
  readonly aAtom: string | null;
  readonly bAtom: string | null;
}

export interface ReadCompound {
  readonly id: string;
  readonly name?: string;
  readonly type?: string;
  readonly parents?: readonly string[];
}

/** A description as read: entities, units, compounds by their ids; links by their ends. */
export interface ReadStructure {
  readonly entities: readonly ReadEntity[];
  readonly units: readonly ReadUnit[];
  readonly links: readonly ReadLink[];
  readonly compounds: readonly ReadCompound[];
  /** The residue table's rows (what a named value is a column of). */
  readonly residueRows: number;
  readonly issues: readonly Issue[];
}

type Code = keyof (typeof StructureIssues)['codes'];

/** Code units' order, the tie-break of a natural order: ids are unique, so every order is total. */
const byCode = (a: string, b: string) => (a < b ? -1 : a > b ? 1 : 0);
export const natural = (a: string, b: string) => a.localeCompare(b, 'en', { numeric: true }) || byCode(a, b);

const isColumn = (x: unknown): x is ArrayLike<unknown> =>
  Array.isArray(x) || (ArrayBuffer.isView(x) && !(x instanceof DataView));
const str = (x: unknown): x is string => typeof x === 'string';
const int = (x: unknown): x is number => typeof x === 'number' && Number.isInteger(x);
const optStr = (x: unknown) => x === undefined || x === null || str(x);
const asOpt = (x: unknown) => (str(x) ? x : undefined);
/** A chain's id, which a unit's id joins to its operator by `-`: no `-` in it, so the join is unambiguous. */
const ChainPattern = /^[\p{L}\p{N}_][\p{L}\p{N}_.+]*$/u;
const OperatorPattern = /^[\p{L}\p{N}_.+-]+$/u;
const OpenWord = /^[a-z][a-z0-9_-]*$/;
const Letters = /^[A-Za-z]*$/;
const EntityTypes = new Set<string>(['polymer', 'non-polymer', 'branched', 'water']);

/** Reads a description within its bounds: what reads, by its identifiers, and what was said. */
export function readStructure(
  d: StructureDescription,
  options: { readonly source: string; readonly water?: boolean; readonly reserved: readonly string[] },
  limits: StructureLimits,
): ReadStructure {
  const issues: Issue[] = [];
  const say = (code: Code, message: string) =>
    issues.push(issueOf(StructureIssues, code, message, { address: options.source }));
  const bad = new Map<string, number>();
  const badRow = (table: string) => bad.set(table, (bad.get(table) ?? 0) + 1);
  const bounded = new Map<string, number>();
  const pastBound = (what: string) => bounded.set(what, (bounded.get(what) ?? 0) + 1);

  /** A table's row count: its columns (required and given optional ones) columns of one length within the bound. */
  const rows = (name: string, table: unknown, required: readonly string[], optional: readonly string[] = []) => {
    if (table === undefined || table === null) {
      if (required.length && ['entities', 'components', 'chains', 'residues'].includes(name))
        say('seq.io.structure.columns', `The description has no ${name} table: read as empty`);
      return 0;
    }
    if (typeof table !== 'object') {
      say('seq.io.structure.columns', `The ${name} table is no table: left out`);
      return 0;
    }
    const t = table as Record<string, unknown>;
    let n = -1;
    for (const c of [...required, ...optional]) {
      const col = t[c];
      if (col === undefined && !required.includes(c)) continue;
      if (!isColumn(col) || (n >= 0 && col.length !== n)) {
        say(
          'seq.io.structure.columns',
          `The ${name} table’s column ${c} is missing, no column or of another length: the table left out`,
        );
        return 0;
      }
      n = col.length;
    }
    if (n > limits.maxRows) {
      say('seq.io.structure.bounds', `The ${name} table has ${n} rows, more than ${limits.maxRows}: left out`);
      return 0;
    }
    return Math.max(0, n);
  };
  const col = (table: unknown, c: string): ArrayLike<unknown> | undefined =>
    (table as Record<string, ArrayLike<unknown> | undefined> | undefined)?.[c];
  const D = (typeof d === 'object' && d !== null ? d : {}) as Partial<Record<keyof StructureDescription, unknown>>;
  if (D !== d) say('seq.io.structure.columns', 'The description is no record: read as empty');

  // entities: unique ids, none a reserved key; all rows of an id given twice left out
  const e = D.entities;
  const nE = rows('entities', e, ['id', 'type'], ['polymerType', 'name', 'letters']);
  const eId = col(e, 'id')!;
  const count = new Map<string, number>();
  for (let i = 0; i < nE; i++) if (str(eId[i])) count.set(eId[i] as string, (count.get(eId[i] as string) ?? 0) + 1);
  const ents: (MutableEntity | undefined)[] = new Array(nE);
  for (let i = 0; i < nE; i++) {
    const id = eId[i];
    const type = col(e, 'type')![i];
    const [pt, name, letters] = ['polymerType', 'name', 'letters'].map((c) => col(e, c)?.[i]);
    if (!str(id) || !str(type) || !EntityTypes.has(type) || !optStr(pt) || !optStr(name) || !optStr(letters)) {
      badRow('entities');
      continue;
    }
    if (!IdPattern.test(id) || options.reserved.includes(id) || count.get(id)! > 1) {
      say('seq.io.structure.id', `The entity ${JSON.stringify(id)} is no id, a reserved key or given twice: left out`);
      continue;
    }
    if (type === 'water' && !options.water) continue;
    if (str(letters) && !Letters.test(letters)) {
      badRow('entities');
      continue;
    }
    if (str(letters) && letters.length > limits.maxLength) {
      pastBound('entity lengths');
      continue;
    }
    ents[i] = {
      id,
      type: type as EntityType,
      ...(str(pt) ? { polymerType: pt } : {}),
      ...(str(name) ? { name } : {}),
      ...(str(letters) ? { letters: letters.toUpperCase() } : {}),
      length: 0,
      components: [],
      identifiers: [],
      atoms: [],
    };
  }
  const entityAt = (x: unknown) => (int(x) && x >= 0 && x < nE ? ents[x] : undefined);
  const named = (x: unknown) => int(x) && x >= 0 && x < nE;

  // components
  const c = D.components;
  const nC = rows('components', c, ['entity', 'position', 'component']);
  for (let i = 0; i < nC; i++) {
    const [en, p, code] = ['entity', 'position', 'component'].map((k) => col(c, k)![i]);
    if (!named(en) || !int(p) || p < 0 || !str(code) || !code) badRow('components');
    else if (p >= limits.maxLength) pastBound('entity lengths');
    else entityAt(en)?.components.push({ position: p, component: code });
  }

  // identifiers and atoms
  const idt = D.identifiers;
  const nI = rows('identifiers', idt, ['entity', 'namespace', 'value']);
  for (let i = 0; i < nI; i++) {
    const [en, ns, v] = ['entity', 'namespace', 'value'].map((k) => col(idt, k)![i]);
    if (!named(en) || !str(ns) || !OpenWord.test(ns) || !str(v)) badRow('identifiers');
    else entityAt(en)?.identifiers.push({ namespace: ns, value: v });
  }
  const at = D.atoms;
  const nA = rows('atoms', at, ['entity', 'residue', 'name', 'element']);
  let hydrogens = 0;
  for (let i = 0; i < nA; i++) {
    const [en, r, name, element] = ['entity', 'residue', 'name', 'element'].map((k) => col(at, k)![i]);
    if (!named(en) || !int(r) || r < 0 || !str(name) || !name || !str(element) || !element) badRow('atoms');
    else if (/^(H|D)$/i.test(element)) hydrogens++;
    else entityAt(en)?.atoms.push({ residue: r, name, element });
  }
  if (hydrogens) say('seq.io.structure.hydrogen', `${hydrogens} hydrogen atoms left out`);

  // each entity's length: its letters', else its components' (a molecule's run of positions from 0)
  for (const x of ents) {
    if (!x) continue;
    x.components.sort((a, b) => a.position - b.position);
    x.identifiers.sort((a, b) => byCode(a.namespace, b.namespace) || byCode(a.value, b.value));
    let length = x.letters?.length ?? (x.components.length ? x.components[x.components.length - 1].position + 1 : 0);
    if (x.type !== 'polymer') {
      length = 0;
      for (const k of x.components) if (k.position === length) length++;
    }
    x.length = length;
    const atoms = x.atoms.filter((a) => a.residue < length);
    if (atoms.length < x.atoms.length) bad.set('atoms', (bad.get('atoms') ?? 0) + x.atoms.length - atoms.length);
    atoms.sort((a, b) => a.residue - b.residue || byCode(a.name, b.name) || byCode(a.element, b.element));
    x.atoms = atoms;
  }

  // units: an id each, of its chain and operator; every row of an id given twice left out
  const ch = D.chains;
  const nU = rows('chains', ch, ['label_asym_id', 'entity'], ['auth_asym_id', 'operator']);
  const unitIds: (string | undefined)[] = new Array(nU);
  const idCount = new Map<string, number>();
  for (let i = 0; i < nU; i++) {
    const [label, en, auth, op] = ['label_asym_id', 'entity', 'auth_asym_id', 'operator'].map((k) => col(ch, k)?.[i]);
    if (!str(label) || !named(en) || !optStr(auth) || !optStr(op)) {
      badRow('chains');
      continue;
    }
    if (!entityAt(en)) continue;
    if (!ChainPattern.test(label) || (str(op) && !OperatorPattern.test(op))) {
      say(
        'seq.io.structure.id',
        `The chain ${JSON.stringify(label)} (operator ${JSON.stringify(op ?? '')}) makes no id: left out`,
      );
      continue;
    }
    const id = str(op) ? `${label}-${op}` : label;
    unitIds[i] = id;
    idCount.set(id, (idCount.get(id) ?? 0) + 1);
  }
  for (const [id, n] of idCount)
    if (n > 1) say('seq.io.structure.id', `The unit ${JSON.stringify(id)} is given ${n} times: left out`);
  const unitRows: number[] = [];
  for (let i = 0; i < nU; i++) if (unitIds[i] !== undefined && idCount.get(unitIds[i]!) === 1) unitRows.push(i);
  unitRows.sort((a, b) => natural(unitIds[a]!, unitIds[b]!));
  if (unitRows.length > limits.maxUnits) {
    say('seq.io.structure.bounds', `${unitRows.length} units, more than ${limits.maxUnits}: the rest left out`);
    unitRows.length = limits.maxUnits;
  }
  const unitOfRow = new Map<number, MutableUnit>();
  const units: MutableUnit[] = unitRows.map((i) => {
    const entity = entityAt(col(ch, 'entity')![i])!;
    const u: MutableUnit = {
      id: unitIds[i]!,
      chain: col(ch, 'label_asym_id')![i] as string,
      ...(str(col(ch, 'auth_asym_id')?.[i]) ? { auth: col(ch, 'auth_asym_id')![i] as string } : {}),
      ...(str(col(ch, 'operator')?.[i]) ? { operator: col(ch, 'operator')![i] as string } : {}),
      entity,
      residues: new Map(),
    };
    unitOfRow.set(i, u);
    return u;
  });

  // residues: a position listed twice in a unit is left out (no row wins by order)
  const r = D.residues;
  const nR = rows('residues', r, ['chain', 'position'], ['auth_seq_id', 'ins_code', 'observed']);
  const twice = new Set<string>();
  for (let i = 0; i < nR; i++) {
    const [chain, p, auth, ins, obs] = ['chain', 'position', 'auth_seq_id', 'ins_code', 'observed'].map(
      (k) => col(r, k)?.[i],
    );
    if (
      !int(chain) ||
      chain < 0 ||
      chain >= nU ||
      !int(p) ||
      p < 0 ||
      !optStr(ins) ||
      (obs !== undefined && typeof obs !== 'number')
    ) {
      badRow('residues');
      continue;
    }
    const u = unitOfRow.get(chain);
    if (!u) continue;
    if (p >= u.entity.length) {
      badRow('residues');
      continue;
    }
    if (u.residues.has(p)) twice.add(`${u.id}\u0000${p}`);
    const label = int(auth) ? `${auth}${asOpt(ins) ?? ''}` : null;
    u.residues.set(p, { observed: obs === undefined || obs !== 0, label, row: i });
  }
  for (const k of twice) {
    const [id, p] = k.split('\u0000');
    units.find((u) => u.id === id)!.residues.delete(Number(p));
    badRow('residues');
  }

  // links between kept units' positions, by their ends
  const l = D.links;
  const nL = rows('links', l, ['a_chain', 'a_position', 'b_chain', 'b_position'], ['type', 'a_atom', 'b_atom']);
  const links: ReadLink[] = [];
  for (let i = 0; i < nL; i++) {
    const [ac, ap, bc, bp, type, aa, ba] = [
      'a_chain',
      'a_position',
      'b_chain',
      'b_position',
      'type',
      'a_atom',
      'b_atom',
    ].map((k) => col(l, k)?.[i]);
    if (
      ![ac, ap, bc, bp].every(int) ||
      (ac as number) < 0 ||
      (ac as number) >= nU ||
      (bc as number) < 0 ||
      (bc as number) >= nU
    ) {
      badRow('links');
      continue;
    }
    if (!optStr(type) || !optStr(aa) || !optStr(ba)) {
      badRow('links');
      continue;
    }
    const a = unitOfRow.get(ac as number);
    const b = unitOfRow.get(bc as number);
    if (!a || !b) continue;
    if (
      (ap as number) < 0 ||
      (bp as number) < 0 ||
      (ap as number) >= a.entity.length ||
      (bp as number) >= b.entity.length
    ) {
      badRow('links');
      continue;
    }
    links.push({
      a: a.id,
      aPosition: ap as number,
      b: b.id,
      bPosition: bp as number,
      type: asOpt(type) ?? null,
      aAtom: asOpt(aa) ?? null,
      bAtom: asOpt(ba) ?? null,
    });
  }
  const unitOrder = new Map(units.map((u, i) => [u.id, i]));
  const nullable = (a: string | null, b: string | null) =>
    byCode(a ?? '', b ?? '') || (a === null ? 0 : 1) - (b === null ? 0 : 1);
  links.sort(
    (x, y) =>
      unitOrder.get(x.a)! - unitOrder.get(y.a)! ||
      x.aPosition - y.aPosition ||
      unitOrder.get(x.b)! - unitOrder.get(y.b)! ||
      x.bPosition - y.bPosition ||
      nullable(x.type, y.type) ||
      nullable(x.aAtom, y.aAtom) ||
      nullable(x.bAtom, y.bAtom),
  );

  // compounds: unique codes; every row of a code given twice left out
  const k = D.compounds;
  const nK = rows('compounds', k, ['id'], ['name', 'type', 'parents']);
  const codes = new Map<string, number>();
  for (let i = 0; i < nK; i++) {
    const id = col(k, 'id')![i];
    if (str(id)) codes.set(id, (codes.get(id) ?? 0) + 1);
  }
  const compounds: ReadCompound[] = [];
  for (let i = 0; i < nK; i++) {
    const [id, name, type, parents] = ['id', 'name', 'type', 'parents'].map((x) => col(k, x)?.[i]);
    const parentsOk = parents === undefined || parents === null || (Array.isArray(parents) && parents.every(str));
    if (!str(id) || !optStr(name) || !optStr(type) || !parentsOk) {
      badRow('compounds');
      continue;
    }
    if (!IdPattern.test(id) || codes.get(id)! > 1) {
      say('seq.io.structure.id', `The compound ${JSON.stringify(id)} is no id or given twice: left out`);
      continue;
    }
    compounds.push({
      id,
      ...(str(name) ? { name } : {}),
      ...(str(type) ? { type } : {}),
      ...(Array.isArray(parents) && parents.length ? { parents: parents as string[] } : {}),
    });
  }
  compounds.sort((a, b) => natural(a.id, b.id));

  for (const [table, n] of bad)
    say(
      'seq.io.structure.rows',
      `${n} ${table} rows do not read (a cell of another type, a row of no table, a position past its unit): left out`,
    );
  for (const [what, n] of bounded)
    say('seq.io.structure.bounds', `${n} rows past the bound on ${what} (${limits.maxLength}): left out`);
  return {
    entities: ents.filter((x): x is MutableEntity => !!x).sort((a, b) => natural(a.id, b.id)),
    units,
    links,
    compounds,
    residueRows: nR,
    issues,
  };
}

type MutableEntity = Omit<ReadEntity, 'components' | 'identifiers' | 'atoms' | 'length'> & {
  length: number;
  components: { position: number; component: string }[];
  identifiers: { namespace: string; value: string }[];
  atoms: { residue: number; name: string; element: string }[];
};

type MutableUnit = Omit<ReadUnit, 'residues'> & {
  residues: Map<number, { observed: boolean; label: string | null; row: number }>;
};
