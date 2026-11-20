/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Issue, IdPattern, issueOf, type Node } from '@molstar/core/repr/tree';
import { type SequenceAlphabet, standardComponent, standardLetter } from '@molstar/core/util/sequence-constants';
import type { DynamicNode, SystemInfo } from '@molstar/seqx-core/data';
import { StructureIssues } from './issues.js';
import type { StructureDescription } from './structure-description.js';

// A structure description into the dynamic nodes molecules are read as (_spec/README.md "Seq io"): a sequence per
// polymer entity, a molecule per other entity, the compounds, and a complex of units with their `auth` numbering,
// observed runs and links. Keys come from the data's identifiers, never from row indices: an entity by its id, a unit
// by its chain and operator, so every adapter of one structure gives equal keys, and a row removed or moved renames
// nothing else.

/** How a description is read: as the mmCIF runtime reads a file (molecules.md, Reading). */
export interface StructureNodesOptions {
  /** The address of the source the nodes are under: what a unit's `of` names its entity through. */
  readonly source: string;
  /** Water entities and their units: left out unless true. */
  readonly water?: boolean;
  /** The complex's id and key: `complex`. */
  readonly complex?: string;
}

export interface StructureNodes {
  readonly nodes: readonly DynamicNode[];
  /** The units' systems (their ids, their lengths): what a source's coverage lists for values by region. */
  readonly systems: readonly SystemInfo[];
  readonly issues: readonly Issue[];
}

/** The compounds' container's id and key (the language's). */
const Compounds = 'compounds';

/** A unit's id: its chain, and its operator after a `-` when it has one; characters no id takes as `_`. */
export function unitId(chain: string, operator: string | undefined): string {
  const id = operator ? `${chain}-${operator}` : chain;
  return id.replace(/[^\p{L}\p{N}_.+-]/gu, '_').replace(/^[.+-]/, '_');
}

/** The address of a source's dynamic node by key: what a unit's `of` names. */
export function nodeRef(source: string, key: string): string {
  return `${source}/(node :key ${key})`;
}

/** A polymer type's alphabet: DNA before RNA (`polydeoxyribonucleotide` holds `ribonucleotide`), else protein. */
function alphabetOf(polymerType: string | undefined): SequenceAlphabet {
  const t = polymerType?.toLowerCase() ?? '';
  if (t.includes('deoxyribonucleotide')) return 'dna';
  if (t.includes('ribonucleotide')) return 'rna';
  return 'protein';
}

const natural = (a: string, b: string) => a.localeCompare(b, 'en', { numeric: true });

interface Entity {
  readonly row: number;
  readonly id: string;
  readonly length: number;
  readonly node: Node;
}

interface Unit {
  readonly id: string;
  readonly chain: number;
  readonly entity: Entity;
}

/**
 * The nodes of a structure description, its units' systems and what reading it said. Tables whose columns differ in
 * length are left out; rows naming no row of another table, or a position past their entity, are left out; each said
 * once a table.
 */
export function structureNodes(d: StructureDescription, options: StructureNodesOptions): StructureNodes {
  const issues: Issue[] = [];
  const say = (code: keyof (typeof StructureIssues)['codes'], message: string) =>
    issues.push(issueOf(StructureIssues, code, message, { address: options.source }));
  const complexId = options.complex ?? 'complex';
  const ok = (table: string, columns: Record<string, ArrayLike<unknown> | undefined>) => {
    const lengths = new Set(
      Object.values(columns)
        .filter((c) => c !== undefined)
        .map((c) => c!.length),
    );
    if (lengths.size <= 1) return true;
    say('seq.io.structure.columns', `The ${table} table’s columns differ in length: left out`);
    return false;
  };
  const e = d.entities;
  const entityRows = ok('entities', { ...e }) ? e.id.length : 0;
  const c = d.chains;
  const chainRows = ok('chains', { ...c }) ? c.label_asym_id.length : 0;
  const comp = d.components;
  const compRows = ok('components', { ...comp }) ? comp.entity.length : 0;
  const r = d.residues;
  const residueRows = ok('residues', { ...r }) ? r.chain.length : 0;

  // the components of each entity, by position (alternatives in their order)
  const byEntity = new Map<number, number[]>();
  let badRows = 0;
  for (let i = 0; i < compRows; i++) {
    const en = comp.entity[i];
    if (en < 0 || en >= entityRows || comp.position[i] < 0) badRows++;
    else {
      let list = byEntity.get(en);
      if (!list) byEntity.set(en, (list = []));
      list.push(i);
    }
  }
  if (badRows) say('seq.io.structure.rows', `${badRows} component rows name no entity or no position: left out`);
  for (const list of byEntity.values()) list.sort((a, b) => comp.position[a] - comp.position[b] || a - b);

  // entities, by id
  const entities = new Map<number, Entity>();
  const taken = new Set([Compounds, complexId]);
  const order = Array.from({ length: entityRows }, (_, i) => i).sort((a, b) => natural(e.id[a], e.id[b]) || a - b);
  let unknownLetters = 0;
  const parents = parentsOf(d);
  for (const i of order) {
    const type = e.type[i];
    if (type === 'water' && !options.water) continue;
    const id = e.id[i];
    if (!IdPattern.test(id) || taken.has(id)) {
      say('seq.io.structure.id', `The entity ${JSON.stringify(id)} is no id, or given twice: left out`);
      continue;
    }
    taken.add(id);
    const rows = byEntity.get(i) ?? [];
    const name = e.name?.[i];
    if (type === 'polymer') {
      const alphabet = alphabetOf(e.polymerType?.[i]);
      let letters = e.letters?.[i]?.toUpperCase();
      if (letters === undefined) {
        const length = rows.length ? comp.position[rows[rows.length - 1]] + 1 : 0;
        const unknown = alphabet === 'protein' ? 'X' : 'N';
        const out = new Array<string>(length).fill(unknown);
        // the first row at a position gives its letter (its alternatives follow)
        for (let k = 0; k < rows.length; k++) {
          if (k > 0 && comp.position[rows[k - 1]] === comp.position[rows[k]]) continue;
          const letter = letterOf(alphabet, comp.component[rows[k]], parents);
          if (letter) out[comp.position[rows[k]]] = letter;
          else unknownLetters++;
        }
        letters = out.join('');
      }
      const modified = modifiedRows(rows, comp, letters, alphabet);
      const params: Record<string, unknown> = { alphabet, letters };
      if (modified) params.modified = modified;
      entities.set(i, { row: i, id, length: letters.length, node: { kind: 'sequence', id, params } as Node });
    } else {
      // a molecule: its residues the first component at each position
      const residues: string[] = [];
      for (const k of rows) if (comp.position[k] === residues.length) residues.push(comp.component[k]);
      const atoms = moleculeAtoms(d, i, residues.length, say);
      const params: Record<string, unknown> = {
        class: type === 'non-polymer' ? (atoms?.name.length === 1 ? 'ion' : 'ligand') : type,
        residues: categoryColumn(residues),
      };
      if (atoms)
        params.atoms = {
          residue: Int32Array.from(atoms.residue),
          name: categoryColumn(atoms.name),
          element: categoryColumn(atoms.element),
        };
      if (name !== undefined) params.name = name;
      const ids = identifiersOf(d, i);
      if (ids) params.identifiers = ids;
      entities.set(i, { row: i, id, length: residues.length, node: { kind: 'molecule', id, params } as Node });
    }
  }
  if (unknownLetters)
    say('seq.io.structure.letters', `${unknownLetters} positions’ components have no standard letter: X or N there`);

  // units, in the chains' order
  const units: Unit[] = [];
  const unitOf = new Map<number, Unit>();
  const unitIds = new Set<string>();
  let badChains = 0;
  for (let i = 0; i < chainRows; i++) {
    const en = c.entity[i];
    if (en < 0 || en >= entityRows) {
      badChains++;
      continue;
    }
    const entity = entities.get(en);
    if (!entity) continue;
    const id = unitId(c.label_asym_id[i], c.operator?.[i]);
    if (!IdPattern.test(id) || unitIds.has(id)) {
      say('seq.io.structure.id', `The unit ${JSON.stringify(id)} is no id, or given twice: left out`);
      continue;
    }
    unitIds.add(id);
    const u = { id, chain: i, entity };
    units.push(u);
    unitOf.set(i, u);
  }
  if (badChains) say('seq.io.structure.rows', `${badChains} chain rows name no entity: left out`);

  // each unit's residues: observed, auth labels
  const observed = new Map<Unit, Uint8Array>();
  const labels = new Map<Unit, (string | null)[]>();
  let badResidues = 0;
  for (let i = 0; i < residueRows; i++) {
    const ch = r.chain[i];
    const u = ch >= 0 && ch < chainRows ? unitOf.get(ch) : undefined;
    const p = r.position[i];
    if (ch < 0 || ch >= chainRows || p < 0 || (u && p >= u.entity.length)) {
      badResidues++;
      continue;
    }
    if (!u) continue;
    let obs = observed.get(u);
    if (!obs) observed.set(u, (obs = new Uint8Array(u.entity.length)));
    obs[p] = r.observed ? r.observed[i] : 1;
    if (r.auth_seq_id) {
      let l = labels.get(u);
      if (!l) labels.set(u, (l = new Array<string | null>(u.entity.length).fill(null)));
      l[p] = `${r.auth_seq_id[i]}${r.ins_code?.[i] ?? ''}`;
    }
  }
  if (badResidues)
    say('seq.io.structure.rows', `${badResidues} residue rows name no chain or a position past it: left out`);

  const children: Node[] = units.map((u) => unitNode(u, d, options.source, observed.get(u), labels.get(u)));
  const links = linksNode(d, unitOf, chainRows, say);
  if (links) children.push(links);

  const nodes: DynamicNode[] = [...entities.values()].map((x) => ({ key: x.id, node: x.node }));
  const compounds = compoundsNode(d, say);
  if (compounds) nodes.push({ key: Compounds, node: compounds });
  if (units.length) {
    const complex: Node = {
      kind: 'complex',
      id: complexId,
      ...(d.name !== undefined ? { params: { name: d.name } } : {}),
      children,
    };
    nodes.push({ key: complexId, node: complex });
  }
  const systems = units.map((u) => ({ name: u.id, length: u.entity.length }));
  return { nodes, systems, issues };
}

/** Each described compound's parents (`MSE`: `MET`). */
function parentsOf(d: StructureDescription): ReadonlyMap<string, readonly string[]> {
  const out = new Map<string, readonly string[]>();
  const t = d.compounds;
  for (let i = 0; i < (t?.id.length ?? 0); i++) {
    const p = t!.parents?.[i];
    if (p?.length && !out.has(t!.id[i])) out.set(t!.id[i], p);
  }
  return out;
}

/** A component's letter: its own standard one, else its parent's (one parent: MSE is `M`), else none. */
function letterOf(alphabet: SequenceAlphabet, component: string, parents: ReadonlyMap<string, readonly string[]>) {
  const own = standardLetter(alphabet, component);
  if (own) return own;
  const p = parents.get(component);
  return p?.length === 1 ? standardLetter(alphabet, p[0]) : undefined;
}

/** The positions whose component is not their letter's standard one, and every row of a position with alternatives. */
function modifiedRows(
  rows: readonly number[],
  comp: StructureDescription['components'],
  letters: string,
  alphabet: SequenceAlphabet,
) {
  const position: number[] = [];
  const component: string[] = [];
  for (let k = 0; k < rows.length; k++) {
    const p = comp.position[rows[k]];
    if (p >= letters.length) continue;
    const several =
      (k > 0 && comp.position[rows[k - 1]] === p) || (k + 1 < rows.length && comp.position[rows[k + 1]] === p);
    if (!several && comp.component[rows[k]] === standardComponent(alphabet, letters[p])) continue;
    position.push(p);
    component.push(comp.component[rows[k]]);
  }
  return position.length ? { position: Int32Array.from(position), component: categoryColumn(component) } : undefined;
}

/** A molecule entity's heavy atoms, by residue (hydrogens left out, said); undefined with none. */
function moleculeAtoms(
  d: StructureDescription,
  entity: number,
  residues: number,
  say: (code: 'seq.io.structure.hydrogen' | 'seq.io.structure.rows', message: string) => void,
) {
  const a = d.atoms;
  if (!a) return undefined;
  const rows: number[] = [];
  let hydrogens = 0;
  let bad = 0;
  for (let i = 0; i < a.entity.length; i++) {
    if (a.entity[i] !== entity) continue;
    if (a.residue[i] < 0 || a.residue[i] >= residues) bad++;
    else if (/^(H|D)$/i.test(a.element[i])) hydrogens++;
    else rows.push(i);
  }
  if (hydrogens) say('seq.io.structure.hydrogen', `${hydrogens} hydrogen atoms left out`);
  if (bad) say('seq.io.structure.rows', `${bad} atom rows name a residue past their molecule: left out`);
  if (!rows.length) return undefined;
  rows.sort((x, y) => a.residue[x] - a.residue[y] || x - y);
  return {
    residue: rows.map((i) => a.residue[i]),
    name: rows.map((i) => a.name[i]),
    element: rows.map((i) => a.element[i]),
  };
}

function identifiersOf(d: StructureDescription, entity: number) {
  const t = d.identifiers;
  if (!t) return undefined;
  const namespace: string[] = [];
  const value: string[] = [];
  for (let i = 0; i < t.entity.length; i++)
    if (t.entity[i] === entity) {
      namespace.push(t.namespace[i]);
      value.push(t.value[i]);
    }
  return namespace.length ? { namespace: categoryColumn(namespace), value } : undefined;
}

/** A unit: its chain and operator, its observed runs when some position was not, its `auth` labels when given. */
function unitNode(
  u: Unit,
  d: StructureDescription,
  source: string,
  observed: Uint8Array | undefined,
  labels: (string | null)[] | undefined,
): Node {
  const c = d.chains;
  const params: Record<string, unknown> = { of: nodeRef(source, u.entity.id), chain: c.label_asym_id[u.chain] };
  const auth = c.auth_asym_id?.[u.chain];
  if (auth !== undefined) params.auth_chain = auth;
  const operator = c.operator?.[u.chain];
  if (operator !== undefined) params.operator = operator;
  const obs = observed ?? new Uint8Array(u.entity.length);
  if (obs.some((x) => x === 0)) {
    const start: number[] = [];
    const end: number[] = [];
    for (let p = 0; p < obs.length; p++) {
      if (!obs[p]) continue;
      if (end.length && end[end.length - 1] === p) end[end.length - 1] = p + 1;
      else {
        start.push(p);
        end.push(p + 1);
      }
    }
    params.observed = { start: Int32Array.from(start), end: Int32Array.from(end) };
  }
  const node: { kind: string; id: string; params: Record<string, unknown>; children?: Node[] } = {
    kind: 'unit',
    id: u.id,
    params,
  };
  if (labels)
    node.children = [
      { kind: 'numbering', id: 'auth', params: { method: 'labels', scheme: 'auth', labels: categoryColumn(labels) } },
    ];
  return node as Node;
}

/** The complex's links between its units' positions, a row a link; undefined with none. */
function linksNode(
  d: StructureDescription,
  unitOf: ReadonlyMap<number, Unit>,
  chainRows: number,
  say: (code: 'seq.io.structure.rows', message: string) => void,
): Node | undefined {
  const l = d.links;
  if (!l) return undefined;
  const cols = { a_unit: [] as string[], a_start: [] as number[], b_unit: [] as string[], b_start: [] as number[] };
  const type: (string | null)[] = [];
  const aAtom: (string | null)[] = [];
  const bAtom: (string | null)[] = [];
  let bad = 0;
  for (let i = 0; i < l.a_chain.length; i++) {
    const [ac, bc, ap, bp] = [l.a_chain[i], l.b_chain[i], l.a_position[i], l.b_position[i]];
    const a = unitOf.get(ac);
    const b = unitOf.get(bc);
    const named = ac >= 0 && ac < chainRows && bc >= 0 && bc < chainRows;
    // a link to a unit left out (water) goes quietly; one naming no chain or a position past its unit is said
    if (!named || (a && b && (ap < 0 || bp < 0 || ap >= a.entity.length || bp >= b.entity.length))) {
      bad++;
      continue;
    }
    if (!a || !b) continue;
    cols.a_unit.push(a.id);
    cols.a_start.push(ap);
    cols.b_unit.push(b.id);
    cols.b_start.push(bp);
    type.push(l.type?.[i] ?? null);
    aAtom.push(l.a_atom?.[i] ?? null);
    bAtom.push(l.b_atom?.[i] ?? null);
  }
  if (bad) say('seq.io.structure.rows', `${bad} link rows name no chain or a position past it: left out`);
  if (!cols.a_unit.length) return undefined;
  const rows: Record<string, unknown> = {
    a_unit: categoryColumn(cols.a_unit),
    a_start: Int32Array.from(cols.a_start),
    a_end: Int32Array.from(cols.a_start, (p) => p + 1),
    b_unit: categoryColumn(cols.b_unit),
    b_start: Int32Array.from(cols.b_start),
    b_end: Int32Array.from(cols.b_start, (p) => p + 1),
  };
  if (l.type) rows.type = categoryColumn(type);
  if (l.a_atom) rows.a_atom = categoryColumn(aAtom);
  if (l.b_atom) rows.b_atom = categoryColumn(bAtom);
  return { kind: 'pairs', id: 'links', params: { rows } } as Node;
}

/** The compounds described, by code (one a code, the first kept); undefined with none. */
function compoundsNode(
  d: StructureDescription,
  say: (code: 'seq.io.structure.id', message: string) => void,
): Node | undefined {
  const t = d.compounds;
  if (!t?.id.length) return undefined;
  const children: Node[] = [];
  const seen = new Set<string>();
  const order = Array.from({ length: t.id.length }, (_, i) => i).sort((a, b) => natural(t.id[a], t.id[b]) || a - b);
  for (const i of order) {
    const id = t.id[i];
    if (!IdPattern.test(id) || seen.has(id)) {
      say('seq.io.structure.id', `The compound ${JSON.stringify(id)} is no id, or given twice: left out`);
      continue;
    }
    seen.add(id);
    const params: Record<string, unknown> = {};
    if (t.name?.[i] !== undefined) params.name = t.name[i];
    if (t.parents?.[i]?.length) params.parents = [...t.parents[i]!];
    if (t.type?.[i] !== undefined) params.class = t.type[i];
    children.push({ kind: 'compound', id, ...(Object.keys(params).length ? { params } : {}) } as Node);
  }
  return children.length ? ({ kind: 'compounds', id: Compounds, children } as Node) : undefined;
}

/**
 * A named value of the description over each unit's positions (NaN where none): what a source serves for a `values`
 * node's `field` on a unit system. Undefined when the description has no such values.
 */
export function unitValues(
  d: StructureDescription,
  field: string,
  systems: readonly SystemInfo[],
  options: Pick<StructureNodesOptions, 'water'> = {},
): ReadonlyMap<string, Float32Array> | undefined {
  const column = d.values?.[field];
  if (!column || column.length !== d.residues.chain.length) return undefined;
  const lengths = new Map(systems.map((s) => [s.name, s.length ?? 0]));
  const out = new Map<string, Float32Array>();
  for (let i = 0; i < column.length; i++) {
    const ch = d.residues.chain[i];
    if (ch < 0 || ch >= d.chains.label_asym_id.length) continue;
    const entity = d.chains.entity[ch];
    if (d.entities.type[entity] === 'water' && !options.water) continue;
    const id = unitId(d.chains.label_asym_id[ch], d.chains.operator?.[ch]);
    const length = lengths.get(id);
    const p = d.residues.position[i];
    if (length === undefined || p < 0 || p >= length) continue;
    let v = out.get(id);
    if (!v) out.set(id, (v = new Float32Array(length).fill(Number.NaN)));
    v[p] = column[i];
  }
  return out;
}
