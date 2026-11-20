/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Issue, issueOf, type Node } from '@molstar/core/repr/tree';
import { type SequenceAlphabet, standardComponent, standardLetter } from '@molstar/core/util/sequence-constants';
import type { DynamicNode, SystemInfo } from '@molstar/seqx-core/data';
import { StructureIssues } from './issues.js';
import type { StructureDescription } from './structure-description.js';
import {
  DefaultStructureLimits,
  type ReadCompound,
  type ReadEntity,
  type ReadLink,
  type ReadUnit,
  readStructure,
  type StructureLimits,
} from './structure-read.js';

// A structure description into the dynamic nodes molecules are read as (_spec/README.md "Seq io"): a sequence per
// polymer entity, a molecule per other entity, the compounds, and a complex of units with their `auth` numbering,
// observed runs and links. Read first (`readStructure`: checked, bounded, keyed and ordered by the data's identifiers),
// so every adapter of one structure gives equal nodes, and a row removed or moved changes nothing else.

/** How a description is read: as the mmCIF runtime reads a file (molecules.md, Reading). */
export interface StructureNodesOptions {
  /** The address of the source the nodes are under: what a unit's `of` names its entity through. */
  readonly source: string;
  /** Water entities and their units: left out unless true. */
  readonly water?: boolean;
  /** The complex's id and key: `complex`. */
  readonly complex?: string;
  /** The bounds it is read within, over `DefaultStructureLimits`. */
  readonly limits?: Partial<StructureLimits>;
}

export interface StructureNodes {
  readonly nodes: readonly DynamicNode[];
  /** The units' systems (their ids, their lengths): what a source's coverage lists for values by region. */
  readonly systems: readonly SystemInfo[];
  readonly issues: readonly Issue[];
  /**
   * A named value of the description over each unit's positions (NaN where none), by the unit's id: what a source
   * serves for a `values` node's `field`. Undefined when the description has no such column a residue row.
   */
  values(field: string): ReadonlyMap<string, Float32Array> | undefined;
}

/** The compounds' container's id and key (the language's). */
const Compounds = 'compounds';

/** A polymer type's alphabet: DNA before RNA (`polydeoxyribonucleotide` holds `ribonucleotide`), else protein. */
function alphabetOf(polymerType: string | undefined): SequenceAlphabet {
  const t = polymerType?.toLowerCase() ?? '';
  if (t.includes('deoxyribonucleotide')) return 'dna';
  if (t.includes('ribonucleotide')) return 'rna';
  return 'protein';
}

/**
 * The nodes of a structure description, its units' systems and what reading it said: what does not read, or is past a
 * bound, is said and left out, never thrown.
 */
export function structureNodes(d: StructureDescription, options: StructureNodesOptions): StructureNodes {
  const complexId = options.complex ?? 'complex';
  const limits = { ...DefaultStructureLimits, ...options.limits };
  const read = readStructure(
    d,
    { source: options.source, water: options.water, reserved: [Compounds, complexId] },
    limits,
  );
  const issues = [...read.issues];
  const parents = new Map(read.compounds.filter((c) => c.parents).map((c) => [c.id, c.parents!]));
  const formulas = new Map(read.compounds.filter((c) => c.formula).map((c) => [c.id, c.formula!]));
  let unknown = 0;
  const nodes: DynamicNode[] = read.entities.map((x) => {
    const made = entityNode(x, parents, formulas);
    unknown += made.unknown;
    return { key: x.id, node: made.node };
  });
  if (unknown)
    issues.push(
      issueOf(
        StructureIssues,
        'seq.io.structure.letters',
        `${unknown} positions’ components have no standard letter: X or N there`,
        {
          address: options.source,
        },
      ),
    );
  const compounds = compoundsNode(read.compounds);
  if (compounds) nodes.push({ key: Compounds, node: compounds });
  if (read.units.length) {
    const children: Node[] = read.units.map((u) => unitNode(u, options.source));
    const links = linksNode(read.links);
    if (links) children.push(links);
    const name = typeof d?.name === 'string' ? { params: { name: d.name } } : {};
    nodes.push({ key: complexId, node: { kind: 'complex', id: complexId, ...name, children } as Node });
  }
  return {
    nodes,
    systems: read.units.map((u) => ({ name: u.id, length: u.entity.length })),
    issues,
    values: (field) => unitValues(d, field, read.units, read.residueRows),
  };
}

/** An entity's node: a polymer a sequence (letters given, else from its components), any other a molecule. */
function entityNode(
  x: ReadEntity,
  parents: ReadonlyMap<string, readonly string[]>,
  formulas: ReadonlyMap<string, string>,
): { node: Node; unknown: number } {
  const comps = x.components;
  if (x.type === 'polymer') {
    const alphabet = alphabetOf(x.polymerType);
    let letters = x.letters;
    let unknown = 0;
    if (letters === undefined) {
      const out = new Array<string>(x.length).fill(alphabet === 'protein' ? 'X' : 'N');
      // the first row at a position gives its letter (its alternatives follow)
      for (let k = 0; k < comps.length; k++) {
        if (k > 0 && comps[k - 1].position === comps[k].position) continue;
        const letter = letterOf(alphabet, comps[k].component, parents);
        if (letter) out[comps[k].position] = letter;
        else unknown++;
      }
      letters = out.join('');
    }
    const params: Record<string, unknown> = { alphabet, letters };
    const modified = modifiedRows(comps, letters, alphabet);
    if (modified) params.modified = modified;
    return { node: { kind: 'sequence', id: x.id, params } as Node, unknown };
  }
  // a molecule: its residues the first component at each position
  const residues: string[] = [];
  for (const k of comps) if (k.position === residues.length) residues.push(k.component);
  const params: Record<string, unknown> = {
    class: x.type === 'non-polymer' ? (isIon(residues, formulas, x.atoms.length) ? 'ion' : 'ligand') : x.type,
    residues: categoryColumn(residues),
  };
  if (x.atoms.length)
    params.atoms = {
      residue: Int32Array.from(x.atoms, (a) => a.residue),
      name: categoryColumn(x.atoms.map((a) => a.name)),
      element: categoryColumn(x.atoms.map((a) => a.element)),
    };
  if (x.name !== undefined) params.name = x.name;
  const links = entityLinksNode(x.links);
  if (x.identifiers.length)
    params.identifiers = {
      namespace: categoryColumn(x.identifiers.map((i) => i.namespace)),
      value: x.identifiers.map((i) => i.value),
    };
  return { node: { kind: 'molecule', id: x.id, params, ...(links ? { children: [links] } : {}) } as Node, unknown: 0 };
}

/**
 * A non-polymer is an ion when its one component is one atom, by its compound's formula (`Zn 2`, `Cl -1`; not
 * `O4 P -3`): its chemistry, not what was modeled (a phosphate with its P alone modeled is a ligand). Without a
 * formula, the atoms described tell (one: an ion).
 */
function isIon(residues: readonly string[], formulas: ReadonlyMap<string, string>, atoms: number): boolean {
  if (residues.length !== 1) return false;
  const formula = formulas.get(residues[0]);
  if (formula === undefined) return atoms === 1;
  let n = 0;
  for (const t of formula.trim().split(/\s+/)) {
    const m = /^([A-Z][a-z]?)(\d*)$/.exec(t);
    if (m) n += m[2] ? Number(m[2]) : 1;
  }
  return n === 1;
}

/** A component's letter: its own standard one, else its parent's (one parent: MSE is `M`), else none. */
function letterOf(alphabet: SequenceAlphabet, component: string, parents: ReadonlyMap<string, readonly string[]>) {
  const own = standardLetter(alphabet, component);
  if (own) return own;
  const p = parents.get(component);
  return p?.length === 1 ? standardLetter(alphabet, p[0]) : undefined;
}

/** The positions whose component is not their letter's standard one, and every row of a position with alternatives. */
function modifiedRows(comps: ReadEntity['components'], letters: string, alphabet: SequenceAlphabet) {
  const position: number[] = [];
  const component: string[] = [];
  for (let k = 0; k < comps.length; k++) {
    const p = comps[k].position;
    if (p >= letters.length) continue;
    const several = (k > 0 && comps[k - 1].position === p) || (k + 1 < comps.length && comps[k + 1].position === p);
    if (!several && comps[k].component === standardComponent(alphabet, letters[p])) continue;
    position.push(p);
    component.push(comps[k].component);
  }
  return position.length ? { position: Int32Array.from(position), component: categoryColumn(component) } : undefined;
}

/** A unit: its chain and operator, its observed runs when some position was not, its `auth` labels when given. */
function unitNode(u: ReadUnit, source: string): Node {
  const params: Record<string, unknown> = { of: `${source}/(node :key ${u.entity.id})`, chain: u.chain };
  if (u.auth !== undefined) params.auth_chain = u.auth;
  if (u.operator !== undefined) params.operator = u.operator;
  const length = u.entity.length;
  let all = true;
  for (let p = 0; p < length && all; p++) all = !!u.residues.get(p)?.observed;
  if (!all) {
    const start: number[] = [];
    const end: number[] = [];
    for (let p = 0; p < length; p++) {
      if (!u.residues.get(p)?.observed) continue;
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
  let labeled = false;
  for (const r of u.residues.values()) labeled ||= r.label !== null;
  if (labeled) {
    const labels = Array.from({ length }, (_, p) => u.residues.get(p)?.label ?? null);
    node.children = [
      { kind: 'numbering', id: 'auth', params: { method: 'labels', scheme: 'auth', labels: categoryColumn(labels) } },
    ];
  }
  return node as Node;
}

/** The complex's links between its units' positions, a row a link; undefined with none. */
function linksNode(links: readonly ReadLink[]): Node | undefined {
  if (!links.length) return undefined;
  const rows: Record<string, unknown> = {
    a_unit: categoryColumn(links.map((l) => l.a)),
    a_start: Int32Array.from(links, (l) => l.aPosition),
    a_end: Int32Array.from(links, (l) => l.aPosition + 1),
    b_unit: categoryColumn(links.map((l) => l.b)),
    b_start: Int32Array.from(links, (l) => l.bPosition),
    b_end: Int32Array.from(links, (l) => l.bPosition + 1),
  };
  if (links.some((l) => l.type !== null)) rows.type = categoryColumn(links.map((l) => l.type));
  if (links.some((l) => l.aAtom !== null)) rows.a_atom = categoryColumn(links.map((l) => l.aAtom));
  if (links.some((l) => l.bAtom !== null)) rows.b_atom = categoryColumn(links.map((l) => l.bAtom));
  return { kind: 'pairs', id: 'links', params: { rows } } as Node;
}

/** A molecule's own links between its positions (a glycan's bonds), a row a link; undefined with none. */
function entityLinksNode(links: ReadEntity['links']): Node | undefined {
  if (!links.length) return undefined;
  const rows: Record<string, unknown> = {
    a_start: Int32Array.from(links, (l) => l.aPosition),
    a_end: Int32Array.from(links, (l) => l.aPosition + 1),
    b_start: Int32Array.from(links, (l) => l.bPosition),
    b_end: Int32Array.from(links, (l) => l.bPosition + 1),
  };
  if (links.some((l) => l.type !== null)) rows.type = categoryColumn(links.map((l) => l.type));
  if (links.some((l) => l.aAtom !== null)) rows.a_atom = categoryColumn(links.map((l) => l.aAtom));
  if (links.some((l) => l.bAtom !== null)) rows.b_atom = categoryColumn(links.map((l) => l.bAtom));
  if (links.some((l) => l.order !== null)) rows.order = categoryColumn(links.map((l) => l.order));
  return { kind: 'pairs', id: 'links', params: { rows } } as Node;
}

/** The compounds described, by code; undefined with none. */
function compoundsNode(compounds: readonly ReadCompound[]): Node | undefined {
  if (!compounds.length) return undefined;
  const children = compounds.map((c) => {
    const params: Record<string, unknown> = {};
    if (c.name !== undefined) params.name = c.name;
    if (c.parents) params.parents = [...c.parents];
    if (c.type !== undefined) params.class = c.type;
    return { kind: 'compound', id: c.id, ...(Object.keys(params).length ? { params } : {}) } as Node;
  });
  return { kind: 'compounds', id: Compounds, children } as Node;
}

/** A named value over each kept unit's positions, from its own residue rows only (NaN where none). */
function unitValues(
  d: StructureDescription,
  field: string,
  units: readonly ReadUnit[],
  residueRows: number,
): ReadonlyMap<string, Float32Array> | undefined {
  const values = d?.values;
  const column = typeof values === 'object' && values !== null ? (values as Record<string, unknown>)[field] : undefined;
  if (!(Array.isArray(column) || (ArrayBuffer.isView(column) && !(column instanceof DataView)))) return undefined;
  const col = column as ArrayLike<unknown>;
  if (col.length !== residueRows) return undefined;
  const out = new Map<string, Float32Array>();
  for (const u of units) {
    const v = new Float32Array(u.entity.length).fill(Number.NaN);
    for (const [p, r] of u.residues) {
      const x = col[r.row];
      if (typeof x === 'number') v[p] = x;
    }
    out.set(u.id, v);
  }
  return out;
}
