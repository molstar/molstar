/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static model's molecules (_spec/README.md "Resolve", Molecules and units): a spec's molecules, sequences'
// modified indices, complexes' units with their copy groups and the `maxUnits` bound, links held on per-row units,
// compounds by code. The model's build calls it with its own placing, caching and refs.

import { categoryColumn, type Issue, issueOf, type Node, printAddress, type Value } from '@molstar/core/repr/tree';
import { inlineSeq, type RangesList, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { type PairsData, readPairs } from './held.js';
import { ResolveIssues } from './issues.js';
import type { AccessPlan, AnnotationEntry, LinkEnds, ResolveConfig, SystemEntry } from './model.js';
import type { CompoundEntry, Molecular } from './molecule-query.js';
import {
  type AtomsEntry,
  type ComplexEntry,
  complexEntry,
  type ModifiedIndex,
  type MoleculeEntry,
  modifiedIndex,
  moleculeEntry,
  type UnitEntry,
  unitEntry,
} from './molecules.js';

/** A node of the spec and its address. */
export interface Placed {
  readonly node: Node;
  readonly address: string;
}

/** An entry kept between resolves by its key. */
export type Cached = <T>(key: string, make: () => T) => T;

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** A child of a complex or a molecule (a unit, links), placed under it. */
export interface Child extends Placed {
  readonly parent: Placed;
}

/** A complex's units and links and a molecule's links, by parent, each addressed under it (and found by address). */
export function placeChildren(data: readonly Placed[], byAddress: Map<string, Placed>): Map<Placed, Child[]> {
  const out = new Map<Placed, Child[]>();
  for (const at of data) {
    if (at.node.kind !== 'complex' && at.node.kind !== 'molecule') continue;
    const list: Child[] = [];
    const ids = at.address.split('/').filter(Boolean);
    for (const c of at.node.children ?? []) {
      if (c.id === undefined || (c.kind !== 'unit' && c.kind !== 'pairs')) continue;
      const placed = { node: c, address: printAddress([...ids, c.id]), parent: at };
      list.push(placed);
      byAddress.set(placed.address, placed);
    }
    out.set(at, list);
  }
  return out;
}

/** Links: a complex's between its units (an end's unit a row's), a molecule's between its residues. */
export function buildLinks(
  children: ReadonlyMap<Placed, readonly Child[]>,
  molecular: Molecular,
  systems: ReadonlyMap<string, SystemEntry>,
  annotations: Map<string, AnnotationEntry>,
  cached: Cached,
  key: (...at: Placed[]) => string,
): void {
  for (const [parent, list] of children) {
    const units = molecular.complexes.get(parent.address)?.units;
    if (parent.node.kind === 'complex' ? !units?.length : !systems.has(parent.address)) continue;
    for (const at of list) {
      if (at.node.kind !== 'pairs') continue;
      const entry = cached(`${at.address}|${key(at, parent)}`, () =>
        linksEntry(at.address, at.node, parent.address, units, systems),
      );
      annotations.set(at.address, entry);
    }
  }
}

/**
 * Molecules (their residues and atoms systems), sequences' modified indices, and complexes' units (systems on their
 * molecules' positions) with their copy groups; a complex past `maxUnits` units is refused (_spec/README.md "Resolve",
 * Molecules and units).
 */
export function buildMolecular(
  data: readonly Placed[],
  children: ReadonlyMap<Placed, readonly Child[]>,
  systems: Map<string, SystemEntry>,
  issues: Issue[],
  config: ResolveConfig,
  cached: Cached,
  key: (...at: Placed[]) => string,
  ref: (at: Placed, key: string, kinds: readonly string[]) => Placed | undefined,
  lengths: ReadonlyMap<string, number>,
): Molecular {
  const molecules = new Map<string, MoleculeEntry>();
  const atoms = new Map<string, AtomsEntry>();
  const modified = new Map<string, ModifiedIndex>();
  const nodes = new Map<string, Node>();
  for (const at of data) {
    if (at.node.kind === 'sequence') {
      const m = cached(`${at.address}|${key(at)}|modified`, () => modifiedIndex(at.node));
      if (m) modified.set(at.address, m);
    }
    if (at.node.kind !== 'molecule') continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    const e = cached(`${at.address}|${key(at)}`, () =>
      moleculeEntry(at.address, at.node, src?.address, lengths.get(at.address)),
    );
    molecules.set(at.address, e.molecule);
    nodes.set(at.address, at.node);
    for (const sys of e.systems) systems.set(sys.address, sys);
    if (e.molecule.atoms) atoms.set(e.molecule.atoms.address, e.molecule.atoms);
  }
  const units = new Map<string, UnitEntry>();
  const complexes = new Map<string, ComplexEntry>();
  for (const c of data) {
    if (c.node.kind !== 'complex') continue;
    const own = children.get(c) ?? [];
    const unitNodes = own.filter((u) => u.node.kind === 'unit');
    if (unitNodes.length > config.maxUnits) {
      issues.push(unitsIssue(c.address, unitNodes.length, config.maxUnits));
      continue;
    }
    const list: UnitEntry[] = [];
    for (const u of unitNodes) {
      const of = ref(u, 'of', ['sequence', 'molecule']);
      const system = of && systems.get(of.address);
      if (!of || !system) continue;
      const e = cached(`${u.address}|${key(u, of)}`, () =>
        unitEntry(u.address, c.address, u.node, { system, molecule: molecules.get(of.address) }),
      );
      list.push(e.unit);
      units.set(u.address, e.unit);
      nodes.set(u.address, u.node);
      for (const sys of e.systems) systems.set(sys.address, sys);
      if (e.atoms) atoms.set(e.atoms.address, e.atoms);
    }
    const links = own.filter((l) => l.node.kind === 'pairs').map((l) => l.address);
    complexes.set(c.address, complexEntry(c.address, c.node, list, links));
  }
  return { molecules, complexes, units, atoms, modified, nodes };
}

/** A complex past `maxUnits` units, refused. */
export function unitsIssue(address: string, count: number, max: number): Issue {
  return issueOf(ResolveIssues, 'seq.resolve.complex.units', `The complex has ${count} units, past ${max}`, {
    address,
    hint: ResolveIssues.codes['seq.resolve.complex.units'].hint,
  });
}

/** A position or a span as a reader says it (from 1). */
const spanText = (s: number, e: number) => (e - s > 1 ? `${s + 1}–${e}` : `${s + 1}`);

/**
 * A link's ends on each system as an inline ranges plan (what a `ranges` layer over links reads, _spec/README.md
 * "Render", Layer types: ranges), built on first ask: a row an end, its id `<key>:a` or `<key>:b`, its label the other
 * end (`→ Glycan 1`), its category the link's kind. Its source address `<links>#<system>`.
 */
function linkEnds(
  address: string,
  d: PairsData,
  systemOf: (k: number, side: 0 | 1) => string,
  systems: ReadonlyMap<string, SystemEntry>,
): LinkEnds | undefined {
  const on = new Set<string>();
  for (let k = 0; k < d.n; k++)
    for (const side of [0, 1] as const) {
      const s = systemOf(k, side);
      if (systems.has(s)) on.add(s);
    }
  if (!on.size) return undefined;
  const plans = new Map<string, AccessPlan>();
  const name = (s: string) => systems.get(s)?.name ?? s;
  const build = (system: string): AccessPlan => {
    const t = { start: [] as number[], end: [] as number[], id: [] as string[], label: [] as string[] };
    const kind: (string | null)[] = [];
    for (let k = 0; k < d.n; k++)
      for (const side of [0, 1] as const) {
        if (systemOf(k, side) !== system) continue;
        const there = systemOf(k, side ? 0 : 1);
        t.start.push(side ? d.bStart[k] : d.aStart[k]);
        t.end.push(side ? d.bEnd[k] : d.aEnd[k]);
        t.id.push(`${d.keys[k]}:${side ? 'b' : 'a'}`);
        t.label.push(`→ ${name(there)} ${side ? spanText(d.aStart[k], d.aEnd[k]) : spanText(d.bStart[k], d.bEnd[k])}`);
        kind.push(d.kind[k] < 0 ? null : d.kinds[d.kind[k]]);
      }
    const s = systems.get(system)!;
    const table: RangesList = {
      form: 'list',
      start: Int32Array.from(t.start),
      end: Int32Array.from(t.end),
      id: categoryColumn(t.id),
      fields: [
        { name: 'label', column: categoryColumn(t.label) },
        { name: 'category', column: categoryColumn(kind) },
      ],
    };
    const value = seqToValue(inlineSeq({ name: system, length: s.length, cyclic: s.cyclic }, { ranges: table }));
    return {
      kind: 'inline',
      source: { address: `${address}#${system}`, format: 'seq', value },
      runtime: SeqRuntime,
      system,
      shape: 'ranges',
    };
  };
  return {
    systems: on,
    plan(system) {
      if (!on.has(system)) return undefined;
      let plan = plans.get(system);
      if (!plan) plans.set(system, (plan = build(system)));
      return plan;
    },
    built: () => plans.values(),
  };
}

/** Links held whole: a complex's on per-row units (`units`, in its order), a molecule's on its residues. */
export function linksEntry(
  address: string,
  node: Node,
  parent: string,
  units?: readonly UnitEntry[],
  systems?: ReadonlyMap<string, SystemEntry>,
): AnnotationEntry {
  const index = units && new Map(units.map((u, i) => [u.id, i]));
  const pairs = readPairs(node, index && ((id) => index.get(id) ?? -1));
  const target = units?.[0]?.address ?? parent;
  const systemOf = (k: number, side: 0 | 1) => {
    const u = side ? pairs.bUnit : pairs.aUnit;
    return units && u ? (units[u[k]]?.address ?? target) : target;
  };
  const ends = systems ? linkEnds(address, pairs, systemOf, systems) : undefined;
  const fields = ['type', ...(pairs.labels ? ['label'] : []), ...Object.keys(pairs.extra ?? {})];
  return {
    address,
    shape: 'pairs',
    target,
    ySystem: target,
    fields,
    pairs,
    ...(units ? { units: units.map((u) => u.address) } : {}),
    ...(ends ? { ends } : {}),
  };
}

export function compoundEntry(n: Node): CompoundEntry {
  const name = p(n, 'name');
  const cls = p(n, 'class');
  const parents = p(n, 'parents');
  return {
    id: n.id ?? '',
    ...(typeof name === 'string' ? { name } : {}),
    ...(typeof cls === 'string' ? { class: cls } : {}),
    ...(Array.isArray(parents) ? { parents: parents.map(String) } : {}),
  };
}
