/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static model's molecules (_spec/README.md "Resolve", Molecules and units): a spec's molecules, sequences'
// modified indices, complexes' units with their copy groups and the `maxUnits` bound, links held on per-row units,
// compounds by code. The model's build calls it with its own placing, caching and refs.

import { type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import { readPairs } from './held.js';
import { ResolveIssues } from './issues.js';
import type { AnnotationEntry, ResolveConfig, SystemEntry } from './model.js';
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

/**
 * Molecules (their residues and atoms systems), sequences' modified indices, and complexes' units (systems on their
 * molecules' positions) with their copy groups; a complex past `maxUnits` units is refused (_spec/README.md "Resolve",
 * Molecules and units).
 */
export function buildMolecular(
  data: readonly Placed[],
  inner: readonly (Placed & { readonly parent: Placed })[],
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
    const children = inner.filter((u) => u.parent === c && u.node.kind === 'unit');
    if (children.length > config.maxUnits) {
      issues.push(unitsIssue(c.address, children.length, config.maxUnits));
      continue;
    }
    const list: UnitEntry[] = [];
    for (const u of children) {
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
    const links = inner.filter((l) => l.parent === c && l.node.kind === 'pairs').map((l) => l.address);
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

/** Links held whole: a complex's on per-row units (`units`, in its order), a molecule's on its residues. */
export function linksEntry(address: string, node: Node, parent: string, units?: readonly UnitEntry[]): AnnotationEntry {
  const index = units && new Map(units.map((u, i) => [u.id, i]));
  const pairs = readPairs(node, index && ((id) => index.get(id) ?? -1));
  const target = units?.[0]?.address ?? parent;
  const fields = ['type', ...(pairs.labels ? ['label'] : []), ...Object.keys(pairs.extra ?? {})];
  return {
    address,
    shape: 'pairs',
    target,
    ySystem: target,
    fields,
    pairs,
    ...(units ? { units: units.map((u) => u.address) } : {}),
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
