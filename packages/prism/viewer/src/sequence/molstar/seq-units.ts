/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A Seq* complex's units as a structure provider made them (seqspec:libraries#the-mol-adapter): their refs and `mvs`
// fields, and the chain instance of a Mol* structure each names.

import {
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  KeyKind,
  type SeqSet,
} from '@molstar/seq-core/model';

const UNIT = CoordinateSystemKind.unit;

import type { Structure } from '@molstar/model/model/structure';
import { type ChainInstance, chainsOf } from './chains.js';

export interface MvsFields {
  label_asym_id?: string;
  auth_asym_id?: string;
  instance_id?: string;
}

/** The complex's units (ids and refs) with their `mvs` fields, in its order. */
export function complexUnits(
  set: SeqSet,
  complex: string,
): { unit: number; ref: string; mvs: MvsFields; first: number }[] {
  const t = { kind: 0, id: 0, row: -1 };
  if (!set.byName(complex, t) || t.kind !== KeyKind.complex) return [];
  const cx = set.complexes[t.id];
  return Array.from(cx?.units ?? [], (u) => {
    const mol = set.molecules[set.units.molecule[u]];
    const label = mol?.numberings.byName(set.strings.find('label') as never);
    return {
      unit: u,
      ref: set.strings.get(set.units.key[u]),
      mvs: (set.units.custom[u]?.doc as { mvs?: MvsFields } | undefined)?.mvs ?? {},
      first: label?.method === 'offset' ? label.first : 1,
    };
  });
}

/** The chain instance a unit's `mvs` fields name: label_asym_id and instance, else the first with its label_asym_id. */
export function chainFor(structure: Structure, mvs: MvsFields): ChainInstance | null {
  const { byInstance, byAsym } = chainLookup(structure);
  return byInstance.get(`${mvs.label_asym_id}@${mvs.instance_id}`) ?? byAsym.get(mvs.label_asym_id ?? '') ?? null;
}

/** The chain instances by `label_asym_id@instance_id` and by the first of each `label_asym_id`, once per index. */
const lookups = new WeakMap<object, { byInstance: Map<string, ChainInstance>; byAsym: Map<string, ChainInstance> }>();
function chainLookup(structure: Structure) {
  const index = chainsOf(structure);
  let l = lookups.get(index);
  if (l) return l;
  l = { byInstance: new Map(), byAsym: new Map() };
  for (const c of index.chains) {
    const key = `${c.labelAsymId}@${c.instanceId}`;
    if (!l.byInstance.has(key)) l.byInstance.set(key, c);
    if (!l.byAsym.has(c.labelAsymId)) l.byAsym.set(c.labelAsymId, c);
  }
  lookups.set(index, l);
  return l;
}

/** A unit's own map onto its molecule: the same indices, over its range. */
export function unitMap(set: SeqSet, u: number): CoordinateMap {
  const start = set.units.rangeStart[u];
  const molecule = { kind: CoordinateSystemKind.molecule, id: set.units.molecule[u], row: -1 };
  return Coords.create({ kind: UNIT, id: u, row: -1 }, molecule, 1, [start, set.units.rangeStop[u] - start, start, 1]);
}

/**
 * The Seq* routes a request from outside the complex takes (a hover on P04637 reaching chain C through its SIFTS
 * mapping): each linked unit onto its molecule, every mapping of those molecules, and the units of the molecules at
 * their other ends. The link registers them with the structure's system (so `detach` drops them); the composer
reads them in a translator of its own.
 */
export function seqRoutes(set: SeqSet, units: ReadonlySet<number>): CoordinateMap[] {
  const molecules = new Set<number>();
  for (const u of units) molecules.add(set.units.molecule[u]);
  const maps: CoordinateMap[] = [];
  const far = new Set<number>();
  const molOf = (s: CoordinateSystem) =>
    s.kind === CoordinateSystemKind.molecule ? s.id : s.kind === UNIT ? set.units.molecule[s.id] : -1;
  for (const m of set.mappings) {
    if (!m) continue;
    const [a, b] = [molOf(m.from), molOf(m.to)];
    if (!molecules.has(a) && !molecules.has(b)) continue;
    maps.push(m.map);
    for (const x of [a, b]) if (x >= 0 && !molecules.has(x)) far.add(x);
  }
  for (let u = 0; u < set.units.count; u++) {
    if (set.units.alive[u] && (units.has(u) || far.has(set.units.molecule[u]))) maps.push(unitMap(set, u));
  }
  return maps;
}
