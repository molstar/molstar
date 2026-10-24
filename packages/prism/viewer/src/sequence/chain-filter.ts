/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, Unit } from '@molstar/model/model/structure';
import type { ChainFilter } from './tracks.js';

// The chain filter (sequence.md §3): which of a structure's chain instances the section shows, as its blocks. The
// chains are named as the Mol* adapter names its units (`unit_name` `label_asym_id`: `@<operator>` off the identity,
// `molstar/complex.ts`), so a block names the unit `<complex>/<name>`. Read from the structure's units alone (no
// Seq* module, no index of its elements): the tree lists them for any structure, shown or not.

export type ChainKind = 'polymer' | 'branched' | 'ligand';

/** A chain instance the sequence can show (waters never). */
export interface SequenceChain {
  /** The unit's name in the sequence's complex. */
  readonly name: string;
  /** "A", or "C [auth A]" where the author's id differs; with its operator off the identity. */
  readonly label: string;
  readonly entity: string;
  /** The entity's description ("Cytochrome P450 3A4"), else its id. */
  readonly entityLabel: string;
  readonly kind: ChainKind;
}

const Kinds: Record<string, ChainKind | undefined> = {
  polymer: 'polymer',
  'non-polymer': 'ligand',
  macrolide: 'ligand',
  branched: 'branched',
};

const cache = new WeakMap<Structure, readonly SequenceChain[]>();

/** The structure's chain instances in its units' order, once each; waters and coarse units left out. */
export function sequenceChains(structure: Structure): readonly SequenceChain[] {
  const hit = cache.get(structure);
  if (hit) return hit;
  const out: SequenceChain[] = [];
  const seen = new Set<string>();
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const h = unit.model.atomicHierarchy;
    const entities = unit.model.entities.data;
    const op = unit.conformation.operator;
    const els = unit.elements;
    // A unit is one chain, or Mol*'s merge of waters and consecutive one-atom chains: each chain once.
    let last = -1;
    for (let i = 0; i < els.length; i++) {
      const c = h.chainAtomSegments.index[els[i]];
      if (c === last) continue;
      last = c;
      const e = h.index.getEntityFromChain(c);
      const kind = Kinds[entities.type.value(e)];
      if (!kind) continue;
      const labelAsym = h.chains.label_asym_id.value(c);
      const authAsym = h.chains.auth_asym_id.value(c);
      const name = op.isIdentity ? labelAsym : `${labelAsym}@${op.instanceId}`;
      if (seen.has(name)) continue;
      seen.add(name);
      const own = authAsym && authAsym !== labelAsym ? `${labelAsym} [auth ${authAsym}]` : labelAsym;
      const entity = entities.id.value(e);
      out.push({
        name,
        label: op.isIdentity ? own : `${own} ${op.instanceId}`,
        entity,
        entityLabel: entities.pdbx_description.value(e)?.[0] || `Entity ${entity}`,
        kind,
      });
    }
  }
  cache.set(structure, out);
  return out;
}

/**
 * The units the section shows under a filter, in the structure's order; undefined: no filter (every polymer unit, then
 * the ligands the 3D view draws: Seq*'s default).
 */
export function filteredUnits(filter: ChainFilter | undefined, chains: readonly SequenceChain[]): string[] | undefined {
  if (!filter || (!filter.polymer && filter.chains === undefined)) return undefined;
  const kept = filter.chains ? new Set(filter.chains) : undefined;
  return chains
    .filter((c) => (!filter.polymer || c.kind === 'polymer') && (!kept || kept.has(c.name)))
    .map((c) => c.name);
}

/** What a filter shows, in a few words ("chains A, B", "polymers", "every copy"); undefined: everything. */
export function describeFilter(filter: ChainFilter | undefined, chains: readonly SequenceChain[]): string | undefined {
  if (!filter) return undefined;
  const parts: string[] = [];
  if (filter.chains !== undefined) {
    // By name alone (the row is one line): the labels, author ids and all, are the inspector's.
    const names = new Set(chains.map((c) => c.name));
    const kept = filter.chains.filter((n) => names.has(n));
    parts.push(
      kept.length === 0
        ? 'no chains'
        : kept.length <= 3
          ? `chain${kept.length > 1 ? 's' : ''} ${kept.join(', ')}`
          : `${kept.length} chains`,
    );
  }
  if (filter.polymer) parts.push('polymers');
  if (filter.each) parts.push('every copy');
  return parts.length ? parts.join(', ') : undefined;
}
