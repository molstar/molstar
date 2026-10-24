/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, Unit } from '@molstar/model/model/structure';
import type { ChainFilter } from './tracks.js';

// The chain filter (sequence.md §3): which of a structure's chains the section shows, as its blocks. The chains are
// named as the Mol* adapter names the complex's units (`unit_name` `label_asym_id`, `molstar/complex.ts`), so a block
// names the unit `<complex>/<name>`. Read from the structure's units alone (no Seq* module, no index of its elements):
// the tree lists them for any structure, shown or not.

export type ChainKind = 'polymer' | 'branched' | 'ligand';

/** A chain the sequence can show (waters never). */
export interface SequenceChain {
  /** The unit's name in the sequence's complex: its `label_asym_id`. */
  readonly name: string;
  /** The author's chain id. */
  readonly auth: string;
  /** "A", or "C [auth A]" where the author's id differs. */
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

const cache = new WeakMap<object, Map<string, readonly SequenceChain[]>>();
const byStructure = new WeakMap<Structure, readonly SequenceChain[]>();

/**
 * The chains of the structure's model, as the sequence's complex has them: the complex is made of the deposited model
 * (one unit a `label_asym_id`, `molstar/complex.ts`), so an assembly's copies are its chains once, and a block names
 * a unit there is. In the complex's order: by entity, then as the chains come. Cached by the model's hierarchy and the
 * chains the structure has (a trajectory's frames share it).
 */
export function sequenceChains(structure: Structure): readonly SequenceChain[] {
  const first = structure.units.find(Unit.isAtomic);
  if (!first) return [];
  const h = first.model.atomicHierarchy;
  let byChains = cache.get(h);
  if (!byChains) cache.set(h, (byChains = new Map()));
  // The chains present, by index (an assembly or a selection may hold fewer than the model).
  // A unit is one chain, its first and last atoms say (O(units) a frame); only one Mol* merged (waters, one-atom
  // chains) is walked.
  const seen = byStructure.get(structure);
  if (seen) return seen;
  const present = new Set<number>();
  const at = h.chainAtomSegments.index;
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit) || unit.model.atomicHierarchy !== h) continue;
    const els = unit.elements;
    const a = at[els[0]];
    const b = at[els[els.length - 1]];
    present.add(a);
    if (a === b) continue;
    let last = a;
    for (let i = 1; i < els.length; i++) {
      const c = at[els[i]];
      if (c !== last) present.add((last = c));
    }
  }
  const order = [...present].sort((a, b) => a - b);
  const key = order.join(',');
  const hit = byChains.get(key);
  if (hit) {
    byStructure.set(structure, hit);
    return hit;
  }
  const entities = first.model.entities;
  const out: { chain: SequenceChain; entity: number }[] = [];
  for (const c of order) {
    const e = h.index.getEntityFromChain(c as never);
    const kind = Kinds[entities.data.type.value(e)];
    if (!kind) continue;
    const name = h.chains.label_asym_id.value(c as never);
    const auth = h.chains.auth_asym_id.value(c as never);
    const entity = entities.data.id.value(e);
    out.push({
      entity: e,
      chain: {
        name,
        auth,
        label: auth && auth !== name ? `${name} [auth ${auth}]` : name,
        entity,
        entityLabel: entities.data.pdbx_description.value(e)?.[0] || `Entity ${entity}`,
        kind,
      },
    });
  }
  // The complex's order: entity, then chain (a stable sort keeps the chains' order within an entity).
  out.sort((a, b) => a.entity - b.entity);
  const list = out.map((x) => x.chain);
  byChains.set(key, list);
  byStructure.set(structure, list);
  return list;
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
