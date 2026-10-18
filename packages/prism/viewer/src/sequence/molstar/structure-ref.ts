/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A structure named as data (seqspec:libraries#the-mol-adapter): an item's `structure` field, a coverage row's, what
// "Show in 3D" loads: its source (a service's entry or a URI), its chains, assembly and model, and the Seq* complex its
// residues are read from.

import type { ServiceEntry } from '@molstar/seq-core/format';
import { KeyKind, type SeqSet } from '@molstar/seq-core/model';
import { MolstarServices } from './registry.js';

/**
 * A structure as data: what an item's JSON field `structure` holds. `complex` names the Seq* complex read from it when
 * that is not the source id in lower case.
 */
export interface StructureRef {
  readonly source: StructureServiceSource | StructureUriSource;
  readonly chains?: readonly string[];
  /** The chains' author IDs (auth_asym_id) in `chains`' order, for what a person reads; absent: `chains`. */
  readonly authChains?: readonly string[];
  readonly assembly?: string | null;
  readonly model?: number;
  readonly complex?: string;
}

export type StructureServiceSource = { readonly service: string; readonly id: string };
export type StructureUriSource = { readonly uri: string; readonly format: string };
const serviceOf = (s: StructureRef['source']): StructureServiceSource | null =>
  'service' in s ? (s as StructureServiceSource) : null;

/** The source's URL and format: a service's (`pdbe`, `afdb`) through its entry, or as written. */
export function structureUrl(ref: StructureRef, services: readonly ServiceEntry[] = MolstarServices) {
  const src = serviceOf(ref.source);
  if (!src) return { url: (ref.source as StructureUriSource).uri, format: (ref.source as StructureUriSource).format };
  const r = services.find((s) => s.name === src.service)?.resolve(src.id, 'structure');
  return r ? { url: r.uri, format: r.format } : null;
}

/** The Seq* complex a structure's residues are read from: `complex`, else its source id in lower case; null: none. */
export function complexOf(set: SeqSet, ref: StructureRef): string | null {
  const t = { kind: 0, id: 0, row: -1 };
  const name = ref.complex ?? serviceOf(ref.source)?.id.toLowerCase() ?? null;
  return name !== null && set.byName(name, t) && t.kind === KeyKind.complex ? name : null;
}
