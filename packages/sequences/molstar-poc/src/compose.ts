/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An MVS state composed from Seq* data (seqspec:libraries#the-mol-adapter): a structure named as data (an item's
// `structure` field, an annotation's `custom.molstar`), its chains, and colored Seq* runs carried onto its residues
// through the document's own residue maps (the units' label numbering, the mappings such as SIFTS), as MVS color
// tables. Pure data: no Mol* structure, no canvas; the viewer half feeds it a view's colors (viewer.ts).

import type { ServiceEntry } from '@molstar/seq-core/format';
import {
  type CoordinateSystem,
  CoordinateSystemKind,
  Extents,
  KeyKind,
  MoleculeKind,
  type SeqSet,
} from '@molstar/seq-core/model';
import { Translator } from '@molstar/seq-core/translate';
import type { MVSData_State } from 'molstar/lib/extensions/mvs/mvs-data';
import { type MolstarDefaults, MvsColorTable, structureMvs } from './mvs.js';
import { MolstarServices } from './registry.js';
import { complexUnits, seqRoutes } from './seq-units.js';

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

/**
 * Seq* runs onto one complex's polymer units, as `label_seq_id` ranges per chain: through a translator holding the
 * units' maps onto their molecules and every mapping of those molecules (a run on P04637 reaches 1TUP's chains C to E
 * by its SIFTS mapping). Made per composition (edge rate).
 */
export class ResidueProjector {
  private readonly t = new Translator();
  private readonly units: { system: CoordinateSystem; chain: string; first: number }[] = [];
  private readonly one = Extents.create(1);
  private readonly out = Extents.create(8);

  constructor(set: SeqSet, complex: string) {
    const ids = new Set<number>();
    const t = { kind: 0, id: 0, row: -1 };
    for (const u of complexUnits(set, complex)) {
      if (!set.byName(u.ref, t) || t.kind !== KeyKind.unit || !u.mvs.label_asym_id) continue;
      if (set.molecules[set.units.molecule[t.id]]?.kind !== MoleculeKind.polymer) continue;
      ids.add(t.id);
      this.units.push({
        system: { kind: CoordinateSystemKind.unit, id: t.id, row: -1 },
        chain: u.mvs.label_asym_id,
        first: u.first,
      });
    }
    for (const m of seqRoutes(set, ids)) this.t.addMap(m);
  }

  /** Each run's residues appended to `table` in its color; a residue a later run reaches is listed again (it wins). */
  add(runs: Extents, colors: readonly string[], table: MvsColorTable): MvsColorTable {
    for (let i = 0; i < runs.count; i++) {
      this.one.count = 0;
      const from = { kind: runs.kind[i], id: runs.id[i], row: runs.row[i] };
      Extents.add(this.one, from, runs.start[i], runs.stop[i]);
      for (const u of this.units) {
        const out = this.out;
        out.count = 0;
        const same = from.kind === u.system.kind && from.id === u.system.id;
        if (same) Extents.add(out, u.system, runs.start[i], runs.stop[i]);
        else this.t.translate(this.one, from, u.system, out);
        for (let k = 0; k < out.count; k++) {
          const last = table.color.length - 1;
          const [beg, end] = [out.start[k] + u.first, out.stop[k] - 1 + u.first];
          if (
            last >= 0 &&
            table.label_asym_id[last] === u.chain &&
            table.color[last] === colors[i] &&
            table.end_label_seq_id[last] + 1 === beg
          ) {
            table.end_label_seq_id[last] = end;
          } else {
            table.label_asym_id.push(u.chain);
            table.beg_label_seq_id.push(beg);
            table.end_label_seq_id.push(end);
            table.color.push(colors[i]);
          }
        }
      }
    }
    return table;
  }
}

export interface ComposeInput {
  readonly set: SeqSet;
  readonly structure: StructureRef;
  /** Colored runs per layer, in order: a later layer wins; `sticks`: its residues as ball and stick too (sites). */
  readonly layers?: readonly {
    readonly runs: Extents;
    readonly colors: readonly string[];
    readonly sticks?: boolean;
  }[];
  readonly selection?: { readonly runs: Extents; readonly label: string; readonly color?: string } | null;
  readonly defaults?: MolstarDefaults;
  readonly services?: readonly ServiceEntry[];
  /** The 3D canvas's background (the host's surface). */
  readonly background?: string;
}

/** The MVS state for a structure and the Seq* colors and selection over it; null: its source resolves to no URL. */
export function composeStructureMvs(
  input: ComposeInput,
): { state: MVSData_State; tables: MvsColorTable[]; sticks: MvsColorTable[] } | null {
  const { set, structure } = input;
  const at = structureUrl(structure, input.services);
  if (!at) return null;
  const complex = complexOf(set, structure);
  const projector = complex ? new ResidueProjector(set, complex) : null;
  const projected = projector
    ? (input.layers ?? [])
        .map((l) => ({ table: projector.add(l.runs, l.colors, MvsColorTable.create()), sticks: !!l.sticks }))
        .filter((t) => MvsColorTable.size(t.table) > 0)
    : [];
  const tables = projected.map((t) => t.table);
  const sticks = projected.filter((t) => t.sticks).map((t) => t.table);
  const sel = input.selection;
  const selection =
    projector && sel
      ? {
          residues: projector.add(sel.runs, Array(sel.runs.count).fill(sel.color ?? '#e8a33d'), MvsColorTable.create()),
          label: sel.label,
        }
      : null;
  const state = structureMvs({
    url: at.url,
    format: at.format as 'bcif',
    model: structure.model,
    assembly: structure.assembly ?? null,
    chains: structure.chains?.map((c) => ({ label_asym_id: c })) ?? null,
    colors: tables,
    sticks,
    selection,
    defaults: input.defaults,
    background: input.background,
  });
  return { state, tables, sticks };
}
