/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's stream description at prepare time (`trackStream`, `lettersStream`) and the draw-side helpers over it
// (`streamedBlock`, `streamStats`) (seqspec:viewer#extension-points).

import { type StreamedData, StreamProp } from '@molstar/seq-core/format';
import { type ComposedPath, type CoordinateSystem, CoordinateSystemKind, type SeqSet } from '@molstar/seq-core/model';
import type { ShownIndex } from './axis.js';
import type { DrawView } from './draw/view.js';
import type { LayerInput } from './registry.js';
import type { StreamLayer } from './stream-contract.js';

/**
 * A layer's streamed annotations (`svs:stream`, one source a layer): its description, per placed block its system, by
 * the block's unit or, for data over an alignment's columns (a served statistic), its alignment.
 */
export function trackStream(input: LayerInput): StreamLayer | null {
  const { set, layer, section } = input;
  let data: StreamedData | undefined;
  let systems: Int32Array | null = null;
  let paths: (ComposedPath | null)[] | null = null;
  let routed = '';
  for (const pl of layer.placements) {
    const a = set.annotations[pl.annotation];
    const d = a?.streamed ? (a.custom.props.get(StreamProp) as StreamedData | undefined) : undefined;
    const blk = section.blocks[pl.block];
    // Through a path (ST7) the system is the path's start, the annotation's target; else the block's own.
    const path = pl.path >= 0 ? layer.paths[pl.path] : null;
    const key = path
      ? systemKey(set, path.from)
      : !blk
        ? -1
        : blk.unit >= 0
          ? set.units.key[blk.unit]
          : (set.alignments[blk.alignment]?.key ?? -1);
    if (!d || (data && d.source !== data.source) || key < 0) continue;
    data = d;
    systems ??= new Int32Array(section.blocks.length).fill(-1);
    systems[pl.block] = d.systems.get(set.strings.get(key)) ?? -1;
    if (path) paths ??= new Array(section.blocks.length).fill(null);
    if (paths) paths[pl.block] = path;
    if (!path) continue;
    // The path's segments in the key: a new map (an edited CDS, a strand flipped) places the held items again.
    routed += `|${pl.block}:${mapHash(path)}`;
  }
  if (!data || !systems) return null;
  const collapse =
    typeof input.params.collapse === 'string' && input.params.collapse !== '' ? input.params.collapse : null;
  return {
    key: `${data.source.ref}#${systems.join()}${routed}${collapse === null ? '' : `|${collapse}`}`,
    rows: 0,
    columns: 0,
    remote: data,
    local: null,
    systems,
    collapse,
    paths,
  };
}

/** A hash of a path's segments (FNV-1a over at, length, toLo, dir; edge-rate). */
function mapHash(p: ComposedPath): string {
  const m = p.map;
  let h = 0x811c9dc5;
  for (let k = 0; k < m.count; k++)
    for (const x of [m.at[k], m.length[k], m.toLo[k], m.dir[k]]) h = Math.imul(h ^ x, 0x01000193) >>> 0;
  return `${m.count}.${m.span}.${h.toString(36)}`;
}

/** A coordinate system's key: a unit's, a molecule's, or an alignment's for its columns; -1: none. */
function systemKey(set: SeqSet, s: CoordinateSystem): number {
  const K = CoordinateSystemKind;
  if (s.kind === K.unit) return set.units.key[s.id];
  if (s.kind === K.molecule) return set.molecules[s.id]?.key ?? -1;
  return s.kind === K.columns ? (set.alignments[s.id]?.key ?? -1) : -1;
}

/**
 * A letters layer's streamed blocks (ST3a): per block the system of its unit's lazy molecule, by its `svs:stream` (one
 * stream a layer); null: none.
 */
export function lettersStream(input: LayerInput): StreamLayer | null {
  const { set, section } = input;
  let data: StreamedData | undefined;
  let systems: Int32Array | null = null;
  for (let b = 0; b < section.blocks.length; b++) {
    const u = section.blocks[b].unit;
    const mol = u >= 0 ? set.molecules[set.units.molecule[u]] : undefined;
    const d = mol?.lazy ? (mol.custom.props.get(StreamProp) as StreamedData | undefined) : undefined;
    if (!mol || !d || (data && d.source !== data.source)) continue;
    data = d;
    systems ??= new Int32Array(section.blocks.length).fill(-1);
    systems[b] = d.systems.get(set.strings.get(mol.key)) ?? -1;
  }
  if (!data || !systems) return null;
  return { key: `${data.source.ref}#${systems.join()}`, rows: 0, columns: 0, remote: data, local: null, systems };
}

/** Draw-side: block b of the drawn layer is streamed (a 1-D layer's block showing one of its systems). */
export function streamedBlock(v: DrawView, b: number): boolean {
  const sys = v.stream?.source.systems;
  return sys != null && b < sys.length && sys[b] >= 0;
}
/** Draw-side: block b's shown [k0, k1) from the layer's regions (`StreamView.values`) into `out`; count 0: none loaded. */
export function streamStats(v: DrawView, b: number, k0: number, k1: number, out: Float64Array) {
  const ax = v.axis;
  const st = v.stream;
  out[0] = out[1] = Number.NaN;
  out[2] = out[3] = 0;
  if (!ax || !st?.values) return;
  const a = ax.indexOf(b, k0 as ShownIndex);
  const c = ax.indexOf(b, (k1 - 1) as ShownIndex);
  out[0] = 1 / (v.cell * v.span[b] * v.dpr);
  if (st.values(b, Math.min(a, c), Math.max(a, c) + 1, out)) return;
  out[0] = out[1] = Number.NaN;
  out[2] = out[3] = 0;
}
