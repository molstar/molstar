/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// BigWig: values (README, "BigWig"). The file's zoom levels are the manifest's levels; a level request resamples the
// zoom records onto the request's bins, level -1 is the data itself.

import { type DataPatch, StreamNotServed, type StreamProvider, type StreamRequest } from '@molstar/seq-core/format';
import { type Chrom, openBbi } from './bbi.js';
import { once, type RangeReader, view } from './ranges.js';

/** A BigWig file as a `values` provider: one coordinate system per chromosome, `levels` the file's zoom reductions. */
export function bigWigProvider(read: RangeReader): StreamProvider {
  const open = once((signal) => openBbi(read, signal));
  return {
    async fetch(r: StreamRequest, signal: AbortSignal): Promise<DataPatch> {
      const bbi = await open(signal);
      if (bbi.kind !== 'bigwig') throw new Error('not a BigWig file (this is a BigBed)');
      const at = {
        source: r.source,
        system: r.system,
        start: r.start,
        stop: r.stop,
        row_start: 0,
        row_stop: 0,
        level: r.level,
      };
      if (r.system === null) {
        const systems = [...bbi.chroms.values()].map((c) => ({
          key: c.name,
          extent: [0, c.size] as [number, number],
          numbering: null,
          label: c.name,
        }));
        const manifest = {
          systems,
          levels: bbi.zooms.map((z) => z.reduction),
          domain: bbi.range ?? ([0, 0] as [number, number]),
        };
        return { ...at, kind: 'manifest', manifest };
      }
      const chrom = bbi.chroms.get(r.system);
      if (!chrom) throw new Error(`BigWig: no chromosome ${r.system}`);
      if (r.level >= bbi.zooms.length)
        throw new StreamNotServed(`BigWig: no level ${r.level} (the file has ${bbi.zooms.length})`);
      if (r.start < 0 || r.stop < r.start) throw new Error(`BigWig: bad region [${r.start}, ${r.stop})`);
      const { start, stop } = r;
      if (r.level < 0) {
        const value = new Float64Array(Math.max(0, Math.min(stop, chrom.size) - start)).fill(Number.NaN);
        for (const b of await bbi.blocks(bbi.dataIndexOffset, chrom, start, stop, signal))
          decode(b, bbi.le, chrom, (s, e, v) =>
            value.fill(v, Math.max(s, start) - start, Math.min(e, stop, chrom.size) - start),
          );
        return { ...at, kind: 'items', items: { value } };
      }
      const zoom = bbi.zooms[r.level];
      const w = zoom.reduction;
      const n = Math.ceil((stop - start) / w);
      const [min, max] = [
        new Float64Array(n).fill(Number.POSITIVE_INFINITY),
        new Float64Array(n).fill(Number.NEGATIVE_INFINITY),
      ];
      const [sum, count] = [new Float64Array(n), new Float64Array(n)];
      for (const b of await bbi.blocks(zoom.indexOffset, chrom, start, stop, signal))
        zoomRecords(b, bbi.le, chrom, (s, e, valid, lo, hi, total) => {
          spread(s, e, start, w, n, (j, share) => {
            count[j] += valid * share;
            sum[j] += total * share;
            min[j] = Math.min(min[j], lo);
            max[j] = Math.max(max[j], hi);
          });
        });
      for (let j = 0; j < n; j++) if (!(count[j] > 0)) min[j] = max[j] = Number.NaN; // no value in the bin
      const mean = sum.map((t, j) => t / count[j]);
      return { ...at, kind: 'bins', bins: { min, max, mean } };
    },
  };
}

/** A data block's items as (start, end, value) in the chromosome `chrom`: bedGraph, variableStep or fixedStep. */
function decode(b: Uint8Array, le: boolean, chrom: Chrom, emit: (start: number, end: number, value: number) => void) {
  const dv = view(b);
  const [id, first, step, span, type, count] = [
    dv.getUint32(0, le),
    dv.getUint32(4, le),
    dv.getUint32(12, le),
    dv.getUint32(16, le),
    b[20],
    dv.getUint16(22, le),
  ];
  if (id !== chrom.id) return;
  for (let i = 0; i < count; i++) {
    if (type === 1) emit(dv.getUint32(24 + 12 * i, le), dv.getUint32(28 + 12 * i, le), dv.getFloat32(32 + 12 * i, le));
    else if (type === 2)
      emit(dv.getUint32(24 + 8 * i, le), dv.getUint32(24 + 8 * i, le) + span, dv.getFloat32(28 + 8 * i, le));
    else if (type === 3) emit(first + i * step, first + i * step + span, dv.getFloat32(24 + 4 * i, le));
    else throw new Error(`BigWig: unknown section type ${type}`);
  }
}

/**
 * A record [s, e) is not on a level's grid (`n` bins of `w` from `start`): `onBin(j, share)` for each bin it overlaps,
 * `share` the part of the record in it (CB-12: BigWig's and BigBed's zoom records both spread so).
 */
export function spread(
  s: number,
  e: number,
  start: number,
  w: number,
  n: number,
  onBin: (bin: number, share: number) => void,
): void {
  for (let j = Math.max(0, Math.floor((s - start) / w)); j < n && start + j * w < e; j++) {
    const share = (Math.min(e, start + (j + 1) * w) - Math.max(s, start + j * w)) / (e - s);
    if (share > 0) onBin(j, share);
  }
}

/** A zoom block's 32-byte records: (start, end, valid bases, min, max, sum) of the chromosome `chrom`. */
export function zoomRecords(
  b: Uint8Array,
  le: boolean,
  chrom: Chrom,
  emit: (start: number, end: number, valid: number, min: number, max: number, sum: number) => void,
) {
  const dv = view(b);
  for (let o = 0; o + 32 <= b.length; o += 32) {
    if (dv.getUint32(o, le) !== chrom.id) continue;
    emit(
      dv.getUint32(o + 4, le),
      dv.getUint32(o + 8, le),
      dv.getUint32(o + 12, le),
      dv.getFloat32(o + 16, le),
      dv.getFloat32(o + 20, le),
      dv.getFloat32(o + 24, le),
    );
  }
}
