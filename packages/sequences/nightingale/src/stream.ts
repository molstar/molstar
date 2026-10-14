/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Streamed values (seqspec:libraries#nightingale): Nightingale cannot stream, so
// `materialize` stays synchronous and a `values` layer over a streamed source describes what to fetch, and `mount`
// fetches it once through the source's provider, at the level whose bins are the finest at least a pixel wide at the
// container's width (else the coarsest). Platform-free: no DOM, no Nightingale.

import { type StreamedData, StreamNotServed, streamRequest, TILE } from '@molstar/seq-core/format';
import type { Issue } from '@molstar/seq-core/model';
import { type AxisBlock, indexAt } from './axis.js';
import { nightingaleIssue } from './issues.js';

/** One placement of a streamed `values` annotation on the axis: its source, the manifest system and block it shows. */
export interface StreamedPart {
  readonly stream: StreamedData;
  /** Its coordinate system's index in the manifest. */
  readonly system: number;
  /** The block it is placed on (directly: its indices are the annotation's). */
  readonly block: AxisBlock;
}
/** A linegraph's streamed data (its `svs` property): the layer's ref for issues, the clamp, the parts. */
export interface StreamedValues {
  readonly ref: string;
  /** The resolved domain values are clamped to; null: none. */
  readonly domain: readonly [number, number] | null;
  readonly parts: readonly StreamedPart[];
}

/** Without a measured width: the level for this many pixels. */
const DEFAULT_WIDTH = 1000;

/** The finest of `widths` (positions a bin, ascending) at least `perPixel` wide; else the coarsest. */
export function streamedLevel(widths: readonly (number | readonly number[])[], axisLength: number, width: number) {
  const perPixel = axisLength / (width > 0 ? width : DEFAULT_WIDTH);
  const k = widths.findIndex((w) => (w as number) >= perPixel);
  return k < 0 ? widths.length - 1 : k;
}

/**
 * The values on an axis of `axisLength` positions: a bin's mean (else the middle of its min and max) on every position
 * it covers, clamped to the domain, null where nothing is served. A failed region leaves its positions null and is an
 * issue; an aborted fetch reports nothing.
 */
export async function fetchStreamed(
  p: StreamedValues,
  axisLength: number,
  width: number,
  signal?: AbortSignal,
): Promise<{ values: (number | null)[]; level: number; issues: Issue[] }> {
  const values: (number | null)[] = new Array(axisLength).fill(null);
  const issues: Issue[] = [];
  // Every tile that failed, said in one issue with the first few (GV-07: one per tile, each its own text).
  const failed: string[] = [];
  const [lo, hi] = p.domain ?? [-Infinity, Infinity];
  let level = -1;
  await Promise.all(
    p.parts.map(async ({ stream: s, system, block: b }) => {
      const { key, extent } = s.manifest.systems[system];
      const lvl = streamedLevel(s.manifest.levels, axisLength, width);
      level = lvl;
      const w = s.manifest.levels[lvl] as number;
      const span = TILE * w;
      // The tiles the block's indices fall in: its two ends would not bound them on a reversed piece or across a
      // cyclic origin (NG-05). One pass over the block, once a mount.
      const tiles = new Set<number>();
      for (let k = 0; k < b.count; k++) {
        const i = indexAt(b, k);
        if (i >= extent[0] && i < extent[1]) tiles.add(Math.floor(i / span));
      }
      const regions = new Map<number, ArrayLike<number>[]>();
      const fetches: Promise<void>[] = [];
      for (const tx of tiles) {
        const [start, stop] = [tx * span, Math.min((tx + 1) * span, extent[1])];
        const at = `${s.source.uri ?? s.source.stream ?? s.source.ref} ${key} [${start}, ${stop})`;
        fetches.push(
          s.provider
            .fetch(
              streamRequest(s.source.ref, key, { start, stop, level: lvl }),
              signal ?? new AbortController().signal,
            )
            .then(
              (r) => {
                const t = r.kind === 'bins' ? r.bins : null;
                if (t?.mean && t.min && t.max) regions.set(tx, [t.min as never, t.max as never, t.mean as never]);
                else failed.push(`${at}: no bins`);
              },
              (e: unknown) => {
                if (signal?.aborted) return;
                const why = e instanceof StreamNotServed ? 'is not served' : e instanceof Error ? e.message : String(e);
                failed.push(`${at}: ${why}`);
              },
            ),
        );
      }
      await Promise.all(fetches);
      for (let k = 0; k < b.count; k++) {
        const i = indexAt(b, k);
        const t = i >= extent[0] && i < extent[1] ? regions.get(Math.floor(i / span)) : undefined;
        if (!t) continue;
        const j = Math.floor((i - Math.floor(i / span) * span) / w);
        const v = Number.isNaN(t[2][j]) ? (t[0][j] + t[1][j]) / 2 : t[2][j];
        if (!Number.isNaN(v)) values[b.offset + k] = Math.min(hi, Math.max(lo, v));
      }
    }),
  );
  if (failed.length > 0) {
    const more = failed.length > 3 ? ` and ${failed.length - 3} more` : '';
    const message = `Nightingale: ${failed.length} tile(s) failed: ${failed.slice(0, 3).join('; ')}${more}`;
    issues.push(nightingaleIssue('svs.stream.load', p.ref, message));
  }
  return { values, level, issues: signal?.aborted ? [] : issues };
}
