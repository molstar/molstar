/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The coordinate systems of a resolved view (CP1b's format part): a block's piece and its layout map, and
// `ResolvedView.systems`, every coordinate system the view draws or links, those on each section's axis, and the
// closure index over the span maps of shown placements and links (seqspec:interaction#channels), built
// once at resolve from placements and links, never from hover, selection or layout.

import { type ComposedPath, type CoordinateSystem, Coords, Mappings, type SeqSet } from '@molstar/seq-core/model';

/** The section's layout coordinate system, which the model does not have: a block's shown positions, from 0. */
const LAYOUT: CoordinateSystem = { kind: -1, id: -1, row: -1 };

/**
 * A piece's layout map: `[start, stop)` of `system` laid out in direction `dir`, shown position 0 at `first`. Across a
 * cyclic origin (`stop <= start`, the domain `[0, length)`) it is two segments, `[start, length)` then `[0, stop)`, in
 * reverse order with `dir` -1; `start === stop` is the whole cyclic system from `start`.
 */
export function layoutOf(system: CoordinateSystem, start: number, stop: number, dir: number, length: number) {
  const wrap = stop <= start;
  const a = wrap ? length - start : stop - start;
  const b = wrap ? stop : 0;
  const quads = dir > 0 ? [start, a, 0, 1, 0, b, a, 1] : [start, a, b, -1, 0, b, 0, -1];
  const first = dir > 0 ? start : stop > 0 ? stop - 1 : length - 1;
  return { map: Coords.create(system, LAYOUT, 1, quads), first };
}

export interface ResolvedSystems {
  /** Every coordinate system the view draws or links. */
  readonly count: number;
  readonly kind: Uint8Array;
  readonly id: Int32Array;
  readonly row: Int32Array;
  /** Its domain. */
  readonly start: Int32Array;
  readonly stop: Int32Array;
  /** CSR: the coordinate systems on each section's axis. */
  readonly sectionStart: Int32Array;
  readonly sectionSystems: Int32Array;
  /** The span maps (span > 1) of shown placements and links, by `to` coordinate system, then support hull. */
  readonly closure: {
    readonly count: number;
    /** Into `paths`. */
    readonly path: Int32Array;
    /** Into the coordinate systems above. */
    readonly from: Int32Array;
    readonly to: Int32Array;
    readonly toStart: Int32Array;
    readonly toStop: Int32Array;
    /** Running maximum of `toStop` within each `to`: a hull search is two binary searches. */
    readonly toMaxStop: Int32Array;
  };
}

/** What `systemsOf` reads of the composed tables. */
export interface SystemsInput {
  readonly sections: readonly { readonly hidden: boolean; readonly blocks: readonly CoordinateSystem[] }[];
  readonly paths: readonly ComposedPath[];
  /** Per placement: its path (-1: direct) and whether its track and section are shown. */
  readonly placements: readonly { readonly path: number; readonly shown: boolean }[];
  readonly links: readonly { readonly path: number; readonly shown: boolean }[];
}

/** `ResolvedView.systems`. */
export function systemsOf(set: SeqSet, v: SystemsInput): ResolvedSystems {
  const list: CoordinateSystem[] = [];
  const index = new Map<string, number>();
  const idOf = (s: CoordinateSystem) => {
    const key = `${s.kind}:${s.id}:${s.row}`;
    let k = index.get(key);
    if (k === undefined) {
      k = list.push(s) - 1;
      index.set(key, k);
    }
    return k;
  };
  const sectionStart = [0];
  const sectionSystems: number[] = [];
  for (const s of v.sections) {
    const own = new Set(s.blocks.map(idOf));
    for (const k of own) sectionSystems.push(k);
    sectionStart.push(sectionSystems.length);
  }
  for (const p of v.placements) if (p.path >= 0) idOf(v.paths[p.path].from);
  for (const l of v.links) if (l.path >= 0) idOf(v.paths[l.path].from);
  // The closure maps: each span map once, however many placements and links use it.
  const closure: number[] = [];
  const seen = new Set<number>();
  const add = (path: number, shown: boolean) => {
    if (path < 0 || !shown || seen.has(path) || v.paths[path].map.span < 2) return;
    seen.add(path);
    const m = v.paths[path].map;
    closure.push(path, idOf(m.from), idOf(m.to), m.toStart, m.toStop);
  };
  for (const p of v.placements) add(p.path, p.shown);
  for (const l of v.links) add(l.path, l.shown);
  const n = closure.length / 5;
  const order = Array.from({ length: n }, (_, i) => 5 * i).sort(
    (a, b) => closure[a + 2] - closure[b + 2] || closure[a + 3] - closure[b + 3],
  );
  const col = (k: number) => Int32Array.from(order, (i) => closure[i + k]);
  const to = col(2);
  const toStop = col(4);
  const toMaxStop = new Int32Array(n);
  for (let i = 0; i < n; i++)
    toMaxStop[i] = i > 0 && to[i - 1] === to[i] ? Math.max(toMaxStop[i - 1], toStop[i]) : toStop[i];
  const domains = list.map((s) => Mappings.domainOf(set, s));
  return {
    count: list.length,
    kind: Uint8Array.from(list, (s) => s.kind),
    id: Int32Array.from(list, (s) => s.id),
    row: Int32Array.from(list, (s) => s.row),
    start: Int32Array.from(domains, (d) => d[0]),
    stop: Int32Array.from(domains, (d) => d[1]),
    sectionStart: Int32Array.from(sectionStart),
    sectionSystems: Int32Array.from(sectionSystems),
    closure: { count: n, path: col(0), from: col(1), to, toStart: col(3), toStop, toMaxStop },
  };
}
