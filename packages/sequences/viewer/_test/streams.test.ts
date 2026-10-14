/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stream manager's region cache (LY4; seqspec:sources#the-viewers-side, seqspec:sources#levels-and-the-zoom-gate),
// headless: the tile-manager block (ST1's migration tests: behavior, not internals) and its reductions; then
// LY4's: the levels under a budget, local compute past one tile and slot -1, the reference server against
// `reduce`, each cell read once, a frame's local work bounded (no stall), content in the key, retries, cooldown
// and one report per outage, cancellation, eviction with the view pinned, and a filtered million-column alignment.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  type Loaders,
  type SeqViewSpecDoc,
  type StreamedData,
  StreamProp,
  TILE,
  type Tile,
} from '@molstar/seq-core/format';
import type { Issue } from '@molstar/seq-core/model';
import { NotServed, type TileTransport, tilesStream } from '@molstar/seq-core/stream/tiles';
import { proceduralRow, tileServer } from '@molstar/seq-core/stream/tiles/server';
import { describe, expect, it } from 'vitest';
import { mismatches } from '../../core/_test/model/mismatches.js';
import { RegionState, SeqViewerConfig } from '../src/index.js';
import type { RowsState } from '../src/rows.js';
import { reduce, type Streams, tileLevels } from '../src/streams.js';
import { layer, opened, preparedRenderers } from './layout-fixtures.js';

const P = { renderers: preparedRenderers() };
const flush = () => new Promise((ok) => setTimeout(ok, 0));
const URI = 'https://example.org/lazy';

/** Loaders over the procedural server whose tile answers wait for `release` (or come at once when `open`). */
function gated(rows: number, columns: number, transport?: TileTransport) {
  const server = tileServer({ rowCount: rows, columnCount: columns });
  const inner = transport ?? server.transport;
  const waiting: (() => void)[] = [];
  const calls: string[] = [];
  const signals: (AbortSignal | undefined)[] = [];
  const g = {
    open: true,
    calls,
    signals,
    release: () => {
      for (const f of waiting.splice(0)) f();
    },
  };
  const loaders: Loaders = {
    load: (r) => Promise.reject(new Error(`no fixture ${r.uri}`)),
    stream: tilesStream((u) =>
      u === URI
        ? (path, signal) => {
            calls.push(path);
            signals.push(signal);
            return g.open || !path.startsWith('tiles/')
              ? inner(path, signal)
              : new Promise((ok, no) => waiting.push(() => inner(path, signal).then(ok, no)));
          }
        : null,
    ),
  };
  return { g, loaders };
}
function lazyDoc(rows: number, columns: number, o: { second?: boolean; filter?: [number, number][] } = {}) {
  const b = createSeqViewSpec({ title: 'lazy' });
  const aln = b.data.alignment({
    ref: 'lazy',
    row_count: rows,
    column_count: columns,
    source: { uri: URI, protocol: 'seqstar-tiles/1' },
  });
  const f = o.filter
    ? b.data.annotation({
        ref: 'kept',
        geometry: 'ranges',
        target: aln,
        value: { start: o.filter.map((r) => r[0]), stop: o.filter.map((r) => r[1]) },
      })
    : undefined;
  const v = b.view({ ref: 'v', label: 'A', style: 'bare', height: 'box', ...(f ? { positions: f } : {}) });
  v.block({ alignment: aln });
  v.track({ label: 'Rows' }).rows();
  if (o.second) b.view({ ref: 'w', label: 'B', style: 'bare' }).block({ alignment: aln }).track({ label: 'R' }).ruler();
  return b.getDocument();
}
/** An alignment by value (the Storybook's large one): a W every 37 columns, gaps at about one in seven. */
function byValue(rows: number, columns: number, salt = 0) {
  const AA = 'ACDEFGHIKLMNPQRSTVWY';
  const lines = ['# STOCKHOLM 1.0'];
  for (let k = 0; k < rows; k++) {
    let r = '';
    for (let c = 0; c < columns; c++) {
      const h = Math.imul((k % 97) * 73856093 + c * 19349663 + (k >> 5) + salt, 2654435761) >>> 0;
      r += c % 37 === 5 ? 'W' : h % 7 === 0 ? '-' : AA[(c * 7 + (h % 5)) % 20];
    }
    lines.push(`R${k}/1-${columns} ${r}`);
  }
  lines.push('//', '');
  return lines.join('\n');
}
function largeDoc(rows: number, columns: number, salt = 0): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'large' });
  const aln = b.data.alignment({
    ref: 'big',
    source: { value: byValue(rows, columns, salt), format: 'stockholm', max_bytes: 1 << 26 },
  });
  const v = b.view({ ref: 'large', label: 'Large', style: 'bare', height: 'box' });
  v.block({ alignment: aln });
  v.track({ label: 'Rows' }).rows();
  return b.getDocument();
}
const tiles = (g: { calls: string[] }) => g.calls.filter((c) => c.startsWith('tiles/'));
const sm = (t: { ctx: { streams: unknown } }) => t.ctx.streams as Streams;
const ref = (t: Awaited<ReturnType<typeof opened>>) => {
  const { k, m } = layer(t.ctx, 'Rows');
  return m.set.strings.get(m.resolved.layers.ref[k]);
};
/** Two frames: the viewport's, then the stream manager's it asks for. */
const frames = (t: { scheduler: { advance(ms: number): void } }) => {
  t.scheduler.advance(16);
  t.scheduler.advance(16);
};
const blank = (): Tile => ({ cls: new Uint8Array(TILE * TILE), share: new Uint8Array(TILE * TILE) });

/** The brute force of a sample: the most frequent residue (ties to the lower code), else the gap; its share of all. */
function brute(
  codes: Uint8Array,
  slot: Int32Array,
  C: number,
  N: number,
  ly: number,
  lx: number,
  y: number,
  x: number,
) {
  const count = new Map<number, number>();
  let n = 0;
  for (let k = y * 2 ** ly; k < Math.min(N, (y + 1) * 2 ** ly); k++)
    for (let c = x * 2 ** lx; slot[k] >= 0 && c < Math.min(C, (x + 1) * 2 ** lx); c++) {
      const code = codes[slot[k] * C + c];
      count.set(code, (count.get(code) ?? 0) + 1);
      n++;
    }
  if (n === 0) return [0, 0];
  const order = (xs: [number, number][]) => xs.sort((a, b) => b[1] - a[1] || a[0] - b[0])[0];
  const residues = [...count].filter(([code]) => code !== 45 && code !== 46);
  const [best, m] = residues.length ? order(residues) : order([...count]);
  return [best, Math.round((m / n) * 255)];
}

describe('the region cache, headless on the in-process server (ported)', () => {
  it('loads only the tiles in view, coarse before fine, and none again for a cached view', async () => {
    const { g, loaders } = gated(40_000, 3000);
    const t = await opened(lazyDoc(40_000, 3000), { ...P, loaders, size: { 'v.main': [1000, 600] } });
    const s = sm(t);
    const vp = t.vp('v.main');
    vp.zoomTo(0, 0, 300);
    s.pass();
    await flush();
    const first = tiles(g);
    expect(first.length).toBeGreaterThan(1);
    expect(first.length).toBeLessThanOrEqual(6);
    const level = (p: string) => p.split('/').slice(1, 3).join('/');
    expect(level(first[0]) !== level(first[first.length - 1])).toBe(true);
    s.pass();
    await flush();
    expect(tiles(g).length).toBe(first.length);
    // A pan by a screen: the new view's tiles only.
    vp.panBy(1000, 0, 'api');
    s.pass();
    await flush();
    const panned = tiles(g).slice(first.length);
    expect(panned.length).toBeGreaterThan(0);
    for (const p of panned) expect(Number(p.split('/')[4]) * 256 * 2 ** Number(p.split('/')[2])).toBeGreaterThan(256);
    // A zoom out: coarser levels, a bounded number of tiles (the sample budget).
    vp.zoomTo(0, 0, 3000);
    vp.setRowScale(0.01, 0);
    t.layout('v.main');
    s.pass();
    await flush();
    const zoomed = tiles(g).slice(first.length + panned.length);
    expect(zoomed.length).toBeGreaterThan(0);
    expect(zoomed.length).toBeLessThanOrEqual(12);
    expect(s.stats.current.inFlight).toBe(0);
  }, 30_000); // the in-process server reduces every cell of 40,000 × 3,000 for the zoomed-out tiles

  it("ST2: a source's alignment takes a margin: the ring around the view's tiles, after them, none at Margin 0", async () => {
    const asked = async (margin: number) => {
      const { g, loaders } = gated(40_000, 3000);
      const t = await opened(lazyDoc(40_000, 3000), { ...P, loaders, size: { 'v.main': [1000, 600] } });
      t.ctx.config.set(SeqViewerConfig.Stream.Margin, margin);
      t.ctx.config.set(SeqViewerConfig.Stream.MaxInFlight, 100);
      t.vp('v.main').zoomTo(0, 0, 300);
      t.vp('v.main').scrollTo(t.vp('v.main').live.current.scrollX, 20_000 * 22, 'api');
      sm(t).pass();
      return tiles(g);
    };
    const [none, ring] = [await asked(0), await asked(1)];
    expect(none.length).toBeGreaterThan(0);
    // The view's tiles first, in the order they had; then tiles around them at the same level.
    expect(ring.slice(0, none.length)).toEqual(none);
    const rest = ring.slice(none.length);
    expect(rest.length).toBeGreaterThan(0);
    expect(new Set(rest.map((p) => p.split('/').slice(1, 3).join('/'))).size).toBe(1);
    expect(new Set(ring).size).toBe(ring.length);
  }, 30_000);

  it("supersession: an older view's request is cancelled, never its queue; eviction keeps the bound", async () => {
    const { g, loaders } = gated(40_000, 3000);
    const t = await opened(lazyDoc(40_000, 3000), { ...P, loaders, size: { 'v.main': [1000, 600] } });
    t.ctx.config.set(SeqViewerConfig.Stream.MaxInFlight, 1);
    t.ctx.config.set(SeqViewerConfig.Stream.CacheBytes, 3 * 2 * TILE * TILE);
    const s = sm(t);
    const vp = t.vp('v.main');
    vp.zoomTo(0, 0, 300);
    g.open = false;
    s.pass();
    expect(tiles(g).length).toBe(1);
    // Far away before anything arrived: the first view's request is aborted, its queue dropped.
    vp.scrollTo(2500 * vp.live.current.cell, 30_000 * 22, 'api');
    s.pass();
    const first = g.calls.findIndex((p) => p.startsWith('tiles/'));
    expect([s.stats.current.cancelled, g.signals[first]?.aborted]).toEqual([1, true]);
    g.release();
    await flush();
    g.open = true;
    for (let i = 0; i < 6; i++) {
      s.pass();
      await flush();
    }
    const after = tiles(g).slice(1);
    expect(after.length).toBeGreaterThan(0);
    for (const p of after) expect(Number(p.split('/')[3])).toBeGreaterThan(0);
    expect(s.stats.current.requests).toBeGreaterThan(3);
    // The view's own regions are pinned, its labels' chunk too (VT-07); at most a tile more than the bound while it fills.
    let labels = 0;
    for (const a of t.ctx.graph.assembled.value?.set.alignments ?? [])
      for (const x of (a?.custom.props.get(StreamProp) as StreamedData | undefined)?.labels.values() ?? [])
        labels += 2 * x.length;
    expect(s.stats.current.bytes).toBeLessThanOrEqual(Math.max(3, after.length) * 2 * TILE * TILE + labels);
  });

  it('pauses while its view is not current (cache kept), closes when its layer leaves and on dispose', async () => {
    const { g, loaders } = gated(5000, 500);
    const doc = lazyDoc(5000, 500, { second: true });
    const t = await opened(doc, { ...P, loaders, size: { 'v.main': [1000, 400], 'w.main': [1000, 400] } });
    const s = sm(t);
    s.pass();
    await flush();
    const n = tiles(g).length;
    expect(n).toBeGreaterThan(0);
    const r = ref(t);
    const view = s.get(r);
    t.ctx.views.setCurrent('w');
    await t.ctx.graph.update();
    s.pass();
    await flush();
    expect([tiles(g).length, s.get(r) === view, (view?.version ?? 0) > 0]).toEqual([n, true, true]);
    t.ctx.views.setCurrent('v');
    await t.ctx.graph.update();
    s.pass();
    await flush();
    expect(tiles(g).length).toBe(n);
    await t.ctx.setDocument({
      ...doc,
      root: { ...doc.root, children: doc.root.children?.filter((c) => c.ref !== 'v') },
    });
    expect(s.get(r)).toBe(null);
    t.ctx.dispose();
    expect(() => s.pass()).not.toThrow();
  });

  it('the pyramid reductions equal a brute force at every level, gappy samples fading', () => {
    const C = 300;
    const N = 70;
    const codes = new Uint8Array(N * C);
    for (let i = 0; i < codes.length; i++)
      codes[i] = i >= 60 * C ? 45 : [45, 45, 45, 46, 65, 67, 68, 69][((i * 7919) % 13) % 8];
    const slot = Int32Array.from({ length: N }, (_, k) => (k * 17) % N);
    const out = blank();
    const m = mismatches();
    for (const [ly, lx] of [
      [0, 0],
      [1, 2],
      [3, 1],
      [6, 8],
    ]) {
      reduce(codes, slot, C, N, ly, lx, 0, 0, 0, TILE, out);
      for (let y = 0; y * 2 ** ly < N && y < TILE; y++)
        for (let x = 0; x * 2 ** lx < C && x < TILE; x++)
          m.check([out.cls[y * TILE + x], out.share[y * TILE + x]], brute(codes, slot, C, N, ly, lx, y, x), [
            ly,
            lx,
            y,
            x,
          ]);
    }
    m.none();
    reduce(codes, slot, C, N, 2, 2, 0, 0, 0, TILE, out);
    expect([out.cls[0] !== 45 && out.cls[0] !== 46, out.share[0] < 128]).toEqual([true, true]);
    const gapRow = Array.from(slot).indexOf(64);
    reduce(codes, slot, C, N, 0, 3, 0, 0, 0, TILE, out);
    expect([out.cls[gapRow * TILE], out.share[gapRow * TILE]]).toEqual([45, 255]);
  });
});

describe('LY4: levels, local compute, the wire format', () => {
  it('a numeric budget coarsens the axis with more samples first, then the other', () => {
    const lv = new Int32Array(4);
    const levels = (cellPx: number, rowPx: number, cols: number, rows: number, budget: number) => {
      tileLevels(Float64Array.of(cellPx, rowPx, cols, rows, budget), lv);
      return Array.from(lv);
    };
    // A sample per device px: 1000 × 100 at 1 px each, no budget.
    expect(levels(1, 1, 1000, 100, Number.POSITIVE_INFINITY)).toEqual([0, 0, 0, 2]);
    // 100,000 samples allowed: x (1000) gives first, to 500 × 100.
    expect(levels(1, 1, 1000, 100, 50_000)).toEqual([0, 1, 0, 2]);
    // Tall: y gives first.
    expect(levels(1, 1, 100, 1000, 50_000)).toEqual([1, 0, 2, 0]);
    // Both, alternating, x first on a tie: 1000 × 1000 under 10,000 samples is 63 × 125 (lx 4, ly 3).
    expect(levels(1, 1, 1000, 1000, 10_000)).toEqual([3, 4, 3, 4]);
    // Below a device px: the finest level whose sample spans one (0.3 px: 4 cells), coarse holding all in a tile.
    expect(levels(0.3, 0.5, 3000, 600, Number.POSITIVE_INFINITY)).toEqual([1, 2, 2, 4]);
  });

  it('reduce past tile (0, 0) and over lanes not held (slot -1)', () => {
    const C = 700;
    const N = 600;
    const codes = new Uint8Array(N * C);
    for (let i = 0; i < codes.length; i++) codes[i] = [45, 65, 67, 68, 69, 70][(i * 7919) % 6];
    const slot = Int32Array.from({ length: N }, (_, k) => (k % 50 === 7 ? -1 : (k * 13) % N));
    const out = blank();
    const m = mismatches();
    for (const [ly, lx, ty, tx] of [
      [1, 0, 1, 0],
      [0, 1, 0, 1],
      [0, 0, 2, 2],
      [1, 1, 1, 1],
    ]) {
      reduce(codes, slot, C, N, ly, lx, ty, tx, 0, TILE, out);
      for (let y = 0; y < TILE; y += 7)
        for (let x = 0; x < TILE; x += 5)
          m.check(
            [out.cls[y * TILE + x], out.share[y * TILE + x]],
            brute(codes, slot, C, N, ly, lx, ty * TILE + y, tx * TILE + x),
            [ly, lx, ty, tx, y, x],
          );
    }
    m.none();
    // A lane not held at level 0: no data.
    reduce(codes, slot, C, N, 0, 0, 0, 0, 0, TILE, out);
    expect([out.cls[7 * TILE], out.share[7 * TILE]]).toEqual([0, 0]);
  });

  it('each cell is read once per sample', () => {
    const C = 300;
    const N = 200;
    const raw = Uint8Array.from({ length: N * C }, (_, i) => [45, 65, 67][i % 3]);
    let reads = 0;
    const codes = new Proxy(raw, {
      get(target, p) {
        if (typeof p === 'string' && /^\d+$/.test(p)) reads++;
        return Reflect.get(target, p);
      },
    }) as Uint8Array;
    const slot = Int32Array.from({ length: N }, (_, k) => k);
    const cells = reduce(codes, slot, C, N, 3, 2, 0, 0, 0, TILE, blank());
    expect([cells, reads]).toEqual([N * C, N * C]);
  });

  it("the reference server's tiles equal `reduce` over the same rows, `.` read as a gap", async () => {
    const [N, C] = [700, 600];
    const rows = Array.from({ length: N }, (_, k) => proceduralRow(k, C).replace(/-(?=A)/g, '.'));
    const server = tileServer({ rowCount: N, columnCount: C, row: (k) => rows[k] });
    const codes = new Uint8Array(N * C);
    for (let k = 0; k < N; k++) for (let c = 0; c < C; c++) codes[k * C + c] = rows[k].charCodeAt(c);
    const slot = Int32Array.from({ length: N }, (_, k) => k);
    const out = blank();
    for (const [ly, lx, ty, tx] of [
      [0, 0, 1, 2],
      [1, 1, 1, 0],
      [2, 3, 0, 0],
    ]) {
      const b = await server.transport(`tiles/${ly}/${lx}/${ty}/${tx}`);
      reduce(codes, slot, C, N, ly, lx, ty, tx, 0, TILE, out);
      // Level 0 keeps `.` where the cell is one; above it a sample's gap class is `-`.
      const cls = Array.from(out.cls, (c) => (c === 46 ? 45 : c));
      expect([Array.from(b.subarray(0, TILE * TILE)), Array.from(b.subarray(TILE * TILE))]).toEqual([
        cls,
        Array.from(out.share),
      ]);
    }
  });

  it('an alignment in memory past one tile per axis: every tile computed in place, as `reduce`', async () => {
    const t = await opened(largeDoc(600, 700), { ...P, size: { 'large.main': [400, 300] } });
    const s = sm(t);
    const vp = t.vp('large.main');
    vp.zoomTo(0, 0, 700);
    vp.setRowScale(0.02, 0);
    t.layout('large.main');
    s.pass();
    const st = layer(t.ctx, 'Rows').state as RowsState;
    const view = s.get(ref(t));
    expect(st.source.local).not.toBe(null);
    expect(s.stats.current.computed).toBeGreaterThan(2);
    let checked = 0;
    for (let ly = 0; ly < 16; ly++)
      for (let lx = 0; lx < 16; lx++)
        for (let ty = 0; ty < 3; ty++)
          for (let tx = 0; tx < 3; tx++) {
            const tile = view?.tile(ly, lx, ty, tx);
            if (!tile) continue;
            const want = blank();
            const { codes, slot } = st.source.local as { codes: Uint8Array; slot: Int32Array };
            reduce(codes, slot, st.columns, st.lanes, ly, lx, ty, tx, 0, TILE, want);
            expect(tile).toEqual(want);
            checked++;
          }
    expect(checked).toBe(s.stats.current.computed);
    expect([view?.state(0, 0, 0, 0), view?.state(0, 0, 0, 9)]).toEqual([RegionState.pending, RegionState.empty]);
  }, 30_000);

  it("a frame's local work is bounded: one tile over a 1,000 × 3,000 alignment takes frames, then lands", async () => {
    const t = await opened(largeDoc(1000, 3000), { ...P, size: { 'large.main': [600, 400] } });
    const s = sm(t);
    const vp = t.vp('large.main');
    vp.zoomTo(0, 0, 3000);
    vp.setRowScale(0.01, 0);
    t.layout('large.main');
    s.pass();
    // The coarse tile holds every cell (3M): more than a frame's work.
    expect(s.stats.current.computed).toBe(0);
    let frames = 1;
    while (s.stats.current.computed === 0 && frames < 10) {
      t.scheduler.advance(16);
      frames++;
    }
    expect([s.stats.current.computed > 0, frames >= 2, frames <= 3]).toEqual([true, true, true]);
  }, 30_000);

  it('an alignment upserted with the same rows and shape but new cells drops its regions', async () => {
    const t = await opened(largeDoc(300, 600), { ...P, size: { 'large.main': [400, 300] } });
    const s = sm(t);
    const vp = t.vp('large.main');
    vp.zoomTo(0, 0, 600);
    vp.setRowScale(0.1, 0);
    t.layout('large.main');
    s.pass();
    const r = ref(t);
    const before = s.get(r);
    const key = before?.source.key;
    expect(s.stats.current.computed).toBeGreaterThan(0);
    await t.ctx.setDocument(largeDoc(300, 600, 1));
    t.layout('large.main');
    s.pass();
    const after = s.get(r);
    const st = layer(t.ctx, 'Rows').state as RowsState;
    expect([after === before, after?.source.key === key]).toEqual([false, false]);
    // What it holds is the new cells' reduction.
    const want = blank();
    const { codes, slot } = st.source.local as { codes: Uint8Array; slot: Int32Array };
    for (let ly = 0; ly < 16; ly++)
      for (let lx = 0; lx < 16; lx++) {
        const tile = after?.tile(ly, lx, 0, 0);
        if (!tile) continue;
        reduce(codes, slot, st.columns, st.lanes, ly, lx, 0, 0, 0, TILE, want);
        expect(tile).toEqual(want);
      }
    // A new column count alone changes the key too.
    await t.ctx.setDocument(largeDoc(300, 601, 1));
    expect(s.get(r)?.source.key === after?.source.key).toBe(false);
  }, 30_000);
});

describe('LY4: failures, cancellation, eviction, filters', () => {
  /** A transport that fails the first `n` asks of each path matching `which` (NotServed with `notServed`). */
  function failing(rows: number, columns: number, n: number, which: RegExp, notServed = false) {
    const server = tileServer({ rowCount: rows, columnCount: columns });
    const seen = new Map<string, number>();
    const transport: TileTransport = (path, signal) => {
      const k = seen.get(path) ?? 0;
      seen.set(path, k + 1);
      if (which.test(path) && k < n)
        return Promise.reject(notServed ? new NotServed(`404 ${path}`) : new Error('down'));
      return server.transport(path, signal);
    };
    return gated(rows, columns, transport);
  }
  async function watched(loaders: Loaders) {
    const t = await opened(lazyDoc(2000, 400), { ...P, loaders, size: { 'v.main': [800, 400] } });
    const issues: Issue[] = [];
    t.ctx.events.issue.subscribe((i) => issues.push(i));
    return { t, issues, s: sm(t) };
  }
  // A streamed region's failure: retried, the rest stands, a warning under its own code (decision 9).
  const loads = (issues: Issue[]) =>
    issues.filter((i) => i.code === 'svs.stream.load' && i.severity === 'warning').length;

  it("`state` at a lane reads that lane's tile row: a failed lower tile under a loaded upper one (VT-08)", async () => {
    const { loaders } = failing(2000, 400, 99, /^tiles\/0\/0\/1\/0$/);
    const t = await opened(lazyDoc(2000, 400), { ...P, loaders, size: { 'v.main': [800, 9000] } });
    const s = sm(t);
    for (let k = 0; k < 4; k++) {
      s.pass();
      await flush();
    }
    expect([s.state(ref(t), 0, 10, 10), s.state(ref(t), 0, 10, 300)]).toEqual([RegionState.loaded, RegionState.failed]);
    t.ctx.dispose();
  });

  it("a click on the failed lower tile's band retries it; a click on the loaded upper band does not (LB-23)", async () => {
    const { driven } = await import('./input-fixtures.js');
    const { g, loaders } = failing(2000, 400, 99, /^tiles\/0\/0\/1\/0$/);
    const t = await driven(lazyDoc(2000, 400), { ...P, loaders, size: { 'v.main': [800, 9000] } });
    const s = t.ctx.streams as Streams;
    for (let k = 0; k < 4; k++) {
      s.pass();
      await flush();
    }
    const view = t.ctx.graph.assembled.value?.view;
    const track = view?.set.strings.get(view.resolved.tracks.ref[0]) as string;
    const asked = () => g.calls.filter((p) => p === 'tiles/0/0/1/0').length;
    const before = asked();
    const at = (row: number) => [t.xOf('v.main', 0, 10), t.yOf('v.main', track, row)] as const;
    // The upper band (row 10, loaded): nothing asked again.
    t.click('v.main', ...at(10));
    t.advance(400);
    s.pass();
    await flush();
    expect(asked()).toBe(before);
    // The lower band (row 300, its tile failed): asked again.
    t.click('v.main', ...at(300));
    t.advance(400);
    s.pass();
    await flush();
    expect(asked()).toBe(before + 1);
    t.ctx.dispose();
  });

  it('four failures: retries at 1, 2 and 4 s, a 30 s cooldown, one report; a pan after it loads the tile', async () => {
    const { g, loaders } = failing(2000, 400, 4, /^tiles\/0\/0\/0\/0$/);
    const { t, issues, s } = await watched(loaders);
    const asked = () => g.calls.filter((p) => p === 'tiles/0/0/0/0').length;
    s.pass();
    await flush();
    expect([asked(), loads(issues)]).toEqual([1, 1]);
    for (const [ms, n] of [
      [1000, 2],
      [2000, 3],
      [4000, 4],
    ]) {
      t.scheduler.advance(ms - 1);
      await flush();
      expect(asked()).toBe(n - 1);
      t.scheduler.advance(1);
      await flush();
      expect(asked()).toBe(n);
    }
    const view = s.get(ref(t));
    expect(view?.state(0, 0, 0, 0)).toBe(RegionState.failed);
    // Cooling down: a pan does not ask; after 30 s it does, and the tile loads.
    t.vp('v.main').panBy(0, 1, 'api');
    frames(t);
    await flush();
    expect(asked()).toBe(4);
    t.scheduler.advance(30_000);
    t.vp('v.main').panBy(0, -1, 'api');
    frames(t);
    await flush();
    expect([asked(), view?.state(0, 0, 0, 0), loads(issues), t.scheduler.pendingTimers]).toEqual([
      5,
      RegionState.loaded,
      1,
      0,
    ]);
  });

  it('Retry ends a cooldown at once', async () => {
    const { g, loaders } = failing(2000, 400, 4, /^tiles\/0\/0\/0\/0$/);
    const { t, s } = await watched(loaders);
    s.pass();
    await flush();
    for (const ms of [1000, 2000, 4000]) {
      t.scheduler.advance(ms);
      await flush();
    }
    expect(s.get(ref(t))?.state(0, 0, 0, 0)).toBe(RegionState.failed);
    expect(s.retry(ref(t))).toBe(1);
    frames(t);
    await flush();
    expect([g.calls.filter((p) => p === 'tiles/0/0/0/0').length, s.get(ref(t))?.state(0, 0, 0, 0)]).toEqual([
      5,
      RegionState.loaded,
    ]);
  });

  it('a label chunk retries alike; NotServed is reported once and never retried, Retry included', async () => {
    const a = failing(2000, 400, 1, /^labels\//);
    const x = await watched(a.loaders);
    x.s.pass();
    await flush();
    x.t.scheduler.advance(1000);
    await flush();
    expect([a.g.calls.filter((p) => p.startsWith('labels/0/')).length, loads(x.issues)]).toEqual([2, 1]);
    expect(x.s.get(ref(x.t))?.label(3)).toBe(`SEQ3/1-${proceduralRow(3, 400).replace(/-/g, '').length}`);
    const b = failing(2000, 400, 99, /^tiles\/0\/0\/0\/0$/, true);
    const y = await watched(b.loaders);
    y.s.pass();
    await flush();
    y.t.scheduler.advance(60_000);
    y.t.vp('v.main').panBy(0, 1, 'api');
    frames(y.t);
    await flush();
    expect([b.g.calls.filter((p) => p === 'tiles/0/0/0/0').length, loads(y.issues)]).toEqual([1, 1]);
    expect(y.s.retry(ref(y.t))).toBe(0);
  });

  it('dispose during a backoff leaves no timer and no frame', async () => {
    const { loaders } = failing(2000, 400, 9, /^tiles\//);
    const { t, s } = await watched(loaders);
    s.pass();
    await flush();
    expect(t.scheduler.pendingTimers).toBeGreaterThan(0);
    t.ctx.dispose();
    expect([t.scheduler.pendingTimers, t.scheduler.pendingFrames]).toEqual([0, 0]);
  });

  it("eviction by bytes: the view's regions pinned past the bound, the least recently needed dropped first", async () => {
    const { g, loaders } = gated(40_000, 3000);
    const t = await opened(lazyDoc(40_000, 3000), { ...P, loaders, size: { 'v.main': [1000, 600] } });
    const tile = 2 * TILE * TILE;
    t.ctx.config.set(SeqViewerConfig.Stream.CacheBytes, tile);
    const s = sm(t);
    const vp = t.vp('v.main');
    vp.zoomTo(0, 0, 600);
    s.pass();
    await flush();
    const view = s.get(ref(t));
    const held = s.stats.current.bytes / tile;
    // More than the bound, all the view's.
    expect(held).toBeGreaterThan(1);
    const first = tiles(g)
      .filter((p) => p.startsWith('tiles/0/0/'))
      .map((p) => p.split('/').slice(3).map(Number));
    vp.scrollTo(2000 * vp.live.current.cell, 20_000 * 22, 'api');
    s.pass();
    await flush();
    // The first view's fine tiles are gone; the new view's all held.
    for (const [ty, tx] of first) expect(view?.state(0, 0, ty, tx)).toBe(RegionState.pending);
    expect(s.stats.current.bytes / tile).toBeLessThanOrEqual(held + 2);
    expect(view?.tile(0, 0, Math.floor(20_000 / TILE), Math.floor(2000 / TILE))).not.toBe(null);
  });

  it("a filtered million-column alignment asks for the kept runs' tiles only, near the view first", async () => {
    const { g, loaders } = gated(100, 1_000_000);
    const doc = lazyDoc(100, 1_000_000, {
      filter: [
        [100, 200],
        [900_000, 900_100],
      ],
    });
    const t = await opened(doc, { ...P, loaders, size: { 'v.main': [1000, 400] } });
    const s = sm(t);
    t.vp('v.main').zoomTo(0, 100, 900_100);
    t.layout('v.main');
    s.pass();
    await flush();
    const asked = tiles(g);
    expect(asked.length).toBeGreaterThan(0);
    expect(asked.length).toBeLessThanOrEqual(8);
    for (const p of asked) {
      const [, , lx, , tx] = p.split('/').map(Number);
      const c0 = tx * TILE * 2 ** lx;
      const c1 = c0 + TILE * 2 ** lx;
      expect((c0 < 200 && c1 > 100) || (c0 < 900_100 && c1 > 900_000)).toBe(true);
    }
  });
});
