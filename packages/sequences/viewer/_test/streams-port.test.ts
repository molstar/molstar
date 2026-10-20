/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1: the stream manager on the port (seqspec:sources#the-viewers-side; seqspec:sources#levels-and-the-zoom-gate),
// headless, over a streamed genome's values (a scripted provider: answers held until released, every request and abort
// recorded). The view's regions come before
// the margin's and coarse before fine; a jump aborts what left the view; a held region is never asked again and a pass
// over a held view changes nothing; eviction under a small `Stream.CacheBytes` keeps the view; pending, loaded and empty
// are told apart; `settled` waits for the view, failures included; a thousand patches leave the document, its history
// and the set untouched; and the values renderer draws the regions (`renderSvg` awaiting them).

import { describe, expect, it } from 'vitest';
import { GRCH38, hg38Streams } from '../../core/_test/format/stream-provider.js';
import { RegionState, SeqViewerConfig } from '../src/index.js';
import { renderSvg } from '../src/svg/index.js';
import {
  CacheBytes,
  GENOME,
  Margin,
  MaxInFlight,
  openGenome,
  type StreamedGenome,
  streamedDoc,
} from './stream-fixtures.js';

const flush = () => new Promise((ok) => setTimeout(ok, 0));
/** Three chromosomes of GRCh38's lengths: the manager's logic at full size, fewer blocks to read. */
const THREE = GRCH38.slice(0, 3);
const CHR2 = 1;

/** The window [start, stop) of block b in view, then one pass. */
function show(t: StreamedGenome, b: number, start: number, stop: number) {
  t.vp(GENOME).zoomTo(b, start, stop);
  t.layout(GENOME);
  t.manager.pass();
}
async function drain(t: StreamedGenome) {
  for (let k = 0; k < 6; k++) {
    t.streams.phylop.release();
    await flush();
    t.manager.pass();
  }
}

describe('requests: the view, then the margin; coarse before fine', () => {
  it('at 10 kb: the coarsest level first, the data in view from the center out, then the margin around it', async () => {
    const t = await openGenome({ gated: true, chromosomes: THREE });
    t.ctx.config.set(MaxInFlight, 100);
    t.streams.phylop.requests.length = 0;
    show(t, CHR2, 100_000_000, 100_010_000);
    const asked = t.streams.phylop.requests.map((r) => [r.system, r.level, r.start, r.stop, r.priority] as const);
    expect(asked[0]).toEqual(['chr2', 5, 0, 242_193_529, 'view']);
    const view = asked.filter((r) => r[1] === -1 && r[4] === 'view');
    const margin = asked.filter((r) => r[4] === 'margin');
    expect(view.map((r) => r[2])).toEqual([24_415 * 4096, 24_414 * 4096, 24_416 * 4096]);
    expect(asked.slice(1, 1 + view.length)).toEqual(view);
    // A view width on each side: [lo - w, lo) and [hi, hi + w), the regions not already asked for the view.
    const v = t.ctx.layout.axis(GENOME);
    const [lo, hi] = [v?.visLo[CHR2] ?? 0, v?.visHi[CHR2] ?? 0];
    const w = hi - lo;
    const cover = (a: number, b: number) =>
      Array.from(
        { length: Math.floor((b - 1) / 4096) - Math.floor(a / 4096) + 1 },
        (_, k) => (Math.floor(a / 4096) + k) * 4096,
      );
    const inView = new Set(cover(lo, hi));
    const want = [...cover(lo - w, lo), ...cover(hi, hi + w)].filter((x) => !inView.has(x));
    expect(margin.map((r) => r[2]).sort()).toEqual(want.sort());
    expect(asked.length).toBe(1 + view.length + margin.length);
    t.ctx.dispose();
  });

  it("no margin when `Stream.Margin` is 0; at most `Stream.MaxInFlight` in flight, the view's first", async () => {
    const t = await openGenome({ gated: true, chromosomes: THREE });
    t.ctx.config.set(Margin, 0);
    t.streams.phylop.requests.length = 0;
    show(t, CHR2, 100_000_000, 100_010_000);
    expect(t.streams.phylop.requests.map((r) => r.priority)).toEqual(['view', 'view', 'view', 'view']);
    t.ctx.config.set(Margin, 1);
    t.ctx.config.set(MaxInFlight, 2);
    show(t, CHR2, 150_000_000, 150_010_000);
    expect(t.manager.stats.current.inFlight).toBe(2);
    t.ctx.dispose();
  });

  it("the whole genome: one coarse region a chromosome, every chromosome's level its own", async () => {
    const t = await openGenome({ gated: true });
    t.ctx.config.set(MaxInFlight, 100);
    t.streams.phylop.requests.length = 0;
    const vp = t.vp(GENOME);
    vp.fit();
    t.layout(GENOME);
    t.manager.pass();
    const asked = t.streams.phylop.requests;
    expect(new Set(asked.map((r) => r.system)).size).toBe(24);
    expect(asked.every((r) => r.level === 5 && r.priority === 'view')).toBe(true);
    expect(asked.length).toBe(24);
    t.ctx.dispose();
  });
});

describe('the cache', () => {
  it('a jump aborts what left the view; a held region is not asked again and a pass over it changes nothing', async () => {
    const t = await openGenome({ gated: true, chromosomes: THREE });
    t.ctx.config.set(MaxInFlight, 100);
    show(t, CHR2, 100_000_000, 100_010_000);
    const first = t.streams.phylop.requests.slice();
    show(t, 2, 50_000_000, 50_010_000);
    // Kept for `Stream.CancelGraceMs` after no pass needs them (VT-02), then aborted.
    expect(t.streams.phylop.aborted.length).toBe(0);
    t.scheduler.advance(t.ctx.config.get(SeqViewerConfig.Stream.CancelGraceMs));
    t.manager.pass();
    expect(t.streams.phylop.aborted.length).toBe(first.length);
    expect(t.manager.stats.current.cancelled).toBe(first.length);
    await drain(t);
    const n = t.streams.phylop.requests.length;
    const version = t.manager.version.value;
    const view = t.manager.get(t.layer)?.version;
    for (let k = 0; k < 3; k++) t.manager.pass();
    expect(t.streams.phylop.requests.length).toBe(n);
    expect([t.manager.version.value, t.manager.get(t.layer)?.version]).toEqual([version, view]);
    expect(t.manager.stats.current.hits).toBeGreaterThan(0);
    t.ctx.dispose();
  });

  it("past a small `Stream.CacheBytes` the regions out of view go, the view's kept", async () => {
    const t = await openGenome({ chromosomes: THREE });
    t.ctx.config.set(CacheBytes, 200_000);
    // The bound is for what the views pin: a prefetch's regions would pin more.
    t.ctx.config.set(SeqViewerConfig.Stream.Prefetch, false);
    for (const at of [10, 60, 110, 160]) {
      show(t, CHR2, at * 1e6, at * 1e6 + 10_000);
      await drain(t);
    }
    expect(t.manager.stats.current.bytes).toBeLessThanOrEqual(200_000 + 3 * 32_768 + 6_144);
    for (const i of [160_000_000, 160_005_000, 160_009_999])
      expect(t.manager.state(t.layer, CHR2, i)).toBe(RegionState.loaded);
    t.ctx.config.set(CacheBytes, 64 * 2 ** 20);
    show(t, CHR2, 10_000_000, 10_010_000);
    expect(t.manager.state(t.layer, CHR2, 10_005_000)).toBe(RegionState.pending);
    t.ctx.dispose();
  });

  it('pending, loaded and empty told apart: off the view pending, outside the source empty, no values empty', async () => {
    const t = await openGenome({
      gated: true,
      chromosomes: THREE,
      answer: (r) => (r.system === 'chr3' ? { ...r, kind: 'items', items: { value: [] } } : null),
    });
    show(t, CHR2, 100_000_000, 100_010_000);
    expect(t.manager.state(t.layer, CHR2, 100_005_000)).toBe(RegionState.pending);
    expect(t.manager.state(t.layer, CHR2, 242_193_529)).toBe(RegionState.empty);
    expect(t.manager.state(t.layer, 2, 5)).toBe(RegionState.pending);
    expect(t.manager.state(t.layer, 7, 5)).toBe(RegionState.empty);
    await drain(t);
    expect(t.manager.state(t.layer, CHR2, 100_005_000)).toBe(RegionState.loaded);
    show(t, 2, 1_000_000, 1_010_000);
    await drain(t);
    expect(t.manager.state(t.layer, 2, 1_005_000)).toBe(RegionState.empty);
    t.ctx.dispose();
  });

  /** Views visited until `n` patches have landed; the document, its history and the set checked unchanged. */
  async function patches(n: number) {
    const t = await openGenome({ chromosomes: THREE });
    const before = {
      base: t.ctx.document.base.value,
      shown: t.ctx.document.shown.value,
      local: t.ctx.document.local.value,
      revision: t.ctx.document.revision.value,
      set: t.ctx.graph.assembled.value?.set,
    };
    for (let k = 0; t.manager.stats.current.coarse + t.manager.stats.current.fine < n; k++) {
      show(t, k % 3, (k * 1_000_003) % 40_000_000, ((k * 1_000_003) % 40_000_000) + 10_000);
      await drain(t);
    }
    expect({
      base: t.ctx.document.base.value,
      shown: t.ctx.document.shown.value,
      local: t.ctx.document.local.value,
      revision: t.ctx.document.revision.value,
      set: t.ctx.graph.assembled.value?.set,
    }).toEqual(before);
    expect(t.ctx.document.shown.value).toBe(before.shown);
    expect(t.ctx.graph.assembled.value?.set).toBe(before.set);
    expect(t.ctx.history.canUndo.value).toBe(false);
    t.ctx.dispose();
  }

  it('two hundred patches leave the document, its history and the set untouched', () => patches(200));

  it('a thousand patches leave the document, its history and the set untouched', { tags: ['heavy'] }, () =>
    patches(1000),
  );
});

describe('`settled`', () => {
  it('waits for every region in view at its level, and for failures too (one report an outage)', async () => {
    const t = await openGenome({ gated: true, chromosomes: THREE });
    show(t, CHR2, 100_000_000, 100_010_000);
    let done = false;
    const settled = t.manager.settled().then(() => {
      done = true;
    });
    for (let k = 0; k < 3; k++) {
      t.scheduler.advance(16);
      await flush();
    }
    expect(done).toBe(false);
    t.streams.phylop.release();
    for (let k = 0; k < 6 && !done; k++) {
      t.scheduler.advance(16);
      await flush();
      t.streams.phylop.release();
    }
    await settled;
    expect(t.manager.state(t.layer, CHR2, 100_005_000)).toBe(RegionState.loaded);
    t.ctx.dispose();

    const issues: string[] = [];
    const f = await openGenome({ chromosomes: THREE, answer: () => new Error('offline') });

    show(f, CHR2, 100_000_000, 100_010_000);
    let failed = false;
    const p = f.manager.settled().then(() => {
      failed = true;
    });
    for (let k = 0; k < 6 && !failed; k++) {
      f.scheduler.advance(16);
      await flush();
    }
    await p;
    expect(f.manager.state(f.layer, CHR2, 100_005_000)).toBe(RegionState.failed);
    expect(f.manager.hasFailed(f.layer)).toBe(true);
    issues.push(...f.ctx.issues.current.value.map((i) => i.code));
    expect(issues).toEqual(['svs.stream.load']);
    f.ctx.dispose();
  });
});

describe('a patch the region cannot take (GI-13)', () => {
  it('fails the region, reported and retried as a failure, never an unhandled rejection that leaves it pending', async () => {
    const f = await openGenome({
      chromosomes: THREE,
      // A provider's answer of another shape: a bins patch with no bins.
      answer: (r) => ({ ...r, kind: 'bins' }) as never,
    });
    show(f, CHR2, 100_000_000, 100_010_000);
    let done = false;
    const p = f.manager.settled().then(() => {
      done = true;
    });
    for (let k = 0; k < 6 && !done; k++) {
      f.scheduler.advance(16);
      await flush();
    }
    await p;
    expect(f.manager.state(f.layer, CHR2, 100_005_000)).toBe(RegionState.failed);
    expect(f.ctx.issues.current.value.map((i) => i.code)).toEqual(['svs.stream.load']);
    f.ctx.dispose();
  });
});

describe('a patch the region cannot take backs off as a failure (LA-04)', () => {
  /** The requests over 40 s of scheduler time with every answer `answer`, and the issues reported. */
  async function over40s(answer: (r: { system: string | null }) => unknown) {
    let calls = 0;
    const f = await openGenome({
      chromosomes: THREE,
      answer: (r) => {
        calls++;
        return answer(r) as never;
      },
    });
    show(f, CHR2, 100_000_000, 100_010_000);
    for (let t = 0; t < 40_000; t += 250) {
      f.scheduler.advance(250);
      await flush();
    }
    const out = { calls, issues: f.ctx.issues.current.value.map((i) => i.code), failed: f.manager.hasFailed(f.layer) };
    f.ctx.dispose();
    return out;
  }

  it('is retried at 1, 2 and 4 s and then after the cooldown, as a rejection is, and its outage is not ended by it', async () => {
    const rejected = await over40s(() => new Error('offline'));
    const malformed = await over40s((r) => ({ ...r, kind: 'bins' }));
    // A few regions (the view and its margin), each asked about six times, not once a second.
    expect(rejected.calls).toBeLessThan(40);
    expect(malformed).toEqual(rejected);
    expect(malformed.issues).toEqual(['svs.stream.load']);
    expect(malformed.failed).toBe(true);
  }, 30_000);

  it('an answer of another shape than asked is a failure too: a values region answered with letters, rows or a manifest (LB-19)', async () => {
    const rejected = await over40s(() => new Error('offline'));
    for (const kind of ['letters', 'rows', 'manifest'])
      expect([
        kind,
        await over40s((r) => ({ ...r, kind, letters: 'ACGT', labels: [], rows: null, manifest: {} })),
      ]).toEqual([kind, rejected]);
  }, 60_000);

  it('a malformed answer during an outage neither ends it nor reports it again (LB-20)', async () => {
    let mode: 'reject' | 'malformed' = 'reject';
    const t = await openGenome({
      chromosomes: THREE,
      answer: (r) =>
        r.level !== -1 ? null : mode === 'reject' ? (new Error('503') as never) : ({ ...r, kind: 'bins' } as never),
    });
    t.ctx.config.set(Margin, 0);
    const loads = () => t.ctx.issues.current.value.filter((i) => i.code === 'svs.stream.load').length;
    const at = async (start: number) => {
      show(t, CHR2, start, start + 10_000);
      for (let k = 0; k < 3; k++) {
        t.manager.pass();
        await flush();
      }
    };
    await at(10_000_000);
    expect(loads()).toBe(1);
    const failedAt = (start: number) => t.manager.state(t.layer, CHR2, start + 5_000);
    expect(failedAt(10_000_000)).toBe(RegionState.failed);
    // The view moves on and its answers are malformed: one outage still, the first view's failed regions kept.
    mode = 'malformed';
    await at(100_000_000);
    expect(loads()).toBe(1);
    expect(failedAt(10_000_000)).toBe(RegionState.failed);
    t.ctx.dispose();
  });
});

describe('a provider whose fetch throws (VT-05)', () => {
  it('fails the region as a rejection would: the pass completes, the slot is freed, it loads once the provider recovers', async () => {
    const { opened, preparedRenderers } = await import('./layout-fixtures.js');
    const streams = hg38Streams({ chromosomes: THREE });
    let throwing = true;
    const loaders = {
      load: (r: { uri?: string }) => Promise.reject(new Error(`no file ${r.uri}`)),
      stream: (src: Parameters<typeof streams.stream>[0]) => {
        const p = streams.stream(src);
        if (!p || src.stream !== 'phylop') return p;
        return {
          fetch: (r: Parameters<typeof p.fetch>[0], signal: AbortSignal) => {
            if (throwing && r.system !== null) throw new Error('bad request');
            return p.fetch(r, signal);
          },
        };
      },
    };
    const t = await opened(streamedDoc(), {
      renderers: preparedRenderers(),
      loaders: loaders as never,
      size: { [GENOME]: [1000, 200] },
    });
    const s = t.ctx.streams as unknown as StreamedGenome['manager'];
    t.vp(GENOME).zoomTo(CHR2, 100_000_000, 100_010_000);
    t.layout(GENOME);
    expect(() => s.pass()).not.toThrow();
    await flush();
    expect(s.stats.current.inFlight).toBe(0);
    throwing = false;
    const m = t.ctx.graph.assembled.value?.view;
    let layer = '';
    for (let k = 0; m && k < m.resolved.layers.count; k++)
      if (m.set.strings.get(m.resolved.layers.type[k]) === 'values')
        layer = m.set.strings.get(m.resolved.layers.ref[k]);
    expect(s.state(layer, CHR2, 100_005_000)).toBe(RegionState.failed);
    expect(t.ctx.issues.current.value.map((i) => i.code)).toContain('svs.stream.load');
    // Its backoff passes, then it is asked again and loads.
    for (let k = 0; k < 8 && s.state(layer, CHR2, 100_005_000) !== RegionState.loaded; k++) {
      t.scheduler.advance(2000);
      s.pass();
      await flush();
    }
    expect(s.state(layer, CHR2, 100_005_000)).toBe(RegionState.loaded);
    t.ctx.dispose();
  });

  it('a fetch that returns a plain patch or nothing, not a promise: taken, or failed, the slot freed (LB-21)', async () => {
    const { opened, preparedRenderers } = await import('./layout-fixtures.js');
    for (const plain of [true, false]) {
      const streams = hg38Streams({ chromosomes: THREE });
      const loaders = {
        load: (r: { uri?: string }) => Promise.reject(new Error(`no file ${r.uri}`)),
        stream: (src: Parameters<typeof streams.stream>[0]) => {
          const p = streams.stream(src);
          if (!p || src.stream !== 'phylop') return p;
          // A JS provider that is not `async`: the patch the scripted one would give, at once, or nothing.
          const patch = (r: Parameters<typeof p.fetch>[0]) => ({
            ...r,
            source: r.source,
            kind: 'bins',
            bins: { mean: [1, 2, 3] },
            row_start: 0,
            row_stop: 0,
          });
          return {
            fetch: (r: Parameters<typeof p.fetch>[0], signal: AbortSignal) =>
              r.system === null ? p.fetch(r, signal) : plain ? patch(r) : undefined,
          };
        },
      };
      const t = await opened(streamedDoc(), {
        renderers: preparedRenderers(),
        loaders: loaders as never,
        size: { [GENOME]: [1000, 200] },
      });
      const s = t.ctx.streams as unknown as StreamedGenome['manager'];
      t.vp(GENOME).zoomTo(CHR2, 100_000_000, 100_010_000);
      t.layout(GENOME);
      expect(() => s.pass()).not.toThrow();
      await flush();
      expect([plain, s.stats.current.inFlight]).toEqual([plain, 0]);
      expect([plain, s.stats.current.errors > 0]).toEqual([plain, !plain]);
      t.ctx.dispose();
    }
  });
});

describe('an outage ends with a load (VT-04)', () => {
  it('the source back: a region loading ends the outage and drops the failed regions out of view; the next is reported', async () => {
    let down = true;
    const t = await openGenome({
      chromosomes: THREE,
      answer: (r) => (down && r.level === -1 ? (new Error('503') as never) : null),
    });
    t.ctx.config.set(Margin, 0);
    const loads = () => t.ctx.issues.current.value.filter((i) => i.code === 'svs.stream.load').length;
    const at = async (start: number) => {
      t.vp(GENOME).zoomTo(CHR2, start, start + 10_000);
      t.layout(GENOME);
      for (let k = 0; k < 3; k++) {
        t.manager.pass();
        await flush();
      }
    };
    await at(10_000_000); // A, while the source is down: one report
    expect([loads(), t.manager.hasFailed(t.layer)]).toEqual([1, true]);
    down = false;
    await at(100_000_000); // B, the source back: loaded, A's failures gone with the outage
    expect(t.manager.state(t.layer, CHR2, 100_005_000)).toBe(RegionState.loaded);
    expect(t.manager.hasFailed(t.layer)).toBe(false);
    down = true;
    await at(200_000_000); // C, down again: a second outage, reported
    expect(t.manager.state(t.layer, CHR2, 200_005_000)).toBe(RegionState.failed);
    expect(loads()).toBe(2);
    t.ctx.dispose();
  });
});

describe('`settled` and the margin (ST2)', () => {
  it("a region the margin asked first in a pass and a view asks after (another block's) is waited for", async () => {
    // Two blocks over one system can share a region: block b's margin can ask it before block b+1's view does. The
    // second ask used to return at once (asked this pass), so the view's region in flight was not counted.
    const t = await openGenome({ gated: true, chromosomes: THREE });
    show(t, CHR2, 100_000_000, 100_010_000);
    const m = t.manager as any;
    const c = m.order[0];
    const r = c.all.find((x: { flight: unknown; group: number }) => x.flight !== null && x.group === 4);
    expect(r).toBeDefined();
    m.clock++;
    m.missing = 0;
    m.want(0, c, r.level, r.ty, r.tx, r.index, 4);
    expect(m.missing).toBe(0);
    m.want(0, c, r.level, r.ty, r.tx, r.index, 1);
    expect(m.missing).toBe(1);
    m.want(0, c, r.level, r.ty, r.tx, r.index, 1);
    expect(m.missing).toBe(1);
    m.queued = 0;
    t.ctx.dispose();
  });
});

describe('drawn from the cache', () => {
  it('the values bars come from the regions: `renderSvg` of a document awaits them', async () => {
    const svg = async (empty: boolean) => {
      const s = hg38Streams({
        chromosomes: THREE,
        answer: empty ? (r) => ({ ...r, kind: 'items', items: { value: [] } }) : undefined,
      });
      const loaders = { load: () => Promise.reject(new Error('none')), stream: s.stream };
      const out = await renderSvg(streamedDoc(), { width: 800, loaders });
      return { out, asked: s.phylop.requests };
    };
    const full = await svg(false);
    const empty = await svg(true);
    expect(full.asked.length).toBeGreaterThan(0);
    expect(full.out.length).toBeGreaterThan(empty.out.length + 1000);
  });
});
