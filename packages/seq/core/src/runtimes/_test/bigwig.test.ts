/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { Asset, DefaultAssetManager, openRanges } from '@molstar/core/util/assets';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { delayedRuntime } from '../../data/async-source.js';
import {
  type AssetReader,
  DefaultOpenLimits,
  type OpenContext,
  type RegionAnswer,
  type ValuesBins,
  type ValuesPositions,
} from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';
import { headlessAssetPolicy, policyAssets } from '../../context/assets.js';
import { binValues } from '../../data/summaries.js';
import { FakeClock } from '../../data/_test/synthetic.js';
import { sourceSuite } from '../../data/_test/source-suite.js';
import { BigWigRuntime, bigWigRuntime } from '../bigwig.js';
import { syntheticSignal, type WriterChrom, type WriterOptions, writeBigWig } from './bigwig-writer.js';

// The bigWig runtime (_spec/README.md "Data", bigWig) on files its test-only writer makes from a seed: both zoom
// levels and level 1, a chromosome of each section type, data at both edges, a gap (an empty region) and a
// chromosome without data. Two hosts read them: core's asset manager by range reads over a fake `fetch` (nothing
// reaches the network), and a bare reader of whole bytes (a test's, a CLI's).

const Zooms = [16, 256];
const Chroms: WriterChrom[] = [
  { name: 'chrA', values: syntheticSignal(1, 10_000, [[2000, 2600]]) },
  { name: 'chrB', values: syntheticSignal(2, 3000, [[0, 100]]), type: 'varStep' },
  { name: 'chrC', values: syntheticSignal(3, 2500, [[1200, 1300]]), type: 'fixedStep' },
  { name: 'chrE', values: new Float32Array(500).fill(Number.NaN) },
];
const Options: WriterOptions = { zooms: Zooms, itemsPerSection: 64, blockSize: 4 };
const Fixture = await writeBigWig(Chroms, Options);
const Url = 'https://x.org/signal.bw';

const limits = DefaultOpenLimits;
const signal = new AbortController().signal;

/** A fake server of `files` by URL, honest about ranges; its requests logged. */
function serve(files: Record<string, Uint8Array>) {
  const log: { url: string; range: string }[] = [];
  vi.stubGlobal('fetch', async (input: RequestInfo | URL, init: RequestInit = {}) => {
    const url = String(input);
    const range = new Headers(init.headers).get('Range') ?? '';
    log.push({ url, range });
    const body = files[url];
    if (!body) return new Response('', { status: 404 });
    const m = /^bytes=(\d+)-(\d+)$/.exec(range);
    if (!m) return new Response(body.slice(), { status: 200 });
    const a = Number(m[1]);
    const b = Math.min(Number(m[2]), body.length - 1);
    if (a >= body.length)
      return new Response('', { status: 416, headers: { 'Content-Range': `bytes */${body.length}` } });
    return new Response(body.slice(a, b + 1), {
      status: 206,
      headers: { 'Content-Range': `bytes ${a}-${b}/${body.length}`, ETag: '"f"' },
    });
  });
  return log;
}

afterEach(() => vi.unstubAllGlobals());

/** Core's manager behind the port, by range reads (pages of 4 KiB: a small file still spans many). */
function rangeContext(manager = new DefaultAssetManager()): OpenContext {
  const assets: AssetReader = {
    readText: () => Promise.reject(new Error('not read as text')),
    readBytes: () => Promise.reject(new Error('read whole')),
    openRange: (uri, s) => openRanges(manager, Asset.Url(uri), { signal: s, pageSize: 4096 }).run(),
  };
  return { assets, signal, limits };
}

/** A bare reader of whole bytes: no range reads. */
function bytesContext(files: Record<string, Uint8Array>): OpenContext {
  return {
    assets: {
      readText: () => Promise.reject(new Error('not read as text')),
      readBytes: async (uri) => {
        if (!files[uri]) throw new Error('not found');
        return files[uri];
      },
    },
    signal,
    limits,
  };
}

const entry = { address: '/data/signal', format: 'bigwig', uri: Url };

async function values(answer: Promise<RegionAnswer>) {
  const a = await answer;
  if (!a.ok) throw new Error(a.issue.message);
  return a.tables.values!;
}

/** What the file holds over `[start, end)`: the generator's values (as float32). */
function expected(chrom: string, start: number, end: number) {
  return Chroms.find((c) => c.name === chrom)!.values.slice(start, end);
}

function expectBins(got: ValuesBins, want: ValuesBins) {
  expect(got.start).toBe(want.start);
  expect(Array.from(got.count)).toEqual(Array.from(want.count));
  for (const k of ['min', 'max', 'mean'] as const)
    for (let i = 0; i < want.count.length; i++) {
      if (Number.isNaN(want[k][i])) expect(got[k][i], `${k} ${i}`).toBeNaN();
      else expect(got[k][i], `${k} ${i}`).toBeCloseTo(want[k][i], 4);
    }
}

describe('the bigWig runtime', () => {
  it('passes the shared source suite through range reads and whole bytes', async () => {
    serve({ [Url]: Fixture });
    await sourceSuite(BigWigRuntime, entry, rangeContext());
    await sourceSuite(BigWigRuntime, entry, bytesContext({ [Url]: Fixture }));
  });

  it('covers the chromosomes with their lengths, values, and level 1 and both zoom levels', async () => {
    serve({ [Url]: Fixture });
    const s = await BigWigRuntime.open(entry, rangeContext());
    expect(s.coverage.systems).toEqual(Chroms.map((c) => ({ name: c.name, length: c.values.length })));
    expect(s.coverage.shapes.every((x) => x.shape === 'values')).toBe(true);
    expect(s.coverage.levels).toEqual([1, ...Zooms]);
    s.dispose();
  });

  it('level 1: every section type read back per position, the gaps missing, both edges held', async () => {
    serve({ [Url]: Fixture });
    const s = await BigWigRuntime.open(entry, rangeContext());
    for (const c of Chroms) {
      const n = c.values.length;
      for (const [a, b] of [
        [0, Math.min(n, 1024)],
        [Math.max(0, n - 1024), n],
        [Math.floor(n / 2) - 7, Math.floor(n / 2) + 900],
      ]) {
        const v = (await values(
          s.read({ system: c.name, start: a, end: b, level: 1, shapes: ['values'] }, signal),
        )) as ValuesPositions;
        expect(v.form).toBe('positions');
        expect(v.start).toBe(a);
        expect(Array.from(v.value)).toEqual(Array.from(expected(c.name, a, b)));
      }
    }
    s.dispose();
  });

  it('a seeded fuzz of regions: level 1 the generator’s values, a zoom level their bins', async () => {
    serve({ [Url]: Fixture });
    const s = await BigWigRuntime.open(entry, rangeContext());
    const rng = new PCG(11);
    for (let i = 0; i < 40; i++) {
      const c = Chroms[rng.int() % Chroms.length];
      const n = c.values.length;
      const a = rng.int() % n;
      const b = Math.min(n, a + 1 + (rng.int() % 1500));
      const v = (await values(
        s.read({ system: c.name, start: a, end: b, level: 1, shapes: ['values'] }, signal),
      )) as ValuesPositions;
      expect(Array.from(v.value), `${c.name} ${a}-${b}`).toEqual(Array.from(expected(c.name, a, b)));
      const L = Zooms[rng.int() % Zooms.length];
      const z = (await values(
        s.read({ system: c.name, start: a, end: b, level: L, shapes: ['values'] }, signal),
      )) as ValuesBins;
      const from = Math.floor(a / L) * L;
      const all = binValues({ form: 'positions', start: 0, value: c.values }, L, n);
      const first = from / L;
      const count = Math.ceil(b / L) - first;
      expectBins(z, {
        form: 'bins',
        start: from,
        min: all.min.subarray(first, first + count),
        max: all.max.subarray(first, first + count),
        mean: all.mean.subarray(first, first + count),
        count: all.count.subarray(first, first + count),
      });
    }
    s.dispose();
  });

  it('an empty region, an empty chromosome and a region past the end answer empty, not failed', async () => {
    serve({ [Url]: Fixture });
    const s = await BigWigRuntime.open(entry, rangeContext());
    const gap = (await values(
      s.read({ system: 'chrA', start: 2000, end: 2600, level: 1, shapes: ['values'] }, signal),
    )) as ValuesPositions;
    expect(gap.value.every(Number.isNaN)).toBe(true);
    const binned = (await values(
      s.read({ system: 'chrA', start: 2048, end: 2560, level: 256, shapes: ['values'] }, signal),
    )) as ValuesBins;
    expect(Array.from(binned.count)).toEqual([0, 0]);
    expect(binned.mean.every(Number.isNaN)).toBe(true);
    const none = (await values(
      s.read({ system: 'chrE', start: 0, end: 500, level: 16, shapes: ['values'] }, signal),
    )) as ValuesBins;
    expect(none.count.every((x) => x === 0)).toBe(true);
    const past = (await values(
      s.read({ system: 'chrB', start: 3000, end: 4000, level: 1, shapes: ['values'] }, signal),
    )) as ValuesPositions;
    expect(past.value.length).toBe(0);
    s.dispose();
  });

  it('a level it does not have, a shape other than values and an unknown chromosome are not served', async () => {
    serve({ [Url]: Fixture });
    const s = await BigWigRuntime.open(entry, rangeContext());
    for (const r of [
      { system: 'chrA', start: 0, end: 64, level: 4, shapes: ['values' as const] },
      { system: 'chrA', start: 0, end: 64, level: 1, shapes: ['ranges' as const] },
      { system: 'chrZ', start: 0, end: 64, level: 1, shapes: ['values' as const] },
    ])
      expect(await s.read(r, signal)).toMatchObject({ ok: false, kind: 'not-served' });
    s.dispose();
  });

  it('reads by range: a region at a zoom level fetches its index and blocks, not the file; a pan back fetches nothing', async () => {
    const log = serve({ [Url]: Fixture });
    const s = await BigWigRuntime.open(entry, rangeContext());
    const opened = log.length;
    await values(s.read({ system: 'chrC', start: 0, end: 1024, level: 1, shapes: ['values'] }, signal));
    const fetched = log
      .slice(opened)
      .map((l) => /bytes=(\d+)-(\d+)/.exec(l.range)!)
      .reduce((n, m) => n + Number(m[2]) + 1 - Number(m[1]), 0);
    expect(log.every((l) => l.range !== '')).toBe(true);
    expect(fetched).toBeLessThan(Fixture.length / 2);
    const n = log.length;
    await values(s.read({ system: 'chrC', start: 0, end: 1024, level: 1, shapes: ['values'] }, signal));
    expect(log.length).toBe(n);
    s.dispose();
  });

  it('under a real host’s policy (policyAssets over core’s manager) it reads by range, never the whole file', async () => {
    const log = serve({ [Url]: Fixture });
    const assets = policyAssets(new DefaultAssetManager(), headlessAssetPolicy({ hosts: ['x.org'] }));
    const s = await BigWigRuntime.open(entry, { assets, signal, limits });
    await values(s.read({ system: 'chrA', start: 0, end: 10_000, level: 256, shapes: ['values'] }, signal));
    expect(log.length).toBeGreaterThan(0);
    expect(log.every((l) => l.range.startsWith('bytes='))).toBe(true);
    s.dispose();
  });

  // the writer's variants: big-endian quick, the other two medium (a file written each, a few hundred ms)
  for (const [name, o, tags] of [
    ['big-endian', { bigEndian: true }, []],
    ['uncompressed', { compress: false }, ['medium']],
    ['two-items-a-node', { blockSize: 2, itemsPerSection: 7 }, ['medium']],
  ] as [string, WriterOptions, string[]][])
    it(`${name} files read the same`, { tags }, async () => {
      const bytes = await writeBigWig(Chroms, { ...Options, ...o });
      const s = await BigWigRuntime.open(entry, bytesContext({ [Url]: bytes }));
      const v = (await values(
        s.read({ system: 'chrB', start: 50, end: 1500, level: 1, shapes: ['values'] }, signal),
      )) as ValuesPositions;
      expect(Array.from(v.value), JSON.stringify(o)).toEqual(Array.from(expected('chrB', 50, 1500)));
      const z = (await values(
        s.read({ system: 'chrA', start: 0, end: 10_000, level: 256, shapes: ['values'] }, signal),
      )) as ValuesBins;
      expectBins(z, binValues({ form: 'positions', start: 0, value: Chroms[0].values }, 256, 10_000));
      s.dispose();
    });

  it('summaries that do not start on the level’s multiples are shared by overlap', async () => {
    // one record over [8, 24) of 16 positions, mean 1: half of it in each 16-wide bin
    const s = await BigWigRuntime.open(entry, bytesContext({ [Url]: await shiftedZoom() }));
    const z = (await values(
      s.read({ system: 'chrS', start: 0, end: 32, level: 16, shapes: ['values'] }, signal),
    )) as ValuesBins;
    expect(Array.from(z.count)).toEqual([8, 8]);
    expect(Array.from(z.mean)).toEqual([1, 1]);
    s.dispose();
  });

  it('an aborted read rejects as aborted; the next read answers', async () => {
    serve({ [Url]: Fixture });
    const s = await BigWigRuntime.open(entry, rangeContext());
    const stop = new AbortController();
    const p = s.read({ system: 'chrA', start: 4096, end: 5120, level: 1, shapes: ['values'] }, stop.signal);
    stop.abort();
    await expect(p).rejects.toMatchObject({ name: 'AbortError' });
    expect((await s.read({ system: 'chrA', start: 4096, end: 5120, level: 1, shapes: ['values'] }, signal)).ok).toBe(
      true,
    );
    s.dispose();
  });

  it('a file that is not a bigWig, one cut short, and a missing one are issues on the source, never wrong values', async () => {
    const open = (bytes: Uint8Array) => BigWigRuntime.open(entry, bytesContext({ [Url]: bytes }));
    await expect(open(new Uint8Array(200))).rejects.toSatisfy(
      (e) => e instanceof SourceOpenError && e.issues[0].code === 'seq.source.malformed',
    );
    await expect(open(Fixture.slice(0, 40))).rejects.toBeInstanceOf(SourceOpenError);
    await expect(BigWigRuntime.open(entry, bytesContext({}))).rejects.toSatisfy(
      (e) => e instanceof SourceOpenError && e.issues[0].code === 'seq.source.unreadable',
    );
    await expect(BigWigRuntime.open({ address: '/d', format: 'bigwig' }, bytesContext({}))).rejects.toSatisfy(
      (e) => e instanceof SourceOpenError && e.issues[0].code === 'seq.source.no_data',
    );
    // cut after its index: the header opens, a region's blocks are past the end
    const s = await open(Fixture.slice(0, Fixture.length - 64));
    const answers = await Promise.all(
      Chroms.map((c) =>
        s.read({ system: c.name, start: 0, end: c.values.length, level: 256, shapes: ['values'] }, signal),
      ),
    );
    expect(answers.some((a) => !a.ok && a.kind === 'malformed')).toBe(true);
    s.dispose();
  });

  it('its synthetic twin: the same answers after a delay on a fake clock', async () => {
    const clock = new FakeClock();
    const twin = delayedRuntime(BigWigRuntime, { clock, delayMs: 40 });
    const s = await twin.open(entry, bytesContext({ [Url]: Fixture }));
    let answered: RegionAnswer | undefined;
    s.read({ system: 'chrA', start: 0, end: 1024, level: 16, shapes: ['values'] }, signal).then((a) => (answered = a));
    await clock.advance(39);
    expect(answered).toBeUndefined();
    await clock.advance(1);
    await vi.waitFor(() => expect(answered?.ok).toBe(true));
    s.dispose();
  });

  // the performance target: a bigWig region at the right zoom level from a local file, under 50 ms to decoded; its
  // 400,000 positions written first (a third of a second): medium unbounded, heavy bounded; the quick lane reads zoom
  // levels from the small fixture above
  describe('the performance target', { tags: ['medium'] }, () =>
    itWithin('a region at the zoom level that fits, from a local file, decoded under 50 ms', async (within) => {
      const length = 400_000;
      const big = await bigFixture(length);
      // a view of the whole chromosome over 1,024 device columns: the coarsest level no wider than a column (EC-170)
      const level = 256;
      const manager = new DefaultAssetManager();
      const file = Asset.File(new File([big as Uint8Array<ArrayBuffer>], 'big.bw'));
      const context: OpenContext = {
        assets: {
          readText: () => Promise.reject(new Error('text')),
          readBytes: () => Promise.reject(new Error('whole')),
          openRange: () => openRanges(manager, file).run(),
        } as AssetReader,
        signal,
        limits,
      };
      const v = (await within(50, async () => {
        const s = await BigWigRuntime.open(entry, context);
        const a = await values(s.read({ system: 'chr1', start: 0, end: length, level, shapes: ['values'] }, signal));
        s.dispose();
        return a;
      })) as ValuesBins;
      expect(v.count.length).toBe(Math.ceil(length / level));
      expect(v.count.reduce((n, x) => n + x, 0)).toBe(length);
    }),
  );
});

// security.md: a file is untrusted input. Each crafted file is refused with a coded issue, before it allocates.
describe('a crafted bigWig is refused by the runtime’s bounds', () => {
  /** The fixture with `patch` applied to a copy (offsets from its header: little-endian). */
  function crafted(patch: (dv: DataView, at: { chromTree: number; fullIndex: number }) => void) {
    const bytes = Fixture.slice();
    const dv = new DataView(bytes.buffer);
    patch(dv, { chromTree: Number(dv.getBigUint64(8, true)), fullIndex: Number(dv.getBigUint64(24, true)) });
    return bytes;
  }
  const whole = { system: 'chrA', start: 0, end: 10_000, level: 1, shapes: ['values' as const] };

  it('an index node of 65,536 items is refused, not read as 2 MB a node', async () => {
    const bytes = crafted((dv, at) => dv.setUint32(at.fullIndex + 4, 65_536, true));
    const s = await BigWigRuntime.open(entry, bytesContext({ [Url]: bytes }));
    const a = await s.read(whole, signal);
    expect(a).toMatchObject({ ok: false, kind: 'malformed', issue: { code: 'seq.source.malformed' } });
    s.dispose();
  });

  it('an index whose children all point at one node (itself): read once a level, refused at the depth bound', async () => {
    const bytes = crafted((dv, at) => {
      const root = at.fullIndex + 48;
      expect(dv.getUint8(root)).toBe(0);
      const count = dv.getUint16(root + 2, true);
      for (let i = 0; i < count; i++) {
        dv.setUint32(root + 4 + 24 * i, 0, true);
        dv.setUint32(root + 4 + 24 * i + 8, 0xffffffff, true);
        dv.setBigUint64(root + 4 + 24 * i + 16, BigInt(root), true);
      }
    });
    const reads: number[] = [];
    const ctx = bytesContext({ [Url]: bytes });
    const s = await BigWigRuntime.open(entry, {
      ...ctx,
      assets: {
        ...ctx.assets,
        openRange: async () => {
          const r = { size: bytes.length };
          return {
            ...r,
            read: async (ranges) => {
              reads.push(ranges.length);
              return ranges.map((x) => bytes.subarray(x.offset, Math.min(bytes.length, x.offset + x.length)));
            },
            dispose() {},
          };
        },
      },
    });
    const n = reads.length;
    expect(await s.read(whole, signal)).toMatchObject({ ok: false, kind: 'malformed' });
    // the header, then one node a level: never the 4^depth a tree of repeated children names
    expect(Math.max(...reads.slice(n))).toBe(1);
    s.dispose();
  });

  it('an offset past 2^53 is malformed, not a RangeError of the range reader', async () => {
    const bytes = crafted((dv) => dv.setUint32(24 + 4, 0xffffffff, true));
    serve({ [Url]: bytes });
    await expect(BigWigRuntime.open(entry, rangeContext())).rejects.toSatisfy(
      (e) => e instanceof SourceOpenError && e.issues[0].code === 'seq.source.malformed',
    );
    // a leaf's block offset past it: the region is malformed
    const leafy = crafted((dv, at) => {
      for (let o = at.fullIndex + 48; ; ) {
        const leaf = dv.getUint8(o) === 1;
        if (leaf) return dv.setUint32(o + 4 + 16 + 4, 0xffffffff, true);
        o = Number(dv.getBigUint64(o + 4 + 16, true));
      }
    });
    serve({ [Url]: leafy });
    const s = await BigWigRuntime.open(entry, rangeContext());
    const answers = await Promise.all(
      Chroms.map((c) =>
        s.read({ system: c.name, start: 0, end: c.values.length, level: 1, shapes: ['values'] }, signal),
      ),
    );
    expect(answers.some((a) => !a.ok && a.kind === 'malformed')).toBe(true);
    s.dispose();
  });

  it('a chromosome tree of 65,535 items of 1,032 bytes a node is refused at open', async () => {
    const bytes = crafted((dv, at) => {
      dv.setUint32(at.chromTree + 4, 65_535, true);
      dv.setUint32(at.chromTree + 8, 1024, true);
    });
    await expect(BigWigRuntime.open(entry, bytesContext({ [Url]: bytes }))).rejects.toSatisfy(
      (e) => e instanceof SourceOpenError && e.issues[0].code === 'seq.source.malformed',
    );
  });

  it('a tree level wider than its bound, and blocks inflating past the region’s bound, are refused', async () => {
    const narrow = await bigWigRuntime({ maxLevelNodes: 2 }).open(entry, bytesContext({ [Url]: Fixture }));
    expect(await narrow.read(whole, signal)).toMatchObject({ ok: false, kind: 'malformed' });
    narrow.dispose();
    const small = await bigWigRuntime({ maxInflatedBytes: 2000 }).open(entry, bytesContext({ [Url]: Fixture }));
    expect(await small.read(whole, signal)).toMatchObject({
      ok: false,
      kind: 'too-large',
      issue: { code: 'seq.source.too_large' },
    });
    // a region within the bound is read
    expect((await small.read({ ...whole, end: 64 }, signal)).ok).toBe(true);
    small.dispose();
  });
});

const bigFixtures = new Map<number, Uint8Array>();
async function bigFixture(length: number) {
  if (!bigFixtures.has(length))
    bigFixtures.set(
      length,
      await writeBigWig([{ name: 'chr1', values: syntheticSignal(9, length) }], {
        zooms: [64, 256, 1024],
        itemsPerSection: 1024,
        blockSize: 64,
      }),
    );
  return bigFixtures.get(length)!;
}

/** A file whose zoom-16 record of `chrS` spans [8, 24): written, then its record moved by 8. */
async function shiftedZoom() {
  const values = new Float32Array(32).fill(Number.NaN);
  values.fill(1, 0, 16);
  const bytes = await writeBigWig([{ name: 'chrS', values }], { zooms: [16], compress: false });
  const dv = new DataView(bytes.buffer);
  // the zoom header's data offset: its count (4 bytes), then the one record (chrom, start, end)
  const data = Number(dv.getBigUint64(64 + 8, true));
  dv.setUint32(data + 4 + 4, 8, true);
  dv.setUint32(data + 4 + 8, 24, true);
  return bytes;
}
