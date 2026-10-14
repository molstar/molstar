/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Records the byte ranges the readers read for _test/cases.ts from the real public files, and the expected answers from
// sources that do not use the readers. Run by hand (`node scripts/record.ts [twobit|bigwig|bigbed|vcf ...]`), never in
// tests: it needs the network and overwrites _test/fixtures.
//   expected.json: UCSC's REST API (api.genome.ucsc.edu getData: sequence, gc5BaseBw, knownGene = gencodeV50.bb) for the
//   2bit, BigWig and BigBed cases, and for the VCF the first/last megabytes of the file inflated with Node's gunzip
//   and filtered by line (no index, no BGZF splitting of ours).

import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { constants, gunzipSync } from 'node:zlib';
import { CASES, type Key, makers, URLS } from '../_test/stream/gmod/cases.ts';
import type { SliceMeta } from '../_test/stream/gmod/slices.ts';
import { fetchRanges, type RangeReader } from '../src/stream/gmod/index.ts';

const fixtures = new URL('../_test/stream/gmod/fixtures/', import.meta.url);
const API = 'https://api.genome.ucsc.edu/getData';
const only = process.argv.slice(2);
const want = (k: string) => only.length === 0 || only.includes(k);
mkdirSync(fixtures, { recursive: true });
const slicesFile = new URL('slices.json', fixtures);
const metas: Record<string, SliceMeta> = existsSync(slicesFile) ? JSON.parse(readFileSync(slicesFile, 'utf8')) : {};
const expectedFile = new URL('expected.json', fixtures);
const expected: Record<string, unknown> = existsSync(expectedFile)
  ? JSON.parse(readFileSync(expectedFile, 'utf8'))
  : {};

const json = async <T>(url: string): Promise<T> => {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return (await r.json()) as T;
};

/** A reader that records every distinct read it passes on. */
async function recording(key: Key, url: string) {
  const head = await fetch(url, { method: 'HEAD' });
  const size = Number(head.headers.get('content-length'));
  const base = fetchRanges(url);
  const got = new Map<string, Uint8Array>();
  const read: RangeReader = async (offset, length, signal) => {
    const k = `${offset}:${length}`;
    if (!got.has(k)) got.set(k, await base(offset, length, signal));
    return got.get(k) as Uint8Array;
  };
  const save = () => {
    const entries = [...got.entries()].map(([k, b]) => ({
      offset: Number(k.split(':')[0]),
      asked: Number(k.split(':')[1]),
      b,
    }));
    entries.sort((a, b) => a.offset - b.offset || a.asked - b.asked);
    writeFileSync(new URL(`${key}.bin`, fixtures), Buffer.concat(entries.map((e) => e.b)));
    metas[key] = {
      url,
      size,
      modified: head.headers.get('last-modified'),
      chunks: entries.map((e) => [e.offset, e.asked, e.b.length]),
    };
    writeFileSync(slicesFile, `${JSON.stringify(metas, null, 1)}\n`);
    console.log(
      key,
      entries.length,
      'ranges',
      entries.reduce((n, e) => n + e.b.length, 0),
      'bytes of',
      size,
    );
  };
  return { read, save, size };
}

const signal = new AbortController().signal;
async function run(key: keyof typeof URLS, read: RangeReader, index?: Uint8Array) {
  const provider =
    key === 'vcf'
      ? makers.vcf(read, index as Uint8Array)
      : (makers[key] as (r: RangeReader) => ReturnType<typeof makers.bigwig>)(read);
  const out: Record<string, unknown> = {};
  for (const [name, r] of CASES[key]) out[name] = await provider.fetch(r, signal);
  return out;
}

if (want('twobit')) {
  const rec = await recording('twobit', URLS.twobit);
  await run('twobit', rec.read);
  rec.save();
  const e: Record<string, unknown> = {};
  for (const [name, r] of CASES.twobit.slice(1)) {
    const end = Math.min(r.stop, r.system === 'chrM' ? 16_569 : r.stop);
    e[name] = (
      await json<{ dna: string }>(`${API}/sequence?genome=hg38;chrom=${r.system};start=${r.start};end=${end}`)
    ).dna;
  }
  expected.twobit = e;
}

if (want('bigwig')) {
  const rec = await recording('bigwig', URLS.bigwig);
  const got = (await run('bigwig', rec.read)) as Record<string, { manifest: { levels: number[] } }>;
  rec.save();
  const e: Record<string, unknown> = { levels: got.manifest.manifest.levels };
  const levels: number[] = got.manifest.manifest.levels;
  for (const [name, r] of CASES.bigwig.slice(1)) {
    const raw = (
      await json<{ gc5BaseBw: { start: number; end: number; value: number }[] }>(
        `${API}/track?genome=hg38;track=gc5BaseBw;chrom=${r.system};start=${r.start};end=${r.stop};maxItemsOutput=-1`,
      )
    ).gc5BaseBw;
    if (r.level < 0) e[name] = raw.map((x) => [x.start, x.end, x.value]);
    else {
      // The raw 5-base windows on the level's grid: per bin the minimum, maximum and mean of the windows' values.
      const w = levels[r.level];
      const bins = Math.ceil((r.stop - r.start) / w);
      const acc = Array.from({ length: bins }, () => ({ min: Infinity, max: -Infinity, sum: 0, n: 0 }));
      for (const x of raw)
        for (let p = x.start; p < x.end; p++) {
          const a = acc[Math.floor((p - r.start) / w)];
          if (!a) continue;
          a.min = Math.min(a.min, x.value);
          a.max = Math.max(a.max, x.value);
          a.sum += x.value;
          a.n++;
        }
      e[name] = acc.map((a) => (a.n ? [a.min, a.max, a.sum / a.n] : null));
    }
  }
  expected.bigwig = e;
}

if (want('bigbed')) {
  const rec = await recording('bigbed', URLS.bigbed);
  await run('bigbed', rec.read);
  rec.save();
  const e: Record<string, unknown> = {};
  for (const [name, r] of CASES.bigbed.slice(1)) {
    const list = (
      await json<{ knownGene: Record<string, string>[] }>(
        `${API}/track?genome=hg38;track=knownGene;chrom=${r.system};start=${r.start};end=${r.stop}`,
      )
    ).knownGene;
    e[name] = list.map((g) => ({
      name: g.name,
      start: g.chromStart,
      end: g.chromEnd,
      strand: g.strand,
      exons: g.chromStarts
        .split(',')
        .filter(Boolean)
        .map((s, k) => [
          Number(g.chromStart) + Number(s),
          Number(g.chromStart) + Number(s) + Number(g.blockSizes.split(',')[k]),
        ]),
    }));
  }
  expected.bigbed = e;
}

if (want('vcf')) {
  const tbi = new Uint8Array(await (await fetch(`${URLS.vcf}.tbi`)).arrayBuffer());
  writeFileSync(new URL('tbi.bin', fixtures), tbi);
  metas.tbi = { url: `${URLS.vcf}.tbi`, size: tbi.length, modified: null, chunks: [[0, tbi.length, tbi.length]] };
  const rec = await recording('vcf', URLS.vcf);
  await run('vcf', rec.read, tbi);
  rec.save();
  // The same lines without the index: from the first (chr 1) and the last (MT) megabytes, as Node inflates them.
  const lines = async (from: number, to: number) => {
    let bytes = await rec_range(from, to);
    // Start at a BGZF block: gzip magic, the BC subfield.
    let at = 0;
    while (
      at < bytes.length &&
      !(bytes[at] === 31 && bytes[at + 1] === 139 && bytes[at + 12] === 66 && bytes[at + 13] === 67)
    )
      at++;
    bytes = bytes.subarray(at);
    return gunzipSync(bytes, { finishFlush: constants.Z_SYNC_FLUSH }).toString('latin1').split('\n');
  };
  const rec_range = fetchRanges(URLS.vcf);
  const size = rec.size;
  const pick = (all: string[], chrom: string, start: number, stop: number) =>
    all
      .filter((l) => l && l[0] !== '#')
      .map((l) => l.split('\t'))
      .filter((f) => f[0] === chrom && Number(f[1]) - 1 < stop && Number(f[1]) - 1 + f[3].length > start)
      .flatMap((f) =>
        f[4].split(',').map((alt) => ({
          start: Number(f[1]) - 1,
          ref: f[3],
          alt,
          name: f[2],
          clnsig: /(?:^|;)CLNSIG=([^;]*)/.exec(f[7])?.[1] ?? '',
        })),
      );
  const [, r1] = CASES.vcf[1];
  const [, r2] = CASES.vcf[2];
  expected.vcf = {
    'chr1 variants': pick(await lines(0, 3_000_000), '1', r1.start, r1.stop),
    'MT variants': pick(await lines(size - 2_000_000, size), 'MT', r2.start, r2.stop),
  };
}

writeFileSync(expectedFile, `${JSON.stringify({ recorded: '2026-10-12', ...expected })}\n`);
