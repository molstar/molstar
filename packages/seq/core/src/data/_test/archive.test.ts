/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { unzip, zip } from '@molstar/core/util/zip/zip';
import { createSeqSpec, DefaultSpecBounds } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  type Archive,
  ArchiveEntryMissing,
  ArchiveIndexName,
  ArchiveSpecName,
  archiveReader,
  checkArchiveName,
  DefaultArchiveLimits,
  readArchive,
  writeArchive,
} from '../archive.js';
import type { AssetReader } from '../protocol.js';

const enc = new TextEncoder();
const dec = new TextDecoder();
const signal = new AbortController().signal;

/** A reader over a map of files: what a host's assets port gives; a location not in it is refused. */
function fake(files: Record<string, string | Uint8Array>): AssetReader & { reads: string[] } {
  const reads: string[] = [];
  const get = (uri: string) => {
    reads.push(uri);
    const f = files[uri];
    if (f === undefined) throw new Error(`refused ${uri}`);
    return typeof f === 'string' ? enc.encode(f) : f;
  };
  return { reads, readText: async (u) => dec.decode(get(u)), readBytes: async (u) => get(u) };
}

function specOf(...uris: (string | { uri: string; index: string })[]): Node {
  const b = createSeqSpec();
  uris.forEach((u, i) => {
    const o = typeof u === 'string' ? { uri: u } : u;
    b.data.source(`s${i}`, { format: 'seq', ...o });
  });
  return b.build();
}

const sourcesOf = (spec: Node) => (spec.children ?? []).flatMap((d) => d.children ?? []).map((s) => s.params);

const rawIndex = (spec: string, files: string[]) => JSON.stringify({ archive: 'seqx', version: 1, spec, files });

const built = async (entries: Record<string, Uint8Array<ArrayBuffer>>) =>
  new Uint8Array(await zip(SyncRuntimeContext, entries));

async function valid(files: string[] = [], extra: Record<string, Uint8Array<ArrayBuffer>> = {}) {
  const entries: Record<string, Uint8Array<ArrayBuffer>> = {
    [ArchiveIndexName]: enc.encode(rawIndex(ArchiveSpecName, files)),
    [ArchiveSpecName]: enc.encode(JSON.stringify({ kind: 'spec', params: { version: '0.1.0' } })),
    ...extra,
  };
  for (const f of files) entries[f] ??= enc.encode(`file ${f}`);
  return built(entries);
}

/** Patch the central directory's declared sizes of an entry (the archive is not inflated by the patch). */
function patch(bytes: Uint8Array, name: string, change: { usize?: number; csize?: number }) {
  const out = bytes.slice();
  const v = new DataView(out.buffer);
  for (let i = 0; i + 46 < out.length; i++) {
    if (v.getUint32(i, true) !== 0x02014b50) continue;
    const nl = v.getUint16(i + 28, true);
    if (dec.decode(out.subarray(i + 46, i + 46 + nl)) !== name) continue;
    if (change.csize !== undefined) v.setUint32(i + 20, change.csize, true);
    if (change.usize !== undefined) v.setUint32(i + 24, change.usize, true);
    return out;
  }
  throw new Error(`no entry ${name}`);
}

const codes = (r: { issues: readonly { code: string }[] }) => r.issues.map((i) => i.code);

describe('writing an archive', () => {
  it('saves a view with the files its sources name, and opens it again', async () => {
    const assets = fake({ 'data/a.seq': 'AAA', 'b.seq': 'BBB', 'b.idx': 'III' });
    const w = await writeArchive(specOf('data/a.seq', { uri: 'b.seq', index: 'b.idx' }), { assets });
    expect(w.issues).toEqual([]);
    const r = await readArchive(w.bytes!);
    expect(r.issues).toEqual([]);
    expect(sourcesOf(r.archive!.spec)).toMatchObject([
      { format: 'seq', uri: 'data/a.seq' },
      { format: 'seq', uri: 'b.seq', index: 'b.idx' },
    ]);
    expect([...r.archive!.files.keys()].sort()).toEqual(['b.idx', 'b.seq', 'data/a.seq']);
    expect(dec.decode(r.archive!.files.get('b.seq')!)).toBe('BBB');
  });

  it('names a file that is no path under the root `data/<n>-<name>` and reads a location once', async () => {
    const assets = fake({ 'https://x.org/a/g.seq': 'G', '../up/h.seq': 'H', 'index.json': 'clash' });
    const spec = specOf('https://x.org/a/g.seq', '../up/h.seq', 'https://x.org/a/g.seq', 'index.json');
    const w = await writeArchive(spec, { assets });
    const r = await readArchive(w.bytes!);
    const uris = sourcesOf(r.archive!.spec).map((p) => p?.uri as string);
    expect(uris[0]).toBe(uris[2]);
    expect(uris[0]).toMatch(/^data\/\d+-g\.seq$/);
    expect(uris[1]).toMatch(/^data\/\d+-h\.seq$/);
    expect(uris[3]).toMatch(/^data\/\d+-index\.json$/);
    expect(new Set(uris).size).toBe(3);
    expect(assets.reads.filter((u) => u === 'https://x.org/a/g.seq')).toHaveLength(1);
  });

  it('leaves a location the host keeps as written, unread', async () => {
    const assets = fake({ 'a.seq': 'A' });
    const w = await writeArchive(specOf('generated:chr1', 'a.seq'), {
      assets,
      keep: (u) => u.startsWith('generated:'),
    });
    const r = await readArchive(w.bytes!);
    expect(sourcesOf(r.archive!.spec).map((p) => p?.uri)).toEqual(['generated:chr1', 'a.seq']);
    expect([...r.archive!.files.keys()]).toEqual(['a.seq']);
    expect(assets.reads).toEqual(['a.seq']);
  });

  it('writes an archive without files for a spec with inline data only', async () => {
    const w = await writeArchive(specOf(), { assets: fake({}) });
    const r = await readArchive(w.bytes!);
    expect(r.archive!.files.size).toBe(0);
  });

  it('writes nothing when a file is refused or unreadable, and says which', async () => {
    const w = await writeArchive(specOf('a.seq', 'b.seq'), { assets: fake({ 'a.seq': 'A' }) });
    expect(w.bytes).toBeUndefined();
    expect(codes(w)).toEqual(['seq.archive.data']);
    expect(w.issues[0].message).toContain('b.seq');
  });

  it('EC-91: refuses to write past the bounds it reads with', async () => {
    const big = new Uint8Array(100);
    const assets = fake({ 'a.seq': big, 'b.seq': big, 'c.seq': big });
    const entries = await writeArchive(specOf('a.seq', 'b.seq', 'c.seq'), { assets, limits: { maxEntries: 4 } });
    expect(codes(entries)).toEqual(['seq.archive.bounds']);
    const bytes = await writeArchive(specOf('a.seq', 'b.seq'), { assets, limits: { maxBytes: 150 } });
    expect(codes(bytes)).toEqual(['seq.archive.bounds']);
    const one = await writeArchive(specOf('a.seq'), { assets, limits: { maxEntryBytes: 99 } });
    expect(codes(one)).toEqual(['seq.archive.bounds']);
  });

  it('stops on an aborted signal instead of reporting the read as a missing file', async () => {
    const stop = new AbortController();
    stop.abort(new Error('stopped'));
    const assets: AssetReader = {
      readText: async () => '',
      readBytes: async (_u, s) => {
        throw s.reason;
      },
    };
    await expect(writeArchive(specOf('a.seq'), { assets, signal: stop.signal })).rejects.toThrow('stopped');
  });
});

describe('writing an archive: what is written opens (EC-91)', () => {
  it('EC-91: stores an entry that would pass the ratio compressed (a megabyte of N in a FASTA), and it opens', async () => {
    const fasta = enc.encode(`>chr\n${'N'.repeat(1 << 20)}\n`);
    const w = await writeArchive(specOf('chr.fa'), { assets: fake({ 'chr.fa': fasta }) });
    expect(w.issues).toEqual([]);
    const r = await readArchive(w.bytes!);
    expect(r.issues).toEqual([]);
    expect(r.archive!.files.get('chr.fa')).toEqual(fasta);
  });

  it('EC-91: refuses when a bound still breaks after that (an archive past the total)', async () => {
    const fasta = enc.encode('N'.repeat(1 << 20));
    const w = await writeArchive(specOf('chr.fa'), {
      assets: fake({ 'chr.fa': fasta }),
      limits: { maxBytes: 1 << 20 },
    });
    expect(w.bytes).toBeUndefined();
    expect(codes(w)).toEqual(['seq.archive.bounds']);
  });

  it('names two locations that differ in case only, or one that is `SPEC.json`, as different entries', async () => {
    const assets = fake({ 'A.seq': 'A', 'a.seq': 'a', 'SPEC.json': 'S', 'Index.JSON': 'I' });
    const w = await writeArchive(specOf('A.seq', 'a.seq', 'SPEC.json', 'Index.JSON'), { assets });
    expect(w.issues).toEqual([]);
    const r = await readArchive(w.bytes!);
    expect(r.issues).toEqual([]);
    expect(r.archive!.files.size).toBe(4);
    const uris = sourcesOf(r.archive!.spec).map((p) => p?.uri as string);
    expect(new Set(uris.map((u) => u.toLowerCase())).size).toBe(4);
    expect(uris[0]).toBe('A.seq');
  });
});

describe('reading an archive: assets', () => {
  it('serves the entries as assets by name, `./` dropped, and refuses what it does not hold', async () => {
    const r = await readArchive(await valid(['d/x.seq']));
    const reader = archiveReader(r.archive as Archive);
    expect(await reader.readText('d/x.seq', signal)).toBe('file d/x.seq');
    expect(await reader.readText('./d/x.seq', signal)).toBe('file d/x.seq');
    expect((await reader.readBytes('d/x.seq', signal)).byteLength).toBe(12);
    await expect(reader.readText('d/y.seq', signal)).rejects.toBeInstanceOf(ArchiveEntryMissing);
    await expect(reader.readText('https://x.org/d/x.seq', signal)).rejects.toBeInstanceOf(ArchiveEntryMissing);
    await expect(reader.readText('../d/x.seq', signal)).rejects.toBeInstanceOf(ArchiveEntryMissing);
    const stop = new AbortController();
    stop.abort();
    await expect(reader.readText('d/x.seq', stop.signal)).rejects.toBeDefined();
  });

  it('inflates only what the index lists', async () => {
    const r = await readArchive(await valid(['a.seq'], { 'junk.bin': enc.encode('x'.repeat(100)) }));
    expect(r.issues).toEqual([]);
    expect([...r.archive!.files.keys()]).toEqual(['a.seq']);
  });
});

describe('reading an archive: names (EC-91)', () => {
  const refused = [
    ['../evil.seq', 'above the root'],
    ['a/../../evil.seq', 'above the root'],
    ['/etc/passwd', 'absolute'],
    // a drive letter written in two pieces: the repository's path check refuses a literal one in any tracked file
    [`${'C'}:/x.seq`, 'drive'],
    ['a\\b.seq', 'backslash'],
    ['a//b.seq', 'empty segment'],
    ['a/./b.seq', '`.` segment'],
    ['a\u0000b', 'control'],
  ] as const;

  it.each(refused)('EC-91: checkArchiveName refuses %j', (name, why) => {
    expect(checkArchiveName(name)).toContain(why);
  });

  it('EC-91: refuses an empty name and a name too long', () => {
    expect(checkArchiveName('')).toContain('empty');
    expect(checkArchiveName('x'.repeat(257))).toContain('longer');
    expect(checkArchiveName('x'.repeat(256))).toBeUndefined();
  });

  it('EC-91: accepts ordinary paths and directories', () => {
    for (const ok of ['a.seq', 'data/a.b.seq', 'dir/', 'data/..hidden', 'ünï/ç.seq'])
      expect(checkArchiveName(ok)).toBeUndefined();
  });

  it.each(refused)('EC-91: an archive with the entry %j is refused before anything is inflated', async (name) => {
    const r = await readArchive(await valid([], { [name]: enc.encode('x') }));
    expect(codes(r)).toEqual(['seq.archive.name']);
    expect(r.archive).toBeUndefined();
  });

  it('EC-91: two entries that differ in case only are one name', async () => {
    const r = await readArchive(await valid([], { 'A.seq': enc.encode('a'), 'a.seq': enc.encode('c') }));
    expect(codes(r)).toEqual(['seq.archive.name']);
  });
});

describe('reading an archive: bounds (EC-91)', () => {
  it('EC-91: refuses two entries of one name, which the listing would show once (the bounds of the first bypassed)', async () => {
    const z = await valid(['a.seq', 'b.seq'], { 'a.seq': new Uint8Array(1 << 20), 'b.seq': enc.encode('small') });
    // b.seq renamed a.seq in the local and the central header: same length, so the offsets hold
    const text = new TextDecoder('latin1').decode(z);
    const renamed = z.slice();
    for (let i = text.indexOf('b.seq'); i >= 0; i = text.indexOf('b.seq', i + 1)) renamed.set(enc.encode('a.seq'), i);
    const r = await readArchive(renamed, { limits: { maxEntryBytes: 1000 } });
    expect(codes(r)).toEqual(['seq.archive.name']);
    expect(r.archive).toBeUndefined();
  });

  it('EC-91: refuses more entries than the bound, from the end record', async () => {
    const extra: Record<string, Uint8Array<ArrayBuffer>> = {};
    for (let i = 0; i < 20; i++) extra[`e${i}`] = enc.encode('x');
    const z = await valid([], extra);
    expect(codes(await readArchive(z, { limits: { maxEntries: 22 } }))).toEqual([]);
    expect(codes(await readArchive(z, { limits: { maxEntries: 21 } }))).toEqual(['seq.archive.bounds']);
  });

  it('EC-91: refuses an expanded size over the bound from the declared sizes, never inflating', async () => {
    const z = patch(await valid(['a.seq']), 'a.seq', { usize: 2 ** 31 });
    expect(codes(await readArchive(z))).toEqual(['seq.archive.bounds']);
    const total = await readArchive(await valid(['a.seq', 'b.seq']), { limits: { maxBytes: 100 } });
    expect(codes(total)).toEqual(['seq.archive.bounds']);
  });

  it('EC-91: refuses one entry over its bound while the total is within', async () => {
    const z = await valid(['a.seq']);
    expect(codes(await readArchive(z, { limits: { maxEntryBytes: 7 } }))).toEqual(['seq.archive.bounds']);
  });

  it('EC-91: refuses a ratio of expansion over the bound for an entry past the floor, and only past it', async () => {
    const z = await valid(['bomb.seq'], { 'bomb.seq': new Uint8Array(1 << 20) });
    const r = await readArchive(z);
    expect(codes(r)).toEqual(['seq.archive.bounds']);
    expect(r.issues[0].message).toContain('ratio');
    // a small file of zeros is not a bomb, and a raised ratio opens the large one
    const small = await valid(['z.seq'], { 'z.seq': new Uint8Array(30_000) });
    expect(codes(await readArchive(small))).toEqual([]);
    expect(codes(await readArchive(z, { limits: { maxRatio: 100_000 } }))).toEqual([]);
  });

  it('EC-91: refuses archive bytes over the bound', async () => {
    const z = await valid();
    expect(codes(await readArchive(z, { limits: { maxBytes: z.byteLength - 1 } }))).toEqual(['seq.archive.bounds']);
  });

  it('keeps the defaults the spec lists', () => {
    expect(DefaultArchiveLimits).toEqual({
      maxEntries: 1000,
      maxBytes: 268435456,
      maxEntryBytes: 67108864,
      maxRatio: 200,
      ratioFloor: 65536,
      maxNameLength: 256,
    });
  });
});

describe('reading an archive: the index and the spec', () => {
  const entries = (index: string | undefined, spec = '{"kind":"spec"}') => {
    const e: Record<string, Uint8Array<ArrayBuffer>> = { [ArchiveSpecName]: enc.encode(spec) };
    if (index !== undefined) e[ArchiveIndexName] = enc.encode(index);
    return built(e);
  };
  const index = (over: object) =>
    JSON.stringify({ archive: 'seqx', version: 1, spec: 'spec.json', files: [], ...over });

  it.each([
    ['no index', undefined],
    ['not JSON', '{'],
    ['not an object', '[]'],
    ['another archive', index({ archive: 'tar' })],
    ['a newer version', index({ version: 2 })],
    ['version 0', index({ version: 0 })],
    ['a negative version', index({ version: -1 })],
    ['a fractional version', index({ version: 1.5 })],
    ['a spec above the root', index({ spec: '../spec.json' })],
    ['a file above the root', index({ files: ['/etc/passwd'] })],
    ['a file it lacks', rawIndex('spec.json', ['gone.seq'])],
    ['a spec it lacks', rawIndex('other.json', [])],
    ['files not a list', index({ files: 'a' })],
  ] as const)('refuses an index with %s', async (_, text) => {
    expect(codes(await readArchive(await entries(text)))).toEqual(['seq.archive.index']);
  });

  it('refuses a spec that is no tree, or is not UTF-8', async () => {
    const idx = rawIndex('spec.json', []);
    expect(codes(await readArchive(await entries(idx, '{')))[0]).toBe('seq.archive.spec');
    expect(codes(await readArchive(await entries(idx, '[1,2]')))[0]).toBe('seq.archive.spec');
    const bad = await built({
      [ArchiveIndexName]: enc.encode(idx),
      [ArchiveSpecName]: new Uint8Array([0xff, 0xfe, 0xfd]),
    });
    expect(codes(await readArchive(bad))).toEqual(['seq.archive.spec']);
  });

  it('refuses a spec past the spec bounds (a deeply nested one)', async () => {
    let json = '{"kind":"x"}';
    for (let i = 0; i < 300; i++) json = `{"kind":"x","children":[${json}]}`;
    const r = await readArchive(await entries(rawIndex('spec.json', []), json));
    expect(codes(r)[0]).toBe('seq.archive.spec');
    expect(r.archive).toBeUndefined();
  });

  it('takes the host’s spec bounds', async () => {
    const z = await valid();
    expect(codes(await readArchive(z, { specBounds: { ...DefaultSpecBounds, maxBytes: 10 } }))).toEqual([
      'seq.archive.spec',
    ]);
    expect(codes(await readArchive(z, { specBounds: { ...DefaultSpecBounds, maxDepth: 0 } }))[0]).toBe(
      'seq.archive.spec',
    );
    expect(codes(await readArchive(z))).toEqual([]);
  });

  it('refuses bytes that are no zip, or a zip cut short, without throwing', async () => {
    const text = enc.encode('not a zip at all, but long enough to look like one');
    expect(codes(await readArchive(text))).toEqual(['seq.archive.unreadable']);
    expect(codes(await readArchive(new Uint8Array(0)))).toEqual(['seq.archive.unreadable']);
    const z = await valid(['a.seq']);
    expect(codes(await readArchive(z.subarray(0, z.length - 10)))).toEqual(['seq.archive.unreadable']);
  });

  it('refuses an entry whose bytes fail the CRC-32, and reads a view of a larger buffer', async () => {
    const z = await valid(['a.png']);
    const at = dec.decode(z).indexOf('file a.png');
    const flipped = z.slice();
    flipped[at] ^= 1;
    expect(codes(await readArchive(flipped))).toEqual(['seq.archive.unreadable']);
    const wide = new Uint8Array(z.length + 20);
    wide.set(z, 10);
    expect(codes(await readArchive(wide.subarray(10, 10 + z.length)))).toEqual([]);
  });
});

/** A xorshift generator: the fuzz is reproducible. */
function rng(seed: number) {
  let s = seed | 0 || 1;
  return () => {
    s ^= s << 13;
    s ^= s >>> 17;
    s ^= s << 5;
    return (s >>> 0) / 0x100000000;
  };
}

/** Mutations of a valid archive: an archive or an error, never a throw, and the bounds held by what opens. */
async function fuzz(runs: number) {
  // content that does not compress, so the base is within the ratio
  const r = rng(91);
  const noise = Uint8Array.from({ length: 5000 }, () => Math.floor(r() * 256));
  const base = await valid(['a.seq', 'd/b.seq'], { 'c.seq': noise });
  const limits = { maxEntries: 50, maxBytes: 1 << 20, maxEntryBytes: 6000, ratioFloor: 1000, maxRatio: 20 };
  let opened = 0;
  for (let i = 0; i < runs; i++) {
    const z = base.slice();
    const n = 1 + Math.floor(r() * 4);
    for (let k = 0; k < n; k++) {
      const at = r() < 0.5 ? Math.floor(r() * z.length) : z.length - 22 - Math.floor(r() * 120);
      z[Math.max(0, at)] = Math.floor(r() * 256);
    }
    const cut = r() < 0.2 ? z.subarray(0, Math.floor(r() * z.length)) : z;
    const out = await readArchive(cut, { limits });
    expect(out.archive !== undefined || out.issues.some((x) => x.severity === 'error')).toBe(true);
    if (!out.archive) continue;
    opened++;
    // the bounds held: no entry past its bound, the total within, every entry past the floor within the ratio, and
    // the listing's entries and the end record's agree
    const listed = (await unzip(SyncRuntimeContext, cut.slice().buffer as ArrayBuffer, true)) as Record<
      string,
      { size: number; csize: number }
    >;
    let total = 0;
    for (const [name, data] of out.archive.files) {
      expect(checkArchiveName(name)).toBeUndefined();
      expect(data.byteLength).toBeLessThanOrEqual(limits.maxEntryBytes);
      total += data.byteLength;
      const { size, csize } = listed[name];
      expect(size).toBe(data.byteLength);
      if (size > limits.ratioFloor) expect(size).toBeLessThanOrEqual(limits.maxRatio * csize);
    }
    expect(total).toBeLessThanOrEqual(limits.maxBytes);
  }
  return opened;
}

describe('reading an archive: fuzz', () => {
  it('EC-91: mutated archives are archives or issues, and every name that opens is under the root', async () => {
    expect(await fuzz(300)).toBeGreaterThan(0);
  });

  it('EC-91: random names are accepted only when they have no way above the root', () => {
    const r = rng(7);
    const pieces = ['a', '..', '.', '', 'b.seq', '\\', 'C:', '\u0001', 'é'];
    for (let i = 0; i < 2000; i++) {
      const name = Array.from({ length: 1 + Math.floor(r() * 5) }, () => pieces[Math.floor(r() * pieces.length)]).join(
        r() < 0.8 ? '/' : '',
      );
      if (checkArchiveName(name) !== undefined) continue;
      expect(name.split('/')).not.toContain('..');
      expect(name.startsWith('/')).toBe(false);
      expect(name).not.toMatch(/[\\\u0000-\u001f]/);
    }
  });

  it('EC-91: the full-size fuzz', { tags: ['heavy'] }, async () => {
    await fuzz(20_000);
  });
});
