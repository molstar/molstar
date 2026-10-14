/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tabix-indexed VCF against the recorded slices of ClinVar (archive 2024-01-07, GRCh38, chromosomes named 1..22, X, Y,
// MT) and its .tbi. The expected variants are lines found without the index: the first and the last 2 to 3 megabytes of
// the file inflated with Node's gunzip and filtered by chromosome and position (scripts/record.ts).

import { readFileSync } from 'node:fs';
import { StreamNotServed } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { vcfProvider } from '../../../src/stream/gmod/index.js';
import { CASES, LENGTHS, makers, req } from './cases.js';
import { sliceReader, whole } from './slices.js';

interface Variant {
  start: number;
  ref: string;
  alt: string;
  name: string;
  clnsig: string;
}
const expected = JSON.parse(readFileSync(new URL('./fixtures/expected.json', import.meta.url), 'utf8')).vcf as Record<
  string,
  Variant[]
>;
const signal = new AbortController().signal;
const tbi = await whole('tbi');
const provider = (maxItems?: number) => makers.vcf(sliceReader('vcf'), tbi, maxItems);

describe('tabix VCF', () => {
  it('lists the indexed chromosomes, with the lengths the host gives', async () => {
    const p = await provider().fetch(req(null), signal);
    if (p.kind !== 'manifest') throw new Error('not a manifest');
    const sizes = Object.fromEntries(p.manifest.systems.map((s) => [s.key, s.extent[1]]));
    expect([sizes[1], sizes.MT]).toEqual([LENGTHS[1], LENGTHS.MT]);
    expect(p.manifest.systems.map((s) => s.key).slice(0, 3)).toEqual(['1', '2', '3']);
    expect(p.manifest.systems.map((s) => s.key)).toContain('NW_009646201.1');
    expect(p.manifest.levels).toEqual([256, 4096, 65536]);
    expect(p.manifest.vocabulary).toEqual({ category: ['SNV', 'MNV', 'insertion', 'deletion', 'structural'] });
  });

  it('without lengths from the host, a chromosome ends where its linear index does', async () => {
    const p = await vcfProvider(sliceReader('vcf'), tbi).fetch(req(null), signal);
    if (p.kind !== 'manifest') throw new Error('not a manifest');
    const sizes = Object.fromEntries(p.manifest.systems.map((s) => [s.key, s.extent[1]]));
    expect(sizes[1]).toBeGreaterThan(LENGTHS[1] - 100_000);
    expect(sizes[1]).toBeLessThanOrEqual(LENGTHS[1]);
    expect(sizes[1] % 16_384).toBe(0);
  });

  for (const [name, r] of CASES.vcf.slice(1))
    it(`${name}: the variants are those found in the file without the index`, async () => {
      const p = await provider().fetch(r, signal);
      if (p.kind !== 'items') throw new Error('not items');
      const t = p.items as Record<string, (string | number)[]>;
      const got = t.id.map((_, j) => ({
        start: t.start[j],
        ref: t.ref[j],
        alt: t.alt[j],
        name: t.name[j],
        clnsig: t.CLNSIG[j],
      }));
      const key = (v: { start: unknown; ref: unknown; alt: unknown; name: unknown }) =>
        `${v.start}:${v.ref}:${v.alt}:${v.name}`;
      expect(got.map(key).sort()).toEqual(expected[name].map(key).sort());
      expect(got.map((v) => v.clnsig).sort()).toEqual(expected[name].map((v) => v.clnsig).sort());
      expect(got.length).toBeGreaterThan(100);
      expect(t.stop.every((stop, j) => stop === (t.start[j] as number) + (t.ref[j] as string).length)).toBe(true);
      expect(new Set(t.id).size).toBe(t.id.length); // an id is one variant: the viewer joins rows of an id
      expect(t.label.every((l) => l !== '')).toBe(true);
    });

  it('types a variant by the lengths of ref and alt', async () => {
    const p = await provider().fetch(CASES.vcf[1][1], signal);
    if (p.kind !== 'items') throw new Error('not items');
    const t = p.items as Record<string, string[]>;
    const type = (j: number) =>
      t.ref[j].length === t.alt[j].length
        ? t.ref[j].length === 1
          ? 'SNV'
          : 'MNV'
        : t.alt[j].length > t.ref[j].length
          ? 'insertion'
          : 'deletion';
    expect(t.category).toEqual(t.id.map((_, j) => type(j)));
    expect(new Set(t.category).size).toBeGreaterThan(1);
  });

  it('answers a data request estimated past maxItems from the index alone (the zoom gate; ST4b); `force` reads it', async () => {
    const [, r] = CASES.vcf[1];
    let reads = 0;
    const base = sliceReader('vcf');
    const counted = (o: number, l: number, sg?: AbortSignal) => {
      reads++;
      return base(o, l, sg);
    };
    const p = makers.vcf(counted, tbi, 100);
    await p.fetch(req(null), signal);
    const before = reads;
    const got = await p.fetch(r, signal); // 908 variants over 25 kb
    if (got.kind !== 'bins') throw new Error('not bins');
    expect(got.estimate).toBe(true);
    const count = got.bins.count as ArrayLike<number>;
    expect([got.level, count.length]).toEqual([0, 98]);
    // No data read: the counts are the index's estimate, about the 908 variants.
    expect(reads).toBe(before);
    const sum = Array.from(count).reduce((a, b) => a + b, 0);
    expect(sum).toBeGreaterThan(908 * 0.6);
    expect(sum).toBeLessThan(908 * 1.4);
    const whole = await p.fetch({ ...r, force: true }, signal);
    expect(whole.kind).toBe('items');
    if (whole.kind !== 'items') throw new Error('not items');
    expect(whole.items.id.length).toBe(908);
  });

  it("gates a region whose chunks pass the request's `max_bytes`, before reading it (ST4b)", async () => {
    const [, r] = CASES.vcf[1];
    const p = provider();
    const small = await p.fetch({ ...r, max_bytes: 1000 }, signal);
    expect(small.kind).toBe('bins');
    const big = await p.fetch({ ...r, max_bytes: 1 << 20 }, signal);
    expect(big.kind).toBe('items');
  });

  it('fails clearly for what it does not serve', async () => {
    const p = provider();
    await expect(p.fetch(req('1', 0, 256, 3), signal)).rejects.toBeInstanceOf(StreamNotServed);
    await expect(p.fetch(req('chr1', 0, 256), signal)).rejects.toThrow('VCF: no chromosome chr1 in the index');
  });

  it('fails naming the range when a region is not in the recorded slices', async () => {
    await expect(provider().fetch({ ...req('17', 43_044_000, 43_060_000), force: true }, signal)).rejects.toThrow(
      /vcf: bytes \d+\.\.\d+ are not in the recorded slices/,
    );
  });
});

describe("a VCF under the document's names (`renamed`, `chrNames`; ST4b)", () => {
  it("lists `1` as `chr1` and `MT` as `chrM`, leaves contigs out, and answers `chr1` as the file's `1`", async () => {
    const { chrNames, renamed } = await import('../../../src/stream/gmod/index.js');
    const p = renamed(provider(), chrNames);
    const m = await p.fetch(req(null), signal);
    if (m.kind !== 'manifest') throw new Error('not a manifest');
    const keys = m.manifest.systems.map((s) => s.key);
    expect(keys.slice(0, 3)).toEqual(['chr1', 'chr2', 'chr3']);
    expect(keys).toContain('chrM');
    expect(keys.some((k) => !k.startsWith('chr'))).toBe(false);
    const [, r] = CASES.vcf.slice(1).find(([, x]) => x.system === '1') ?? [];
    if (!r) throw new Error('no case on 1');
    const direct = await provider().fetch(r, signal);
    const mapped = await p.fetch({ ...r, system: 'chr1' }, signal);
    expect(mapped.system).toBe('chr1');
    expect((mapped as { items: unknown }).items).toEqual((direct as { items: unknown }).items);
    await expect(p.fetch({ ...r, system: 'chr_none' }, signal)).rejects.toBeInstanceOf(StreamNotServed);
  });
});
