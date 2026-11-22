/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { DefaultOpenLimits, type OpenedSource } from '../../data/protocol.js';
import { GenbankCauses } from '../genbank-issues.js';
import { GenbankRuntime, type GenbankSource, genbankRuntime, parseGenbank } from '../genbank.js';
import { memoryAssets } from './bed-fixtures.js';
import { generated, letters, record } from './genbank-fixtures.js';

// The GenBank runtime's bounds (_spec/README.md "Data", GenBank; security.md: a file is untrusted input): each refused
// before what it bounds is held (the file's: `seq.source.too_large`; a feature's: the feature left out, said), and a
// file costs its own length: the parser's searches stay in their line, counted here as the characters they scan (no
// clock; the wall-clock bounds are `genbank-time.test.ts`, heavy).

afterEach(() => vi.restoreAllMocks());

const signal = new AbortController().signal;

async function outcome(text: string, over = {}, limits = DefaultOpenLimits): Promise<[string[], GenbankSource?]> {
  try {
    const s = (await genbankRuntime(over).open(
      { address: '/g', format: 'genbank', uri: 'f' },
      {
        assets: memoryAssets({ f: text }),
        signal,
        limits,
      },
    )) as GenbankSource;
    return [['opened', ...(s.issues ?? []).map((i) => i.code)], s];
  } catch (e) {
    return [((e as { issues: { code: string }[] }).issues ?? []).map((i) => i.code)];
  }
}

const one = (qualifiers: (readonly [string, string?])[], location = '1..10') =>
  record('b', letters(100), [
    { key: 'misc_feature', location, qualifiers },
    { key: 'gene', location: '20..30' },
  ]);

async function features(s: OpenedSource | undefined) {
  const a = await s!.read({ system: 'b', start: 0, end: 100, level: 1, shapes: ['ranges'] }, signal);
  return a.ok && a.tables.ranges?.form === 'list' ? a.tables.ranges.start.length : -1;
}

describe('bounds of the file: refused before it is held (`seq.source.too_large`)', () => {
  it('its bytes (`limits.maxBytes`), before a line is read', async () => {
    const [codes] = await outcome(generated(10), {}, { ...DefaultOpenLimits, maxBytes: 1000 });
    expect(codes).toEqual(['seq.source.too_large']);
  });

  it('its records (`maxRecords`), counted as they end', async () => {
    const text = Array.from({ length: 4 }, (_, k) => record(`r${k}`, 'acgt', [])).join('');
    expect((await outcome(text, { maxRecords: 3 }))[0]).toEqual(['seq.source.too_large']);
    expect((await outcome(text, { maxRecords: 4 }))[0][0]).toBe('opened');
  });

  it('its features (`limits.maxRows`), counted as they end', async () => {
    const limits = { ...DefaultOpenLimits, maxRows: 99 };
    expect((await outcome(generated(100), {}, limits))[0]).toEqual(['seq.source.too_large']);
    expect((await outcome(generated(99), {}, limits))[0][0]).toBe('opened');
  });
});

describe('bounds of a feature: the feature left out, said; the file read', () => {
  it('its qualifiers (`maxQualifiers`)', async () => {
    const qs = Array.from({ length: 5 }, (_, k) => [`q${k}`, '"v"'] as const);
    const [codes, s] = await outcome(one(qs), { maxQualifiers: 4 });
    expect(codes).toEqual(['opened', 'seq.source.genbank.feature']);
    expect(await features(s)).toBe(1);
    expect(await features((await outcome(one(qs), { maxQualifiers: 5 }))[1])).toBe(2);
  });

  it('a qualifier’s text (`maxQualifierBytes`): one line, many lines, and a repeated name together', async () => {
    const long = `"${'x'.repeat(2000)}"`;
    const wrapped = `"${Array.from({ length: 100 }, () => 'y'.repeat(30)).join('\n                     ')}"`;
    for (const qs of [
      [['note', long]],
      [['note', wrapped]],
      Array.from({ length: 80 }, () => ['note', '"zzzzzzzzzzzzzzz"'] as const),
    ]) {
      const [codes, s] = await outcome(one(qs as never), { maxQualifierBytes: 1000 });
      expect(codes).toEqual(['opened', 'seq.source.genbank.feature']);
      expect(await features(s)).toBe(1);
    }
  });

  it('a location’s pieces (`maxLocationPieces`) and its nesting (16)', async () => {
    const pieces = `join(${Array.from({ length: 6 }, (_, k) => `${k * 3 + 1}..${k * 3 + 2}`).join(',')})`;
    expect((await outcome(one([], pieces), { maxLocationPieces: 5 }))[0]).toEqual([
      'opened',
      'seq.source.genbank.location',
    ]);
    expect((await outcome(one([], pieces), { maxLocationPieces: 6 }))[0]).toEqual(['opened']);
    const deep = `${'complement('.repeat(17)}1..5${')'.repeat(17)}`;
    expect((await outcome(one([], deep)))[0]).toEqual(['opened', 'seq.source.genbank.location']);
  });

  it('a location’s text past 64 characters a piece is left out before it is held', async () => {
    const wide = `join(${'1..2,'.repeat(200)}1..2)`;
    const [codes] = await outcome(one([], wide), { maxLocationPieces: 10 });
    expect(codes).toEqual(['opened', 'seq.source.genbank.location']);
  });

  it('a quote never closed: the indented lines after it its text, its feature left out at the next line from the first column, the letters read', async () => {
    const text = one([['note', '"never closed']]);
    const [codes, s] = await outcome(text);
    expect(codes).toEqual(['opened', 'seq.source.genbank.feature']);
    expect(s!.coverage.systems[0].length).toBe(100);
    expect(await features(s)).toBe(0);
  });

  it('`/translation` is never held: one of 100,000 characters passes a 1,000-character bound', async () => {
    const t = `"${Array.from({ length: 2000 }, () => 'M'.repeat(50)).join('\n                     ')}"`;
    const [codes, s] = await outcome(one([['translation', t]]), { maxQualifierBytes: 1000 });
    expect(codes).toEqual(['opened', 'seq.source.genbank.translation']);
    expect(await features(s)).toBe(2);
  });
});

/** The characters a search scans: from its start to what it found (inclusive), or to the end. */
const scanned = (length: number, found: number, from = 0) => (found < 0 ? length : found + 1) - Math.max(0, from);

describe('a file costs its own length', () => {
  const parse = (text: string) =>
    parseGenbank(text, new GenbankCauses(), {
      maxRecords: 1000,
      maxFeatures: 1e6,
      maxQualifiers: 256,
      maxQualifierBytes: 65536,
      maxLocationPieces: 10_000,
    });

  it('`indexOf` scans at most the text once: lines of qualifiers without `=`, quotes left open, long ORIGINs', () => {
    const flags = Array.from({ length: 3000 }, (_, k) => `                     /flag${k}`).join('\n');
    const open = `                     /note="${Array.from({ length: 3000 }, () => 'w /x=y').join('\n                     ')}`;
    const text = [
      'LOCUS       s   600 bp    DNA     linear',
      'FEATURES             Location/Qualifiers',
      '     misc_feature    1..10',
      flags,
      '     misc_feature    1..10',
      open,
      ...[
        'ORIGIN',
        ...Array.from({ length: 10 }, (_, k) => `${String(k * 60 + 1).padStart(9)} ${'acgtacgtac '.repeat(6)}`),
      ],
      '//',
      '',
    ].join('\n');
    const indexOf = String.prototype.indexOf;
    let total = 0;
    vi.spyOn(String.prototype, 'indexOf').mockImplementation(function (this: string, v: string, from?: number) {
      const r = indexOf.call(this, v, from);
      if (this.length === text.length) total += scanned(this.length, r, from);
      return r;
    });
    const f = parse(text);
    vi.restoreAllMocks();
    expect(f.ok && f.records[0].letters?.length).toBe(600);
    expect(total).toBeLessThanOrEqual(text.length);
  });

  it('a DEFINITION of 100,000 lines stops growing at its bound (each line would copy it again)', () => {
    const lines = Array.from({ length: 100_000 }, () => `            ${'d'.repeat(60)}`).join('\n');
    const text = `LOCUS       s   4 bp    DNA     linear\nDEFINITION  x\n${lines}\nORIGIN\n        1 acgt\n//\n`;
    const slice = String.prototype.slice;
    let copied = 0;
    vi.spyOn(String.prototype, 'slice').mockImplementation(function (this: string, a?: number, b?: number) {
      const r = slice.call(this, a, b);
      if (this.length > 1000 && this.length !== text.length) copied += r.length;
      return r;
    });
    const f = parse(text);
    vi.restoreAllMocks();
    expect(f.ok && f.records[0].definition?.length).toBe(65536);
    // the definition copied while it grows to the bound, then never again: about bound² / (2 × 61), not lines × bound
    expect(copied).toBeLessThan(65536 * 1200);
  });
});

describe('the runtime with its defaults', () => {
  it('reads the default file of every form', async () => {
    const s = await GenbankRuntime.open(
      { address: '/g', format: 'genbank', value: generated(5) as never },
      {
        assets: memoryAssets({}),
        signal,
        limits: DefaultOpenLimits,
      },
    );
    expect(s.coverage.systems[0].name).toBe('gen');
  });
});
