/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { gzipSync } from 'node:zlib';
import { columnValues } from '@molstar/core/repr/tree';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { checkCoverage } from '../../data/seq-format.js';
import { checkAnswer } from '../../data/tables.js';
import { DefaultOpenLimits, type OpenedSource, type RangesList, type RegionAnswer } from '../../data/protocol.js';
import { GenbankRuntime, type GenbankSource, genbankRuntime, parseLocation } from '../genbank.js';
import { memoryAssets } from './bed-fixtures.js';
import { generated, letters, plasmidGenbank, record } from './genbank-fixtures.js';

// The GenBank runtime (_spec/README.md "Data", GenBank): the location grammar (every form, nesting, partial ends, a
// join over a circular record's origin, each refusal per feature); records as sequences (letters, length, cyclic,
// alphabet, aliases) with their features as ranges; multi-record files; the name index (EC-228); the examples'
// generated plasmid record; fuzzed and truncated files without a throw (EC-40, EC-92); issues one a cause with its
// count and places (EC-86); ids `<key>@<first>-<last>` with `#2` on a repeat (EC-196); the format's base fixed (EC-195);
// a location past the record's length refused, never clipped into a wrong place (EC-51).

const signal = new AbortController().signal;
const limits = DefaultOpenLimits;

async function open(text: string | Uint8Array, extra: Record<string, unknown> = {}, runtime = GenbankRuntime) {
  return (await runtime.open({ address: '/data/gb', format: 'genbank', uri: 'f', ...extra } as never, {
    assets: memoryAssets({ f: text }),
    signal,
    limits,
  })) as GenbankSource;
}

async function read(s: OpenedSource, system: string, start: number, end: number, shapes = ['letters', 'ranges']) {
  return s.read({ system, start, end, level: 1, shapes: shapes as never }, signal);
}

function list(a: RegionAnswer): RangesList {
  if (!a.ok || a.tables.ranges?.form !== 'list') throw new Error('no list');
  return a.tables.ranges;
}

const field = (r: RangesList, name: string) => columnValues(r.fields!.find((f) => f.name === name)!.column);
const codes = (s: OpenedSource) => (s.issues ?? []).map((i) => i.code);
/** The codes of an open's refusal, or of the opened source's issues. */
async function outcome(p: Promise<OpenedSource>): Promise<string[]> {
  try {
    const s = await p;
    return ['opened', ...codes(s)];
  } catch (e) {
    const issues = (e as { issues?: { code: string }[] }).issues;
    if (!issues) throw e;
    return issues.map((i) => i.code);
  }
}

describe('the location grammar', () => {
  const at = (text: string, length = 1000, circular = false) => parseLocation(text, length, circular, 10_000);

  it('a range, a single base, `a.b` as its hull; 1-based inclusive to 0-based half-open (EC-195)', () => {
    expect(at('10..20')).toMatchObject({ ok: true, pieces: [9, 20], start: 9, end: 20, strand: 1 });
    expect(at('7')).toMatchObject({ ok: true, pieces: [6, 7] });
    expect(at('102.110')).toMatchObject({ ok: true, pieces: [101, 110], approximate: true });
  });

  it('`complement`, `join`, `order` and their nesting as INSDC writes it', () => {
    expect(at('complement(10..20)')).toMatchObject({ ok: true, pieces: [9, 20], strand: -1 });
    expect(at('join(1..5,10..20)')).toMatchObject({ ok: true, pieces: [0, 5, 9, 20], start: 0, end: 20, strand: 1 });
    expect(at('order(30..40,1..5)')).toMatchObject({ ok: true, pieces: [0, 5, 29, 40] });
    expect(at('complement(join(1..5,10..20))')).toMatchObject({ ok: true, pieces: [0, 5, 9, 20], strand: -1 });
    expect(at('join(complement(10..20),complement(1..5))')).toMatchObject({
      ok: true,
      pieces: [0, 5, 9, 20],
      strand: -1,
    });
    expect(at('join(complement(1..5),10..20)')).toMatchObject({ ok: true, strand: 0 });
    expect(at('complement(complement(1..5))')).toMatchObject({ ok: true, strand: 1 });
    // pieces overlapping or touching are one (a slippage's shared base)
    expect(at('join(1..10,10..20)')).toMatchObject({ ok: true, pieces: [0, 20] });
  });

  it('partial ends kept: `<` and `>`', () => {
    expect(at('<1..200')).toMatchObject({ ok: true, partialStart: true, partialEnd: false });
    expect(at('1..>200')).toMatchObject({ ok: true, partialStart: false, partialEnd: true });
    expect(at('complement(<5..>50)')).toMatchObject({ ok: true, partialStart: true, partialEnd: true, strand: -1 });
    expect(at('join(<1..10,20..>30)')).toMatchObject({ ok: true, partialStart: true, partialEnd: true });
  });

  it('a join over a circular record’s origin is one item across it, both strands', () => {
    expect(at('join(2600..2686,1..120)', 2686, true)).toMatchObject({
      ok: true,
      pieces: [2599, 2806],
      start: 2599,
      end: 2806,
      wraps: true,
    });
    for (const t of ['complement(join(2600..2686,1..120))', 'join(complement(1..120),complement(2600..2686))'])
      expect(at(t, 2686, true)).toMatchObject({ ok: true, pieces: [2599, 2806], strand: -1, wraps: true });
    // pieces apart over the origin stay apart
    expect(at('join(2600..2680,10..120)', 2686, true)).toMatchObject({ ok: true, pieces: [2599, 2680, 2695, 2806] });
    // on a linear record the same text is two pieces where they are
    expect(at('join(2600..2686,1..120)', 2686, false)).toMatchObject({ ok: true, pieces: [0, 120, 2599, 2686] });
    // more than one lap is refused
    expect(at('join(2600..2686,1..2686,1..10)', 2686, true)).toMatchObject({ ok: false });
  });

  it('refused per feature: remote references, `a^b` sites, past the length (EC-51), malformed, too deep, too many pieces', () => {
    const why = (t: string, length = 1000) => {
      const r = at(t, length);
      return r.ok ? 'read' : r.why;
    };
    expect(why('J00194.1:100..202')).toMatch(/remote/);
    expect(why('join(1..10,J00194.1:100..202)')).toMatch(/remote/);
    expect(why('123^124')).toMatch(/between two bases/);
    expect(why('1..1001')).toMatch(/past the record/);
    expect(why('20..10')).toMatch(/before its start/);
    expect(why('0..10')).toMatch(/before 1/);
    expect(why('99999999999..1')).toMatch(/2\^31/);
    for (const bad of ['', 'join(1..5', 'join(1..5,)', 'complement(1..5', '1..5)', 'gap(100)', '1..', '..5', 'join()'])
      expect(why(bad), bad).not.toBe('read');
    const deep = `${'complement('.repeat(17)}1..5${')'.repeat(17)}`;
    expect(why(deep)).toMatch(/deeper than 16/);
    expect(why(`${'complement('.repeat(15)}1..5${')'.repeat(15)}`)).toBe('read');
    const many = `join(${Array.from({ length: 11 }, (_, k) => `${k * 2 + 1}..${k * 2 + 1}`).join(',')})`;
    expect(parseLocation(many, 1000, false, 10)).toMatchObject({ ok: false });
    expect(parseLocation(many, 1000, false, 11)).toMatchObject({ ok: true });
  });

  it('a nesting 20,000 deep is refused at 16, no stack overflow', () => {
    const t = `${'join('.repeat(20_000)}1..5${')'.repeat(20_000)}`;
    expect(at(t)).toMatchObject({ ok: false });
  });
});

describe('records and their features', () => {
  const text = record(
    'demo',
    letters(300),
    [
      {
        key: 'gene',
        location: '10..100',
        qualifiers: [
          ['gene', '"abcA"'],
          ['locus_tag', '"DEM_1"'],
        ],
      },
      {
        key: 'CDS',
        location: 'complement(join(10..40,60..100))',
        qualifiers: [
          ['gene', '"abcA"'],
          ['product', '"a protein, with a ""quoted"" word"'],
          ['note', '"one"'],
          ['note', '"two"'],
          ['translation', '"MKV\n                     LLA"'],
          ['pseudo'],
        ],
      },
      { key: 'misc_feature', location: '<1..>20', qualifiers: [['label', 'site']] },
      { key: 'misc_feature', location: '<1..>20' },
      { key: 'repeat_region', location: '200^201' },
      { key: 'gene', location: 'J00194.1:100..202' },
      { key: 'gene', location: '250..400' },
    ],
    { accession: 'AB000001' },
  );

  it('a record is a sequence: LOCUS’s name, its letters’ length, its alphabet, its accession an alias', async () => {
    const s = await open(text);
    expect(checkCoverage(s.coverage)).toBeUndefined();
    expect(s.coverage.systems).toEqual([{ name: 'demo', length: 300, alphabet: 'dna', aliases: ['AB000001.1'] }]);
    expect(s.coverage.shapes.map((x) => x.shape)).toEqual(['letters', 'ranges']);
    expect(s.records).toEqual([
      {
        name: 'demo',
        accession: 'AB000001.1',
        definition: 'A generated record demo.',
        organism: 'synthetic construct',
        molecule: 'DNA',
        length: 300,
        circular: false,
        features: 4,
      },
    ]);
    const a = await read(s, 'demo', 5, 15, ['letters']);
    expect(a.ok && a.tables.letters).toEqual({ form: 'letters', start: 5, letters: letters(300).slice(5, 15) });
  });

  it('features as ranges: key the category, label by its rule, qualifiers fields, strand, partial ends, ids (EC-196)', async () => {
    const s = await open(text);
    const a = await read(s, 'demo', 0, 300);
    const request = { system: 'demo', start: 0, end: 300, level: 1, shapes: ['letters', 'ranges'] } as const;
    expect(checkAnswer(a, request)).toBeUndefined();
    const r = list(a);
    const ids = new Set(columnValues(r.id!));
    expect(ids).toEqual(new Set(['gene@10-100', 'CDS@10-100', 'misc_feature@1-20', 'misc_feature@1-20#2']));
    const byId = (id: string, name: string) => field(r, name)[columnValues(r.id!).indexOf(id)];
    expect(byId('gene@10-100', 'label')).toBe('abcA');
    expect(byId('gene@10-100', 'category')).toBe('gene');
    expect(byId('misc_feature@1-20', 'label')).toBe('site');
    expect(byId('misc_feature@1-20#2', 'label')).toBe('misc_feature');
    expect(byId('misc_feature@1-20', 'partial_start')).toBe(1);
    expect(byId('misc_feature@1-20', 'partial_end')).toBe(1);
    expect(byId('gene@10-100', 'partial_end')).toBe(0);
    expect(byId('CDS@10-100', 'product')).toBe('a protein, with a "quoted" word');
    expect(byId('CDS@10-100', 'note')).toBe('one; two');
    expect(byId('CDS@10-100', 'pseudo')).toBe('true');
    expect(r.fields!.map((f) => f.name)).not.toContain('translation');
    // the CDS's two pieces, its extent on each, its strand
    const cds = [...columnValues(r.id!).keys()].filter((k) => columnValues(r.id!)[k] === 'CDS@10-100');
    expect(cds.map((k) => [r.start[k], r.end[k], r.itemStart![k], r.itemEnd![k], r.strand![k]])).toEqual([
      [9, 40, 9, 100, -1],
      [59, 100, 9, 100, -1],
    ]);
  });

  it('each refusal said once a cause with its places (EC-86), never the file; `/translation` counted', async () => {
    const s = await open(text);
    const issues = s.issues ?? [];
    const location = issues.filter((i) => i.code === 'seq.source.genbank.location');
    expect(location.map((i) => i.message).sort()).toEqual(
      [
        expect.stringMatching(/a position past the record’s length.*: 1 \(line \d+\)/),
        expect.stringMatching(/a remote reference.*: 1 \(line \d+\)/),
        expect.stringMatching(/a site between two bases.*: 1 \(line \d+\)/),
      ].sort((a, b) => String(a).localeCompare(String(b))) as never,
    );
    expect(codes(s)).toContain('seq.source.genbank.translation');
    expect(issues.find((i) => i.code === 'seq.source.genbank.translation')!.severity).toBe('info');
  });

  it('`base` and `ends` not read, said (EC-195)', async () => {
    const s = await open(text, { base: 0, ends: 'exclusive' });
    expect(codes(s)).toContain('seq.source.genbank.transform');
    expect(list(await read(s, 'demo', 0, 300)).start).toContain(9);
  });

  it('header fields not read counted once; letters that disagree with LOCUS: the letters win, said', async () => {
    const t = record('demo', letters(120), [], { length: 100 }).replace(
      'SOURCE',
      'KEYWORDS    a.\nCOMMENT     c\nSOURCE',
    );
    const s = await open(t);
    expect(s.coverage.systems[0].length).toBe(120);
    const header = (s.issues ?? []).filter((i) => i.code === 'seq.source.genbank.header');
    expect(header).toHaveLength(1);
    expect(header[0].message).toMatch(/: 2 \(KEYWORDS, COMMENT\)/);
    expect(codes(s)).toContain('seq.source.genbank.length');
  });
});

describe('circular records, molecules and several records', () => {
  it('`circular` a cyclic system; a join over its origin one piece past the length, found from either end', async () => {
    const t = record('ring', letters(500), [{ key: 'CDS', location: 'complement(join(450..500,1..30))' }], {
      circular: true,
    });
    const s = await open(t);
    expect(s.coverage.systems[0]).toMatchObject({ name: 'ring', length: 500, cyclic: true });
    const r = list(await read(s, 'ring', 460, 470));
    expect([r.start[0], r.end[0], r.itemStart![0], r.itemEnd![0], r.strand![0]]).toEqual([449, 530, 449, 530, -1]);
    expect(columnValues(r.id!)).toEqual(['CDS@450-530']);
  });

  it('the molecule its alphabet: RNA, mRNA, a protein’s `aa` (strand none); a protein’s letters past its alphabet refuse the file', async () => {
    const rna = await open(record('r', 'acgu', [], { molecule: 'mRNA' }));
    expect(rna.coverage.systems[0].alphabet).toBe('rna');
    expect(rna.records[0].molecule).toBe('mRNA');
    const p = await open(record('p', 'mkvlaa', [{ key: 'Region', location: '2..4' }], { protein: true }));
    expect(p.coverage.systems[0].alphabet).toBe('protein');
    expect(list(await read(p, 'p', 0, 6)).strand![0]).toBe(0);
    expect(await outcome(open(record('p', 'mk1la', [], { protein: true }).replace('mk1la', 'mk%la')))).toEqual([
      'seq.source.malformed',
    ]);
  });

  it('several records: a sequence each, features on their own; a repeated name refused; a record without ORIGIN ranges only', async () => {
    const t = [
      record('one', letters(100, 1), [{ key: 'gene', location: '1..10', qualifiers: [['gene', '"x"']] }]),
      record('two', letters(200, 2), [{ key: 'gene', location: '5..50', qualifiers: [['gene', '"x"']] }]),
      record('one', letters(50, 3), []),
      record('bare', '', [{ key: 'gene', location: '1..10' }], { length: 900 }).replace('ORIGIN\n', ''),
    ].join('');
    const s = await open(t);
    expect(s.coverage.systems.map((x) => [x.name, x.length])).toEqual([
      ['one', 100],
      ['two', 200],
      ['bare', 900],
    ]);
    expect(s.coverage.shapes.filter((x) => x.system === 'bare').map((x) => x.shape)).toEqual(['ranges']);
    expect(codes(s)).toContain('seq.source.genbank.record');
    expect(columnValues(list(await read(s, 'two', 0, 200)).id!)).toEqual(['gene@5-50']);
    const hits = await s.findNames!('x', { max: 10 }, signal);
    expect(hits.map((h) => [h.system, h.id])).toEqual([
      ['one', 'gene@1-10'],
      ['two', 'gene@5-50'],
    ]);
  });

  it('a CONTIG assembly of other records is not read, said; the records beside it are', async () => {
    const contig = record('asm', '', [], { length: 5000 }).replace('ORIGIN\n', 'CONTIG      join(A1.1:1..5000)\n');
    const s = await open(`${contig}${record('ok', letters(10), [])}`);
    expect(s.coverage.systems.map((x) => x.name)).toEqual(['ok']);
    expect((s.issues ?? []).find((i) => i.code === 'seq.source.genbank.record')!.message).toMatch(/CONTIG/);
  });
});

describe('refused, the file', () => {
  it('no LOCUS line; a file cut short before any `//`; a truncated last record said, the others kept', async () => {
    expect(await outcome(open('>fasta\nacgt\n'))).toEqual(['seq.source.malformed']);
    const one = record('one', letters(100), []);
    expect(await outcome(open(one.slice(0, one.length - 4)))).toEqual(['seq.source.malformed']);
    const two = `${one}${record('two', letters(100), [])}`;
    const s = await open(two.slice(0, two.length - 4));
    expect(s.coverage.systems.map((x) => x.name)).toEqual(['one']);
    expect((s.issues ?? []).find((i) => i.code === 'seq.source.genbank.record')!.message).toMatch(/cut short/);
  });

  it('gzip read by its magic, never its extension', async () => {
    const s = await open(new Uint8Array(gzipSync(record('z', letters(30), []))));
    expect(s.coverage.systems[0].name).toBe('z');
  });
});

describe('the name index (EC-228)', () => {
  it('`/label`, `/gene`, `/locus_tag`, the LOCUS name and the accession; exact first, then by prefix', async () => {
    const t = record(
      'pX',
      letters(200),
      [
        {
          key: 'CDS',
          location: '1..30',
          qualifiers: [
            ['gene', '"lacZ"'],
            ['label', '"lacZ alpha"'],
          ],
        },
        {
          key: 'CDS',
          location: '40..90',
          qualifiers: [
            ['gene', '"lacI"'],
            ['locus_tag', '"PX_002"'],
          ],
        },
      ],
      { accession: 'ZZ123' },
    );
    const s = await open(t);
    expect(s.coverage.names).toEqual(['label', 'gene', 'locus_tag', 'locus', 'accession']);
    const find = async (q: string, prefix = false) =>
      (await s.findNames!(q, { max: 20, prefix }, signal)).map((h) => [h.field, h.id ?? '', h.label, h.start, h.end]);
    expect(await find('LACZ')).toEqual([['gene', 'CDS@1-30', 'lacZ alpha', 0, 30]]);
    expect(await find('lac', true)).toEqual([
      ['gene', 'CDS@40-90', 'lacI', 39, 90],
      ['gene', 'CDS@1-30', 'lacZ alpha', 0, 30],
      ['label', 'CDS@1-30', 'lacZ alpha', 0, 30],
    ]);
    expect(await find('px_002')).toEqual([['locus_tag', 'CDS@40-90', 'lacI', 39, 90]]);
    expect(await find('px')).toEqual([['locus', '', 'pX', 0, 200]]);
    expect(await find('zz123.1')).toEqual([['accession', '', 'pX', 0, 200]]);
    expect(await s.findNames!('lac', { max: 1, prefix: true }, signal)).toHaveLength(1);
  });
});

describe('the examples’ generated plasmid record', () => {
  it('pUC19c turned: cyclic, 2,686 letters as read (lowercase), bla across the origin, lacZα complemented', async () => {
    const s = await open(plasmidGenbank());
    expect(s.coverage.systems).toEqual([{ name: 'pUC19c', length: 2686, alphabet: 'dna', cyclic: true }]);
    expect(s.records[0]).toMatchObject({ organism: 'Cloning vector pUC19c', circular: true, features: 6 });
    expect(codes(s).sort()).toEqual(['seq.source.genbank.header', 'seq.source.genbank.translation']);
    const a = await read(s, 'pUC19c', 0, 2686);
    const letters = a.ok ? (a.tables.letters?.letters as string) : '';
    expect(letters).toMatch(/^[acgt]{2686}$/);
    const r = list(a);
    const ids = columnValues(r.id!);
    const row = (label: string) => ids.indexOf(ids[field(r, 'label').indexOf(label)]);
    const bla = row('bla (AmpR)');
    expect([r.start[bla], r.end[bla], r.strand![bla], ids[bla]]).toEqual([2311, 3172, -1, 'CDS@2312-3172']);
    const lacZ = row('lacZα');
    expect([r.start[lacZ], r.end[lacZ], r.strand![lacZ]]).toEqual([831, 1155, -1]);
    const primer = row('primer site');
    expect([r.start[primer], r.end[primer], r.strand![primer]]).toEqual([2656, 2711, 1]);
    expect((await s.findNames!('bla', { max: 5 }, signal)).map((h) => h.id)).toEqual(['CDS@2312-3172']);
  });
});

describe('fuzzed and truncated files (EC-40, EC-92)', () => {
  it('bytes flipped, files cut: every outcome a source whose answers check, or a refusal; never a throw', async () => {
    const rng = new PCG(11);
    const base = new TextEncoder().encode(`${generated(40)}${plasmidGenbank()}`);
    for (let k = 0; k < 60; k++) {
      const b = base.slice(0, Math.floor(rng.float() * base.length));
      for (let m = 0; m < 30; m++) b[Math.floor(rng.float() * b.length)] = Math.floor(rng.float() * 256);
      // structural characters too: quotes, parentheses, slashes, newlines
      for (let m = 0; m < 30; m++) b[Math.floor(rng.float() * b.length)] = '"()/\n.,<>^:'.charCodeAt(m % 11);
      let s: OpenedSource;
      try {
        s = await open(b);
      } catch (e) {
        expect((e as { issues?: unknown[] }).issues?.length).toBeGreaterThan(0);
        continue;
      }
      expect(checkCoverage(s.coverage)).toBeUndefined();
      for (const sys of s.coverage.systems) {
        const request = {
          system: sys.name,
          start: 0,
          end: sys.length ?? 1000,
          level: 1,
          shapes: ['letters', 'ranges'],
        };
        const a = await s.read(request as never, signal);
        expect(a.ok).toBe(true);
        expect(checkAnswer(a, request as never)).toBeUndefined();
      }
    }
  });
});

describe('size', () => {
  // the quick representative of 10,000 features under 100 ms (`genbank-time.test.ts`, heavy): no clock here
  it('1,000 features parsed, each an item', async () => {
    const s = await open(generated(1000));
    const r = list(await read(s, 'gen', 0, s.coverage.systems[0].length!, ['ranges']));
    expect(new Set(columnValues(r.id!)).size).toBe(1000);
  });

  it('options bound what they name', async () => {
    const rt = genbankRuntime({ maxRecords: 1 });
    expect(await outcome(open(`${record('a', 'ac', [])}${record('b', 'ac', [])}`, {}, rt))).toEqual([
      'seq.source.too_large',
    ]);
  });
});
