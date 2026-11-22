/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node, Value } from '@molstar/core/repr/tree';
import { itWithin } from '@molstar/test-harness/timing';
import { beforeAll, describe, expect, it, vi } from 'vitest';
import { assembly, colinearBlocks, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { includePattern } from '@molstar/seqx-core/lang/include-pattern';
import { locationOfResult, resolveAddressIn } from '@molstar/seqx-core/query';
import { resolveView } from '@molstar/seqx-core/view';
import { blockNodes } from '../../view/blocks.js';
import { assemblyNames, includeFilter, lookupName, naturalOrder } from '../assembly-names.js';
import { loadGenome } from '../genome.js';
import { resolveSpec, sourceLengths } from '../model.js';
import { overlayQueryModel } from '../query-model.js';
import { EmptyOverlay } from '../overlay.js';

// Assemblies (_spec/README.md "Resolve", Assemblies): the names rule, the order, include, derived sequences from inline
// rows and a source's coverage, overrides, the issues, the genome forms and a map between assemblies. The genome chunk
// is loaded first (its own test, `genome-lazy.test.ts`, holds the pending state before it lands).

beforeAll(async () => {
  await loadGenome();
});

const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);

const hg38 = [
  { name: 'chr2', length: 2000 },
  { name: 'chr10', length: 1000 },
  { name: 'chr1', length: 3000, aliases: ['NC_000001.11'] },
  { name: 'chrUn_KI270302v1', length: 20 },
  { name: 'chrM', length: 16 },
  { name: 'chrX', length: 1500 },
  { name: 'HLA-A*01:01:01:01', length: 30 },
];

function spec(extra: Record<string, unknown> = {}, body?: (b: ReturnType<typeof createSeqSpec>) => void): Node {
  const b = createSeqSpec();
  b.data.assembly('hg38', { ...assembly('GRCh38', hg38, { aliases: ['hg38'] }), ...extra } as never, (a) =>
    a.sequence('mito', { name: 'MT', cyclic: true, label: 'mitochondrion' } as never),
  );
  body?.(b);
  return b.build();
}

const A = '/data/hg38';

/** A location with plain arrays, to compare. */
const plain = (l: { parts: readonly { system: string; ranges: ArrayLike<number> }[]; rects: readonly unknown[] }) => ({
  parts: l.parts.map((p) => ({ system: p.system, ranges: Array.from(p.ranges) })),
  rects: [...l.rects],
});

describe('the names rule (EC-189)', () => {
  const seqs = [
    { name: 'chr1', length: 10, aliases: ['1x'] },
    { name: 'MT', length: 16569 },
    { name: 'chrM', length: 16571 },
    { name: '2', length: 5, aliases: ['1x'] },
  ];
  const names = assemblyNames(seqs);

  it('reads a name exactly, then its aliases, then with the `chr` prefix added or dropped', () => {
    expect([lookupName(names, 'chr1'), lookupName(names, '1'), lookupName(names, 'chr2')]).toEqual([0, 0, 3]);
    expect(lookupName(assemblyNames(seqs, 'exact'), '1')).toBeUndefined();
  });

  it('never joins `M` and `MT`, and drops an alias that names two sequences (EC-189)', () => {
    expect([lookupName(names, 'MT'), lookupName(names, 'M'), lookupName(names, 'chrMT')]).toEqual([1, 2, 1]);
    expect(lookupName(names, '1x')).toBeUndefined();
    expect(names.dropped).toEqual(['1x']);
    // an alias table joining hg19's chrM with GRCh37's MT (different lengths) is refused
    const joined = assemblyNames(seqs, 'chr', [{ name: 'chrM', aliases: ['MT'] }]);
    expect(joined.dropped).toContain('MT');
    expect(lookupName(joined, 'MT')).toBe(1);
  });

  it('says an alias that is another sequence’s name dropped, whatever the lengths (EC-189)', () => {
    const r = assemblyNames([
      { name: 'chr1', length: 10, aliases: ['x'] },
      { name: 'x', length: 10 },
      { name: 'chr2', aliases: ['chr1'] },
    ]);
    expect(r.dropped).toEqual(['x', 'chr1']);
    expect([lookupName(r, 'x'), lookupName(r, 'chr1')]).toEqual([1, 0]);
  });

  it('never maps an alias to two sequences, over random alias tables (property)', () => {
    const rng = new PCG(7);
    for (let round = 0; round < 50; round++) {
      const n = 2 + (rng.int() % 6);
      const list = Array.from({ length: n }, (_, i) => ({
        name: `s${i}`,
        length: 1 + (rng.int() % 3),
        aliases: Array.from({ length: rng.int() % 3 }, () => `a${rng.int() % 4}`),
      }));
      const r = assemblyNames(list);
      for (let k = 0; k < 4; k++) {
        const claimed = list.filter((s) => s.aliases.includes(`a${k}`)).map((s) => s.name);
        const i = lookupName(r, `a${k}`);
        if (new Set(claimed).size > 1) expect(i).toBeUndefined();
        else if (claimed.length) expect(list[i!].name).toBe(claimed[0]);
      }
    }
  });
});

describe('order and include', () => {
  it('sorts `natural`: numbers, X, Y, W, Z, M, MT, then the rest as listed', () => {
    const names = ['chrUn', 'MT', 'chr10', 'chrY', '2', 'chrM', 'chrX', 'chr1', 'HLA-A'];
    expect(naturalOrder(names).map((i) => names[i])).toEqual([
      'chr1',
      '2',
      'chr10',
      'chrX',
      'chrY',
      'chrM',
      'MT',
      'chrUn',
      'HLA-A',
    ]);
  });

  it('keeps names through the names rule, patterns and a length floor', () => {
    const list = hg38;
    const names = assemblyNames(list);
    const kept = (v: Parameters<typeof includeFilter>[0]) =>
      list.filter((s, i) => includeFilter(v, names)!(i, s)).map((s) => s.name);
    expect(kept(['1', 'X'])).toEqual(['chr1', 'chrX']);
    expect(kept(['^chr[0-9]+$'])).toEqual(['chr2', 'chr10', 'chr1']);
    expect(kept({ min_length: 1500 })).toEqual(['chr2', 'chr1', 'chrX']);
    expect(includeFilter(undefined, names)).toBeUndefined();
  });
});

describe('include patterns', () => {
  it('take the design’s forms and refuse groups, alternation and counted repeats', () => {
    const names = ['chr1', 'chr22', 'chrX', 'chrM', 'chrUn_KI270302v1', 'chr1_KI270706v1_random', 'HLA-A*01:01'];
    const kept = (pattern: string) => names.filter(includePattern(pattern)!);
    expect(kept('^chr[0-9XYM]+$')).toEqual(['chr1', 'chr22', 'chrX', 'chrM']);
    expect(kept('^chrUn_')).toEqual(['chrUn_KI270302v1']);
    expect(kept('^chr\\d+_.*_random$')).toEqual(['chr1_KI270706v1_random']);
    expect(kept('^HLA-A\\*')).toEqual(['HLA-A*01:01']);
    for (const bad of ['^(a+)+$', '^a|b', '^a{2,}', '^(chr)1', 'chr1', '^[z-a]'])
      expect(includePattern(bad), bad).toBeUndefined();
  });

  it('match 100,000 untrusted names in a table, never a backtracking search (no ReDoS)', () => {
    const list = Array.from({ length: 100_000 }, (_, i) => ({ name: `${'a'.repeat(40)}${i}!`, length: 1 }));
    const names = assemblyNames(list.slice(0, 10));
    // a nested quantifier is no pattern (the language check refuses it): read as a name, matching none
    expect(list.filter((x, i) => includeFilter(['^(a+)+$'], names)!(i, x)).length).toBe(0);
    const keep = includeFilter(['^a*a*a*a*a*a*b$'], names)!;
    expect(list.filter((x, i) => keep(i, x)).length).toBe(0);
  });
});

describe('an assembly’s derived sequences', () => {
  it('derives a system each from inline rows, in `natural` order, an override by name (EC-187, EC-188, EC-200)', () => {
    const m = resolveSpec(spec());
    const a = m.assemblies.get(A)!;
    expect(a.names).toEqual(['GRCh38', 'hg38']);
    expect(a.systems).toEqual([
      '/data/hg38/(sequence :name chr1)',
      '/data/hg38/(sequence :name chr2)',
      '/data/hg38/(sequence :name chr10)',
      '/data/hg38/(sequence :name chrX)',
      '/data/hg38/(sequence :name chrM)',
      '/data/hg38/(sequence :name chrUn_KI270302v1)',
      '/data/hg38/(sequence :name HLA-A*01:01:01:01)',
    ]);
    // `MT` reaches chrM through the prefix rule? never: `M` only. The override names `MT`, which this assembly lacks:
    // dropped, and said
    const chrM = m.systems.get('/data/hg38/(sequence :name chrM)')!;
    expect([chrM.length, chrM.cyclic, chrM.alphabet]).toEqual([16, false, 'dna']);
    expect(codes(m.issues)).toEqual(['seq.resolve.assembly.override']);
  });

  it('drops a sequence named twice and says so, the first kept', () => {
    const b = createSeqSpec();
    b.data.assembly(
      'g',
      assembly('toy', [
        { name: 'chr1', length: 10 },
        { name: 'chr1', length: 99 },
      ]),
    );
    const m = resolveSpec(b.build());
    expect(m.assemblies.get('/data/g')!.systems).toEqual(['/data/g/(sequence :name chr1)']);
    expect(m.systems.get('/data/g/(sequence :name chr1)')!.length).toBe(10);
    expect(codes(m.issues)).toEqual(['seq.resolve.assembly.duplicate']);
  });

  it('says an override naming no sequence of the assembly, dropped (`seq.resolve.assembly.override`)', () => {
    const b = createSeqSpec();
    b.data.assembly('hg38', assembly('GRCh38', hg38), (a) => {
      a.sequence('mito', { name: 'M', cyclic: true });
      a.sequence('typo', { name: 'chrQ', label: 'nothing' });
    });
    const m = resolveSpec(b.build());
    expect(codes(m.issues)).toEqual(['seq.resolve.assembly.override']);
    expect(m.issues[0].message).toContain('chrQ');
    expect(m.issues[0].address).toBe('/data/hg38');
    // the others apply, the system list is the derived sequences'
    expect(m.systems.get('/data/hg38/(sequence :name chrM)')!.cyclic).toBe(true);
    expect(m.assemblies.get('/data/hg38')!.systems).toHaveLength(hg38.length);
  });

  it('applies an override found through the names rule', () => {
    const b = createSeqSpec();
    b.data.assembly('hg38', assembly('GRCh38', hg38), (a) =>
      a.sequence('mito', { name: 'M', cyclic: true, label: 'mitochondrion' }),
    );
    const s = resolveSpec(b.build()).systems.get('/data/hg38/(sequence :name chrM)')!;
    expect([s.cyclic, s.label]).toEqual([true, 'mitochondrion']);
  });

  it('keeps the whole assembly in its order before `include`: the `system` domain (EC-191, its model part)', () => {
    const m = resolveSpec(spec({ include: ['chr2'], order: 'length' }));
    const a = m.assemblies.get(A)!;
    expect(a.systems).toEqual(['/data/hg38/(sequence :name chr2)']);
    expect(a.order.map((i) => a.sequences[i].name)).toEqual([
      'chr1',
      'chr2',
      'chrX',
      'chr10',
      'HLA-A*01:01:01:01',
      'chrUn_KI270302v1',
      'chrM',
    ]);
  });

  it('keeps `include`, orders by `length` or a list', () => {
    const m = resolveSpec(spec({ include: { min_length: 1500 }, order: 'length' }));
    expect(m.assemblies.get(A)!.systems.map((s) => s.split(':name ')[1])).toEqual(['chr1)', 'chr2)', 'chrX)']);
    const listed = resolveSpec(spec({ order: ['X', '10'], include: ['X', '10', '1'] }));
    expect(listed.assemblies.get(A)!.systems.map((s) => s.split(':name ')[1])).toEqual(['chrX)', 'chr10)', 'chr1)']);
  });

  it('is pending until its source opens, then joins its coverage through the names rule (EC-112, EC-186, EC-188)', () => {
    const b = createSeqSpec();
    b.data.source('sizes', { format: 'chrom-sizes', uri: 'hg38.chrom.sizes' });
    b.data.assembly('hg38', { name: 'GRCh38', source: '/data/sizes', sequences: { name: ['chr1'], length: [3000] } });
    const root = b.build();
    const before = resolveSpec(root);
    expect(before.assemblies.get(A)!.pending).toBe(true);
    expect(before.assemblies.get(A)!.systems).toEqual([]);
    const coverage = {
      systems: [
        { name: '1', length: 3000 },
        { name: '2', length: 2000 },
      ],
      shapes: [],
      levels: [1],
    };
    const lengths = sourceLengths(before, (s) => (s === '/data/sizes' ? coverage : undefined));
    const after = resolveSpec(root, { lengths });
    const a = after.assemblies.get(A)!;
    expect(a.pending).toBe(false);
    expect(a.systems).toEqual(['/data/hg38/(sequence :name chr1)', '/data/hg38/(sequence :name "2")']);
    // read by region by the name its source lists it by
    expect(after.systems.get(a.systems[0])!.name).toBe('1');
    expect(after.hash).not.toBe(before.hash);
  });

  it('is unavailable when a sequence has no length, said with its names (EC-188)', () => {
    const b = createSeqSpec();
    b.data.source('fa', { format: 'fasta', uri: 'g.fa' });
    b.data.assembly('g', { name: 'toy', source: '/data/fa' });
    const root = b.build();
    const m0 = resolveSpec(root);
    const coverage = { systems: [{ name: 'c1' }, { name: 'c2', length: 4 }], shapes: [], levels: [1] };
    const m = resolveSpec(root, { lengths: sourceLengths(m0, () => coverage) });
    expect(m.assemblies.get('/data/g')!.systems).toEqual([]);
    expect(codes(m.issues)).toContain('seq.resolve.assembly.lengths');
  });

  it('says a source over it that names another assembly, and inline rows naming no sequence (EC-187)', () => {
    const root = spec({}, (b) => {
      b.data.source('bed', { format: 'bed', uri: 'x.bed', assembly: 'hg19' });
      b.data.ranges('peaks', { target: A, source: '/data/bed' });
      b.data.ranges('inline', {
        target: A,
        rows: { start: [0, 1], end: [1, 2], target: ['1', 'chr99'] } as never,
      });
    });
    const said = codes(resolveSpec(root).issues).filter((c) => c.startsWith('seq.resolve.assembly'));
    // (the fixture's `MT` override names no sequence of it: its own code, first)
    expect(said).toEqual([
      'seq.resolve.assembly.override',
      'seq.resolve.assembly.mismatch',
      'seq.resolve.assembly.name',
    ]);
  });
});

describe('resolve at genome scale', () => {
  itWithin('an assembly of 100,000 sequences kept to 25 by `include`, linear', (within) => {
    const names = Array.from({ length: 100_000 }, (_, i) => (i < 25 ? `chr${i + 1}` : `scaffold_${i}`));
    const b = createSeqSpec();
    b.data.assembly(
      'big',
      assembly(
        'big',
        names.map((name, i) => ({ name, length: 1000 + i })),
        { include: ['^chr[0-9]+$'] },
      ),
    );
    const root = b.build();
    resolveSpec(root); // warm
    // one sample bounded at 50 ms (the design's representative)
    const m = within(50, () => resolveSpec(root));
    expect(m.assemblies.get('/data/big')!.systems.length).toBe(25);
  });
});

describe('rows over an assembly with `include`', () => {
  it('a row on a sequence `include` leaves out is no unknown name, by its name or its alias alike', () => {
    const root = spec({ include: ['chr2'] }, (b) =>
      b.data.ranges('rows', {
        target: A,
        rows: { start: [0, 0, 0], end: [1, 1, 1], target: ['chr1', '1', 'chr99'] } as never,
      }),
    );
    const said = resolveSpec(root).issues.filter((i) => i.code === 'seq.resolve.assembly.name');
    expect(said.map((i) => i.message)).toEqual(['Rows name 1 sequences the assembly lacks (chr99)']);
  });
});

describe('the genome forms', () => {
  const root = spec();
  // resolved inside the tests: the genome chunk loads before them (`beforeAll`)
  const at = (address: string) => resolveAddressIn(root, overlayQueryModel(resolveSpec(root), EmptyOverlay), address);
  const loc = (address: string) => {
    const r = at(address);
    return r.result
      ? plain(locationOfResult(root, overlayQueryModel(resolveSpec(root), EmptyOverlay), r.result).location)
      : r.issues.map((i) => i.code);
  };

  it('names a derived sequence by any of its names, printing the assembly’s (EC-200)', () => {
    const r = at('/data/hg38/(sequence :name "1")');
    expect(r.status).toBe('resolved');
    expect(r.result?.type === 'nodes' && r.result.nodes[0].address).toBeTruthy();
    expect(loc('/data/hg38/(sequence :name NC_000001.11)/(all)')).toEqual({
      parts: [{ system: '/data/hg38/(sequence :name chr1)', ranges: [0, 3000] }],
      rects: [],
    });
    expect(loc('/data/hg38/(sequence :name "HLA-A*01:01:01:01")/(all)')).toEqual({
      parts: [{ system: '/data/hg38/(sequence :name HLA-A*01:01:01:01)', ranges: [0, 30] }],
      rects: [],
    });
  });

  it('reads a locus as readers type it: the conformance table (EC-188, EC-189, EC-200)', () => {
    const table: [string, unknown][] = [
      ['chr1', [0, 3000]],
      ['1', [0, 3000]],
      ['chr1:101', [100, 101]],
      ['chr1:1,001-2,000', [1000, 2000]],
      ['X:1-10', [0, 10]],
      ['HLA-A*01:01:01:01', [0, 30]],
      ['HLA-A*01:01:01:01:2-3', [1, 3]],
      ['chr1:2,000-9,999', [1999, 3000]],
    ];
    for (const [text, ranges] of table) {
      const l = loc(`/data/hg38/(locus :text "${text}")`) as { parts?: { ranges: number[] }[] };
      expect(l.parts?.[0].ranges, text).toEqual(ranges);
    }
    expect(loc('/data/hg38/(locus :text "chr99:1-2")')).toContain('seq.query.locus.name');
    expect(loc('/data/hg38/(locus :text "MT")')).toContain('seq.query.locus.name');
    expect(loc('/data/hg38/(locus :text "chr1:x-y")')).toContain('seq.query.position.invalid');
  });

  it('says an unknown name whole, never a shorter one at a position, and never reads a bare number as a name', () => {
    const said = (address: string) => at(address).issues.map((i) => `${i.code}: ${i.message}`);
    expect(said('/data/hg38/(locus :text "HLA-A*01:01:01:02")')).toEqual([
      'seq.query.locus.name: No sequence of the assembly is named HLA-A*01:01:01:02',
    ]);
    // `01` reads as the number 1 (chr1's alias): refused, its text lost
    expect(codes(at('/data/hg38/(sequence :name 01)').issues)).toEqual(['seq.query.form.param']);
    expect(at('/data/hg38/(sequence :name "01")').issues.map((i) => i.code)).toEqual(['seq.query.locus.name']);
  });

  it('takes `(block :name …)` under a section of assemblies alone', () => {
    const r = spec({}, (b) => {
      const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
      b.views.view('v', (v) => v.section('plain', { of: s }, (sec) => sec.block('b', { of: s })));
    });
    const l = resolveAddressIn(r, overlayQueryModel(resolveSpec(r), EmptyOverlay), '/views/v/plain/(block :name s)');
    expect(codes(l.issues)).toEqual(['seq.query.form.prefix']);
  });

  it('names a derived block of a section of assemblies by its sequence’s name (EC-200)', () => {
    const r = spec({}, (b) => b.views.view('v', (v) => genome(v, 'genome', A)));
    const m = resolveSpec(r);
    const l = resolveAddressIn(r, overlayQueryModel(m, EmptyOverlay), '/views/v/genome/(block :name M)');
    expect(l.result && plain(locationOfResult(r, overlayQueryModel(m, EmptyOverlay), l.result).location).parts).toEqual(
      [{ system: '/data/hg38/(sequence :name chrM)', ranges: [0, 16] }],
    );
  });
});

describe('unread params', () => {
  it('are said, one info a node, naming the lane that reads each', () => {
    const b = createSeqSpec();
    const g = b.data.assembly('g', { ...assembly('g', hg38), names: 'ucsc' } as never);
    b.data.ranges('over', { target: g, rows: { start: [0], end: [1], target: ['chr1'] } } as never);
    b.views.view('v', (v) => {
      v.section('a', { of: g, member_order: 'x' } as never, (s) =>
        s.block('b', { of: '/data/g/(sequence :name chr1)', direction: 'reversed', scale: 2 } as never),
      );
      v.section('c', { of: g });
      v.bridge('br', { a: '/views/v/a', b: '/views/v/c' } as never);
    });
    const said = resolveSpec(b.build())
      .issues.filter((i) => i.code === 'seq.resolve.param.unread')
      .map((i) => `${i.address} ${i.message}`);
    expect(said).toEqual([
      '/data/g `assembly.names`: read by lane GF',
      '/views/v/a `section.member_order`: read by lane LA',
    ]);
    expect(resolveSpec(spec()).issues.filter((i) => i.code === 'seq.resolve.param.unread')).toEqual([]);
  });
});

describe('sections of assemblies', () => {
  it('derive a block a sequence, two assemblies in turn, keyed by name (EC-200, EC-210)', () => {
    const b = createSeqSpec();
    b.data.assembly(
      'peach',
      assembly('peach', [
        { name: 'Pp01', length: 10 },
        { name: 'Pp02', length: 5 },
      ]),
    );
    b.data.assembly('grape', assembly('grape', [{ name: 'chr1', length: 8 }]));
    b.views.view('v', (v) => genome(v, 'two', ['/data/peach', '/data/grape']));
    const view = resolveView(resolveSpec(b.build()), '/views/v')!;
    const s = view.sections[0];
    expect(s.blocks?.map((x) => x.address)).toEqual([
      '/views/v/two/(block :name Pp01)',
      '/views/v/two/(block :name Pp02)',
      '/views/v/two/(block :name chr1)',
    ]);
    expect(s.axis.length).toBe(23);
  });

  it('derive their blocks once a build', () => {
    const root = spec({}, (b) => b.views.view('v', (v) => genome(v, 'g', A)));
    const model = resolveSpec(root);
    const include = vi.spyOn(model.assemblies.get(A)!, 'include');
    expect(resolveView(model, '/views/v')!.sections[0].blocks?.length).toBe(7);
    expect(include).toHaveBeenCalledTimes(1);
  });

  it('show their `include` and say past `maxAssemblyBlocks` (EC-210)', () => {
    const root = spec({}, (b) => b.views.view('v', (v) => genome(v, 'g', A, { include: ['1', '2', 'X'] })));
    expect(resolveView(resolveSpec(root), '/views/v')!.sections[0].blocks?.length).toBe(3);
    const capped = resolveView(resolveSpec(root, { config: { maxAssemblyBlocks: 2 } }), '/views/v')!;
    expect(capped.sections[0].blocks?.length).toBe(2);
    expect(codes(capped.issues)).toContain('seq.view.section.blocks');
  });
});

describe('a block’s direction', () => {
  /** `block.direction` flipped: `reversed` against the section, else back to the section’s. */
  const flip = (n: Node): Node => ({
    ...n,
    params: { ...n.params, direction: n.params?.direction === 'reversed' ? 'section' : 'reversed' },
  });

  // What the model reads of an override is the derived block's params (`blockNodes`); the view's block axes read them
  // (`view/_test/genome-blocks.test.ts`), so the property is on the params.
  it('an override flipped twice derives the same blocks, flipped once only its own block changes (EC-206, property)', () => {
    const rng = new PCG(3);
    const model = resolveSpec(spec());
    const section = (children: Node[]): Node => ({ kind: 'section', id: 'g', params: { of: A }, children });
    // a block's direction read as the model reads it: none is the section's
    const derived = (children: Node[]) =>
      blockNodes(model, section(children), '/views/v/g').map(
        (n): Record<string, Value> => ({ direction: 'section', ...n.params }),
      );
    const shown = derived([]).map((p) => p.of as string);
    for (let round = 0; round < 20; round++) {
      const overrides: Node[] = [];
      for (let i = 0; i < 1 + (rng.int() % 4); i++)
        overrides.push({
          kind: 'block',
          id: `o${i}`,
          params: {
            of: shown[rng.int() % shown.length],
            ...(rng.int() % 2 ? { direction: 'reversed' } : {}),
          },
        });
      const k = rng.int() % overrides.length;
      const base = derived(overrides);
      expect(derived(overrides.map((c, i) => (i === k ? flip(flip(c)) : c)))).toEqual(base);
      const once = derived(overrides.map((c, i) => (i === k ? flip(c) : c)));
      const target = overrides[k].params!.of;
      // the last override of a sequence wins: flipping an earlier one of the same sequence changes nothing
      const wins = overrides.map((c) => c.params!.of).lastIndexOf(target) === k;
      const changed = once.flatMap((params, i) => (JSON.stringify(params) === JSON.stringify(base[i]) ? [] : [i]));
      expect(changed.length).toBe(wins ? 1 : 0);
      for (const i of changed) expect({ ...once[i], direction: 0 }).toEqual({ ...base[i], direction: 0 });
    }
  });
});

describe('a map between assemblies', () => {
  function synteny(): Node {
    const b = createSeqSpec();
    b.data.assembly(
      'a',
      assembly('A', [
        { name: 'chr1', length: 100 },
        { name: 'chr2', length: 100 },
      ]),
    );
    b.data.assembly(
      'b',
      assembly('B', [
        { name: 'chr1', length: 100 },
        { name: 'chr2', length: 100 },
      ]),
    );
    b.data.map('ab', {
      from: '/data/a',
      to: '/data/b',
      method: 'segments',
      segments: {
        from: [0, 50, 0],
        to: [10, 99, 0],
        length: [20, 30, 40],
        strand: ['+', '-', '+'],
        from_name: ['chr1', 'chr1', 'chr2'],
        to_name: ['chr2', 'chr2', 'chr1'],
      },
    });
    return b.build();
  }

  it('is one map over every sequence pair: a path between two sequences is one step of it (EC-207)', () => {
    const m = resolveSpec(synteny());
    const a1 = '/data/a/(sequence :name chr1)';
    const b2 = '/data/b/(sequence :name chr2)';
    const r = m.path(a1, b2);
    expect(r.issue).toBeUndefined();
    expect(r.path?.edges.map((e) => e.map)).toEqual(['/data/ab']);
    expect(m.path('/data/a/(sequence :name chr2)', '/data/b/(sequence :name chr1)').path).toBeDefined();
    expect(m.path(a1, '/data/b/(sequence :name chr1)').path).toBeUndefined();
  });

  it('carries a location across it, each run in its own direction (EC-206, EC-212)', () => {
    const m = resolveSpec(synteny());
    const p = m.project(
      { parts: [{ system: '/data/a/(sequence :name chr1)', ranges: new Int32Array([0, 5, 50, 52]) }], rects: [] },
      '/data/b/(sequence :name chr2)',
    );
    expect(p.location.parts.map((x) => [x.system, Array.from(x.ranges)])).toEqual([
      ['/data/b/(sequence :name chr2)', [10, 15, 98, 100]],
    ]);
  });
});

describe('the builder’s genome helpers', () => {
  it('writes a genome section with a ruler in bases, and an ideogram navigator given bands', () => {
    const root = spec({}, (b) => {
      const bands = b.data.ranges('bands', { target: A, rows: { start: [0], end: [10] } });
      b.views.view('v', (v) => genome(v, 'g', A, { bands }));
    });
    const view = root.children!.find((c) => c.kind === 'views')!.children![0];
    expect(view.children!.map((c) => c.kind)).toEqual(['section', 'navigators']);
    expect(view.children![0].children![0].children![0].params).toMatchObject({ type: 'ruler', units: 'bases' });
  });

  it('writes one sequence’s blocks in another’s order, short runs merged, reversed where a run is', () => {
    const runs = {
      from: [500_000, 0, 200_000, 260_000],
      to: [0, 1_000_000, 700_000, 400_000],
      length: [300_000, 200_000, 60_000, 150_000],
      strand: ['+', '+', '-', '+'],
    };
    expect(colinearBlocks(runs, '/data/b')).toEqual([
      { of: '/data/b', range: [640_001, 1_200_000] },
      { of: '/data/b', range: [400_000, 550_000] },
      { of: '/data/b', range: [0, 300_000] },
    ]);
    expect(colinearBlocks(runs, '/data/b', { minLength: 10 })[1]).toEqual({
      of: '/data/b',
      range: [640_001, 700_001],
      direction: 'reversed',
    });
  });
});
