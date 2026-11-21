/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { itWithin } from '@molstar/test-harness/timing';
import { beforeAll, describe, expect, it } from 'vitest';
import { assembly, colinearBlocks, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { locationOfResult, resolveAddressIn } from '@molstar/seqx-core/query';
import { resolveView } from '@molstar/seqx-core/view';
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
    // `MT` reaches chrM through the prefix rule? never: `M` only. The override names `MT`, which this assembly lacks
    const chrM = m.systems.get('/data/hg38/(sequence :name chrM)')!;
    expect([chrM.length, chrM.cyclic, chrM.alphabet]).toEqual([16, false, 'dna']);
    expect(codes(m.issues)).toEqual([]);
  });

  it('applies an override found through the names rule', () => {
    const b = createSeqSpec();
    b.data.assembly('hg38', assembly('GRCh38', hg38), (a) =>
      a.sequence('mito', { name: 'M', cyclic: true, label: 'mitochondrion' }),
    );
    const s = resolveSpec(b.build()).systems.get('/data/hg38/(sequence :name chrM)')!;
    expect([s.cyclic, s.label]).toEqual([true, 'mitochondrion']);
  });

  it('keeps `include`, orders by `length` or a list', () => {
    const m = resolveSpec(spec({ include: { min_length: 1500 }, order: 'length' }));
    expect(m.assemblies.get(A)!.systems.map((s) => s.split(':name ')[1])).toEqual(['chr1)', 'chr2)', 'chrX)']);
    const listed = resolveSpec(spec({ order: ['X', '10'], include: ['X', '10', '1'] }));
    expect(listed.assemblies.get(A)!.systems.map((s) => s.split(':name ')[1])).toEqual(['chrX)', 'chr10)', 'chr1)']);
  });

  it('is pending until its source opens, then joins its coverage through the names rule (EC-112, EC-188)', () => {
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
    expect(said).toEqual(['seq.resolve.assembly.mismatch', 'seq.resolve.assembly.name']);
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
    const r = at('/data/hg38/(sequence :name 1)');
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

  it('names a derived block of a section of assemblies by its sequence’s name (EC-200)', () => {
    const r = spec({}, (b) => b.views.view('v', (v) => genome(v, 'genome', A)));
    const m = resolveSpec(r);
    const l = resolveAddressIn(r, overlayQueryModel(m, EmptyOverlay), '/views/v/genome/(block :name M)');
    expect(l.result && plain(locationOfResult(r, overlayQueryModel(m, EmptyOverlay), l.result).location).parts).toEqual(
      [{ system: '/data/hg38/(sequence :name chrM)', ranges: [0, 16] }],
    );
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

  it('reversed twice is the identity at the model level, over random blocks (property)', () => {
    const rng = new PCG(3);
    for (let round = 0; round < 20; round++) {
      const b = createSeqSpec();
      const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
      b.views.view('v', (v) =>
        v.section('sec', {}, (sec) => {
          for (let i = 0; i < 1 + (rng.int() % 4); i++) {
            const start = rng.int() % 20;
            sec.block(`b${i}`, {
              of: s,
              range: [start, start + 1 + (rng.int() % 20)],
              ...(rng.int() % 2 ? { direction: 'reversed' as const } : {}),
            });
          }
        }),
      );
      const root = b.build();
      const k = rng.int() % 4;
      const [data, views] = root.children!;
      const view = views.children![0];
      const section = view.children![0];
      const flipped = { ...section, children: section.children!.map((c, i) => (i === k ? flip(flip(c)) : c)) };
      const twice: Node = { ...root, children: [data, { ...views, children: [{ ...view, children: [flipped] }] }] };
      const blocks = (r: Node) => JSON.stringify(resolveView(resolveSpec(r), '/views/v')!.sections[0].blocks);
      expect(blocks(twice)).toBe(blocks(root));
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

  it('carries a location across it, each run in its own direction (EC-212)', () => {
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
