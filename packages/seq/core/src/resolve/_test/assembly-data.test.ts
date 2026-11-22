/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import {
  type Coverage,
  DataStore,
  openSeq,
  type RangesList,
  readSeq,
  type SeqData,
  seqToValue,
} from '@molstar/seqx-core/data';
import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { beforeAll, describe, expect, it } from 'vitest';
import { loadGenome } from '../genome.js';
import { resolveSpec, sourceLengths } from '../model.js';

// Data over an assembly read by sequence (_spec/README.md "Resolve", Assemblies: data over an assembly): a property
// over generated assemblies (names with and without `chr`, aliases, `include`, `order`): each derived sequence's plan
// reads exactly the rows naming it by the names rule, so a region read over the assembly is the per-sequence reads
// joined; a source's plans name each sequence as the source lists it (EC-187, EC-189).

beforeAll(async () => {
  await loadGenome();
});

const signal = new AbortController().signal;

interface Generated {
  readonly sequences: { name: string; length: number; aliases: string[] }[];
  readonly include?: string[];
  readonly order?: 'length' | 'source';
  /** Rows: a start, an end and the name a row's `target` gives (a name, an alias, the other prefix, or unknown). */
  readonly rows: { start: number; end: number; target: string; seq: number }[];
}

function generate(seed: number): Generated {
  const rng = new PCG(seed);
  const int = (n: number) => rng.int() % n;
  const n = 2 + int(8);
  const sequences = Array.from({ length: n }, (_, i) => ({
    name: int(2) ? `chr${i + 1}` : `${i + 1}`,
    length: 1000 + int(9000),
    aliases: int(3) ? [] : [`NC_${1000 + i}`],
  }));
  const include = int(2) ? sequences.filter(() => int(3) > 0).map((s) => s.name) : undefined;
  const order = (['length', 'source', undefined] as const)[int(3)];
  const rows = Array.from({ length: 5 + int(60) }, () => {
    const seq = int(n + 1);
    if (seq === n) return { start: 0, end: 1, target: 'chrUnknown', seq: -1 };
    const s = sequences[seq];
    const bare = s.name.replace(/^chr/, '');
    const names = [s.name, bare === s.name ? `chr${bare}` : bare, ...s.aliases];
    const start = int(s.length - 10);
    return { start, end: start + 1 + int(9), target: names[int(names.length)], seq };
  });
  return { sequences, ...(include ? { include } : {}), ...(order ? { order } : {}), rows };
}

function specOf(g: Generated, source?: SeqData) {
  const b = createSeqSpec();
  const a = b.data.assembly(
    'g',
    assembly('G', g.sequences, {
      ...(g.include ? { include: g.include } : {}),
      ...(g.order ? { order: g.order } : {}),
    } as never),
  );
  const src = source && b.data.source('src', { format: 'seq', value: seqToValue(source) as never });
  const rows = {
    start: g.rows.map((r) => r.start),
    end: g.rows.map((r) => r.end),
    target: g.rows.map((r) => r.target),
  };
  b.data.ranges('bed', { target: a, ...(src ? { source: src } : { rows }) } as never);
  return b.build();
}

/** A plan's rows read whole (an inline plan's own source), as `start-end` sorted. */
async function readPlan(plan: { kind: string; source: unknown; system: string }, length: number): Promise<string[]> {
  if (plan.kind !== 'inline') throw new Error('an inline plan');
  const value = (plan.source as { value: unknown }).value;
  const opened = openSeq(readSeq(value as never).data!);
  const a = await opened.read({ system: plan.system, start: 0, end: length, level: 1, shapes: ['ranges'] }, signal);
  if (!a.ok) throw new Error(a.issue.message);
  const t = a.tables.ranges as RangesList;
  return Array.from(t.start, (s, i) => `${s}-${t.end[i]}`).sort();
}

describe('data over an assembly, read by sequence', () => {
  it('inline rows: each sequence’s plan holds exactly the rows naming it (names, aliases, the prefix), over generated assemblies with `include` and `order` (EC-189)', async () => {
    for (let seed = 1; seed <= 40; seed++) {
      const g = generate(seed);
      const m = resolveSpec(specOf(g));
      const a = m.assemblies.get('/data/g')!;
      const e = m.annotations.get('/data/bed')!;
      expect(e.target).toBe('/data/g');
      const joined: string[] = [];
      for (let i = 0; i < g.sequences.length; i++) {
        const system = a.systemAt[a.find(g.sequences[i].name)!];
        const want = g.rows
          .filter((r) => r.seq === i)
          .map((r) => `${r.start}-${r.end}`)
          .sort();
        const plan = system === undefined ? undefined : e.bySystem!.get(system);
        // a sequence `include` leaves out reads nothing; one no row names has no plan
        if (system === undefined || !want.length) {
          expect(plan, `seed ${seed} ${g.sequences[i].name}`).toBeUndefined();
          continue;
        }
        const got = await readPlan(plan!, g.sequences[i].length);
        expect(got, `seed ${seed} ${g.sequences[i].name}`).toEqual(want);
        joined.push(...got);
      }
      // the reads joined are every row whose sequence is kept, once each
      const kept = g.rows.filter((r) => r.seq >= 0 && a.systemAt[a.find(g.sequences[r.seq].name)!] !== undefined);
      expect(joined.length).toBe(kept.length);
      // a row naming no sequence: counted, said once
      const unknown = m.issues.filter((i) => i.code === 'seq.resolve.assembly.name');
      expect(unknown.length).toBe(g.rows.some((r) => r.seq < 0) ? 1 : 0);
    }
  });

  it('a source: each sequence read under the name the source lists it by, after it opens (EC-187, EC-189)', () => {
    for (let seed = 1; seed <= 20; seed++) {
      const g = generate(seed);
      // the source lists each sequence by its other spelling (`1` for `chr1`, and back), and one name of its own
      const listed = g.sequences.map((s) => (s.name.startsWith('chr') ? s.name.slice(3) : `chr${s.name}`));
      const coverage: Coverage = {
        systems: [
          ...listed.map((name, i) => ({ name, length: g.sequences[i].length })),
          { name: 'scaffold', length: 5 },
        ],
        shapes: [...listed, 'scaffold'].map((system) => ({ system, shape: 'ranges' as const })),
        levels: [1],
      };
      const source: SeqData = { coverage, pairs: [] };
      const root = specOf(g, source);
      const before = resolveSpec(root);
      const e0 = before.annotations.get('/data/bed')!;
      // before the source opens: the printed names
      const a0 = before.assemblies.get('/data/g')!;
      for (let i = 0; i < a0.sequences.length; i++) {
        const system = a0.systemAt[i];
        if (system !== undefined) expect(e0.bySystem!.get(system)?.system).toBe(a0.sequences[i].name);
      }
      const lengths = sourceLengths(before, (s) => (s === '/data/src' ? coverage : undefined));
      const m = resolveSpec(root, { lengths });
      const a = m.assemblies.get('/data/g')!;
      const e = m.annotations.get('/data/bed')!;
      for (let i = 0; i < g.sequences.length; i++) {
        const system = a.systemAt[a.find(g.sequences[i].name)!];
        if (system === undefined) expect(e.bySystem!.size).toBeLessThan(g.sequences.length);
        else expect(e.bySystem!.get(system)?.system).toBe(listed[i]);
      }
      expect(m.issues.filter((i) => i.code === 'seq.resolve.assembly.name').map((i) => i.message)).toEqual([
        'The source names 1 sequences the assembly lacks (scaffold)',
      ]);
    }
  });

  it('over two assemblies, a listed name neither has is said, once; inline values, which name no sequence, say they are not read', () => {
    const b = createSeqSpec();
    const a = b.data.assembly('a', assembly('A', [{ name: 'chrA1', length: 100 }]));
    const g = b.data.assembly('b', assembly('B', [{ name: 'chrB1', length: 100 }]));
    const coverage: Coverage = {
      systems: ['chrA1', 'chrB1', 'chrZ'].map((name) => ({ name, length: 100 })),
      shapes: ['chrA1', 'chrB1', 'chrZ'].map((system) => ({ system, shape: 'pairs' as const })),
      levels: [1],
    };
    const src = b.data.source('src', { format: 'seq', value: seqToValue({ coverage, pairs: [] }) as never });
    b.data.pairs('syn', { a, b: g, source: src });
    b.data.values('vals', { target: a, values: new Float32Array(100) });
    const root = b.build();
    const lengths = sourceLengths(resolveSpec(root), (s) => (s === '/data/src' ? coverage : undefined));
    const said = resolveSpec(root, { lengths })
      .issues.filter((i) => i.code.startsWith('seq.resolve.assembly.'))
      .map((i) => `${i.address} ${i.code}: ${i.message}`);
    expect(said).toEqual([
      '/data/syn seq.resolve.assembly.name: The source names 1 sequences neither assembly has (chrZ)',
      '/data/vals seq.resolve.assembly.inline: inline values over an assembly is not read yet: name a source, or target a sequence',
    ]);
  });

  it('a region read over the assembly through the store is the per-sequence reads joined, in the assembly’s order (with `M` and `MT`, EC-189)', async () => {
    const cases: Generated[] = [
      ...Array.from({ length: 12 }, (_, k) => generate(100 + k)),
      {
        // `M` is `chrM` by the prefix; `MT` is its own sequence, never `M`'s
        sequences: [
          { name: 'chrM', length: 16_571, aliases: [] },
          { name: 'MT', length: 16_569, aliases: [] },
          { name: 'chr1', length: 3000, aliases: [] },
        ],
        rows: [
          { start: 10, end: 20, target: 'M', seq: 0 },
          { start: 30, end: 40, target: 'chrM', seq: 0 },
          { start: 50, end: 60, target: 'MT', seq: 1 },
          { start: 70, end: 80, target: '1', seq: 2 },
        ],
      },
    ];
    for (const g of cases) {
      const m = resolveSpec(specOf(g));
      const a = m.assemblies.get('/data/g')!;
      const e = m.annotations.get('/data/bed')!;
      const store = new DataStore();
      const plan = e.access;
      if (!plan || plan.kind !== 'inline') {
        expect(g.rows.every((r) => r.seq < 0 || a.systemAt[a.find(g.sequences[r.seq].name)!] === undefined)).toBe(true);
        continue;
      }
      store.attach(plan.source.address, openSeq(readSeq(plan.source.value as never).data!));
      // the region over the assembly: every derived sequence whole, in the assembly's order
      const wants = a.systems.flatMap((system) => {
        const p = e.bySystem!.get(system);
        const length = m.systems.get(system)!.length!;
        return p
          ? [
              {
                source: p.kind === 'inline' ? p.source.address : p.source,
                system: p.system,
                start: 0,
                end: length,
                level: 1,
              },
            ]
          : [];
      });
      store.want(wants);
      await store.whenSettled();
      const joined: string[] = [];
      for (const w of wants) {
        // a row once (a tile holds every row meeting it: one across a tile's edge is in two), by its row
        const rows = new Map<number, [number, number]>();
        for (const k of store.keysOf(w)) {
          const entry = store.get(k);
          if (entry?.status !== 'loaded') continue;
          const t = entry.tables.ranges as RangesList;
          for (let i = 0; i < t.start.length; i++) rows.set(t.row![i], [t.start[i], t.end[i]]);
        }
        joined.push(
          ...[...rows.values()].sort((x, y) => x[0] - y[0] || x[1] - y[1]).map(([st, en]) => `${w.system}:${st}-${en}`),
        );
      }
      // the per-sequence rows, joined in the assembly's order
      const expected = a.systems.flatMap((system) => {
        const i = a.systemAt.indexOf(system);
        return g.rows
          .filter((r) => r.seq >= 0 && a.find(g.sequences[r.seq].name) === i)
          .sort((x, y) => x.start - y.start || x.end - y.end)
          .map((r) => `${system}:${r.start}-${r.end}`);
      });
      expect(joined).toEqual(expected);
    }
  });
});

describe('inline pairs, grids and matrices over an assembly', () => {
  const toy = (b: ReturnType<typeof createSeqSpec>) =>
    b.data.assembly(
      'g',
      assembly('G', [
        { name: 'chr1', length: 4 },
        { name: 'chr2', length: 2 },
      ]),
    );
  const plans = (root: ReturnType<ReturnType<typeof createSeqSpec>['build']>, address: string) => {
    const m = resolveSpec(root);
    const a = m.annotations.get(address)!;
    return { m, a, systems: [...(a.bySystem?.keys() ?? [])].map((s) => s.replace(/^.*:name |\)$/g, '')) };
  };

  it('a side on a sequence is that system: pairs between a sequence and an assembly hold both’s units', () => {
    const b = createSeqSpec();
    const g = toy(b);
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
    b.data.pairs('p', {
      a: s,
      b: g,
      rows: { a_start: [0, 2], a_end: [1, 3], b_start: [0, 1], b_end: [1, 2], b_name: ['2', 'chr1'] },
    });
    const { a } = plans(b.build(), '/data/p');
    expect(a.target).toBe('/data/s');
    // held, each end on its unit: the sequence first, then the assembly's sequences in order
    expect(a.units?.map((u) => u.replace(/^.*:name |\)$/g, ''))).toEqual(['/data/s', 'chr1', 'chr2']);
    expect([...a.pairs!.aUnit!]).toEqual([0, 0]);
    expect([...a.pairs!.bUnit!]).toEqual([2, 1]);
  });

  it('a grid’s sequence pair of the wrong count is left out and said; a matrix’s states its inline names', () => {
    const b = createSeqSpec();
    const g = toy(b);
    b.data.grid('h', {
      x: g,
      values: [1, 2, 3, 4, 5, 6, 7, 8, 9],
      x_name: ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr2'],
      y_name: ['chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2'],
    });
    b.data.matrix('m', { target: g, values: [1, 2, 3, 4], states: ['x', 'y'], name: ['chr2', 'chr2'] });
    const root = b.build();
    const { m } = plans(root, '/data/h');
    expect(m.issues.filter((i) => i.code === 'seq.resolve.grid.size').map((i) => i.message)).toEqual([
      "A sequence pair's grid has 1 of 4 values: left out",
    ]);
    expect(m.systems.get('/data/m/(states)')?.length).toBe(2);
    expect(plans(root, '/data/m').systems).toEqual(['chr2']);
  });

  const said = (m: ReturnType<typeof resolveSpec>) =>
    m.issues.filter((i) => i.code.startsWith('seq.resolve.assembly.')).map((i) => `${i.code}: ${i.message}`);
  const pair = (n: number) => ({
    a_start: Array.from({ length: n }, () => 0),
    a_end: Array.from({ length: n }, () => 1),
    b_start: Array.from({ length: n }, () => 0),
    b_end: Array.from({ length: n }, () => 1),
  });

  it('ends past a sequence’s length: left out, counted and said once (EC-189, EC-208)', () => {
    const b = createSeqSpec();
    const g = toy(b);
    b.data.pairs('p', {
      a: g,
      rows: { ...pair(3), b_end: [1, 3, 5], a_name: ['chr1', 'chr1', 'chr1'], b_name: ['chr2', 'chr2', 'chr2'] },
    });
    const { m, a } = plans(b.build(), '/data/p');
    expect(a.pairs!.n).toBe(1);
    expect(said(m)).toEqual(['seq.resolve.assembly.bounds: 2 rows are past their sequence’s length: left out']);
  });

  it('rows naming no sequence (null, a short column, no string) and items past a length: counted and said once (EC-189)', () => {
    const b = createSeqSpec();
    const g = toy(b);
    b.data.pairs('p', {
      a: g,
      rows: { ...pair(3), a_name: ['chr1', null, 'chr2'], b_name: ['chr1', 'chr2', null] } as never,
    });
    b.data.matrix('m', {
      target: g,
      values: [1, 2, 3, 4, 5, 6],
      states: ['x'],
      name: ['chr2', 'chr2', 'chr2', 'chr1', 7, null] as never,
    });
    const root = b.build();
    const m = resolveSpec(root);
    expect(m.annotations.get('/data/p')!.pairs!.n).toBe(1);
    expect(said(m)).toEqual([
      'seq.resolve.assembly.name: 2 rows name no sequence: left out',
      'seq.resolve.assembly.name: 2 items name no sequence: left out',
      'seq.resolve.assembly.bounds: 1 items are past their sequence’s length: left out',
    ]);
  });

  it('no name is unknown while the assembly is pending (its source not open): said once it lists its sequences', () => {
    const b = createSeqSpec();
    const sizes = b.data.source('sizes', { format: 'chrom-sizes', value: 'chr1\t100\n' as never });
    const g = b.data.assembly('g', { name: 'G', source: sizes });
    b.data.pairs('p', { a: g, rows: { ...pair(1), a_name: ['chrZ'], b_name: ['chr1'] } });
    const root = b.build();
    expect(said(resolveSpec(root))).toEqual([]);
  });

  it('two assemblies’ sequences in order, each assembly’s whole before the next (a grid’s systems)', () => {
    const b = createSeqSpec();
    const three = (n: string) => [1, 2, 3].map((i) => ({ name: `${n}${i}`, length: 1 }));
    const a = b.data.assembly('a', assembly('A', three('a'), { order: 'source' }));
    const c = b.data.assembly('c', assembly('C', three('c'), { order: 'source' }));
    const xs = ['a3', 'a1', 'a2', 'a3', 'a1', 'a2'];
    const ys = ['c2', 'c3', 'c1', 'c1', 'c2', 'c3'];
    b.data.grid('h', { x: a, y: c, values: xs.map(() => 1), x_name: xs, y_name: ys });
    const { systems } = plans(b.build(), '/data/h');
    expect(systems).toEqual(['a1', 'a2', 'a3', 'c1', 'c2', 'c3']);
  });

  it('a grid over 3,400 one-base scaffolds serves the pairs with values alone; a pair past the bound refused unread', () => {
    const b = createSeqSpec();
    const n = 3400;
    const seqs = Array.from({ length: n }, (_, i) => ({ name: `scaffold_${i}`, length: 1 }));
    const g = b.data.assembly('g', assembly('S', [...seqs, { name: 'chrBig', length: 5000 }], { order: 'source' }));
    const names = [...seqs.map((q) => q.name), 'chrBig'];
    b.data.grid('h', { x: g, values: names.map(() => 1), x_name: names, y_name: names });
    const t = performance.now();
    const m = resolveSpec(b.build());
    const a = m.annotations.get('/data/h')!;
    const value = (a.access as { source: { value: unknown } }).source.value;
    const data = readSeq(value as never).data!;
    expect(data.pairs.length).toBe(n);
    expect(performance.now() - t).toBeLessThan(2000);
    expect(m.issues.filter((i) => i.code === 'seq.resolve.grid.size').map((i) => i.message)).toEqual([
      "A sequence pair's grid of 25000000 values is past 16777216",
    ]);
  });
});
