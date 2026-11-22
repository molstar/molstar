/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { binPairs, type PairsList, type SeqData, type SeqPair, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import type { RenderConfigOverrides } from '../config.js';
import type { BlockPlacement } from '../layout/blocks.js';
import { recordDrawLists } from '../record.js';
import { hitView } from '../tracks-hit.js';
import { settle } from './helpers.js';

// Inline pairs, grids and matrices over an assembly (_spec/README.md "Resolve", Assemblies: data over an assembly):
// each row names its sequences (`a_name`/`b_name`, `x_name`/`y_name`, `name`) through the names rule, and draws as the
// same rows given per sequence by a source.

beforeAll(async () => {
  await loadGenome();
});

const Size = { width: 900, height: 240 };

interface Seq {
  readonly name: string;
  readonly length: number;
  readonly aliases?: readonly string[];
}

/** Pairs as rows naming their sequences (by any name or alias the assembly knows them by). */
type Rows = {
  readonly a_start: number[];
  readonly a_end: number[];
  readonly b_start: number[];
  readonly b_end: number[];
  readonly a_name: string[];
  readonly b_name: string[];
  readonly type?: string[];
};

/**
 * A section of `seqs` with a links track, its data inline (`rows`) or a `seq` source of the same rows given per
 * sequence (each listed by its own name, every pair in the read of each sequence an end is on).
 */
export function linksSpec(
  seqs: readonly Seq[],
  rows: Rows,
  inline: boolean,
  options: { include?: string[]; order?: 'source' | 'length' | 'natural'; layout?: 'rings' } = {},
): Node {
  const b = createSeqSpec({ title: 'links over an assembly' });
  const g = b.data.assembly('toy', assembly('Toy', seqs, { order: options.order ?? 'source' }));
  let links;
  if (inline) links = b.data.pairs('links', { a: g, rows });
  else {
    const src = b.data.source('src', { format: 'seq', value: seqToValue(perSequence(seqs, rows)) as never });
    links = b.data.pairs('links', { a: g, source: src });
  }
  b.views.view('v', (v) =>
    genome(
      v,
      'genome',
      g,
      { ...(options.include ? { include: options.include } : {}), ...(options.layout ? { layout: 'rings' } : {}) },
      (s) => s.track('links').layer('links', { type: 'pairs', data: links }),
    ),
  );
  return b.build();
}

/** A sequence's own name for any of its names (or none). */
function canonical(seqs: readonly Seq[], name: string): string | undefined {
  const bare = (x: string) => x.replace(/^chr/i, '');
  return seqs.find((s) => s.name === name || s.aliases?.includes(name) || bare(s.name) === bare(name))?.name;
}

/** The rows as a `seq` source: a pairs list a sequence, the rows with an end on it, named by its own name. */
function perSequence(seqs: readonly Seq[], rows: Rows): SeqData {
  const pairs: SeqPair[] = [];
  for (const s of seqs) {
    const at: number[] = [];
    for (let i = 0; i < rows.a_start.length; i++) {
      const a = canonical(seqs, rows.a_name[i]);
      const b = canonical(seqs, rows.b_name[i]);
      if (a && b && (a === s.name || b === s.name)) at.push(i);
    }
    pairs.push({
      request: { system: s.name, start: 0, end: s.length, level: 1, shapes: ['pairs'] },
      answer: {
        ok: true,
        tables: {
          pairs: {
            form: 'list',
            aSystem: categoryColumn(at.map((i) => canonical(seqs, rows.a_name[i]) ?? null)),
            bSystem: categoryColumn(at.map((i) => canonical(seqs, rows.b_name[i]) ?? null)),
            aStart: Int32Array.from(at, (i) => rows.a_start[i]),
            aEnd: Int32Array.from(at, (i) => rows.a_end[i]),
            bStart: Int32Array.from(at, (i) => rows.b_start[i]),
            bEnd: Int32Array.from(at, (i) => rows.b_end[i]),
            ...(rows.type ? { type: categoryColumn(at.map((i) => rows.type![i])) } : {}),
            row: Int32Array.from(at),
          },
        },
      },
    });
  }
  return {
    coverage: {
      systems: seqs.map((s) => ({ name: s.name, length: s.length })),
      shapes: seqs.map((s) => ({ system: s.name, shape: 'pairs' as const, ...(rows.type ? { binBy: 'type' } : {}) })),
      levels: [1],
    },
    pairs,
  };
}

/**
 * The data pass's marks of a frame, as recorded: the geometry, never the keys, empty clips dropped; the pairs' strokes
 * uncoloured (a source's
 * pairs are coloured by kind apart from held and inline ones: a drift the report names, not this lane's).
 */
async function marks(spec: Node, config: RenderConfigOverrides = {}): Promise<string> {
  const s = await settle(spec, [SeqRuntime], config);
  const f = await s.frame({ size: Size, dpr: 1 });
  let out = recordDrawLists(f.lists, ['data']).replace(/(path [^\n]* stroke) #[0-9a-f]{8}/g, '$1');
  // the arcs in one order (held and read pairs are gathered in their own)
  for (const w of ['path', 'rect'])
    out = out.replace(
      new RegExp(`(?:${w} [^\\n]*\\n)+`, 'g'),
      (run) => `${run.split('\n').filter(Boolean).sort().join('\n')}\n`,
    );
  // a clip with nothing in it draws nothing: a source's block it lists with no rows (an inline node plans none)
  for (let was = ''; was !== out; ) {
    was = out;
    out = out.replace(/clip [^\n]*\nunclip\n/g, '');
  }
  return out;
}

function seeded(seed: number) {
  let x = seed >>> 0 || 1;
  return () => {
    x ^= x << 13;
    x ^= x >>> 17;
    x ^= x << 5;
    return (x >>> 0) / 4294967296;
  };
}

export const Toy: Seq[] = [
  { name: 'chr1', length: 3000 },
  { name: 'chr2', length: 2000 },
  { name: 'chr3', length: 1000 },
];

describe('inline pairs over an assembly', () => {
  const rows: Rows = {
    a_start: [100, 500, 200],
    a_end: [200, 600, 300],
    b_start: [2500, 100, 800],
    b_end: [2600, 200, 900],
    // `2` is chr2 by the names rule; chr99 the assembly lacks
    a_name: ['chr1', 'chr1', 'chr99'],
    b_name: ['chr1', '2', 'chr3'],
  };

  it('each end on its sequence’s block: a pair between two sequences is a chord between their blocks', async () => {
    const s = await settle(linksSpec(Toy, rows, true));
    const f = await s.frame({ size: Size, dpr: 1 });
    const p = f.layout.sections[0].placement as BlockPlacement;
    // the pair across chr1 and chr2 is hit from inside the gap between them, above both blocks' strips
    const strip = (name: string) => {
      for (let i = 0; i < p.strips.n; i++) if (p.strips.axes[i].system.endsWith(`:name ${name})`)) return i;
      return -1;
    };
    expect(strip('chr1')).toBeGreaterThanOrEqual(0);
    const keys = new Set<string>();
    for (let x = 0; x < Size.width; x += 2)
      for (let y = 0; y < Size.height; y += 2) {
        const h = hitView(f.layout, x, y);
        if (h?.element && h.layer.endsWith('/links')) keys.add(String(h.element.key));
      }
    // rows 0 and 1 drawn (keyed by their rows, held as a complex's links are); row 2 names no sequence of the assembly
    expect([...keys].sort()).toEqual(['0', '1']);
    const said = s.model.issues.filter((i) => i.code.startsWith('seq.resolve.assembly'));
    expect(said.map((i) => i.message)).toEqual(['Rows name 1 sequences the assembly lacks (chr99)']);
  });

  it('draw as the same pairs given per sequence by a source', async () => {
    expect(await marks(linksSpec(Toy, rows, true))).toEqual(await marks(linksSpec(Toy, rows, false)));
  });

  it('a property: over generated assemblies (aliases, `include`, `order`), inline rows draw as the rows per sequence', async () => {
    for (let seed = 1; seed <= 4; seed++) {
      const r = seeded(seed * 7919);
      const n = 2 + Math.floor(r() * 4);
      const seqs: Seq[] = Array.from({ length: n }, (_, i) => ({
        name: i === n - 1 && r() < 0.5 ? 'chrM' : `chr${i + 1}`,
        length: 500 + Math.floor(r() * 4000),
        ...(r() < 0.5 ? { aliases: [`NC_${i + 1}`] } : {}),
      }));
      const nameOf = (s: Seq) => {
        const k = r();
        return k < 0.3 && s.aliases ? s.aliases[0] : k < 0.6 ? s.name.replace(/^chr/, '') : s.name;
      };
      const rows: Rows = { a_start: [], a_end: [], b_start: [], b_end: [], a_name: [], b_name: [], type: [] };
      for (let i = 0; i < 5 + Math.floor(r() * 20); i++) {
        const a = seqs[Math.floor(r() * n)];
        const b = seqs[Math.floor(r() * n)];
        const as = Math.floor(r() * (a.length - 50));
        const bs = Math.floor(r() * (b.length - 50));
        rows.a_start.push(as);
        rows.a_end.push(as + 1 + Math.floor(r() * 40));
        rows.b_start.push(bs);
        rows.b_end.push(bs + 1 + Math.floor(r() * 40));
        rows.a_name.push(nameOf(a));
        rows.b_name.push(nameOf(b));
        rows.type!.push(r() < 0.5 ? 'link' : 'contact');
      }
      const include = r() < 0.5 ? seqs.filter(() => r() < 0.7).map((s) => s.name) : undefined;
      const order = (['source', 'length', 'natural'] as const)[Math.floor(r() * 3)];
      const options = {
        ...(include?.length ? { include } : {}),
        order,
        ...(r() < 0.3 ? { layout: 'rings' as const } : {}),
      };
      const inline = await marks(linksSpec(seqs, rows, true, options));
      expect(inline, `seed ${seed}`).toEqual(await marks(linksSpec(seqs, rows, false, options)));
      expect(inline.length).toBeGreaterThan(0);
    }
  });
});

describe('inline pairs over an assembly past the gate', () => {
  // region-reads Decided 7: held and BEDPE pairs over a genome are one picture, binned by the pairs level
  const r = seeded(3);
  const rows: Rows = { a_start: [], a_end: [], b_start: [], b_end: [], a_name: [], b_name: [] };
  // chromosomes of megabases, the ends clustered at a few hot spots: the bins halve the rows (levels kept)
  const Big: Seq[] = Toy.map((t) => ({ name: t.name, length: t.length * 1000 }));
  const spot = (q: Seq) => Math.floor((Math.floor(r() * 4) + 0.5) * (q.length / 4) + r() * 2000);
  for (let i = 0; i < 400; i++) {
    const a = Big[Math.floor(r() * 3)];
    const b = Big[Math.floor(r() * 3)];
    const as = spot(a);
    const bs = spot(b);
    rows.a_start.push(as);
    rows.a_end.push(as + 5);
    rows.b_start.push(bs);
    rows.b_end.push(bs + 5);
    rows.a_name.push(a.name);
    rows.b_name.push(b.name);
  }
  const Gate = { pairs: { maxArcs: 200 } };

  it('bin by the pairs level: an arc a bin pair of the rows at the level read, never by device column', async () => {
    const s = await settle(linksSpec(Big, rows, true), [SeqRuntime], Gate);
    const f = await s.frame({ size: Size, dpr: 1 });
    const keys = new Set<string>();
    const levels = new Set<number>();
    for (let x = 0; x < Size.width; x += 3)
      for (let y = 0; y < Size.height; y += 3) {
        const h = hitView(f.layout, x, y);
        if (!h?.element || !h.layer.endsWith('/links')) continue;
        keys.add(String(h.element.key));
        if (h.bin !== undefined) levels.add(h.bin);
      }
    expect(keys.size).toBeGreaterThan(0);
    for (const k of keys) expect(k.startsWith('bin:'), k).toBe(true);
    expect(levels.size).toBe(1);
    const level = [...levels][0];
    expect(level).toBeGreaterThan(1);
    // the bin pairs of the rows at that level, on their sequences (the protocol's binning, computed apart)
    const list: PairsList = {
      form: 'list',
      aSystem: categoryColumn(rows.a_name),
      bSystem: categoryColumn(rows.b_name),
      aStart: Int32Array.from(rows.a_start),
      aEnd: Int32Array.from(rows.a_end),
      bStart: Int32Array.from(rows.b_start),
      bEnd: Int32Array.from(rows.b_end),
    };
    const bins = binPairs(list, level, { systems: true });
    expect(bins.count.length).toBeGreaterThan(Gate.pairs.maxArcs / 10);
    const arcs = recordDrawLists(f.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('path ') && l.includes(' C '));
    // an arc a bin pair (a bin reaching past its chromosome's end is a caret, a file's as well: TODO)
    expect(arcs.length).toBeGreaterThan(0);
    expect(arcs.length).toBeLessThanOrEqual(bins.count.length);
    // never by device column: no density strip, no `seq.render.pairs.limit`
    expect(f.layout.issues.map((i) => i.code)).not.toContain('seq.render.pairs.limit');
    // held and the same pairs from a file read at one level: one picture
    expect(await marks(linksSpec(Big, rows, true), Gate)).toEqual(await marks(linksSpec(Big, rows, false), Gate));
  });
});

const Small: Seq[] = [
  { name: 'chr1', length: 30 },
  { name: 'chr2', length: 20 },
  { name: 'chr3', length: 10 },
];

/** A matrix of `states` over `Small`: inline (`name` an item, in the rows' order) or a source's items per sequence. */
function matrixSpec(items: { name: string; values: number[] }[], states: number, inline: boolean): Node {
  const b = createSeqSpec({ title: 'a matrix over an assembly' });
  const g = b.data.assembly('toy', assembly('Toy', Small, { order: 'source' }));
  const names = Array.from({ length: states }, (_, i) => `s${i}`);
  let m;
  if (inline)
    m = b.data.matrix('m', {
      target: g,
      values: items.flatMap((i) => i.values),
      states: names,
      name: items.map((i) => i.name),
    });
  else {
    // the sequences with items, as the inline rows name them
    const named = Small.filter((q) => items.some((i) => canonical(Small, i.name) === q.name));
    const data: SeqData = {
      coverage: {
        systems: named.map((q) => ({ name: q.name, length: q.length })),
        shapes: named.map((q) => ({ system: q.name, shape: 'matrix' as const, states: names })),
        levels: [1],
      },
      pairs: named.map((q) => {
        const own = items.filter((i) => canonical(Small, i.name) === q.name);
        return {
          request: { system: q.name, start: 0, end: q.length, level: 1, shapes: ['matrix' as const] },
          answer: {
            ok: true as const,
            tables: {
              matrix: {
                form: 'items' as const,
                start: Int32Array.from(own, (_, j) => j),
                end: Int32Array.from(own, (_, j) => j + 1),
                states,
                value: Float32Array.from(own.flatMap((i) => i.values)),
              },
            },
          },
        };
      }),
    };
    const src = b.data.source('src', { format: 'seq', value: seqToValue(data) as never });
    m = b.data.matrix('m', { target: g, source: src });
  }
  b.views.view('v', (v) =>
    genome(v, 'genome', g, {}, (s) => s.track('m', { height: 60 }).layer('m', { type: 'matrix', data: m })),
  );
  return b.build();
}

describe('an inline matrix over an assembly', () => {
  it('a sequence’s items its positions in order (`name` an item), drawn as the same items given per sequence', async () => {
    const r = seeded(11);
    const items: { name: string; values: number[] }[] = [];
    // chr1 and chr3 by any of their names, interleaved; chr2 none
    for (let i = 0; i < 18; i++)
      items.push({ name: i % 3 === 2 ? '3' : i % 2 ? 'chr1' : '1', values: [r(), r(), r()] });
    const inline = await marks(matrixSpec(items, 3, true));
    expect(inline).toEqual(await marks(matrixSpec(items, 3, false)));
    expect(inline.split('\n').filter((l) => l.startsWith('rect') || l.startsWith('image')).length).toBeGreaterThan(0);
  });
});

const Empty = {
  form: 'sparse' as const,
  x: 0,
  y: 0,
  xBin: new Int32Array(0),
  yBin: new Int32Array(0),
  value: new Float32Array(0),
  fill: Number.NaN,
};

/** A grid over `Small` both ways: inline (`x_name`/`y_name` a value) or a source's grid per sequence pair. */
function gridSpec(cells: { x: string; y: string; values: number[] }[], inline: boolean): Node {
  const b = createSeqSpec({ title: 'a grid over an assembly' });
  const g = b.data.assembly('toy', assembly('Toy', Small, { order: 'source' }));
  let grid;
  if (inline)
    grid = b.data.grid('h', {
      x: g,
      values: cells.flatMap((c) => c.values),
      x_name: cells.flatMap((c) => c.values.map(() => c.x)),
      y_name: cells.flatMap((c) => c.values.map(() => c.y)),
    });
  else {
    const xs = [...new Set(cells.map((c) => canonical(Small, c.x)!))];
    const ys = [...new Set(cells.map((c) => canonical(Small, c.y)!))];
    const length = (n: string) => Small.find((q) => q.name === n)!.length;
    const pairs: SeqPair[] = [];
    for (const x of xs)
      for (const y of ys) {
        const c = cells.find((k) => canonical(Small, k.x) === x && canonical(Small, k.y) === y);
        const nx = length(x);
        const ny = length(y);
        pairs.push({
          request: { system: x, start: 0, end: nx, level: 1, shapes: ['grid'], y: { system: y, start: 0, end: ny } },
          answer: {
            ok: true,
            tables: { grid: c ? { form: 'dense', x: 0, y: 0, nx, ny, value: Float32Array.from(c.values) } : Empty },
          },
        });
      }
    const data: SeqData = {
      coverage: {
        systems: [...new Set([...xs, ...ys])].map((name) => ({ name, length: length(name) })),
        shapes: xs.map((system) => ({ system, shape: 'grid' as const, y: ys })),
        levels: [1],
      },
      pairs,
    };
    const src = b.data.source('src', { format: 'seq', value: seqToValue(data) as never });
    grid = b.data.grid('h', { x: g, source: src });
  }
  b.views.view('v', (v) => {
    const top = genome(v, 'genome', g);
    v.section('plot', { layout: 'grid', follow: top.address }).track('plot', {}, (t) => {
      t.layer('h', { type: 'grid', data: grid });
    });
  });
  return b.build();
}

describe('an inline grid over an assembly', () => {
  it('a sequence pair’s values one grid (`x_name`/`y_name` a value), drawn as the same grids given per pair', async () => {
    const r = seeded(5);
    const cell = (x: string, y: string, nx: number, ny: number) => ({
      x,
      y,
      values: Array.from({ length: nx * ny }, () => Math.round(r() * 100)),
    });
    const cells = [cell('chr1', 'chr1', 30, 30), cell('1', '2', 30, 20), cell('chr3', 'chr2', 10, 20)];
    const inline = await marks(gridSpec(cells, true));
    expect(inline).toEqual(await marks(gridSpec(cells, false)));
    expect(inline).toContain('image');
  });
});
