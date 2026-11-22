/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { type Coverage, type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import type { BlockPlacement } from '../layout/blocks.js';
import type { Frame } from '../renderer.js';
import { hitView } from '../tracks-hit.js';
import { recordDrawLists } from '../record.js';
import { toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';
import { settleNavigator } from './navigator-fixtures.js';

// Data over an assembly read by sequence (_spec/README.md "Resolve", Assemblies: data over an assembly): each block of
// a section of the assembly reads its own sequence's rows, inline (split by `target`) or from a source by region.

beforeAll(async () => {
  await loadGenome();
});

const Toy = [
  { name: 'chr1', length: 3000 },
  { name: 'chr2', length: 2000 },
  { name: 'chr3', length: 1000 },
];

const Size = { width: 900, height: 200 };

/** The keys hit along each strip's band `band`, by strip (its sequence's name). */
function keysByStrip(f: Frame, band: number): Record<string, string[]> {
  const sec = f.layout.sections[0];
  const p = sec.placement as BlockPlacement;
  const s = p.strips;
  const out: Record<string, string[]> = {};
  for (let i = 0; i < s.n; i++) {
    const a = s.axes[i];
    const keys = new Set<string>();
    const y0 = s.bandY[i * s.bands + band];
    for (let x = a.x0; x < a.x0 + a.width; x += 1)
      for (let y = y0; y < y0 + s.bandH[i * s.bands + band]; y += 3) {
        const h = hitView(f.layout, x, y);
        if (h?.element && h.layer.endsWith('/bed')) keys.add(String(h.element.key));
      }
    out[a.system.replace(/^.*:name /, '').replace(/\)$/, '')] = [...keys].sort();
  }
  return out;
}

/** A toy assembly's section with a `bed` ranges layer over `rows` (inline) or `source` (a seq source's value). */
function bedSpec(over: { rows?: Record<string, unknown>; source?: SeqData; include?: string[] }): Node {
  const b = createSeqSpec({ title: 'bed over an assembly' });
  const g = b.data.assembly('toy', assembly('Toy', Toy, { order: 'source' } as never));
  const src = over.source && b.data.source('src', { format: 'seq', value: seqToValue(over.source) as never });
  const bed = b.data.ranges('bed', { target: g, ...(over.rows ? { rows: over.rows } : { source: src }) } as never);
  b.views.view('v', (v) =>
    genome(v, 'genome', g, over.include ? { include: over.include } : {}, (s) =>
      s.track('bed').layer('bed', { type: 'ranges', data: bed }),
    ),
  );
  return b.build();
}

/** A source listing chr1 and chr3 as `1` and `3`, and a name the assembly lacks. */
const Listed: SeqData = {
  coverage: {
    systems: [
      { name: '1', length: 3000 },
      { name: '3', length: 1000 },
      { name: 'chrUn', length: 50 },
    ],
    shapes: ['1', '3', 'chrUn'].map((system) => ({ system, shape: 'ranges' as const })),
    levels: [1],
  },
  pairs: [
    { system: '1', end: 3000, id: ['p', 'q'], start: [10, 1000], stop: [200, 2000] },
    { system: '3', end: 1000, id: ['r'], start: [500], stop: [900] },
    { system: 'chrUn', end: 50, id: ['u'], start: [0], stop: [10] },
  ].map((x) => ({
    request: { system: x.system, start: 0, end: x.end, level: 1, shapes: ['ranges' as const] },
    answer: {
      ok: true as const,
      tables: {
        ranges: {
          form: 'list' as const,
          start: Int32Array.from(x.start),
          end: Int32Array.from(x.stop),
          id: categoryColumn(x.id),
        },
      },
    },
  })),
};

describe('ranges over an assembly', () => {
  it('inline rows are split by their `target` through the names rule, each drawn on its block (EC-189)', async () => {
    const s = await settle(
      bedSpec({
        rows: {
          start: [100, 500, 200, 10, 0],
          end: [400, 900, 600, 300, 10],
          id: ['a', 'b', 'c', 'd', 'x'],
          target: ['chr1', '1', '2', 'chr3', 'chr99'],
        },
      }),
    );
    const f = await s.frame({ size: Size });
    expect(keysByStrip(f, 1)).toEqual({ chr1: ['a', 'b'], chr2: ['c'], chr3: ['d'] });
    // a row naming no sequence is counted, once
    const said = s.model.issues.filter((i) => i.code === 'seq.resolve.assembly.name');
    expect(said.map((i) => i.message)).toEqual(['Rows name 1 sequences the assembly lacks (chr99)']);
    expect(s.model.issues.filter((i) => i.code === 'seq.resolve.param.unread')).toEqual([]);
  });

  it('a source is read by region on each block, its names mapped both ways (`1` is `chr1`, EC-187)', async () => {
    const s = await settle(bedSpec({ source: Listed }), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    expect(keysByStrip(f, 1)).toEqual({ chr1: ['p', 'q'], chr2: [], chr3: ['r'] });
    const said = s.model.issues.filter((i) => i.code === 'seq.resolve.assembly.name');
    expect(said.map((i) => i.message)).toEqual(['The source names 1 sequences the assembly lacks (chrUn)']);
  });

  it('the layer’s wants over every block are bounded together (`maxTilesPerLayer` summed over blocks): past it, zoom in', async () => {
    const wide = await settle(bedSpec({ source: Listed }), [SeqRuntime], {}, {}, { maxTilesPerLayer: 2 });
    const f = await wide.frame({ size: Size });
    // two blocks read (chr1, chr3), a tile each at least: past two together
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.zoom_in');
    expect(keysByStrip(f, 1)).toEqual({ chr1: [], chr2: [], chr3: [] });
    const fits = await settle(bedSpec({ source: Listed }), [SeqRuntime], {}, {}, { maxTilesPerLayer: 64 });
    expect(keysByStrip(await fits.frame({ size: Size }), 1)).toEqual({ chr1: ['p', 'q'], chr2: [], chr3: ['r'] });
  });

  it('a sequence the section’s `include` leaves out reads nothing', async () => {
    const s = await settle(
      bedSpec({
        rows: { start: [100, 200], end: [400, 600], id: ['a', 'c'], target: ['chr1', 'chr2'] },
        include: ['chr2', 'chr3'],
      }),
    );
    const f = await s.frame({ size: Size });
    expect(keysByStrip(f, 1)).toEqual({ chr2: ['c'], chr3: [] });
  });
});

describe('the ideogram over an assembly', () => {
  /** A toy assembly's section with cytobands over every chromosome (or `only` those) and an ideogram navigator. */
  function ideogram(only?: string): Node {
    const rows = {
      start: [] as number[],
      end: [] as number[],
      id: [] as string[],
      stain: [] as string[],
      target: [] as string[],
    };
    for (const c of Toy) {
      if (only && c.name !== only) continue;
      for (let i = 0; i < 4; i++) {
        rows.start.push((c.length * i) / 4);
        rows.end.push((c.length * (i + 1)) / 4);
        rows.id.push(`${c.name}p${i}`);
        rows.stain.push(['gneg', 'gpos50', 'acen', 'gpos100'][i]);
        // the bands file's names without the prefix: `1` is `chr1`
        rows.target.push(c.name.replace('chr', ''));
      }
    }
    const b = createSeqSpec({ title: 'ideogram over an assembly' });
    const g = b.data.assembly('toy', assembly('Toy', Toy, { order: 'source' } as never));
    const bands = b.data.ranges('bands', { target: g, rows } as never);
    b.views.view('v', (v) => genome(v, 'genome', g, { bands }));
    return b.build();
  }

  const rects = async (spec: Node) => {
    const s = await settleNavigator(spec);
    const { nav } = await s.frames({ size: Size }, { size: { width: 900, height: 24 } });
    return recordDrawLists(nav.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('rect ')).length;
  };

  it('reads its bands by sequence: every chromosome’s bands drawn (GC3’s follow-up)', async () => {
    const all = await rects(ideogram());
    const one = await rects(ideogram('chr1'));
    expect(one).toBeGreaterThan(0);
    expect(all).toBeGreaterThanOrEqual(3 * one - 2);
  });
});

describe('data over 3,000 contigs at the fit', () => {
  // medium: 3,000 blocks resolved and laid out (seconds); the quick lane's representatives are the 3-block cases above
  it('the small-blocks rule holds: no read on a contig under a device pixel, the layers within their bounds', {
    tags: ['medium'],
  }, async () => {
    const sequences = toyGenome(24, 2976, 3);
    const names = sequences.map((q) => q.name);
    // a feature and a link on every sequence, as one `seq` source over all of them
    const coverage: Coverage = {
      systems: sequences.map((q) => ({ name: q.name, length: q.length })),
      shapes: names.flatMap((system) => [
        { system, shape: 'ranges' as const },
        { system, shape: 'pairs' as const },
      ]),
      levels: [1, 65_536],
    };
    // each sequence's rows whole at both levels (a source's coarse level: the fit reads it, never level 1)
    const pairs: SeqData['pairs'] = [1, 65_536].flatMap((level) =>
      sequences.map((q, i) => {
        const other = names[(i * 7) % names.length];
        return {
          request: { system: q.name, start: 0, end: q.length, level, shapes: ['ranges', 'pairs'] },
          answer: {
            ok: true,
            tables: {
              ranges: {
                form: 'list',
                start: Int32Array.from([0]),
                end: Int32Array.from([10]),
                id: categoryColumn([`${q.name}.f`]),
              },
              pairs: {
                form: 'list',
                aSystem: categoryColumn([q.name]),
                bSystem: categoryColumn([other]),
                aStart: Int32Array.from([0]),
                aEnd: Int32Array.from([10]),
                bStart: Int32Array.from([0]),
                bEnd: Int32Array.from([10]),
              },
            },
          },
        };
      }),
    );
    const b = createSeqSpec({ title: '3,000 contigs' });
    const g = b.data.assembly('toy', assembly('Toy', sequences, { order: 'source' } as never));
    const src = b.data.source('src', { format: 'seq', value: seqToValue({ coverage, pairs }) as never });
    const bed = b.data.ranges('bed', { target: g, source: src } as never);
    const links = b.data.pairs('links', { a: g, source: src } as never);
    b.views.view('v', (v) =>
      genome(v, 'genome', g, {}, (s) => {
        s.track('bed').layer('bed', { type: 'ranges', data: bed });
        s.track('links').layer('links', { type: 'pairs', data: links });
      }),
    );
    const s = await settle(b.build(), [SeqRuntime]);
    const f = await s.frame({ size: Size, dpr: 1 });
    const p = f.layout.sections[0].placement as BlockPlacement;
    const read = new Set<string>();
    for (let i = 0; i < p.strips.n; i++) read.add(p.strips.axes[i].system.replace(/^.*:name /, '').replace(/\)$/, ''));
    // most contigs are under a device pixel: no strip, so no read
    expect(p.strips.n).toBeLessThan(500);
    const asked = new Set(f.wants.map((w) => w.system));
    for (const name of asked) expect(read.has(name), name).toBe(true);
    expect(asked.size).toBeLessThanOrEqual(p.strips.n);
    expect(f.layout.issues.map((i) => i.code)).not.toContain('seq.render.zoom_in');
    expect(s.store.issues.value.map((i) => i.code)).not.toContain('seq.store.too_many_tiles');
  });
});
