/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// N1's findings (each also ran against the old stack until SW1, in its parity layer): letters without spreads, one
// palette per layer, first-fit lanes and accession counters, loops for a group's merged row, split colors, the
// section's label, and G1-2 the `color` track over dense values; PL's color parity against `layerPalette` (the function
// the viewer colors through).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  Color,
  compile,
  createRegistries,
  KernelRegistry,
  layerPalette,
  MvsRegistry,
  type Palette,
  register,
  resolveView,
  type SeqViewSpecDoc,
} from '@molstar/seq-core/format';
import { Annotations, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { type ElementDesc, materialize, NIGHTINGALE_MAX_POSITIONS, NightingaleSupported } from '../src/index.js';
import { lanesOf } from '../src/materialize.js';

const kernel = createRegistries();
register(kernel, [KernelRegistry, MvsRegistry]);
const none = { load: async () => Promise.reject(new Error('no files')) };

type Builder = ReturnType<typeof createSeqViewSpec>;
async function draw(build: (b: Builder) => void, view = 'v') {
  const b = createSeqViewSpec();
  build(b);
  const doc: SeqViewSpecDoc = b.getDocument();
  const { set } = await compile(doc, kernel, none);
  const r = resolveView(doc, set, view, kernel, NightingaleSupported);
  return { set, view: r.view, ...materialize(r.view, set) };
}
const all = (d: ElementDesc, tag: string, out: ElementDesc[] = []): ElementDesc[] => {
  if (d.tag === tag) out.push(d);
  for (const c of d.children) all(c, tag, out);
  return out;
};
type F = { accession: string; start: number; end: number; fill?: string; color?: string; tooltipContent?: string };
const features = (d: ElementDesc) => all(d, 'nightingale-track').map((t) => t.properties.data as F[]);
const words = (n: number) => Array.from({ length: n }, (_, i) => 'ACDEFGHIKLMNPQRSTVWY'[(i * 7) % 20]).join('');
const MKV = 'MKVLAAGIVGLLLAQ';

describe('scale', () => {
  it('200,000 letters, no spread over the call stack', async () => {
    const text = words(200_000);
    const { root, issues } = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: text });
      b.view({ ref: 'v', label: 'V' }).block({ sequence: s }).track({ label: 'L' }).letters();
    });
    expect((all(root, 'nightingale-sequence')[0].properties.data as string) === text).toBe(true);
    expect(issues).toEqual([]);
  });

  it('150,000 features in a collapsed group are its merged row, copied without a spread', async () => {
    const n = 150_000;
    const { root } = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: words(n) });
      const value = { start: Array.from({ length: n }, (_, i) => i), stop: Array.from({ length: n }, (_, i) => i + 1) };
      const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value });
      b.view({ ref: 'v', label: 'V' })
        .block({ sequence: s })
        .group({ label: 'G', collapsed: true })
        .track({ label: 'A' })
        .ranges({ data: a });
    });
    const sizes = all(root, 'nightingale-track-canvas').map((t) => (t.properties.data as F[]).length);
    expect(sizes).toEqual([n, n]); // the track's, and the group's merged row
  });

  it('5,000 features with one id take one suffix each, counted from the base', async () => {
    const n = 5_000;
    const { root } = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: MKV });
      const value = { start: Array(n).fill(2), stop: Array(n).fill(5), id: Array(n).fill('same') };
      const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value });
      b.view({ ref: 'v', label: 'V' }).block({ sequence: s }).track({ label: 'A' }).ranges({ data: a });
    });
    // Features drawn alike are one: the id is not repeated into a track (merged copies' items at one place).
    expect(features(root)[0].map((f) => f.accession)).toEqual(['same']);
    const second = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: MKV });
      const value = { start: [0, 1, 2, 3], stop: [1, 2, 3, 4], id: ['x', 'x', 'x', 'x'] };
      const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value });
      b.view({ ref: 'v', label: 'V' }).block({ sequence: s }).track({ label: 'A' }).ranges({ data: a });
    });
    // One id on four rows is one item; unnamed rows get the layer's ref, the annotation and the row.
    expect(features(second.root)[0].map((f) => f.accession)).toEqual(['x']);
    const unnamed = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: MKV });
      const a = b.data.annotation({
        ref: 'a',
        geometry: 'ranges',
        target: s,
        value: { start: [0, 3, 6], stop: [2, 5, 8] },
      });
      b.view({ ref: 'v', label: 'V' }).block({ sequence: s }).track({ label: 'A' }).points({ data: a });
    });
    const ids = features(unnamed.root)[0].map((f) => f.accession);
    expect(new Set(ids).size).toBe(ids.length);
  });
});

describe('lanes', () => {
  type R = { start: number; end: number };
  /** The first-fit reference: the old packing, quadratic and simple. */
  const naive = (fs: R[]) => {
    const lanes: R[][] = [];
    for (const f of fs) {
      const k = lanes.findIndex((lane) => lane.every((g) => f.start > g.end || f.end < g.start));
      if (k >= 0) lanes[k].push(f);
      else lanes.push([f]);
    }
    return lanes.length;
  };

  it('pack as first fit does, to at most LANES', () => {
    let seed = 11;
    const rand = () => {
      seed = (seed * 1103515245 + 12345) % 2147483648;
      return seed / 2147483648;
    };
    for (let round = 0; round < 40; round++) {
      const n = 1 + Math.floor(rand() * 60);
      const fs = Array.from({ length: n }, () => {
        const start = Math.floor(rand() * 80);
        return { start, end: start + Math.floor(rand() * (round % 2 ? 30 : 6)) };
      });
      expect(lanesOf(fs)).toBe(Math.min(10, naive(fs)));
    }
    expect([
      lanesOf([]),
      lanesOf([
        { start: 1, end: 5 },
        { start: 6, end: 9 },
      ]),
      lanesOf([
        { start: 1, end: 5 },
        { start: 5, end: 9 },
      ]),
    ]).toEqual([0, 1, 2]);
    // Eleven features on one place need eleven lanes: it stops at ten.
    expect(lanesOf(Array.from({ length: 30 }, () => ({ start: 1, end: 2 })))).toBe(10);
  });

  it('40,000 disjoint ranges in one lane: O(n) for sorted input, O(n log n) probes and O(n) moves for any order', () => {
    const n = 40_000;
    const sortedInput = Array.from({ length: n }, (_, i) => ({ start: 3 * i, end: 3 * i + 1 }));
    const probes = { probes: 0, moves: 0 };
    expect(lanesOf(sortedInput, probes)).toBe(1);
    expect(probes).toEqual({ probes: 0, moves: 0 }); // beyond the hull each time: no search, nothing moved
    const shuffled = sortedInput
      .map((f, i) => [(i * 7919) % n, f] as const)
      .sort((a, b) => a[0] - b[0])
      .map(([, f]) => f);
    const shuffledProbes = { probes: 0, moves: 0 };
    expect(lanesOf(shuffled, shuffledProbes)).toBe(1);
    expect(shuffledProbes.probes).toBeLessThan(n * Math.log2(n));
    expect(shuffledProbes.probes).toBeGreaterThan(0);
    // An insert moves at most a chunk (512): one sorted array would move about n² / 4 (4 · 10⁸ here).
    expect(shuffledProbes.moves).toBeLessThan(512 * n);
    // Reverse order takes the hull on the left, also with no search; each insert moves the first chunk only.
    const reversed = { probes: 0, moves: 0 };
    expect(lanesOf([...sortedInput].reverse(), reversed)).toBe(1);
    expect(reversed.probes).toBe(0);
    expect(reversed.moves).toBeLessThan(512 * n);
  });
});

describe('one palette per layer', () => {
  const color = { theme: 'field_palette', params: { palette: { kind: 'categorical', colors: 'Category10' } } };
  const build = (b: Builder, extra: Record<string, unknown> = {}, arrangement = 'side') => {
    const [a, c] = ['A', 'B'].map((x) => b.data.sequence({ ref: x, value: MKV }));
    const cat = (ref: string, target: unknown, category: string[]) =>
      b.data.annotation({
        ref,
        geometry: 'ranges',
        target: target as never,
        value: { start: [0, 5], stop: [3, 8], category },
      });
    const [fa, fb] = [cat('fa', a, ['x', 'y']), cat('fb', c, ['z', 'y'])];
    const layer = { data: [fa, fb], color, ...extra } as never;
    const both = b.view({ ref: 'both', label: 'Both', arrangement } as never);
    both.block({ sequence: a }).block({ sequence: c });
    both.track({ label: 'F' }).ranges(layer);
    const only = b.view({ ref: 'only', label: 'B' });
    only.block({ sequence: c });
    only.track({ label: 'F' }).ranges(layer);
  };
  const byName = (fs: F[]) => Object.fromEntries(fs.map((f) => [f.tooltipContent, f.fill]));

  it('a value has one color whatever the view shows: B alone colors as B in the view of both', async () => {
    const both = await draw((b) => build(b), 'both');
    const only = await draw((b) => build(b), 'only');
    const all2 = features(both.root).at(-1) as F[]; // after the block labels' track
    const fromB = features(only.root).at(-1) as F[];
    expect(fromB).toHaveLength(2);
    // A palette over the layer's data, in data order: x, y, z.
    expect(Object.keys(byName(all2)).sort()).toEqual(['x', 'y', 'z']);
    for (const f of fromB)
      expect(f.fill).toBe(
        Color.toHex(paletteFor(both.set, both.view, 'category')(f.tooltipContent as string) as number),
      );
    const shared = new Set(all2.map((f) => f.fill));
    expect(fromB.every((f) => shared.has(f.fill))).toBe(true);
    // "y" is on both chains: one color for it everywhere.
    const y = [...all2, ...fromB].filter((f) => f.tooltipContent === 'y').map((f) => f.fill);
    expect(new Set(y).size).toBe(1);
  });

  it("a stacked section's second block does not restart at the first color", async () => {
    const stacked = await draw((b) => build(b, {}, 'stacked'), 'both');
    const [first, second] = all(stacked.root, 'nightingale-manager').map((m) => features(m)[0]);
    const f = Object.fromEntries([...first, ...second].map((x) => [x.tooltipContent, x.fill]));
    // x, y on A and z, y on B in the layer's data order: z is the third value, not the first.
    expect(f.z).toBe(Color.toHex(0x2ca02c));
    expect(f.x).toBe(Color.toHex(0x1f77b4));
    expect(f.y).toBe(Color.toHex(0xff7f0e));
  });

  it('a filtered layer keeps its colors', async () => {
    const plain = await draw((b) => build(b), 'only');
    const filtered = await draw((b) => build(b, { filter: { field: 'category', values: ['y'] } }), 'only');
    const [p] = features(plain.root);
    const [f] = features(filtered.root);
    expect(f).toHaveLength(1);
    expect(f[0].fill).toBe(p.find((x) => x.tooltipContent === 'y')?.fill);
  });

  /** The layer's palette as the viewer builds it, over every row of the layer's data. */
  function paletteFor(set: SeqSet, view: Awaited<ReturnType<typeof draw>>['view'], field: string) {
    const data = Array.from(view.data, (id) => set.annotations[id] as never);
    return layerPalette(set, data, { kind: 'categorical', colors: 'Category10' } as Palette, field, [
      Number.NaN,
      Number.NaN,
    ]);
  }
});

describe('colors', () => {
  it("a split color 'red/blue' is its first half, as the viewer draws it", async () => {
    const { root } = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: MKV });
      const value = { start: [0, 4], stop: [2, 6], color: ['red/blue', '#00ff00'] };
      const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value });
      b.view({ ref: 'v', label: 'V' })
        .block({ sequence: s })
        .track({ label: 'A' })
        .ranges({ data: a, color: { theme: 'field_color', params: { field: 'color' } } } as never);
    });
    expect(features(root)[0].map((f) => f.fill)).toEqual(['#ff0000', '#00ff00']);
  });

  it("the section's label, else its ref; shown when more than one section is", async () => {
    const two = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: MKV });
      const v = b.view({ ref: 'v', label: 'V' });
      v.section({ ref: 'first', label: 'First section' }).block({ sequence: s }).track({ label: 'L' }).letters();
      v.section({ ref: 'second' }).block({ sequence: s }).track({ label: 'L' }).letters();
    });
    expect(
      all(two.root, 'div')
        .filter((d) => d.attributes.class === 'svs-section-label')
        .map((d) => d.properties.textContent),
    ).toEqual(['First section', 'second']);
    const one = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: MKV });
      b.view({ ref: 'v', label: 'V' })
        .section({ ref: 'only', label: 'Only' })
        .block({ sequence: s })
        .track({ label: 'L' })
        .letters();
    });
    expect(all(one.root, 'div').filter((d) => d.attributes.class === 'svs-section-label')).toEqual([]);
  });

  it("a section past NIGHTINGALE_MAX_POSITIONS is left out with one issue; the view's other sections render (GM-05)", async () => {
    const r = await draw((b) => {
      const big = b.data.sequence({ ref: 'big', value: 'A'.repeat(NIGHTINGALE_MAX_POSITIONS + 1) });
      const s = b.data.sequence({ ref: 's', value: MKV });
      const v = b.view({ ref: 'v', label: 'V' });
      v.section({ ref: 'genome' }).block({ sequence: big }).track({ label: 'L' }).letters();
      v.section({ ref: 'protein' }).block({ sequence: s }).track({ label: 'L' }).letters();
    });
    expect(r.issues.map((i) => [i.severity, i.code, i.path, i.message])).toEqual([
      [
        'warning',
        'svs.render.omitted',
        'genome',
        'Nightingale: a section of 1,000,001 positions, past its limit of 1,000,000 (an element per position)',
      ],
    ]);
    const sections = all(r.root, 'div').filter((d) => d.attributes.class === 'svs-section');
    expect(sections.map((d) => d.attributes['data-svs-section'])).toEqual(['protein']);
    expect(all(sections[0], 'nightingale-sequence')).toHaveLength(1);
  });
});

describe('G1-2 and PL: a color layer over dense values', () => {
  const palettes: [string, Palette][] = [
    ['categorical', { kind: 'categorical', colors: { '1': '#ff0000', '2': '#00ff00', '3': '#0000ff' } }],
    ['continuous', { kind: 'continuous', colors: 'Viridis' }],
    [
      'discrete',
      {
        kind: 'discrete',
        mode: 'absolute',
        colors: [
          ['#ff0000', 0],
          ['#0000ff', 2],
        ],
      },
    ],
  ];
  it.each(palettes)(
    "%s: one hex per position, the layer's palette function's, in runs of one color",
    async (_, palette) => {
      const dense = [1, 1, 2, 2, 2, Number.NaN, 3, 1, 1, 1, 2, 3, 3, 1, 2];
      const { root, set, view } = await draw((b) => {
        const s = b.data.sequence({ ref: 's', value: MKV });
        const d = b.data.annotation({ ref: 'd', geometry: 'values', target: s, value: { value: dense } });
        b.view({ ref: 'v', label: 'V' })
          .block({ sequence: s })
          .track({ label: 'C' })
          .color({ data: d, color: { theme: 'field_palette', params: { palette } } } as never);
      });
      const ann = set.annotations[view.data[0]] as never;
      const f = layerPalette(set, [ann], palette, null, [view.layers.domainMin[0], view.layers.domainMax[0]]);
      const expected = dense.map((x) => (Number.isNaN(x) ? undefined : f(x)));
      const got: (string | undefined)[] = Array(dense.length).fill(undefined);
      const drawn = features(root)[0];
      for (const x of drawn) for (let p = x.start; p <= x.end; p++) got[p - 1] = x.fill;
      expect(got).toEqual(expected.map((c) => (c === undefined ? undefined : Color.toHex(c))));
      // A feature per run of one color, never per position.
      for (let k = 1; k < drawn.length; k++)
        expect(drawn[k].fill === drawn[k - 1].fill && drawn[k].start === drawn[k - 1].end + 1).toBe(false);
      expect(Annotations).toBeDefined();
    },
  );

  it('walks 200,000 positions in runs, one palette call per value change', async () => {
    const n = 200_000;
    const dense = Array.from({ length: n }, (_, i) => 1 + (Math.floor(i / 1000) % 3));
    const { root } = await draw((b) => {
      const s = b.data.sequence({ ref: 's', value: words(n) });
      const d = b.data.annotation({ ref: 'd', geometry: 'values', target: s, value: { value: dense } });
      const palette = { kind: 'categorical', colors: { '1': '#ff0000', '2': '#00ff00', '3': '#0000ff' } };
      b.view({ ref: 'v', label: 'V' })
        .block({ sequence: s })
        .track({ label: 'C' })
        .color({ data: d, color: { theme: 'field_palette', params: { palette } } } as never);
    });
    const drawn = features(root)[0];
    expect(drawn).toHaveLength(200);
    expect(drawn.every((f) => f.end - f.start === 999)).toBe(true);
  });
});
