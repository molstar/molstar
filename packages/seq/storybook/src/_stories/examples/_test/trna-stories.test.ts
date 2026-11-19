/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { AlignmentFileRuntime } from '@molstar/seqx-core/data/alignment-readers';
import { theme } from '@molstar/seqx-core/lang';
import {
  DarkTokens,
  type GridLay,
  LightTokens,
  type PairsLay,
  type RangesLay,
  recordDrawLists,
  type RowsLay,
  viewLegends,
} from '@molstar/seqx-core/render';
import { TrnaViews } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { DerivedData } from '@molstar/seqx-examples/trna/derived-data';
import { without } from '@molstar/seqx-examples/trna/trna';
import type { UiHistory } from '@molstar/seqx-viewer/ui-history';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

// The tRNA dataset's views, headless, drawn against the data: 1EHZ's anchors in the arcs' lay (the four stems, the six
// tertiary pairs, the kinds), its modified residues' ranges, the dot plot's two triangles, both themes in the SVG, the
// seed row and the seed; and the editing views, each a patch undone (a kind of pairs hidden, the covariation's palette).
// The interaction (hover, click, the crosshair and the rectangle) is in `pairs-grids-stories.test.ts`.

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const trna = (id: string) => TrnaViews.find((v) => v.id === id)!;
const arcsOf = (v: SeqViewer, section = 0) => v.frame!.layout.sections[section].sectionLays[0].lay as PairsLay;
const keysOf = (lay: PairsLay) => {
  const keys = lay.prep!.data!.keys as string[];
  return Array.from({ length: lay.n }, (_, k) => keys[lay.pair[k]]).sort();
};
/** A pair id `i-j` from 1-based positions. */
const ids = (pairs: readonly (readonly [number, number])[]) => pairs.map(([i, j]) => `${i}-${j}`);

describe('tRNA: 1EHZ’s arcs against its anchors', () => {
  it('the stems nest (none crosses another); of the tertiary pairs four cross the stems, 26–44 crosses nothing, 54–58 only tertiary pairs', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const keys = keysOf(arcsOf(viewer));
    const pairs = keys.map((k) => k.split('-').map(Number) as [number, number]);
    const cross = (p: number[], q: number[]) =>
      (p[0] < q[0] && q[0] < p[1] && p[1] < q[1]) || (q[0] < p[0] && p[0] < q[1] && q[1] < p[1]);
    const tertiary = new Set(['8-14', '15-48', '18-55', '19-56', '26-44', '54-58']);
    const crossing = (k: string) => pairs.filter((q) => cross(k.split('-').map(Number), q)).map((q) => q.join('-'));
    // stems (and the loop’s 32–38) cross no stem: they cross only the tertiary pairs
    for (const k of keys.filter((k) => !tertiary.has(k)))
      for (const c of crossing(k)) expect(tertiary.has(c)).toBe(true);
    expect(crossing('26-44')).toEqual([]);
    expect(crossing('54-58').sort()).toEqual(['18-55', '19-56']);
    for (const k of ['8-14', '15-48', '18-55', '19-56']) expect(crossing(k).some((c) => !tertiary.has(c))).toBe(true);
  });

  it('the four stems nested, the six tertiary pairs among them, every pair drawn once', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const keys = keysOf(arcsOf(viewer));
    expect(keys).toHaveLength(27);
    // the acceptor stem 1–7 : 66–72, the D stem 10–13 : 22–25, the anticodon stem 27–30 : 40–43, the T stem 49–53 : 61–65
    const stem = (a: number, b: number, n: number, ...skip: number[]) =>
      Array.from({ length: n }, (_, k) => [a + k, b - k] as const).filter(([i]) => !skip.includes(i));
    for (const [a, b, n] of [
      [1, 72, 7],
      [10, 25, 4],
      [27, 43, 4],
      [49, 65, 5],
    ] as const)
      for (const id of ids(stem(a, b, n))) expect(keys).toContain(id);
    // the entry's table lacks the anticodon stem's 31–39, and has the loop's 32–38 instead
    expect(keys).not.toContain('31-39');
    expect(keys).toContain('32-38');
    for (const id of ids([
      [8, 14],
      [15, 48],
      [18, 55],
      [19, 56],
      [26, 44],
      [54, 58],
    ]))
      expect(keys).toContain(id);
  });

  it('the kinds: one wobble (4–69), the tertiary and loop pairs non-canonical but 19–56', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const d = arcsOf(viewer).prep!.data!;
    const kind = (id: string) => d.kinds[d.kind[(d.keys as string[]).indexOf(id)]];
    expect(kind('4-69')).toBe('wobble');
    expect(kind('15-48')).toBe('non-canonical');
    expect(kind('19-56')).toBe('canonical');
    expect(kind('1-72')).toBe('canonical');
  });

  it('the modified residues: 14 ranges over the letters, each its code, the anticodon’s Gm34 among them', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const section = viewer.frame!.layout.sections[0];
    const track = section.tracks.find((t) => t.track.address.endsWith('/modified'))!;
    const lay = track.layers[0].lay as RangesLay;
    expect(lay.keys).toHaveLength(14);
    expect([...lay.keys]).toContain('OMG34');
    expect([...lay.keys]).toContain('H2U16');
    expect([...lay.keys]).toContain('PSU55');
  });

  for (const t of ['light', 'dark'] as const)
    it(`${t}: the SVG draws every arc, with no issue, and differs from the other theme’s`, async () => {
      const ctx = createSeqContext({
        registry: [{ kind: 'format-runtime', runtime: SeqRuntime }],
        ports: { assets: () => nodeFixtureAssets() },
      });
      const r = await ctx.renderSvg(trna('arcs').build(), {
        size: { width: 720, height: 320 },
        tokens: t === 'dark' ? DarkTokens : LightTokens,
      });
      expect(r.issues).toEqual([]);
      const paths = recordDrawLists(r.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('path'));
      expect(paths).toHaveLength(27);
      const other = await ctx.renderSvg(trna('arcs').build(), {
        size: { width: 720, height: 320 },
        tokens: t === 'light' ? DarkTokens : LightTokens,
      });
      expect(other.svg).not.toBe(r.svg);
      ctx.retire();
    });
});

describe('tRNA: the arcs’ legend and the modified residues’ codes', () => {
  it('the legend lists the pair kinds under the arcs’ track; zoomed in, each modified residue is labeled by its code', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const arcs = viewLegends(viewer.frame!.layout).legends.find((l) => l.title === 'Base pairs');
    expect(arcs?.entries.kind === 'swatches' && arcs.entries.items.map((i) => i.label)).toEqual([
      'canonical',
      'wobble',
      'non-canonical',
    ]);
    // at fit a residue's mark is narrower than its code: none labeled; zoomed in, each code fits its mark
    const modified = () =>
      viewer!.frame!.layout.sections[0].tracks.find((t) => t.track.label === 'Modified')!.layers[0].lay as RangesLay;
    expect(Array.from({ length: modified().m }, (_, j) => modified().mlabel[j]).every((w) => !(w > 0))).toBe(true);
    viewer.view.setRange('/views/main/trna', 0, 14);
    await vi.advanceTimersByTimeAsync(50);
    viewer.draw();
    const lay = modified();
    const labelled = Array.from({ length: lay.m }, (_, j) =>
      lay.mlabel[j] > 0 ? lay.labels[lay.mitem[j]] : '',
    ).filter(Boolean);
    expect(labelled.slice(0, 3)).toEqual(['2MG', 'H2U', 'H2U']);
  });
});

describe('tRNA: the dot plot’s two triangles', () => {
  it('the covariation only above the diagonal, exactly the cells computed; the structure’s 27 pairs below', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    // its letters read from the file, their length its once open
    const trnaSystem = viewer.spec.model.value.systems.get('/data/trna');
    expect([trnaSystem?.source, trnaSystem?.letters, trnaSystem?.length]).toEqual(['/data/ehz', undefined, 76]);
    const section = viewer.frame!.layout.sections[1];
    expect(viewer.frame!.layout.issues).toEqual([]);
    const grid = section.lays[0][0].lay as GridLay;
    // 76 positions a side: one tile, a pixel a cell; a cell has a color where the covariation was computed
    expect(grid.images).toHaveLength(1);
    const im = grid.images[0];
    expect([im.width, im.height]).toEqual([76, 76]);
    const cells = new Set<string>();
    for (let y = 0; y < 76; y++)
      for (let x = 0; x < 76; x++) if (im.data[4 * (y * 76 + x) + 3] > 0) cells.add(`${y}-${x}`);
    // `upper`: y before x, and exactly the entries computed (their first position down, their second across)
    expect([...cells].every((c) => Number(c.split('-')[0]) < Number(c.split('-')[1]))).toBe(true);
    expect(cells).toEqual(new Set(DerivedData.covariation.map(([i, j]) => `${i}-${j}`)));
    const dots = section.lays[0][1].lay as PairsLay;
    // a dot a pair, the structure’s below the diagonal
    expect(dots.mode).toBe('dots');
    expect(dots.nd).toBe(27);
  });
});

describe('tRNA: the dot plot in both themes', () => {
  /** The WCAG contrast ratio of two colors given as [r, g, b] bytes. */
  const contrast = (a: ArrayLike<number>, b: ArrayLike<number>) => {
    const lum = (c: ArrayLike<number>) => {
      const ch = [c[0], c[1], c[2]].map((v) => {
        const x = v / 255;
        return x <= 0.03928 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4;
      });
      return 0.2126 * ch[0] + 0.7152 * ch[1] + 0.0722 * ch[2];
    };
    const [hi, lo] = [lum(a), lum(b)].sort((x, y) => y - x);
    return (hi + 0.05) / (lo + 0.05);
  };
  const rgb = (c: number) => [(c >>> 24) & 255, (c >>> 16) & 255, (c >>> 8) & 255];

  it('the lowest covariation is as quiet as the track under it in the light theme as in the dark (the ramp’s low end), the highest stands out', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const v = viewer;
    const entries = DerivedData.covariation as unknown as readonly (readonly number[])[];
    const lowest = entries.reduce((m, e) => (e[2] < m[2] ? e : m));
    const highest = entries.reduce((m, e) => (e[2] > m[2] ? e : m));
    const pixel = (e: readonly number[]) => {
      const im = (v.frame!.layout.sections[1].lays[0][0].lay as GridLay).images[0];
      const k = 4 * (e[0] * im.width + e[1]);
      return [im.data[k], im.data[k + 1], im.data[k + 2]];
    };
    for (const tokens of [LightTokens, DarkTokens]) {
      v.setTokens(tokens);
      await vi.advanceTimersByTimeAsync(1000);
      v.draw();
      const surface = rgb(tokens.colors['seq-track']);
      expect(contrast(pixel(lowest), surface), 'the lowest').toBeLessThan(1.6);
      expect(contrast(pixel(highest), surface), 'the highest').toBeGreaterThan(3);
    }
  });
});

describe('tRNA: the seed row and the seed', () => {
  it('the seed row: its letters and SS_cons pairs projected onto it, nested only', async () => {
    viewer = await openExample(trna('seed-row'), { gene: '' });
    const lay = arcsOf(viewer);
    expect(lay.n).toBe(21);
    // nothing crosses: the largest nested set of pairs is all of them
    expect(viewer.frame!.layout.issues).toEqual([]);
  });

  it('the seed draws in both themes, settled, with no issue', async () => {
    for (const t of ['light', 'dark'] as const) {
      const ctx = createSeqContext({
        registry: [
          { kind: 'format-runtime', runtime: SeqRuntime },
          { kind: 'format-runtime', runtime: AlignmentFileRuntime },
        ],
        ports: { assets: () => nodeFixtureAssets() },
      });
      const r = await ctx.renderSvg(trna('seed').build(), {
        size: { width: 720, height: 420 },
        tokens: t === 'dark' ? DarkTokens : LightTokens,
      });
      expect(r.issues).toEqual([]);
      expect(r.settled).toBe(true);
      ctx.retire();
    }
  });

  it('the seed viewer: the rows layer holds 954 rows over 118 columns; its legend lists A, C, G, U', async () => {
    viewer = await openExample(trna('seed'), { gene: '' });
    const lay = viewer.frame!.layout.sections[0].tracks.flatMap((t) => t.layers).find((l) => l.type.name === 'rows')!
      .lay as RowsLay;
    expect(lay.total).toBe(954);
    expect(lay.columns).toBe(118);
    const legend = viewLegends(viewer.frame!.layout).legends.find((l) => l.entries.kind === 'swatches')!;
    const items = legend.entries.kind === 'swatches' ? legend.entries.items.map((i) => i.label) : [];
    expect(items).toContain('U');
    expect(items).not.toContain('T');
  });
});

describe('tRNA: the editing views, each a patch undone', () => {
  const arcs = ['views', 'main', 'trna', 'pairs', 'arcs'];

  it('hide the tertiary pairs: a filter by that field in the arcs’ `data`, one patch; undo shows them again', async () => {
    viewer = await openExample(trna('arcs'), { gene: '', history: true });
    const v = viewer;
    expect(arcsOf(v).n).toBe(27);
    const o = await v.spec.patch('Hide the tertiary pairs', {
      ops: [{ op: 'set', at: arcs, key: 'data', value: without('/data/pairs', 'tertiary', 'yes') }],
    });
    expect(o.ok).toBe(true);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    expect(keysOf(arcsOf(v))).toHaveLength(21);
    for (const k of ['8-14', '15-48', '18-55', '19-56', '26-44', '54-58']) expect(keysOf(arcsOf(v))).not.toContain(k);
    expect(keysOf(arcsOf(v))).toContain('32-38');
    expect((await (v.spec.history as UiHistory).stack.undo()).ok).toBe(true);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    expect(arcsOf(v).n).toBe(27);
  });

  it('change the covariation’s palette: a theme patch on the grid layer, one step; undo restores it', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '', history: true });
    const v = viewer;
    const layer = ['views', 'main', 'dots', 'plot', 'covariation'];
    const param = () => {
      let n = v.spec.root.value;
      for (const id of layer) n = n.children!.find((c) => c.id === id)!;
      return JSON.stringify(n.params?.color);
    };
    const before = param();
    // the cells' colors as drawn: the tile's opaque pixels
    const fills = () => {
      const im = (v.frame!.layout.sections[1].lays[0][0].lay as GridLay).images[0];
      const out: number[] = [];
      for (let k = 0; k < im.data.length && out.length < 40; k += 4)
        if (im.data[k + 3] > 0) out.push(im.data[k], im.data[k + 1], im.data[k + 2]);
      return out;
    };
    const drawn = fills();
    const key = JSON.stringify(viewLegends(v.frame!.layout).legends.map((l) => l.entries));
    expect(drawn.length).toBeGreaterThan(0);
    const o = await v.spec.patch('Recolor the covariation', {
      ops: [
        {
          op: 'set',
          at: layer,
          key: 'color',
          value: theme.field('value', { domain: [0, 1.2], palette: theme.continuous({ colors: 'blues' }) }),
        },
      ],
    });
    expect(o.ok).toBe(true);
    expect(param()).not.toBe(before);
    // it redraws: the same cells in other colors, the legend’s ramp too
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    expect(fills()).not.toEqual(drawn);
    expect(JSON.stringify(viewLegends(v.frame!.layout).legends.map((l) => l.entries))).not.toBe(key);
    expect((await (v.spec.history as UiHistory).stack.undo()).ok).toBe(true);
    expect(param()).toBe(before);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    expect(fills()).toEqual(drawn);
  });
});
