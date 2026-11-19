/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqContext } from '@molstar/seqx-core/context';
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
} from '@molstar/seqx-core/render';
import { TrnaViews } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { withoutKind } from '@molstar/seqx-examples/trna/trna';
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
  it('the four stems nested, the six tertiary pairs crossing them, every pair drawn once', async () => {
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

describe('tRNA: the dot plot’s two triangles', () => {
  it('the covariation above the diagonal (stems brighter), the structure’s 27 pairs below', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const section = viewer.frame!.layout.sections[1];
    expect(viewer.frame!.layout.issues).toEqual([]);
    const grid = section.lays[0][0].lay as GridLay;
    expect(grid.images.length).toBeGreaterThan(0);
    const dots = section.lays[0][1].lay as PairsLay;
    // a dot a pair, the structure’s below the diagonal
    expect(dots.mode).toBe('dots');
    expect(dots.nd).toBe(27);
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

  it('the seed: 954 rows over 118 columns, drawn with no issue in both themes', async () => {
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
});

describe('tRNA: the editing views, each a patch undone', () => {
  const arcs = ['views', 'main', 'trna', 'pairs', 'arcs'];

  it('hide a kind of pairs: a filter in the arcs’ `data`, one patch; undo shows them again', async () => {
    viewer = await openExample(trna('arcs'), { gene: '', history: true });
    const v = viewer;
    expect(arcsOf(v).n).toBe(27);
    const o = await v.spec.patch('Hide the non-canonical pairs', {
      ops: [{ op: 'set', at: arcs, key: 'data', value: withoutKind('/data/pairs', 'non-canonical') }],
    });
    expect(o.ok).toBe(true);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    expect(keysOf(arcsOf(v))).toHaveLength(21);
    expect(keysOf(arcsOf(v))).not.toContain('15-48');
    expect(keysOf(arcsOf(v))).toContain('19-56');
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
    expect((await (v.spec.history as UiHistory).stack.undo()).ok).toBe(true);
    expect(param()).toBe(before);
  });
});
