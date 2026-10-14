/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The render prop `unobservedGray` (seqspec:viewer#the-render-backend): a unit's unobserved positions drawn in the skin's
// `unobserved` token, in plain text, in tinted text (`wrappedColors`) and in the letters' cells, zoomed in and zoomed
// out (the cells' classes keep them apart from the letters, and the prop is part of the classes' identity), opaque.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { SimpleParamsSchema } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { SkinToken } from '../src/draw/backend.js';
import { KernelRenderers, type LayerRenderer, type RenderProps, themeProvider } from '../src/index.js';
import { columns, drawn, entryOf, record, rects, texts } from './draw-fixtures.js';
import { opened } from './layout-fixtures.js';

const GRAY = SkinToken.unobserved;
const INK = SkinToken.ink;

/** A unit of `n` residues of `value` repeated, observed at [lo, hi) only: the rest is in the sequence, not the model. */
function doc(o: { layout?: 'text'; n?: number; value?: string; lo?: number; hi?: number; theme?: string } = {}) {
  const n = o.n ?? 40;
  const value = (o.value ?? 'MKVLATGHEK').repeat(Math.ceil(n / (o.value ?? 'MKVLATGHEK').length)).slice(0, n);
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'chain', value });
  b.data.unit({ ref: 'A', sequence: 'chain', observed: [[o.lo ?? 10, o.hi ?? 30]] } as never);
  const v = b.view({ ref: 'v', label: 'V', ...(o.layout ? { layout: o.layout } : {}) });
  v.block({ unit: 'A' } as never);
  const track = v.track({ label: 'Letters' });
  if (o.theme === 'none') track.letters();
  else if (o.theme)
    track.letters({ color: { theme: o.theme, fallback: { theme: 'uniform', params: { color: '#336699' } } } } as never);
  else track.letters({ color: { theme: 'letter', params: { colors: 'Hydrophobicity' } } } as never);
  return b.getDocument();
}
const TEXT: Partial<RenderProps> = { wrappedGroupLabels: 'none' };

describe('unobserved positions (unobservedGray)', () => {
  it("are the skin's unobserved token in tinted text (wrappedColors), the observed in their own colors", async () => {
    const t = await drawn(doc({ layout: 'text' }), {
      width: 900,
      height: 200,
      render: { ...TEXT, wrappedColors: true, unobservedGray: true },
    });
    const gray = texts(record(t.ctx, 'Letters').rec).map((d) => d.color === GRAY);
    expect(gray).toHaveLength(40);
    expect(gray.slice(0, 10).every(Boolean)).toBe(true);
    expect(gray.slice(10, 30).some(Boolean)).toBe(false);
    expect(gray.slice(30).every(Boolean)).toBe(true);
  });

  it('are the token in plain text too (no wrappedColors), the observed letters in the ink', async () => {
    const t = await drawn(doc({ layout: 'text' }), {
      width: 900,
      height: 200,
      render: { ...TEXT, unobservedGray: true },
    });
    const colors = texts(record(t.ctx, 'Letters').rec).map((d) => d.color);
    expect(colors).toHaveLength(40);
    expect(colors.slice(0, 10).every((c) => c === GRAY)).toBe(true);
    expect(colors.slice(10, 30).every((c) => c === INK)).toBe(true);
    expect(colors.slice(30).every((c) => c === GRAY)).toBe(true);
  });

  it('are as any other position without the prop, in text and in cells', async () => {
    for (const layout of ['text', undefined] as const) {
      const t = await drawn(doc({ layout }), { width: 900, height: 200, render: { ...TEXT, wrappedColors: true } });
      const { rec } = record(t.ctx, 'Letters');
      expect(texts(rec).some((d) => d.color === GRAY)).toBe(false);
      expect(rects(rec).some((r) => r[4] === GRAY)).toBe(false);
    }
  });

  it('are gray cells zoomed in (linear), opaque, their letters in the ink; the observed cells in the theme', async () => {
    const t = await drawn(doc(), { width: 900, height: 200, render: { unobservedGray: true } });
    const { rec } = record(t.ctx, 'Letters');
    const cells = rects(rec);
    expect(cells).toHaveLength(40);
    const gray = cells.map((r) => r[4] === GRAY);
    expect(gray.slice(0, 10).every(Boolean)).toBe(true);
    expect(gray.slice(10, 30).some(Boolean)).toBe(false);
    expect(gray.slice(30).every(Boolean)).toBe(true);
    expect(cells.filter((r) => r[4] === GRAY).every((r) => r[5] === 255)).toBe(true);
    const letters = texts(rec);
    expect(letters.slice(0, 10).every((d) => d.color === INK)).toBe(true);
  });

  it('are their own class zoomed out: a run of one letter across the observed boundary is not one color', async () => {
    // 4,000 residues of one letter in 900 px (a few positions a device column): the classes are keyed by letter, and
    // without a key of their own the unobserved run would take the first color for the whole block.
    const d = doc({ n: 4000, value: 'A', lo: 1000, hi: 3000 });
    const t = await drawn(d, { width: 900, height: 200, render: { unobservedGray: true } });
    const e = entryOf(t.ctx, 'Letters');
    const y = Math.floor(e.top + e.height / 2);
    const { rec } = record(t.ctx, 'Letters');
    const l = t.layout('v.main');
    const [x0, x1] = [Math.ceil(l.dataLeft), Math.floor(l.dataRight)];
    const cols = columns(rec, x0, x1, y);
    const w = x1 - x0;
    // The first and last fifth are unobserved, the middle observed.
    expect(cols.slice(2, Math.floor(w / 5)).every((c) => c === GRAY)).toBe(true);
    expect(cols.slice(Math.ceil((2 * w) / 5), Math.floor((3 * w) / 5)).every((c) => c >= 0)).toBe(true);
    expect(cols.slice(Math.ceil((4 * w) / 5), w - 2).every((c) => c === GRAY)).toBe(true);
  });

  it('rebuild the zoomed-out classes when the prop changes (it is part of their identity)', async () => {
    const d = doc({ n: 4000, value: 'A', lo: 1000, hi: 3000 });
    const t = await drawn(d, { width: 900, height: 200, render: { unobservedGray: false } });
    const e = entryOf(t.ctx, 'Letters');
    const y = Math.floor(e.top + e.height / 2);
    const l = t.layout('v.main');
    const at = () => columns(record(t.ctx, 'Letters').rec, Math.ceil(l.dataLeft), Math.ceil(l.dataLeft) + 20, y);
    expect(at().some((c) => c === GRAY)).toBe(false);
    t.ctx.setRenderProps({ unobservedGray: true });
    await t.ctx.graph.update();
    expect(at().every((c) => c === GRAY)).toBe(true);
    t.ctx.setRenderProps({ unobservedGray: false });
    await t.ctx.graph.update();
    expect(at().some((c) => c === GRAY)).toBe(false);
  });

  it("are opaque under a theme with an alpha: an unobserved cell never reads the previous position's", async () => {
    // A theme whose alpha depends on the position: 100 below 20, 200 from it.
    const alphaTheme = themeProvider('test:alpha', ['color'], 0, () => ({
      value: () => 0x336699,
      alpha: (loc) => (loc.index < 20 ? 100 : 200),
    }));
    const t = await opened(doc({ theme: 'test:alpha' }), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      registry: [
        {
          svs: { themes: [{ name: 'test:alpha', params: { color: SimpleParamsSchema({}) }, fallback: null }] },
          themes: [alphaTheme],
        },
      ],
      render: { unobservedGray: true },
    });
    t.ctx.viewports.get('v.main')?.setSize(900, 200, 1);
    await t.ctx.graph.update();
    const cells = rects(record(t.ctx, 'Letters').rec);
    expect(cells).toHaveLength(40);
    expect(cells.map((r) => r[5])).toEqual([
      ...new Array(10).fill(255),
      ...new Array(10).fill(100),
      ...new Array(10).fill(200),
      ...new Array(10).fill(255),
    ]);
  });
});
