/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The values layer's reference lines (_spec/README.md "Render", Layer types, `values`; EC-217): a hairline at a value
// through the layer's own scale, over the data and under the hover, never hit; outside a closed domain not drawn and
// said once, inside an open one part of it; its label at the well's right end; the accessible words.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { contrastRgba } from '@molstar/seqx-core/theme/index';
import { describe, expect, it } from 'vitest';
import { DarkTokens, hitView, LightTokens, type ValuesLay } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { ValuesType } from '../values.js';
import { holdsReference, loadReference, ReferenceHooks } from '../values-reference-hooks.js';
import { BigRuntime } from './fixtures.js';
import { settle } from './helpers.js';
import { layerLay, readPaths, readRects } from './small-helpers.js';

const Size = { width: 400, height: 80 };
const LineInk = LightTokens.colors['msui-text'];

type Line = { at: number; label?: string; color?: string; style?: string };

function spec(
  values: readonly number[],
  layer: Record<string, unknown> = {},
  log?: readonly [number | null, number | null],
  track: Record<string, unknown> = {},
) {
  const b = createSeqSpec({ title: 'reference' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'A'.repeat(values.length) });
  const v = b.data.values('v', { target: s, values: Float32Array.from(values) });
  const scale = log === undefined ? undefined : b.views.scale('log', { type: 'log', domain: log } as never);
  b.views.view('main', (view) =>
    view
      .section('s', { of: s })
      .track('v', track as never)
      .layer('values', { type: 'values', data: v, ...(scale ? { scale } : {}), ...layer } as never),
  );
  return b.build();
}

async function frame(
  values: readonly number[],
  layer: Record<string, unknown> = {},
  opts: { height?: number; log?: readonly [number | null, number | null]; width?: number } = {},
) {
  const s = await settle(spec(values, layer, opts.log), [], { heights: { values: opts.height ?? 60 } } as never);
  const f = await s.frame({ size: { width: opts.width ?? Size.width, height: Size.height }, settled: true });
  const lay = layerLay<ValuesLay>(f, 0);
  return { s, f, lay };
}

/** The line a layer draws at `i`: its solid rect's middle, else its dashed path's y. */
const lineYs = (f: Awaited<ReturnType<typeof frame>>['f']) => [
  ...readRects(f.lists)
    .filter((r) => r.fill === LineInk && r.w > 100)
    .map((r) => r.y + r.h / 2),
  ...readPaths(f.lists)
    .filter((p) => !p.closed && p.stroke === LineInk && p.dash > 0)
    .map((p) => p.xy[1]),
];

/** Where a value sits in the well: the layer's box, its domain and its headroom (`Fill`). */
const heightOf = (lay: ValuesLay, v: number) => {
  const [lo, hi] = lay.domain!;
  const t = lay.yLog ? (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo)) : (v - lo) / (hi - lo);
  return lay.box.y + lay.box.height - t * 0.9 * lay.box.height;
};

const Positions = Array.from({ length: 10 }, (_, i) => i + 1);

describe('a reference line is where a value is', () => {
  it('on bars, a line and points alike: the bar’s top and the line at its value coincide', async () => {
    for (const encoding of ['bars', 'line', 'points']) {
      const { f, lay } = await frame(Positions, { encoding, domain: [0, 10], reference: [{ at: 7, style: 'solid' }] });
      const [y] = lineYs(f);
      // a hairline snaps to a device pixel: within half of one
      expect(Math.abs(y - heightOf(lay, 7)), encoding).toBeLessThanOrEqual(0.5 + 1e-6);
    }
    const { f } = await frame(Positions, { domain: [0, 10], reference: [{ at: 7, style: 'solid' }] });
    const bar = readRects(f.lists).find(
      (r) =>
        (r.fill & 255) === 0x99 &&
        Math.abs(r.h - (heightOf(layerLay<ValuesLay>(f, 0), 0) - heightOf(layerLay<ValuesLay>(f, 0), 7))) < 0.6,
    )!;
    expect(bar).toBeDefined();
  });

  it('a line at the domain’s top sits where a value there sits, the headroom above it', async () => {
    const { f, lay } = await frame(Positions, { domain: [0, 10], reference: [{ at: 10, style: 'solid' }] });
    const [y] = lineYs(f);
    expect(Math.abs(y - (lay.box.y + 0.1 * lay.box.height))).toBeLessThanOrEqual(0.5 + 1e-6);
  });

  it('on a log scale, at its value: halfway up the span between the decades 1 and 10^4 for 100', async () => {
    const { f, lay } = await frame(
      [1, 10, 100, 1000, 10000],
      { reference: [{ at: 100, style: 'solid' }] },
      { log: [1, 10000] },
    );
    const [y] = lineYs(f);
    expect(Math.abs(y - heightOf(lay, 100))).toBeLessThanOrEqual(0.5 + 1e-6);
    expect(Math.abs(y - (lay.box.y + lay.box.height - 0.5 * 0.9 * lay.box.height))).toBeLessThanOrEqual(0.5 + 1e-6);
  });

  it('on a log scale a line at or under zero is not on the scale: none, and no said outside (the layer’s own info says why)', async () => {
    for (const at of [0, -3]) {
      const { f, lay } = await frame([1, 10, 100], { reference: [{ at, style: 'solid' }] }, { log: [1, 100] });
      expect(lineYs(f)).toEqual([]);
      expect(lay.ref.shown[0]).toBe(0);
      expect(f.layout.issues.filter((i) => i.code === 'seq.render.reference.outside')).toEqual([]);
    }
  });

  it('two lines: two hairlines, each its own', async () => {
    const { f, lay } = await frame(Positions, {
      domain: [0, 10],
      reference: [
        { at: 7, style: 'solid' },
        { at: 3, style: 'solid' },
      ],
    });
    const ys = lineYs(f).sort((a, b) => a - b);
    expect(ys).toHaveLength(2);
    expect(Math.abs(ys[0] - heightOf(lay, 7))).toBeLessThanOrEqual(0.5 + 1e-6);
    expect(Math.abs(ys[1] - heightOf(lay, 3))).toBeLessThanOrEqual(0.5 + 1e-6);
  });
});

describe('a line outside the domain, and an open domain (EC-217)', () => {
  it('a closed domain’s outside line is not drawn, said once with its index and value', async () => {
    const { f, lay } = await frame(Positions, {
      domain: [0, 6],
      reference: [
        { at: 3, style: 'solid' },
        { at: 7.3, style: 'solid' },
        { at: 9, style: 'solid' },
        { at: -1, style: 'solid' },
      ],
    });
    expect(lineYs(f)).toHaveLength(1);
    expect([...lay.ref.shown]).toEqual([1, 0, 0, 0]);
    const said = f.layout.issues.filter((i) => i.code === 'seq.render.reference.outside');
    expect(said).toHaveLength(1);
    expect(said[0].severity).toBe('info');
    expect(said[0].message).toContain('1');
    expect(said[0].message).toContain('7.3');
  });

  it('an open domain includes the lines: data to 6 and a line at 7.3, the domain [0, 7.3] and the line drawn', async () => {
    const { f, lay } = await frame([1, 2, 3, 4, 5, 6], { reference: [{ at: 7.3, style: 'solid' }] });
    expect(lay.domain).toEqual([0, 7.3]);
    expect(lineYs(f)).toHaveLength(1);
    expect(f.layout.issues.filter((i) => i.code === 'seq.render.reference.outside')).toEqual([]);
    // a line under the data is part of an open low end too
    const low = await frame([1, 2, 3], { reference: [{ at: -2, style: 'solid' }] });
    expect(low.lay.domain).toEqual([-2, 3]);
    expect(lineYs(low.f)).toHaveLength(1);
  });

  it('a plot with no value over the line shows the line over the data’s top: the absence is seen', async () => {
    const { f, lay } = await frame([1, 2, 3], { reference: [{ at: 7.3, style: 'solid' }] });
    const top = readRects(f.lists)
      .filter((r) => (r.fill & 255) === 0x99)
      .map((r) => r.y);
    expect(Math.min(...top)).toBeGreaterThan(lineYs(f)[0]);
    expect(lay.domain![1]).toBe(7.3);
  });

  it('a log scale’s open end takes the lines too; a closed end of a scale node does not move', async () => {
    const open = await frame([1, 10], { reference: [{ at: 1000, style: 'solid' }] }, { log: [null, null] });
    expect(open.lay.domain![1]).toBe(1000);
    expect(lineYs(open.f)).toHaveLength(1);
    const closed = await frame([1, 10], { reference: [{ at: 1000, style: 'solid' }] }, { log: [1, 100] });
    expect(lineYs(closed.f)).toEqual([]);
    expect(closed.f.layout.issues.filter((i) => i.code === 'seq.render.reference.outside')).toHaveLength(1);
  });
});

describe('the label', () => {
  const Many: Line[] = [{ at: 7.3, label: 'genome-wide' }];

  it('sits at the well’s right end, above the line, 10 px, over a patch of the well’s ground', async () => {
    const { f, lay } = await frame(Positions, { domain: [0, 10], reference: Many });
    const r = lay.ref;
    expect(r.labelShown[0]).toBe(1);
    const [x, y, w, h] = [...r.patch.slice(0, 4)];
    expect(x + w).toBeCloseTo(lay.box.x + lay.box.width, 3);
    expect(y + h).toBeCloseTo(r.top[0], 3);
    expect(h).toBe(12);
    expect(r.font).toBe('10px sans');
    const patch = readRects(f.lists, 'data').find((q) => Math.abs(q.x - x) < 0.01 && Math.abs(q.y - y) < 0.01);
    expect(patch?.fill).toBe(LightTokens.colors['seq-track']);
    expect(r.text[0]).toBe('genome-wide · 7.3');
  });

  it('goes below the line when the line is within a label height of the well’s top', async () => {
    const { lay } = await frame(Positions, { domain: [0, 10], reference: [{ at: 10, label: 'top' }] });
    const r = lay.ref;
    expect(r.labelShown[0]).toBe(1);
    // under the hairline's bottom edge, never over its lower half
    expect(r.patch[1]).toBeCloseTo(r.top[0] + r.width, 3);
    expect(r.patch[1]).toBeGreaterThanOrEqual(r.y[0] + r.width / 2 - 1e-6);
  });

  it('is not drawn in a well under 40 px, the line is', async () => {
    const short = await frame(Positions, { domain: [0, 10], reference: Many }, { height: 30 });
    expect(short.lay.box.height).toBeLessThan(40);
    expect(short.lay.ref.labelShown[0]).toBe(0);
    expect(short.lay.ref.shown[0]).toBe(1);
    expect(readRects(short.f.lists, 'data').filter((q) => q.fill === LightTokens.colors['seq-track']).length).toBe(0);
    const tall = await frame(Positions, { domain: [0, 10], reference: Many }, { height: 44 });
    expect(tall.lay.box.height).toBeGreaterThanOrEqual(40);
    expect(tall.lay.ref.labelShown[0]).toBe(1);
  });

  it('of two overlapping labels the later is dropped, its line stays', async () => {
    const { lay } = await frame(Positions, {
      domain: [0, 10],
      reference: [
        { at: 7.3, label: 'genome-wide' },
        { at: 7.0, label: 'suggestive' },
      ],
    });
    expect([...lay.ref.shown.slice(0, 2)]).toEqual([1, 1]);
    expect([...lay.ref.labelShown.slice(0, 2)]).toEqual([1, 0]);
    // far apart: both
    const apart = await frame(Positions, {
      domain: [0, 10],
      reference: [
        { at: 9, label: 'genome-wide' },
        { at: 3, label: 'suggestive' },
      ],
    });
    expect([...apart.lay.ref.labelShown.slice(0, 2)]).toEqual([1, 1]);
  });

  it('a line without a label draws its value alone, the value as the legends write one', async () => {
    const plain = await frame(Positions, { domain: [0, 10], reference: [{ at: 5 }, { at: 8.123456 }] });
    expect([...plain.lay.ref.labelShown.slice(0, 2)]).toEqual([1, 1]);
    expect(plain.lay.ref.text.slice(0, 2)).toEqual(['5', '8.123']);
  });

  it('the SVG export draws the same text', async () => {
    const r = await renderSvg(spec(Positions, { domain: [0, 10], reference: Many }, undefined, { height: 60 }), {
      size: Size,
    });
    expect(r.svg).toContain('genome-wide · 7.3');
  });

  it('a long label is cut to the well’s right third with an ellipsis', async () => {
    const cut = await frame(
      Positions,
      { domain: [0, 10], reference: [{ at: 5, label: 'a very long threshold label that cannot fit' }] },
      { width: 240 },
    );
    const r = cut.lay.ref;
    expect(r.text[0].endsWith('…')).toBe(true);
    expect(r.patch[2]).toBeLessThanOrEqual(240 / 3 + 1e-6);
  });
});

describe('the label and the hover', () => {
  it('the label is drawn with the data pass, before the hover and the cursor: a hovered dot under it stays seen', async () => {
    const { s, f, lay } = await frame(Positions, {
      domain: [0, 10],
      encoding: 'points',
      reference: [{ at: 9, label: 'genome-wide' }],
    });
    const [x, y, w, h] = [...lay.ref.patch.slice(0, 4)];
    expect(readRects(f.lists, 'labels')).toEqual([]);
    expect(readRects(f.lists, 'data').some((q) => Math.abs(q.x - x) < 0.01 && Math.abs(q.y - y) < 0.01)).toBe(true);
    // the dot at the right end, under the label's patch: hovering it draws its ring in the hover pass, after the patch
    const dotX = lay.x1[9] - 1;
    const hit = hitView(f.layout, dotX, y + h / 2);
    expect(hit).toBeDefined();
    expect(dotX).toBeGreaterThan(x);
    expect(w).toBeGreaterThan(0);
    expect(s).toBeDefined();
  });
});

describe('marks, export, hit', () => {
  it('a solid line is a rect, a dashed and a dotted one a dashed path, over the data (after it in the pass)', async () => {
    const { f } = await frame(Positions, {
      domain: [0, 10],
      reference: [{ at: 2, style: 'solid' }, { at: 5 }, { at: 8, style: 'dotted' }],
    });
    const rects = readRects(f.lists).filter((r) => r.fill === LineInk && r.w > 100);
    expect(rects).toHaveLength(1);
    expect(rects[0].h).toBeCloseTo(1, 5);
    const dashed = readPaths(f.lists).filter((p) => p.stroke === LineInk && p.dash > 0);
    expect(dashed.map((p) => p.dash)).toEqual([4, 1.5]);
    expect(dashed.every((p) => p.strokeWidth === 1 && p.xy.length === 4)).toBe(true);
    // after the bars: the last rect of the data pass before the line is a bar
    const all = readRects(f.lists);
    expect(all.findIndex((r) => r.fill === LineInk && r.w > 100)).toBeGreaterThan(
      all.findIndex((r) => (r.fill & 255) === 0x99),
    );
  });

  it('the SVG export draws the dashed line with its dash array, the solid as a rect', async () => {
    const r = await renderSvg(
      spec(Positions, { domain: [0, 10], reference: [{ at: 5 }, { at: 8, style: 'dotted' }] }),
      { size: Size },
    );
    expect(r.svg).toMatch(/stroke-dasharray="4"/);
    expect(r.svg).toMatch(/stroke-dasharray="1.5"/);
    const solid = await renderSvg(spec(Positions, { domain: [0, 10], reference: [{ at: 5, style: 'solid' }] }), {
      size: Size,
    });
    expect(solid.svg).not.toMatch(/stroke-dasharray/);
  });

  it('is never hit: a point on the line hits the layer as the data does', async () => {
    const { f, lay } = await frame(Positions, { domain: [0, 10], reference: [{ at: 5, label: 'mid' }] });
    const at = hitView(f.layout, lay.box.x + 200, lay.ref.y[0]);
    expect(at?.layer).toBeDefined();
    const bare = await frame(Positions, { domain: [0, 10] });
    expect(hitView(bare.f.layout, bare.lay.box.x + 200, bare.lay.box.y + 10)?.layer).toBe(at?.layer);
  });

  it('reads aloud: "reference line at 7.3 (genome-wide)" a drawn line, none for an outside one', async () => {
    const { lay } = await frame(Positions, {
      domain: [0, 8],
      reference: [{ at: 7.3, label: 'genome-wide' }, { at: 5 }, { at: 9, label: 'far' }],
    });
    expect(ValuesType.summary!(lay)).toBe('reference line at 7.3 (genome-wide), reference line at 5');
    const none = await frame(Positions, { domain: [0, 8] });
    expect(ValuesType.summary!(none.lay)).toBeUndefined();
  });
});

describe('light and dark', () => {
  it('the default line (`msui-text`) reads against the well’s ground at 3:1 in both themes, and against `seq-value` in light', () => {
    for (const t of [LightTokens, DarkTokens])
      expect(contrastRgba(t.colors['msui-text'], t.colors['seq-track'])).toBeGreaterThanOrEqual(3);
    expect(contrastRgba(LightTokens.colors['msui-text'], LightTokens.colors['seq-value'])).toBeGreaterThanOrEqual(3);
    // dark: no ink is 3:1 against both a near-black ground and the light series blue (1.75:1 for `msui-text`)
  });

  it('the line is drawn in `msui-text`, the label’s ink stays muted', async () => {
    const { f, lay } = await frame(Positions, { domain: [0, 10], reference: [{ at: 7, label: 'x' }] });
    expect(lineYs(f)).toHaveLength(1);
    expect(lay.ref.ink).toBe(LightTokens.colors['msui-text-muted']);
  });
});

describe('a shared scale', () => {
  it('takes the lines of a layer with no value in view: two tracks on one open scale agree', async () => {
    const b = createSeqSpec({ title: 'shared' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const v = b.data.values('signal', { target: chr, source: src });
    const scale = b.views.scale('shared', { type: 'linear' } as never);
    b.views.view('main', (view) => {
      const s = view.section('chr', { of: chr });
      s.track('a').layer('a', { type: 'values', data: v, scale } as never);
      s.track('b').layer('b', { type: 'values', data: v, scale, reference: [{ at: 500 }] } as never);
    });
    const s = await settle(b.build(), [BigRuntime]);
    // the first frames: nothing read yet (no value in view), the lines the only extent
    const input = { size: { width: 400, height: 120 }, viewport: { start: 0, end: 50 }, settled: true } as const;
    const f = s.renderer.frame(input);
    const [a, bLay] = [layerLay<ValuesLay>(f, 0), layerLay<ValuesLay>(f, 1)];
    expect(a.n).toBe(0);
    expect(bLay.domain![1]).toBeGreaterThanOrEqual(500);
    expect(a.domain![1]).toBe(bLay.domain![1]);
  });
});

describe('the lines’ chunk', () => {
  it('is wanted by a section with a layer that has lines, not by one without; a model of one waits for it', async () => {
    const section = (node: ReturnType<typeof spec>) =>
      node.children!.find((c) => c.kind === 'views')!.children![0].children![0];
    expect(holdsReference(section(spec(Positions, { reference: [{ at: 1 }] })))).toBe(true);
    expect(holdsReference(section(spec(Positions, { reference: [] })))).toBe(false);
    expect(holdsReference(section(spec(Positions)))).toBe(false);
    await frame(Positions, { reference: [{ at: 1 }] });
    // the model's open awaited it: the lines drawn at the first frame
    expect(ReferenceHooks.chunk!.drawReference).toBeTypeOf('function');
    await loadReference();
    expect(ReferenceHooks.chunk!.layoutReference).toBeTypeOf('function');
  });
});
