/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The section layout interface's contract with followers and layers: a follower reads its leader's viewport as the
// leader's layout reads it (a window over a cyclic origin wraps); a `section` layer is laid out once a section, with the
// placement, over the page; a layout's label sink reaches a layer; a band's own scroll and row zoom reach its layers
// and the key compares them by value.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  bundledTextMetrics,
  type LabelSink,
  type LayerType,
  type LayoutContext,
  LightTokens,
  layoutView,
  type SectionLayoutType,
  type StripPlacement,
  sameWindow,
  TracksLayout,
} from '../index.js';
import { BuiltInLayerTypes } from '../built-ins.js';
import { P53 } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 640, height: 300 };

/** A layout that places as tracks but over the whole axis, its viewport a window followers read (as a ring's). */
const Whole: SectionLayoutType<StripPlacement> = {
  ...TracksLayout,
  name: 'whole',
  followWindow: (_s, a, v) => v ?? { start: 0, end: a.total },
};

function cyclicSpec() {
  const b = createSeqSpec({ title: 'cyclic' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true });
  b.views.view('v', (v) => {
    const lead = v.section('ring', { of: c });
    lead.track('l').layer('l', { type: 'letters' });
    v.section('detail', { of: c, follow: lead }).track('l').layer('l', { type: 'letters' });
  });
  return b.build();
}

describe('the layout contract', () => {
  it('a follower follows its leader’s viewport as the leader’s layout reads it: a window over the origin wraps', async () => {
    const s = await settle(cyclicSpec());
    const layout = layoutView({
      model: s.model,
      view: s.view,
      size: Size,
      viewports: { '/views/v/ring': { start: 90, end: 110 } },
      tokens: LightTokens,
      text: bundledTextMetrics(),
      data: s.store,
      layouts: { tracks: Whole },
    });
    const [ring, detail] = layout.sections;
    // the leader places its whole axis; the follower its window, 20 positions over the origin
    expect(ring.axis!.u1 - ring.axis!.u0).toBe(100);
    expect(detail.axis!.total).toBe(20);
    expect(detail.axis!.pieces).toEqual([90, 100, 0, 10]);
  });

  it('a section layer is laid out once a section with the placement, over the section’s box on the page', async () => {
    const b = createSeqSpec({ title: 'section scope' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
    b.views.view('v', (v) => {
      const s = v.section('a', { of: p });
      s.track('l').layer('l', { type: 'letters' });
      s.track('pairs').layer('pairs', { type: 'ruler' });
    });
    const s = await settle(b.build());
    const seen: LayoutContext[] = [];
    const pairs: LayerType<null, null> = {
      name: 'pairs',
      scope: 'section',
      prepare: () => null,
      layout(_p, ctx) {
        seen.push(ctx);
        return null;
      },
      draw(_l, pass, out) {
        if (pass === 'data') out.fill(1, 2, 3, 4, 0xff0000ff);
      },
      hit: () => undefined,
    };
    const layout = layoutView({
      model: s.model,
      view: s.view,
      size: Size,
      tokens: LightTokens,
      text: bundledTextMetrics(),
      data: s.store,
      types: { ...BuiltInLayerTypes, ruler: pairs as LayerType<any, any> },
    });
    const sec = layout.sections[0];
    expect(seen.length).toBe(1);
    expect(seen[0].place).toBe(sec.placement);
    expect(seen[0].box).toEqual(sec.placement.box);
    expect(sec.sectionLays.length).toBe(1);
    expect(sec.lays.flat().some((l) => l.type === pairs)).toBe(false);
  });

  it('a layout’s label sink reaches its layers; a band’s scroll and row zoom reach its layers', async () => {
    const b = createSeqSpec({ title: 'labels' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
    b.views.view('v', (v) => v.section('a', { of: p }).track('l').layer('l', { type: 'letters' }));
    const s = await settle(b.build());
    const labels: string[] = [];
    const sink: LabelSink = { label: (text) => labels.push(text) };
    const seen: LayoutContext[] = [];
    const probe: LayerType<null, null> = {
      name: 'probe',
      prepare: () => null,
      layout(_p, ctx) {
        seen.push(ctx);
        if (ctx.labels && ctx.labels !== 'inside') ctx.labels.label('TAD', 0, 61, { owner: 'o', key: 'k' });
        return null;
      },
      draw: () => {},
      hit: () => undefined,
    };
    const Sinking: SectionLayoutType<StripPlacement> = {
      ...TracksLayout,
      place(section, ctx, viewport, prev) {
        return TracksLayout.place(section, { ...ctx, layStrip: (pl, i) => ctx.layStrip(pl, i, sink) }, viewport, prev);
      },
    };
    layoutView({
      model: s.model,
      view: s.view,
      size: Size,
      viewport: { start: 0, end: 393, bands: { '/views/v/a/l': { bandTop: 3, rowScale: 12 } } },
      tokens: LightTokens,
      text: bundledTextMetrics(),
      data: s.store,
      types: { ...BuiltInLayerTypes, letters: probe as LayerType<any, any> },
      layouts: { tracks: Sinking },
    });
    expect(labels).toEqual(['TAD']);
    expect(seen[0].bandTop).toBe(3);
    expect(seen[0].rowScale).toBe(12);
  });

  it('the key compares a viewport’s bands by value', () => {
    const a = { start: 0, end: 10, bands: { t: { bandTop: 2, rowScale: 8 } } };
    expect(sameWindow(a, { start: 0, end: 10, bands: { t: { bandTop: 2, rowScale: 8 } } })).toBe(true);
    expect(sameWindow(a, { start: 0, end: 10, bands: { t: { bandTop: 3, rowScale: 8 } } })).toBe(false);
    expect(sameWindow(a, { start: 0, end: 10 })).toBe(false);
  });
});
