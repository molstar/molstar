/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section placed nowhere that says why (an issue at its address: its `of` names no system) is one row naming its
// issue, never a blank well per track (_spec/README.md "Render", Tracks layout, Empty tracks; the empty-well rule); a
// section still waiting for its source says nothing (pending). The canvas's draw lists and the SVG say the same.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { BuiltInLayerTypes, BuiltInLayouts } from '../built-ins.js';
import {
  bundledTextMetrics,
  DefaultRenderConfig,
  LightTokens,
  recordDrawLists,
  Renderer,
  renderConfig,
} from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 400 };
const Issue = 'The section shows no system of the spec';

function spec(): Node {
  const b = createSeqSpec({ title: 'section issue' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(50) });
  const here = b.data.ranges('here', { target: s, rows: { start: [1], end: [10], id: ['a'] } });
  b.views.view('main', (v) => {
    v.section('ok', { of: s }).track('t', { label: 'Here' }).layer('l', { type: 'ranges', data: here });
    const lost = v.section('lost', { of: '/data/nothing' });
    lost.track('a', { label: 'Genes' }).layer('l', { type: 'ranges', data: here });
    lost.track('b', { label: 'Sites' }).layer('l', { type: 'ranges', data: here });
    lost.track('c', { label: 'Ruler' }).layer('r', { type: 'ruler' });
  });
  return b.build();
}

describe('a section placed nowhere with an issue', () => {
  it('is one row of the issue, its other tracks not drawn as blank wells, in the draw lists and in the SVG', async () => {
    const st = await settle(spec());
    const f = await st.frame({ size: Size, settled: true });
    const lost = f.layout.sections[1];
    expect(st.view.issues.map((i) => [i.code, i.address])).toEqual([['seq.view.section.system', lost.section.address]]);
    expect(lost.tracks.map((t) => t.box.height)).toEqual([DefaultRenderConfig.track.emptyPx]);
    expect(lost.empties).toEqual(['note']);
    expect(recordDrawLists(f.lists)).toContain(`text "${Issue}"`);
    expect(f.layout.sections[0].note).toBeUndefined();
    const svg = (await renderSvg(spec(), { size: Size })).svg;
    expect(svg.split(`>${Issue}<`).length - 1).toBe(1);
    expect(svg).not.toContain('>Sites<');
  });

  it('only its placement issue collapses it: another issue at its address leaves the section as it was', async () => {
    const st = await settle(spec());
    const lost = st.view.sections[1];
    const other = { ...st.view.issues[0], code: 'seq.view.block.unit', severity: 'warning' as const };
    const view = { ...st.view, issues: [other] };
    const renderer = new Renderer({
      types: BuiltInLayerTypes,
      layouts: BuiltInLayouts,
      model: st.model,
      view,
      data: st.store,
      tokens: LightTokens,
      text: bundledTextMetrics(),
      config: renderConfig({}),
    });
    const f = renderer.frame({ size: Size, settled: true });
    expect(f.layout.sections[1].section.address).toBe(lost.address);
    expect(f.layout.sections[1].note).toBeUndefined();
    expect(f.layout.sections[1].tracks).toHaveLength(3);
  });
});
