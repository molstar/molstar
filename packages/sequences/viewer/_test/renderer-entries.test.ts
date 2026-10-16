/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LY3b (seqspec:viewer#extension-points): the kernel renderers' headless entries in `DefaultRegistry`: their members
// per type, hit locations through the layout's hit test (letters `of` and `row`, ranges in lanes, pairs' arcs and ends,
// a matrix's state band), legends through each channel's theme, `residue_form`, and the pairs geometry the draw,
// overlay and
// hit share.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { case3, case6, case14 } from '../../core/_test/format/cases.js';
import {
  createSeqViewer,
  DefaultRenderProps,
  Hit,
  HoverKind,
  KernelRenderers,
  type LayerRenderer,
  type SeqViewerContext,
  type Theme,
} from '../src/index.js';
import { DefaultSeqViewerSpec } from '../src/default-spec.js';
import type { PairsState } from '../src/prepare.js';
import { Arcs } from '../src/renderer-helpers.js';
import type { MatrixState } from '../src/rows.js';
import type { LegendSink } from '../src/themes.js';
import { layer, opened } from './layout-fixtures.js';
import { allTypes, rangesDoc } from './prepare-fixtures.js';

const K = { renderers: [...KernelRenderers] };

/** The top, height and center y of a track's first entry. */
function band(ctx: SeqViewerContext, section: string, track: number) {
  const l = ctx.layout.get(section)?.current;
  const s = ctx.viewports.get(section)?.live.current;
  if (!l || !s) throw new Error('no layout');
  for (let e = 0; e < l.entryCount; e++) {
    if (l.entryTrack[e] !== track) continue;
    const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - s.scrollY;
    return { top, height: l.entryHeight[e], y: top + l.entryHeight[e] / 2 };
  }
  throw new Error('no entry');
}
const byName = (name: string) => KernelRenderers.find((r) => r.name === name) as LayerRenderer;

describe('the entries', () => {
  it('DefaultRegistry registers one per kernel type and seqv:map, each paired with its format entry', () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    const types = ctx.registry.svs.list('layerTypes').map((e) => e.name);
    expect(ctx.registry.list('renderers').map((r) => r.name)).toEqual(KernelRenderers.map((r) => r.name));
    expect(types.filter((t) => !ctx.registry.get('renderers', t))).toEqual([]);
    ctx.dispose();
  });

  it('members per type: geometries, stability, sizes, miniatures, content', () => {
    const row = (r: LayerRenderer) => [
      r.name,
      r.geometries.join(' '),
      typeof r.stable === 'function' ? 'fn' : r.stable,
      r.minHeight,
      r.miniature ?? true,
      r.content ?? true,
    ];
    expect(KernelRenderers.map(row)).toEqual([
      ['letters', 'linear wrapped circular', 'fn', 6, false, true],
      ['ranges', 'linear wrapped circular', false, 6, true, false],
      ['points', 'linear wrapped circular', false, 6, true, false],
      ['emphasis', 'linear wrapped circular', true, 6, true, true],
      ['variants', 'linear circular', false, 6, true, false],
      ['secondary', 'linear wrapped circular', false, 6, true, false],
      ['values', 'linear wrapped circular', false, 6, true, false],
      ['color', 'linear wrapped circular', true, 6, true, false],
      ['ruler', 'linear circular', false, 6, false, false],
      ['pairs', 'linear circular', false, 6, false, false],
      ['rows', 'linear wrapped', false, 1, false, true],
      ['logo', 'linear circular', false, 8, true, false],
      ['matrix', 'linear', false, 4, true, false],
      ['dot_plot', 'linear', false, 24, false, false],
      ['frame', 'linear wrapped circular', false, 6, false, false],
      ['seqv:map', 'linear', false, 4, false, false],
    ]);
    const h = (name: string, s: unknown = null) => byName(name).height?.(s as never, {}, DefaultRenderProps);
    expect([h('letters'), h('pairs'), h('rows'), h('matrix', { states: 3 }), h('matrix', { states: 40 })]).toEqual([
      22,
      48,
      DefaultRenderProps.rowHeight,
      Math.max(24, 3 * Math.min(18, DefaultRenderProps.rowHeight)),
      160,
    ]);
    expect([byName('variants').splitHeight, byName('matrix').bands?.count({ states: 5 } as never)]).toEqual([8, 5]);
  });

  it("letters: stable, but not with `of` or a ligand's one labeled cell", async () => {
    const stable = byName('letters').stable as (s: unknown, p: unknown) => boolean;
    const t = await opened(case6(), K);
    expect(stable(layer(t.ctx, 'Residues').state, {})).toBe(false);
    const u = await opened(case3(), K);
    const m = u.ctx.graph.assembled.value?.view;
    expect(m?.prepared.some((p) => p.renderer === 'rows')).toBe(true);
    const c14 = await opened(case14(), K);
    expect([stable(layer(c14.ctx, 'Protein').state, {}), stable(layer(c14.ctx, 'DNA').state, {})]).toEqual([
      false,
      true,
    ]);
  });
});

describe('hit locations', () => {
  it("case 14: a residue over its codon gives the protein's location and system; a base, the gene's", async () => {
    const t = await opened(case14(), K);
    const s = t.ctx.viewports.sections.value[0];
    t.vp(s).setView(10, 0, 'api');
    const p = layer(t.ctx, 'Protein');
    const dna = layer(t.ctx, 'DNA');
    const h = Hit.create();
    const at = (track: number, i: number) => {
      t.ctx.layout.hitTest(s, t.xOf(s, 0, i), band(t.ctx, s, track).y, h);
      return [h.kind, p.m.set.strings.get(p.m.set.units.key[h.seq.unit]), h.seq.index, h.system.id === h.seq.unit];
    };
    expect(at(p.t, 13)).toEqual([HoverKind.seq, 'protein/unit', 0, true]);
    expect(at(p.t, 77)).toEqual([HoverKind.seq, 'protein/unit', 11, true]);
    expect(at(p.t, 45)).toEqual([HoverKind.seq, 'protein/unit', 11, true]);
    expect(at(dna.t, 77)).toEqual([HoverKind.seq, 'gene/unit', 77, true]);
  });

  it("ranges packed in lanes in prepare: a hit in each lane's band gives its item, the gap between none", async () => {
    // [10, 30) and [15, 20) overlap: lanes 0 and 1.
    const t = await opened(
      rangesDoc([
        [10, 30],
        [15, 20],
      ]),
      K,
    );
    const s = t.ctx.viewports.sections.value[0];
    const f = layer(t.ctx, 'Features');
    const b = band(t.ctx, s, f.t);
    const h = Hit.create();
    const at = (y: number) => {
      t.ctx.layout.hitTest(s, t.xOf(s, 0, 17), y, h);
      return h.kind === HoverKind.item ? h.item : -1;
    };
    // Lanes of (height - 4) / 2 px from 2 px down, a px of gap above and below each feature.
    const lh = (b.height - 4) / 2;
    expect([at(b.top + 2 + lh / 2), at(b.top + 2 + lh * 1.5)]).toEqual([0, 1]);
    expect([h.system.id, h.system.row]).toEqual([h.annotation, 1]);
    expect(at(b.top + 2 + lh - 0.5)).toBe(-1);
    expect(h.kind).toBe(HoverKind.seq);
    // Points and values: the item at the position; values the position.
    t.ctx.layout.hitTest(s, t.xOf(s, 0, 12), band(t.ctx, s, layer(t.ctx, 'Sites').t).y, h);
    expect([h.kind, h.item]).toEqual([HoverKind.item, 0]);
    t.ctx.layout.hitTest(s, t.xOf(s, 0, 12), band(t.ctx, s, layer(t.ctx, 'Values').t).y, h);
    expect([h.kind, h.seq.index, h.layer]).toEqual([HoverKind.seq, 12, layer(t.ctx, 'Values').k]);
  });

  it('pairs: an arc hit gives its pair (the topmost drawn), an end its pair and side; hover arcs hit once drawn', async () => {
    for (const arcs of ['always', 'hover'] as const) {
      const t = await opened(allTypes({ length: 120, arcs }), K);
      const s = t.ctx.viewports.sections.value[0];
      t.vp(s).setSize(1400, 600, 1);
      t.vp(s).setView(4, 0, 'api');
      const c = layer(t.ctx, 'Contacts');
      const st = c.state as PairsState;
      const v = t.ctx.layout.axis(s);
      if (!v) throw new Error('no axis');
      const b = band(t.ctx, s, c.t);
      v.top = b.top;
      v.height = b.height;
      v.facing = 1;
      // A pair within chain A: its apex, by the shared geometry.
      let k = 0;
      while (st.block[2 * k] !== 0 || st.block[2 * k + 1] !== 0) k++;
      const p = c.m.resolved.layers.params[c.k];
      Arcs.frame(v, st, p);
      Arcs.ends(v, st, k);
      Arcs.apex(st);
      const x = (st.xs[0] + st.xs[1]) / 2;
      const y = st.frame[0] + st.frame[1] * st.frame[9];
      const h = Hit.create();
      t.ctx.layout.hitTest(s, x, y, h);
      if (arcs === 'hover') {
        expect(h.kind).not.toBe(HoverKind.item);
        st.drawn[k] = 1;
        t.ctx.layout.hitTest(s, x, y, h);
      }
      expect([h.kind, h.annotation >= 0, h.side]).toEqual([HoverKind.item, true, -1]);
      expect(st.annotation[k]).toBe(h.annotation);
      // An end: on the ends' line, at its residue.
      t.ctx.layout.hitTest(s, st.xs[0], st.frame[0], h);
      expect([h.kind, h.side]).toEqual([HoverKind.item, 0]);
    }
  });

  it("letters with `row`: an alignment row's residue off a gap, its row system; a gap hits the column", async () => {
    const b = createSeqViewSpec({ title: 'Row' });
    const aln = b.data.alignment({ ref: 'aln', value: ['AC-DE', 'ACKDE'] });
    const v = b.view({ ref: 'v', label: 'V' }).block({ alignment: aln });
    v.track({ label: 'Row 0' }).letters({ row: 0 });
    const t = await opened(b.getDocument(), K);
    const s = t.ctx.viewports.sections.value[0];
    t.vp(s).setView(10, 0, 'api');
    const r = layer(t.ctx, 'Row 0');
    const h = Hit.create();
    const at = (i: number) => {
      t.ctx.layout.hitTest(s, t.xOf(s, 0, i), band(t.ctx, s, r.t).y, h);
      return [h.kind, h.row.row, h.row.index, h.system.row];
    };
    expect(at(3)).toEqual([HoverKind.row, 0, 2, 0]);
    expect(at(2)[0]).toBe(HoverKind.column);
  });

  it('matrix: the state band under the point, from the band function the draw uses', async () => {
    const d = createSeqViewSpec({ title: 'Matrix' });
    const a = d.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR' });
    const value = {
      index: [0, 1, 2],
      value: [
        [0.1, 0.5, 0.9],
        [1, 0, 0.2],
        [0.3, 0.3, 0.3],
      ],
    };
    const m = d.data.annotation({ ref: 'm', geometry: 'matrix', target: a, states: ['x', 'y', 'z'], value } as never);
    d.view({ ref: 'v', label: 'V' }).block({ sequence: a }).track({ label: 'Fingerprint' }).matrix({ data: m });
    const t = await opened(d.getDocument(), K);
    const s = t.ctx.viewports.sections.value[0];
    const fp = layer(t.ctx, 'Fingerprint');
    const st = fp.state as MatrixState;
    expect(st.states).toBeGreaterThan(1);
    const b = band(t.ctx, s, fp.t);
    const h = Hit.create();
    const states: number[] = [];
    for (let k = 0; k < st.states; k++) {
      t.ctx.layout.hitTest(s, t.xOf(s, 0, 0), b.top + ((k + 0.5) * b.height) / st.states, h);
      states.push(h.state);
    }
    expect(states).toEqual(Array.from({ length: st.states }, (_, k) => k));
    expect([h.kind, h.layer]).toEqual([HoverKind.seq, fp.k]);
  });
});

describe('legends, residue forms, the pairs geometry', () => {
  it("a ranges layer's legend lists its categories in palette order with their colors, through its themes", async () => {
    const t = await opened(allTypes({ length: 300 }), K);
    const f = layer(t.ctx, 'Features');
    const items: [string, number][] = [];
    const sink: LegendSink = { item: (label, color) => items.push([label, color]), ramp: () => {}, size: () => {} };
    const pl = f.m.prepared[f.k];
    byName('ranges').legend?.(f.m.resolved.layers.params[f.k], pl.themes as (Theme | null)[], sink, pl.state);
    expect(items.map((x) => x[0])).toEqual(['c0', 'c1', 'c2', 'c3']);
    expect(new Set(items.map((x) => x[1])).size).toBe(4);
  });

  it("a values layer's color over values: a continuous ramp in its legend", async () => {
    const b = createSeqViewSpec({ title: 'Ramp' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR' });
    const d = b.data.annotation({
      ref: 'd',
      geometry: 'values',
      target: a,
      value: { index: [0, 5, 9], value: [0, 5, 20] },
    });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: a }).track({ label: 'Heat' }).color({ data: d });
    const t = await opened(b.getDocument(), K);
    const c = layer(t.ctx, 'Heat');
    const ramps: [number, number, number][] = [];
    const sink: LegendSink = {
      item: () => {},
      ramp: (stops, lo, hi) => ramps.push([stops.length, lo, hi]),
      size: () => {},
    };
    byName('color').legend?.(null, c.m.prepared[c.k].themes as (Theme | null)[], sink, c.state);
    expect(ramps).toEqual([[32, 0, 20]]);
  });

  it('nested arcs nest: heights proportional to span, the widest in view filling the band', async () => {
    const t = await opened(allTypes({ length: 120 }), K);
    const st = layer(t.ctx, 'Contacts').state as PairsState;
    st.frame[3] = 40;
    st.frame[2] = 40;
    const y = (a: number, b: number, x: number) => {
      st.xs[0] = a;
      st.xs[1] = b;
      Arcs.apex(st);
      const u = (x - a) / (b - a);
      return 4 * st.frame[9] * u * (1 - u);
    };
    for (let x = 10; x <= 30; x++) expect(y(10, 30, x)).toBeLessThan(y(0, 40, x));
  });

  it('case 6: contacts below the letters face up toward them (ends at the top, arcs below); `side` wins', async () => {
    const t = await opened(case6(), K);
    const s = t.ctx.viewports.sections.value[0];
    const c = layer(t.ctx, 'Contacts (STI)');
    const v = t.ctx.layout.axis(s);
    const ax = v?.axis;
    if (!v || !ax) throw new Error('no axis');
    expect(ax.facing[c.k]).toBe(1);
    const b = band(t.ctx, s, c.t);
    [v.top, v.height, v.facing] = [b.top, b.height, ax.facing[c.k]];
    const st = c.state as PairsState;
    const p = c.m.resolved.layers.params[c.k];
    Arcs.frame(v, st, p);
    expect([st.frame[0] - b.top < b.height / 2, st.frame[1]]).toEqual([true, 1]);
    Arcs.frame(v, st, { ...p, side: 'above' });
    expect([st.frame[0] - b.top > b.height / 2, st.frame[1]]).toEqual([true, -1]);
  });
});
