/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The viewer half, headless (seqspec:viewer#extension-points). Case 6's poses through `InteractionsViewerRegistry`: the
// fingerprint renderer prepares the KLIFS-style cells from the runtime property (poses 3 and 4 lose T315's and M318's
// hydrogen bond), its bands sized and labeled as a matrix's, its hit naming (ligand, position, the stripe's kind)
// through its `label` member, the legend of kinds; the shared scenes' systems, mixed data and the absent property
//, the track floor, a declared geometry with no draw; the `interactions:kind` themes against their
// kernel fallbacks; without the viewer half, the `matrix` fallback reported; the viewer core importing nothing from
// here. The draws: draws.test.ts; the hit's heap guard: fingerprint.heap.ts.

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import type { DocNode, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { case6 } from '../../../core/_test/format/cases.js';
import { legendGroups } from '../../src/chrome.js';
import { GlyphShape } from '../../src/draw/backend.js';
import { Dirty, SectionComposer } from '../../src/draw/index.js';
import {
  Hit,
  HoverKind,
  KernelRenderers,
  type LayerRenderer,
  type SeqViewerContext,
  ThemeLocation,
} from '../../src/index.js';
import {
  Fingerprints,
  InteractionsRegistry,
  interactionColors,
  interactionGlyphs,
  Kind,
  KindStyle,
  Kinds,
} from '../../src/interactions/index.js';
import {
  FINGERPRINT,
  FingerprintRenderer,
  type FingerprintState,
  InteractionsViewerRegistry,
} from '../../src/interactions/viewer/index.js';
import type { MatrixState } from '../../src/rows.js';
import { LayerFrame, RecordingBackend } from '../../src/testing/index.js';
import { opened } from '../layout-fixtures.js';
import {
  absentScene,
  alignmentScene,
  layerOf,
  manyScene,
  mappingScene,
  mixedScene,
  open,
  type Scene,
  unitScene,
} from './fixtures.js';

const SECTION = 'poses-view.main';
const T315 = 314;
const M318 = 317;
const LIGANDS = ['1iep/C', 'Vina pose 1', 'Vina pose 2', 'Vina pose 3', 'Vina pose 4'];

const poses = (registry = [InteractionsViewerRegistry]) =>
  opened(case6(), { renderers: [...KernelRenderers], view: 'poses-view', registry });
const stateOf = (ctx: SeqViewerContext) => layerOf(ctx).state as FingerprintState;

/** The point (x of a block index, y of a fraction down the track's entry) tested, and what hovering it labels. */
function probe(t: Awaited<ReturnType<typeof open>>, scene: Scene) {
  const l = t.layout(scene.section);
  const e = Array.from(l.entryTrack.subarray(0, l.entryCount)).indexOf(t.trackIndex(scene.track));
  const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e];
  const height = l.entryHeight[e];
  const hit = Hit.create();
  return {
    top,
    height,
    hit: (block: number, index: number, y: number) => {
      t.ctx.layout.hitTest(scene.section, t.xOf(scene.section, block, index), y, hit);
      return hit;
    },
    label: (block: number, index: number, y: number) => {
      t.ctx.layout.hitTest(scene.section, t.xOf(scene.section, block, index), y, hit);
      t.ctx.hover.setFromHit(hit, 'pointer');
      return t.ctx.labels.compute(t.ctx.hover.live.current);
    },
  };
}

describe("case 6's poses through InteractionsViewerRegistry", () => {
  it('prepares the KLIFS-style cells from the runtime property: a stripe per kind, the hydrogen bonds poses 3 and 4 lose', async () => {
    const t = await poses();
    const codes = t.ctx.issues.current.value.map((i) => i.code);
    expect(codes.filter((c) => c.startsWith('svs.render') || c.startsWith('seqv.'))).toEqual([]);
    const s = stateOf(t.ctx);
    const f = s.fp[0];
    if (!f) throw new Error('no fingerprint');
    expect(s.fallback).toBeNull();
    // The kinds present, in stripe order (hydrogen bond, π stacking, cation–π, hydrophobic, the rest): a stripe each.
    expect(Array.from(s.kinds)).toEqual([Kind.hbond, Kind.pi_stacking, Kind.cation_pi, Kind.hydrophobic]);
    expect(s.matrix.labels).toEqual(LIGANDS);
    // Every contacted pocket position has its column; nothing else does; the contacted list is ascending.
    expect(Array.from(s.colOf[0]).filter((c) => c >= 0)).toEqual(Array.from({ length: f.columns.length }, (_, c) => c));
    expect(Array.from(s.contacted[0])).toEqual([...s.contacted[0]].sort((a, b) => a - b));
    expect(s.contacted[0].length).toBe(f.columns.length);
    const kinds = (at: number) => Array.from(LIGANDS, (_, q) => f.cells[q * f.columns.length + s.colOf[0][at]]);
    const hb = 1 << Kind.hbond;
    expect(kinds(T315)).toEqual([hb, hb, hb, 0, 1 << Kind.hydrophobic]);
    expect(kinds(M318)).toEqual([hb, hb, hb, 0, 0]);
    // Against the crystal, poses 3 and 4 lose the hydrogen bond at both.
    for (const pose of [3, 4]) {
      const lost: number[] = [];
      Fingerprints.diff(f, 0, pose, (p, l) => l & hb && lost.push(p), null);
      expect(lost).toEqual(expect.arrayContaining([T315, M318]));
    }
    // A band per ligand, as a matrix's (five of 18 px) under the pocket's numbers and letters.
    expect(t.entries(SECTION).map(([, ref, , , h]) => [ref, h])).toEqual([
      ['poses-view.ruler', 26],
      ['poses-view.pocket', 22],
      ['poses-view.fingerprint', 90],
    ]);
  });

  it("hits (ligand, position, the stripe's kind) and labels it; the legend lists the kinds", async () => {
    const t = await poses();
    const { k } = layerOf(t.ctx);
    const top = t.yOf(SECTION, 'poses-view.fingerprint') - 45;
    /** At DPR 1, 88 px for five rows 18 px apart (16 and a gap of 2), four stripes of 4 px each. */
    const label = (index: number, band: number, stripe: number, part = -1) => {
      const hit = Hit.create();
      const y = top + 1 + band * 18 + stripe * 4 + 2;
      expect(t.ctx.layout.hitTest(SECTION, t.xOf(SECTION, 0, index), y, hit)).toBe(true);
      expect([hit.kind, hit.layer, hit.seq.index, hit.state, hit.part]).toEqual([HoverKind.seq, k, index, band, part]);
      t.ctx.hover.setFromHit(hit, 'pointer');
      return t.ctx.labels.compute(t.ctx.hover.live.current);
    };
    expect(label(T315, 0, 0, Kind.hbond)?.title).toBe('1iep/C · T315 · hbond');
    // Over a stripe whose kind the cell lacks (its ground): the kinds there.
    expect(label(T315, 0, 3)?.title).toBe('1iep/C · T315 · hbond');
    expect(label(T315, 3, 0)?.title).toBe('Vina pose 3 · T315 · no contact');
    expect(label(T315, 4, 3, Kind.hydrophobic)?.title).toBe('Vina pose 4 · T315 · hydrophobic');
    const m318 = label(M318, 4, 1);
    expect(m318?.title).toBe('Vina pose 4 · M318 · no contact');
    expect(m318?.lines).toContainEqual({ key: 'ligand', value: 'Vina pose 4' });
    // Moving between bands at one position is a new hovered element (the state is part of it), as between stripes.
    let v0 = t.ctx.hover.live.version.value;
    label(M318, 0, 1);
    expect(t.ctx.hover.live.version.value).toBe(v0 + 1);
    label(M318, 0, 0, Kind.hbond);
    expect(t.ctx.hover.live.version.value).toBe(v0 + 2);
    v0 = t.ctx.hover.live.version.value;
    label(M318, 0, 0, Kind.hbond);
    expect(t.ctx.hover.live.version.value).toBe(v0);

    const legend = legendGroups(t.ctx).find((g) => g.track === 'Fingerprint');
    const item = (k: keyof typeof KindStyle) => ({
      label: KindStyle[k][0],
      color: KindStyle[k][1],
      glyph: GlyphShape[KindStyle[k][2] as keyof typeof GlyphShape],
    });
    expect(legend?.items).toEqual([item('hbond'), item('pi_stacking'), item('cation_pi'), item('hydrophobic')]);
    // Every block has the property: no ramp.
    expect(legend?.ramps).toEqual([]);
  });

  it("the interactions:kind themes color and mark case 6's contact arcs by kind, as their kernel fallbacks do", async () => {
    // Case 6's chains view with the extension's theme references on its arcs.
    const map = (n: DocNode): DocNode => ({
      ...n,
      params:
        n.kind === 'layer' && n.params?.type === 'pairs'
          ? { ...n.params, color: interactionColors(), glyph: interactionGlyphs() }
          : n.params,
      children: n.children?.map(map),
    });
    const doc: SeqViewSpecDoc = { ...case6(), root: map(case6().root) };
    const t = await opened(doc, { renderers: [...KernelRenderers], registry: [InteractionsViewerRegistry] });
    expect(t.ctx.issues.current.value.filter((i) => i.code === 'svs.render.theme')).toEqual([]);
    const m = t.ctx.graph.assembled.value?.view;
    const k = m?.prepared.findIndex((p) => p.renderer === 'pairs') ?? -1;
    if (!m || k < 0) throw new Error('no pairs layer');
    const a = m.set.annotations[m.resolved.data[m.resolved.layers.dataStart[k]]];
    const [color, glyph] = [m.prepared[k].themes[0], m.prepared[k].themes[1]];
    const loc = ThemeLocation.create();
    loc.annotation = a?.id ?? -1;
    for (let r = 0; a && r < a.count; r++) {
      loc.item = r;
      const kind = m.set.strings.get(a.category?.[r] ?? -1) as keyof typeof KindStyle;
      expect([color?.value(loc), glyph?.value(loc)]).toEqual([
        Number.parseInt(KindStyle[kind][1].slice(1), 16),
        ['circle', 'diamond', 'triangle', 'square', 'pin', 'star', 'cross', 'tick'].indexOf(KindStyle[kind][2]),
      ]);
    }
  });
});

describe('one provider serves color and glyph; a missing category is `other` everywhere', () => {
  it("the theme's color, glyph and legend per row of the contacts equal the kernel fallbacks', null categories as `other`", async () => {
    const scene = unitScene({ arcs: true });
    const withExt = await open(scene, { registry: [InteractionsViewerRegistry] });
    const without = await open(scene, { registry: [{ svs: InteractionsRegistry }] });
    // The arcs are the scene's second view.
    for (const t of [withExt, without]) {
      t.ctx.viewports.get('w.main')?.setSize(1000, 400, 1);
      t.ctx.views.setCurrent('w');
      await t.ctx.graph.update();
    }
    const themes = (ctx: SeqViewerContext) => {
      const m = ctx.graph.assembled.value?.view;
      const k = m?.prepared.findIndex((p) => p.renderer === 'pairs') ?? -1;
      if (!m || k < 0) throw new Error('no pairs layer');
      return { m, themes: m.prepared[k].themes, annotation: m.resolved.data[m.resolved.layers.dataStart[k]] };
    };
    const a = themes(withExt.ctx);
    const b = themes(without.ctx);
    // Without the extension each theme falls back (reported), to the kernel's field_palette and field_map.
    expect(without.ctx.issues.current.value.filter((i) => i.code === 'svs.render.theme').length).toBeGreaterThan(0);
    const rows = withExt.ctx.graph.assembled.value?.set.annotations[a.annotation]?.count ?? 0;
    expect(rows).toBeGreaterThan(8);
    const loc = ThemeLocation.create();
    const cats = new Set<string>();
    for (let r = 0; r < rows; r++) {
      loc.annotation = a.annotation;
      loc.item = loc.row = r;
      const set = a.m.set;
      cats.add(set.strings.get(set.annotations[a.annotation]?.category?.[r] ?? -1));
      const other = Object.assign(ThemeLocation.create(), { annotation: b.annotation, item: r, row: r });
      expect([a.themes[0]?.value(loc), a.themes[1]?.value(loc)]).toEqual([
        b.themes[0]?.value(other),
        b.themes[1]?.value(other),
      ]);
    }
    // The document has rows of a kind this vocabulary does not know and rows with none: both `other`.
    expect(cats.has('weird')).toBe(true);
    expect(cats.has('')).toBe(true);
    const other = Object.assign(ThemeLocation.create(), { category: Kind.other });
    expect(a.themes[0]?.value(other)).toBe(Number.parseInt(KindStyle.other[1].slice(1), 16));
    withExt.ctx.dispose();
    without.ctx.dispose();
  });
});

describe("the shared scenes' systems", () => {
  it("a unit system: the fingerprint of one chain's contacts, hit and labeled in the unit's numbering", async () => {
    const scene = unitScene();
    const t = await open(scene);
    const p = probe(t, scene);
    // 72 px for four rows 18 apart; five stripes (hbond, π, cation–π, hydrophobic, other) of 3 px a row of 16.
    expect(p.height).toBe(72);
    const h = p.hit(0, 7, p.top + 1 + 0 * 18 + 2);
    expect([h.kind, h.state, h.part]).toEqual([HoverKind.seq, 0, Kind.hbond]);
    // Pose 3 lost the bond at residue 7; a row's gap belongs to the row above it.
    const lost = p.label(0, 7, p.top + 1 + 2 * 18 + 2);
    expect(lost?.title).toBe('L pose 3 · K8 · no contact');
    expect(p.hit(0, 7, p.top + 1 + 16 + 1).state).toBe(0);
    // Rows with a kind of no vocabulary entry, and none: `other` both.
    const odd = p.label(0, 7, p.top + 1 + 3 * 18 + 14);
    expect(odd?.title).toMatch(/^L pose 4 · K8 · (other|no contact)$/);
    t.ctx.dispose();
  });

  it('an alignment system: the hover is a column, titled `column N`, and the row is the pose', async () => {
    const scene = alignmentScene();
    const t = await open(scene);
    const p = probe(t, scene);
    const h = p.hit(0, 10, p.top + 1 + 2);
    expect([h.kind, h.state, h.index]).toEqual([HoverKind.column, 0, 10]);
    expect(p.label(0, 10, p.top + 1 + 2)?.title).toMatch(/^L pose 1 · column 11 · /);
    // Pose 3 contacts the second chain's residue 10, which sits in column 11.
    expect(p.label(0, 11, p.top + 1 + 2 * 18 + 2)?.title).toMatch(/^L pose 3 · column 12 · hbond$/);
    t.ctx.dispose();
  });

  it("a block reached through a mapping: the fingerprint's positions land on the mapped block's", async () => {
    const scene = mappingScene();
    const t = await open(scene);
    const s = stateOf(t.ctx);
    // System positions 8, 12, 20 and 30 are the block's 3, 7, 15 and 25.
    expect(Array.from(s.contacted[0])).toEqual([3, 7, 15, 25]);
    const p = probe(t, scene);
    expect(p.hit(0, 3, p.top + 1 + 2).state).toBe(0);
    expect(p.hit(0, 4, p.top + 1 + 2).part).toBe(-1);
    t.ctx.dispose();
  });
});

describe('rows and sizing', () => {
  it('20 poses × 8 kinds: 18 px rows of 16, every stripe 2 px; a track `height` hint raised to a device pixel a row', async () => {
    const t = await open(manyScene());
    expect(t.entries('v.main').at(-1)?.[4]).toBe(360);
    const s = stateOf(t.ctx);
    expect(s.kinds.length).toBe(8);
    // No 160 cap: 20 rows × 18.
    expect(FingerprintRenderer.height?.(s, {}, t.ctx.layout.props('v.main'))).toBe(360);
    expect(FingerprintRenderer.minTrackHeight?.(s, {})).toBe(22);
    t.ctx.dispose();
    // A hint of one row gives a row's height for 40 rows: raised to 42 (S + 2); without the floor it stays one row's.
    const rowHeight = t.ctx.layout.props('v.main').rowHeight;
    const u = await open(manyScene(40, 3, 1));
    expect(u.entries('v.main').at(-1)?.[4]).toBe(42);
    const bare = { ...FingerprintRenderer, minTrackHeight: undefined } as LayerRenderer;
    const v = await open(manyScene(40, 3, 1), { registry: [{ ...InteractionsViewerRegistry, renderers: [bare] }] });
    expect(v.entries('v.main').at(-1)?.[4]).toBe(rowHeight);
    // A taller hint wins as before.
    const w = await open(manyScene(20, 8, 40));
    expect(w.entries('v.main').at(-1)?.[4]).toBe(40 * rowHeight);
    for (const x of [u, v, w]) x.ctx.dispose();
  });

  it("a track without a hint is not raised, and few poses keep the matrix's 24 px floor", async () => {
    const t = await open(unitScene({ contacts: [[[3, 'hbond']]] }));
    expect(t.entries('v.main').at(-1)?.[4]).toBe(24);
    t.ctx.dispose();
  });
});

describe('mixed data and the absent property', () => {
  it("each fingerprint bounds its own rows, and the block without one is hit as a matrix's bands", async () => {
    const scene = mixedScene();
    const t = await open(scene);
    const s = stateOf(t.ctx);
    expect(Array.from(s.rowsOf)).toEqual([5, 3, 0]);
    expect(s.fallback?.rows.map((r) => r !== null)).toEqual([false, false, true]);
    const p = probe(t, scene);
    const row = (q: number) => p.top + 1 + q * 18 + 8;
    // The fingerprint of 5 rows hits row 4; the cut one has no row 3: nothing of the layer's own there.
    expect(p.hit(0, 3, row(4)).state).toBe(4);
    expect(p.hit(1, 4, row(2)).state).toBe(2);
    const past = p.hit(1, 4, row(3));
    expect([past.state, past.part]).toEqual([-1, -1]);
    // The matrix block: its band under the point, at any position of its block.
    expect(p.hit(2, 0, row(4)).state).toBe(4);
    t.ctx.dispose();
  });

  it("labels: the value line of the block without the property, the fingerprint's title elsewhere", async () => {
    const scene = mixedScene();
    const t = await open(scene);
    const p = probe(t, scene);
    const row = (q: number) => p.top + 1 + q * 18 + 8;
    const plain = p.label(2, 3, row(1));
    // Plain matrix state values: (3 + 1) % 3 → 1 + 1 = 2.
    expect(plain?.lines).toContainEqual({ key: 'a pose 2', value: '2' });
    expect(plain?.title).not.toMatch(/·/);
    expect(p.label(0, 3, row(0))?.title).toMatch(/^a pose 1 · /);
    // Nothing of the fingerprint's own past its cut rows.
    expect(p.label(1, 4, row(4))?.title ?? '').not.toMatch(/·/);
    t.ctx.dispose();
  });

  it("the value line of a block with fewer states than the layer's bands: read with its own stride, none past them", async () => {
    // A resolved layer's matrices share their states; a host's own `LayerInput` need not: the plain block (C) made a
    // matrix of 3 states under the layer's 5 bands, its values 100 + 10 · row + state.
    const scene = mixedScene();
    const t = await open(scene);
    const s = stateOf(t.ctx);
    const values = Float64Array.from({ length: 4 * 3 }, (_, i) => 100 + 10 * Math.floor(i / 3) + (i % 3));
    for (const m of [s.matrix, s.fallback as MatrixState]) {
      const ann = m.ann[2] as NonNullable<MatrixState['ann'][number]>;
      (m.ann as unknown[])[2] = { ...ann, matrix: { ...ann.matrix, values } };
      m.stride[2] = 3;
    }
    expect([s.matrix.states, s.fallback?.rows[2]?.[4]]).toEqual([5, 1]);
    const p = probe(t, scene);
    const row = (q: number) => p.top + 1 + q * 18 + 8;
    const line = (q: number) =>
      (p.label(2, 4, row(q))?.lines ?? []).filter((l) => l.key === `a pose ${q + 1}`).map((l) => l.value);
    // Position 4 is the matrix's row 1.
    expect([0, 1, 2].map(line)).toEqual([['110'], ['111'], ['112']]);
    expect([3, 4].map(line)).toEqual([[], []]);
    t.ctx.dispose();
  });

  it("the legend lists the kinds and, as there are blocks without the property, the counts' ramp, captioned", async () => {
    const t = await open(mixedScene());
    const g = legendGroups(t.ctx).find((x) => x.track === 'FP');
    expect(g?.items.map((i) => i.label)).toEqual(['hydrogen bond', 'salt bridge', 'π stacking', 'hydrophobic contact']);
    expect(g?.ramps.map((r) => [r.caption, r.lo])).toEqual([['contacts (number of kinds)', 0]]);
    t.ctx.dispose();
  });

  it("the property absent: the layer is its matrix's counts in Blues, the legend the ramp alone (never empty)", async () => {
    const t = await open(absentScene());
    const s = stateOf(t.ctx);
    expect([s.kinds.length, s.fallback?.rows.length, s.matrix.states]).toEqual([0, 1, 4]);
    const g = legendGroups(t.ctx).find((x) => x.track === 'FP');
    expect([g?.items.length, g?.ramps.length, g?.ramps[0]?.caption]).toEqual([0, 1, 'contacts (number of kinds)']);
    // Blues from 0 (light) to the largest count (dark).
    const [lo, hi] = [g?.ramps[0]?.stops[0].color, g?.ramps[0]?.stops.at(-1)?.color];
    expect(lo).not.toBe(hi);
    t.ctx.dispose();
  });
});

describe('a declared geometry with no draw', () => {
  it('is reported once over several frames, as svs.render.fallback (a warning), and drawn as bars', async () => {
    const bare = { ...FingerprintRenderer, draw: undefined, tiny: undefined } as LayerRenderer;
    const t = await opened(case6(), {
      renderers: [...KernelRenderers],
      view: 'poses-view',
      registry: [{ ...InteractionsViewerRegistry, renderers: [bare] }],
    });
    for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 400, 1);
    await t.ctx.graph.update();
    const composer = new SectionComposer(t.ctx, SECTION);
    let bars = 0;
    for (let i = 0; i < 4; i++) {
      const r = new RecordingBackend();
      composer.frame(r, i === 0 ? Dirty.tracks : Dirty.tracks);
      composer.redraw();
      bars = Array.from(r.rectItems.color.subarray(0, r.rectItems.n)).filter((c) => c === -6).length;
    }
    // Drawn otherwise than written: a fallback, not a failure (GV-10).
    const said = t.ctx.issues.current.value.filter((i) => i.code === 'svs.render.fallback');
    expect(said.map((e) => [e.severity, e.message])).toEqual([
      ['warning', expect.stringContaining('declares linear and has no draw')],
    ]);
    expect(bars).toBeGreaterThan(0);
    t.ctx.dispose();
  });
});

describe('in a miniature', () => {
  it("`seqv:map` draws the fingerprint's tracks in lanes: the kinds' colors are in the map's image", async () => {
    const scene = unitScene({ map: true });
    const t = await open(scene);
    const m = t.ctx.graph.assembled.value?.view;
    const k = m?.prepared.findIndex((p) => p.renderer === 'seqv:map') ?? -1;
    const entry = t.ctx.registry.get('renderers', 'seqv:map');
    if (!entry) throw new Error('no map');
    const f = new LayerFrame(t.ctx, entry, scene.section, k);
    f.draw(f.direct);
    const layer = f.direct.images.content.find((c): c is RecordingBackend => c instanceof RecordingBackend);
    const colors = new Set(Array.from(layer?.rectItems.color.subarray(0, layer.rectItems.n) ?? []));
    expect(t.ctx.issues.current.value.filter((i) => i.code === 'seqv.renderer.error')).toEqual([]);
    for (const kind of [Kind.hbond, Kind.hydrophobic])
      expect(colors.has(Number.parseInt(KindStyle[Kinds[kind]][1].slice(1), 16))).toBe(true);
    t.ctx.dispose();
  });
});

describe('without the viewer half', () => {
  it('the format half alone: the layer is drawn as its `matrix` fallback, reported once', async () => {
    const t = await poses([{ svs: InteractionsRegistry }] as never);
    expect(t.ctx.graph.assembled.value?.view?.prepared.map((p) => p.renderer)).toContain('matrix');
    expect(layerOf(t.ctx).k).toBe(-1);
    const fallback = t.ctx.issues.current.value.filter((i) => i.code === 'svs.render.fallback');
    expect(fallback.map((i) => i.message)).toEqual([expect.stringContaining(FINGERPRINT)]);
  });

  it('the viewer core imports nothing from this package (and Biome refuses it: layers.grit)', () => {
    const root = new URL('../../src/', import.meta.url);
    const files: string[] = [];
    const walk = (dir: string) => {
      for (const name of readdirSync(dir)) {
        const p = `${dir}/${name}`;
        if (statSync(p).isDirectory()) {
          if (name !== 'interactions') walk(p);
        } else if (/\.(ts|tsx)$/.test(name)) files.push(p);
      }
    };
    walk(fileURLToPath(root).replace(/[\\/]$/, ''));
    expect(files.length).toBeGreaterThan(40);
    const imports = /(?:from|import)\s*\(?\s*["']([^"']+)["']/g;
    const offending = files.flatMap((f) =>
      [...readFileSync(f, 'utf8').matchAll(imports)]
        .map((m) => m[1])
        .filter((s) => s.startsWith('@molstar/seq-viewer/interactions') || /(^|\/)interactions\//.test(s))
        .map((s) => `${f}: ${s}`),
    );
    expect(offending).toEqual([]);
    const grit = readFileSync(new URL('../../../layers.grit', import.meta.url), 'utf8');
    expect(grit.includes('@molstar/seq-viewer/interactions')).toBe(true);
  });
});
