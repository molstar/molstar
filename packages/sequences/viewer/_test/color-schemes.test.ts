/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Color schemes in the viewer (seqspec:interaction#color-schemes-at-runtime): the pickers (a legend group's title, a
// track's menu) list only the schemes that fit a layer's color; the runtime override recolors through the theme stage,
// survives a patch and resets on a new document; the action and the command; a dark skin's variant; the scheme's
// legend.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { SeqSets } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import { dataSpan, editRange, legendGroups, resetRange, trackMenu } from '../src/chrome.js';
import { colorSchemePickers, pickedTheme } from '../src/color-schemes.js';
import { KernelRenderers, SeqViewerConfig, type SeqViewerContext, type Theme } from '../src/index.js';
import { layer, opened } from './layout-fixtures.js';
import { at, context } from './theme-fixtures.js';

const code = (ch: string) => ch.charCodeAt(0);
const color = (ctx: SeqViewerContext, label: string) => {
  const r = layer(ctx, label);
  return { theme: r.m.prepared[r.k].themes[0] as Theme, ref: r.m.set.strings.get(r.m.resolved.layers.ref[r.k]), r };
};
const options = (ctx: SeqViewerContext, label: string) =>
  colorSchemePickers(ctx, layer(ctx, label).t).map((p) => p.options.map((o) => o.value));
const settle = (ctx: SeqViewerContext) => ctx.graph.update();

function dna() {
  const b = createSeqViewSpec({ title: 'DNA' });
  const s = b.data.sequence({ ref: 's', alphabet: 'dna', value: 'ACGTACGGTCAGTTAC' });
  b.view({ ref: 'v', label: 'V' }).block({ sequence: s }).track({ label: 'DNA' }).letters();
  return b.getDocument();
}

describe('the scheme theme', () => {
  it("is the theme its scheme stands for: the same values; another subject's falls back, reported", () => {
    const set = SeqSets.build({
      molecules: [{ key: 'p', residues: 'ACDEFGHIKLMNPQRSTVWY' }],
      units: [{ key: 'A', molecule: 'p' }],
    }).set;
    const c = context(set, { type: 'letters' });
    const scheme = c.color('scheme', { name: 'taylor', overrides: { W: '#123456' } });
    const letter = c.color('letter', { colors: 'Taylor' });
    for (const ch of 'ARNDCQEGHILKMFPSTYV')
      expect(scheme.value(at({ block: 0, letter: code(ch) }))).toBe(letter.value(at({ block: 0, letter: code(ch) })));
    expect(scheme.value(at({ block: 0, letter: code('W') }))).toBe(0x123456);
    expect(c.issues).toEqual([]);
    const off = c.color('scheme', { name: 'viridis' });
    expect(off.value(at({ block: 0, letter: code('K') }))).toBe(0x6464ff); // `letter`: Zappo for protein
    expect(c.issues.map((i) => i.code)).toEqual(['svs.render.theme']);
  });
});

describe('a picked scheme that stops fitting (VS-15)', () => {
  it("draws the layer's own color, not `letter`", () => {
    const set = SeqSets.build({
      molecules: [{ key: 'p', residues: 'ACDEFGHIKLMNPQRSTVWY' }],
      units: [{ key: 'A', molecule: 'p' }],
    }).set;
    const c = context(set, { type: 'letters' });
    const own = { theme: 'uniform', params: { color: '#123456' } };
    const ref = pickedTheme(own, 'viridis');
    const picked = c.color(ref.theme, ref.params as Record<string, unknown>, ref.fallback ?? undefined);
    expect(picked.value(at({ block: 0, letter: code('K') }))).toBe(0x123456);
  });
});

describe('the picker', () => {
  it("lists only the schemes that fit: a protein alignment's rows, its points, its values; a DNA sequence", async () => {
    const t = await opened(case3());
    const [rows, points] = options(t.ctx, 'Members');
    expect(rows).toEqual([null, 'clustal_x', 'zappo', 'taylor', 'hydrophobicity', 'percent_identity']);
    expect(points).toEqual([null, 'category10', 'okabe_ito', 'tableau10']);
    expect(options(t.ctx, 'Conservation')).toEqual([[null, 'viridis', 'cividis', 'magma', 'blues', 'rd_bu']]);
    const labels = colorSchemePickers(t.ctx, layer(t.ctx, 'Members').t)[0].options.map((o) => o.label);
    expect(labels).toContain('Percent identity (colorblind-safe)');
    // The members' track lists nothing: no legend group; its two pickers are in its menu, each named by its layer type.
    expect(legendGroups(t.ctx).find((g) => g.ref === 'family.members')).toBeUndefined();
    const drill = { layer: null, open: () => {} };
    expect(trackMenu(t.ctx, 'family.main', 'family.members', drill).map((i) => i.label)).toContain('Colors, rows ▸');
    expect(trackMenu(t.ctx, 'family.main', 'family.members', drill).map((i) => i.label)).toContain('Colors, points ▸');
    t.ctx.dispose();
    const d = await opened(dna());
    expect(options(d.ctx, 'DNA')).toEqual([[null, 'nucleotide', 'purine_pyrimidine', 'nucleotide_cb']]);
    d.ctx.config.set(SeqViewerConfig.UI.ColorSchemes, false);
    await d.ctx.setDocument(dna()); // a new view, read again under the config
    expect(legendGroups(d.ctx)).toEqual([]);
    d.ctx.dispose();
  });
});

describe('the override', () => {
  it('recolors through the theme stage; survives a patch; a new document resets it', async () => {
    const t = await opened(case3());
    const { ref } = color(t.ctx, 'Members');
    const before = color(t.ctx, 'Members').theme.value(at({ block: 0, column: 3, letter: code('R') }));
    expect(t.ctx.actions.run('layer.colorScheme', { layer: ref, scheme: 'taylor' })).toBe(true);
    await settle(t.ctx);
    const taylor = color(t.ctx, 'Members');
    expect(taylor.theme.value(at({ block: 0, column: 3, letter: code('R') }))).toBe(0x0000ff);
    expect(taylor.theme.value(at({ block: 0, column: 9, letter: code('A') }))).toBe(0xccff00);
    expect(colorSchemePickers(t.ctx, taylor.r.t)[0].value).toBe('taylor');
    // A patch elsewhere keeps it.
    const r = await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'family', params: { label: 'F' } }],
    });
    expect(r.ok).toBe(true);
    expect(color(t.ctx, 'Members').theme.value(at({ block: 0, column: 3, letter: code('R') }))).toBe(0x0000ff);
    expect(t.ctx.ui.colorSchemes.value.get(ref)).toBe('taylor');
    // Back to the document's own.
    t.ctx.actions.run('layer.colorScheme', { layer: ref, scheme: null });
    await settle(t.ctx);
    expect(color(t.ctx, 'Members').theme.value(at({ block: 0, column: 3, letter: code('R') }))).toBe(before);
    // A new document clears every override.
    t.ctx.actions.run('layer.colorScheme', { layer: ref, scheme: 'zappo' });
    await t.ctx.setDocument(case3());
    expect(t.ctx.ui.colorSchemes.value.size).toBe(0);
    t.ctx.dispose();
  });

  it('the command: a link switches it, a scheme that does not fit is refused, write: true patches the document', async () => {
    const t = await opened(case3());
    const { ref } = color(t.ctx, 'Members');
    const run = (params: Record<string, unknown>) =>
      t.ctx.commands.run({ command: 'seqv:colorScheme', params: { layer: ref, ...params } }, 'host:links');
    expect(await run({ scheme: 'zappo' })).toBe(true);
    await settle(t.ctx);
    expect(color(t.ctx, 'Members').theme.value(at({ block: 0, column: 3, letter: code('K') }))).toBe(0x6464ff);
    await run({ scheme: 'viridis' }); // runs, but the action refuses: letters are not values
    expect(t.ctx.ui.colorSchemes.value.get(ref)).toBe('zappo');
    expect(await t.ctx.commands.run({ command: 'seqv:colorScheme', params: { scheme: 'zappo' } }, 'api')).toBe(false);
    await run({ scheme: 'hydrophobicity', write: true });
    await settle(t.ctx);
    const node = (n: { ref?: string; params?: Record<string, unknown>; children?: unknown[] }): unknown =>
      n.ref === ref ? n.params?.color : (n.children ?? []).map((c) => node(c as never)).find(Boolean);
    expect(node(t.ctx.document.shown.value?.root as never)).toEqual({
      theme: 'scheme',
      params: { name: 'hydrophobicity', field: null, domain: null },
    });
    expect(t.ctx.ui.colorSchemes.value.has(ref)).toBe(false);
    expect(color(t.ctx, 'Members').theme.value(at({ block: 0, column: 3, letter: code('I') }))).toBe(0xff0000);
    t.ctx.dispose();
  });

  it("a picked scheme that reads the residues runs again after a residue edit, the document's color reading none (VS-17)", async () => {
    const b = createSeqViewSpec({ title: 'Seq' });
    const s = b.data.sequence({ ref: 's', alphabet: 'protein', value: 'MKVLAAGWKR' });
    b.view({ ref: 'v', label: 'V' })
      .block({ sequence: s })
      .track({ label: 'Seq' })
      .letters({ color: { theme: 'letter' } } as never);
    const t = await opened(b.getDocument(), { renderers: [...KernelRenderers] });
    const { ref } = color(t.ctx, 'Seq');
    const schemes = options(t.ctx, 'Seq').flat();
    const scheme = schemes.find((x) => x !== 'letter' && x !== null) as string;
    expect(t.ctx.actions.run('layer.colorScheme', { layer: ref, scheme })).toBe(true);
    await settle(t.ctx);
    const before = color(t.ctx, 'Seq').theme;
    await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 's', params: { value: 'MKVLAAGWKW' } }],
    });
    await settle(t.ctx);
    expect(color(t.ctx, 'Seq').theme).not.toBe(before);
    t.ctx.dispose();
  });

  it('a dark skin draws the dark variant; the field and domain carry over', async () => {
    const t = await opened(case3());
    const { ref } = color(t.ctx, 'Conservation');
    t.ctx.actions.run('layer.colorScheme', { layer: ref, scheme: 'blues' });
    await settle(t.ctx);
    const hi = color(t.ctx, 'Conservation').r.m.resolved.layers.domainMax[color(t.ctx, 'Conservation').r.k];
    expect(color(t.ctx, 'Conservation').theme.value(at({ value: hi }))).toBe(0x08306b);
    t.ctx.setSkin({ ...t.ctx.skin.value, dark: true });
    await settle(t.ctx);
    expect(color(t.ctx, 'Conservation').theme.value(at({ value: hi }))).toBe(0xf7fbff);
    t.ctx.dispose();
  });

  it("a scheme's labeled groups are its legend", async () => {
    const d = await opened(dna(), { renderers: [...KernelRenderers] });
    const { ref } = color(d.ctx, 'DNA');
    d.ctx.actions.run('layer.colorScheme', { layer: ref, scheme: 'purine_pyrimidine' });
    await settle(d.ctx);
    const g = legendGroups(d.ctx).find((x) => x.track === 'DNA');
    expect(g?.items.map((i) => [i.label, i.color])).toEqual([
      ['purine', '#ff83fa'],
      ['pyrimidine', '#40e0d0'],
    ]);
    expect([g?.colors?.value, g?.colors?.label]).toEqual(['purine_pyrimidine', 'Purine/pyrimidine']);
    d.ctx.dispose();
  });
});

/** A color track over values 0..20 (the document's own scale: sequential from 0); `palette`: a `field_palette`'s. */
function heat(palette?: Record<string, unknown>) {
  const b = createSeqViewSpec({ title: 'Heat' });
  const s = b.data.sequence({ ref: 's', value: 'MKTAYIAKQRQISFVKSHFS' });
  const d = b.data.annotation({
    ref: 'd',
    geometry: 'values',
    target: s,
    value: { index: [0, 5, 9], value: [0, 5, 20] },
  });
  const own = palette ? { color: { theme: 'field_palette', params: { field: 'value', palette } } } : {};
  b.view({ ref: 'v', label: 'V' })
    .block({ sequence: s })
    .track({ label: 'Heat' })
    .color({ data: d, ...own } as never);
  return b.getDocument();
}
describe("the color range (the legend's editor)", () => {
  const group = (ctx: SeqViewerContext) => legendGroups(ctx).find((g) => g.track === 'Heat')?.range;
  const value = (ctx: SeqViewerContext, x: number) => color(ctx, 'Heat').theme.value(at({ value: x }));

  it("is wider or narrower than the data, the ends' colors beyond it; a bad edit is refused; reset", async () => {
    const t = await opened(heat(), { renderers: [...KernelRenderers] });
    const r0 = group(t.ctx);
    if (!r0) throw new Error('no range');
    expect([r0.edited, r0.lo, r0.hi, r0.dataLo, r0.dataHi]).toEqual([false, r0.lo, r0.hi, r0.lo, r0.hi]);
    const own = [value(t.ctx, r0.lo), value(t.ctx, r0.hi)];
    // Wider: 0..20 becomes 0..40 (and the data's span is marked inside the ramp); beyond it the end color, not none.
    expect(editRange(t.ctx, r0, 'hi', '40')).toBe(true);
    await settle(t.ctx);
    const wide = group(t.ctx);
    expect([wide?.edited, wide?.lo, wide?.hi, wide?.dataHi]).toEqual([true, r0.lo, 40, r0.dataHi]);
    const span = wide && dataSpan(wide);
    expect(span?.[0]).toBe(0);
    expect(span?.[1]).toBeCloseTo((r0.dataHi - r0.lo) / (40 - r0.lo));
    expect(value(t.ctx, r0.hi)).not.toBe(own[1]);
    expect([value(t.ctx, 40), value(t.ctx, 400)]).toEqual([own[1], own[1]]);
    expect(value(t.ctx, -40)).toBe(own[0]);
    // Narrower than the data: the data beyond the ends takes the end colors.
    if (!wide) throw new Error('no range');
    expect(editRange(t.ctx, wide, 'lo', String(r0.hi / 2))).toBe(true);
    await settle(t.ctx);
    expect([value(t.ctx, r0.lo), value(t.ctx, r0.hi / 2)]).toEqual([own[0], own[0]]);
    // Refused: not a number, the ends out of order.
    const now = group(t.ctx);
    if (!now) throw new Error('no range');
    expect([editRange(t.ctx, now, 'hi', 'abc'), editRange(t.ctx, now, 'hi', String(now.lo))]).toEqual([false, false]);
    expect(editRange(t.ctx, now, 'lo', '')).toBe(false);
    expect(t.ctx.ui.colorRanges.value.size).toBe(1);
    resetRange(t.ctx, now);
    await settle(t.ctx);
    expect([group(t.ctx)?.edited, value(t.ctx, r0.hi)]).toEqual([false, own[1]]);
    t.ctx.dispose();
  });

  it("a discrete palette's ends' colors beyond the range; none offered on an absolute palette, whose ramp spans its checkpoints (VS-14)", async () => {
    const rgb = ['#ff0000', '#00ff00', '#0000ff'];
    const d = await opened(heat({ kind: 'discrete', colors: rgb }), { renderers: [...KernelRenderers] });
    const r = group(d.ctx);
    if (!r) throw new Error('no range');
    expect(editRange(d.ctx, r, 'hi', '10')).toBe(true);
    await settle(d.ctx);
    expect([value(d.ctx, 20), value(d.ctx, 5), value(d.ctx, -5)]).toEqual([0x0000ff, 0x00ff00, 0xff0000]);
    d.ctx.dispose();
    const absolute = {
      kind: 'discrete',
      mode: 'absolute',
      colors: [
        ['#ff0000', 0],
        ['#0000ff', 10],
      ],
    };
    const a = await opened(heat(absolute), { renderers: [...KernelRenderers] });
    expect(group(a.ctx)).toBeNull();
    a.ctx.dispose();
    // Features scored 10..60 under an absolute ramp over 0..100: the legend's ramp is the palette's, not 0..1.
    const b = createSeqViewSpec({ title: 'Scores' });
    const s = b.data.sequence({ ref: 's', value: 'MKTAYIAKQRQISFVKSHFS' });
    const f = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: s,
      value: { start: [0, 5], stop: [4, 9], score: [10, 60] },
    });
    const ramp = {
      kind: 'continuous',
      mode: 'absolute',
      colors: [
        ['#ff0000', 0],
        ['#0000ff', 100],
      ],
    };
    const color = { theme: 'field_palette', params: { field: 'score', palette: ramp } };
    b.view({ ref: 'v', label: 'V' })
      .block({ sequence: s })
      .track({ label: 'Scores' })
      .ranges({ data: f, color } as never);
    const c = await opened(b.getDocument(), { renderers: [...KernelRenderers] });
    const g = legendGroups(c.ctx).find((x) => x.track === 'Scores');
    expect(g?.range).toBeNull();
    expect([g?.ramps[0]?.lo, g?.ramps[0]?.hi]).toEqual([0, 100]);
    c.ctx.dispose();
  });

  it("the command; a patch keeps it; a new document clears it; it sets a scheme's domain too", async () => {
    const t = await opened(heat(), { renderers: [...KernelRenderers] });
    const { ref } = color(t.ctx, 'Heat');
    const run = (range: unknown) =>
      t.ctx.commands.run({ command: 'seqv:colorRange', params: { layer: ref, range } }, 'host:links');
    expect(await run([0, 3])).toBe(true);
    await run([2, 1]); // refused: the action says no
    expect(t.ctx.ui.colorRanges.value.get(ref)).toEqual([0, 3]);
    t.ctx.actions.run('layer.colorScheme', { layer: ref, scheme: 'viridis' });
    await settle(t.ctx);
    expect([group(t.ctx)?.lo, group(t.ctx)?.hi]).toEqual([0, 3]);
    expect(value(t.ctx, 30)).toBe(value(t.ctx, 3));
    await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'v', params: { label: 'F' } }],
    });
    expect(t.ctx.ui.colorRanges.value.get(ref)).toEqual([0, 3]);
    await run(null);
    expect(t.ctx.ui.colorRanges.value.size).toBe(0);
    await run([0, 3]);
    await t.ctx.setDocument(heat());
    expect(t.ctx.ui.colorRanges.value.size).toBe(0);
    t.ctx.dispose();
  });
});
