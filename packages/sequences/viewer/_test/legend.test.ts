/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The legend (seqspec:interaction#the-legend; UX1): a group per track with something to list, keyed `section/trackRef`,
// tracks with identical legends merged into one (titled by their common name and count); every layer's channels in it
// (colors, glyphs, a ramp for a continuous palette, sample sizes); at most `UI.LegendItems` items then "+n more", at
// most `UI.LegendGroups` groups then "+n more", expanded through `ui.legendOpen`; unnamed items' one color left out.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { case9 } from '../../core/_test/format/cases.js';
import {
  type ColorsDrill,
  groupTitle,
  legendGroups,
  legendShown,
  legendShownGroups,
  schemeItems,
  trackMenu,
} from '../src/chrome.js';
import {
  createSeqViewer,
  KernelRenderers,
  type LayerRenderer,
  SeqViewerConfig,
  type SeqViewerContext,
} from '../src/index.js';
import { headerLegend } from '../src/legend.js';
import { StarsFormat, StarsRenderer } from '../src/testing/index.js';
import { SizeUnit } from '../src/themes.js';
import { starsDoc } from './conformance-fixtures.js';
import { opened } from './layout-fixtures.js';
import { distinctDoc, kindsDoc } from './legend-fixtures.js';

const K = { renderers: [...KernelRenderers] };
const STARS = { renderers: [...KernelRenderers, StarsRenderer as LayerRenderer], registry: [{ svs: StarsFormat }] };
const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};
const A = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ';
/** The groups without the color scheme pickers (color-schemes.test.ts has them): what the layers list. */
const listed = (ctx: SeqViewerContext) => {
  ctx.config.set(SeqViewerConfig.UI.ColorSchemes, false);
  return legendGroups(ctx);
};

/** A ranges annotation over `a` named by its labels (`n` ranges, a category each when given). */
function ranges(b: ReturnType<typeof createSeqViewSpec>, ref: string, a: unknown, labels: string[], by = 'label') {
  const start = labels.map((_, k) => k);
  const value = { start, stop: start.map((x) => x + 1), [by]: labels };
  return b.data.annotation({ ref, geometry: 'ranges', target: a as never, value } as never);
}

describe('legend groups', () => {
  it('leave out items with no category, id or label, and a track with nothing else; labels stay (ported)', async () => {
    const b = createSeqViewSpec({ title: 'Legend' });
    const a = b.data.sequence({ ref: 'a', label: 'Chain A', value: A });
    const at = (start: number[], more = {}) => ({ start, stop: start.map((x) => x + 1), ...more });
    const unnamed = b.data.annotation({ ref: 'u', geometry: 'ranges', target: a, value: at([2, 9]) });
    const labels = { label: ['Catalytic H', 'Ligand contact'] };
    const named = b.data.annotation({ ref: 'n', geometry: 'ranges', target: a, value: at([4, 12], labels) });
    const v = b.view({ ref: 'v', label: 'Legend' });
    v.block({ sequence: a });
    v.track({ label: 'Unnamed' }).points({ data: unnamed });
    v.track({ label: 'Named' }).points({ data: named });
    const t = await opened(b.getDocument(), K);
    const groups = listed(t.ctx);
    expect(groups.map((g) => [g.track, g.items.map((i) => i.label)])).toEqual([
      ['Named', ['Catalytic H', 'Ligand contact']],
    ]);
    t.ctx.dispose();
  });

  it('a palette shorter than the values repeats; the legend lists every value with its repeated color (ported)', async () => {
    const b = createSeqViewSpec({ title: 'Repeat' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: A });
    const label = Array.from({ length: 12 }, (_, k) => `Ig ${k + 1}`);
    const d = ranges(b, 'd', a, label);
    const v = b.view({ ref: 'v', label: 'Repeat' });
    v.block({ sequence: a });
    v.track({ label: 'Domains' }).ranges({ data: d });
    const t = await opened(b.getDocument(), K);
    const [g] = listed(t.ctx);
    expect(g.items.map((i) => i.label)).toEqual(label);
    expect(g.items.map((i) => i.color).slice(9)).toEqual(['#17becf', '#1f77b4', '#ff7f0e']);
    t.ctx.dispose();
  });

  it('are one per (section, track), keyed section/trackRef: two tracks with different legends keep two groups', async () => {
    const t = await opened(starsDoc({ ranges: true }), STARS);
    const groups = listed(t.ctx);
    // Stars and Ranges are colored alike over the same data, but only Stars lists a sample size: not merged.
    expect(groups.map((g) => [g.section, g.track])).toEqual([
      ['v.main', 'Stars'],
      ['v.main', 'Ranges'],
    ]);
    expect(groups.map((g) => g.key)).toEqual(groups.map((g) => `${g.section}/${g.ref}`));
    expect(new Set(groups.map((g) => g.key)).size).toBe(2);
    expect(groups[0].ramps).toEqual(groups[1].ramps);
    t.ctx.dispose();
  });

  it('hiding one track removes its group alone: the others keep theirs, and their keys', async () => {
    const t = await opened(starsDoc({ ranges: true }), STARS);
    const [stars, ranges] = listed(t.ctx);
    // One list per assembled view: read again, the same object; a new view, a new list.
    expect(listed(t.ctx)[0]).toBe(stars);
    expect(t.ctx.actions.run('track.hide', stars.ref)).toBe(true);
    await flush();
    expect(listed(t.ctx).map((g) => g.key)).toEqual([ranges.key]);
    expect(t.ctx.actions.run('track.show', stars.ref)).toBe(true);
    await flush();
    expect(listed(t.ctx).map((g) => g.key)).toEqual([stars.key, ranges.key]);
    t.ctx.dispose();
  });

  it("list every channel of a track's layers in one group, each item once", async () => {
    const b = createSeqViewSpec({ title: 'Two layers' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: A });
    const first = ranges(b, 'r', a, ['alpha', 'beta']);
    const second = ranges(b, 'p', a, ['beta', 'gamma']);
    const v = b.view({ ref: 'v', label: 'V' });
    v.block({ sequence: a });
    const track = v.track({ ref: 'v.t', label: 'Both' });
    track.ranges({ data: first });
    track.points({ data: second });
    const t = await opened(b.getDocument(), K);
    const groups = listed(t.ctx);
    expect(groups.map((g) => g.ref)).toEqual(['v.t']);
    // `beta` has the palette's second color in the first layer and its first in the second: two different swatches.
    const labels = groups[0].items.map((i) => i.label);
    expect(labels.filter((l) => l === 'alpha' || l === 'gamma')).toEqual(['alpha', 'gamma']);
    expect(new Set(groups[0].items.map((i) => `${i.label}${i.color}`)).size).toBe(groups[0].items.length);
    t.ctx.dispose();
  });

  it('show a ramp for a continuous palette over its domain, and a sample size: stars sized and colored by affinity', async () => {
    const t = await opened(starsDoc(), STARS);
    const [g] = listed(t.ctx);
    expect(g.items).toEqual([]);
    expect(g.ramps).toHaveLength(1);
    const [ramp] = g.ramps;
    expect([ramp.stops.length, ramp.lo, ramp.hi]).toEqual([16, 6, 9]);
    // pKd 6 is the palette's bottom (black), 9 its top (white): the ramp's ends.
    expect([ramp.stops[0], ramp.stops.at(-1)]).toEqual([
      { t: 0, color: '#000000' },
      { t: 1, color: '#ffffff' },
    ]);
    expect(g.sizes).toEqual([
      { label: 'affinity (pKd)', value: 10, unit: SizeUnit.px, channel: null, field: null, color: null },
    ]);
    t.ctx.dispose();
  });

  it("show a discrete palette's bands as hard steps: two stops at each edge (pLDDT's bands)", async () => {
    const b = createSeqViewSpec({ title: 'Bands' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR' });
    const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { index: [0, 9], value: [30, 95] } });
    const palette = {
      kind: 'discrete',
      mode: 'absolute',
      colors: [
        ['#FF7D45', 0],
        ['#FFDB13', 50],
        ['#65CBF3', 70],
        ['#0053D6', 90],
      ],
    };
    b.view({ ref: 'v', label: 'V' })
      .block({ sequence: a })
      .track({ label: 'pLDDT' })
      .values({
        data: d,
        domain: [0, 100],
        color: { theme: 'field_palette', params: { field: 'value', palette } },
      } as never);
    const t = await opened(b.getDocument(), K);
    const [ramp] = listed(t.ctx)[0].ramps;
    const pct = (x: number) => Math.round(x * 1000) / 10;
    expect([ramp.lo, ramp.hi]).toEqual([0, 100]);
    expect(ramp.stops.map((s) => [pct(s.t), s.color])).toEqual([
      [0, '#ff7d45'],
      [50, '#ff7d45'],
      [50, '#ffdb13'],
      [70, '#ffdb13'],
      [70, '#65cbf3'],
      [90, '#65cbf3'],
      [90, '#0053d6'],
      [100, '#0053d6'],
    ]);
    t.ctx.dispose();
  });

  it("show a ramp for a numeric color layer over values, the domain's two ends (GenomeScale's case)", async () => {
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
    const [g] = listed(t.ctx);
    expect([g.track, g.items.length, g.ramps.map((r) => [r.lo, r.hi])]).toEqual(['Heat', 0, [[0, 20]]]);
    t.ctx.dispose();
  });

  it("EN1: layers in one domain_group list one ramp, over the union, and follow a member's data", async () => {
    const make = async (group: boolean) => {
      const b = createSeqViewSpec({ title: 'Shared' });
      const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR' });
      const column = (ref: string, hi: number) =>
        b.data.annotation({ ref, geometry: 'values', target: a, value: { index: [0, 5, 9], value: [0, hi / 2, hi] } });
      const sec = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
      const more = group ? { domain_group: 'g' } : {};
      sec.track({ label: 'Small' }).color({ data: column('small', 1), ...more });
      sec.track({ label: 'Big' }).color({ data: column('big', 8), ...more });
      return opened(b.getDocument(), K);
    };
    const ramps = (t: Awaited<ReturnType<typeof make>>) => listed(t.ctx).map((g) => g.ramps.map((r) => [r.lo, r.hi]));
    const alone = await make(false);
    expect(ramps(alone)).toEqual([[[0, 1]], [[0, 8]]]);
    alone.ctx.dispose();
    const t = await make(true);
    expect(ramps(t)).toEqual([[[0, 8]]]);
    // A member's data changes: the other layer's themes follow.
    const patch = {
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'big', params: { value: { index: [0, 9], value: [0, 16] } } }],
    };
    expect((await t.ctx.applyPatch(patch as never)).ok).toBe(true);
    expect(ramps(t)).toEqual([[[0, 16]]]);
    t.ctx.dispose();
  });

  it("have no group for tracks with nothing to list (case 9's letters and ruler), and no legend without a view", async () => {
    const t = await opened(case9(), K);
    expect(listed(t.ctx).every((g) => g.items.length + g.ramps.length + g.sizes.length > 0)).toBe(true);
    t.ctx.dispose();
    expect(legendGroups(createSeqViewer())).toEqual([]);
  });
});

describe('"+n more"', () => {
  async function many(n: number) {
    const b = createSeqViewSpec({ title: 'Many' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: A });
    const d = ranges(
      b,
      'd',
      a,
      Array.from({ length: n }, (_, k) => `cat${k}`),
    );
    const v = b.view({ ref: 'v', label: 'V' });
    v.block({ sequence: a });
    v.track({ ref: 'v.t', label: 'Many' }).ranges({ data: d });
    return opened(b.getDocument(), K);
  }
  const shown = (ctx: SeqViewerContext) => {
    const [g] = listed(ctx);
    return [g.items.length, legendShown(ctx, g).items.length, legendShown(ctx, g).more];
  };

  it('a group shows `UI.LegendItems` items and the rest as a count, until it is opened through `ui.legendOpen`', async () => {
    const t = await many(30);
    expect(shown(t.ctx)).toEqual([30, 24, 6]);
    const [g] = listed(t.ctx);
    t.ctx.ui.toggleLegend(g.key);
    expect([...t.ctx.ui.legendOpen.value]).toEqual(['v.main/v.t']);
    expect(shown(t.ctx)).toEqual([30, 30, 0]);
    t.ctx.ui.toggleLegend(g.key);
    expect(shown(t.ctx)).toEqual([30, 24, 6]);
    t.ctx.dispose();
  });

  it("a track header's tooltip lists its first six and the rest as a count, made once per track and view", async () => {
    const t = await many(30);
    const view = t.ctx.graph.assembled.value?.view;
    if (!view) throw new Error('no view');
    const track = view.resolved.tracks.ref.indexOf(view.set.strings.find('v.t'));
    const h = headerLegend(t.ctx.registry, view, track);
    expect([h?.items.map((i) => i.label), h?.more]).toEqual([['cat0', 'cat1', 'cat2', 'cat3', 'cat4', 'cat5'], 24]);
    expect(headerLegend(t.ctx.registry, view, track)).toBe(h);
    t.ctx.dispose();
  });

  it("follows the config's cap, and a group at the cap has no control", async () => {
    const t = await many(30);
    t.ctx.config.set(SeqViewerConfig.UI.LegendItems, 30);
    expect(shown(t.ctx)).toEqual([30, 30, 0]);
    t.ctx.config.set(SeqViewerConfig.UI.LegendItems, 10);
    expect(shown(t.ctx)).toEqual([30, 10, 20]);
    t.ctx.dispose();
  });

  it('an opened group is dropped from `ui.legendOpen` when its section leaves the document, after the update', async () => {
    const t = await many(30);
    t.ctx.ui.toggleLegend('v.main/v.t');
    t.ctx.ui.toggleLegend('gone.main/x');
    expect(t.ctx.ui.legendOpen.value.size).toBe(2);
    await t.ctx.setDocument(case9());
    await flush();
    expect([...t.ctx.ui.legendOpen.value]).toEqual([]);
    t.ctx.dispose();
  });
});

/** A track menu's owner, as both UIs keep it: the submenu's layer, and the menu opened again. */
function owner() {
  const d = { layer: null as string | null, opened: 0 } as { layer: string | null; opened: number } & ColorsDrill;
  d.open = (l) => {
    d.layer = l;
    d.opened++;
  };
  return d;
}

describe('merged groups and the color schemes', () => {
  it("25 tracks of three kinds give three groups, each titled by its tracks' common name and count", async () => {
    const t = await opened(kindsDoc(), K);
    const groups = legendGroups(t.ctx);
    expect(groups.map((g) => [g.track, g.refs.length, g.items.map((i) => i.label), g.colors?.layers.length])).toEqual([
      ['Genes (5)', 5, ['gene'], 5],
      ['Transcripts (10)', 10, ['UTR', 'exon'], 10],
      ['Domains (10)', 10, ['domain'], 10],
    ]);
    expect(groups.map((g) => g.key)).toEqual(['v.main/gene1', 'v.main/tx1.1', 'v.main/dom1.1']);
    // Read again, the same list: computed once per assembled view.
    expect(legendGroups(t.ctx)).toBe(groups);
    t.ctx.dispose();
  });

  it('titles: common leading words, else trailing ones, else the layer type; a single track its label', () => {
    expect(groupTitle(['Gene G1', 'Gene G2'], 'ranges')).toBe('Genes (2)');
    expect(groupTitle(['G1.1 domain', 'G2.1 domain'], 'ranges')).toBe('Domains (2)');
    expect(groupTitle(['Alpha', 'Beta', 'Gamma'], 'ranges')).toBe('Ranges (3)');
    expect(groupTitle(['Members', 'Members'], 'rows')).toBe('Members (2)');
    expect(groupTitle(['SH2 1', 'SH2 2'], 'ranges')).toBe('SH2 (2)');
    expect(groupTitle(['Heat'], 'color')).toBe('Heat');
  });

  it("a track with nothing to list has no group, and its menu offers Colors: a submenu of its layer's schemes", async () => {
    const t = await opened(kindsDoc(true), K);
    expect(legendGroups(t.ctx).map((g) => g.track)).toEqual(['Transcripts (10)', 'Domains (10)']);
    const d = owner();
    const labels = () => trackMenu(t.ctx, 'v.main', 'gene2', d).map((i) => i.label);
    expect(labels()).toContain('Colors ▸');
    // Without an owner to keep the submenu (a host calling the helper alone), no Colors item.
    expect(trackMenu(t.ctx, 'v.main', 'gene2').map((i) => i.label)).not.toContain('Colors ▸');
    trackMenu(t.ctx, 'v.main', 'gene2', d)
      .find((i) => i.label === 'Colors ▸')
      ?.run();
    expect([d.layer, d.opened]).toEqual(['l.gene2', 1]);
    const sub = trackMenu(t.ctx, 'v.main', 'gene2', d);
    expect(sub[0].label).toBe('‹ Colors');
    expect(sub.slice(1).map((i) => [i.label, i.checked])).toEqual([
      ['document', true],
      ['Category 10', false],
      ['Okabe–Ito (colorblind-safe)', false],
      ['Tableau 10', false],
    ]);
    sub.find((i) => i.key === 'scheme:okabe_ito')?.run();
    expect([...t.ctx.ui.colorSchemes.value]).toEqual([['l.gene2', 'okabe_ito']]);
    await t.ctx.graph.update();
    expect(trackMenu(t.ctx, 'v.main', 'gene2', d).find((i) => i.checked)?.key).toBe('scheme:okabe_ito');
    sub[0].run();
    expect(d.layer).toBe(null);
    // `UI.ColorSchemes` off: no Colors item, and no picker in the legend.
    t.ctx.config.set(SeqViewerConfig.UI.ColorSchemes, false);
    expect(labels()).not.toContain('Colors ▸');
    t.ctx.config.set(SeqViewerConfig.UI.TrackActions, ['hide']);
    t.ctx.config.set(SeqViewerConfig.UI.ColorSchemes, true);
    expect(labels()).toEqual(['Hide track']);
    t.ctx.dispose();
  });

  it("a merged group's switch applies to every track in it; one track's own pick splits it off", async () => {
    const t = await opened(kindsDoc(), K);
    const tx = legendGroups(t.ctx)[1];
    const colors = tx.colors as NonNullable<typeof tx.colors>;
    schemeItems(t.ctx, colors.layers, colors)
      .find((i) => i.key === 'scheme:tableau10')
      ?.run();
    expect(t.ctx.ui.colorSchemes.value.size).toBe(10);
    await t.ctx.graph.update();
    const after = legendGroups(t.ctx);
    expect(after.map((g) => [g.track, g.colors?.value, g.colors?.label])).toEqual([
      ['Genes (5)', null, null],
      ['Transcripts (10)', 'tableau10', 'Tableau 10'],
      ['Domains (10)', null, null],
    ]);
    expect(after[1].key).toBe(tx.key);
    // From one transcript's menu, its own scheme: that track has its own group, in its first track's order.
    const d = owner();
    d.open('l.tx2.1');
    trackMenu(t.ctx, 'v.main', 'tx2.1', d)
      .find((i) => i.key === 'scheme:okabe_ito')
      ?.run();
    await t.ctx.graph.update();
    expect(legendGroups(t.ctx).map((g) => [g.track, g.refs.length])).toEqual([
      ['Genes (5)', 5],
      ['Transcripts (9)', 9],
      ['Domains (10)', 10],
      ['Transcript G2.1', 1],
    ]);
    // `UI.ColorSchemes` off: no picker, so the groups merge by their content alone.
    t.ctx.config.set(SeqViewerConfig.UI.ColorSchemes, false);
    await t.ctx.setDocument(kindsDoc());
    expect(legendGroups(t.ctx).map((g) => [g.track, g.colors])).toEqual([
      ['Genes (5)', null],
      ['Transcripts (10)', null],
      ['Domains (10)', null],
    ]);
    t.ctx.dispose();
  });

  it('the legend shows `UI.LegendGroups` groups and the rest as "+n more", until opened through `ui.legendOpen`', async () => {
    const t = await opened(distinctDoc(12), K);
    const all = legendGroups(t.ctx);
    const shown = () => {
      const s = legendShownGroups(t.ctx, all);
      return [s.groups.length, s.more, s.toggle];
    };
    expect(all).toHaveLength(12);
    expect(shown()).toEqual([8, 4, true]);
    const { key } = legendShownGroups(t.ctx, all);
    expect(key).toBe('v.main/*');
    t.ctx.ui.toggleLegend(key);
    expect(shown()).toEqual([12, 0, true]);
    // The key is the section's: it stays while the section does.
    await flush();
    expect(t.ctx.ui.legendOpen.value.has(key)).toBe(true);
    t.ctx.ui.toggleLegend(key);
    t.ctx.config.set(SeqViewerConfig.UI.LegendGroups, 12);
    expect(shown()).toEqual([12, 0, false]);
    t.ctx.dispose();
  });
});

describe("a values layer's y scale", () => {
  it("starts at the group's given domain, not at 0; at 0 for the data's own span", async () => {
    const make = async (given: boolean) => {
      const b = createSeqViewSpec({ title: 'Axis' });
      const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR' });
      const column = (ref: string, lo: number) =>
        b.data.annotation({ ref, geometry: 'values', target: a, value: { index: [0, 9], value: [lo, lo + 10] } });
      const sec = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
      const g = { domain_group: 'g' };
      sec
        .track({ label: 'One' })
        .values({ data: column('one', 6), ...g, ...(given ? { domain: [5, 30] } : {}) } as never);
      sec.track({ label: 'Two' }).values({ data: column('two', 8), ...g } as never);
      const t = await opened(b.getDocument(), K);
      const view = t.ctx.graph.assembled.value?.view;
      const lows = view?.prepared.map((p) => (p.state as { lo: number } | null)?.lo);
      t.ctx.dispose();
      return lows;
    };
    expect(await make(true)).toEqual([5, 5]);
    expect(await make(false)).toEqual([0, 0]);
  });
});

describe('size legends (EN2b)', () => {
  const scale = (field: string, range: [number, number], more: Record<string, unknown> = {}) => ({
    theme: 'field_scale',
    params: { field, range, ...more },
  });
  function doc() {
    const b = createSeqViewSpec({ title: 'Sizes' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: A });
    const count = [1, 4, 9];
    const pairs = b.data.annotation({
      ref: 'pp',
      geometry: 'pairs',
      target: a,
      value: { start: [0, 2, 4], b_target: count.map(() => 'a'), b_start: [20, 22, 24], count },
    } as never);
    const items = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: a,
      value: { start: [0, 10], stop: [8, 18], kind: ['weak', 'strong'], score: [2, 8] },
    } as never);
    const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
    v.track({ label: 'Contacts' }).pairs({
      data: pairs,
      width: scale('count', [1, 8], { scale: 'sqrt' }),
      color: {
        theme: 'field_palette',
        params: { field: 'count', palette: { kind: 'continuous', colors: ['#000000', '#ffffff'] } },
      },
    } as never);
    v.track({ label: 'Kinds' }).ranges({
      data: items,
      thickness: scale('kind', [0, 0], { map: { weak: 0.3, strong: 1 } }),
    } as never);
    v.track({ label: 'Ruled' }).points({
      data: items,
      radius: {
        theme: 'rules',
        params: { rules: [{ when: [{ measure: { field: 'score' }, min: 5 }], size: 6, label: 'high' }] },
      },
    } as never);
    v.track({ label: 'Uniform' }).points({ data: items, radius: { theme: 'uniform', params: { size: 4 } } } as never);
    return b.getDocument();
  }
  it("a field_scale lists its domain's ends and middle through its scale, colored by a ramp over the same field", async () => {
    const t = await opened(doc(), K);
    const groups = listed(t.ctx);
    const contacts = groups.find((g) => g.track === 'Contacts');
    expect(contacts?.ramps).toHaveLength(1);
    const sizes = contacts?.sizes ?? [];
    expect(sizes.map((x) => [x.label, x.channel, x.field])).toEqual([
      ['1', 'width', 'count'],
      ['5', 'width', 'count'],
      ['9', 'width', 'count'],
    ]);
    // sqrt: the middle of the domain is sqrt(0.5) of the range.
    expect(sizes.map((x) => Number(x.value.toFixed(3)))).toEqual([1, Number((1 + Math.SQRT1_2 * 7).toFixed(3)), 8]);
    expect(sizes.every((x) => x.unit === SizeUnit.px)).toBe(true);
    // Merged with the ramp over count (black to white): the ends' colors.
    expect([sizes[0].color, sizes[2].color]).toEqual(['#000000', '#ffffff']);
    t.ctx.dispose();
  });
  it("a string field's map lists its entries; a rule its label; a uniform size nothing", async () => {
    const t = await opened(doc(), K);
    const groups = listed(t.ctx);
    const kinds = groups.find((g) => g.track === 'Kinds')?.sizes ?? [];
    expect(kinds.map((x) => [x.label, x.value, x.unit, x.channel, x.color])).toEqual([
      ['weak', 0.3, SizeUnit.band, 'thickness', null],
      ['strong', 1, SizeUnit.band, 'thickness', null],
    ]);
    const ruled = groups.find((g) => g.track === 'Ruled')?.sizes ?? [];
    expect(ruled.map((x) => [x.label, x.value, x.channel, x.field])).toEqual([['high', 6, 'radius', null]]);
    expect(groups.some((g) => g.track === 'Uniform')).toBe(false);
    t.ctx.dispose();
  });
});
