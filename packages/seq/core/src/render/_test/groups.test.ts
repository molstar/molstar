/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's track groups in the tracks layout (_spec/README.md "Render", Track groups): an open group's header line,
// its tracks indented, the wider gaps; a folded group's one row (its summary in each layer kind's form: values strips,
// the ranges family in one lane), an authored summary shown only while folded (EC-103), `none` the header alone; the
// header's hit; text drawing every group open; a folded member's layers laid only in their summary form.

import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { groupsSpec } from '../../view/_test/groups-spec.js';
import { BuiltInLayerTypes, BuiltInLayouts } from '../built-ins.js';
import { DefaultRenderConfig } from '../config.js';
import { exportSvg } from '../export.js';
import { bundledTextMetrics, LightTokens, layoutView, type ViewLayout } from '../index.js';
import { lanesHeight } from '../sizing.js';
import { settle } from './helpers.js';

const Size = { width: 640, height: 600 };
const Section = '/views/v/p53';

async function lay(spec: Node, groups?: readonly string[]): Promise<ViewLayout> {
  const s = await settle(spec);
  const viewports = groups ? { [Section]: { start: 0, end: 53, groups } } : undefined;
  // the rows loaded: a layout alone reads only what the store holds
  await s.frame({ size: Size, ...(viewports ? { viewports } : {}) });
  return layoutView({
    model: s.model,
    view: s.view,
    size: Size,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    data: s.store,
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    ...(viewports ? { viewports } : {}),
  });
}

const bands = (l: ViewLayout) => l.sections[0].tracks.map((t) => t.track.address.slice(Section.length + 1));
const config = DefaultRenderConfig;

describe('track groups in the tracks layout', () => {
  it('opens with each group’s `collapsed`: the header lines, indented labels, wider gaps; a folded group’s summary bands', async () => {
    const l = await lay(groupsSpec());
    // Variants folded at opening with an authored summary: only it is drawn there (EC-103)
    expect(bands(l)).toEqual([
      'letters',
      'features/domains',
      'features/motifs',
      'variants/density',
      'signals/a',
      'signals/b',
    ]);
    const s = l.sections[0];
    expect(s.groups!.map((g) => [g.group.address.slice(Section.length + 1), g.folded, g.first, g.last])).toEqual([
      ['features', false, 1, 2],
      ['variants', true, 3, 3],
      ['signals', false, 4, 5],
    ]);
    const [letters, domains, motifs, density] = s.tracks;
    // the header line above Features' first well, after the group's wider gap
    const features = s.groups![0];
    expect(features.box.y).toBe(domains.well.y - config.header.groupPx);
    expect(features.box.y - (letters.well.y + letters.well.height)).toBe(config.groupGap);
    // within a group: `trackGap`; a group's tracks indented; a folded group's band leaves its label to the header
    expect(motifs.well.y - (domains.well.y + domains.well.height)).toBe(config.trackGap);
    expect([letters.indent, domains.indent, motifs.indent]).toEqual([undefined, 12, 12]);
    expect(density.labelText).toBe('');
    expect(s.groups![1].count).toBe('2 tracks');
    expect(s.groups![1].label).toBe('Variants');
  });

  it('EC-103: a group’s authored summary track is drawn only while it is folded', async () => {
    const open = await lay(groupsSpec(), []);
    expect(bands(open)).toEqual([
      'letters',
      'features/domains',
      'features/motifs',
      'variants/clinvar',
      'variants/gnomad',
      'signals/a',
      'signals/b',
    ]);
    const folded = await lay(groupsSpec(), ['variants']);
    expect(bands(folded)).toContain('variants/density');
    expect(bands(folded)).not.toContain('variants/clinvar');
  });

  it('`auto`: values as a strip a member (under the lane, none when no member draws in one); the ranges family in one shared lane; a folded member lays no open form', async () => {
    const l = await lay(groupsSpec({ variants: { summary: 'auto' } }), ['features', 'variants', 'signals']);
    expect(bands(l)).toEqual([
      'letters',
      'features#summary',
      'variants#summary',
      'variants#strip0',
      // values alone: their strips with no empty lane above them
      'signals#strip0',
      'signals#strip1',
    ]);
    const s = l.sections[0];
    const [, featuresLane, variantsLane, densityStrip, a, b] = s.tracks;
    // the lane: every lane-form layer of the group (domains and motifs), keyed `#summary`
    expect(featuresLane.track.layers.map((x) => x.address)).toEqual([
      '/views/v/p53/features/domains/l#summary',
      '/views/v/p53/features/motifs/l#summary',
    ]);
    expect(variantsLane.track.layers.map((x) => x.type)).toEqual(['variants', 'variants']);
    // strips: `stripPx` high, 1 px apart
    expect(densityStrip.box.height).toBe(config.group.summary.stripPx);
    expect(b.well.y - (a.well.y + a.well.height)).toBe(1);
    expect(s.groups!.map((g) => g.count)).toEqual(['2 tracks', '3 tracks', '2 tracks']);
    // no open form of a folded member laid
    expect(s.lays.flat().some((x) => x.key === '/views/v/p53/features/domains/l')).toBe(false);
    // the folded row's header box spans its bands, padded at its top to the header's line (strips alone)
    const sig = s.groups![2];
    expect(sig.box.height).toBe(config.header.groupPx);
    expect(sig.box.y + sig.box.height).toBe(b.well.y + b.well.height);
  });

  it('`lanes: one`: a ranges track in one lane with its labels off, folded (the summary lane) and open', async () => {
    const lanesOf = (l: ViewLayout, track: string) => {
      const t = l.sections[0].tracks.find((x) => x.track.address === `${Section}/${track}`)!;
      return t.layers[0].lay as unknown as { laneCount: number; one: boolean; box: { height: number } };
    };
    // open, packed: its labels inside
    const packed = await lay(groupsSpec(), []);
    expect(lanesOf(packed, 'features/domains').one).toBe(false);
    // open, `lanes: one` on the track: one lane, no labels, the track one lane high
    const spec = groupsSpec();
    const domains = spec.children![1].children![0].children![0].children![1].children![0].children![0];
    domains.params = { ...domains.params, lanes: 'one' };
    const one = await lay(spec, []);
    const d = lanesOf(one, 'features/domains');
    expect([d.laneCount, d.one]).toEqual([1, true]);
    expect(one.sections[0].tracks[1].box.height).toBe(lanesHeight(1, config));
    // its labels drawn packed, none in one lane
    const svg = async (n: Node) =>
      (
        await exportSvg(n, {
          size: { width: 640, height: 'auto' },
          groups: { [Section]: [] },
          types: BuiltInLayerTypes,
          layouts: BuiltInLayouts,
        })
      ).svg;
    expect(await svg(groupsSpec())).toContain('>PRD<');
    expect(await svg(spec)).not.toContain('>PRD<');
    // folded: the summary lane lays every member in one lane, no labels
    const folded = await lay(groupsSpec(), ['features']);
    const lane = folded.sections[0].tracks[1];
    for (const l of lane.layers) {
      const r = l.lay as unknown as { laneCount: number; one: boolean };
      expect([r.laneCount, r.one]).toEqual([1, true]);
    }
  });

  it('`none`: a folded group draws its header line alone', async () => {
    const l = await lay(groupsSpec({ variants: { summary: 'none' } }), ['variants']);
    expect(bands(l)).toContain('variants#header');
    const t = l.sections[0].tracks.find((x) => x.track.address.endsWith('#header'))!;
    expect(t.box.height).toBe(config.header.groupPx);
    expect(t.track.layers).toEqual([]);
  });

  it('the header’s hit: `group-header`, the group’s address; the label column only', async () => {
    const l = await lay(groupsSpec(), []);
    const s = l.sections[0];
    const h = s.groups![0];
    const hit = s.layout.hitChrome!(s, h.box.x + 20, h.box.y + h.box.height / 2, undefined);
    expect(hit).toMatchObject({ target: 'group-header', layer: `${Section}/features`, section: Section });
    expect(s.layout.hitChrome!(s, 400, h.box.y + 2, undefined)).toBeUndefined();
    // folded: its row's label column is the group's header, not the summary row's label
    const f = (await lay(groupsSpec(), ['features'])).sections[0];
    const fh = f.groups![0];
    expect(f.layout.hitChrome!(f, fh.box.x + 20, fh.box.y + fh.box.height / 2, undefined)).toMatchObject({
      target: 'group-header',
      layer: `${Section}/features`,
    });
  });

  it('a text section draws every group open with no header line, said once when some are folded', async () => {
    const spec = groupsSpec();
    const sec = spec.children![1].children![0].children![0];
    sec.params = { ...sec.params, layout: 'text' };
    const s = await settle(spec);
    const l = layoutView({
      model: s.model,
      view: s.view,
      size: Size,
      tokens: LightTokens,
      text: bundledTextMetrics(),
      data: s.store,
      types: BuiltInLayerTypes,
      layouts: BuiltInLayouts,
    });
    expect(l.sections[0].groups).toBeUndefined();
    expect(bands(l)).not.toContain('variants/density');
    expect(l.issues.map((i) => i.code)).toContain('seq.render.groups.layout');
  });

  it('a scrolled section: a group header outside its scrolled area is neither drawn nor hit', async () => {
    const s = await settle(groupsSpec(), [], { section: { maxPx: 120 } });
    const frame = (top: number) =>
      s.frame({ size: Size, settled: true, viewports: { [Section]: { start: 0, end: 53, top, groups: [] } } });
    const at0 = (await frame(0)).layout.sections[0];
    const features = at0.groups!.find((g) => g.group.address.endsWith('/features'))!;
    expect(at0.placement.clip).toBe(true);
    const f = await frame(200);
    const sec = f.layout.sections[0];
    const p = sec.placement;
    // every header left is wholly in the scrolled area; Features' scrolled above it is gone
    for (const g of sec.groups ?? []) {
      expect(g.box.y).toBeGreaterThanOrEqual(p.scrollY!);
      expect(g.box.y + g.box.height).toBeLessThanOrEqual(p.scrollEnd || p.box.y + p.box.height);
    }
    expect(sec.groups?.some((g) => g.group.address.endsWith('/features')) ?? false).toBe(false);
    const hit = sec.layout.hitChrome!(sec, features.box.x + 20, p.scrollY! + 2, undefined);
    expect(hit?.target).not.toBe('group-header');
  });

  it('an export draws the folds it is given (`groups`), else each group’s `collapsed`', async () => {
    const svg = (groups?: Record<string, readonly string[]>) =>
      exportSvg(groupsSpec(), { size: { width: 640, height: 'auto' }, ...(groups ? { groups } : {}) });
    const opening = await svg();
    const open = await svg({ [Section]: [] });
    expect(open.svg).toContain('ClinVar');
    expect(opening.svg).not.toContain('ClinVar');
    expect(opening.svg).toContain('Variants');
    // the chevron is chrome: a figure draws none (the canvas does)
    expect(opening.svg).not.toMatch(/[▸▾]/);
    expect((await lay(groupsSpec(), [])).sections[0].groups!.every((g) => g.chevron)).toBe(true);
  });
});
