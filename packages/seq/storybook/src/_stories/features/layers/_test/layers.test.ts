/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  DarkTokens,
  hitView,
  type RangesLay,
  recordDrawLists,
  type RulerLay,
  type SecondaryLay,
  type VariantsLay,
} from '@molstar/seqx-core/render';
import { LayerViews } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { GeneFile, Secondary, Variants } from '@molstar/seqx-examples/p53/data';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// The small layers' stories, headless on a fake clock, each composed as its story composes it (a viewer of the view and
// its SVG export): the canvas and the SVG draw the same lists in both themes, and what each story says to look at is
// there. The layers' own tests are `core/src/render/_test`.

const Protein = '/views/main/protein';
const Gene = '/views/main/gene';

let viewer: SeqViewer | undefined;
let gene = '';

beforeAll(async () => {
  gene = await nodeFixtureAssets().reader.readText(GeneFile, new AbortController().signal);
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open(id: string): Promise<SeqViewer> {
  viewer = await openExample(LayerViews.find((v) => v.id === id)!, { gene, delayMs: 0 });
  return viewer;
}

/** The viewer's SVG export, its asynchronous work let go on the fake clock. */
async function svgOf(v: SeqViewer) {
  const p = v.renderSvg({ idPrefix: 'x-' });
  await vi.advanceTimersByTimeAsync(400);
  return p;
}

const tracks = (v: SeqViewer, section: number) => v.frame!.layout.sections[section].tracks;

describe('the small layers’ stories, headless', () => {
  for (const view of LayerViews) {
    it(`${view.title}: the canvas and the SVG draw the same lists, in the light and the dark theme`, async () => {
      const v = await open(view.id);
      expect(v.layoutIssues.value).toEqual([]);
      // the rails off: an export has none, nor their row under the ruler
      v.setRailsShown(false);
      v.draw();
      const light = await svgOf(v);
      expect(light.issues).toEqual([]);
      for (const pass of ['data', 'labels'] as const)
        expect(recordDrawLists(light.lists, [pass])).toBe(recordDrawLists(v.frame!.lists, [pass]));
      // the theme switch: the viewer takes the dark table, its export follows
      v.setTokens(DarkTokens);
      v.draw();
      const dark = await svgOf(v);
      expect(dark.issues).toEqual([]);
      // the same marks in other colors
      expect(recordDrawLists(dark.lists, ['data'])).not.toBe(recordDrawLists(light.lists, ['data']));
      expect(dark.lists.passes.data.size).toBe(light.lists.passes.data.size);
    });
  }

  it('Ruler: the protein’s labels rise, the gene’s fall, and each label is the position under its tick', async () => {
    const v = await open('ruler');
    const labels = (section: number) => {
      const lay = tracks(v, section)[0].layers[0].lay as RulerLay;
      const out: { x: number; label: string }[] = [];
      for (let j = 0; j < lay.n; j++) if (lay.shown[j]) out.push({ x: lay.x[j], label: lay.labels[j] });
      return { lay, out: out.sort((a, b) => a.x - b.x) };
    };
    const protein = labels(0);
    const geneLabels = labels(1);
    expect(protein.out.length).toBeGreaterThan(2);
    expect(geneLabels.out.length).toBeGreaterThan(2);
    const nums = (o: typeof protein.out) => o.map((l) => Number(l.label));
    expect(nums(protein.out)).toEqual([...nums(protein.out)].sort((a, b) => a - b));
    expect(nums(geneLabels.out)).toEqual([...nums(geneLabels.out)].sort((a, b) => b - a));
    for (const { lay, out } of [protein, geneLabels])
      for (const l of out) expect(String(lay.axis.positionAt(l.x) + 1)).toBe(l.label);
  });

  it('Secondary structure: each strand has its arrowhead, on the protein and on the gene; each mark is hit', async () => {
    const v = await open('secondary');
    v.view.setRange(Protein, 90, 240);
    v.draw();
    for (const section of [0, 1]) {
      const track = tracks(v, section).find((t) => t.layers[0].type.name === 'secondary')!;
      const lay = track.layers[0].lay as SecondaryLay;
      const inView = Secondary.id.filter((_, i) => Secondary.end[i] > 90 && Secondary.start[i] < 240);
      const keys = new Set<string | number>();
      for (let j = 0; j < lay.m; j++) keys.add(lay.rows.keys[lay.mitem[j]]);
      expect([...keys].sort()).toEqual([...inView].sort());
      for (let j = 0; j < lay.m; j++) {
        const x = (lay.mx0[j] + lay.mx1[j]) / 2;
        expect(hitView(v.frame!.layout, x, track.box.y + track.box.height / 2)?.element?.key).toBe(
          lay.rows.keys[lay.mitem[j]],
        );
      }
      // the strands wholly in view point right on both sections (the gene reads along the protein)
      for (let j = 0; j < lay.m; j++) {
        const i = Secondary.id.indexOf(String(lay.rows.keys[lay.mitem[j]]));
        const whole = Secondary.start[i] >= 90 && Secondary.end[i] <= 240;
        if (Secondary.category[i] === 'strand' && whole) expect(lay.mtip[j]).toBe(1);
        else expect(lay.mtip[j]).toBe(0);
      }
    }
  });

  it('Variants: each variant in view has a head hit as itself, those of one site stacked; zoomed out, a strip', async () => {
    const v = await open('variants');
    v.view.setRange(Protein, 160, 320);
    v.draw();
    const track = tracks(v, 0).find((t) => t.layers[0].type.name === 'variants')!;
    const lay = track.layers[0].lay as VariantsLay;
    const inView = Variants.id.filter((_, i) => Variants.start[i] >= 160 && Variants.start[i] < 320);
    const heads: { key: string | number; x: number; y: number }[] = [];
    for (let j = 0; j < lay.m; j++) heads.push({ key: lay.rows.keys[lay.mrow[j]], x: lay.mx[j], y: lay.my[j] });
    expect(heads.map((h) => h.key).sort()).toEqual([...inView].sort());
    for (const h of heads) expect(hitView(v.frame!.layout, h.x, h.y)?.element?.key).toBe(h.key);
    const r248 = heads.filter((h) => h.key === 'R248Q' || h.key === 'R248W');
    expect(new Set(r248.map((h) => h.x)).size).toBe(1);
    expect(new Set(r248.map((h) => h.y)).size).toBe(2);
    // far out (the gene’s window follows): the heads give way to the count strip
    v.view.setRange(Gene, 0, 1500);
    v.draw();
    const gtrack = tracks(v, 1).find((t) => t.layers[0].type.name === 'variants')!;
    expect((gtrack.layers[0].lay as VariantsLay).density).toBeDefined();
  });

  it('Values as a line and an area: three tracks of one signal; the line and the area break at the missing residues', async () => {
    const v = await open('values-line');
    const t = tracks(v, 0);
    expect(t.map((x) => x.track.label)).toEqual(['', 'Bars', 'Line', 'Area']);
    const svg = await svgOf(v);
    // the line’s two paths and the area’s two polygons and two paths, as stroked and filled paths
    expect(svg.svg.match(/<path /g)!.length).toBeGreaterThanOrEqual(6);
  });

  it('Ranges as lines and arrows: each style drawn, the dash anchored to its domain over a pan', async () => {
    const v = await open('ranges-lines');
    const lays = tracks(v, 0)
      .filter((t) => t.layers[0].type.name === 'ranges')
      .map((t) => t.layers[0].lay as RangesLay);
    expect(lays.map((l) => l.mode)).toEqual(lays.map(() => 'items'));
    expect(lays.map((l) => l.style)).toEqual([1, 2, 3, 4, 5, 6, 7, 0, 1]);
    expect(lays[7].arrow).toBe(1);
    expect(lays[8].arrow).toBe(3);
    // panned and zoomed in, the dashed and dotted domains are still drawn as dashed paths (the phase is the core's test)
    for (const [from, to] of [
      [150, 250],
      [200, 300],
    ]) {
      v.view.setRange(Protein, from, to);
      v.draw();
      const dashes = new Set<number>();
      const paths = v.frame!.lists.passes.data.paths;
      for (let i = 0; i < paths.n; i++) if (paths.dash[i] > 0) dashes.add(paths.dash[i]);
      expect([...dashes].sort()).toEqual([2, 6]);
    }
  });
  it('Tracks in their wells: the padding inside each well, the wells squircles their data is clipped to', async () => {
    const v = await open('wells');
    const [letters, domains, values] = tracks(v, 0);
    expect(letters.box.y - letters.well.y).toBe(8);
    expect(letters.well.y + letters.well.height).toBe(letters.box.y + letters.box.height);
    expect(domains.box.y - domains.well.y).toBe(6);
    expect(domains.well.height - domains.box.height).toBe(12);
    expect(values.well).toEqual(values.box);
    // the wells and their clips squircles; nothing hit in the padding, the data hit below it
    const background = recordDrawLists(v.frame!.lists, ['background']).split('\n').slice(3, 6);
    expect(background.every((l) => / squircle [\d.]+$/.test(l))).toBe(true);
    const clips = recordDrawLists(v.frame!.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('clip '));
    expect(clips.length).toBe(3);
    expect(clips.every((l) => / squircle [\d.]+$/.test(l))).toBe(true);
    expect(hitView(v.frame!.layout, letters.box.x + 20, letters.well.y + 4)).toBeUndefined();
    expect(hitView(v.frame!.layout, letters.box.x + 20, letters.box.y + 4)?.position).toBeDefined();
    const svg = await svgOf(v);
    expect(svg.svg.match(/<clipPath id="[^"]*"><path /g)?.length).toBe(3);
  });
});
