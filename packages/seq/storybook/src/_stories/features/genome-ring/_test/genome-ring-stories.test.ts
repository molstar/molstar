/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PassWriter } from '@molstar/seqx-core/extension';
import { DarkTokens, DrawLists, LightTokens, viewLegends } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import type { ExampleView } from '@molstar/seqx-examples';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { layerAria } from '@molstar/seqx-viewer/layer-aria';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { describerFor } from '@molstar/seqx-viewer/describe';
import { emphasisHeld, Kind, RibbonsType } from '@molstar/seqx-ext-rondo';
import { hoverStrongestStem, ribbonsOf, ringStoryView, strongestStem } from '../views.js';
import { bandsLanded } from './ring-helpers.js';

// Features/Genome Ring, headless on a fake clock as the stories compose them (a viewer with `rondoViewer`, or none for
// the fallback's, its SVG export light and dark): a twin a default of Design › Genome ring, chr4's ring stepped through
// its bins landing and then its items.

let viewer: SeqViewer | undefined;

// the genome chunk landed first: its import is real I/O, never on the fake clock
beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const Rondo = [rondoViewer] as const;

async function open(view: ExampleView, extensions: typeof Rondo | 'none' = Rondo, wait = 3000) {
  viewer?.retire();
  viewer = await openExample(view, {
    delayMs: 0,
    wait,
    fixtures: true,
    ...(extensions !== 'none' ? { extensions } : {}),
  });
  viewer.draw();
  return viewer;
}

/** The view's SVG export, its clock run (a promise settled by the timers the export waits on). */
async function svgOf(v: SeqViewer, dark: boolean, navigators?: 'include') {
  const r = v.renderSvg({
    idPrefix: 'r-',
    tokens: dark ? DarkTokens : LightTokens,
    ...(navigators ? { navigators } : {}),
  });
  await vi.advanceTimersByTimeAsync(400);
  return r;
}

/** The view's SVG export light and dark (the story's SVG follows the toolbar's theme): settled, no error, two pictures. */
async function exported(v: SeqViewer) {
  const light = await svgOf(v, false);
  const dark = await svgOf(v, true);
  expect(light.settled).toBe(true);
  expect(light.issues.filter((i) => i.severity === 'error')).toEqual([]);
  expect(dark.svg).not.toBe(light.svg);
  return { light, dark };
}

const errors = (v: SeqViewer) => [
  ...v.store.issues.value.filter((i) => i.severity === 'error'),
  ...v.layoutIssues.value.filter((i) => i.severity === 'error'),
];

const descriptions = (v: SeqViewer) => layerAria(v.frame!.layout).map((a) => a.description);

describe('Genome ring stories, headless', () => {
  it('1. the genome at the fit: 23 chromosomes’ stems, ribbons wide by count; light and dark', async () => {
    const v = await open(ringStoryView('ring'));
    const l = ribbonsOf(v)!;
    expect(errors(v)).toEqual([]);
    expect(l.items).toBe(false);
    expect(new Set(Array.from(l.stemSystem.subarray(0, l.ns))).size).toBe(23);
    expect(l.nr).toBeGreaterThan(100);
    expect(descriptions(v).some((d) => /ribbons, 23 chromosomes; strongest/.test(d))).toBe(true);
    const { light } = await exported(v);
    expect(light.issues.filter((i) => i.code === 'seq.view.layer.extension')).toEqual([]);
  });

  it('2. EC-235: the same spec without the extension: core’s chords, the issue once with the version said', async () => {
    const v = await open(ringStoryView('ring'), 'none');
    expect(ribbonsOf(v)).toBeUndefined();
    const types = v.frame!.layout.sections.flatMap((s) => s.sectionLays.map((x) => x.type.name));
    expect(types).toContain('pairs');
    const { light } = await exported(v);
    const said = light.issues.filter((i) => i.code === 'seq.view.layer.extension');
    expect(said.map((i) => i.message)).toEqual(['1 layer(s) need `rondo` 0.1: drawn as `pairs`']);
  });

  it('3. EC-169: chr4’s ring on a fake clock: the hatch, the bins landing as ribbons, then its links as items', async () => {
    viewer = await openExample(ringStoryView('ring-chromosome'), {
      delayMs: 150,
      wait: 0,
      fixtures: true,
      extensions: Rondo,
    });
    const seen: string[] = [];
    for (let t = 0; t < 80 && seen.at(-1) !== 'items'; t++) {
      viewer.draw();
      const l = ribbonsOf(viewer);
      const now = !l ? 'absent' : l.empty ? l.empty : l.items ? 'items' : 'bins';
      if (seen.at(-1) !== now) seen.push(now);
      await vi.advanceTimersByTimeAsync(50);
    }
    // pending (never "none" while a tile is wanted), then bins, then items: nothing else between
    expect(seen).toEqual(['pending', 'bins', 'items']);
    const l = ribbonsOf(viewer)!;
    expect(l.nr).toBeGreaterThan(0);
    // its outward links as stubs inward, labelled with their partners
    expect(l.nb).toBeGreaterThan(0);
    expect(l.stubLabel.slice(0, l.nb).some((x) => /^chr\w+/.test(x))).toBe(true);
    await exported(viewer);
  });

  it('4. a pair (chr4 × chr7): its ribbons between the two, the rest stubs', async () => {
    const v = await open(ringStoryView('ring-pair'));
    const l = ribbonsOf(v)!;
    expect(errors(v)).toEqual([]);
    const shown = new Set<string>();
    for (let s = 0; s < l.ns; s++) if (l.stemPlaced[s]) shown.add(l.names[l.stemSystem[s]]);
    expect([...shown].sort()).toEqual(['chr4', 'chr7']);
    expect(l.nr).toBeGreaterThan(0);
    expect(l.nb).toBeGreaterThan(0);
    await exported(v);
  });

  it('5. the square beside: the export with the navigators draws more than without', async () => {
    const v = await open(ringStoryView('ring'));
    const without = await svgOf(v, false);
    const svg = await svgOf(v, false, 'include');
    expect(svg.settled).toBe(true);
    expect(svg.svg.length).toBeGreaterThan(without.svg.length);
  });

  it('6. hovering a stem: the stem the hover, its ribbons lit in the hover pass', async () => {
    const v = await open(ringStoryView('ring'));
    const l = ribbonsOf(v)!;
    const { stem } = strongestStem(l);
    hoverStrongestStem(v);
    expect(v.emphasis.hover?.element?.key).toBe(l.stemKey[stem]);
    expect(v.frame!.lists.passes.hover.paths.n).toBeGreaterThan(1);
  });

  it('7. EC-210: two assemblies on one ring: every ribbon from a peach stem to a grape stem', async () => {
    const v = await open(ringStoryView('ring-two-assemblies'));
    const l = ribbonsOf(v)!;
    expect(errors(v)).toEqual([]);
    expect(l.nr).toBeGreaterThan(0);
    const peach = (stem: number) => l.systems[l.stemSystem[stem]].includes('/peach/');
    for (let r = 0; r < l.nr; r++) expect(peach(l.ribA[r])).not.toBe(peach(l.ribB[r]));
    await exported(v);
  });

  it('8. crowded (3,000 contigs): stems and ribbons drawn, no error', async () => {
    const v = await open(ringStoryView('ring-crowded'));
    const l = ribbonsOf(v)!;
    expect(errors(v)).toEqual([]);
    expect(l.ns).toBeGreaterThan(0);
    expect(l.nr).toBeGreaterThan(0);
    await exported(v);
  });

  it('9. pending, failed, nothing: the hatch while the links wait; the issue once when failed; “no links”', async () => {
    let v = await open(ringStoryView('ring-pending'), Rondo, 1000);
    expect(ribbonsOf(v)!.empty).toBe('pending');
    expect(descriptions(v)).toContain('links loading');
    v = await open(ringStoryView('ring-failed'), Rondo, 10_000);
    expect(v.store.issues.value.filter((i) => i.severity === 'error').length).toBe(1);
    // a source that never opened is failed on the store (its issue says why): no ribbon, no hatch waiting, said
    const failed = ribbonsOf(v)!;
    expect(failed.nr).toBe(0);
    expect(failed.empty).toBe('failed');
    expect(failed.reader.pending ?? []).toEqual([]);
    expect(descriptions(v)).toContain('links failed');
    v = await open(ringStoryView('ring-no-links'));
    expect(ribbonsOf(v)!.empty).toBe('none');
    expect(descriptions(v)).toContain('no links');
  });

  it('10. dark: tokens only, the ribbons’ hues the scheme’s derived dark', async () => {
    const v = await open(ringStoryView('ring'));
    const { light, dark } = await exported(v);
    const fills = (svg: string) => new Set([...svg.matchAll(/fill="(#[0-9a-fA-F]{6})"/g)].map((m) => m[1]));
    const lf = fills(light.svg);
    const df = fills(dark.svg);
    expect(lf.size).toBeGreaterThan(10);
    let shared = 0;
    for (const f of lf) if (df.has(f)) shared++;
    expect(shared).toBeLessThan(lf.size / 2);
  });

  it('the showcase settles with no issue: the cytobands read for the assembly’s chromosomes alone', async () => {
    const said = (v: SeqViewer) =>
      [...v.modelIssues.value, ...v.layoutIssues.value, ...v.store.issues.value].map((i) => `${i.code}: ${i.message}`);
    for (const id of ['ring', 'ring-chromosome', 'ring-pair', 'ring-two-assemblies', 'ring-crowded']) {
      const v = await open(ringStoryView(id), Rondo, 5000);
      // the cytobands' rows in (the issue said when their source lists its sequences)
      if (id === 'ring' || id === 'ring-chromosome' || id === 'ring-pair') await bandsLanded(v);
      await vi.advanceTimersByTimeAsync(2000);
      v.draw();
      expect(said(v), id).toEqual([]);
    }
  });

  it('the hue ring’s legend lists the assembly’s chromosomes alone, never the stains its bands carry', async () => {
    const v = await open(ringStoryView('ring'));
    await bandsLanded(v);
    const legend = viewLegends(v.frame!.layout, Number.POSITIVE_INFINITY).legends.find(
      (l) => l.title === 'Chromosomes',
    )!;
    const labels = legend.entries.kind === 'swatches' ? legend.entries.items.map((i) => i.label) : [];
    // the assembly's 24 (its hue domain, taken before `include`: chrY keeps its hue when it joins)
    expect(labels).toEqual([...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY']);
  });

  it('hovering a band: the card names it and its stain; the ring lights exactly the ribbons with an end in it', async () => {
    for (const id of ['ring', 'ring-chromosome']) {
      const v = await open(ringStoryView(id));
      await bandsLanded(v);
      const l = ribbonsOf(v)!;
      const f = l.frame;
      const d = describerFor({ ...v, layout: () => v.frame?.layout });
      // three bands under the cytobands' track: a sweep round the ring at its radius
      const seen = new Set<string>();
      for (let a = 0.05; a < 2 * Math.PI && seen.size < 3; a += 0.11) {
        const x = f.cx + (f.r + 24) * Math.sin(a);
        const y = f.cy - (f.r + 24) * Math.cos(a);
        const h = v.hit(x, y);
        if (!h?.element || !h.layer.endsWith('/bands/bands') || seen.has(String(h.element.key))) continue;
        seen.add(String(h.element.key));
        v.input.move(x, y);
        v.draw();
        const e = v.emphasis.hover!;
        // the emphasis is the band: its extent the band's, on its chromosome
        const [, range] = String(h.element.key).split('@');
        const [start, end] = range.split('-').map(Number);
        expect(Array.from(e.extent!.parts[0].ranges), `${id} ${h.element.key}`).toEqual([start, end]);
        // the lit polygons: those whose stem (a ribbon's either end) has a position in the band, the hole washed
        const touches = (s: number) =>
          l.systems[l.stemSystem[s]] === h.system && l.stemStart[s] < end && l.stemEnd[s] > start;
        let lit = 0;
        for (let i = 0; i < l.polys.n; i++) {
          const k = l.polyKind[i];
          const r = l.polyRef[i];
          if (
            k === Kind.ribbon
              ? touches(l.ribA[r]) || touches(l.ribB[r])
              : k === Kind.stub
                ? touches(l.stubStem[r])
                : touches(r)
          )
            lit++;
        }
        // the ribbons' own hover pass: the wash, then each lit polygon
        const lists = new DrawLists();
        RibbonsType.draw(l, 'hover', lists.writer('hover') as PassWriter, v.emphasis);
        expect(lists.passes.hover.paths.n, `${id} ${h.element.key}`).toBe(1 + lit);
        expect(lit).toBeLessThan(l.polys.n - 1);
        // the card's words name the band and its stain
        const card = d.describe(v.hover.current)!;
        expect(card.title).toBe(h.state);
        expect(card.text.startsWith(`${h.state}, Stain: `)).toBe(true);
      }
      expect(seen.size, id).toBe(3);
      // nothing held past the draw
      expect(emphasisHeld()).toBe(false);
    }
  });
});
