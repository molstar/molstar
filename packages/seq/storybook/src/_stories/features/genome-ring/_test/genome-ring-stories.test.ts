/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import type { ExampleView } from '@molstar/seqx-examples';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { layerAria } from '@molstar/seqx-viewer/layer-aria';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { hoverStrongestStem, ribbonsOf, ringStoryView, strongestStem } from '../views.js';

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
    // a source that never opened has no coverage, so core's tile reads call its tiles pending: no ribbon, the hatch
    // stays (the issue says why)
    expect(ribbonsOf(v)!.nr).toBe(0);
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
});
