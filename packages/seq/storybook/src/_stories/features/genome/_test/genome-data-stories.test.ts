/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type BlockPlacement, DarkTokens, hitView, LightTokens, recordDrawLists } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample, turnUntil } from '../../../examples/_test/draw-checks.js';
import { GenomeViews } from '../views.js';

// Features/Genome Coordinates, data over an assembly, headless on a fake clock as the stories compose them (a viewer of
// the view, its navigator, its SVG export light and dark): a BED over the toy assembly's rows, hg38's cytobands over
// every chromosome as a bands track and as the ideogram navigator (EC-187, EC-189), and data written in the spec over
// an assembly by its name columns: links on a chord ring (`a_name`, `b_name`), a matrix over two chromosomes (`name`).

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

const view = (id: string) => GenomeViews.find((v) => v.id === id)!;

async function open(id: string) {
  viewer = await openExample(view(id), { delayMs: 0, wait: 2000 });
  await vi.advanceTimersByTimeAsync(2000);
  viewer.draw();
  return viewer;
}

async function exported(v: SeqViewer, dark: boolean) {
  const r = v.renderSvg({ idPrefix: 'gd-', tokens: dark ? DarkTokens : LightTokens, navigators: 'include' } as never);
  await vi.advanceTimersByTimeAsync(2000);
  return r;
}

/** The element keys hit along a strip's band, a few px apart, on each of its lanes (every 4 px down it). */
function keysOn(v: SeqViewer, strip: number, band: number): Set<string> {
  const p = v.frame!.layout.sections[0].placement as BlockPlacement;
  const s = p.strips;
  const a = s.axes[strip];
  const top = s.bandY[strip * s.bands + band];
  const height = s.bandH[strip * s.bands + band];
  const keys = new Set<string>();
  for (let y = top + 2; y < top + height; y += 4)
    for (let x = a.x0; x < a.x0 + a.width; x += 0.5) {
      const h = hitView(v.frame!.layout, x, y);
      if (h?.element) keys.add(String(h.element.key));
    }
  return keys;
}

describe('Genome coordinates: data over an assembly, headless', () => {
  it('a BED over the assembly (EC-189): each chromosome’s block draws the file’s own features, named without `chr`', async () => {
    const v = await open('bed-over-rows');
    const p = v.frame!.layout.sections[0].placement as BlockPlacement;
    expect(p.strips.n).toBe(24);
    // the features band is the second (the ruler first): every block's features its own chromosome's
    for (const i of [0, 9, 22]) {
      const name = p.labels[i].replace(/^chr/, '');
      const keys = [...keysOn(v, i, 1)];
      expect(keys.length, p.labels[i]).toBeGreaterThan(3);
      expect(keys.every((k) => k.startsWith(`${name}.f`))).toBe(true);
    }
    expect(v.view.model.value.issues).toEqual([]);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect([svg.settled, svg.issues]).toEqual([true, []]);
    }
  });

  it('cytobands over the genome: a bands track on every chromosome and the ideogram over all, the file’s other contigs said once', async () => {
    const v = await open('cytobands');
    // the recorded file is read from the disk (real I/O, off the fake clock): moved until its names are listed
    await turnUntil(() => v.view.model.value.issues.length > 0, 200);
    await vi.advanceTimersByTimeAsync(2000);
    v.draw();
    const p = v.frame!.layout.sections[0].placement as BlockPlacement;
    // chrM's 16.6 kb is under a device pixel at the fit: no strip
    expect(p.strips.n).toBe(24);
    expect(v.view.model.value.issues.map((i) => i.code)).toEqual(['seq.resolve.assembly.name']);
    expect(v.view.model.value.issues[0].message).toMatch(/^The source names 686 sequences the assembly lacks/);
    // the ideogram: every chromosome's bands, in the navigator's picture
    const svg = await exported(v, false);
    expect(svg.settled).toBe(true);
    expect(svg.issues.map((i) => i.code)).toEqual(['seq.resolve.assembly.name']);
    const dark = await exported(v, true);
    expect(dark.svg).not.toBe(svg.svg);
    expect(
      recordDrawLists(v.frame!.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('rect ')).length,
    ).toBeGreaterThan(500);
  });

  it('links written in the spec on a ring (`a_name`, `b_name`): chords between chromosomes, every row read', async () => {
    const v = await open('inline-links-ring');
    expect(v.view.model.value.issues).toEqual([]);
    const lists = recordDrawLists(v.frame!.lists, ['data']);
    const chords = lists.split('\n').filter((l) => l.startsWith('path ') && l.includes(' C '));
    // sixty links, a few within one chromosome: every one drawn (none past the gate at the fit)
    expect(chords.length).toBe(60);
    // a link between two chromosomes has its ends on two units
    const keys = new Set<string>();
    const { width, height } = v.frame!.layout;
    for (let x = 0; x < width; x += 4)
      for (let y = 0; y < height; y += 4) {
        const h = hitView(v.frame!.layout, x, y);
        if (h?.element && h.layer.endsWith('/links')) keys.add(String(h.element.key));
      }
    // held as a complex's links are, each end on its chromosome's unit: the links hit by their rows
    const pairs = v.view.model.value.annotations.get('/data/links')!.pairs!;
    expect(pairs.n).toBe(60);
    const across = [...keys].filter((k) => pairs.aUnit![Number(k)] !== pairs.bUnit![Number(k)]);
    expect(across.length).toBeGreaterThan(10);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect([svg.settled, svg.issues]).toEqual([true, []]);
      expect(svg.svg.match(/<path/g)?.length ?? 0).toBeGreaterThanOrEqual(60);
    }
  });

  it('a matrix written in the spec over two chromosomes (`name`): each block draws its own items', async () => {
    const v = await open('inline-matrix');
    expect(v.view.model.value.issues).toEqual([]);
    const p = v.frame!.layout.sections[0].placement as BlockPlacement;
    expect(p.labels).toEqual(['chrI', 'chrII']);
    // the matrix band (after the ruler) holds cells over each block: a rect a run of items and sample
    const s = p.strips;
    const y0 = s.bandY[1];
    const y1 = y0 + s.bandH[1];
    const rects = recordDrawLists(v.frame!.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .map((l) => l.split(' ').slice(1, 3).map(Number));
    for (const i of [0, 1]) {
      const a = s.axes[i];
      const cells = rects.filter(([x, y]) => x >= a.x0 && x < a.x0 + a.width && y >= y0 && y < y1);
      expect(cells.length, p.labels[i]).toBeGreaterThan(20);
    }
    const light = await exported(v, false);
    const dark = await exported(v, true);
    expect([light.settled, light.issues, dark.settled, dark.issues]).toEqual([true, [], true, []]);
    expect(dark.svg).not.toBe(light.svg);
  });
});
