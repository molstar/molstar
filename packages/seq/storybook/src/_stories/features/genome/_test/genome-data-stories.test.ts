/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type BlockPlacement, DarkTokens, hitView, LightTokens, recordDrawLists } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { GenomeViews } from '../views.js';

// Features/Genome Coordinates, data over an assembly, headless on a fake clock as the stories compose them (a viewer of
// the view, its navigator, its SVG export light and dark): a BED over the toy assembly's rows, and hg38's cytobands over
// every chromosome as a bands track and as the ideogram navigator (EC-187, EC-189).

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

/** The element keys hit along a strip's band, a few px apart. */
function keysOn(v: SeqViewer, strip: number, band: number): Set<string> {
  const p = v.frame!.layout.sections[0].placement as BlockPlacement;
  const s = p.strips;
  const a = s.axes[strip];
  const y = s.bandY[strip * s.bands + band] + s.bandH[strip * s.bands + band] / 2;
  const keys = new Set<string>();
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
    for (let i = 0; i < 100 && !v.view.model.value.issues.length; i++) {
      await new Promise((r) => setImmediate(r));
      await vi.advanceTimersByTimeAsync(200);
    }
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
});
