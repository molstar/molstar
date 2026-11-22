/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { go, openEmpty, ReadableWindow } from '../openings.js';
import { advance, exported, openRegion } from './region-twin.js';

// Features/Region reads/Zoom in, headless on a fake clock (the story's composition): a matrix served at its items'
// level alone, past what a view asks the store for (its band the waiting hatch, nothing read, no wait resolving it),
// and links with nothing in the stretch in view (no arc, no placeholder); zoomed in the matrix is read (EC-169, EC-208).

let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

interface MatrixLay {
  readonly pend: { readonly length: number };
  readonly items: boolean;
  readonly tables: { readonly length: number };
  readonly n: number;
}

interface PairsLay {
  readonly n: number;
  readonly carets: { readonly length: number };
  readonly densityN: number;
}

const matrix = (v: SeqViewer) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === 'Genotypes')!.layers[0].lay as unknown as MatrixLay;
const links = (v: SeqViewer) => v.frame!.layout.sections[0].sectionLays[0].lay as unknown as PairsLay;

describe('Features/Region reads/Zoom in, headless', () => {
  it('the genotypes past what a view asks: the waiting hatch, nothing read however long; links with nothing in view draw nothing', async () => {
    const v = (viewer = await openRegion('zoom-in'));
    await advance(v, 2000);
    openEmpty(v);
    for (const t of [500, 5000, 30_000]) {
      await advance(v, t);
      // items at 300 kb: more tiles than a want is planned for (`maxTilesPerWant`): none is read, the band stays hatched
      expect(matrix(v).items).toBe(true);
      expect(matrix(v).tables.length).toBe(0);
      expect(matrix(v).pend.length).toBeGreaterThan(0);
      expect(matrix(v).n).toBe(0);
    }
    // the links: no end in 400 to 700 kb, no arc, no caret, no density strip, nothing pending
    expect(links(v).n).toBe(0);
    expect(links(v).carets.length).toBe(0);
    expect(links(v).densityN).toBe(0);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues).toEqual([]);
    }
  });

  it('zoomed in to 20 kb the genotypes are read and drawn', async () => {
    const v = (viewer = await openRegion('zoom-in'));
    await advance(v, 2000);
    openEmpty(v);
    await advance(v, 3000);
    expect(matrix(v).tables.length).toBe(0);
    go(v, 'chr', ReadableWindow);
    await advance(v, 5000);
    expect(matrix(v).pend.length).toBe(0);
    expect(matrix(v).tables.length).toBeGreaterThan(0);
    expect(matrix(v).n).toBeGreaterThan(0);
    // the links panned to their own stretch: the arcs are there
    go(v, 'chr', [0, 300_000]);
    await advance(v, 5000);
    expect(links(v).n).toBeGreaterThan(0);
  });
});
