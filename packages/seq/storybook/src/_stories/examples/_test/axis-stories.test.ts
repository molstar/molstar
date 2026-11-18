/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The axis's stories, headless on a fake clock, each composed as its story composes it (a viewer of the view and its
// SVG export): the protein and its gene on one axis drawn by the canvas and by `renderSvg` from one spec, their hits
// and the export's `data-*` agreeing; a follower on a reversed collapsed axis moved by its leader and leading it back.

import { hitView, type LettersLay, recordDrawLists } from '@molstar/seqx-core/render';
import { positionAt } from '@molstar/seqx-core/view';
import { P53Views } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { GeneFile } from '@molstar/seqx-examples/p53/data';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

const Gene = '/views/main/gene';
const Protein = '/views/main/protein';

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
  viewer = await openExample(P53Views.find((v) => v.id === id)!, { gene, delayMs: 300 });
  return viewer;
}

/** The viewer's SVG export, its asynchronous work let go on the fake clock. */
async function svgOf(v: SeqViewer) {
  const p = v.renderSvg({ idPrefix: 'x-' });
  await vi.advanceTimersByTimeAsync(400);
  return p;
}

describe('Protein and gene on one axis, headless', () => {
  it('the canvas and the SVG draw the same lists; a residue’s hit and the export’s data-* agree', async () => {
    const v = await open('one-axis');
    v.view.setRange(Gene, 450, 480);
    v.draw();
    const svg = await svgOf(v);
    expect(svg.issues).toEqual([]);
    for (const pass of ['data', 'labels'] as const)
      expect(recordDrawLists(svg.lists, [pass])).toBe(recordDrawLists(v.frame!.lists, [pass]));
    expect(svg.svg).toContain(`data-section="${Gene}" data-direction="reversed" data-pieces="coding"`);
    // a residue under the pointer: its element on p53, its codon on the gene, the cell in both pictures
    const sec = v.frame!.layout.sections[0];
    const lay = sec.tracks[0].layers[0].lay as LettersLay;
    expect(lay.n).toBeGreaterThan(0);
    const j = Math.floor(lay.n / 2);
    const h = hitView(v.frame!.layout, (lay.x0[j] + lay.x1[j]) / 2, lay.box.y + 4)!;
    expect(h.extent?.parts[0].system).toBe('/data/p53');
    expect(h.extent?.parts[0].ranges[0]).toBe(lay.pos[j]);
    expect(h.hull![1] - h.hull![0]).toBe(3);
    // the export's rect of that cell is the canvas's
    const rect = `rect ${lay.x0[j]} ${lay.box.y} ${lay.x1[j] - lay.x0[j]} ${lay.box.height}`;
    expect(recordDrawLists(svg.lists, ['data'])).toContain(rect);
    // a click selects the residue; the SVG shows the selection the canvas shows
    v.selection.select(h, v.view.resolved.value!.sections[0].axis);
    v.draw();
    const selected = await svgOf(v);
    expect(recordDrawLists(selected.lists, ['selection'])).toBe(recordDrawLists(v.frame!.lists, ['selection']));
  });
});

describe('Protein on its gene (a follower on a reversed collapsed axis), headless', () => {
  it('moved by its leader the same way on screen, leading it back', async () => {
    const v = await open('protein-on-gene');
    const axis = v.view.resolved.value!.sections[1].axis;
    expect([axis.dir, axis.pieces]).toEqual([-1, 'coding']);
    v.view.setRange(Protein, 100, 120);
    const g0 = v.view.get(Gene)!;
    // three units a residue: the coding axis's offsets are three times the protein's positions
    expect(g0.start).toBeCloseTo(300, 6);
    expect(g0.end).toBeCloseTo(360, 6);
    v.view.pan(Protein, 5);
    expect(v.view.get(Gene)!.start).toBeCloseTo(315, 6);
    // the gene leads back: its pan moves the protein the same way, a third as far
    v.view.pan(Gene, -15);
    expect(v.view.get(Protein)!.start).toBeCloseTo(100, 6);
    expect(v.view.get(Gene)!.start).toBeCloseTo(300, 6);
    // the gene's window in positions: read down (a reversed axis), never over a hidden intron
    const g = v.view.get(Gene)!;
    const a = positionAt(axis, g.start, 'start')!;
    const b = positionAt(axis, g.end, 'end')!;
    expect(a).toBeGreaterThan(b);
  });
});
