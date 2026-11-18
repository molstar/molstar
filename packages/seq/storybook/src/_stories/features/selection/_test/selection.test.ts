/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { recordDrawLists } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { GeneFile } from '@molstar/seqx-examples/p53/data';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  cellPoint,
  clickAt,
  composeViewer,
  describeSelection,
  dragFrom,
  type SelectionLook,
  selectionSvg,
} from '../compose.js';

// The Selection story, headless on a fake clock, composed as the story composes it: on both p53 views a click, a drag
// and the cursor through the viewer's input; the selection drawn in the canvas's lists and in the SVG export beside it
// alike, and projected onto the gene's section in both. The performance representative: a drag-select step and a cursor
// step (input, then the frame) under 4 ms.

const tracks = P53Views.find((v) => v.id === 'tracks')!;
const onGene = P53Views.find((v) => v.id === 'protein-on-gene')!;

let gene = '';
beforeAll(async () => {
  gene = await nodeFixtureAssets().reader.readText(GeneFile, new AbortController().signal);
});

let viewer: SeqViewer | undefined;

// a frame tick asked for by an earlier file on timers since dropped would never run: this file's own host function
const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** `clock`: the bound of the performance case reads the real `performance`, its timers faked alone. */
async function make(view = tracks, clock = true, look?: SelectionLook) {
  vi.useFakeTimers({
    toFake: clock ? ['setTimeout', 'clearTimeout', 'Date', 'performance'] : ['setTimeout', 'clearTimeout'],
  });
  const v = composeViewer(view, {
    ports: { assets: () => ({ reader: { readText: async () => gene } }) },
    ...(look ? { look } : {}),
  });
  viewer = v;
  v.input.mac = false;
  v.resize(720);
  await v.load(view.build());
  await vi.advanceTimersByTimeAsync(400);
  v.draw();
  return v;
}

/** The selection pass's rects within a section's band, as recorded. */
function band(v: SeqViewer, section: number): string[] {
  const s = v.frame!.layout.sections[section];
  return recordDrawLists(v.frame!.lists, ['selection'])
    .split('\n')
    .filter((l) => l.startsWith('rect ') && Number(l.split(' ')[2]) >= s.y && Number(l.split(' ')[2]) < s.y + s.height);
}

describe('Selection, headless', () => {
  it('p53 tracks: a click on a domain selects it; the SVG beside the canvas draws the same selection', async () => {
    const v = await make();
    v.view.setRange('/views/main/protein', 90, 130);
    await vi.advanceTimersByTimeAsync(200);
    v.draw();
    // the DNA-binding domain covers 100, in whichever lane it took
    const [x, top] = cellPoint(v, 0, 1, 100);
    const y = [0, 1, 2].map((lane) => top + lane * 12).find((ly) => v.hit(x, ly)?.element?.key === 'DBD')!;
    clickAt(v, [x, y]);
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['DBD']);
    expect(describeSelection(v)).toContain('(item :id DBD)');
    v.draw();
    const svg = await selectionSvg(v, 't-');
    expect(recordDrawLists(svg.lists, ['selection'])).toBe(recordDrawLists(v.frame!.lists, ['selection']));
    expect(band(v, 0).length).toBeGreaterThan(0);
  });

  it('protein on its gene: a residue selected on the protein shows on the gene, lighter, in the canvas and the SVG', async () => {
    const v = await make(onGene);
    v.view.setRange('/views/main/protein', 170, 180);
    await vi.advanceTimersByTimeAsync(400);
    v.draw();
    clickAt(v, cellPoint(v, 0, 0, 174));
    v.draw();
    expect(band(v, 0).some((l) => l.includes('#8f4c062e'))).toBe(true);
    expect(band(v, 1).some((l) => l.includes('#8f4c0617'))).toBe(true);
    const svg = await selectionSvg(v, 'g-');
    expect(recordDrawLists(svg.lists, ['selection'])).toBe(recordDrawLists(v.frame!.lists, ['selection']));
  });

  it('protein on its gene: a drag on the gene selects nucleotides, projected onto the protein', async () => {
    const v = await make(onGene);
    v.view.setRange('/views/main/protein', 170, 180);
    await vi.advanceTimersByTimeAsync(400);
    v.draw();
    const s = v.frame!.layout.sections[1];
    const y = s.tracks[0].box.y + 6;
    const vis = s.axis!.visible([]);
    const a = Math.floor(vis[0]) + 3;
    dragFrom(v, [cellPoint(v, 1, 0, a)[0], y], [cellPoint(v, 1, 0, a + 8)[0], y]);
    const sel = v.selection.value.value.shown;
    expect(sel.parts.map((p) => p.system)).toEqual(['/data/tp53']);
    v.draw();
    expect(band(v, 0).some((l) => l.includes('#8f4c0617'))).toBe(true);
    expect(band(v, 1).some((l) => l.includes('#8f4c062e'))).toBe(true);
  });

  it('the cursor walks the protein with the arrow keys, and Shift extends the selection from it', async () => {
    const v = await make();
    v.input.key({ key: 'ArrowRight' });
    v.input.key({ key: 'ArrowRight' });
    v.input.key({ key: 'ArrowRight', shiftKey: true });
    v.input.key({ key: 'ArrowRight', shiftKey: true });
    expect(v.selection.cursor.value?.position).toBe(3);
    expect([...v.selection.value.value.shown.parts[0].ranges]).toEqual([1, 4]);
    v.input.key({ key: 'Escape' });
    expect(v.selection.value.value.shown.parts).toEqual([]);
  });
});

describe('Selection: the performance representative', () => {
  // the tint over the data (the default) redraws the emphasis alone; under it, a selection's change draws the base too
  for (const place of ['over', 'under'] as const)
    itWithin(
      `a drag-select step and a cursor step (the input, then the frame) under 4 ms each, the tint ${place}`,
      async (within) => {
        const v = await make(onGene, false, { place });
        v.view.setRange('/views/main/protein', 150, 250);
        await vi.advanceTimersByTimeAsync(400);
        v.draw();
        const [x0, y] = cellPoint(v, 0, 0, 160);
        const [x1] = cellPoint(v, 0, 0, 240);
        v.input.move(x0, y);
        v.input.down(x0, y);
        let i = 0;
        const dragStep = () => {
          i++;
          v.input.move(x0 + ((x1 - x0) * (i % 40)) / 40 + 4, y);
          v.draw();
        };
        const cursorStep = () => {
          v.input.key({ key: i++ % 20 < 10 ? 'ArrowRight' : 'ArrowLeft' });
          v.draw();
        };
        // warm-ups
        for (let k = 0; k < 40; k++) dragStep();
        const n = 10;
        const t0 = performance.now();
        for (let k = 0; k < n; k++) within(4, dragStep);
        console.info(`drag-select step (${place}): ${((performance.now() - t0) / n).toFixed(2)} ms`);
        v.input.up();
        v.selection.setCursor({ section: '/views/main/protein', system: '/data/p53', position: 200 });
        for (let k = 0; k < 40; k++) cursorStep();
        const t1 = performance.now();
        for (let k = 0; k < n; k++) within(4, cursorStep);
        console.info(`cursor step (${place}): ${((performance.now() - t1) / n).toFixed(2)} ms`);
        expect(v.selection.value.value.shown.parts.length).toBeGreaterThan(0);
      },
    );
});

it('the status line names nothing before a selection', async () => {
  const v = await make();
  expect(describeSelection(v)).toBe('selection: — · cursor: —');
});
