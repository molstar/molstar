/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text stories, headless on a fake clock, composed as the story composes them (a viewer of the view and its SVG
// export): the canvas shows the lines in view, the export every line; a hit on a letter through the placement names its
// position; the letters drawn are the protein's.

import { hitView, recordDrawLists, type TextPlacement } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { P53 } from '@molstar/seqx-examples/p53/data';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const lettersOf = (record: string) =>
  record
    .split('\n')
    .filter((l) => l.startsWith('glyphs'))
    .map((l) =>
      l
        .replace(/^glyphs "[^"]*" [\d.]+ \S+ /, '')
        .split(', ')
        .map((x) => x[0])
        .join(''),
    )
    .join('');

describe('Text, headless', () => {
  it('the canvas shows lines of the protein in reading order; a letter’s hit names its position; the export every line', async () => {
    viewer = await openExample(P53Views.find((v) => v.id === 'text')!, { gene: '' });
    const v = viewer;
    const sec = v.frame!.layout.sections[0];
    const p = sec.placement as TextPlacement;
    expect(p.strips.n).toBeGreaterThan(1);
    const shown = lettersOf(recordDrawLists(v.frame!.lists, ['data']));
    expect(shown.length).toBeGreaterThan(0);
    expect(P53.startsWith(shown)).toBe(true);
    // the second line's fifth letter
    const a = p.strips.axes[1];
    const x = a.xOfUnit(a.u0 + 4.5);
    const y = p.strips.top[1] + p.strips.bandY[p.strips.bands] + 4;
    expect(hitView(v.frame!.layout, x, y)?.position).toBe(a.u0 + 4);
    const svg = v.renderSvg({ idPrefix: 'x-' });
    await vi.advanceTimersByTimeAsync(400);
    const r = await svg;
    expect(r.issues).toEqual([]);
    expect(lettersOf(recordDrawLists(r.lists, ['data']))).toBe(P53);
  });

  it('the wheel scrolls the text’s lines only while it overflows its box, and is given back at its ends', async () => {
    viewer = await openExample(P53Views.find((v) => v.id === 'text')!, { gene: '' });
    const v = viewer;
    const at = () => {
      const sec = v.frame!.layout.sections[0];
      return { p: sec.placement as TextPlacement, y: sec.y + 20 };
    };
    // wide: every line fits its box, the page keeps the wheel
    v.resize(1400);
    v.draw();
    expect(at().p.contentHeight).toBeLessThanOrEqual(at().p.box.height);
    expect(v.input.wheel(400, at().y, { deltaX: 0, deltaY: 40 })).toBe(false);
    // narrow: the text overflows; down scrolls its lines, up at its top gives the wheel back
    v.resize(360);
    v.draw();
    expect(at().p.contentHeight).toBeGreaterThan(at().p.box.height);
    expect(v.input.wheel(200, at().y, { deltaX: 0, deltaY: -40 })).toBe(false);
    expect(v.input.wheel(200, at().y, { deltaX: 0, deltaY: 120 })).toBe(true);
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    expect(at().p.anchor).toBeGreaterThan(0);
    expect(at().p.strips.axes[0].u0).toBeGreaterThan(0);
  });
});
