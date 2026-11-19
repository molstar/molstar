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
import { openExample } from './draw-checks.js';

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
});
