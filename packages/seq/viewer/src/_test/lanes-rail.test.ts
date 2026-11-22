/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The lanes rail by gesture (core's _spec/README.md "Render", Track heights; the rows rail's gesture, in lanes): a band
// whose lanes overflow has the rows rail itself at its right edge: its window dragged scrolls the lanes (the stored
// `bandTop`, session state: no history step), a press beside it centres the lanes there, Escape cancels, Page Up and
// Page Down scroll with the pointer in the band, and the pan-only rail takes a grip as its body.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { LanesRail, type RangesLay, type Rail } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Tall = '/views/v/s/tall';

function spec(): Node {
  const b = createSeqSpec({ title: 'lanes' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  const rows = {
    start: Array.from({ length: 40 }, (_, i) => i),
    end: Array.from({ length: 40 }, () => 150),
    id: Array.from({ length: 40 }, (_, i) => `r${i}`),
  };
  const tall = b.data.ranges('tall', { target: s, rows });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s, heights: 'fixed' } as never);
    sec.track('tall', { height: 240 } as never).layer('l', { type: 'ranges', data: tall });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make() {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, config: [[viewerItem('view.jumpMs'), 0]] });
  viewer = v;
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  v.draw();
  v.input.mac = false;
  return v;
}

const rail = (v: SeqViewer): Rail => v.frame!.layout.sections[0].lanesRails![0];
const lay = (v: SeqViewer) => {
  for (const t of v.frame!.layout.sections[0].tracks)
    for (const l of t.layers) if (l.key.startsWith(Tall)) return l.lay as RangesLay;
  throw new Error('no lanes');
};
const yOf = (r: Rail, u: number) => r.x + (u / r.total) * r.width;
const top = (v: SeqViewer) => v.view.band(Tall)?.bandTop ?? 0;
const settle = (v: SeqViewer) => {
  frames!.flush();
  v.draw();
};

describe('the lanes rail by gesture', () => {
  it('is the rows rail: a target with its window the lanes in view; a grip is the window’s body', async () => {
    const v = await make();
    const r = rail(v);
    expect(r.active).toBe(true);
    expect(r.end - r.start).toBeCloseTo(lay(v).room, 6);
    const x = r.top + 6;
    const h = v.hit(x, yOf(r, 1));
    expect(h?.navigator).toBe(LanesRail);
    expect(h?.part).toBe('body');
    expect(v.hit(x, yOf(r, r.end + 5))?.part).toBe('outside');
    // the window's very start is its grip on the rows rail, body here
    expect(v.hit(x, yOf(r, 0) + 0.5)?.part).toBe('body');
  });

  it('a drag on its window scrolls the lanes (clamped, session state: no step); Escape puts them back', async () => {
    const v = await make();
    const r = rail(v);
    const x = r.top + 6;
    const y0 = yOf(r, 1);
    expect(v.input.down(x, y0)).toBe(true);
    v.input.move(x, y0 + (6 / r.total) * r.width);
    expect(top(v)).toBeCloseTo(6, 6);
    expect(v.input.cursor).toBe('grabbing');
    v.input.move(x, yOf(r, 500));
    expect(top(v)).toBeCloseTo(lay(v).scrollMax, 6);
    v.input.up();
    expect(v.view.history.value.back).toBe(0);
    settle(v);
    expect(lay(v).top).toBeCloseTo(lay(v).scrollMax, 6);
    // Escape: back where the press found it
    const r2 = rail(v);
    const y1 = yOf(r2, r2.start + 1);
    v.input.down(x, y1);
    v.input.move(x, y1 - 40);
    expect(top(v)).toBeLessThan(lay(v).scrollMax);
    v.input.key({ key: 'Escape' });
    expect(top(v)).toBeCloseTo(lay(v).scrollMax, 6);
  });

  it('a press beside the window centres the lanes there', async () => {
    const v = await make();
    const r = rail(v);
    const x = r.top + 6;
    const y = yOf(r, 25);
    expect(v.hit(x, y)?.part).toBe('outside');
    v.input.down(x, y);
    v.input.up();
    v.input.click(x, y, { detail: 1 });
    expect(top(v)).toBeCloseTo(25 - lay(v).room / 2, 6);
    expect(v.view.history.value.back).toBe(0);
  });

  it('Page Down and Page Up with the pointer in the band', async () => {
    const v = await make();
    const l = lay(v);
    const x = l.box.x + 20;
    const y = l.box.y + 30;
    v.input.move(x, y);
    const page = Math.floor(l.room) - 1;
    expect(v.input.key({ key: 'PageDown' })).toBe(true);
    expect(top(v)).toBe(page);
    // clamped at the last lane
    v.input.key({ key: 'PageDown' });
    expect(top(v)).toBe(Math.min(l.scrollMax, 2 * page));
    v.input.key({ key: 'PageUp' });
    expect(top(v)).toBe(Math.min(l.scrollMax, 2 * page) - page);
    v.input.key({ key: 'PageDown' });
    v.input.key({ key: 'PageDown' });
    // at the end the key is the page's
    expect(top(v)).toBe(l.scrollMax);
    expect(v.input.key({ key: 'PageDown' })).toBe(false);
    v.view.setBand(Tall, { bandTop: 0 });
    expect(v.input.key({ key: 'PageUp' })).toBe(false);
  });
});
