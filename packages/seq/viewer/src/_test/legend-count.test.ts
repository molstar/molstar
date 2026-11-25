// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesLay } from '@molstar/seqx-core/render';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { legendRows } from '../legends.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { stubCanvasContext } from './canvas-stub.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
beforeEach(stubCanvasContext);
afterEach(() => {
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
});

const Chromosomes = '/views/main/p/chroms';

function spec(n: number) {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'dna', letters: 'A'.repeat(n * 10) });
  const d = b.data.ranges('d', {
    target: p,
    rows: {
      start: Array.from({ length: n }, (_, i) => i * 10),
      end: Array.from({ length: n }, (_, i) => i * 10 + 8),
      id: Array.from({ length: n }, (_, i) => `i${i}`),
      category: Array.from({ length: n }, (_, i) => `chr${i + 1}`),
    },
  });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('chroms', { label: 'Chromosomes' }).layer('d', {
      type: 'ranges',
      data: d,
      color: theme.field('category', { list: false }),
    });
  });
  return b.build();
}

async function rest() {
  for (let k = 0; k < 3; k++) {
    await viewer!.settle();
    frames!.timers.fire();
    frames!.flush();
  }
}

async function load(n: number) {
  frames = installFrames();
  viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      [viewerItem('view.motion'), 'cut'],
    ],
  });
  await viewer.load(spec(n));
  viewer.resize(500);
  await viewer.settle();
  frames.flush();
}

describe('a legend’s values in view', () => {
  it('lists the values drawn, never those only read (out of view)', async () => {
    await load(30);
    await rest();
    expect(legendRows(viewer!, { of: Chromosomes })[0].more).toBe(18);
    // focused on one item: the window holds a few of the 30 items, the legend lists what it draws
    viewer!.runCommand('view.focus', { section: '/views/main/p', focus: { start: 0, end: 30 } });
    await rest();
    const row = legendRows(viewer!, { of: Chromosomes })[0];
    const entries = row.entries.length + row.more;
    // the values of the marks the layer draws (its hit marks), read off the layout
    const drawn = new Set<string>();
    for (const sec of viewer!.frame!.layout.sections)
      for (const t of sec.tracks)
        for (const l of t.layers) {
          const lay = l.lay as RangesLay;
          if (l.type.name !== 'ranges') continue;
          for (let k = 0; k < lay.hn; k++) {
            const v = lay.painted[lay.hitem[k]];
            if (v !== undefined) drawn.add(v);
          }
        }
    expect(entries).toBe(drawn.size);
    expect(entries).toBeLessThan(30);
    expect(entries).toBeGreaterThan(0);
  });
});
