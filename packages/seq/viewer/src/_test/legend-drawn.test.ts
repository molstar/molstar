// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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
const N = 30;

function spec(kind: 'variants' | 'bands') {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'dna', letters: 'A'.repeat(N * 10) });
  const stains = ['gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100', 'gvar'];
  const rows =
    kind === 'bands'
      ? {
          start: Array.from({ length: N }, (_, i) => i * 10),
          end: Array.from({ length: N }, (_, i) => i * 10 + 8),
          id: Array.from({ length: N }, (_, i) => `i${i}`),
          // the first ten bands one stain, the rest the others: the left of the view holds one
          stain: Array.from({ length: N }, (_, i) => (i < 10 ? 'gneg' : stains[1 + (i % 5)])),
        }
      : {
          start: Array.from({ length: N }, (_, i) => i * 10),
          end: Array.from({ length: N }, (_, i) => i * 10 + 1),
          id: Array.from({ length: N }, (_, i) => `i${i}`),
          category: Array.from({ length: N }, (_, i) => `chr${i + 1}`),
        };
  const d = b.data.ranges('d', { target: p, rows: rows as never });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('chroms', { label: 'Chromosomes' }).layer(
      'd',
      kind === 'bands'
        ? { type: 'bands', data: d, color: theme.field('stain', { list: false }) }
        : { type: 'variants', data: d, color: theme.field('category', { list: false }) },
    );
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

async function load(kind: 'variants' | 'bands') {
  frames = installFrames();
  viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      [viewerItem('view.motion'), 'cut'],
    ],
  });
  await viewer.load(spec(kind));
  viewer.resize(500);
  await viewer.settle();
  frames.flush();
  await rest();
}

const total = () => {
  const r = legendRows(viewer!, { of: Chromosomes })[0];
  return r.entries.length + r.more;
};

describe('a legend lists the values drawn, never those only read', () => {
  it('variants: a view over 3 of 30 rows lists their values', async () => {
    await load('variants');
    expect(total()).toBe(30);
    viewer!.runCommand('view.focus', { section: '/views/main/p', focus: { start: 0, end: 30 } });
    await rest();
    expect(total()).toBeLessThan(30);
    expect(total()).toBeGreaterThan(0);
  });

  it('bands: a view over the first bands lists their stains only', async () => {
    await load('bands');
    const all = total();
    viewer!.runCommand('view.focus', { section: '/views/main/p', focus: { start: 0, end: 60 } });
    await rest();
    expect(total()).toBeLessThan(all);
    expect(total()).toBeGreaterThan(0);
  });
});
