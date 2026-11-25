// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { BridgeView } from '../models/bridge-hooks.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// The bands' ARIA groups in a mounted viewer (_spec/README.md "Bridges", ARIA words): set when a frame settles, a
// group a bridge labelled by its `label`, its words the ribbons in view.

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

function spec() {
  const b = createSeqSpec({ title: 'bridge' });
  const x = b.data.sequence('x', { alphabet: 'dna', letters: 'ACGT'.repeat(50) });
  const y = b.data.sequence('y', { alphabet: 'dna', letters: 'ACGT'.repeat(50) });
  const p = b.data.pairs('p', {
    a: x,
    b: y,
    rows: { a_start: [0, 50], a_end: [40, 90], b_start: [10, 100], b_end: [50, 140] },
  });
  b.views.view('v', (v) => {
    v.section('peach', { of: x }).track('t').layer('l', { type: 'letters' });
    v.bridge('br', { a: '/views/v/peach', b: '/views/v/grape', label: 'Peach to grape' } as never, (br) =>
      br.layer('ribbons', { type: 'pairs', data: p } as never),
    );
    v.section('grape', { of: y }).track('t').layer('l', { type: 'letters' });
  });
  return b.build();
}

describe('a bridge’s ARIA group', () => {
  it('is a group labelled by the bridge’s label, its words the ribbons in view', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    await viewer.load(spec());
    viewer.resize(800);
    await viewer.settle();
    frames.flush();
    viewer.draw();
    // what the mount calls when the frame settles (`mount.ts`)
    const container = document.createElement('div');
    BridgeView.parts!.aria(viewer, container);
    const groups = [...container.querySelectorAll('[role="group"]')];
    expect(groups.map((g) => [g.getAttribute('aria-label'), g.textContent])).toEqual([
      ['Peach to grape', '2 ribbons between peach and grape in view'],
    ]);
  });
});
