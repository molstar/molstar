// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { lostReleases, stubPointerCapture } from '@molstar/prism-ui/testing/gestures';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import type { PaeMatrix } from '../pae/data.js';
import { type PaeBox, PaePlot } from '../pae/plot.js';

// The plot's box drag attaches once (prism AGENTS.md, gestures; `drag`, dom/scoped.ts): pressed, moved and left without
// a release N times, it leaves no listener behind, and the next drag makes a box.

beforeAll(stubPointerCapture);
afterEach(() => {
  vi.restoreAllMocks();
  document.body.innerHTML = '';
});

/** 4 × 4, one chain. */
const matrix: PaeMatrix = {
  name: 'PAE',
  size: 4,
  values: new Float32Array(16).map((_, k) => Math.abs((k >> 2) - (k & 3))),
  max: 3,
  residues: [1, 2, 3, 4].map((seq) => ({ chain: 'A', seq, insCode: '', name: 'ALA' })),
};

describe('the PAE plot', () => {
  it('a box drag: no listener left by lost releases, and the next drag makes a box', () => {
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockReturnValue(null);
    const boxes: (PaeBox | undefined)[] = [];
    const host = document.createElement('div');
    document.body.append(host);
    const dispose = render(() => <PaePlot matrix={matrix} onBox={(b) => boxes.push(b)} />, host);
    const plot = host.querySelector<HTMLElement>('[role="img"]')!;
    // jsdom lays nothing out: the plot's box as a browser would have it, 100 px a cell.
    plot.getBoundingClientRect = () => ({ left: 0, top: 0, width: 400, height: 400 }) as DOMRect;
    const r = lostReleases({
      target: () => plot,
      from: { x: 50, y: 50 },
      to: { x: 250, y: 350 },
      reset: () => (boxes.length = 0),
      worked: () => boxes.length === 1 && JSON.stringify(boxes[0]) === '{"i0":0,"i1":3,"j0":0,"j1":2}',
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
    dispose();
  });
});
