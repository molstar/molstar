/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PlacePoint, type TextPlacement } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, describe, expect, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// The text viewer's performance representatives (performance.md, text layout): on the p53 text figure in a box that
// scrolls, a hover (the input's move through the placement, the hover's redraw) under 2 ms; a drag-select's step (the
// unit under the pointer in reading order, the selection's redraw) and a cursor's key step a line down (the layout's
// step, the reveal, the redraw) under 4 ms each. The bounds are the heavy lane's (`itWithin`, best of three); the
// quick lane runs them unbounded. The clock of the bound is the real one: the timers alone are faked.

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open() {
  viewer = await openExample(P53Views.find((v) => v.id === 'text')!, { gene: '', realClock: true });
  viewer.input.mac = false;
  viewer.resize(360);
  viewer.draw();
  return viewer;
}

/** Where a display unit is drawn on the page. */
function at(v: SeqViewer, u: number) {
  const p = v.frame!.layout.sections[0].placement as TextPlacement;
  const pt = new PlacePoint();
  p.place(u, pt);
  return { x: p.strips.axes[pt.strip].xOfUnit(u + 0.5), y: pt.y, upl: p.geometry!.unitsPerLine };
}

describe('the text viewer’s time', () => {
  itWithin('a hover over a text', async (within) => {
    const v = await open();
    const a = at(v, 3);
    const b = at(v, a.upl + 7);
    v.input.move(a.x, a.y);
    v.draw();
    const t0 = performance.now();
    within(2, () => {
      v.input.move(b.x, b.y);
      v.draw();
    });
    console.info(`text hover: ${(performance.now() - t0).toFixed(3)} ms`);
    expect(v.hover.current.position).toBe(a.upl + 7);
  });

  itWithin('a drag-select step across lines', async (within) => {
    const v = await open();
    const a = at(v, 3);
    const b = at(v, 2 * a.upl + 5);
    v.input.down(a.x, a.y);
    v.input.move(a.x + 20, a.y);
    v.draw();
    const t0 = performance.now();
    within(4, () => {
      v.input.move(b.x, b.y);
      v.draw();
    });
    console.info(`text drag-select step: ${(performance.now() - t0).toFixed(3)} ms`);
    v.input.up();
    expect([...v.selection.value.value.shown.parts[0].ranges]).toEqual([3, 2 * a.upl + 6]);
  });

  itWithin('a cursor’s key step a line down', async (within) => {
    const v = await open();
    const a = at(v, 3);
    v.input.click(a.x, a.y);
    v.draw();
    const t0 = performance.now();
    within(4, () => {
      v.input.key({ key: 'ArrowDown' });
      v.draw();
    });
    console.info(`text key step: ${(performance.now() - t0).toFixed(3)} ms`);
    expect(v.selection.cursor.value?.position).toBe(3 + a.upl);
  });
});
