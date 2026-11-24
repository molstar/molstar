/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a drag-select re-runs, a drag kind each (the viewer's _spec, "SelectionModel and HighlightModel", Drags), on the
// views the author dragged on, composed as the stories compose them with the SVG twin beside: over many moves and a few
// frames, the pointer is hit once a frame (never a move; a frame with the pointer held where it was hit, not again),
// the selection written at most once a frame, the hover neither looked up nor changed, nothing laid out again, the SVG
// twin not exported; on release the last moves land once, the hover comes back where the pointer rests and the twin is
// exported once (a task after the release's frame), and the selection is the one a drag straight to the last point
// makes.

// @module-tag medium

import { loadGenome } from '@molstar/seqx-core/resolve';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { followSvg } from '../../../follow-svg.js';
import { DragKinds, openDragView } from './drag-kinds.js';

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  vi.useRealTimers();
});

const Frames = 8;
const MovesAFrame = 5;

describe('a drag-select re-runs its work once a frame, bounded by the screen', () => {
  for (const kind of DragKinds)
    it(kind.name, async () => {
      const v = await openDragView(kind.view);
      v.input.mac = false;
      const host = { innerHTML: '' } as HTMLElement;
      const twin = followSvg(v, host);
      const path = kind.path(v, Frames * MovesAFrame + 3);
      const hit = vi.spyOn(v, 'hit');
      const svg = vi.spyOn(v, 'renderSvg');
      let writes = 0;
      const off = v.selection.value.listen(() => writes++);
      expect(v.input.down(...path[0])).toBe(true);
      v.input.move(...path[1]);
      v.draw();
      expect(
        v.selection.value.value.shown.parts.length + v.selection.value.value.location.rects.length,
      ).toBeGreaterThan(0);
      hit.mockClear();
      svg.mockClear();
      writes = 0;
      const hover = v.hover.version.value;
      const layouts = v.counts.layouts;
      for (let f = 0; f < Frames; f++) {
        for (let m = 0; m < MovesAFrame; m++) v.input.move(...path[2 + f * MovesAFrame + m]);
        // nothing between the frames: no hit, no write
        expect(hit).toHaveBeenCalledTimes(f);
        v.draw();
      }
      expect(hit).toHaveBeenCalledTimes(Frames);
      expect(writes).toBeGreaterThan(0);
      expect(writes).toBeLessThanOrEqual(Frames);
      expect(v.hover.version.value).toBe(hover);
      expect(v.hover.current.owner).toBeUndefined();
      expect(v.counts.layouts - layouts).toBe(0);
      expect(svg).not.toHaveBeenCalled();
      // a frame with the pointer held where the last one applied it: nothing hit again
      v.input.move(...path[path.length - 2]);
      v.draw();
      expect(hit).toHaveBeenCalledTimes(Frames);
      // the last moves, with no frame after them: the release applies them, then hovers where the pointer rests
      const last = path[path.length - 1];
      v.input.move(last[0] + 0.5, last[1] + 0.5);
      v.input.move(...last);
      v.input.up();
      expect(hit).toHaveBeenCalledTimes(Frames + 2);
      v.draw();
      // the twin's export a task after that frame (follow-svg.ts), once
      expect(svg).not.toHaveBeenCalled();
      await vi.advanceTimersByTimeAsync(0);
      expect(svg).toHaveBeenCalledTimes(1);
      const dragged = v.selection.value.value.location;
      off();
      // the release's selection is a drag straight to the last point's
      v.selection.clear();
      v.input.down(...path[0]);
      v.input.move(...last);
      v.input.up();
      expect(v.selection.value.value.location).toEqual(dragged);
      twin.dispose();
      v.retire();
    });
});
