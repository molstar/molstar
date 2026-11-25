/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick allocation tier (`test:perf`) of the hover label (stage §3.1): beside the pointer, a move and the frame that
// places it (one read, one write of cached offsets) allocate nothing, nor does the layer's box read again after a resize
// or a scroll; the hide and show swap neither. The label and its
// layer are stand-ins (a DOM's rect and style writes are the browser's, not this code's).

import { quickCase } from '@molstar/test-harness/quick';
import { expect, it } from 'vitest';
import { createRoot, DomHost } from '../../dom/index.js';
import { type FollowedLabel, type FollowedLayer, HoverFollower } from '../hover-label.js';

/** The page's one frame tick, run by the steps. */
let tick: ((t: number) => void) | undefined;
DomHost.requestAnimationFrame = (f) => {
  tick = f;
  return 1;
};
let now = 0;
const runFrame = () => {
  const f = tick;
  tick = undefined;
  now += 16;
  f?.(now);
};

/** Writes, counted (so nothing is optimized away). */
let writes = 0;
const rect = { left: 10, top: 20 };
const layer: FollowedLayer = { getBoundingClientRect: () => rect, clientWidth: 1200, clientHeight: 800 };
const label: FollowedLabel = {
  offsetWidth: 140,
  offsetHeight: 40,
  classList: {
    toggle: () => {
      writes++;
      return true;
    },
  },
  style: {
    setProperty: () => {
      writes++;
    },
  },
  setAttribute: () => {
    writes++;
  },
  removeAttribute: () => {
    writes++;
  },
};

const cases: Record<string, () => (i: number) => void> = {
  /** Beside the pointer: a move a frame, and the frame placing it (flipped at the edge now and then). */
  cursorMove() {
    let follower!: HoverFollower;
    createRoot(() => {
      follower = new HoverFollower(label, () => layer, false, 'out');
    });
    follower.setContent(true);
    return (i) => {
      follower.move(100 + (i % 1100), 100 + (i % 700), false, true);
      runFrame();
    };
  },
  /** Beside the pointer, the stage resized or scrolled now and then: the layer's box read again in the next frame. */
  cursorResized() {
    let follower!: HoverFollower;
    createRoot(() => {
      follower = new HoverFollower(label, () => layer, false, 'out');
    });
    follower.setContent(true);
    return (i) => {
      if ((i & 15) === 0) follower.invalidate();
      follower.move(100 + (i % 1100), 100 + (i % 700), false, true);
      runFrame();
    };
  },
  /** In the corner: a move a frame, held and released now and then (hidden and shown at once). */
  cornerMove() {
    let follower!: HoverFollower;
    createRoot(() => {
      follower = new HoverFollower(label, () => layer, true, 'out');
    });
    follower.setContent(true);
    return (i) => {
      follower.move(100 + (i % 500), 100, (i & 15) === 0, true);
      runFrame();
    };
  },
};

for (const [name, make] of Object.entries(cases)) {
  it(name, ({ task }) => {
    quickCase(task, `hover label: ${name}`, make());
    expect(writes).toBeGreaterThan(0);
  });
}
