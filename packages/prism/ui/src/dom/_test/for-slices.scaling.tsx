// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Scheduler } from '@molstar/core/task/util/scheduler';
import { afterAll, beforeAll, it, vi } from 'vitest';
import { checkLinear, scaling } from '../../_test/scaling.js';
import { observableState } from '../../reactive/observable-state.js';
import { For, timeSlicing } from '../control-flow.js';
import { render } from '../render.js';
import { flush } from '../scope.js';

// A long, non-virtualized list's first render over many time slices grows linearly (`test:scaling`,
// ui/src/_test/scaling.ts): each slice makes its own new entries, not a reconcile of the whole list again (a walk of
// n items per slice: n² / slice for the whole).

// The continuations, run here at once (Mol*'s `setImmediate` would run them a task later).
const continuations: (() => void)[] = [];
beforeAll(() => {
  timeSlicing(true);
  vi.spyOn(Scheduler, 'setImmediate').mockImplementation((fn: () => void) => {
    continuations.push(fn);
    return 0;
  });
});
afterAll(() => {
  timeSlicing('auto');
  vi.restoreAllMocks();
});

// A clock that moves 20 ms a reading inside the run: a slice makes its first entries (64) and is spent, so n entries
// are n / 64 slices at any machine's speed. The real clock is back when the run ends (the harness times it).
const realNow = performance.now.bind(performance);
let virtual = 0;
function withSlicedClock(fn: () => void) {
  performance.now = () => (virtual += 20);
  try {
    fn();
  } finally {
    performance.now = realNow;
  }
}

it('a long list rendered over many slices, each making only its new entries', () => {
  const s = scaling(
    {
      setup: (n) => {
        const items = observableState<readonly number[]>([]);
        const box = document.createElement('div');
        document.body.append(box);
        const dispose = render(() => <For each={items}>{(i) => <div>{`${i}`}</div>}</For>, box);
        flush();
        return { n, items, box, dispose };
      },
      run: ({ n, items, box }) => {
        withSlicedClock(() => {
          items.set(Array.from({ length: n }, (_, i) => i));
          flush();
          for (let i = 0; i < continuations.length; i++) continuations[i]();
          continuations.length = 0;
        });
        if (box.childElementCount !== n) throw new Error(`rendered ${box.childElementCount} of ${n}`);
      },
      teardown: ({ box, dispose }) => {
        dispose();
        box.remove();
      },
    },
    4000,
  );
  checkLinear('For: a long list over many slices', s);
});
