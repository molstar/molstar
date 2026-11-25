/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories' SVG twin waits for the canvas (follow-svg.ts): a change's export runs in a task of its own after the frame
// that draws the change, never in the input's task nor in the microtasks after its frame. Composed as the Copies story
// composes it (the 200 copies on the host's stack, the twin following), frames queued as a page's are and flushed by
// hand: a split, then a merge, each to drawn with its microtasks; the twin lands a task after the next frame, once. Heavy: the split
// to drawn under 16 ms with its microtasks (it was 25 to 50 ms while the twin exported at once: a second layout).

import { ComplexViews } from '@molstar/seqx-examples';
import { copiesScale } from '@molstar/seqx-examples/complex/blocks';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { DomHost } from '@molstar/ui/dom';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, describe, expect, vi } from 'vitest';
import { composeEditing, type EditingSession } from '../../features/editing/compose.js';
import { type FollowedSvg, followSvg } from '../../follow-svg.js';
import { setCopies } from '../complex/compose.js';

const raf = DomHost.requestAnimationFrame;
const frames: ((t: number) => void)[] = [];
let session: EditingSession | undefined;
let twin: FollowedSvg | undefined;

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  frames.length = 0;
  twin?.dispose();
  twin = undefined;
  session?.dispose();
  session = undefined;
});

/** The page's frame: the queued callbacks run (a callback that asks for another waits for the next). */
const frame = () => {
  for (const f of frames.splice(0)) f(performance.now());
};
/** The microtasks after a frame drained (an awaited chain of them), as the measurements drain them. */
const drain = async () => {
  for (let i = 0; i < 30; i++) await null;
};
/** A macrotask. */
const task = () => new Promise((r) => setTimeout(r, 0));

describe('the SVG twin waits for the canvas', () => {
  itWithin(
    'the 200 copies split and merged: each drawn with its microtasks before the twin exports, under 16 ms',
    async (within) => {
      DomHost.requestAnimationFrame = (fn) => void frames.push(fn);
      const base = ComplexViews.find((v) => v.id === 'blocks')!;
      const build = () => copiesScale(200);
      const s = composeEditing(
        { ...base, id: 'copies-scale', build },
        { ports: { assets: () => nodeFixtureAssets() }, delayMs: 0 },
      );
      session = s;
      const v = s.viewer;
      v.resize(1440, 900);
      const host = { innerHTML: '' } as unknown as HTMLElement;
      twin = followSvg(v, host, { prefix: 't' });
      const exports = vi.spyOn(v, 'renderSvg');
      await v.load(build());
      await v.settle();
      for (let i = 0; i < 20; i++) {
        await new Promise((r) => setTimeout(r, 20));
        frame();
      }
      await twin.draw();
      // warm: a split and a merge, each landed
      for (const copies of ['each', 'merge'] as const) {
        void setCopies(s, 'copies', copies);
        frame();
        await drain();
        await task();
        await new Promise((r) => setTimeout(r, 20));
        frame();
      }
      for (const copies of ['each', 'merge'] as const) {
        exports.mockClear();
        const before = host.innerHTML;
        // the input, its frame and the microtasks after it: the edit drawn, the twin not started
        await within(16, async () => {
          void setCopies(s, 'copies', copies);
          frame();
          await drain();
        });
        expect(v.frame!.layout.sections[0].section.blocks!.length).toBe(copies === 'each' ? 200 : 1);
        expect(exports).not.toHaveBeenCalled();
        // the edit's step lands in those microtasks (its revision): the twin waits for the frame after it, then a task
        frame();
        expect(exports).not.toHaveBeenCalled();
        await task();
        expect(exports).toHaveBeenCalledTimes(1);
        await exports.mock.results[0].value;
        await drain();
        expect(host.innerHTML).not.toBe(before);
        await new Promise((r) => setTimeout(r, 20));
        frame();
        await task();
      }
      expect(exports).toHaveBeenCalledTimes(1);
    },
    60_000,
  );
});
