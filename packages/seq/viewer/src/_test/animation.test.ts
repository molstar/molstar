/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The animation port's default adapter (_spec/README.md "Ports", Animation): one frame callback at a time, however
// requests are made from inside its tick (a plan's `done` starting the next, a `progress` requesting another).

import { State } from '@molstar/core/reactive';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { plan } from '@molstar/seqx-core/story';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { frameAnimation } from '../animation.js';
import type { SeqAnimationParticipant } from '../ports/index.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
afterEach(() => {
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

function participant(onProgress: (t: number) => void = () => {}): SeqAnimationParticipant {
  return {
    progress: (_p, t) => onProgress(t),
    show: () => {},
    clear: () => {},
    draw: () => {},
    drawn: new State(0),
    durationOf: (p) => p.duration,
  };
}

describe('frameAnimation', () => {
  it('a request made in a plan’s `done` asks for one tick, not two', () => {
    vi.useFakeTimers({ toFake: ['Date', 'performance'] });
    frames = installFrames();
    const root = createSeqSpec().build();
    const p = plan({ root }, { root }, { duration: 100 });
    const anim = frameAnimation()(participant());
    let chained = false;
    anim.request(p, {
      duration: 100,
      done: () => {
        chained = true;
        anim.request(p, { duration: 100 });
      },
    });
    vi.advanceTimersByTime(150);
    // run the one queued tick: the plan ends, its `done` requests the next
    const before = frames.callbacks;
    expect(before).toBe(1);
    frames.flushOnce();
    expect(chained).toBe(true);
    expect(frames.callbacks).toBe(1);
    anim.dispose?.();
  });
});
