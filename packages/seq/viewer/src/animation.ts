/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The animation port's adapters (_spec/README.md "Ports", Animation), an entry of its own
// (`@molstar/seqx-viewer/animation`) for the composition roots: the page's frames (`frameAnimation`, the default: a
// plain viewer animates with no host code) and a fake clock stepped by hand (`headlessAnimation`: tests and export).

import { cancelFrameCallback, DomHost, requestFrameCallback } from '@molstar/ui/dom';
import type { TransitionPlan } from '@molstar/seqx-core/story';
import type {
  AnimationRequest,
  AnimationRequestOptions,
  SeqAnimation,
  SeqAnimationParticipant,
  SeqAnimationPort,
} from '@molstar/seqx-viewer/ports';

interface Running {
  readonly plan: TransitionPlan;
  readonly options: AnimationRequestOptions;
  readonly t0: number;
  live: boolean;
}

/** The time a request reached, [0, 1]. */
const timeOf = (r: Running, now: number) => {
  const d = r.options.duration;
  return d > 0 ? Math.min(1, Math.max(0, (now - r.t0) / d)) : 1;
};

/**
 * The default adapter: requested plans run on the page's frames (`@molstar/ui/dom`'s frame callbacks, on `DomHost`:
 * its clock is `DomHost.now`), timed from the request, the last frame at `t` 1.
 */
export function frameAnimation(): SeqAnimationPort {
  return (participant: SeqAnimationParticipant): SeqAnimation => {
    const running = new Set<Running>();
    let frame: number | undefined;
    const tick = () => {
      frame = undefined;
      const now = DomHost.now();
      for (const r of [...running]) {
        if (!r.live) continue;
        const t = timeOf(r, now);
        participant.progress(r.plan, t);
        if (t >= 1) {
          r.live = false;
          running.delete(r);
          r.options.done?.();
        }
      }
      // a request made from a `progress` or a `done` above already asked for the next tick
      if (running.size) frame ??= requestFrameCallback(tick);
    };
    return {
      request(plan, options): AnimationRequest {
        const r: Running = { plan, options, t0: DomHost.now(), live: true };
        running.add(r);
        frame ??= requestFrameCallback(tick);
        return {
          cancel() {
            r.live = false;
            running.delete(r);
          },
        };
      },
      dispose() {
        running.clear();
        if (frame !== undefined) cancelFrameCallback(frame);
        frame = undefined;
      },
    };
  };
}

/** A fake clock the application steps: what a headless twin and an export run on. */
export interface HeadlessClock {
  /** The clock's time, ms. */
  readonly now: number;
  /** Plans running. */
  readonly running: number;
  /** Moves the clock by `ms`: each running plan is shown at its time and drawn. */
  tick(ms: number): void;
  /** Ticks by `stepMs` until no plan runs (at most `maxSteps`); returns the ticks taken. */
  runAll(stepMs: number, maxSteps?: number): number;
}

/** The headless adapter: a clock that moves only when ticked, each tick drawn (`participant.draw`). */
export function headlessAnimation(): { readonly port: SeqAnimationPort; readonly clock: HeadlessClock } {
  let now = 0;
  const running = new Set<Running>();
  const owners = new Map<Running, SeqAnimationParticipant>();
  const clock: HeadlessClock = {
    get now() {
      return now;
    },
    get running() {
      return running.size;
    },
    tick(ms) {
      now += ms;
      const drawn = new Set<SeqAnimationParticipant>();
      for (const r of [...running]) {
        if (!r.live) continue;
        const p = owners.get(r)!;
        const t = timeOf(r, now);
        p.progress(r.plan, t);
        drawn.add(p);
        if (t >= 1) {
          r.live = false;
          running.delete(r);
          owners.delete(r);
          r.options.done?.();
        }
      }
      for (const p of drawn) p.draw();
    },
    runAll(stepMs, maxSteps = 10_000) {
      let n = 0;
      while (running.size && n < maxSteps) {
        clock.tick(stepMs);
        n++;
      }
      return n;
    },
  };
  const port: SeqAnimationPort = (participant) => {
    return {
      request(plan, options) {
        const r: Running = { plan, options, t0: now, live: true };
        running.add(r);
        owners.set(r, participant);
        return {
          cancel() {
            r.live = false;
            running.delete(r);
            owners.delete(r);
          },
        };
      },
      dispose() {
        for (const [r, p] of owners)
          if (p === participant) {
            running.delete(r);
            owners.delete(r);
          }
      },
    };
  };
  return { port, clock };
}
