/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewerApi } from '../api/area.js';

// Changes at pointer rate as undo steps (a drag, keys held, a color field dragged across): the changes in quick
// succession are one `api.group`, the latest winning per key; what is not applied yet is shown already (`overlay`).
// While the gesture is held the viewer gets previews spaced out (`throttleMs`, at least as long as the last took: an
// isosurface or a direct volume is heavy, and the pointer must stay smooth); let go, the latest applies at once. A gesture ends on `settle` (the pointer let go) or after a moment's quiet.
// The map levels and the transfer function use it (sources/volume-levels.ts, sources/volume-transfer.ts).

export interface GestureSteps<C extends object> {
  /** A change to `key`, merged into what is pending for it; starts a step when none is open. */
  change(key: string, change: C): void;
  /** The changes so far are one step; the next starts another. */
  settle(): void;
  /** What is shown for a key over what the state says: the changes of steps still in progress, merged. */
  overlay(key: string): C | undefined;
}

/** A step in progress: its changes not applied yet, the ones shown until it ends, and its loop's wake-up. */
interface Step<C> {
  readonly pending: Map<string, C>;
  readonly shown: Map<string, C>;
  wake?: () => void;
  closed: boolean;
}

/** After this long without a change, the changes so far are one step. */
const QuietMs = 600;
/** Previews while a gesture is held: at most this often. */
const ThrottleMs = 250;

/** Waits for a change or the gesture's end (`wake`), or `ms` if given. */
function sleep(step: { wake?: () => void }, ms?: number): Promise<void> {
  return new Promise<void>((resolve) => {
    const timer = ms === undefined ? undefined : setTimeout(done, ms);
    function done() {
      clearTimeout(timer);
      if (step.wake === done) step.wake = undefined;
      resolve();
    }
    step.wake = done;
  });
}

export function gestureSteps<C extends object>(
  api: ViewerApi,
  options: {
    /** The step's label in the history. */
    readonly label: string;
    /** Applies a key's merged change, as part of the step (`api` is the group's). */
    readonly apply: (api: ViewerApi, key: string, change: C) => Promise<unknown>;
    /** What is shown changed (a change came, or a step ended and its overlay went). */
    readonly changed: () => void;
    readonly quietMs?: number;
    /** While a gesture is held, previews at most this often (and no more often than the last one took). */
    readonly throttleMs?: number;
  },
): GestureSteps<C> {
  const quietMs = options.quietMs ?? QuietMs;
  const throttleMs = options.throttleMs ?? ThrottleMs;
  const steps: Step<C>[] = [];
  let open: Step<C> | undefined;
  let quiet: ReturnType<typeof setTimeout> | undefined;

  const run = (step: Step<C>) =>
    api
      .group(
        options.label,
        async (g) => {
          // When the last preview landed and how long it took: the next waits at least as long (and `throttleMs`), so
          // a heavy update leaves the page time to follow the pointer.
          let last = Number.NEGATIVE_INFINITY;
          let took = 0;
          for (;;) {
            const next = step.pending.entries().next();
            if (!next.done) {
              // While the gesture is held, previews are spaced out; let go, the latest goes at once (what waited for
              // its turn is replaced by it, merged per key).
              const wait = step.closed ? 0 : Math.max(throttleMs, took) - (performance.now() - last);
              if (wait > 0) {
                await sleep(step, wait);
                continue;
              }
              const [key, change] = next.value;
              step.pending.delete(key);
              const t0 = performance.now();
              await options.apply(g, key, change);
              last = performance.now();
              took = last - t0;
              continue;
            }
            if (step.closed) return;
            await sleep(step);
          }
        },
        { onFailure: 'continue' },
      )
      .catch(() => undefined)
      .finally(() => {
        steps.splice(steps.indexOf(step), 1);
        options.changed();
      });

  const settle = () => {
    clearTimeout(quiet);
    const step = open;
    open = undefined;
    if (!step) return;
    step.closed = true;
    step.wake?.();
  };

  return {
    change(key, change) {
      let step = open;
      if (!step) {
        step = open = { pending: new Map(), shown: new Map(), closed: false };
        steps.push(step);
        void run(step);
      }
      step.pending.set(key, { ...step.pending.get(key), ...change });
      step.shown.set(key, { ...step.shown.get(key), ...change });
      step.wake?.();
      options.changed();
      clearTimeout(quiet);
      quiet = setTimeout(settle, quietMs);
    },
    settle,
    overlay(key) {
      let out: C | undefined;
      for (const step of steps) {
        const c = step.shown.get(key);
        if (c) out = { ...out, ...c } as C;
      }
      return out;
    },
  };
}
