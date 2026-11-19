/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { abortError, type FormatRuntime, type OpenedSource, type RegionRequest } from './protocol.js';
import type { SeqData } from './seq-format.js';
import { openSeq } from './seq-runtime.js';
import { RealClock, type StoreClock } from './store.js';

// Asynchronous sources (_spec/README.md "Data", Asynchronous sources): a source whose every answer comes after a delay
// on a clock, the store's `StoreClock` (real timers by default, a fake one in tests, so a test never waits). It is the
// synthetic by-region source of tests and stories: `seq` data answered as a remote source would answer it.

export interface DelayOptions {
  /** The clock the delay runs on: real timers by default. */
  readonly clock?: StoreClock;
  /** Each answer's delay in ms, or by its request (0: as soon as the inner source answers, no timer). */
  readonly delayMs?: number | ((request: RegionRequest) => number);
}

/** A source answering as `inner` does, each answer after `delayMs` on `clock`; an abort rejects at once. */
export function delayedSource(inner: OpenedSource, options: DelayOptions = {}): OpenedSource {
  const clock = options.clock ?? RealClock;
  const delayOf = options.delayMs ?? 0;
  return {
    coverage: inner.coverage,
    ...(inner.issues && { issues: inner.issues }),
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      const delay = typeof delayOf === 'number' ? delayOf : delayOf(request);
      if (delay <= 0) return inner.read(request, signal);
      return new Promise((resolve, reject) => {
        const onAbort = () => {
          clock.clearTimeout(timer);
          reject(abortError());
        };
        const timer = clock.setTimeout(() => {
          signal.removeEventListener('abort', onAbort);
          inner.read(request, signal).then(resolve, reject);
        }, delay);
        signal.addEventListener('abort', onAbort, { once: true });
      });
    },
    dispose() {
      inner.dispose();
    },
  };
}

/** The synthetic asynchronous source: `seq` data (its coverage and request/answer pairs) answered after a delay. */
export function asyncSeqSource(data: SeqData, options: DelayOptions = {}): OpenedSource {
  return delayedSource(openSeq(data), options);
}

/**
 * A runtime that opens what `inner` opens and delays every answer: the application's choice of format names
 * (`formats`, the inner's by default), so a spec names the delayed one (`seq-async`) apart from the immediate one.
 */
export function delayedRuntime(
  inner: FormatRuntime,
  options: DelayOptions & { readonly name?: string; readonly formats?: readonly string[] } = {},
): FormatRuntime {
  return {
    name: options.name ?? `${inner.name}-delayed`,
    formats: options.formats ?? inner.formats,
    async open(source, ctx) {
      return delayedSource(await inner.open(source, ctx), options);
    },
  };
}
