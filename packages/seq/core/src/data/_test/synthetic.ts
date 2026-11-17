/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import { abortError, type Coverage, type OpenedSource, type RegionAnswer, type RegionRequest } from '../protocol.js';
import type { StoreClock } from '../store.js';

// Test doubles: a fake clock (timers run by `advance`, never real ones) and a synthetic asynchronous source.

interface Timer {
  at: number;
  seq: number;
  fn: () => void;
}

/** Lets promise callbacks run: the store's `then` after a source resolves. */
export async function flush(): Promise<void> {
  for (let i = 0; i < 8; i++) await Promise.resolve();
}

export class FakeClock implements StoreClock {
  private _now = 0;
  private _seq = 0;
  private readonly _timers = new Map<number, Timer>();

  now(): number {
    return this._now;
  }

  setTimeout(fn: () => void, ms: number): unknown {
    const seq = ++this._seq;
    this._timers.set(seq, { at: this._now + ms, seq, fn });
    return seq;
  }

  clearTimeout(handle: unknown): void {
    this._timers.delete(handle as number);
  }

  get pending(): number {
    return this._timers.size;
  }

  /** Runs the timers due within `ms`, in time order, letting promises settle after each. */
  async advance(ms: number): Promise<void> {
    const until = this._now + ms;
    await flush();
    for (;;) {
      let next: Timer | undefined;
      for (const t of this._timers.values())
        if (t.at <= until && (!next || t.at < next.at || (t.at === next.at && t.seq < next.seq))) next = t;
      if (!next) break;
      this._timers.delete(next.seq);
      this._now = next.at;
      next.fn();
      await flush();
    }
    this._now = until;
  }
}

export interface SyntheticRead {
  readonly request: RegionRequest;
  aborted: boolean;
  answered: boolean;
}

/**
 * A source answering after `delayMs` on the fake clock: by `answer(request)`, by default a ranges table with one
 * range per 100 positions of the region. It records its reads and honors aborts.
 */
export class SyntheticSource implements OpenedSource {
  readonly reads: SyntheticRead[] = [];
  disposed = false;

  constructor(
    readonly coverage: Coverage,
    private readonly clock: FakeClock,
    readonly options: { delayMs?: number; answer?: (r: RegionRequest, n: number) => RegionAnswer | Error } = {},
  ) {}

  read(request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer> {
    const read: SyntheticRead = { request, aborted: false, answered: false };
    const n = this.reads.push(read);
    return new Promise((resolve, reject) => {
      const timer = this.clock.setTimeout(() => {
        read.answered = true;
        const a = (this.options.answer ?? rangesAnswer)(request, n);
        if (a instanceof Error) reject(a);
        else resolve(a);
      }, this.options.delayMs ?? 100);
      signal.addEventListener('abort', () => {
        read.aborted = true;
        this.clock.clearTimeout(timer);
        reject(abortError());
      });
    });
  }

  dispose(): void {
    this.disposed = true;
  }

  /** Reads still waiting for their answer. */
  get open(): SyntheticRead[] {
    return this.reads.filter((r) => !r.aborted && !r.answered);
  }
}

export function rangesAnswer(r: RegionRequest): RegionAnswer {
  const starts: number[] = [];
  for (let p = Math.ceil(r.start / 100) * 100; p < r.end; p += 100) starts.push(p);
  return {
    ok: true,
    tables: {
      ranges: {
        form: 'rows',
        start: Int32Array.from(starts),
        end: Int32Array.from(starts, (s) => s + 10),
        id: categoryColumn(starts.map((s) => `r${s}`)),
      },
    },
  };
}

/** One system `chr1` of `length`, ranges, at levels 1 and 10. */
export function syntheticCoverage(length = 1_000_000): Coverage {
  return {
    systems: [{ name: 'chr1', length, aliases: ['1'] }],
    shapes: [{ system: 'chr1', shape: 'ranges' }],
    levels: [1, 10],
  };
}
