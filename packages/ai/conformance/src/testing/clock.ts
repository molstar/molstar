/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Clock } from '@molstar/ai-harness';

interface Sleeper {
  readonly at: number;
  readonly order: number;
  resolve(): void;
  reject(reason: unknown): void;
  cleanup(): void;
}

/** Lets pending promises settle: a number of rounds of the event loop, then a timer's turn. */
export async function settle(rounds = 12): Promise<void> {
  for (let i = 0; i < rounds; i += 1) await new Promise<void>((resolve) => setImmediate(resolve));
  await new Promise<void>((resolve) => setTimeout(resolve, 0));
}

/**
 * Virtual time. Nothing waits for real: `sleep` parks until time is advanced, either explicitly
 * (`advance`) or by `runUntil`, which moves time to the next sleeper whenever nothing else can run.
 */
export class SimulatedClock implements Clock {
  private time: number;
  private sleepers: Sleeper[] = [];
  private counter = 0;
  /** Calls to `now` and `sleep`: work in progress, so `runUntil` does not move time under it. */
  private activity = 0;

  constructor(start = 1_000_000) {
    this.time = start;
  }

  now(): number {
    this.activity += 1;
    return this.time;
  }

  sleep(ms: number, signal?: AbortSignal): Promise<void> {
    this.activity += 1;
    return new Promise<void>((resolve, reject) => {
      if (signal?.aborted) {
        reject(signal.reason);
        return;
      }
      const onAbort = () => {
        this.sleepers = this.sleepers.filter((s) => s !== sleeper);
        reject(signal?.reason);
      };
      const sleeper: Sleeper = {
        at: this.time + Math.max(0, ms),
        order: this.counter++,
        resolve,
        reject,
        cleanup: () => signal?.removeEventListener('abort', onAbort),
      };
      signal?.addEventListener('abort', onAbort, { once: true });
      this.sleepers.push(sleeper);
      this.sleepers.sort((a, b) => a.at - b.at || a.order - b.order);
    });
  }

  /** The time of the next sleeper, if any. */
  next(): number | undefined {
    return this.sleepers[0]?.at;
  }

  /** Moves time forward by `ms`, waking every sleeper due on the way, in order. */
  async advance(ms: number): Promise<void> {
    const target = this.time + ms;
    for (;;) {
      const first = this.sleepers[0];
      if (first === undefined || first.at > target) break;
      this.sleepers.shift();
      this.time = Math.max(this.time, first.at);
      first.cleanup();
      first.resolve();
      await settle();
    }
    this.time = target;
  }

  /**
   * Runs `promise` to completion in virtual time: settles what can run, and when nothing can, moves
   * time to the next sleeper. Fails if nothing can run and nothing sleeps (a deadlock).
   */
  async runUntil<T>(promise: Promise<T>, maxSteps = 100_000): Promise<T> {
    let done = false;
    let value: T | undefined;
    let error: unknown;
    let failed = false;
    promise.then(
      (result) => {
        done = true;
        value = result;
      },
      (reason: unknown) => {
        done = true;
        failed = true;
        error = reason;
      },
    );
    for (let step = 0; step < maxSteps; step += 1) {
      // Quiescence: settle until a whole settle passes with no clock activity.
      for (let quiet = 0, rounds = 0; quiet < 2 && rounds < 10_000 && !done; rounds += 1) {
        const seen = this.activity;
        await settle();
        quiet = this.activity === seen ? quiet + 1 : 0;
      }
      if (done) break;
      const first = this.sleepers[0];
      if (first === undefined) throw new Error('SimulatedClock.runUntil: deadlock (nothing runs, nothing sleeps)');
      await this.advance(first.at - this.time);
    }
    if (!done) throw new Error('SimulatedClock.runUntil: too many steps');
    if (failed) throw error;
    return value as T;
  }
}
