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

/**
 * Hashes in flight. WebCrypto resolves off the event loop, so the clock cannot tell a runner that is waiting for one
 * from a runner with nothing left to do, and would move virtual time under it (a lease renewal, a deadlock that is
 * none). Every digest the harness takes goes through `crypto.subtle.digest`; this counts the ones not yet finished.
 */
let digestsInFlight = 0;
const tracked = Symbol.for('molstar.ai-conformance.digest-tracked');
const subtle = globalThis.crypto?.subtle;
if (subtle && !(subtle.digest as unknown as Record<symbol, unknown>)[tracked]) {
  const digest = subtle.digest.bind(subtle) as (...args: unknown[]) => Promise<ArrayBuffer>;
  const counted = async (...args: unknown[]) => {
    digestsInFlight += 1;
    try {
      return await digest(...args);
    } finally {
      digestsInFlight -= 1;
    }
  };
  (counted as unknown as Record<symbol, unknown>)[tracked] = true;
  (subtle as unknown as { digest: unknown }).digest = counted;
}

/**
 * Lets pending promises (including WebCrypto's, which resolve off the microtask queue) settle: a number of rounds of the
 * event loop, and then for as long as a hash is in flight (its continuation may start the next one).
 */
export async function settle(rounds = 12): Promise<void> {
  for (let i = 0; i < rounds; i += 1) await new Promise<void>((resolve) => setImmediate(resolve));
  await new Promise<void>((resolve) => setTimeout(resolve, 0));
  for (let waited = 0; waited < 100_000; ) {
    for (; digestsInFlight > 0 && waited < 100_000; waited += 1)
      await new Promise<void>((resolve) => setImmediate(resolve));
    for (let i = 0; i < 3; i += 1) await new Promise<void>((resolve) => setImmediate(resolve));
    if (digestsInFlight === 0) break;
  }
}

/** Real time `runUntil` waits for work still in flight before it calls an idle clock a deadlock. */
const DeadlockGraceMs = 5_000;

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
  /** Work the clock cannot see (storage in flight); `runUntil` waits while it is busy. */
  private busy = 0;

  /** Brackets work the clock cannot see, such as a storage operation. */
  async track<T>(work: () => Promise<T>): Promise<T> {
    this.busy += 1;
    try {
      return await work();
    } finally {
      this.busy -= 1;
      this.activity += 1;
    }
  }

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

  /**
   * Whether something is parked for at least `ms` more of virtual time: a step in a long wait (a slow model or
   * tool), not the runner's own short housekeeping sleeps. What a case that crashes mid-step waits for.
   */
  parked(ms: number): boolean {
    return this.sleepers.some((sleeper) => sleeper.at - this.time >= ms);
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
      // Quiescence: settle until a whole settle passes with no clock activity (WebCrypto work, for
      // one, resolves off the microtask queue and may outlast a single settle).
      for (let quiet = 0, rounds = 0; quiet < 2 && rounds < 10_000 && !done; rounds += 1) {
        const seen = this.activity;
        await settle();
        quiet = this.activity === seen && this.busy === 0 && digestsInFlight === 0 ? quiet + 1 : 0;
      }
      if (done) break;
      // Nothing runs and nothing sleeps: a deadlock, or work that resolves off the event loop (WebCrypto on a loaded
      // machine) that has not finished yet. Give it real time, restarted by every sign of life, before saying so.
      let first = this.sleepers[0];
      for (
        let waited = 0, seen = this.activity;
        first === undefined && !done && waited < DeadlockGraceMs;
        waited += 5
      ) {
        await new Promise<void>((resolve) => setTimeout(resolve, 5));
        await settle();
        if (this.activity !== seen) {
          seen = this.activity;
          waited = 0;
        }
        first = this.sleepers[0];
      }
      if (done) break;
      if (first === undefined) throw new Error('SimulatedClock.runUntil: deadlock (nothing runs, nothing sleeps)');
      await this.advance(first.at - this.time);
    }
    if (!done) throw new Error('SimulatedClock.runUntil: too many steps');
    if (failed) throw error;
    return value as T;
  }
}
