/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { HarnessObservable } from '../contracts/observe.js';

/** A hot subject: no replay; a subscriber that throws is isolated (harness.md §13). */
export class Subject<T> implements HarnessObservable<T> {
  private readonly subscribers = new Set<{ next(value: T): void; error?(error: unknown): void; done?(): void }>();

  subscribe(onNext: (value: T) => void, onError?: (error: unknown) => void, onComplete?: () => void): () => void {
    const subscriber = {
      next: onNext,
      ...(onError !== undefined ? { error: onError } : {}),
      ...(onComplete !== undefined ? { done: onComplete } : {}),
    };
    this.subscribers.add(subscriber);
    return () => {
      this.subscribers.delete(subscriber);
    };
  }

  next(value: T): void {
    for (const subscriber of [...this.subscribers]) {
      try {
        subscriber.next(value);
      } catch (error) {
        try {
          subscriber.error?.(error);
        } catch {
          // An observer cannot affect execution.
        }
      }
    }
  }

  complete(): void {
    for (const subscriber of [...this.subscribers]) {
      try {
        subscriber.done?.();
      } catch {
        // Ignored, as above.
      }
    }
    this.subscribers.clear();
  }
}
