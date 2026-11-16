/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export const ReactiveConfig: {
  /**
   * Told of an error nothing handles: a listener or a subscriber's `next` that throws, an `error` with no handler, a
   * `derive` that throws while its listeners are told. Called at once. Unset, the error is rethrown in a `setTimeout`
   * (as rxjs), so it surfaces as uncaught without stopping the producer or the other listeners.
   */
  onUnhandledError: ((err: any) => void) | null;
} = { onUnhandledError: null };

export function reportUnhandledError(err: any): void {
  const handler = ReactiveConfig.onUnhandledError;
  if (handler) {
    handler(err);
    return;
  }
  setTimeout(() => {
    throw err;
  });
}
