/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createPendingWork } from '../pending-work.js';

function deferred() {
  let resolve!: () => void;
  let reject!: (e: unknown) => void;
  const promise = new Promise<void>((res, rej) => {
    resolve = res;
    reject = rej;
  });
  return { promise, resolve, reject };
}

describe('pending work', () => {
  it('is ready without tracked work', async () => {
    const pending = createPendingWork();
    expect(pending.isReady()).toBe(true);
    await expect(pending.whenReady()).resolves.toBeUndefined();
  });

  it('waits for all tracked work, including work tracked while waiting', async () => {
    const pending = createPendingWork();
    const a = deferred();
    const b = deferred();
    pending.track(a.promise, 'a');
    let ready = false;
    const waiting = pending.whenReady().then(() => {
      ready = true;
    });

    pending.track(b.promise, 'b');
    a.resolve();
    await a.promise;
    await Promise.resolve();
    expect(ready).toBe(false);
    expect(pending.isReady()).toBe(false);

    b.resolve();
    await waiting;
    expect(ready).toBe(true);
    expect(pending.isReady()).toBe(true);
  });

  it('emits settled per tracked promise', async () => {
    const pending = createPendingWork();
    let count = 0;
    pending.settled.subscribe(() => count++);
    const a = deferred();
    const b = deferred();
    pending.track(a.promise, 'a');
    pending.track(b.promise, 'b');
    a.resolve();
    b.resolve();
    await pending.whenReady();
    expect(count).toBe(2);
  });

  it('reports rejections and keeps going', async () => {
    const error = jest.spyOn(console, 'error').mockImplementation(() => {});
    const pending = createPendingWork();
    const a = deferred();
    pending.track(a.promise, 'a').catch(() => {});
    a.reject(new Error('failed'));
    await pending.whenReady();
    expect(pending.isReady()).toBe(true);
    expect(error).toHaveBeenCalledWith("'a' failed", expect.any(Error));
    error.mockRestore();
  });

  it('counts tracked work, also after it settled', async () => {
    const pending = createPendingWork();
    expect(pending.trackedCount).toBe(0);
    const a = deferred();
    pending.track(a.promise, 'a');
    a.resolve();
    await pending.whenReady();
    expect(pending.isReady()).toBe(true);
    expect(pending.trackedCount).toBe(1);
  });

  it('releases waiters on clear', async () => {
    const pending = createPendingWork();
    pending.track(new Promise(() => {}), 'never');
    const waiting = pending.whenReady();
    pending.clear();
    await expect(waiting).resolves.toBeUndefined();
    expect(pending.isReady()).toBe(true);
  });
});
