/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DomHost } from '@molstar/ui/dom';
import { vi } from 'vitest';

// A fake frame host: `DomHost`'s frames are collected and run by `flush`, its timers are counted and run by `fire`, so a
// test drives a viewer frame by frame and can tell what a retired viewer left behind.

export interface Timers {
  /** Timers set and neither fired nor cleared. */
  readonly live: number;
  /** Runs the live timers. */
  fire(): void;
}

export interface Frames {
  readonly timers: Timers;
  /** Frames asked for and not yet run. */
  readonly pending: number;
  /** Runs the asked frames, and those they ask for, until none is left. */
  flush(): void;
  restore(): void;
}

export function installFrames(): Frames {
  let queue: (() => void)[] = [];
  const live = new Map<number, () => void>();
  let id = 0;
  const raf = vi.spyOn(DomHost, 'requestAnimationFrame').mockImplementation((fn) => {
    queue.push(() => fn(DomHost.now()));
  });
  const set = vi.spyOn(DomHost, 'setTimeout').mockImplementation((fn) => {
    live.set(++id, fn);
    return id;
  });
  const clear = vi.spyOn(DomHost, 'clearTimeout').mockImplementation((h) => {
    live.delete(h as number);
  });
  const timers: Timers = {
    get live() {
      return live.size;
    },
    fire() {
      for (const [k, fn] of [...live]) {
        live.delete(k);
        fn();
      }
    },
  };
  return {
    timers,
    get pending() {
      return queue.length;
    },
    flush() {
      for (let i = 0; i < 20 && queue.length; i++) {
        const run = queue;
        queue = [];
        for (const f of run) f();
      }
    },
    restore() {
      raf.mockRestore();
      set.mockRestore();
      clear.mockRestore();
    },
  };
}
