/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Schedulers (seqspec:viewer#the-context): every timer and frame goes through one. The root touches no DOM global, so
// the rAF scheduler takes its host (a window) as a port; the manual one fires due timers and frames in `advance`.

export interface Scheduler {
  now(): number;
  requestFrame(cb: (t: number) => void): number;
  cancelFrame(id: number): void;
  setTimer(cb: () => void, ms: number): number;
  clearTimer(id: number): void;
}

export interface ManualScheduler extends Scheduler {
  /** Moves the clock: fires due timers in due order (timers set meanwhile included), then the frames requested before. */
  advance(ms: number): void;
  readonly pendingTimers: number;
  readonly pendingFrames: number;
}

interface ManualTimer {
  readonly id: number;
  readonly due: number;
  readonly cb: () => void;
}

export function manualScheduler(start = 0): ManualScheduler {
  let time = start;
  let nextId = 1;
  let timers: ManualTimer[] = [];
  let frames = new Map<number, (t: number) => void>();
  return {
    now: () => time,
    requestFrame(cb) {
      const id = nextId++;
      frames.set(id, cb);
      return id;
    },
    cancelFrame: (id) => void frames.delete(id),
    setTimer(cb, ms) {
      const id = nextId++;
      timers.push({ id, due: time + Math.max(0, ms), cb });
      return id;
    },
    clearTimer(id) {
      timers = timers.filter((t) => t.id !== id);
    },
    advance(ms) {
      const end = time + Math.max(0, ms);
      for (;;) {
        let first: ManualTimer | null = null;
        for (const t of timers) if (t.due <= end && (!first || t.due < first.due)) first = t;
        if (!first) break;
        const fired = first;
        timers = timers.filter((t) => t !== fired);
        time = Math.max(time, fired.due);
        fired.cb();
      }
      time = end;
      const due = frames;
      frames = new Map();
      for (const cb of due.values()) cb(time);
    },
    get pendingTimers() {
      return timers.length;
    },
    get pendingFrames() {
      return frames.size;
    },
  };
}

/** What the rAF scheduler needs from its host; a browser `window` (with `performance.now`) fits. */
export interface FrameHost {
  requestAnimationFrame(cb: (t: number) => void): number;
  cancelAnimationFrame(id: number): void;
  setTimeout(cb: () => void, ms: number): number;
  clearTimeout(id: number): void;
  now(): number;
}

export function rafScheduler(host: FrameHost): Scheduler {
  return {
    now: () => host.now(),
    requestFrame: (cb) => host.requestAnimationFrame(cb),
    cancelFrame: (id) => host.cancelAnimationFrame(id),
    setTimer: (cb, ms) => host.setTimeout(cb, ms),
    clearTimer: (id) => host.clearTimeout(id),
  };
}

/**
 * The context's view of its port: records what is outstanding, so `clearAll` (the context's `dispose`) clears every
 * timer and frame, and refuses new ones after (they return -1). Each callback's wrapper is made once and kept
 *: a frame loop or a timer with a stable callback allocates nothing.
 */
export interface TrackedScheduler extends Scheduler {
  readonly outstanding: number;
  readonly closed: boolean;
  clearAll(): void;
}

type Callback = (t: number) => void;

export function trackedScheduler(port: Scheduler): TrackedScheduler {
  // What is outstanding, few at a time, swap-removed: the port's id, the wrapper, when it fires (a frame 0, a timer its
  // due time) and a sequence. A wrapper that fires settles its own entry that fires first.
  const ids: number[] = [];
  const wraps: (Callback | null)[] = [];
  const dues: number[] = [];
  const seqs: number[] = [];
  const timers: boolean[] = [];
  let n = 0;
  let seq = 0;
  let closed = false;
  const frameWraps = new WeakMap<Callback, Callback>();
  const timerWraps = new WeakMap<() => void, Callback>();
  const remove = (k: number) => {
    n--;
    ids[k] = ids[n];
    wraps[k] = wraps[n];
    dues[k] = dues[n];
    seqs[k] = seqs[n];
    timers[k] = timers[n];
    wraps[n] = null;
  };
  const settle = (w: Callback) => {
    let at = -1;
    for (let k = 0; k < n; k++)
      if (wraps[k] === w && (at < 0 || dues[k] < dues[at] || (dues[k] === dues[at] && seqs[k] < seqs[at]))) at = k;
    if (at >= 0) remove(at);
  };
  const track = (id: number, w: Callback, due: number, timer: boolean) => {
    ids[n] = id;
    wraps[n] = w;
    dues[n] = due;
    seqs[n] = seq++;
    timers[n] = timer;
    n++;
    return id;
  };
  const untrack = (id: number) => {
    for (let k = 0; k < n; k++) if (ids[k] === id) return remove(k);
  };
  // Made outside the request methods: a closure there over `cb` would allocate a context on every call.
  const wrapFrame = (cb: Callback) => {
    const f: Callback = (t) => {
      settle(f);
      if (!closed) cb(t);
    };
    frameWraps.set(cb, f);
    return f;
  };
  const wrapTimer = (cb: () => void) => {
    const f: Callback = () => {
      settle(f);
      if (!closed) cb();
    };
    timerWraps.set(cb, f);
    return f;
  };
  return {
    now: () => port.now(),
    requestFrame(cb) {
      if (closed) return -1;
      const w = frameWraps.get(cb) ?? wrapFrame(cb);
      return track(port.requestFrame(w), w, 0, false);
    },
    cancelFrame(id) {
      if (id < 0) return;
      untrack(id);
      port.cancelFrame(id);
    },
    setTimer(cb, ms) {
      if (closed) return -1;
      const w = timerWraps.get(cb) ?? wrapTimer(cb);
      return track(port.setTimer(w as () => void, ms), w, port.now() + Math.max(0, ms), true);
    },
    clearTimer(id) {
      if (id < 0) return;
      untrack(id);
      port.clearTimer(id);
    },
    get outstanding() {
      return n;
    },
    get closed() {
      return closed;
    },
    clearAll() {
      closed = true;
      for (let k = 0; k < n; k++) {
        if (timers[k]) port.clearTimer(ids[k]);
        else port.cancelFrame(ids[k]);
      }
      while (n > 0) remove(n - 1);
    },
  };
}
