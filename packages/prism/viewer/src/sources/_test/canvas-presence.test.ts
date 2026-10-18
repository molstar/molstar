/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// @vitest-environment jsdom

import { Subject } from 'rxjs';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import type { PluginContext } from '@molstar/plugin/context';
import { mountCanvas, unmountCanvas } from '../canvas-host.js';
import {
  CanvasContexts,
  followCanvasPresence,
  holdCanvas,
  pausesOutOfView,
  unfollowCanvasPresence,
} from '../canvas-presence.js';

// The canvas' presence (canvas-presence.ts) with fake plugins: a fake canvas context whose release and restore report
// a task later (as browsers do; a restore can be told to be refused; with `asyncLoss`, as WebGPU, the context reads as
// lost only then too), a canvas that records commits, one IntersectionObserver for the page driven by hand.

type Entry = { target: Element; isIntersecting: boolean };
let observerCallback: ((entries: Entry[]) => void) | undefined;
let observersMade = 0;
const observed = new Set<Element>();
class FakeIntersectionObserver {
  constructor(cb: (entries: Entry[]) => void) {
    observerCallback = cb;
    observersMade++;
  }
  observe(el: Element) {
    observed.add(el);
  }
  unobserve(el: Element) {
    observed.delete(el);
  }
  disconnect() {
    observed.clear();
  }
}

const tick = () => new Promise<void>((r) => setTimeout(r, 0));

/**
 * The viewer's plugin makes its context in `mountAsync`, as the plugin does: at once (WebGL's `getContext`), but the
 * plugin has it only a task later (the backend chosen asynchronously). `contexts.peak`: the most live at once.
 */
function mountLater(
  v: ReturnType<typeof fakeViewer>,
  all: ReturnType<typeof fakeViewer>[],
  contexts: { peak: number },
) {
  const plugin = v.plugin as unknown as Record<string, unknown>;
  const { canvas3dContext, canvas3d } = plugin;
  plugin.canvas3dContext = undefined;
  plugin.canvas3d = undefined;
  const state = { created: false, mounts: 0 };
  plugin.mountAsync = async () => {
    state.mounts++;
    state.created = true;
    const live = all.filter((o) => o.created() && !o.isLost()).length;
    contexts.peak = Math.max(contexts.peak, live);
    await tick();
    plugin.canvas3dContext = canvas3dContext;
    plugin.canvas3d = canvas3d;
    return true;
  };
  plugin.unmount = () => {};
  v.created = () => state.created;
  return state;
}

function fakeViewer(inView: boolean, asyncLoss = false) {
  const root = document.createElement('div');
  root.getBoundingClientRect = () =>
    (inView
      ? { width: 100, height: 100, top: 0, left: 0, bottom: 100, right: 100 }
      : { width: 100, height: 100, top: 5000, left: 0, bottom: 5100, right: 100 }) as DOMRect;
  let lost = false;
  let released = false;
  let lostEventSent = false;
  const contextLost = new Subject<number>();
  const contextRestored = new Subject<number>();
  const control = { refuseRestore: false, holdRestore: false, restores: 0 };
  let restoring: ((restored: boolean) => void) | undefined;
  const finishRestore = () => {
    lost = false;
    lostEventSent = false;
    contextRestored.next(0);
    restoring?.(true);
    restoring = undefined;
  };
  const context = {
    graphics: {
      get isContextLost() {
        return lost;
      },
      xr: { isPresenting: () => false },
    },
    contextLost,
    contextRestored,
    get isReleased() {
      return released;
    },
    release() {
      if (lost || released) return false;
      released = true;
      if (!asyncLoss) lost = true;
      setTimeout(() => {
        lost = true;
        lostEventSent = true;
        contextLost.next(0);
      });
      return true;
    },
    restore() {
      if (!lostEventSent || control.refuseRestore) return Promise.reject(new Error('restore not allowed'));
      released = false;
      control.restores++;
      return new Promise<boolean>((resolve) => {
        restoring = resolve;
        if (!control.holdRestore) setTimeout(finishRestore);
      });
    },
  };
  const loop = { isAnimating: true, starts: 0, stops: 0 };
  const scene = { pending: 0, commits: 0, inTransition: false };
  const plugin = {
    layout: { root },
    canvas3dContext: context,
    canvas3d: {
      camera: {
        transition: {
          get inTransition() {
            return scene.inTransition;
          },
          complete() {
            scene.inTransition = false;
          },
        },
      },
      commit() {
        scene.commits++;
        scene.pending = 0;
      },
      requestDraw() {},
    },
    animationLoop: {
      get isAnimating() {
        return loop.isAnimating;
      },
      start() {
        loop.isAnimating = true;
        loop.starts++;
      },
      stop() {
        loop.isAnimating = false;
        loop.stops++;
      },
    },
  } as unknown as PluginContext;
  return {
    plugin,
    loop,
    root,
    scene,
    control,
    finishRestore,
    isLost: () => lost,
    /** Whether its context was made (see `mountLater`). */
    created: () => true,
    show: (on: boolean) => observed.has(root) && observerCallback?.([{ target: root, isIntersecting: on }]),
  };
}

describe('canvas presence', () => {
  const made: ReturnType<typeof fakeViewer>[] = [];
  const make = (inView: boolean, asyncLoss = false) => {
    const v = fakeViewer(inView, asyncLoss);
    followCanvasPresence(v.plugin);
    made.push(v);
    return v;
  };
  beforeEach(() => {
    (globalThis as any).IntersectionObserver = FakeIntersectionObserver;
    observersMade = 0;
    CanvasContexts.budget = 3;
    CanvasContexts.max = 4;
  });
  afterEach(() => {
    for (const v of made.splice(0)) unfollowCanvasPresence(v.plugin);
    CanvasContexts.budget = 12;
    CanvasContexts.max = 15;
  });

  it('stops the loop out of view and starts it in view', async () => {
    const v = make(true);
    v.show(false);
    await tick();
    expect(v.loop.isAnimating).toBe(false);
    v.show(true);
    await tick();
    expect(v.loop.isAnimating).toBe(true);
  });

  it('keeps live contexts within the budget, releasing those out of view first', async () => {
    const shown = [make(true), make(true)];
    const hidden = [make(false), make(false), make(false)];
    await tick();
    await tick();
    expect(CanvasContexts.stats().live).toBe(3);
    for (const v of shown) expect(v.isLost()).toBe(false);
    expect(hidden.filter((v) => v.isLost()).length).toBe(2);
  });

  it('restores a context when its viewer comes into view, without passing the budget', async () => {
    const a = make(true);
    const b = make(false);
    const c = make(false);
    const d = make(false);
    await tick();
    await tick();
    const lostOnes = [b, c, d].filter((v) => v.isLost());
    expect(lostOnes.length).toBeGreaterThan(0);
    const back = lostOnes[0];
    a.show(false);
    back.show(true);
    await tick();
    await tick();
    await tick();
    expect(back.isLost()).toBe(false);
    expect(CanvasContexts.stats().live).toBeLessThanOrEqual(CanvasContexts.budget);
    expect(CanvasContexts.stats().dark).toBe(0);
  });

  it('never restores before the browser said the context was lost', async () => {
    const a = make(true);
    const b = make(true);
    const c = make(true);
    const d = make(false);
    await tick();
    await tick();
    // d was released at mount; show it at once, then hide another: the restore waits for d's lost event.
    expect(d.isLost()).toBe(true);
    c.show(false);
    d.show(true);
    await tick();
    await tick();
    await tick();
    expect(d.isLost()).toBe(false);
    expect([a, b].every((v) => !v.isLost())).toBe(true);
  });

  it('counts a context released before its loss is reported (WebGPU): not released twice, restored in view', async () => {
    const [, , c] = [make(true, true), make(true, true), make(true, true)];
    const d = make(false, true);
    // d released at mount, its loss not reported yet; a viewer followed now reconciles in this task.
    const e = make(false, true);
    expect(CanvasContexts.stats().live).toBeLessThanOrEqual(CanvasContexts.budget);
    await tick();
    await tick();
    expect(d.isLost()).toBe(true);
    expect(e.isLost()).toBe(true);
    c.show(false);
    d.show(true);
    for (let i = 0; i < 4; i++) await tick();
    expect(d.isLost()).toBe(false);
    expect(CanvasContexts.stats().live).toBeLessThanOrEqual(CanvasContexts.budget);
  });

  it('lets viewers in view use up to max when more than the budget show', async () => {
    for (let i = 0; i < 5; i++) make(false);
    await tick();
    await tick();
    for (const v of made) v.show(true);
    for (let i = 0; i < 6; i++) await tick();
    const s = CanvasContexts.stats();
    expect(s.live).toBe(CanvasContexts.max);
    expect(s.dark).toBe(1);
  });

  it('leaves a held canvas alone: no release, no loop changes, and restores it first if it was released', async () => {
    const a = make(true);
    const b = make(true);
    const c = make(true);
    const d = make(false);
    await tick();
    await tick();
    expect(d.isLost()).toBe(true);
    // An export of d (scrolled away): it gets its context back before the work starts.
    const release = await holdCanvas(d.plugin);
    expect(d.isLost()).toBe(false);
    const stopsBefore = d.loop.stops;
    const startsBefore = d.loop.starts;
    // While held, scrolling the others in and out neither releases d nor touches its loop.
    a.show(false);
    b.show(false);
    c.show(false);
    for (let i = 0; i < 4; i++) await tick();
    expect(d.isLost()).toBe(false);
    expect(d.loop.stops).toBe(stopsBefore);
    expect(d.loop.starts).toBe(startsBefore);
    release();
    for (let i = 0; i < 3; i++) await tick();
    // Released and out of view: its loop stops again.
    expect(d.loop.isAnimating).toBe(false);
  });

  it("follows every canvas with one observer, so a frame's changes come in one callback", async () => {
    const a = make(true);
    const b = make(true);
    make(false);
    expect(observersMade).toBe(1);
    observerCallback?.([
      { target: a.root, isIntersecting: false },
      { target: b.root, isIntersecting: false },
    ]);
    await tick();
    expect(a.loop.isAnimating).toBe(false);
    expect(b.loop.isAnimating).toBe(false);
    expect(CanvasContexts.stats().inView).toBe(0);
  });

  it('keeps the budget when many viewers are followed in one task, before any microtask', () => {
    for (let i = 0; i < 10; i++) {
      make(false);
      expect(CanvasContexts.stats().live).toBeLessThanOrEqual(CanvasContexts.budget);
    }
    // In view past max: those followed last go dark instead of the browser dropping one.
    for (let i = 0; i < 6; i++) {
      make(true);
      expect(CanvasContexts.stats().live).toBeLessThanOrEqual(CanvasContexts.max);
    }
    expect(made.slice(10, 14).every((v) => !v.isLost())).toBe(true);
    expect(made.slice(14).every((v) => v.isLost())).toBe(true);
  });

  it('gives a dark viewer in view the slot another frees', async () => {
    for (let i = 0; i < 5; i++) make(true);
    await tick();
    const dark = made[4];
    expect(dark.isLost()).toBe(true);
    made[0].show(false);
    for (let i = 0; i < 4; i++) await tick();
    expect(dark.isLost()).toBe(false);
    expect(CanvasContexts.stats().dark).toBe(0);
  });

  it('commits what is pending before held work renders a viewer out of view, and ends its transition', async () => {
    const v = make(false);
    await tick();
    expect(v.loop.isAnimating).toBe(false);
    // A structure loaded while it was scrolled away: added to the scene, not committed (no tick runs).
    v.scene.pending = 1;
    v.scene.inTransition = true;
    const release = await holdCanvas(v.plugin);
    expect(v.scene.commits).toBe(1);
    expect(v.scene.pending).toBe(0);
    expect(v.scene.inTransition).toBe(false);
    release();
  });

  it('puts the loop back out of view after a screenshot started it', async () => {
    const v = make(false);
    await tick();
    const release = await holdCanvas(v.plugin);
    // ViewportScreenshotHelper.draw ends with animationLoop.start({ immediate: true }).
    v.plugin.animationLoop.start();
    release();
    await tick();
    expect(v.loop.isAnimating).toBe(false);
  });

  it('waits for a restore already under way', async () => {
    const a = make(true);
    const b = make(true);
    const c = make(true);
    const d = make(false);
    await tick();
    await tick();
    expect(d.isLost()).toBe(true);
    // d comes into view: its restore is asked for and not done yet.
    d.control.holdRestore = true;
    c.show(false);
    d.show(true);
    await tick();
    expect(d.control.restores).toBe(1);
    expect(d.isLost()).toBe(true);
    let held = false;
    const hold = holdCanvas(d.plugin).then((release) => {
      held = true;
      return release;
    });
    await tick();
    expect(held).toBe(false);
    d.finishRestore();
    const release = await hold;
    expect(d.isLost()).toBe(false);
    release();
    expect([a, b].every((v) => !v.isLost())).toBe(true);
  });

  it('lets a waiting hold go when its viewer is disposed, and gives the slot to a dark viewer in view', async () => {
    for (let i = 0; i < 3; i++) make(true);
    const d = make(false);
    await tick();
    await tick();
    d.control.holdRestore = true;
    const hold = holdCanvas(d.plugin);
    await tick();
    // d's restore counts as live: a viewer shown now is past max and goes dark.
    const dark = make(true);
    expect(dark.isLost()).toBe(true);
    unfollowCanvasPresence(d.plugin);
    made.splice(made.indexOf(d), 1);
    const release = await hold;
    release();
    for (let i = 0; i < 4; i++) await tick();
    expect(dark.isLost()).toBe(false);
    expect(CanvasContexts.stats().live).toBeLessThanOrEqual(CanvasContexts.max);
  });

  it('stops waiting when the hold is aborted, and fails when the restore is refused', async () => {
    for (let i = 0; i < 3; i++) make(true);
    const d = make(false);
    await tick();
    await tick();
    expect(d.isLost()).toBe(true);
    d.control.holdRestore = true;
    const stop = new AbortController();
    const hold = holdCanvas(d.plugin, { signal: stop.signal });
    await tick();
    stop.abort(new Error('Stopped'));
    await expect(hold).rejects.toThrow('Stopped');
    d.finishRestore();
    await tick();
    // Nothing is held after the abort: out of view, its loop stays stopped.
    expect(d.loop.isAnimating).toBe(false);
    for (let i = 0; i < 3; i++) await tick();
    expect(d.isLost()).toBe(true);
    d.control.holdRestore = false;
    d.control.refuseRestore = true;
    await expect(holdCanvas(d.plugin)).rejects.toThrow('restore not allowed');
    // An aborted signal fails at once.
    await expect(holdCanvas(d.plugin, { signal: AbortSignal.abort(new Error('Already')) })).rejects.toThrow('Already');
  });

  it('pauses out of view unless automated or told otherwise', () => {
    expect(pausesOutOfView('always')).toBe(true);
    expect(pausesOutOfView('never')).toBe(false);
    const nav = globalThis.navigator as { webdriver?: boolean };
    const was = Object.getOwnPropertyDescriptor(nav, 'webdriver');
    Object.defineProperty(nav, 'webdriver', { value: true, configurable: true });
    expect(pausesOutOfView('auto')).toBe(false);
    Object.defineProperty(nav, 'webdriver', { value: false, configurable: true });
    expect(pausesOutOfView('auto')).toBe(true);
    if (was) Object.defineProperty(nav, 'webdriver', was);
    else delete (nav as { webdriver?: boolean }).webdriver;
  });

  it('mounts viewers together within the limit while their contexts are made asynchronously', async () => {
    const contexts = { peak: 0 };
    const all: ReturnType<typeof fakeViewer>[] = [];
    for (let i = 0; i < 8; i++) {
      const v = fakeViewer(true);
      made.push(v);
      all.push(v);
      mountLater(v, all, contexts);
    }
    const mounted = all.map((v) => mountCanvas(v.plugin, v.root, { pauseOutOfView: 'always' }));
    expect(await Promise.all(mounted)).toEqual(all.map(() => true));
    for (let i = 0; i < 3; i++) await tick();
    // One past max at most (below the browser's limit); then those followed last go dark.
    expect(contexts.peak).toBeLessThanOrEqual(CanvasContexts.max + 1);
    expect(CanvasContexts.stats().live).toBe(CanvasContexts.max);
    expect(CanvasContexts.stats().dark).toBe(all.length - CanvasContexts.max);
  });

  it('lets a mount waiting for a slot go when its viewer is unmounted meanwhile', async () => {
    CanvasContexts.budget = 1;
    CanvasContexts.max = 1;
    const contexts = { peak: 0 };
    const all: ReturnType<typeof fakeViewer>[] = [];
    const states = [0, 1, 2].map(() => {
      const v = fakeViewer(true);
      made.push(v);
      all.push(v);
      return mountLater(v, all, contexts);
    });
    const mounted = all.map((v) => mountCanvas(v.plugin, v.root, { pauseOutOfView: 'always' }));
    unmountCanvas(all[2].plugin);
    expect(await Promise.all(mounted)).toEqual([true, true, false]);
    expect(states[2].mounts).toBe(0);
    expect(contexts.peak).toBeLessThanOrEqual(2);
  });
});
