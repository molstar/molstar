/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { Subscription } from 'rxjs';

// A viewer's canvas works only while it can be seen (styling.md §5: what is not shown does no work; README decision 29:
// several viewers on a page). An IntersectionObserver on the canvas' host decides (one for the page's viewers, so a
// frame's changes come in one callback): out of view (scrolled away, a hidden tab, a collapsed panel: anything that
// leaves it no area on screen), the plugin's animation loop stops without drawing; in view again, it starts and draws once
// what changed meanwhile.
//
// WebGL contexts are the page's: a browser keeps only so many live (Chrome 16) and drops the oldest when one more is
// made, and a context it dropped cannot be restored. So the page's viewers share a budget (`CanvasContexts.budget`,
// below the browser's limit) and the page never reaches the limit: past the budget, the contexts of viewers out of view
// are released (`WEBGL_lose_context.loseContext`), the longest unseen first; a viewer coming into view gets its context
// back (`restoreContext`; Mol*'s Canvas3D restores its resources on `webglcontextrestored`), up to `max` when more
// viewers than the budget are in view at once; past `max`, the viewers in view followed last give theirs up until a slot
// frees. The changes of a frame are reconciled together, releases before restores, and a context being restored counts
// as live until it is: a restore is asynchronous. A viewer just followed is reconciled at once, in the task that made its
// context: viewers mounted together must not pass the browser's limit before a microtask runs. A viewer in an XR session
// is left alone.

interface Waiter {
  resolve(): void;
  reject(e: unknown): void;
}

interface Presence {
  readonly plugin: PluginContext;
  readonly root: HTMLElement;
  readonly lose: WEBGL_lose_context | null;
  /** The order it was followed in: past `max`, those in view followed last give their contexts up first. */
  readonly order: number;
  visible: boolean;
  /** When it was last in view (ms), for the release order; 0: never. */
  seen: number;
  /** Released by the pool (not lost on its own): only those can be restored. */
  released: boolean;
  /** A restore asked for and not done yet: counted as live. */
  restoring: boolean;
  /** The browser has told it the context is lost (`webglcontextlost`): a restore before that is refused. */
  lostSeen: boolean;
  /** Held by work that renders on purpose (an export, a thumbnail): its loop and context are left alone. */
  holds: number;
  /** Waiting for its context to come back (`holdCanvas`). */
  waiters: Waiter[];
  subscriptions: Subscription[];
}

const viewers = new Set<Presence>();
const byPlugin = new WeakMap<PluginContext, Presence>();
const byRoot = new Map<Element, Presence>();
let followed = 0;

/**
 * The page's WebGL contexts: `budget` kept live (viewers out of view release theirs past it), `max` the most ever live,
 * below the browser's limit (Chrome's is 16; past it the browser drops a context that cannot come back), which viewers
 * in view may use when more than `budget` of them show at once.
 */
export const CanvasContexts = {
  budget: 12,
  max: 15,
  /** The viewers followed, those in view, the live contexts, and those in view without one (released or lost). */
  stats() {
    let inView = 0;
    let live = 0;
    let dark = 0;
    for (const p of viewers) {
      if (p.visible) inView++;
      if (isLive(p)) live++;
      if (p.visible && !isLive(p)) dark++;
    }
    return { viewers: viewers.size, inView, live, dark };
  },
};

function isLost(p: Presence): boolean {
  return !!p.plugin.canvas3d?.webgl.isContextLost;
}

/** Holds a context now or soon (a restore under way). */
function isLive(p: Presence): boolean {
  return p.restoring || !isLost(p);
}

function inXr(p: Presence): boolean {
  return !!p.plugin.canvas3d?.webgl.xr.session;
}

let scheduled = false;
function schedule(): void {
  if (scheduled) return;
  scheduled = true;
  queueMicrotask(reconcile);
}

function release(p: Presence): void {
  p.released = true;
  p.lose!.loseContext();
}

/** The frame's changes together: loops stopped out of view, releases to make room, then restores and loops started. */
function reconcile(): void {
  scheduled = false;
  for (const p of viewers) {
    if (p.visible || p.holds || inXr(p)) continue;
    const loop = p.plugin.animationLoop;
    if (loop.isAnimating) loop.stop({ noDraw: true });
  }
  const wanted = [...viewers].filter(
    (p) => (p.visible || p.holds > 0) && isLost(p) && p.released && p.lostSeen && !p.restoring && p.lose,
  );
  let live = 0;
  for (const p of viewers) if (isLive(p)) live++;
  const releasable = (p: Presence) => !p.holds && !p.restoring && !isLost(p) && !!p.lose && !inXr(p);
  const unseen = [...viewers].filter((p) => !p.visible && releasable(p)).sort((a, b) => a.seen - b.seen);
  for (const p of unseen) {
    if (live + wanted.length <= CanvasContexts.budget) break;
    release(p);
    live--;
  }
  // More in view than the page may hold (many viewers shown or mounted at once): those followed last go dark until a
  // slot frees, rather than the browser dropping one for good.
  if (live > CanvasContexts.max) {
    const shown = [...viewers].filter((p) => p.visible && releasable(p)).sort((a, b) => b.order - a.order);
    for (const p of shown) {
      if (live <= CanvasContexts.max) break;
      release(p);
      live--;
    }
  }
  for (const p of wanted) {
    // A held canvas gets its context whatever the count: the work holding it needs it now.
    if (live >= CanvasContexts.max && !p.holds) continue;
    p.released = false;
    p.restoring = true;
    live++;
    try {
      p.lose!.restoreContext();
    } catch (e) {
      // Refused: it stays released, to be asked again on the next change; work waiting for it is told.
      p.released = true;
      p.restoring = false;
      live--;
      for (const w of p.waiters.splice(0)) w.reject(e);
    }
  }
  for (const p of viewers) {
    if (!p.visible || p.holds || isLost(p)) continue;
    const loop = p.plugin.animationLoop;
    if (!loop.isAnimating) {
      loop.start();
      p.plugin.canvas3d?.requestDraw();
    }
  }
}

function inViewNow(root: HTMLElement): boolean {
  const rect = root.getBoundingClientRect();
  const view = root.ownerDocument.defaultView;
  if (!view) return true;
  return (
    rect.width > 0 &&
    rect.height > 0 &&
    rect.bottom > 0 &&
    rect.right > 0 &&
    rect.top < view.innerHeight &&
    rect.left < view.innerWidth
  );
}

/** Resolves once the context is back; rejects when the restore is refused or `signal` aborts. */
function untilLive(p: Presence, signal: AbortSignal | undefined): Promise<void> {
  return new Promise<void>((resolve, reject) => {
    const abort = () => {
      const i = p.waiters.indexOf(waiter);
      if (i >= 0) p.waiters.splice(i, 1);
      reject(signal!.reason);
    };
    const waiter: Waiter = {
      resolve() {
        signal?.removeEventListener('abort', abort);
        resolve();
      },
      reject(e) {
        signal?.removeEventListener('abort', abort);
        reject(e);
      },
    };
    signal?.addEventListener('abort', abort, { once: true });
    p.waiters.push(waiter);
    schedule();
  });
}

/**
 * What changed while the loop was stopped is made current before held work renders: Canvas3D commits added
 * representations and resolves camera resets only in its tick, which a stopped loop does not run, and an image pass
 * draws what is committed. A camera transition the stopped loop left under way is ended at its target.
 */
function settle(p: Presence): void {
  const canvas = p.plugin.canvas3d;
  if (!canvas || isLost(p)) return;
  canvas.commit(true);
  if (!p.plugin.animationLoop.isAnimating) canvas.camera.transition.complete();
}

/**
 * Work that renders on purpose (an image or video export, a snapshot's thumbnail, a screenshot) holds the canvas: its
 * loop and context are left alone until the returned release, and a released context (or one being restored) is back
 * first. What is pending in the scene is committed before this resolves, so a viewer out of view renders what it holds.
 * The release puts the loop back to what being in view says, whoever started it meanwhile (the screenshot helper starts
 * it after each image). A canvas not followed (headless, or pausing off) resolves at once; one whose viewer is disposed
 * while waiting resolves too. Rejects when `signal` aborts or the browser refuses the restore (nothing stays held).
 */
export async function holdCanvas(plugin: PluginContext, options: { signal?: AbortSignal } = {}): Promise<() => void> {
  options.signal?.throwIfAborted();
  const p = byPlugin.get(plugin);
  if (!p) return () => {};
  p.holds++;
  let done = false;
  const release = () => {
    if (done) return;
    done = true;
    p.holds--;
    schedule();
  };
  if (isLost(p) && (p.released || p.restoring)) {
    try {
      await untilLive(p, options.signal);
    } catch (e) {
      release();
      throw e;
    }
  }
  if (byPlugin.get(plugin) === p) settle(p);
  return release;
}

/**
 * Whether a page pauses canvases out of view (`options.ui.pauseOutOfView`): `auto` does, except in an automated browser
 * (`navigator.webdriver`: Playwright, Puppeteer, a capture service), where a page is rendered on purpose whatever is in
 * view; `always` and `never` force it.
 */
export function pausesOutOfView(mode: 'auto' | 'always' | 'never'): boolean {
  if (mode === 'never') return false;
  if (mode === 'always') return true;
  return !(typeof navigator !== 'undefined' && navigator.webdriver);
}

let observer: IntersectionObserver | undefined;
function observe(root: HTMLElement): void {
  if (typeof IntersectionObserver === 'undefined') return;
  observer ??= new IntersectionObserver((entries) => {
    // Entries in time order: a target's last one says where it is now.
    for (const e of entries) {
      const p = byRoot.get(e.target);
      if (!p || e.isIntersecting === p.visible) continue;
      p.visible = e.isIntersecting;
      if (p.visible) p.seen = performance.now();
    }
    schedule();
  });
  observer.observe(root);
}

/**
 * Starts following the canvas: pauses it out of view and puts its context in the page's budget, at once (call it in the
 * task that made the context). Following it again does nothing.
 */
export function followCanvasPresence(plugin: PluginContext): void {
  if (byPlugin.has(plugin)) return;
  const root = plugin.layout.root;
  const webgl = plugin.canvas3d?.webgl;
  if (!root || !webgl) return;
  // In view now? One read at mount: the observer answers a frame later, and a page making many viewers at once must not
  // pass the browser's limit meanwhile.
  const visible = inViewNow(root);
  const p: Presence = {
    plugin,
    root,
    lose: webgl.gl.getExtension('WEBGL_lose_context'),
    order: followed++,
    visible,
    seen: visible ? performance.now() : 0,
    released: false,
    restoring: false,
    lostSeen: false,
    holds: 0,
    waiters: [],
    subscriptions: [],
  };
  viewers.add(p);
  byPlugin.set(plugin, p);
  byRoot.set(root, p);
  const lost = plugin.canvas3dContext?.contextLost;
  if (lost) {
    p.subscriptions.push(
      lost.subscribe(() => {
        p.lostSeen = true;
        schedule();
      }),
    );
  }
  p.subscriptions.push(
    webgl.contextRestored.subscribe(() => {
      p.restoring = false;
      p.lostSeen = false;
      for (const w of p.waiters.splice(0)) w.resolve();
      plugin.canvas3d?.requestDraw();
      schedule();
    }),
  );
  observe(root);
  reconcile();
}

/** Stops following the canvas: work waiting for its context is let go, and the slot it frees goes to another viewer. */
export function unfollowCanvasPresence(plugin: PluginContext): void {
  const p = byPlugin.get(plugin);
  if (!p) return;
  observer?.unobserve(p.root);
  for (const s of p.subscriptions) s.unsubscribe();
  viewers.delete(p);
  byPlugin.delete(plugin);
  byRoot.delete(p.root);
  if (!viewers.size) {
    observer?.disconnect();
    observer = undefined;
  }
  for (const w of p.waiters.splice(0)) w.resolve();
  schedule();
}
