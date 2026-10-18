/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec2 } from '@molstar/core/math/linear-algebra';
import { now } from '@molstar/core/util/now';
import type { PluginContext } from '@molstar/plugin/context';
import { admitCanvas, followCanvasPresence, pausesOutOfView, unfollowCanvasPresence } from './canvas-presence.js';

// The canvas in a page's viewer. Resizes (README decision 28): the plugin waits out a 50 ms debounce and a 100 ms
// throttle, and sizing the canvas element clears it, so the frame after is blank until the next draw. Here the host's
// size is observed (a ResizeObserver: once per frame while it changes, after layout, before paint), the canvas is sized
// and drawn in the same callback, so the cleared canvas never shows; a size it already has is not set again (setting it
// clears too). The plugin's own trailing resize goes through the same path. Canvas3D is not touched (decision 28).

const observers = new WeakMap<PluginContext, ResizeObserver>();
/** The mount under way or done, per plugin: an unmount while a mount waits for a slot cancels it. */
const mounts = new WeakMap<PluginContext, object>();

/** Mounts the plugin's canvas in `host` (a page's viewer); false when the canvas cannot start (no WebGL). */
export async function mountCanvas(
  plugin: PluginContext,
  host: HTMLElement,
  o: { pauseOutOfView?: 'auto' | 'always' | 'never' } = {},
): Promise<boolean> {
  const pauses = pausesOutOfView(o.pauseOutOfView ?? 'auto');
  const token = {};
  mounts.set(plugin, token);
  // The context is made asynchronously (the backend is chosen first): the mount waits for a slot in the page's budget and
  // its context counts as live until it is followed, so viewers mounted together never pass the browser's limit (README
  // decision 29).
  const admitted = pauses ? await admitCanvas() : undefined;
  try {
    // Unmounted while waiting for a slot.
    if (mounts.get(plugin) !== token) return false;
    const ok = !!(await plugin.mountAsync(host));
    if (ok) {
      smoothResizes(plugin);
      longPress(plugin);
      if (pauses) followCanvasPresence(plugin);
    } else unfollowCanvasPresence(plugin);
    return ok;
  } finally {
    admitted?.();
  }
}

export function unmountCanvas(plugin: PluginContext): void {
  mounts.delete(plugin);
  unfollowCanvasPresence(plugin);
  observers.get(plugin)?.disconnect();
  observers.delete(plugin);
  presses.get(plugin)?.();
  presses.delete(plugin);
  plugin.unmount();
}

function smoothResizes(plugin: PluginContext): void {
  const root = plugin.layout.root;
  if (!root || typeof ResizeObserver === 'undefined') return;
  const resize = plugin.handleResize;
  let last = '';
  const handle = () => {
    const context = plugin.canvas3dContext;
    if (!context?.canvas) return;
    const size = `${root.offsetWidth}x${root.offsetHeight}@${window.devicePixelRatio}x${context.pixelScale}`;
    if (size === last) return;
    last = size;
    resize();
    // Drawn now, before this frame paints.
    plugin.canvas3d?.tick(now());
  };
  // The plugin calls it as a property (its debounced resize, a layout change): ours stands in.
  (plugin as { handleResize: () => void }).handleResize = handle;
  const observer = new ResizeObserver(handle);
  observer.observe(root);
  observers.set(plugin, observer);
}

const presses = new WeakMap<PluginContext, () => void>();
/** How long a touch is held to be a long press (ms), and how far it may move (px). */
export const PressMs = 450;
export const PressSlop = 8;

/**
 * Touch (accessibility.md): what hover gives a mouse, a long press gives a finger: the label of what is under it (the
 * plugin's highlight, so the hover label shows it), until the next touch.
 */
export function longPress(plugin: PluginContext): void {
  // On the canvas' parent, for the canvas of the moment: a backend fallback replaces the canvas element in place.
  const host = plugin.canvas3dContext?.canvas?.parentElement;
  if (!host) return;
  const onCanvas = (e: Event) => e.target === plugin.canvas3dContext?.canvas;
  let timer: ReturnType<typeof setTimeout> | undefined;
  let x = 0;
  let y = 0;
  const cancel = () => {
    if (timer !== undefined) clearTimeout(timer);
    timer = undefined;
  };
  const show = async () => {
    timer = undefined;
    const canvas = plugin.canvas3dContext?.canvas;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    // read asynchronously on every backend (WebGPU has no synchronous read)
    const picked = await plugin.canvas3d?.asyncIdentify(Vec2.create(x - rect.left, y - rect.top))?.wait();
    const loci = picked && plugin.canvas3d?.getLoci(picked.id).loci;
    if (loci) plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
  };
  const down = (e: PointerEvent) => {
    if (e.pointerType !== 'touch' || !onCanvas(e)) return;
    cancel();
    x = e.clientX;
    y = e.clientY;
    timer = setTimeout(show, PressMs);
  };
  const move = (e: PointerEvent) => {
    if (timer !== undefined && Math.hypot(e.clientX - x, e.clientY - y) > PressSlop) cancel();
  };
  host.addEventListener('pointerdown', down);
  host.addEventListener('pointermove', move);
  host.addEventListener('pointerup', cancel);
  host.addEventListener('pointercancel', cancel);
  presses.set(plugin, () => {
    cancel();
    host.removeEventListener('pointerdown', down);
    host.removeEventListener('pointermove', move);
    host.removeEventListener('pointerup', cancel);
    host.removeEventListener('pointercancel', cancel);
  });
}
