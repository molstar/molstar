/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { now } from '@molstar/core/util/now';
import type { PluginContext } from '@molstar/plugin/context';

// The canvas in a page's viewer. Resizes (README §3.5, step 2): the plugin waits out a 50 ms debounce and a 100 ms
// throttle, and sizing the canvas element clears it, so the frame after is blank until the next draw. Here the host's
// size is observed (a ResizeObserver: once per frame while it changes, after layout, before paint), the canvas is sized
// and drawn in the same callback, so the cleared canvas never shows; a size it already has is not set again (setting it
// clears too). The plugin's own trailing resize goes through the same path. Canvas3D is not touched (decision 28).

const observers = new WeakMap<PluginContext, ResizeObserver>();

/** Mounts the plugin's canvas in `host` (a page's viewer); false when the canvas cannot start (no WebGL). */
export async function mountCanvas(plugin: PluginContext, host: HTMLElement): Promise<boolean> {
  const ok = !!(await plugin.mountAsync(host));
  if (ok) smoothResizes(plugin);
  return ok;
}

export function unmountCanvas(plugin: PluginContext): void {
  observers.get(plugin)?.disconnect();
  observers.delete(plugin);
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
