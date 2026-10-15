/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';

/** Mounts the plugin's canvas in `host` (a page's viewer); false when the canvas cannot start (no WebGL). */
export async function mountCanvas(plugin: PluginContext, host: HTMLElement): Promise<boolean> {
  return !!(await plugin.mountAsync(host));
}

export function unmountCanvas(plugin: PluginContext): void {
  plugin.unmount();
}
