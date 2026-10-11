/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { onCleanup } from '../../dom/scope.js';
import type { ViewRenderer } from '../../dom/view.js';

/** Mounts a contributed view (any technology) into a container owned by the current scope. */
export function ViewHost<Ctx>(props: { renderer: ViewRenderer<Ctx>; ctx: Ctx; class?: string }) {
  const host = (<div class={props.class} style="display: contents" />) as HTMLElement;
  const dispose = props.renderer.mount(host, props.ctx);
  onCleanup(dispose);
  return host;
}
