/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The stories' SVG beside a canvas: the viewer's export of the spec it shows, drawn again when the view comes to rest,
// the theme changes, or any state the story names moves. A later draw drops an earlier one still in flight; after
// `dispose` none lands. SVG ids are prefixed per story instance, so SVGs in one page never share a clip or a mask.

let instances = 0;

export interface FollowedSvg {
  /** Draws now (and returns when it has landed or been dropped). */
  draw(): Promise<void>;
  /** Stops following: the listeners go, an export in flight is dropped. */
  dispose(): void;
}

export interface Listenable {
  listen(fn: () => void): () => void;
}

/** Keeps `host` showing the viewer's SVG; `also` are further states whose change draws it again. */
export function followSvg(
  viewer: SeqViewer,
  host: HTMLElement,
  options: { readonly prefix?: string; readonly also?: readonly Listenable[] } = {},
): FollowedSvg {
  const idPrefix = `${options.prefix ?? 'sv'}${++instances}-`;
  let alive = true;
  let latest = 0;
  const draw = async () => {
    const mine = ++latest;
    const r = await viewer.renderSvg({ idPrefix });
    if (!alive || mine !== latest) return;
    // our own output: every string from data or the spec is escaped by the backend
    host.innerHTML = r.svg;
  };
  const offs = [viewer.rested, viewer.tokensVersion, ...(options.also ?? [])].map((state) =>
    state.listen(() => void draw()),
  );
  return {
    draw,
    dispose() {
      alive = false;
      for (const off of offs) off();
    },
  };
}
