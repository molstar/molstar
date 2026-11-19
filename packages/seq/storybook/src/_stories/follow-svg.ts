/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RenderSvgOptions, RenderSvgResult } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The stories' SVG beside a canvas: the viewer's export of the spec it shows, drawn again when the view comes to rest,
// the theme changes, the spec is edited (its revision: an edit or an undo, at once, not a rest later, so the twin is never seen a step
// behind the canvas) or a story step is shown (a preview or an animation frame waits for the rest), a clade is folded or opened, or any state the story names moves. A later draw drops an earlier one still in flight; after
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

/** Each band's folds, in the bands' order. */
function foldsOf(viewer: SeqViewer): readonly (readonly string[] | undefined)[] {
  return [...viewer.view.bands.value.values()].map((b) => b.folds);
}

/** Keeps `host` showing the viewer's SVG; `also` are further states whose change draws it again. */
export function followSvg(
  viewer: SeqViewer,
  host: HTMLElement,
  options: {
    readonly prefix?: string;
    readonly also?: readonly Listenable[];
    /** The export's own options (`navigators`). */
    readonly svg?: Partial<RenderSvgOptions>;
    /** Called with each export that lands (a dropped one never is). */
    readonly onDraw?: (r: RenderSvgResult) => void;
  } = {},
): FollowedSvg {
  const idPrefix = `${options.prefix ?? 'sv'}${++instances}-`;
  let alive = true;
  let latest = 0;
  const draw = async () => {
    const mine = ++latest;
    const r = await viewer.renderSvg({ idPrefix, ...options.svg });
    if (!alive || mine !== latest) return;
    // our own output: every string from data or the spec is escaped by the backend
    host.innerHTML = r.svg;
    options.onDraw?.(r);
  };
  // the folds are session state: only a change of them draws (a scroll moves `bands` too; the rest follows it)
  let folds = foldsOf(viewer);
  const onBands = () => {
    const next = foldsOf(viewer);
    if (next.length === folds.length && next.every((f, i) => f === folds[i])) return;
    folds = next;
    void draw();
  };
  const offs = [
    viewer.rested,
    viewer.tokensVersion,
    viewer.spec.revision,
    viewer.spec.step,
    ...(options.also ?? []),
  ].map((state) => state.listen(() => void draw()));
  offs.push(viewer.view.bands.listen(onBands));
  return {
    draw,
    dispose() {
      alive = false;
      for (const off of offs) off();
    },
  };
}
