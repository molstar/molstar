/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RenderSvgOptions, RenderSvgResult } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom/host';

// The stories' SVG beside a canvas: the viewer's export of the spec it shows, drawn again when the view comes to rest,
// the theme changes, the spec is edited (its revision: an edit or an undo, the frame after it, not a rest later, so the
// twin is never seen a step behind the canvas) or a story step is shown (a preview or an animation frame waits for the
// rest), a clade or a track group is folded or opened (by the reader, or restored by Back), the selection is set or
// cleared, a highlight channel (a search's hits) is set or cleared, or any state the story names moves. A later draw
// drops an earlier one still in flight; after `dispose` none lands. A selection set during a drag is drawn after the
// first frame after its release, never per move: an export is a whole picture (half a second on the 4,000 × 4,000 dot
// plot), and a drag's work is bounded by the screen (the viewer's _spec, "SelectionModel and HighlightModel", Drags).
// SVG ids are prefixed per story instance, so SVGs in one page never share a clip or a mask.
//
// A change's export waits for the viewer's frame that draws the change, then runs in a task of its own: never in the
// input's task, nor in the microtasks after its frame. Why: the export is a second resolve, layout and draw of the whole
// view (and an `innerHTML` of its picture), 6 to 25 ms on the 200 copies, the ring and a search's Next (2026-11-24);
// run at once it doubled the edit's time to drawn. Changes before the export starts are one export.

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
    request();
  };
  // a change waits for the frame that draws it (asked for), then the export runs in a task of its own, once however
  // many changes came before it starts
  let due = false;
  let task: unknown;
  const run = () => {
    task = undefined;
    if (alive) void draw();
  };
  const request = () => {
    if (!alive) return;
    due = true;
    viewer.schedule();
  };
  // the rest's frame is drawn already: its export's task at once
  const later = () => {
    if (alive && task === undefined) task = DomHost.setTimeout(run, 0);
  };
  // a selection set while a drag runs: drawn once, after the first frame after the release
  let afterDrag = false;
  const onSelection = () => {
    if (viewer.input.dragging) afterDrag = true;
    else request();
  };
  const onFrame = () => {
    if (afterDrag && !viewer.input.dragging) {
      afterDrag = false;
      due = true;
    }
    if (!due || task !== undefined) return;
    due = false;
    task = DomHost.setTimeout(run, 0);
  };
  const offs = [
    viewer.tokensVersion,
    viewer.spec.revision,
    viewer.spec.step,
    viewer.highlight.version,
    viewer.view.groupFolds,
    ...(options.also ?? []),
  ].map((state) => state.listen(request));
  offs.push(
    viewer.rested.listen(later),
    viewer.view.bands.listen(onBands),
    viewer.selection.value.listen(onSelection),
    viewer.drawn.listen(onFrame),
  );
  return {
    draw,
    dispose() {
      alive = false;
      if (task !== undefined) DomHost.clearTimeout(task);
      for (const off of offs) off();
    },
  };
}
