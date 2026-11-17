/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { bind, frames, observeResize, on, requireScope, runInScope } from '../dom/index.js';
import { classes } from './_generated/hover-label.module.css.js';
import type { HoverAt, HoverLabelModel, HoverLine, HoverRun } from './stage-model.js';

// The hover label (stage §3.1): what is under the pointer, in the stage's corner or beside the pointer. Its content is the
// model's lines (text only); its place beside the pointer is the label's size read and one write a frame, allocating
// nothing (`_test/hover-label.quick.ts`); the layer's box is read once and again only after a resize, a scroll or the
// pointer coming in.

/** The label's distance from the pointer (px). */
const Offset = 12;

/** A line as an element: plain text (the first is the title), or formatted runs. Text only (`textContent`). */
function lineElement(line: HoverLine, index: number): HTMLElement {
  const el = document.createElement('span');
  if (typeof line === 'string') {
    el.className = index === 0 ? classes.primary : classes.secondary;
    el.textContent = line;
    return el;
  }
  el.className =
    line.role === 'title'
      ? classes.primary
      : line.role === 'item'
        ? `${classes.secondary} ${classes.item}`
        : classes.secondary;
  for (const run of line.runs) el.append(runNode(run));
  return el;
}

function runNode(run: HoverRun): Node {
  let node: Node = document.createTextNode(run.text);
  const wrap = (tag: string, className?: string) => {
    const el = document.createElement(tag);
    if (className) el.className = className;
    el.append(node);
    node = el;
  };
  if (run.code) wrap('code', classes.code);
  if (run.sub) wrap('sub');
  if (run.sup) wrap('sup');
  if (run.small) wrap('small');
  if (run.italic) wrap('em');
  if (run.bold) wrap('strong');
  if (run.link) wrap('span', classes.link);
  return node;
}

/** What the follower writes to: the label (an interface, so its per-frame path is measured with stand-ins). */
export interface FollowedLabel {
  readonly offsetWidth: number;
  readonly offsetHeight: number;
  readonly classList: { toggle(name: string, force?: boolean): boolean };
  readonly style: { setProperty(name: string, value: string): void };
  setAttribute(name: string, value: string): void;
  removeAttribute(name: string): void;
}

/** What the label is placed in (its offset parent): read once, and again after `HoverFollower.invalidate`. */
export interface FollowedLayer {
  getBoundingClientRect(): { readonly left: number; readonly top: number };
  readonly clientWidth: number;
  readonly clientHeight: number;
}

/** `${n}px` for whole pixels, made once each: a frame's write allocates no string. */
const PxStrings: string[] = [];
const PxCached = 4096;
function px(n: number): string {
  if (n < 0 || n >= PxCached) return `${n}px`;
  let s = PxStrings[n];
  if (s === undefined) PxStrings[n] = s = `${n}px`;
  return s;
}

/**
 * Whether and where the label shows (stage §3.1): at once when the pointer comes with content (in the corner; beside the
 * pointer in the frame that places it); a drag (a button held) or no content hides it at once; the pointer leaving
 * waits a frame (out and back, or across a small element, keeps it). Beside the pointer it is placed once a frame: the
 * label's size read, then the two offsets written (custom properties of cached strings), flipped near the edges; the
 * layer and its box are cached until `invalidate` (a resize, a scroll: no rect allocated a frame). Made in a scope: its
 * frames are the scope's.
 */
export class HoverFollower {
  #x = 0;
  #y = 0;
  #inside = false;
  #dragging = false;
  #content = false;
  #shown = false;
  #left = -1;
  #top = -1;
  /** The layer's box, cached (`#cached`) until `invalidate`. */
  #cached = false;
  #layerLeft = 0;
  #layerTop = 0;
  #layerWidth = 0;
  #layerHeight = 0;
  readonly #label: FollowedLabel;
  readonly #layer: () => FollowedLayer | null;
  readonly #corner: boolean;
  readonly #out: string;
  readonly #placing: { request(): void; cancel(): void };
  readonly #leaving: { request(): void; cancel(): void };

  constructor(label: FollowedLabel, layer: () => FollowedLayer | null, corner: boolean, outClass: string) {
    this.#label = label;
    this.#layer = layer;
    this.#corner = corner;
    this.#out = outClass;
    this.#placing = frames(() => this.#place());
    this.#leaving = frames(() => {
      if (!this.#visible()) this.#show(false);
    });
  }

  get shown(): boolean {
    return this.#shown;
  }

  /** The pointer at (x, y) (client px); `inside`: over what the label is for; `dragging`: a button held. */
  move(x: number, y: number, dragging: boolean, inside: boolean): void {
    this.#x = x;
    this.#y = y;
    this.#dragging = dragging;
    this.#inside = inside;
    this.#update();
  }

  leave(): void {
    this.#inside = false;
    this.#update();
  }

  /** The layer moved or resized (a resize, a scroll, the pointer coming in): its box is read again in the next frame. */
  invalidate(): void {
    this.#cached = false;
  }

  setContent(has: boolean): void {
    this.#content = has;
    this.#update();
  }

  #visible(): boolean {
    return this.#inside && !this.#dragging && this.#content;
  }

  #update(): void {
    if (this.#visible()) {
      this.#leaving.cancel();
      if (this.#corner) this.#show(true);
      else this.#placing.request();
    } else if (this.#dragging || !this.#content) {
      this.#leaving.cancel();
      this.#show(false);
    } else if (this.#shown) this.#leaving.request();
  }

  #show(on: boolean): void {
    if (this.#shown === on) return;
    this.#shown = on;
    this.#label.classList.toggle(this.#out, !on);
    // Out of the accessibility tree while hidden, whatever the CSS (not loaded, or still fading).
    if (on) this.#label.removeAttribute('aria-hidden');
    else this.#label.setAttribute('aria-hidden', 'true');
  }

  /** The layer's box, read when not cached; false while there is no layer (not laid out). */
  #measure(): boolean {
    if (this.#cached) return true;
    const layer = this.#layer();
    if (!layer) return false;
    const rect = layer.getBoundingClientRect();
    this.#layerLeft = rect.left;
    this.#layerTop = rect.top;
    this.#layerWidth = layer.clientWidth;
    this.#layerHeight = layer.clientHeight;
    this.#cached = true;
    return true;
  }

  /** One read (the label's size; the layer's box when not cached), then one write (its offsets, only when they moved). */
  #place(): void {
    if (!this.#visible()) return;
    if (this.#measure()) {
      const x = this.#x - this.#layerLeft;
      const y = this.#y - this.#layerTop;
      const w = this.#label.offsetWidth;
      const h = this.#label.offsetHeight;
      const left = Math.max(0, Math.round(x + Offset + w > this.#layerWidth ? x - Offset - w : x + Offset));
      const top = Math.max(0, Math.round(y + Offset + h > this.#layerHeight ? y - Offset - h : y + Offset));
      if (left !== this.#left) {
        this.#left = left;
        this.#label.style.setProperty('--msui-hover-x', px(left));
      }
      if (top !== this.#top) {
        this.#top = top;
        this.#label.style.setProperty('--msui-hover-y', px(top));
      }
    }
    this.#show(true);
  }
}

export interface HoverLabelProps {
  model: HoverLabelModel;
  /** Where the pointer is followed (the stage); its listeners attach once it is in the page. */
  area: () => HTMLElement | null;
  /** What the pointer must be over for the label to show (the canvas, the band); default anything in `area`. */
  over?: (target: Element) => boolean;
  /** `corner` (default): put at its offset parent's bottom right (`--msui-hover-label-*`); `cursor`: by the pointer. */
  placement?: HoverAt;
  class?: string;
}

/**
 * The hover label (stage §3.1): the model's lines, by the pointer or in the corner; hidden while a button is held, out
 * of the accessibility tree while hidden. An emptied label stays a frame (a fast move from one item to the next swaps
 * its text in place, no blink), then fades out with its last text.
 */
export function HoverLabel(props: HoverLabelProps) {
  const corner = (props.placement ?? 'corner') === 'corner';
  // The lines in a body of their own: a long label is cut at the body's height and fades there (its edge untouched).
  const body = (<div class={classes.body} />) as HTMLElement;
  const label = (
    <div
      class={`${classes.label} ${classes.out} ${corner ? classes.corner : classes.cursor}${props.class ? ` ${props.class}` : ''}`}
      role="tooltip"
      aria-hidden="true"
      data-placement={corner ? 'corner' : 'cursor'}
    >
      {body}
    </div>
  ) as HTMLElement;
  const follower = new HoverFollower(label, () => label.offsetParent as HTMLElement | null, corner, classes.out);

  // What shows: new lines at once; emptied, only if still empty a frame later (the swap above).
  const shown = new State<readonly HoverLine[]>([]);
  const hiding = frames(() => shown.set([]));
  bind(props.model.lines, (lines) => {
    if (lines.length) {
      hiding.cancel();
      shown.set(lines);
    } else if (shown.value.length) hiding.request();
  });
  bind(shown, (lines) => {
    // Emptied, the last text stays for the fade.
    if (lines.length) body.replaceChildren(...lines.map(lineElement));
    follower.setContent(lines.length > 0);
  });

  // The pointer, on the area, found once the label is in the page: its listeners attached then, once, in this scope
  // (detached while paused).
  const scope = requireScope('HoverLabel');
  queueMicrotask(() => {
    if (scope.disposed) return;
    const area = props.area();
    if (!area) return;
    const over = props.over;
    runInScope(scope, () => {
      on<PointerEvent>(area, 'pointermove', (e) =>
        follower.move(e.clientX, e.clientY, e.buttons !== 0, !over || over(e.target as Element)),
      );
      on(area, 'pointerleave', () => follower.leave());
      if (corner) return;
      // The layer's box, cached by the follower: read again after the stage resizes, anything scrolls, or the pointer
      // comes in (a stage moved without resizing).
      const invalidate = () => follower.invalidate();
      on(area, 'pointerenter', invalidate);
      observeResize(area, invalidate);
      on(area.ownerDocument, 'scroll', invalidate, { capture: true, passive: true });
    });
  });
  return label;
}
