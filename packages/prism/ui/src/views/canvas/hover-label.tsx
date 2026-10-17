/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, requireScope, runInScope } from '../../dom/scope.js';
import { frames, on } from '../../dom/scoped.js';
import type { HoverLabelModel, HoverLine, HoverRun } from '../../model/playback.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/hover-label.module.css.js';

const Offset = 14;

/** A line as elements: plain text (the first is the title), or formatted runs. Text only (`textContent`), no markup. */
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

/**
 * The hover label next to the cursor. Its content comes from the model; its position follows the pointer over the
 * canvas directly (one transform write per frame, no state or allocation per move), flipped near the edges and hidden
 * while a button is held (rotating, dragging) or the pointer is outside. `area` is the canvas itself (ShellView marks
 * it `data-msp-canvas`), not what floats over it: over the panel or the chrome there is no label, though the hovered
 * tree row highlights its part on the canvas.
 */
export function HoverLabel(props: {
  model: HoverLabelModel;
  area: () => HTMLElement | null;
  /**
   * What the pointer must be over for the label to show (the canvas, the sequence band): over anything else in the area
   * (a panel floating over it, whose rows highlight atoms and have their own tooltips) there is none. Default: all.
   */
  over?: (target: Element) => boolean;
  /**
   * Where it shows: `corner` (the default), at the area's bottom right, above the corner row (the toolbar: the area's
   * `--msp-corner-room`), where it stays put and covers nothing near the pointer; `cursor`, beside the pointer.
   */
  placement?: 'corner' | 'cursor';
}) {
  const corner = (props.placement ?? 'corner') === 'corner';
  const label = (
    <div
      class={`${classes.label} ${classes.out}${corner ? ` ${classes.corner}` : ''}`}
      role="tooltip"
      aria-hidden="true"
    />
  ) as HTMLElement;

  // What shows: a new label at once; an empty one only if it is still empty a frame later, so a pointer moving fast from
  // one atom to the next (Mol* empties the label between them) swaps the text in place instead of blinking it off.
  // Frames on the page's one tick (`frames`: asked for from bindings and listeners, cancelled with the scope).
  const shown = observableState<readonly HoverLine[]>([]);
  const hiding = frames(() => shown.set([]));
  bind(props.model.lines, (lines) => {
    if (lines.length) {
      hiding.cancel();
      shown.set(lines);
    } else if (shown.value.length) hiding.request();
  });

  let x = 0;
  let y = 0;
  let inside = false;
  let dragging = false;
  let hasContent = false;

  // Hiding waits a frame as emptying does: a fast move crossing a small element over the canvas (or leaving it for a
  // moment) keeps the label rather than blinking it; a drag (rotating) hides it at once. Hidden is the `out` class: it
  // appears at once and fades out (CSS, `--msp-motion-fast`), holding its last text while it fades.
  let isShown = false;
  const show = (on: boolean) => {
    if (isShown === on) return;
    isShown = on;
    label.classList.toggle(classes.out, !on);
    // Out of the accessibility tree while hidden, whatever the CSS (it may not be loaded, or still fading).
    if (on) label.removeAttribute('aria-hidden');
    else label.setAttribute('aria-hidden', 'true');
  };
  const place = () => {
    const visible = inside && !dragging && hasContent;
    if (visible || dragging || !hasContent) {
      leaving.cancel();
      show(visible);
    } else if (isShown) leaving.request();
    if (!visible || corner) return;
    const area = label.offsetParent as HTMLElement | null;
    if (!area) return;
    // The pointer in the layer the label is placed in (the canvas may sit elsewhere in it).
    const rect = area.getBoundingClientRect();
    const px = x - rect.left;
    const py = y - rect.top;
    const w = label.offsetWidth;
    const h = label.offsetHeight;
    const left = px + Offset + w > area.clientWidth ? px - Offset - w : px + Offset;
    const top = py + Offset + h > area.clientHeight ? py - Offset - h : py + Offset;
    label.style.transform = `translate(${Math.max(0, left)}px, ${Math.max(0, top)}px)`;
  };
  const leaving = frames(() => {
    if (!(inside && !dragging && hasContent)) show(false);
  });
  // Placed once a frame, asked for from the listeners and the lines' binding.
  const placing = frames(place);
  const schedule = () => placing.request();

  bind(shown.asReadonly, (lines) => {
    // Emptied, the last text stays for the fade out.
    if (lines.length) label.replaceChildren(...lines.map(lineElement));
    hasContent = lines.length > 0;
    schedule();
  });

  // The pointer is tracked on the canvas area (the element the label is positioned in), found once it is in the page:
  // its listeners attached then, once, in this scope (dom/scoped.ts: detached while it is paused).
  const scope = requireScope('HoverLabel');
  queueMicrotask(() => {
    if (scope.disposed) return;
    const area = props.area();
    if (!area) return;
    const over = props.over;
    runInScope(scope, () => {
      on<PointerEvent>(area, 'pointermove', (e) => {
        x = e.clientX;
        y = e.clientY;
        inside = !over || over(e.target as Element);
        dragging = e.buttons !== 0;
        schedule();
      });
      on(area, 'pointerleave', () => {
        inside = false;
        schedule();
      });
    });
  });
  return label;
}
