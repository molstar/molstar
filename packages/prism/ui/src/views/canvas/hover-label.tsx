/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, onCleanup } from '../../dom/scope.js';
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
    <div class={corner ? `${classes.label} ${classes.corner}` : classes.label} hidden={true} role="tooltip" />
  ) as HTMLElement;

  // What shows: a new label at once; an empty one only if it is still empty a frame later, so a pointer moving fast from
  // one atom to the next (Mol* empties the label between them) swaps the text in place instead of blinking it off.
  const shown = observableState<readonly HoverLine[]>([]);
  let hideFrame = 0;
  bind(props.model.lines, (lines) => {
    if (lines.length) {
      if (hideFrame) cancelAnimationFrame(hideFrame);
      hideFrame = 0;
      shown.set(lines);
    } else if (!hideFrame && shown.value.length) {
      hideFrame = requestAnimationFrame(() => {
        hideFrame = 0;
        shown.set([]);
      });
    }
  });
  onCleanup(() => {
    if (hideFrame) cancelAnimationFrame(hideFrame);
    if (leaveFrame) cancelAnimationFrame(leaveFrame);
  });

  let x = 0;
  let y = 0;
  let inside = false;
  let dragging = false;
  let hasContent = false;
  let frame = 0;

  // Hiding waits a frame as emptying does: a fast move crossing a small element over the canvas (or leaving it for a
  // moment) keeps the label rather than blinking it; a drag (rotating) hides it at once.
  let leaveFrame = 0;
  const place = () => {
    frame = 0;
    const visible = inside && !dragging && hasContent;
    if (visible || dragging || !hasContent) {
      if (leaveFrame) cancelAnimationFrame(leaveFrame);
      leaveFrame = 0;
      if (label.hidden === visible) label.hidden = !visible;
    } else if (!label.hidden && !leaveFrame) {
      leaveFrame = requestAnimationFrame(() => {
        leaveFrame = 0;
        if (!(inside && !dragging && hasContent)) label.hidden = true;
      });
    }
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
  const schedule = () => {
    if (!frame) frame = requestAnimationFrame(place);
  };

  bind(shown.asReadonly, (lines) => {
    label.replaceChildren(...lines.map(lineElement));
    hasContent = lines.length > 0;
    schedule();
  });

  // The pointer is tracked on the canvas area (the element the label is positioned in).
  queueMicrotask(() => {
    const area = props.area();
    if (!area) return;
    const over = props.over;
    const onMove = (e: PointerEvent) => {
      x = e.clientX;
      y = e.clientY;
      inside = !over || over(e.target as Element);
      dragging = e.buttons !== 0;
      schedule();
    };
    const onLeave = () => {
      inside = false;
      schedule();
    };
    area.addEventListener('pointermove', onMove);
    area.addEventListener('pointerleave', onLeave);
    onCleanupArea = () => {
      area.removeEventListener('pointermove', onMove);
      area.removeEventListener('pointerleave', onLeave);
    };
  });
  let onCleanupArea: () => void = () => {};
  onCleanup(() => {
    onCleanupArea();
    if (frame) cancelAnimationFrame(frame);
  });
  return label;
}
