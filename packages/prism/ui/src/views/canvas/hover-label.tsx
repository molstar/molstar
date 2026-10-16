/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { For } from '../../dom/control-flow.js';
import { bind, onCleanup } from '../../dom/scope.js';
import type { HoverLabelModel } from '../../model/playback.js';
import { classes } from './_generated/hover-label.module.css.js';

const Offset = 14;

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
}) {
  const label = (
    <div class={classes.label} hidden={true} role="tooltip">
      <For each={props.model.lines}>
        {(line, index) => (
          <span class={index.pipe(map((i) => (i === 0 ? classes.primary : classes.secondary)))}>{line}</span>
        )}
      </For>
    </div>
  ) as HTMLElement;

  let x = 0;
  let y = 0;
  let inside = false;
  let dragging = false;
  let hasContent = false;
  let frame = 0;

  const place = () => {
    frame = 0;
    const visible = inside && !dragging && hasContent;
    if (label.hidden === visible) label.hidden = !visible;
    if (!visible) return;
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

  bind(props.model.lines, (lines) => {
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
