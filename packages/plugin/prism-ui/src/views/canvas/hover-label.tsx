/** @jsxImportSource @molstar/plugin-prism-ui/dom */

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
 * while a button is held (rotating, dragging) or the pointer is outside.
 */
export function HoverLabel(props: { model: HoverLabelModel; area: () => HTMLElement | null }) {
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
    const w = label.offsetWidth;
    const h = label.offsetHeight;
    const left = x + Offset + w > area.clientWidth ? x - Offset - w : x + Offset;
    const top = y + Offset + h > area.clientHeight ? y - Offset - h : y + Offset;
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
    const onMove = (e: PointerEvent) => {
      const rect = area.getBoundingClientRect();
      x = e.clientX - rect.left;
      y = e.clientY - rect.top;
      inside = true;
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
