/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { For } from '../../dom/control-flow.js';
import { bind, requireScope, runInScope } from '../../dom/scope.js';
import { frames, on } from '../../dom/scoped.js';
import type { HoverLabelModel } from '../../model/playback.js';
import { derive } from '../../reactive/state.js';
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
          <span class={derive(index, (i) => (i === 0 ? classes.primary : classes.secondary))}>{line}</span>
        )}
      </For>
    </div>
  ) as HTMLElement;

  let x = 0;
  let y = 0;
  let inside = false;
  let dragging = false;
  let hasContent = false;

  const place = () => {
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
  // Asked for from the pointer's listeners and the lines' binding, outside the scope (where `frame()` cannot be): a
  // frame callback by hand, in the page's one frame tick, cancelled with the scope.
  // Placed once a frame on the page's one tick, asked for from the listeners and the lines' binding (`frames`).
  const placing = frames(place);
  const schedule = () => {
    placing.request();
  };

  bind(props.model.lines, (lines) => {
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
