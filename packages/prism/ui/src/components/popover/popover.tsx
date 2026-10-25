/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from '../../dom/jsx-types.js';
import { eventTarget, isInside } from '../../dom/events.js';
import type { MaybeReactive } from '../../dom/scope.js';
import { frame, on } from '../../dom/scoped.js';
import { Floating } from '../select/select.js';
import { classes } from './_generated/popover.module.css.js';

/** Whether an event happened in a floating layer other than `own` (a select's options opened from inside it). */
function inOtherLayer(e: Event, own: HTMLElement | undefined): boolean {
  const path = typeof e.composedPath === 'function' ? e.composedPath() : [];
  return path.some((n) => n instanceof HTMLElement && n !== own && n.hasAttribute('popover') && !own?.contains(n));
}

/**
 * A panel of controls in a popover from its button (a header's sliders): raised and placed as a menu is
 * (`MenuPopover`), a dialog rather than a menu. Escape (the focus back on the button) or a press outside it and its
 * button closes it; a press or Escape in a layer it opened (a select's options) is that layer's, as is an Escape with
 * the focus on another surface. The focus goes to its first control once it shows.
 */
export function Popover(props: {
  anchor: HTMLElement;
  onClose(refocus?: boolean): void;
  children: Child;
  /** A heading at its top, which also names it. */
  title?: MaybeReactive<string>;
  'aria-label'?: string;
  /** Its edge that lines up with the anchor's (default `end`: a button at the right). */
  align?: 'start' | 'end';
}) {
  let panel: HTMLElement | undefined;
  const layer = () => panel?.closest<HTMLElement>('[popover]') ?? undefined;
  // Escape closes it before the panel it is in (the document's handlers run in capture order); one meant for a layer
  // opened from it (the focus there) is that layer's.
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented) return;
    if (!isInside(e, panel) && inOtherLayer(e, layer())) return;
    // Escape meant for another surface (a dialog, the command menu: the focus there) is not this one's
    const target = eventTarget(e);
    const elsewhere =
      target instanceof Element &&
      target !== document.body &&
      target !== document.documentElement &&
      !isInside(e, panel) &&
      !isInside(e, props.anchor);
    if (elsewhere) return;
    e.preventDefault();
    props.onClose(true);
  };
  const onDown = (e: PointerEvent) => {
    if (isInside(e, panel) || isInside(e, props.anchor) || inOtherLayer(e, layer())) return;
    props.onClose();
  };
  on(document, 'keydown', onKey, { capture: true });
  on(document, 'pointerdown', onDown, { capture: true });
  return (
    <Floating anchor={props.anchor} align={props.align ?? 'end'}>
      <div
        class={classes.popover}
        role="dialog"
        aria-label={props['aria-label'] ?? props.title}
        ref={(el: HTMLElement) => {
          panel = el;
          frame(() => el.querySelector<HTMLElement>('button:not(:disabled), input, select, [tabindex="0"]')?.focus());
        }}
      >
        {props.title !== undefined ? <div class={classes.title}>{props.title}</div> : null}
        {props.children}
      </div>
    </Floating>
  );
}
