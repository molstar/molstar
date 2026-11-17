/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Child, frame, type MaybeReactive } from '../../dom/index.js';
import { Floating } from '../select/select.js';
import { classes } from './_generated/popover.module.css.js';
import { dismissable } from './dismiss.js';

/**
 * A panel of controls in a popover from its button (a header's settings): a raised surface placed as a menu is
 * (`MenuPopover`), a dialog rather than a menu. Escape (the focus back on the button) or a press outside it and its
 * button closes it; a press or Escape in a layer it opened (a select's options) is that layer's, as is an Escape with
 * the focus on another surface. The focus goes to its first control once it shows. Its sides pad by the surface's pad
 * (styles §5): what is in it (a form, rows) pads the rest of the inset.
 */
export function Popover(props: {
  anchor: HTMLElement;
  onClose(refocus?: boolean): void;
  children: Child;
  /** A heading at its top (a header, styles §4), which also names it. */
  title?: MaybeReactive<string>;
  'aria-label'?: string;
  /** Its edge that lines up with the anchor's (default `end`: a button at the right), where it has no column. */
  align?: 'start' | 'end';
  /**
   * The column it covers: its left edge and width. Default: the surface its anchor is in (`.msui-surface`), so it is
   * never wider than the panel it opens from; `null` for none (as wide as its content, at most 340px, by the anchor).
   */
  column?: HTMLElement | null;
}) {
  let panel: HTMLElement | undefined;
  // One rule with the other overlays (dismiss.ts): a layer opened from it (a select's options) takes its own presses and
  // keys first; Escape with the focus on another surface is that surface's. A dialog: Tab stays in it.
  dismissable({
    layer: () => panel,
    anchor: props.anchor,
    close: (reason) => (reason === 'outside' ? props.onClose() : props.onClose(true)),
  });
  return (
    <Floating
      anchor={props.anchor}
      align={props.align ?? 'end'}
      column={
        props.column === undefined
          ? (props.anchor.closest<HTMLElement>('.msui-surface') ?? undefined)
          : (props.column ?? undefined)
      }
    >
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
