/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isState, type ReadonlyState, State } from '@molstar/core/reactive';
import {
  activity,
  effect,
  type Child,
  type MaybeReactive,
  Motion,
  onCleanup,
  Show,
  setAttr,
  toStream,
  uniqueId,
  watch,
} from '../../dom/index.js';
import { held } from '../../state.js';
import { Icon } from '../icon/icon.js';
import { classes } from './_generated/expandable.module.css.js';

export interface ExpandableProps {
  /** The header's text (one line: it fades past its width). */
  title: Child;
  /**
   * Whether it is open. A `State` is set when the header is pressed (the host keeps it, e.g. to remember it); a plain
   * value is how it starts (it keeps its own); a read-only state or a stream is followed (`onToggle` updates its source).
   * Default: closed.
   */
  open?: State<boolean> | MaybeReactive<boolean>;
  /** Told the new state when the header is pressed (after a `State` or a plain value took it). */
  onToggle?(open: boolean): void;
  /** What trails the title in the header, before the chevron (a count, a status): not part of the button. */
  trailing?: Child;
  /**
   * `well` (default): an open panel is a well (`--msui-surface-sunken`) holding its header and body, the header hovers
   * as a row; `flat`: no backgrounds, a header a row high (inside a surface that already is one: an inspector).
   */
  variant?: 'well' | 'flat';
  /** The heading level of the header (an ARIA heading around its button, WAI-ARIA accordion). Default 3. */
  level?: number;
  /**
   * The body's content. Rendered the first time it opens (`lazy`, the default) or at once; then kept while it closes,
   * inert and paused (`activity`: its bindings hold no subscriptions), so what lives in it (a scroll, a draft) stays.
   */
  children: () => Child;
  /** Render the body at once, not on first open. */
  eager?: boolean;
  class?: string;
}

/**
 * An expandable panel: a header (its title, what trails it, a chevron that turns) that opens and closes a body under it.
 * Opened, it is a quiet well (`variant="well"`) holding both, the header its top; closed, only the header shows, a row
 * that hovers. Opening grows the body in, closing plays it back (styles §7.1, `Motion.close`) before it hides. The body
 * stays mounted once rendered, inert and paused while closed (dom §6). Stack several in a grid (`ExpandableList`): each
 * keeps its own state. WAI-ARIA disclosure in an accordion's markup: a button in a heading, `aria-expanded`,
 * `aria-controls` the body, a region named by the header.
 */
export function Expandable(props: ExpandableProps) {
  // A writable state or a plain value takes the toggle here; a read-only source is the caller's to update.
  const own =
    props.open instanceof State
      ? (props.open as State<boolean>)
      : !isState(props.open) && !toStream(props.open)
        ? new State(!!props.open)
        : undefined;
  const source: ReadonlyState<boolean> = own ?? held<boolean>(props.open as MaybeReactive<boolean>, false);
  // What the header says (its chevron turns at once), ahead of the body while it closes.
  const expanded = new State(source.value);
  // Whether the body shows: off only once a closing played.
  const shown = new State(source.value);
  // Rendered once it first opened (or at once, `eager`).
  const rendered = new State(!!props.eager || source.value);
  // Opened from here on: the body is revealed (styles §7.1); the first render stays still.
  const animate = new State(false);
  const id = uniqueId('msui-expandable');
  let body!: HTMLElement;
  let header!: HTMLElement;
  let closing: Animation | undefined;

  const stopClosing = () => {
    closing?.cancel();
    closing = undefined;
  };
  onCleanup(stopClosing);

  // The panel as a function of its source: open at once, close after the closing played (where it plays).
  const apply = (open: boolean) => {
    if (open === expanded.value) return;
    animate.set(true);
    stopClosing();
    expanded.set(open);
    if (open) {
      rendered.set(true);
      shown.set(true);
      return;
    }
    const anim = shown.value ? Motion.close?.(body, 'grow') : undefined;
    if (!anim) {
      shown.set(false);
      return;
    }
    closing = anim;
    anim.finished.then(
      () => {
        if (closing !== anim) return;
        closing = undefined;
        shown.set(false);
        anim.cancel();
      },
      () => {},
    );
  };
  // Paused with its scope (dom §6): resumed, it follows the value it missed.
  effect(source, apply);

  const toggle = () => {
    const next = !expanded.value;
    if (own) own.set(next);
    else apply(next);
    props.onToggle?.(next);
  };

  const root = (
    <section
      class={props.class ? `${classes.expandable} ${props.class}` : classes.expandable}
      data-variant={props.variant ?? 'well'}
    >
      <div class={classes.head} role="heading" aria-level={String(props.level ?? 3)}>
        <button
          type="button"
          id={`${id}-header`}
          class={classes.header}
          aria-controls={`${id}-body`}
          ref={(el: HTMLElement) => (header = el)}
          onClick={toggle}
        >
          <span class={classes.title}>{props.title}</span>
          {props.trailing != null ? <span class={classes.trailing}>{props.trailing}</span> : null}
          <Icon name="chevron-down" class={classes.chevron} />
        </button>
      </div>
      <div
        class={classes.body}
        id={`${id}-body`}
        role="region"
        aria-labelledby={`${id}-header`}
        ref={(el: HTMLElement) => (body = el)}
      >
        {/* Kept once rendered; closed, inert and paused (`activity`). Its `hidden` is written outside it (below), or it
            would pause with it. */}
        <Show when={rendered}>{() => activity(shown.asReadonly, props.children)}</Show>
      </div>
    </section>
  ) as HTMLElement;

  watch([expanded, shown, animate], () => {
    setAttr(header, 'aria-expanded', String(expanded.value));
    setAttr(root, 'data-open', shown.value);
    setAttr(body, 'hidden', !shown.value);
    // Revealed each time it shows again (an animation restarts when its element leaves `display: none`).
    setAttr(body, 'data-reveal', animate.value);
  });
  return root;
}

/** Expandable panels stacked: a unit apart, each its own. */
export function ExpandableList(props: { children: Child; class?: string }) {
  return <div class={props.class ? `${classes.list} ${props.class}` : classes.list}>{props.children}</div>;
}
