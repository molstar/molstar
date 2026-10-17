/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { activeElement, eventTarget, isInside } from '../../dom/events.js';
import { bind, type MaybeReactive, onCleanup, toStream } from '../../dom/scope.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import { Icon } from '../icon/icon.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/select.module.css.js';

import { derive } from '../../reactive/state.js';
import { requestFrameCallback } from '../../reactive/frame-state.js';
/** An option: its value, its label, and the group it is listed under (Mol*'s option categories), if any. */
export type SelectOption = readonly [value: string, label: string, group?: string];

export interface SelectProps {
  options: readonly SelectOption[];
  value?: MaybeReactive<string>;
  onChange?(value: string): void;
  size?: 'sm' | 'md';
  class?: string;
  title?: string;
  'aria-label'?: string;
  /**
   * The options float over what is below instead of pushing it down (a short settings row, where moving the layout
   * would be jarring); in forms they open inline, under the row, as Mol*'s do.
   */
  over?: boolean;
}

/**
 * A choice that expands inline (Mol*'s action menu): the button shows the current option; opened, the options list
 * right under it, the current one marked. Looks and works the same on every platform; arrows move, Enter picks, Esc
 * closes. Forms put the list under the whole row (`SelectButton` + `OptionList`).
 */
export function Select(props: SelectProps) {
  const open = observableState(false);
  const value = toValue(props.value);
  let anchor!: HTMLElement;
  return (
    <span
      class={`${classes.select} ${props.class ?? ''}`}
      ref={(el) => {
        anchor = el;
      }}
    >
      <SelectButton
        label={value.pipe(map((v) => props.options.find((o) => o[0] === v)?.[1] ?? ''))}
        open={open}
        size={props.size}
        title={props.title}
        aria-label={props['aria-label']}
      />
      <Show when={open}>
        {() => {
          const list = (
            <OptionList
              options={props.options}
              value={value}
              onPick={(v) => {
                open.set(false);
                props.onChange?.(v);
              }}
              onClose={() => open.set(false)}
              floating={props.over}
              aria-label={props['aria-label'] ?? props.title}
            />
          );
          return props.over ? <Floating anchor={anchor}>{list}</Floating> : list;
        }}
      </Show>
    </span>
  );
}

/** The closed choice: the current option's label and a chevron; toggles `open`. */
export function SelectButton(props: {
  label: Observable<string>;
  open: ObservableState<boolean>;
  size?: 'sm' | 'md';
  title?: string;
  'aria-label'?: string;
}) {
  return (
    <button
      type="button"
      class={classes.button}
      data-size={props.size ?? 'md'}
      aria-haspopup="listbox"
      aria-expanded={props.open}
      aria-label={props['aria-label']}
      title={props.title}
      onClick={() => props.open.set(!props.open.value)}
    >
      <span class={classes.label}>{props.label}</span>
      <Icon name="chevron-down" size={14} class={classes.chevron} />
    </button>
  );
}

/**
 * The options, inline: the current one marked; focused on open (the current option), arrows move, Enter or a click
 * picks, Esc closes, and so does a click elsewhere. Options with a group are listed under it as v5 lists Mol*'s
 * categories: each group a fold, the current option's open, the others closed; options without one come first.
 */
export function OptionList(props: {
  options: readonly SelectOption[];
  value: Observable<string | undefined>;
  onPick(value: string): void;
  onClose(): void;
  /** Floating over the page (a popover): the wheel stays in it. Inline (default), at its ends it scrolls the panel. */
  floating?: boolean;
  /** The list's accessible name (a listbox has one): what the choice is, as its button's. */
  'aria-label'?: string;
}) {
  // The current value when opened: its group starts open.
  let current: string | undefined;
  bind(props.value, (v) => {
    current = v;
  });
  // Each option's button by its value (the render order is the loose options, then each group's).
  const buttons = new Map<string, HTMLButtonElement>();
  // What the arrows move over: the options and the group heads that show (not inside a closed group).
  const navigable = () =>
    [...list.querySelectorAll<HTMLElement>('[role="option"], [data-group-head]')].filter((b) => !b.closest('[hidden]'));
  const focus = (i: number) => {
    const all = navigable();
    all[Math.max(0, Math.min(all.length - 1, i))]?.focus();
  };
  const option = ([value, label]: SelectOption) => {
    const button = (
      <button
        type="button"
        role="option"
        class={classes.option}
        aria-selected={props.value.pipe(
          map((v) => v === value),
          distinctUntilChanged(),
        )}
        onClick={() => props.onPick(value)}
      >
        <span>{label}</span>
        <Icon name="check" size={14} class={classes.check} />
      </button>
    ) as HTMLButtonElement;
    buttons.set(value, button);
    return button;
  };
  const groups = new Map<string, SelectOption[]>();
  const loose: SelectOption[] = [];
  for (const o of props.options) {
    if (!o[2]) loose.push(o);
    else if (groups.has(o[2])) groups.get(o[2])!.push(o);
    else groups.set(o[2], [o]);
  }
  const group = (name: string, options: readonly SelectOption[]) => {
    const open = observableState(options.some((o) => o[0] === current));
    return (
      <div class={classes.group} role="group" aria-label={name}>
        <button
          type="button"
          class={classes.groupHead}
          data-group-head=""
          aria-expanded={open}
          onClick={() => open.set(!open.value)}
        >
          <Icon name="chevron-right" size={13} class={classes.groupChevron} />
          <span>{name}</span>
        </button>
        <div class={classes.groupBody} hidden={derive(open.asReadonly, (o) => !o)}>
          {options.map(option)}
        </div>
      </div>
    );
  };
  const list = (
    <ScrollArea
      class={classes.list}
      contentClass={classes.listContent}
      chain={!props.floating}
      role="listbox"
      aria-label={props['aria-label']}
      onKeyDown={(e) => {
        const all = navigable();
        const i = all.indexOf(activeElement(list) as HTMLElement);
        if (e.key === 'ArrowDown') focus(i + 1);
        else if (e.key === 'ArrowUp') focus(i - 1);
        else if (e.key === 'Home') focus(0);
        else if (e.key === 'End') focus(all.length - 1);
        else if (e.key === 'Escape') props.onClose();
        else return;
        e.preventDefault();
        e.stopPropagation();
      }}
    >
      {loose.map(option)}
      {[...groups].map(([name, options]) => group(name, options))}
    </ScrollArea>
  ) as HTMLElement;
  // Focus the current option once placed (else the first that shows); close on a press outside the list (its scrollbar
  // too) and its button.
  queueMicrotask(() => {
    const at = current === undefined ? undefined : buttons.get(current);
    if (at && !at.closest('[hidden]')) at.focus();
    else focus(0);
  });
  const onDown = (e: PointerEvent) => {
    const target = eventTarget(e) as Element;
    if (isInside(e, list) || target.closest?.('[aria-haspopup="listbox"][aria-expanded="true"]')) return;
    props.onClose();
  };
  document.addEventListener('pointerdown', onDown, true);
  onCleanup(() => document.removeEventListener('pointerdown', onDown, true));
  return list;
}

function toValue(value: MaybeReactive<string> | undefined): Observable<string | undefined> {
  const stream = toStream(value as MaybeReactive<string | undefined>);
  if (stream) return stream as Observable<string | undefined>;
  const state = observableState<string | undefined>(value as string | undefined);
  return state.asReadonly;
}

/**
 * What floats over the layout, under its anchor and at least as wide (above it when there is no room below): in the top layer
 * (a manual popover), so no panel's clipping or fade cuts it; it follows the anchor as things scroll or resize.
 */
/**
 * A menu in a popover from its button (a row's ⋯): raised, as a select's options; Escape (the focus back on the button)
 * or a press outside it and its button closes it; the focus goes to its first enabled entry once it shows.
 */
export function MenuPopover(props: {
  anchor: HTMLElement;
  onClose(refocus?: boolean): void;
  children: Child;
  'aria-label'?: string;
}) {
  let menu: HTMLElement | undefined;
  // Escape closes the menu first, before the panel it is in (the document's handlers run in capture order).
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented) return;
    e.preventDefault();
    props.onClose(true);
  };
  const onDown = (e: PointerEvent) => {
    if (isInside(e, menu) || isInside(e, props.anchor)) return;
    props.onClose();
  };
  document.addEventListener('keydown', onKey, true);
  document.addEventListener('pointerdown', onDown, true);
  onCleanup(() => {
    document.removeEventListener('keydown', onKey, true);
    document.removeEventListener('pointerdown', onDown, true);
  });
  return (
    <Floating anchor={props.anchor} align="end">
      <div
        class={classes.menuPopover}
        role="menu"
        aria-label={props['aria-label']}
        ref={(el: HTMLElement) => {
          menu = el;
          // Into it once it shows (the popover opens on the next microtask).
          requestFrameCallback(() => el.querySelector<HTMLElement>('button:not(:disabled)')?.focus());
        }}
      >
        {props.children}
      </div>
    </Floating>
  );
}

export function Floating(props: {
  anchor: HTMLElement;
  children: Child;
  /** Its edge that lines up with the anchor's: `start` (left, the default) or `end` (a menu from a button at the right). */
  align?: 'start' | 'end';
}) {
  const el = (
    <div class={classes.floating} attr:popover="manual">
      {props.children}
    </div>
  ) as HTMLElement;
  // At least as wide as the anchor, as wide as its options need: written only when the anchor's width changes, so
  // a scroll or resize reads the layout once and then writes (no layout forced after a write).
  let minWidth = -1;
  const place = () => {
    const r = props.anchor.getBoundingClientRect();
    if (r.width !== minWidth) {
      minWidth = r.width;
      el.style.minWidth = `${r.width}px`;
    }
    const h = el.offsetHeight;
    const w = el.offsetWidth;
    const below = window.innerHeight - r.bottom;
    const top = below >= h + 8 || below >= r.top ? r.bottom + 2 : r.top - h - 2;
    // Kept on screen.
    const left = props.align === 'end' ? r.right - w : r.left;
    el.style.left = `${Math.max(8, Math.min(left, window.innerWidth - 8 - w))}px`;
    el.style.top = `${Math.max(4, top)}px`;
  };
  // A scroll of its own list does not move the anchor (the listener captures every element's scroll).
  const onScroll = (e: Event) => {
    if (!isInside(e, el)) place();
  };
  queueMicrotask(() => {
    if (!el.isConnected) return;
    el.showPopover?.();
    place();
    // The list focused its current option while still hidden (that does nothing): now that it shows, again.
    (
      el.querySelector<HTMLElement>('[aria-selected="true"]') ?? el.querySelector<HTMLElement>('[role=option]')
    )?.focus();
  });
  window.addEventListener('scroll', onScroll, true);
  window.addEventListener('resize', place);
  onCleanup(() => {
    window.removeEventListener('scroll', onScroll, true);
    window.removeEventListener('resize', place);
  });
  return el;
}
