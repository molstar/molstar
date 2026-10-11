/** @jsxImportSource @molstar/plugin-prism-ui/dom */

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
import { bind, type MaybeReactive, onCleanup, toStream } from '../../dom/scope.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import { Icon } from '../icon/icon.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/select.module.css.js';

export type SelectOption = readonly [value: string, label: string];

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
 * picks, Esc closes, and so does a click elsewhere.
 */
export function OptionList(props: {
  options: readonly SelectOption[];
  value: Observable<string | undefined>;
  onPick(value: string): void;
  onClose(): void;
}) {
  const buttons: HTMLButtonElement[] = [];
  const focus = (i: number) => buttons[Math.max(0, Math.min(buttons.length - 1, i))]?.focus();
  const list = (
    <ScrollArea
      class={classes.list}
      contentClass={classes.listContent}
      role="listbox"
      onKeyDown={(e) => {
        const i = buttons.indexOf(document.activeElement as HTMLButtonElement);
        if (e.key === 'ArrowDown') focus(i + 1);
        else if (e.key === 'ArrowUp') focus(i - 1);
        else if (e.key === 'Home') focus(0);
        else if (e.key === 'End') focus(buttons.length - 1);
        else if (e.key === 'Escape') props.onClose();
        else return;
        e.preventDefault();
        e.stopPropagation();
      }}
    >
      {props.options.map(([value, label]) => {
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
        buttons.push(button);
        return button;
      })}
    </ScrollArea>
  ) as HTMLElement;
  // Focus the current option once placed; close on a press outside the list (its scrollbar too) and its button.
  let current: string | undefined;
  bind(props.value, (v) => {
    current = v;
  });
  queueMicrotask(() =>
    focus(
      Math.max(
        0,
        props.options.findIndex((o) => o[0] === current),
      ),
    ),
  );
  const onDown = (e: PointerEvent) => {
    const target = e.target as Element;
    if (list.contains(target) || target.closest?.('[aria-haspopup="listbox"][aria-expanded="true"]')) return;
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
function Floating(props: { anchor: HTMLElement; children: Child }) {
  const el = (
    <div class={classes.floating} attr:popover="manual">
      {props.children}
    </div>
  ) as HTMLElement;
  const place = () => {
    const r = props.anchor.getBoundingClientRect();
    const h = el.offsetHeight;
    const below = window.innerHeight - r.bottom;
    const top = below >= h + 8 || below >= r.top ? r.bottom + 2 : r.top - h - 2;
    // At least as wide as the anchor, as wide as its options need; kept on screen.
    el.style.minWidth = `${r.width}px`;
    el.style.left = `${Math.max(8, Math.min(r.left, window.innerWidth - 8 - el.offsetWidth))}px`;
    el.style.top = `${Math.max(4, top)}px`;
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
  window.addEventListener('scroll', place, true);
  window.addEventListener('resize', place);
  onCleanup(() => {
    window.removeEventListener('scroll', place, true);
    window.removeEventListener('resize', place);
  });
  return el;
}
