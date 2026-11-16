/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, KeyedSelector, type ReadonlyState, type State } from '@molstar/core/reactive';
import {
  activeElement,
  type Child,
  eventTarget,
  frame,
  isInside,
  type MaybeReactive,
  observeResize,
  on,
  Show,
  setAttr,
  watch,
} from '../../dom/index.js';
import { held, state } from '../../state.js';
import { Icon } from '../icon/icon.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/select.module.css.js';

/** An option: its value, its label, and the group it is listed under (Mol*'s option categories), if any. */
export type SelectOption = readonly [value: string, label: string, group?: string];

export interface SelectProps {
  options: readonly SelectOption[];
  value?: MaybeReactive<string>;
  onChange?(value: string): void;
  /** As `Button`'s sizes: `md` (default) a row's height (24px), on its own; `sm` a control (20px), in a form row. */
  size?: 'sm' | 'md';
  class?: string;
  /** Tooltip (`data-tip`); the accessible name when there is no `aria-label`. */
  title?: string;
  'aria-label'?: string;
  /**
   * Where the options open. `true` (default): over the field, the current option on the field's text, in the top layer:
   * the layout never moves. `false`: inline, under the field, pushing what follows down (a list meant to stay in the
   * flow of a panel).
   */
  over?: boolean;
}

/**
 * A choice: the button shows the current option; opened, the options list over it (the current one on the button's
 * text, marked), or inline under it (`over={false}`). Looks and works the same on every platform; arrows move, Enter
 * picks, Esc closes. A form composes its own (`SelectButton` + `OptionList` in a `Floating` with `placement="over"`).
 */
export function Select(props: SelectProps) {
  const open = state(false);
  const over = props.over ?? true;
  const value = held<string | undefined>(props.value, undefined);
  const label = derive(value, (v) => props.options.find((o) => o[0] === v)?.[1] ?? '');
  let anchor!: HTMLElement;
  return (
    <span
      class={props.class ? `${classes.select} ${props.class}` : classes.select}
      ref={(el) => {
        anchor = el;
      }}
    >
      <SelectButton
        label={label}
        open={open}
        size={props.size}
        title={props.title}
        aria-label={props['aria-label'] ?? props.title}
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
              floating={over}
              aria-label={props['aria-label'] ?? props.title}
            />
          );
          return over ? (
            <Floating anchor={anchor} placement="over">
              {list}
            </Floating>
          ) : (
            list
          );
        }}
      </Show>
    </span>
  );
}

/**
 * The closed choice: the current option's label and a chevron; toggles `open`. The label is a string or a state of one,
 * or a node the caller writes (a text node a form row's watch sets).
 */
export function SelectButton(props: {
  label: MaybeReactive<string> | Node;
  open: State<boolean>;
  size?: 'sm' | 'md';
  title?: string;
  'aria-label'?: string;
}) {
  const button = (
    <button
      type="button"
      class={classes.button}
      data-size={props.size ?? 'md'}
      aria-haspopup="listbox"
      aria-label={props['aria-label']}
      data-tip={props.title}
      onClick={() => props.open.set(!props.open.value)}
    >
      <span class={classes.label}>{props.label as Child}</span>
      <Icon name="chevron-down" class={classes.chevron} />
    </button>
  ) as HTMLButtonElement;
  watch([props.open], () => setAttr(button, 'aria-expanded', String(props.open.value)));
  return button;
}

/**
 * The options, inline: the current one marked; focused on open (the current option), arrows move, Enter or a click
 * picks, Esc closes, and so does a press elsewhere. Options with a group are listed under it: each group a fold, the
 * current option's open, the others closed; options without one come first.
 */
export function OptionList(props: {
  options: readonly SelectOption[];
  value: ReadonlyState<string | undefined>;
  onPick(value: string): void;
  onClose(): void;
  /** Floating over the page (a popover): the wheel stays in it. Inline (default), at its ends it scrolls the panel. */
  floating?: boolean;
  /** The list's accessible name (a listbox has one): what the choice is, as its button's. */
  'aria-label'?: string;
}) {
  // The current value when opened: its group starts open, its option takes the focus.
  const current = props.value.value;
  // The marked option: one selector over the value, so a change rewrites two options, not all.
  const selected = new KeyedSelector(props.value);
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
      <button type="button" role="option" class={classes.option} onClick={() => props.onPick(value)}>
        <span class={classes.optionLabel}>{label}</span>
        <Icon name="check" class={classes.check} />
      </button>
    ) as HTMLButtonElement;
    const on = selected.key(value);
    watch([on], () => setAttr(button, 'aria-selected', String(on.value)));
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
    const open = state(options.some((o) => o[0] === current));
    const head = (
      <button type="button" class={classes.groupHead} data-group-head="" onClick={() => open.set(!open.value)}>
        <span>{name}</span>
        <Icon name="chevron-down" class={classes.groupChevron} />
      </button>
    ) as HTMLButtonElement;
    const body = (<div class={classes.groupBody}>{options.map(option)}</div>) as HTMLElement;
    watch([open], () => {
      setAttr(head, 'aria-expanded', String(open.value));
      body.hidden = !open.value;
    });
    return (
      <div class={classes.group} role="group" aria-label={name}>
        {head}
        {body}
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
  on(document, 'pointerdown', onDown, { capture: true });
  return list;
}

/**
 * A menu in a popover from its button (a row's ⋯): raised, as a select's floating options; Escape (the focus back on
 * the button) or a press outside it and its button closes it; the focus goes to its first enabled entry once it shows.
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
  on(document, 'keydown', onKey, { capture: true });
  on(document, 'pointerdown', onDown, { capture: true });
  return (
    <Floating anchor={props.anchor} align="end">
      <div
        class={classes.menuPopover}
        role="menu"
        aria-label={props['aria-label']}
        ref={(el: HTMLElement) => {
          menu = el;
          // Into it once it shows (the popover opens on the next microtask).
          frame(() => el.querySelector<HTMLElement>('button:not(:disabled)')?.focus());
        }}
      >
        {props.children}
      </div>
    </Floating>
  );
}

/**
 * What floats over the layout in the top layer (a manual popover), so no panel's clipping or fade cuts it, at least as
 * wide as its anchor; it follows the anchor as things scroll or resize. `below` (default): under the anchor (above it
 * when there is no room below). `over`: a select's options over their field, the current option (`aria-selected`)
 * centered on the anchor so its text lies on the field's (the option rows write at the control inset, as the field);
 * set once, when it opens: it then keeps that place by the anchor (a fold opened in it grows it downward). `above`:
 * over the anchor (under it when there is no room above). `left`, `right`: beside the anchor (a toolbar's menu, a
 * submenu), flipped to the other side where it would leave the window; not widened to the anchor.
 */
export function Floating(props: {
  anchor: HTMLElement;
  children: Child;
  /**
   * Its edge that lines up with the anchor's: `start` (left, or top beside it; the default) or `end` (a menu from a
   * button at the right).
   */
  align?: 'start' | 'end';
  placement?: 'below' | 'over' | 'above' | 'left' | 'right';
  /** Its distance from the anchor, px (default 2): `below`, `above`, `left`, `right`. */
  gap?: number;
  /**
   * How far its aligned edge lies past the anchor's, outward, px (default 0): a menu whose rows' icons line up with
   * the button's (the menu's padding).
   */
  shift?: number;
  /**
   * The column it belongs to (`below`, `above`): it takes the column's left edge and width instead of lining up with
   * the anchor, so a panel opened from a header covers its panel and nothing past it.
   */
  column?: HTMLElement;
}) {
  const placement = props.placement ?? 'below';
  const beside = placement === 'left' || placement === 'right';
  const gap = props.gap ?? 2;
  const shift = props.shift ?? 0;
  const el = (
    <div
      class={classes.floating}
      attr:popover="manual"
      data-placement={placement}
      data-column={props.column ? '' : undefined}
    >
      {props.children}
    </div>
  ) as HTMLElement;
  // At least as wide as the anchor, as wide as its options need: written only when the anchor's width changes, so a
  // scroll or resize reads the layout once and then writes (no layout forced after a write).
  let minWidth = -1;
  let columnWidth = -1;
  // `over`: its top from the anchor's, measured once on opening.
  let overOffset: number | undefined;
  const place = () => {
    const r = props.anchor.getBoundingClientRect();
    const column = !beside && props.column?.isConnected ? props.column.getBoundingClientRect() : undefined;
    if (column && column.width !== columnWidth) {
      columnWidth = column.width;
      el.style.width = `${column.width}px`;
    } else if (!column && !beside && r.width !== minWidth) {
      minWidth = r.width;
      el.style.minWidth = `${r.width}px`;
    }
    const h = el.offsetHeight;
    const w = el.offsetWidth;
    let top: number;
    let left: number;
    if (beside) {
      // Its side, or the other where it would leave the window and the other has more room.
      const right = window.innerWidth - r.right;
      const onRight = placement === 'right' ? right >= w + 8 || right >= r.left : !(r.left >= w + 8 || r.left >= right);
      left = onRight ? r.right + gap : r.left - gap - w;
      top = props.align === 'end' ? r.bottom - h + shift : r.top - shift;
    } else {
      if (placement === 'over') {
        if (overOffset === undefined) overOffset = measureOver(el, r);
        top = r.top + overOffset;
      } else {
        const below = window.innerHeight - r.bottom;
        const fitsBelow = below >= h + 8 || below >= r.top;
        const fitsAbove = r.top >= h + 8 || r.top >= below;
        top = (placement === 'above' ? !fitsAbove : fitsBelow) ? r.bottom + gap : r.top - h - gap;
      }
      left = column ? column.left : props.align === 'end' ? r.right - w + shift : r.left - shift;
    }
    // Kept on screen: one that grew past the bottom (a submenu opened in place) moves up, over its anchor if it must.
    el.style.left = `${Math.max(8, Math.min(left, window.innerWidth - 8 - w))}px`;
    el.style.top = `${Math.max(4, Math.min(top, window.innerHeight - 8 - h))}px`;
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
    (el.querySelector<HTMLElement>('[aria-selected="true"]') ?? el.querySelector<HTMLElement>('[role=option]'))?.focus({
      preventScroll: true,
    });
  });
  on(window, 'scroll', onScroll, { capture: true });
  on(window, 'resize', place);
  // Its content grew or shrank (a menu's entry opened its sub-entries): placed again, after that layout.
  observeResize(el, () => {
    if (el.isConnected) place();
  });
  return el;
}

/**
 * Where an `over` list goes, from its anchor's top: the current option (else the first) scrolled to the middle of its
 * list, then centered on the anchor. Read once, on opening (the layout of a list just shown).
 */
function measureOver(el: HTMLElement, anchor: DOMRect): number {
  const option =
    el.querySelector<HTMLElement>('[aria-selected="true"]') ?? el.querySelector<HTMLElement>('[role=option]');
  if (!option) return anchor.height + 2;
  // Its scrolling viewport (a long list): the option in its middle, as far as the list scrolls.
  let viewport: HTMLElement | undefined;
  for (let p = option.parentElement; p && p !== el; p = p.parentElement) {
    if (p.scrollHeight > p.clientHeight + 1) {
      viewport = p;
      const vr = p.getBoundingClientRect();
      const or = option.getBoundingClientRect();
      p.scrollTop += or.top - vr.top - (vr.height - or.height) / 2;
      break;
    }
  }
  const o = option.getBoundingClientRect();
  const r = el.getBoundingClientRect();
  let offset = anchor.height / 2 - (o.top - r.top + o.height / 2);
  // Past the window's edge: the list scrolls instead, as far as it can, so the option stays on the field (`place`
  // keeps it on screen, and what the scroll could not take moves it off the field).
  if (viewport) {
    const top = anchor.top + offset;
    const excess = top < 4 ? top - 4 : Math.max(0, top + r.height - (window.innerHeight - 8));
    if (excess !== 0) {
      const before = viewport.scrollTop;
      viewport.scrollTop = before - excess;
      offset += viewport.scrollTop - before;
    }
  }
  return offset;
}
