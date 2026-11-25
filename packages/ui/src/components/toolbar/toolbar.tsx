/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, isState, KeyedSelector, type ReadonlyState } from '@molstar/core/reactive';
import {
  deepActiveElement,
  focusFirstItem,
  For,
  menuKeys,
  on,
  onActive,
  onCleanup,
  type Restartable,
  restartable,
  rovingFocus,
  Show,
  setAttr,
  type WatchSource,
  watch,
} from '../../dom/index.js';
import { held, state } from '../../state.js';
import { Icon, type IconName } from '../icon/icon.js';
import { dismissable } from '../popover/dismiss.js';
import { Floating } from '../select/select.js';
import { classes } from './_generated/toolbar.module.css.js';

// A toolbar (styles §4.1, Toolbars): groups of square icon buttons, each group a block of its own (a surface), the gaps
// between blocks the grouping. WAI-ARIA's toolbar pattern (one Tab stop, the arrows move along it) and its menu button
// pattern (a button's menu opens from the keyboard at its first item, on hover with a mouse, on a tap on touch).

/** What a button looks like now, for a button whose icon and title change (a theme switch, a mode). */
export interface ToolbarFace {
  icon: IconName;
  title: string;
}

/** A row of a button's menu. */
export interface ToolbarMenuRow {
  id: string;
  label: string;
  /** A 12px mark in the row's icon slot (styles §4.1: a menu's icons mark commands); the slot is kept without one. */
  icon?: IconName;
  /** Quiet text at the row's end (a key, a value). */
  hint?: string;
  /**
   * The rows' grouping: consecutive rows with the same `group` are a group; a hairline separates groups, or the
   * group's heading (`heading` on its first row) does.
   */
  group?: string;
  /** A heading over this row's group (read on the group's first row): xs muted text on the menu's text column. */
  heading?: string;
  /**
   * A toggle (`menuitemcheckbox`) or, with `choice`, one of a single choice (`menuitemradio`, the chosen one also in the
   * accent): both show a check mark while true, so the state is never told by colour alone.
   */
  pressed?: ReadonlyState<boolean>;
  choice?: boolean;
  /**
   * Not usable now: greyed, not run, but reached by the arrows (APG menus); a string is the reason, in its tooltip and
   * its accessible description.
   */
  disabled?: boolean | string;
  run?(): void;
  /** A submenu: its rows open beside the menu (hover, a click, →; ← or Esc closes it). */
  items?: readonly ToolbarMenuRow[];
}

/** A toolbar button. */
export interface ToolbarItem {
  id: string;
  icon: IconName;
  /** Its accessible name and tooltip. */
  title: string;
  /** Its icon and title when they change; `icon` and `title` are then what it starts with only in the types. */
  face?: ReadonlyState<ToolbarFace>;
  /** A short text under the icon (`captions: 'below'`) or beside it (`'beside'`): a key, a word. */
  caption?: string;
  /** A toggle: `aria-pressed`, shown pressed while true. */
  pressed?: ReadonlyState<boolean>;
  /**
   * Not usable now (it keeps its place in the arrows' order, `aria-disabled`); a string is the reason, in its tooltip
   * and its accessible description (`aria-description`).
   */
  disabled?: boolean | string | ReadonlyState<boolean | string>;
  /** Whether it shows (default true); a block whose buttons are all hidden goes too. */
  visible?: boolean | ReadonlyState<boolean>;
  /** What a click does; `anchor` is the button (to place what it opens). Without it a button with a menu opens it. */
  run?(anchor: HTMLElement): void;
  /**
   * Its menu: it opens on hover (a mouse), on a tap (touch) and from the keyboard; the button shows a corner mark in
   * place of a tooltip. With `run` too, a click runs and a tap opens the menu with the button's own action first.
   */
  menu?: readonly ToolbarMenuRow[] | ReadonlyState<readonly ToolbarMenuRow[]>;
}

/** A block of buttons. */
export interface ToolbarGroup {
  id: string;
  /** The block's accessible name (`role="group"`). */
  label?: string;
  items: readonly ToolbarItem[];
}

/** Where menus open: above or below a horizontal bar, left or right of a vertical one (flipped at the window edge). */
export type ToolbarMenuSide = 'above' | 'below' | 'left' | 'right';

export interface ToolbarProps {
  groups: readonly ToolbarGroup[] | ReadonlyState<readonly ToolbarGroup[]>;
  /** `horizontal` (default): ← → move, ↓ ↑ open a menu. `vertical`: ↓ ↑ move, → ← open a menu. */
  orientation?: 'horizontal' | 'vertical';
  /** Default: `below` a horizontal bar, `right` of a vertical one. */
  menuSide?: ToolbarMenuSide;
  /** Where a caption goes: `below` the icon (default in a vertical bar) or `beside` it (in a horizontal one). */
  captions?: 'below' | 'beside';
  /**
   * How groups are told apart: `blocks` (default), each group a surface of its own, the gaps between them the grouping;
   * `separators`, one surface for the whole bar, a hairline between groups.
   */
  grouping?: 'blocks' | 'separators';
  /** Blocks that float over content (a canvas): raised, with a hairline and a shadow. */
  elevated?: boolean;
  'aria-label'?: string;
  class?: string;
}

/**
 * A button's menu led by its own action (a tap's menu on touch; a folded toolbar's submenu, stage §3.2): a row running
 * `run` under the button's face as it is now (id `${item.id}:run`), unless the item has no action or a row of the menu
 * already is it (by id: the item's, or its own row's). A row that only shares its label stays (a label is not an
 * identity: two actions may read the same).
 */
export function withOwnAction(
  item: ToolbarItem,
  rows: readonly ToolbarMenuRow[],
  anchor: () => HTMLElement,
): readonly ToolbarMenuRow[] {
  if (!item.run) return rows;
  const id = `${item.id}:run`;
  for (const r of rows) if (r.id === id || r.id === item.id) return rows;
  const face = item.face?.value;
  const own: ToolbarMenuRow = {
    id,
    label: face?.title ?? item.title,
    icon: face?.icon ?? item.icon,
    group: '\u0000own',
    ...(item.pressed ? { pressed: item.pressed } : {}),
    run: () => item.run?.(anchor()),
  };
  return [own, ...rows];
}

/** How long a menu stays after the pointer leaves its button or itself: time to travel between them. */
export const MenuCloseDelayMs = 160;
/** A menu's distance from its button (px): a block's padding (`--msui-surface-pad`) and a unit clear of the block. */
const MenuGap = 8;
/** A menu's own padding (px, `--msui-surface-pad`): the menu shifts by it so its icon column is the button's icon's. */
const MenuPad = 4;

/** The open menu: whose, from where, and how it opened. */
interface OpenMenu {
  readonly item: ToolbarItem;
  readonly button: HTMLElement;
  /** Opened from the keyboard: the focus goes to its first item. */
  readonly focus: boolean;
  /** Opened by a tap: the button's own action is its first row. */
  readonly touch: boolean;
}

/**
 * Which menu is open (one per toolbar) and its close timer. A button follows its own key of the open menu
 * (`KeyedSelector`): opening a menu writes to the button that opened it and the one that closed, not to every button.
 */
class ToolbarMenus {
  readonly open = state<OpenMenu | undefined>(undefined, Object.is);
  readonly shown = new KeyedSelector(derive(this.open, (o) => o?.item));
  /** The open menu's outer element (the menu and its submenus), while open. */
  el: HTMLElement | undefined;
  /** The open menu's first level (`role="menu"`). */
  menu: HTMLElement | undefined;
  private readonly timer: Restartable;

  constructor() {
    // Scoped (dom §6): held while the toolbar is paused, dropped with it.
    this.timer = restartable(() => this.hide());
  }

  show(item: ToolbarItem, button: HTMLElement, focus: boolean, touch = false) {
    this.keep();
    const current = this.open.value;
    if (current?.item === item && current.button === button) {
      // Open already (by hover): from the keyboard, into it.
      if (focus && this.menu) focusFirstItem(this.menu);
      return;
    }
    this.open.set({ item, button, focus, touch });
  }

  keep() {
    this.timer.cancel();
  }

  hideSoon() {
    this.timer.start(MenuCloseDelayMs);
  }

  hide() {
    this.keep();
    if (this.open.value) this.open.set(undefined);
  }

  /** Whether the focus is in the open menu (opened from the keyboard: the pointer leaving does not close it). */
  hasFocus(): boolean {
    const active = deepActiveElement();
    return !!active && !!this.el?.contains(active);
  }
}

/**
 * A toolbar: `groups` of buttons, each a block (a surface), in a row or a column. One Tab stop: the arrows along it move
 * between its buttons (roving `tabindex`, across blocks). A button with a menu opens it on hover (a mouse; it closes
 * `MenuCloseDelayMs` after the pointer leaves the button and the menu), on a tap (touch), and from the keyboard (the
 * arrows across the bar, Enter or Space on a button without its own action: at the first item). Esc closes it and the
 * focus goes back to the button. Tooltips (`data-tip`) need `installTooltips` on a root around it.
 */
export function Toolbar(props: ToolbarProps) {
  const orientation = props.orientation ?? 'horizontal';
  const side = props.menuSide ?? (orientation === 'horizontal' ? 'below' : 'right');
  const captions = props.captions ?? (orientation === 'horizontal' ? 'beside' : 'below');
  const menus = new ToolbarMenus();
  const ctx: ButtonContext = {
    menus,
    orientation,
    tipSide: side === 'above' ? 'top' : side === 'left' ? 'left' : 'bottom',
    items: new WeakMap(),
    pointer: 'mouse',
  };
  const root = (
    <div
      class={props.class ? `${classes.toolbar} ${props.class}` : classes.toolbar}
      role="toolbar"
      aria-orientation={orientation}
      aria-label={props['aria-label']}
      data-orientation={orientation}
      data-captions={captions}
      data-grouping={props.grouping ?? 'blocks'}
      data-elevated={props.elevated ? '' : undefined}
      onPointerDown={(e: PointerEvent) => {
        ctx.pointer = e.pointerType || 'mouse';
      }}
    >
      <For each={props.groups} by={(g) => g.id}>
        {(group) => (
          <div class={classes.group} role="group" aria-label={group.label}>
            {group.items.map((item) => (
              <ToolbarButton item={item} ctx={ctx} />
            ))}
          </div>
        )}
      </For>
      <Show when={menus.open}>{(open) => <ToolbarMenu open={open} menus={menus} side={side} />}</Show>
    </div>
  ) as HTMLElement;
  // One Tab stop (dom §7, `rovingFocus`): the focused button, else the last focused that still shows, else the first;
  // the arrows along the bar move between buttons (wrapping; Home, End; ← → swapped right to left), across its blocks.
  // A disabled button keeps its place (`aria-disabled`: its reason is read).
  rovingFocus(root, { items: '[data-toolbar-button]', axis: orientation });
  // Hover: two listeners for the whole bar (pointer over and out bubble), not two per button.
  on(root, 'pointerover', (e: PointerEvent) => {
    if (e.pointerType !== 'mouse') return;
    const target = e.target as Element;
    if (menus.el?.contains(target)) return menus.keep();
    const button = target.closest?.('[data-toolbar-button]') as HTMLElement | null;
    const item = button && ctx.items.get(button);
    if (!item?.menu || !button || isDisabled(item)) return;
    menus.show(item, button, false);
  });
  on(root, 'pointerout', (e: PointerEvent) => {
    const open = menus.open.value;
    if (e.pointerType !== 'mouse' || !open) return;
    const to = e.relatedTarget as Node | null;
    if (to && (menus.el?.contains(to) || open.button.contains(to))) return;
    if (!menus.hasFocus()) menus.hideSoon();
  });
  // Paused (hidden): no menu stays open over it.
  onActive((active) => {
    if (!active) menus.hide();
  });
  return root;
}

/** What the buttons of one toolbar share. */
interface ButtonContext {
  readonly menus: ToolbarMenus;
  readonly orientation: 'horizontal' | 'vertical';
  readonly tipSide: 'top' | 'bottom' | 'left';
  /** Each button's item (the hover listener's lookup). */
  readonly items: WeakMap<Element, ToolbarItem>;
  /** The last press's pointer type (`touch`: a tap opens a menu rather than running the button). */
  pointer: string;
}

function isDisabled(item: ToolbarItem): boolean {
  const d = item.disabled;
  return !!(isState(d) ? (d as ReadonlyState<boolean | string>).value : d);
}

function ToolbarButton(props: { item: ToolbarItem; ctx: ButtonContext }) {
  const { item, ctx } = props;
  const { menus } = ctx;
  const hasMenu = item.menu !== undefined;
  const face = item.face;
  const pressed = item.pressed;
  const disabled = isState(item.disabled) ? (item.disabled as ReadonlyState<boolean | string>) : undefined;
  const visible = isState(item.visible) ? (item.visible as ReadonlyState<boolean>) : undefined;
  const expanded = hasMenu ? menus.shown.key(item) : undefined;
  // The icon: fixed, or the face's (replaced when its name changes).
  const slot = (<span class={classes.face}>{face ? null : <Icon name={item.icon} />}</span>) as HTMLElement;
  let shownIcon: IconName | undefined;
  const button = (
    <button
      type="button"
      class={classes.button}
      tabIndex={-1}
      data-toolbar-button=""
      data-captioned={item.caption ? '' : undefined}
      data-tip-side={ctx.tipSide}
      aria-haspopup={hasMenu ? 'menu' : undefined}
      onClick={(e: MouseEvent) => {
        const tap = ctx.pointer !== 'mouse';
        ctx.pointer = 'mouse';
        if (isDisabled(item)) return;
        const fromKeys = e.detail === 0;
        if (hasMenu) {
          // A tap: no hover to show the menu, so the tap opens it (or closes it).
          if (tap && !fromKeys) {
            if (expanded!.value) menus.hide();
            else menus.show(item, button, false, !!item.run);
            return;
          }
          // A button that only holds a menu opens it (from the keyboard, at its first item).
          if (!item.run) return menus.show(item, button, fromKeys);
        }
        menus.hide();
        item.run?.(button);
      }}
      onKeyDown={
        hasMenu
          ? (e: KeyboardEvent) => {
              // Escape on the button with its menu open is the menu's (dismiss.ts).
              const opens =
                ctx.orientation === 'horizontal'
                  ? e.key === 'ArrowDown' || e.key === 'ArrowUp'
                  : e.key === 'ArrowRight' || e.key === 'ArrowLeft';
              if (!opens || isDisabled(item)) return;
              e.preventDefault();
              menus.show(item, button, true);
            }
          : undefined
      }
      onFocusOut={
        hasMenu
          ? (e: FocusEvent) => {
              const to = e.relatedTarget as Node | null;
              if (expanded!.value && !(to && menus.el?.contains(to))) menus.hideSoon();
            }
          : undefined
      }
    >
      {slot}
      {item.caption ? (
        <span class={classes.caption} aria-hidden="true">
          {item.caption}
        </span>
      ) : null}
      {hasMenu ? <span class={classes.mark} aria-hidden="true" /> : null}
    </button>
  ) as HTMLButtonElement;
  ctx.items.set(button, item);
  // What changes (its face, pressed, open menu, disabled, shown): one watch for all; a button with none has none.
  const render = () => {
    const f = face?.value;
    const title = f ? f.title : item.title;
    if (f && f.icon !== shownIcon) {
      shownIcon = f.icon;
      slot.replaceChildren(Icon({ name: f.icon }));
    }
    const d = disabled ? disabled.value : item.disabled;
    const reason = typeof d === 'string' && d ? d : undefined;
    setAttr(button, 'aria-label', title);
    // A button with a menu shows the menu instead of a tooltip (both would open in the same spot); disabled, its reason.
    setAttr(button, 'data-tip', reason ? `${title}: ${reason}` : hasMenu && !d ? undefined : title);
    setAttr(button, 'aria-disabled', d ? 'true' : undefined);
    // The reason, for assistive technology too (the tooltip is only seen).
    setAttr(button, 'aria-description', reason);
    if (pressed) setAttr(button, 'aria-pressed', String(pressed.value));
    if (expanded) setAttr(button, 'aria-expanded', String(expanded.value));
    if (visible) setAttr(button, 'hidden', !visible.value);
    else if (item.visible === false) setAttr(button, 'hidden', true);
  };
  const sources: WatchSource[] = [];
  if (face) sources.push(face);
  if (pressed) sources.push(pressed);
  if (expanded) sources.push(expanded);
  if (disabled) sources.push(disabled);
  if (visible) sources.push(visible);
  if (sources.length) watch(sources, render);
  else render();
  return button;
}

/** The open menu of a toolbar: its first level by its button, its submenus beside it. */
function ToolbarMenu(props: { open: OpenMenu; menus: ToolbarMenus; side: ToolbarMenuSide }) {
  const { open, menus } = props;
  const { item, button } = open;
  const rows = held(item.menu, [] as readonly ToolbarMenuRow[]);
  // Opened by a tap on a button with its own action: that action first (the tap did not run it).
  const anchor = () => button;
  // Closed while it had the focus (Esc, a choice): the focus goes back to its button, not to the page.
  let hadFocus = false;
  const vertical = props.side === 'left' || props.side === 'right';
  let first: HTMLElement | undefined;
  const el = (
    <div
      class={classes.menuRoot}
      onFocusIn={() => {
        hadFocus = true;
        menus.keep();
      }}
      onFocusOut={(e: FocusEvent) => {
        const to = e.relatedTarget as Node | null;
        if (!to || !(el.contains(to) || button.contains(to))) menus.hideSoon();
      }}
    >
      <Show when={rows}>
        {(list) => (
          <MenuLevel
            rows={open.touch ? withOwnAction(item, list, anchor) : list}
            label={item.face?.value.title ?? item.title}
            menus={menus}
            back={vertical ? () => menus.hide() : undefined}
            ref={(m) => {
              first = m;
              menus.menu = m;
            }}
          />
        )}
      </Show>
    </div>
  ) as HTMLElement;
  menus.el = el;
  const floating = (
    <Floating anchor={button} placement={props.side} gap={MenuGap} shift={vertical ? 0 : MenuPad}>
      {el}
    </Floating>
  );
  // After the floating shows (its own microtask, queued as it was made): the focus into the menu.
  if (open.focus) queueMicrotask(() => first && focusFirstItem(first));
  // A press outside it and its button, Escape (in it, on its button) or Tab closes it, the focus back on the button
  // (dismiss.ts); a submenu takes its own Escape first.
  dismissable({ layer: () => el, anchor: button, tab: true, close: () => menus.hide(), refocus: () => button });
  onCleanup(() => {
    if (menus.el === el) menus.el = undefined;
    if (menus.menu === first) menus.menu = undefined;
    if (!hadFocus) return;
    const now = deepActiveElement();
    if (!now || el.contains(now)) button.focus({ preventScroll: true });
  });
  return floating;
}

/**
 * One level of a menu: its rows (`role="menu"`), and the submenu open beside it, if any. ↓ ↑ Home End move (wrapping),
 * → or a click opens a row's submenu (at its first item from the keyboard), hover opens it, hovering another row closes
 * it; the arrow back toward what it opened from (`back`: ← for a level on the right, → for one flipped to the left)
 * and Esc in a submenu close it, Esc in the first level closes the menu, Tab anywhere closes it all (dismiss.ts).
 */
function MenuLevel(props: {
  rows: readonly ToolbarMenuRow[];
  label: string;
  menus: ToolbarMenus;
  /** The arrow back here: to the row it opened from (a submenu), or to the button (a menu beside a vertical bar). */
  back?: () => void;
  ref?: (menu: HTMLElement) => void;
}) {
  const { menus } = props;
  const sub = state<{ row: ToolbarMenuRow; el: HTMLElement; focus: boolean } | undefined>(undefined, Object.is);
  // Which row's submenu is open, by row: opening one writes to two rows.
  const subOpen = new KeyedSelector(derive(sub, (s) => s?.row));
  const rowOf = new WeakMap<Element, ToolbarMenuRow>();
  // The open submenu's level (its ref): the focus into it from the keyboard.
  let subMenu: HTMLElement | undefined;
  const openSub = (row: ToolbarMenuRow, el: HTMLElement, focus: boolean) => {
    if (sub.value?.row === row) {
      if (focus && subMenu) focusFirstItem(subMenu);
      return;
    }
    sub.set({ row, el, focus });
  };
  const closeSub = () => {
    if (sub.value) sub.set(undefined);
  };
  const run = (row: ToolbarMenuRow, el: HTMLElement, fromKeys: boolean) => {
    if (row.disabled) return;
    if (row.items) return openSub(row, el, fromKeys);
    menus.hide();
    row.run?.();
  };
  // The arrow back: toward the side it opened from, as placed (`Floating` flips a level that would leave the window).
  const backKey = () => {
    const side = menu.closest<HTMLElement>('[popover]')?.dataset.side;
    return side === 'right' ? 'ArrowLeft' : side === 'left' ? 'ArrowRight' : undefined;
  };
  // A trailing column for the check or the chevron, only in a menu where some row has one.
  const ends = props.rows.some((r) => r.items || r.pressed);
  const menu = (
    <div
      class={classes.menu}
      role="menu"
      aria-label={props.label}
      data-ends={ends ? '' : undefined}
      onKeyDown={(e: KeyboardEvent) => {
        if (e.defaultPrevented) return;
        const active = deepActiveElement();
        const row = active && rowOf.get(active);
        if (e.key === 'ArrowRight' && row?.items && !row.disabled) {
          e.preventDefault();
          openSub(row, active!, true);
        } else if (props.back && e.key === backKey()) {
          e.preventDefault();
          props.back();
        } else menuKeys(e);
      }}
      onClick={(e: MouseEvent) => {
        const el = (e.target as Element).closest?.('[data-menu-row]') as HTMLElement | null;
        const row = el && rowOf.get(el);
        if (row && el) run(row, el, e.detail === 0);
      }}
    >
      {props.rows.map((row, i, all) => {
        const starts = i > 0 && row.group !== all[i - 1].group;
        const el = MenuRow({ row, open: row.items ? subOpen.key(row) : undefined });
        rowOf.set(el, row);
        return (
          <>
            {starts && !row.heading ? <div class={classes.divider} role="separator" /> : null}
            {row.heading && (i === 0 || starts) ? (
              <div class={classes.heading} role="presentation" data-first={i === 0 ? '' : undefined}>
                {row.heading}
              </div>
            ) : null}
            {el}
          </>
        );
      })}
    </div>
  ) as HTMLElement;
  props.ref?.(menu);
  // Hover over a row: its submenu opens, another row's closes it.
  on(menu, 'pointerover', (e: PointerEvent) => {
    if (e.pointerType !== 'mouse') return;
    const el = (e.target as Element).closest?.('[data-menu-row]') as HTMLElement | null;
    const row = el && rowOf.get(el);
    if (!row || !el) return;
    if (row.items && !row.disabled) openSub(row, el, false);
    else closeSub();
  });
  return (
    <div class={classes.level}>
      {menu}
      <Show when={sub}>
        {(s) => {
          const back = () => {
            closeSub();
            s.el.focus({ preventScroll: true });
          };
          const floating = (
            <Floating anchor={s.el} placement="right" gap={MenuGap} shift={MenuPad}>
              <MenuLevel
                rows={s.row.items ?? []}
                label={s.row.label}
                menus={menus}
                back={back}
                ref={(m) => {
                  subMenu = m;
                }}
              />
            </Floating>
          ) as HTMLElement;
          // A layer of its own, inside the menu's (dismiss.ts): Esc here closes this level back to its row, Tab the whole
          // menu (the focus on the button); a press outside is the menu's.
          dismissable({
            layer: () => floating,
            anchor: s.el,
            outside: false,
            tab: true,
            close: (reason) => (reason === 'tab' ? menus.hide() : back()),
            refocus: (reason) => (reason === 'tab' ? menus.open.value?.button : undefined),
          });
          if (s.focus) queueMicrotask(() => subMenu && focusFirstItem(subMenu));
          return floating;
        }}
      </Show>
    </div>
  );
}

/** A menu row (styles §5): the icon slot at the inset, the label on the text column, a hint, the check or chevron. */
function MenuRow(props: { row: ToolbarMenuRow; open: ReadonlyState<boolean> | undefined }): HTMLElement {
  const { row } = props;
  const toggle = row.pressed;
  const role = toggle ? (row.choice ? 'menuitemradio' : 'menuitemcheckbox') : 'menuitem';
  // A disabled row's reason: its tooltip, and its description (it takes the focus, so it is read).
  const reason = typeof row.disabled === 'string' && row.disabled ? row.disabled : undefined;
  const el = (
    <button
      type="button"
      role={role}
      tabIndex={-1}
      class={classes.row}
      data-menu-row=""
      aria-haspopup={row.items ? 'menu' : undefined}
      aria-disabled={row.disabled ? 'true' : undefined}
      aria-description={reason}
      data-tip={reason}
    >
      <span class={classes.slot}>{row.icon ? <Icon name={row.icon} /> : null}</span>
      <span class={classes.label}>{row.label}</span>
      <span class={classes.hint}>{row.hint ?? ''}</span>
      {row.items ? (
        <span class={classes.end}>
          <Icon name="chevron-right" class={classes.chevron} />
        </span>
      ) : toggle ? (
        <span class={classes.end}>
          <Icon name="check" class={classes.check} />
        </span>
      ) : null}
    </button>
  ) as HTMLElement;
  const { open } = props;
  if (toggle) watch([toggle], () => setAttr(el, 'aria-checked', String(toggle.value)));
  if (open) watch([open], () => setAttr(el, 'aria-expanded', String(open.value)));
  return el;
}
