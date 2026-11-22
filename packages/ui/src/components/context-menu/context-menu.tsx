/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRoot, DomHost, isInside, menuKeys, on, rememberFocus } from '../../dom/index.js';
import { Icon, type IconName } from '../icon/icon.js';
import { dismissable } from '../popover/dismiss.js';
import { classes } from './_generated/context-menu.module.css.js';

/** One entry: what it is called, what it does; an icon is optional. */
export interface ContextMenuItem {
  id: string;
  label: string;
  icon?: IconName;
  /** Why it cannot be used now: greyed, not run, the reason in its tooltip and its description (still focusable). */
  disabled?: string;
  /** Its key, as text (`⌘K`), quiet at the row's end. */
  shortcut?: string;
  /**
   * One of a single choice (`menuitemradio`): the chosen one (true) in the accent with a check at the row's end. Set on
   * every item of the choice.
   */
  checked?: boolean;
  run(): void;
  /**
   * A submenu: a flyout beside this row (right of it, left of the menu where there is no room), the menu staying open
   * and the row `aria-expanded`. A click, →, Enter, Space or the pointer resting on the row (`restMs`) opens it; ← or
   * Escape in it closes it alone (the focus back on the row), a choice in it closes the whole menu. `run` is not called
   * for such a row.
   */
  items?: readonly ContextMenuItem[];
  /** The row took the focus or the pointer (a preview of its choice): each time, in the flyouts too. */
  onActive?(): void;
}

/** How a menu behaves beyond its items. */
export interface ContextMenuOptions {
  /**
   * How long the pointer rests on a submenu's row before its flyout opens, and on another row of its level before the
   * flyout closes (the time to travel into it: arriving in the flyout cancels it); default 300 ms. 0 turns resting
   * off: a click or a key opens a flyout, and it stays until another opens or it is closed. The tolerance for a
   * diagonal move into a flyout is this time window, not a geometric triangle: crossing the rows between the row and
   * the flyout inside it keeps the flyout, taking longer closes it.
   */
  readonly restMs?: number;
  /** A flyout closed alone (← or Escape, or the rest on another row), the menu still open. */
  readonly onSubmenuClose?: () => void;
}

/** The default rest, ms. */
const DefaultRestMs = 300;

/** Where a menu goes: its corner (x, y), and where it goes instead when it would leave the window. */
interface MenuPlace {
  /** The menu's left edge, or its right edge with `alignEnd`. */
  readonly x: number;
  /** Its top. */
  readonly y: number;
  /** Its right edge at `x` (it opens to the left of `x`). */
  readonly alignEnd?: boolean;
  /** Its bottom here instead of its top at `y` where it would leave the window's bottom (default: `y`). */
  readonly bottom?: number;
}

/**
 * Opens a menu at a point of the viewport (client px), in `root` (a UI root, so its tokens apply): below and to the
 * right of the point, flipped to the other side where it would leave the window (so near the bottom edge it opens
 * upward). It measures once, on opening. Closes on a choice, on Escape or Tab, on a press elsewhere (dismiss.ts) and when
 * the window resizes or scrolls; `close` ends it early. Returns that close. The arrows move between its items (no Tab
 * stops; a disabled one is reached, its reason read); closed, the focus goes back to where it was. Its rows are list rows (styles §4, §5): an icon slot, kept empty in a menu where some item has an
 * icon, the label, the shortcut. An item with `items` opens a flyout beside its row, the menu staying open.
 */
export function openContextMenu(
  root: HTMLElement,
  at: { x: number; y: number },
  items: readonly ContextMenuItem[],
  onClose?: () => void,
  options?: ContextMenuOptions,
): () => void {
  return openMenu(root, at, items, onClose, options);
}

/**
 * Opens a menu beside a button (a rail's, at a window's edge): on its `side`, its top at the button's top, or its bottom
 * at the button's bottom where it would leave the window (a button low on the rail opens upward); a unit clear of the
 * button. Otherwise as `openContextMenu`.
 */
export function openMenuBeside(
  root: HTMLElement,
  anchor: Element,
  items: readonly ContextMenuItem[],
  side: 'left' | 'right' = 'left',
  onClose?: () => void,
): () => void {
  const r = anchor.getBoundingClientRect();
  // A unit (styles §3) clear of the button.
  const gap = 4;
  return openMenu(
    root,
    side === 'left'
      ? { x: r.left - gap, y: r.top, alignEnd: true, bottom: r.bottom }
      : { x: r.right + gap, y: r.top, bottom: r.bottom },
    items,
    onClose,
  );
}

/** An open flyout: its element, the row it opened from, the menu level that row is in. */
interface Flyout {
  readonly el: HTMLElement;
  readonly row: HTMLElement;
  readonly dispose: () => void;
}

/** Ids for the flyouts and their rows. */
let flyoutIds = 0;

const firstRow = (el: HTMLElement) =>
  el.querySelector<HTMLElement>('[role^="menuitem"]:not([aria-disabled="true"])') ??
  el.querySelector<HTMLElement>('[role^="menuitem"]');

function openMenu(
  root: HTMLElement,
  at: MenuPlace,
  items: readonly ContextMenuItem[],
  onClose?: () => void,
  options: ContextMenuOptions = {},
): () => void {
  let closed = false;
  let dispose = () => {};
  const restMs = options.restMs ?? DefaultRestMs;
  const restore = rememberFocus();
  // The open flyouts: `flyouts[d]` hangs from a row of level `d` (the menu is level 0).
  const flyouts: Flyout[] = [];
  let timer: unknown;
  let timerRow: HTMLElement | undefined;
  const stopTimer = () => {
    if (timer !== undefined) DomHost.clearTimeout(timer);
    timer = timerRow = undefined;
  };
  /** Closes the flyouts of levels `from` and deeper, the innermost first. */
  const closeFlyouts = (from: number, tell = true) => {
    const was = flyouts.length;
    while (flyouts.length > from) {
      const f = flyouts.pop()!;
      f.row.setAttribute('aria-expanded', 'false');
      f.dispose();
      f.el.remove();
    }
    if (tell && was > flyouts.length) options.onSubmenuClose?.();
  };
  const close = () => {
    if (closed) return;
    closed = true;
    stopTimer();
    // the focus in a flyout (a choice there, Tab) is lost with it: put back all the same
    const active = menu.ownerDocument.activeElement;
    const hadFocus = !!active && (menu.contains(active) || flyouts.some((f) => f.el.contains(active)));
    closeFlyouts(0, false);
    // Its listeners go with its scope.
    dispose();
    restore(hadFocus ? menu.ownerDocument.body : menu);
    menu.remove();
    onClose?.();
  };

  /** One level's element for `list`; `depth` 0 is the menu, 1 its flyout, and so on. */
  const level = (list: readonly ContextMenuItem[], depth: number): HTMLElement => {
    const icons = list.some((i) => i.icon);
    const checks = list.some((i) => i.checked !== undefined);
    const rowOf = new Map<Element, ContextMenuItem>();
    const choose = (item: ContextMenuItem, row: HTMLElement, fromKeys: boolean) => {
      if (item.disabled) return;
      if (item.items?.length) return openFlyout(item, row, depth, fromKeys);
      close();
      item.run();
    };
    const el = (
      <div
        class={classes.menu}
        role="menu"
        data-icons={icons ? '' : undefined}
        onKeyDown={(e: KeyboardEvent) => {
          if (e.defaultPrevented) return;
          const row = (e.target as Element).closest?.('[role^="menuitem"]') ?? null;
          const item = row && rowOf.get(row);
          if (e.key === 'ArrowRight' && item?.items?.length && !item.disabled) {
            e.preventDefault();
            openFlyout(item, row as HTMLElement, depth, true);
          } else if (e.key === 'ArrowLeft' && depth > 0) {
            // back to the row this flyout came from: the flyout alone closed
            e.preventDefault();
            const from = flyouts[depth - 1]?.row;
            closeFlyouts(depth - 1);
            from?.focus({ preventScroll: true });
          } else menuKeys(e);
        }}
        onContextMenu={(e: Event) => e.preventDefault()}
        onPointerEnter={stopTimer}
        // off the level (into its flyout, whose own enter takes over, or away): no rest runs after
        onPointerLeave={stopTimer}
      >
        {list.map((item) => {
          let row!: HTMLElement;
          const rest = () => {
            if (restMs <= 0) return;
            const open = flyouts[depth];
            if (open?.row === row) return stopTimer();
            if (timerRow === row) return;
            stopTimer();
            const sub = !item.disabled && item.items?.length ? item : undefined;
            // a row with a submenu opens it on rest; another row closes this level's flyout on rest (a diagonal
            // move into the flyout cancels it)
            if (!sub && !open) return;
            timerRow = row;
            timer = DomHost.setTimeout(() => {
              timer = timerRow = undefined;
              if (sub) openFlyout(sub, row, depth, false);
              else closeFlyouts(depth);
            }, restMs);
          };
          row = (
            <button
              type="button"
              role={item.checked === undefined ? 'menuitem' : 'menuitemradio'}
              aria-checked={item.checked === undefined ? undefined : String(item.checked)}
              aria-haspopup={item.items?.length ? 'menu' : undefined}
              aria-expanded={item.items?.length ? 'false' : undefined}
              tabIndex={-1}
              class={classes.row}
              aria-disabled={item.disabled ? 'true' : undefined}
              aria-description={item.disabled || undefined}
              data-tip={item.disabled}
              onClick={(e: MouseEvent) => choose(item, row, e.detail === 0)}
              onFocus={() => item.onActive?.()}
              onPointerOver={() => {
                if (!item.disabled) item.onActive?.();
                rest();
              }}
              // a menu that opened under a still pointer gets no pointerover for its row: the first move is the rest's start
              onPointerMove={rest}
            >
              {icons ? <span class={classes.slot}>{item.icon ? <Icon name={item.icon} /> : null}</span> : null}
              <span class={classes.label}>{item.label}</span>
              {item.shortcut ? (
                <kbd class={classes.keys}>{item.shortcut}</kbd>
              ) : item.items?.length ? (
                // a submenu's arrow at the row's end, as menus have it (the label never carries it)
                <span class={classes.check} data-submenu-arrow="">
                  <Icon name="chevron-right" />
                </span>
              ) : checks ? (
                <span class={classes.check}>{item.checked ? <Icon name="check" /> : null}</span>
              ) : (
                <span />
              )}
            </button>
          ) as HTMLElement;
          rowOf.set(row, item);
          return row;
        })}
      </div>
    ) as HTMLElement;
    return el;
  };

  const openFlyout = (item: ContextMenuItem, row: HTMLElement, depth: number, focus: boolean) => {
    stopTimer();
    const open = flyouts[depth];
    if (open?.row === row) {
      if (focus) firstRow(open.el)?.focus({ preventScroll: true });
      return;
    }
    closeFlyouts(depth);
    const parent = row.closest<HTMLElement>('[role="menu"]')!;
    let stop = () => {};
    let el!: HTMLElement;
    createRoot((d) => {
      stop = d;
      el = level(item.items!, depth + 1);
      // A layer of its own inside the menu's (dismiss.ts): Escape closes this flyout back to its row (anywhere in the
      // menu: the innermost layer takes it), a press elsewhere is the menu's own, a press on the row its click's.
      dismissable({
        layer: () => el,
        anchor: row,
        outside: false,
        // Tab in a flyout closes the whole menu (its rows are no Tab stops; the focus goes back where it was)
        tab: true,
        close: (reason) => (reason === 'tab' ? close() : closeFlyouts(depth)),
        refocus: () => row,
        ownsKey: (e) => isInside(e, parent),
      });
    });
    // a long list scrolls inside the window
    el.style.maxHeight = `${window.innerHeight - 8}px`;
    el.style.overflowY = 'auto';
    root.append(el);
    // the row names what it opens, the flyout what it came from (ARIA's menu button pattern)
    row.id ||= `msui-cm-row-${++flyoutIds}`;
    el.id = `msui-cm-sub-${flyoutIds}`;
    el.setAttribute('aria-labelledby', row.id);
    row.setAttribute('aria-controls', el.id);
    row.setAttribute('aria-expanded', 'true');
    flyouts[depth] = { el, row, dispose: stop };
    // Beside the row, the menu staying: right of the parent, its top at the row's (less the menu's padding); left of
    // the parent where there is no room; kept inside the window.
    const r = row.getBoundingClientRect();
    const p = parent.getBoundingClientRect();
    const w = el.offsetWidth;
    const h = el.offsetHeight;
    const gap = 2;
    const right = p.right + gap;
    const x = right + w > window.innerWidth - 4 ? Math.max(4, p.left - gap - w) : right;
    const top = r.top - 6;
    const y = top + h > window.innerHeight - 4 ? Math.max(4, Math.min(top, r.bottom + 6 - h)) : Math.max(4, top);
    el.style.transform = `translate(${Math.round(x)}px, ${Math.round(y)}px)`;
    if (focus) firstRow(el)?.focus({ preventScroll: true });
  };

  let menu!: HTMLElement;
  createRoot((d) => {
    dispose = d;
    menu = level(items, 0);
    // In the menu's own scope (dom §6): removed when it closes. Escape in another UI root on the page is not its; Tab
    // closes it, the focus back where it was (and on from there). Its flyouts are layers opened after it: a press in
    // one is not outside.
    dismissable({
      layer: () => menu,
      tab: true,
      close,
      ownsKey: (e) => isInside(e, menu) || ownsKey(root, e),
    });
    on(window, 'resize', close);
    on(window, 'scroll', close, { capture: true });
  });
  root.append(menu);
  // Measured once, now: the flip at the window's edges.
  const w = menu.offsetWidth;
  const h = menu.offsetHeight;
  const left = at.alignEnd ? at.x - w : at.x;
  const x = at.alignEnd
    ? left < 4
      ? Math.min(window.innerWidth - w - 4, at.x)
      : left
    : left + w > window.innerWidth - 4
      ? Math.max(4, at.x - w)
      : left;
  const y = at.y + h > window.innerHeight - 4 ? Math.max(4, (at.bottom ?? at.y) - h) : at.y;
  menu.style.transform = `translate(${Math.round(x)}px, ${Math.round(y)}px)`;
  // The chosen one of a choice first, else the first usable.
  (menu.querySelector<HTMLElement>('[aria-checked="true"]:not([aria-disabled="true"])') ?? firstRow(menu))?.focus({
    preventScroll: true,
  });
  return close;
}

/**
 * Whether a key pressed outside the menu is its root's: the focus is in the root's UI (its shadow host, or the UI root
 * it is in), or nothing on the page has the focus and the pointer is over it. Outside any UI root, every key is.
 */
function ownsKey(root: HTMLElement, e: KeyboardEvent): boolean {
  const node = root.getRootNode();
  const host =
    typeof ShadowRoot !== 'undefined' && node instanceof ShadowRoot
      ? node.host
      : (root.closest('.msui-root') ?? undefined);
  if (!host) return true;
  if (isInside(e, host)) return true;
  const doc = host.ownerDocument;
  const focused = doc.activeElement;
  if (focused && focused !== doc.body && focused !== doc.documentElement) return false;
  return host.matches(':hover');
}
