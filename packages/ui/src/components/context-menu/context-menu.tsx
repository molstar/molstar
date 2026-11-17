/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRoot, isInside, menuKeys, on, rememberFocus } from '../../dom/index.js';
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
}

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
 * icon, the label, the shortcut.
 */
export function openContextMenu(
  root: HTMLElement,
  at: { x: number; y: number },
  items: readonly ContextMenuItem[],
  onClose?: () => void,
): () => void {
  return openMenu(root, at, items, onClose);
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

function openMenu(
  root: HTMLElement,
  at: MenuPlace,
  items: readonly ContextMenuItem[],
  onClose?: () => void,
): () => void {
  let closed = false;
  let dispose = () => {};
  const restore = rememberFocus();
  const close = () => {
    if (closed) return;
    closed = true;
    // Its listeners go with its scope.
    dispose();
    restore(menu);
    menu.remove();
    onClose?.();
  };
  const icons = items.some((i) => i.icon);
  const checks = items.some((i) => i.checked !== undefined);
  let menu!: HTMLElement;
  createRoot((d) => {
    dispose = d;
    menu = (
      <div
        class={classes.menu}
        role="menu"
        data-icons={icons ? '' : undefined}
        onKeyDown={menuKeys}
        onContextMenu={(e: Event) => e.preventDefault()}
      >
        {items.map((item) => (
          <button
            type="button"
            role={item.checked === undefined ? 'menuitem' : 'menuitemradio'}
            aria-checked={item.checked === undefined ? undefined : String(item.checked)}
            tabIndex={-1}
            class={classes.row}
            aria-disabled={item.disabled ? 'true' : undefined}
            aria-description={item.disabled || undefined}
            data-tip={item.disabled}
            onClick={() => {
              if (item.disabled) return;
              close();
              item.run();
            }}
          >
            {icons ? <span class={classes.slot}>{item.icon ? <Icon name={item.icon} /> : null}</span> : null}
            <span class={classes.label}>{item.label}</span>
            {item.shortcut ? (
              <kbd class={classes.keys}>{item.shortcut}</kbd>
            ) : checks ? (
              <span class={classes.check}>{item.checked ? <Icon name="check" /> : null}</span>
            ) : (
              <span />
            )}
          </button>
        ))}
      </div>
    ) as HTMLElement;
    // In the menu's own scope (dom §6): removed when it closes. Escape in another UI root on the page is not its; Tab
    // closes it, the focus back where it was (and on from there).
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
  (
    menu.querySelector<HTMLElement>('[aria-checked="true"]:not([aria-disabled="true"])') ??
    menu.querySelector<HTMLElement>('[role^="menuitem"]:not([aria-disabled="true"])')
  )?.focus({ preventScroll: true });
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
