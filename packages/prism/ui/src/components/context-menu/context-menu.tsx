/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isInside, ownsKey } from '../../dom/events.js';
import { menuKeys, rememberFocus } from '../../dom/focus.js';
import { createRoot } from '../../dom/scope.js';
import { Icon, type IconName } from '../icon/icon.js';
import { classes } from './_generated/context-menu.module.css.js';

/** One entry: what it is called, what it does; an icon (shown as the canvas menus show theirs) is optional. */
export interface ContextMenuItem {
  id: string;
  label: string;
  icon?: IconName;
  /** Why it cannot be used now (greyed). */
  disabled?: string;
  run(): void;
}

/**
 * Opens a menu at a point of the viewport (client px), in `root` (the app's root, so its tokens apply): below and to the
 * right of the point, flipped to the other side where it would leave the window (so near the bottom edge it opens upward).
 * It measures once, on opening. Closes on a choice, on Escape, on a press elsewhere and when the window resizes or
 * scrolls; `close` ends it early. Returns that close. The arrows move between its items; closed, the focus goes back to
 * where it was (accessibility.md).
 */
export function openContextMenu(
  root: HTMLElement,
  at: { x: number; y: number },
  items: readonly ContextMenuItem[],
  onClose?: () => void,
): () => void {
  let closed = false;
  let dispose = () => {};
  const restore = rememberFocus();
  const close = () => {
    if (closed) return;
    closed = true;
    document.removeEventListener('pointerdown', outside, true);
    document.removeEventListener('keydown', onKey, true);
    window.removeEventListener('resize', close);
    window.removeEventListener('scroll', close, true);
    dispose();
    restore(menu);
    menu.remove();
    onClose?.();
  };
  const outside = (e: PointerEvent) => {
    if (!isInside(e, menu)) close();
  };
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape') return;
    // Escape in another viewer on the page is not this menu's.
    if (!isInside(e, menu) && !ownsKey(root, e)) return;
    e.preventDefault();
    e.stopPropagation();
    close();
  };
  let menu!: HTMLElement;
  createRoot((d) => {
    dispose = d;
    menu = (
      <div class={classes.menu} role="menu" onKeyDown={menuKeys} onContextMenu={(e: Event) => e.preventDefault()}>
        {items.map((item) => (
          <button
            type="button"
            role="menuitem"
            class={classes.row}
            disabled={!!item.disabled}
            title={item.disabled}
            onClick={() => {
              close();
              item.run();
            }}
          >
            {item.icon ? <Icon name={item.icon} size={15} /> : <span />}
            <span>{item.label}</span>
          </button>
        ))}
      </div>
    ) as HTMLElement;
  });
  root.append(menu);
  // Measured once, now: the flip at the window's edges.
  const w = menu.offsetWidth;
  const h = menu.offsetHeight;
  const x = at.x + w > window.innerWidth - 4 ? Math.max(4, at.x - w) : at.x;
  const y = at.y + h > window.innerHeight - 4 ? Math.max(4, at.y - h) : at.y;
  menu.style.transform = `translate(${Math.round(x)}px, ${Math.round(y)}px)`;
  document.addEventListener('pointerdown', outside, true);
  document.addEventListener('keydown', onKey, true);
  window.addEventListener('resize', close);
  window.addEventListener('scroll', close, true);
  (menu.querySelector('button:not(:disabled)') as HTMLElement | null)?.focus({ preventScroll: true });
  return close;
}
