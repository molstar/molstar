/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ContextMenuItem, openContextMenu } from '@molstar/ui/components/context-menu/context-menu';
import { DomHost } from '@molstar/ui/dom/host';
import type { HeaderMenu, HeaderMenuItem, HeaderMenuPlace, HeaderMenuPortFactory } from './header-entries.js';

// The header menu's default adapter (_spec/README.md "Ports", Header menu): `@molstar/ui`'s context menu in the canvas's
// element, at the place in the canvas. Its keyboard is the menu's (the first item focused, the arrows, Enter, Esc, a
// press outside), and closed it gives the focus back to the canvas. An item with a submenu opens it in its place (the
// context menu has no nested menus) by a click, by → and by the pointer resting on it (`submenuRestMs`); ← in a
// submenu so opened goes back to its menu, onto its row. A scheme previews while its row is focused or hovered, and a
// menu closed without a choice ends the preview in nothing. Its events are its own: the view under it never gets them (`overlay-events.ts`).

/** @internal `@molstar/ui`'s context menu over the canvas: one at a time (the lazy chunk `headerMenuOverlay()` loads). */
export function contextMenuPort(): HeaderMenuPortFactory {
  return (container) => {
    if (!container) return { show() {}, hide() {} };
    let close: (() => void) | undefined;
    // the menu a submenu opened from its row came from (← goes back to it); the pointer's rest on a submenu's row
    let parent: { menu: HeaderMenu; place: HeaderMenuPlace } | undefined;
    let rest: unknown;
    let restRow: HeaderMenuItem | undefined;
    const stopRest = () => {
      if (rest !== undefined) DomHost.clearTimeout(rest);
      rest = restRow = undefined;
    };
    return {
      show(menu, place, behavior) {
        stopRest();
        const from = parent;
        parent = undefined;
        close?.();
        const r = container.getBoundingClientRect();
        const at = { x: r.left + place.x, y: r.top + place.y };
        const flat: HeaderMenuItem[] = menu.groups.flat();
        const items: ContextMenuItem[] = flat.map((it, i) => ({
          id: `${it.name}:${i}`,
          label: it.items ? `${it.label} ▸` : it.kind === 'state' && it.text ? `${it.label} (${it.text})` : it.label,
          ...(it.icon ? { icon: it.icon } : {}),
          ...(it.enabled ? {} : { disabled: it.why ?? 'Not available' }),
          ...(it.key ? { shortcut: it.key } : {}),
          ...(it.checked !== undefined ? { checked: it.checked } : {}),
          run: () => (it.items ? sub(it) : void behavior.runItem(it)),
        }));
        // in this menu's place, as a click on its row: this menu closed first (the headers' state then follows the new one)
        const sub = (it: HeaderMenuItem) => {
          close?.();
          parent = { menu, place };
          if (!behavior.openMenu(menu.target, place, it.name)) parent = undefined;
        };
        const before = container.lastElementChild;
        let done = false;
        close = openContextMenu(container, at, items, () => {
          if (done) return;
          done = true;
          close = undefined;
          // however it closed (Esc, a press outside, a choice): no rest opens a submenu after it
          stopRest();
          parent = undefined;
          behavior.closeMenu();
        });
        // a row focused or hovered previews its choice (a scheme), when its command has a preview
        const el = container.lastElementChild;
        if (el && el !== before) {
          const rows = [...el.querySelectorAll<HTMLElement>('[role^="menuitem"]')];
          const rowOf = (e: Event) => {
            const i = rows.indexOf((e.target as HTMLElement).closest('[role^="menuitem"]') as HTMLElement);
            return i < 0 ? undefined : flat[i];
          };
          const onRow = (e: Event) => {
            const it = rowOf(e);
            if (it && !it.items && it.enabled) behavior.preview(it);
          };
          el.addEventListener('focusin', onRow);
          el.addEventListener('pointerover', onRow);
          // the pointer resting on a submenu's row opens it; moving on to another row (or off the menu) does not
          // (a move onto a child of the same row keeps its rest)
          el.addEventListener('pointerover', (e) => {
            const it = rowOf(e);
            if (it && it === restRow) return;
            stopRest();
            if (!it?.items || !it.enabled) return;
            restRow = it;
            rest = DomHost.setTimeout(() => {
              rest = restRow = undefined;
              sub(it);
            }, behavior.submenuRestMs);
          });
          el.addEventListener('pointerleave', stopRest);
          // → on a submenu's row opens it, ← in a submenu goes back to its menu
          el.addEventListener('keydown', (e: Event) => {
            const key = (e as KeyboardEvent).key;
            if (e.defaultPrevented) return;
            const it = rowOf(e);
            if (key === 'ArrowRight' && it?.items && it.enabled) {
              e.preventDefault();
              sub(it);
            } else if (key === 'ArrowLeft' && from) {
              e.preventDefault();
              close?.();
              behavior.openMenu(from.menu.target, from.place);
            }
          });
        }
      },
      hide() {
        stopRest();
        parent = undefined;
        const c = close;
        close = undefined;
        c?.();
      },
      dispose() {
        this.hide();
      },
    };
  };
}
