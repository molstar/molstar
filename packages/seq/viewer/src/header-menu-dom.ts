/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ContextMenuItem, openContextMenu } from '@molstar/ui/components/context-menu/context-menu';
import type { HeaderMenuItem, HeaderMenuPortFactory } from './header-entries.js';

// The header menu's default adapter (_spec/README.md "Ports", Header menu): `@molstar/ui`'s context menu in the canvas's
// element, at the place in the canvas. Its keyboard is the menu's (the first item focused, the arrows, Enter, Esc, a
// press outside), and closed it gives the focus back to the canvas. An item with a submenu opens it in its place (the
// context menu has no nested menus); a scheme previews while its row is focused or hovered, and a menu closed without a
// choice ends the preview in nothing.

/** @internal `@molstar/ui`'s context menu over the canvas: one at a time (the lazy chunk `headerMenuOverlay()` loads). */
export function contextMenuPort(): HeaderMenuPortFactory {
  return (container) => {
    if (!container) return { show() {}, hide() {} };
    let close: (() => void) | undefined;
    return {
      show(menu, place, behavior) {
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
          run: () => {
            if (it.items) behavior.openMenu(menu.target, place, it.name);
            else void behavior.runItem(it);
          },
        }));
        const before = container.lastElementChild;
        let done = false;
        close = openContextMenu(container, at, items, () => {
          if (done) return;
          done = true;
          close = undefined;
          behavior.closeMenu();
        });
        // a row focused or hovered previews its choice (a scheme), when its command has a preview
        const el = container.lastElementChild;
        if (el && el !== before) {
          const rows = [...el.querySelectorAll<HTMLElement>('[role^="menuitem"]')];
          const onRow = (e: Event) => {
            const i = rows.indexOf((e.target as HTMLElement).closest('[role^="menuitem"]') as HTMLElement);
            const it = flat[i];
            if (it && !it.items && it.enabled) behavior.preview(it);
          };
          el.addEventListener('focusin', onRow);
          el.addEventListener('pointerover', onRow);
        }
      },
      hide() {
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
