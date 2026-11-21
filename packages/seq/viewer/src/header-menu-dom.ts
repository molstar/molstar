/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ContextMenuItem, openContextMenu } from '@molstar/ui/components/context-menu/context-menu';
import type { HeaderMenuItem, HeaderMenuPortFactory } from './header-entries.js';

// The header menu's default adapter (_spec/README.md "Ports", Header menu): `@molstar/ui`'s context menu in the canvas's
// element, at the place in the canvas. Its keyboard is the menu's (the first item focused, the arrows, Enter, Esc, a
// press outside), and closed it gives the focus back to the canvas. An item with a submenu opens it as a flyout beside
// its row (the menu staying open, the row marked): by a click, →, Enter, Space or the pointer resting on it
// (`submenuRestMs`); ← or Esc in the flyout closes it alone, a choice closes all. A scheme previews while its row is
// focused or hovered, and a menu closed without a choice (a flyout too) ends the preview in nothing. Its events are its
// own: the view under it never gets them (`overlay-events.ts`).

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
        const itemOf = (it: HeaderMenuItem, i: number): ContextMenuItem => ({
          id: `${it.name}:${i}`,
          label: it.items ? `${it.label} ▸` : it.kind === 'state' && it.text ? `${it.label} (${it.text})` : it.label,
          ...(it.icon ? { icon: it.icon } : {}),
          ...(it.enabled ? {} : { disabled: it.why ?? 'Not available' }),
          ...(it.key ? { shortcut: it.key } : {}),
          ...(it.checked !== undefined ? { checked: it.checked } : {}),
          ...(it.items ? { items: it.items.map(itemOf) } : {}),
          // a row focused or hovered previews its choice (a scheme), when its command has a preview
          onActive: () => {
            if (!it.items && it.enabled) behavior.preview(it);
          },
          run: () => void behavior.runItem(it),
        });
        const items = menu.groups.flat().map(itemOf);
        let done = false;
        close = openContextMenu(
          container,
          at,
          items,
          () => {
            if (done) return;
            done = true;
            close = undefined;
            behavior.closeMenu();
          },
          // a flyout closed alone ends its preview; the rest of the menu stays
          { restMs: behavior.submenuRestMs, onSubmenuClose: () => behavior.cancelPreview() },
        );
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
