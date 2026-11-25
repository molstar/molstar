/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  HeaderMenuPort,
  HeaderMenuPortFactory,
  HeaderPartsPort,
  HeaderPartsPortFactory,
} from './header-entries.js';

// The headers' DOM adapters (_spec/README.md "Ports", Header menu and Header parts): what the default spec names is a
// thin port; its DOM (`@molstar/ui`'s context menu, buttons and icons: `header-menu-dom.ts`, `header-parts-dom.ts`) is
// a lazy chunk loaded at the menu's first open and at a host's first parts, as the mmCIF runtime is. Why: the hover
// path's slots and chips are drawn on the canvas; the menu's DOM is needed only once someone asks for it.

/** `@molstar/ui`'s context menu over the canvas, loaded at its first open. */
export function headerMenuOverlay(): HeaderMenuPortFactory {
  return (container) => {
    if (!container) return { show() {}, hide() {} };
    let port: HeaderMenuPort | undefined;
    let loading: Promise<HeaderMenuPort> | undefined;
    // the latest show waits for the chunk; a hide or a dispose before it lands drops it
    let asked = 0;
    let gone = false;
    return {
      show(menu, place, behavior) {
        const n = ++asked;
        if (port) return port.show(menu, place, behavior);
        loading ??= import('./header-menu-dom.js').then((m) => (port = m.contextMenuPort()(container)));
        void loading.then((p) => {
          if (n === asked && !gone) p.show(menu, place, behavior);
        });
      },
      hide() {
        asked++;
        port?.hide();
      },
      loaded: () => loading?.then(() => undefined),
      dispose() {
        gone = true;
        port?.dispose?.();
      },
    };
  };
}

/** `@molstar/ui` buttons over the headers' boxes, loaded at their first update. */
export function headerPartsDom(): HeaderPartsPortFactory {
  return (container, behavior) => {
    let port: HeaderPartsPort | undefined;
    let loading: Promise<HeaderPartsPort> | undefined;
    let last: Parameters<HeaderPartsPort['update']>[0] | undefined;
    let gone = false;
    const load = () =>
      (loading ??= import('./header-parts-dom.js').then((m) => (port = m.partsDomPort()(container, behavior))));
    return {
      update(parts) {
        last = parts;
        if (port) return port.update(parts);
        void load().then((p) => {
          if (!gone && last) p.update(last);
        });
      },
      loaded: () => load().then(() => undefined),
      focus(target) {
        port?.focus?.(target);
      },
      dispose() {
        gone = true;
        port?.dispose?.();
      },
    };
  };
}
