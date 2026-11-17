/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { on, requireScope } from '../../dom/index.js';

// A theme switch around an element (what a view that reads computed colors, a canvas's or a stage's ink, reads again):
// a `data-msui-theme` set on it or an ancestor (styles §1), or the system's light or dark (`prefers-color-scheme`,
// what `auto` and `light-dark()` follow). One attribute filter and one media query: nothing runs between switches.

/**
 * Calls `fn` when the theme around `el` may have changed, in the current scope: while it is active (a switch while
 * paused is told as it resumes), ended with it.
 */
export function onThemeChange(el: Element, fn: () => void): void {
  const scope = requireScope('onThemeChange');
  const doc = el.ownerDocument;
  const view = doc.defaultView;
  const scheme = typeof view?.matchMedia === 'function' ? view.matchMedia('(prefers-color-scheme: dark)') : undefined;
  if (scheme && typeof scheme.addEventListener === 'function') on(scheme, 'change', fn);
  if (typeof MutationObserver === 'undefined') return;
  // Only a theme on `el` or around it (its own descendants' themes, a stage's ink, are not its theme).
  const observer = new MutationObserver((records) => {
    for (const r of records) {
      if (r.target.contains(el)) {
        fn();
        return;
      }
    }
  });
  let paused = false;
  scope.onActive((active) => {
    if (!active) {
      paused = true;
      observer.disconnect();
      return;
    }
    observer.observe(doc.documentElement, { attributes: true, attributeFilter: ['data-msui-theme'], subtree: true });
    if (paused) {
      paused = false;
      fn();
    }
  });
  scope.onCleanup(() => observer.disconnect());
}
