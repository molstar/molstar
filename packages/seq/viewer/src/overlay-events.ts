/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Whose a DOM event in the mounted canvas is (_spec/README.md "Input", Overlays): the view's when it happened on the
// container itself or on one of its canvases; an overlay's (the header menu, the header parts, the keys panel, a host's
// part attached to the container) when it happened in any other element inside it. Why: the overlays live in the
// container (their tokens, their place over the canvas), so their events bubble to its listeners; the view acts only on
// its own, or a click on a menu item selects the residue under it.

/** @internal Whether `e` happened on the view itself (`container` or one of its own canvases), not in an overlay in it. */
export function isViewEvent(container: HTMLElement, e: Event): boolean {
  const t = e.target as Node | null;
  return t === container || (!!t && t.parentNode === container && t.nodeName === 'CANVAS');
}
