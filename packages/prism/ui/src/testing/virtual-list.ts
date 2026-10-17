/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// For tests of what a `VirtualList` shows: its rows are placed by their own transforms and recycled ones stay in the
// DOM hidden when put away (view-models.md §3), so the DOM order is not the list's.

/** The row element an element is in (the list's positioned wrapper), or undefined. */
function rowOf(el: Element): HTMLElement | undefined {
  return el.closest<HTMLElement>('[style*="translateY"]') ?? undefined;
}

/** A row's offset (its `translateY`, px). */
function offset(row: HTMLElement | undefined): number {
  const m = row ? /translateY\((-?[\d.]+)px\)/.exec(row.style.transform) : null;
  return m ? Number.parseFloat(m[1]) : 0;
}

/** The elements matching `selector` under `root` that a virtual list shows, in the list's order (not hidden). */
export function shownInOrder<E extends Element = HTMLElement>(root: ParentNode, selector: string): E[] {
  const found = [...root.querySelectorAll<E>(selector)].filter((el) => !el.closest('[hidden]'));
  return found
    .map((el, i) => ({ el, i, y: offset(rowOf(el)) }))
    .sort((a, b) => a.y - b.y || a.i - b.i)
    .map((e) => e.el);
}
