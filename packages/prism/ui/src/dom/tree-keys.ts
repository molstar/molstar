/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyboard of an ARIA tree (accessibility.md, WAI-ARIA's tree pattern): a tree is reached by Tab once, then the
// arrows move between its items (Down, Up), open and enter (Right), close and leave (Left), Home and End jump, Enter
// and Space do what a click on the row does. Items are the `role="treeitem"` elements under the tree; each item's row
// is its `[data-tree-row]` element (clicked on Enter) and its open/close button its `[data-tree-toggle]`. Focus roves:
// the focused item has `tabindex="0"`, the others `-1`. A tree only binds `onKeyDown={treeKeys}` and
// `onFocusIn={treeFocus}`.

const Item = '[role="treeitem"]';

function items(tree: HTMLElement): HTMLElement[] {
  return [...tree.querySelectorAll<HTMLElement>(Item)];
}

function focusItem(tree: HTMLElement, item: HTMLElement | undefined): void {
  if (!item) return;
  for (const i of items(tree)) i.tabIndex = i === item ? 0 : -1;
  item.focus();
}

/** Tab into the tree: the first item, or the one picked (`aria-selected`), takes the focus. */
export function treeFocus(e: FocusEvent): void {
  const tree = e.currentTarget as HTMLElement;
  if (e.target !== tree) return;
  const all = items(tree);
  focusItem(tree, all.find((i) => i.getAttribute('aria-selected') === 'true') ?? all[0]);
}

/** The arrows, Home, End, Enter and Space on an item of the tree. */
export function treeKeys(e: KeyboardEvent): void {
  const tree = e.currentTarget as HTMLElement;
  const item = (e.target as HTMLElement).closest<HTMLElement>(Item);
  // Keys on a button inside a row (the eye, the menu) are the button's.
  if (!item || e.target !== item) return;
  const all = items(tree);
  const at = all.indexOf(item);
  const expanded = item.getAttribute('aria-expanded');
  const row = item.querySelector<HTMLElement>(':scope > [data-tree-row]');
  const toggle = row?.querySelector<HTMLElement>('[data-tree-toggle]');
  let handled = true;
  switch (e.key) {
    case 'ArrowDown':
      focusItem(tree, all[at + 1]);
      break;
    case 'ArrowUp':
      focusItem(tree, all[at - 1]);
      break;
    case 'Home':
      focusItem(tree, all[0]);
      break;
    case 'End':
      focusItem(tree, all[all.length - 1]);
      break;
    case 'ArrowRight':
      if (expanded === 'false') toggle?.click();
      else if (expanded === 'true') focusItem(tree, item.querySelector<HTMLElement>(Item) ?? undefined);
      break;
    case 'ArrowLeft':
      if (expanded === 'true') toggle?.click();
      else focusItem(tree, item.parentElement?.closest<HTMLElement>(Item) ?? undefined);
      break;
    case 'Enter':
    case ' ':
      row?.click();
      break;
    default:
      handled = false;
  }
  if (handled) e.preventDefault();
}
