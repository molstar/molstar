/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { onActive, onCleanup } from './scope.js';
import { on } from './scoped.js';

// The keyboard's focus around what opens and closes (dom §7): a dialog keeps Tab inside it and gives the
// focus back to what opened it; a menu moves between its items by the arrows (WAI-ARIA's menu pattern: a disabled item
// is focusable, so its reason is read); a tab list or a toolbar is one Tab stop whose items the arrows move between
// (roving `tabindex`; right-to-left text flips ← →). Across shadow roots: a viewer's shell lives in one.

const Focusable =
  'a[href], button:not([disabled]), input:not([disabled]):not([type="hidden"]), select:not([disabled]), ' +
  'textarea:not([disabled]), [tabindex]:not([tabindex="-1"]), [contenteditable="true"]';

const MenuItem = '[role="menuitem"], [role="menuitemcheckbox"], [role="menuitemradio"]';

/** The focused element of the page, through shadow roots; undefined when nothing (or only the body) has it. */
export function deepActiveElement(doc: Document = document): HTMLElement | undefined {
  let active: Element | null = doc.activeElement;
  while (active?.shadowRoot?.activeElement) active = active.shadowRoot.activeElement;
  if (!active || active === doc.body || active === doc.documentElement) return undefined;
  // A shadow host whose focused element went (removed): the focus is lost, not on the host.
  if (active.shadowRoot && !active.shadowRoot.activeElement && (active as HTMLElement).tabIndex < 0) return undefined;
  return active as HTMLElement;
}

/** What can take the focus by Tab inside `root`, in order (nothing hidden or inert; no layout read). */
export function focusables(root: Element): HTMLElement[] {
  return [...root.querySelectorAll<HTMLElement>(Focusable)].filter((el) => !el.closest('[hidden], [inert]'));
}

/** Keeps Tab inside `root` (a dialog, a card): past the last focusable to the first, and back with Shift. */
export function trapTab(e: KeyboardEvent, root: HTMLElement): void {
  if (e.key !== 'Tab' || e.defaultPrevented) return;
  const all = focusables(root);
  if (!all.length) return;
  const active = deepActiveElement(root.ownerDocument);
  const inside = !!active && root.contains(active);
  if (e.shiftKey && (!inside || active === all[0])) {
    e.preventDefault();
    all[all.length - 1].focus();
  } else if (!e.shiftKey && (!inside || active === all[all.length - 1])) {
    e.preventDefault();
    all[0].focus();
  }
}

/**
 * Remembers what has the focus now (call it as a dialog or menu opens); the returned function gives it back (call it
 * as it closes), or to `fallback` when that is gone. It does nothing when the focus has moved on to something else
 * meanwhile: only a focus that is lost (on the body) or still inside `within` (the closing element) goes back.
 */
export function rememberFocus(fallback?: () => HTMLElement | null | undefined): (within?: Element | null) => void {
  const before = deepActiveElement();
  return (within) => {
    const now = deepActiveElement();
    if (now && !(within && within.contains(now))) return;
    const target = before?.isConnected ? before : fallback?.();
    if (target?.isConnected) target.focus({ preventScroll: true });
  };
}

/** Whether `el` was focused from the keyboard (its focus ring shows), not by a click. */
export function focusedByKeyboard(el: Element | null | undefined): boolean {
  if (!el) return false;
  try {
    return el.matches(':focus-visible');
  } catch {
    // An environment without :focus-visible (jsdom): only a focused element counts.
    return deepActiveElement(el.ownerDocument) === el;
  }
}

/**
 * A menu's items the arrows move over, in order: not hidden, not `disabled` (it cannot take the focus); an
 * `aria-disabled` one is (APG menus: focusable, not activated, so its reason is reachable).
 */
export function menuItems(menu: Element): HTMLElement[] {
  return [...menu.querySelectorAll<HTMLElement>(MenuItem)].filter(
    (i) => !i.hasAttribute('disabled') && !i.closest('[hidden]'),
  );
}

const usable = (i: Element) => i.getAttribute('aria-disabled') !== 'true';

/**
 * Focuses a menu's first usable item (or its checked one, `checked`), else its first; true when there was one. Opening
 * on a usable item: Enter there runs something.
 */
export function focusFirstItem(menu: Element, checked = false): boolean {
  const items = menuItems(menu);
  const item =
    (checked ? items.find((i) => i.getAttribute('aria-checked') === 'true') : undefined) ??
    items.find(usable) ??
    items[0];
  item?.focus();
  return !!item;
}

/**
 * The arrows on a `role="menu"`: ↓ and ↑ move between its items (wrapping; disabled ones too), Home and End jump. Bind
 * as `onKeyDown={menuKeys}` on the menu.
 */
export function menuKeys(e: KeyboardEvent): void {
  if (e.defaultPrevented) return;
  const menu = e.currentTarget as HTMLElement;
  const items = menuItems(menu);
  if (!items.length) return;
  const at = items.indexOf(deepActiveElement(menu.ownerDocument) as HTMLElement);
  const n = items.length;
  let next: HTMLElement | undefined;
  switch (e.key) {
    case 'ArrowDown':
      next = items[at < 0 ? 0 : (at + 1) % n];
      break;
    case 'ArrowUp':
      next = items[at < 0 ? n - 1 : (at - 1 + n) % n];
      break;
    case 'Home':
      next = items[0];
      break;
    case 'End':
      next = items[n - 1];
      break;
    default:
      return;
  }
  e.preventDefault();
  next.focus();
}

/** Whether `el`'s text runs right to left (`dir="rtl"`, or the CSS direction): ← and → swap along a row. */
export function isRtl(el: Element): boolean {
  const dir = el.closest('[dir]')?.getAttribute('dir')?.toLowerCase();
  if (dir === 'rtl') return true;
  if (dir === 'ltr') return false;
  return typeof getComputedStyle === 'function' && getComputedStyle(el).direction === 'rtl';
}

export interface RovingFocusOptions {
  /** The items (a selector, under the root); one inside `[hidden]` is skipped. */
  items: string;
  /**
   * The arrows that move: `horizontal` ← → (swapped right to left), `vertical` ↓ ↑, `both` (a tab list) all four.
   * Home and End jump.
   */
  axis: 'horizontal' | 'vertical' | 'both';
  /**
   * The Tab stop when no item has the focus, before the first: a tab list's selected tab. Default: the last focused
   * that still shows.
   */
  preferred?(items: readonly HTMLElement[]): HTMLElement | undefined;
  /** Attributes whose change can move the Tab stop (besides items added, removed or `hidden`). */
  attributes?: readonly string[];
}

/**
 * A group of items as one Tab stop (roving `tabindex`: the focused item, else `preferred`, else the first) whose items
 * the arrows move between (wrapping). Enter or Space activates the focused item (a button). In a scope: its listeners
 * are `on`'s, its `MutationObserver` (items added, removed, hidden) is connected while the scope is active; the returned
 * function ends it sooner.
 */
export function rovingFocus(root: HTMLElement, options: RovingFocusOptions): () => void {
  const items = () => [...root.querySelectorAll<HTMLElement>(options.items)].filter((t) => !t.closest('[hidden]'));
  let last: HTMLElement | undefined;
  const preferred =
    options.preferred ?? ((all: readonly HTMLElement[]) => (last && all.includes(last) ? last : undefined));
  const settle = () => {
    const all = items();
    const active = deepActiveElement(root.ownerDocument);
    const current = all.find((t) => t === active) ?? preferred(all) ?? all[0];
    last = current;
    for (const t of all) {
      const v = t === current ? 0 : -1;
      if (t.tabIndex !== v) t.tabIndex = v;
    }
  };
  const onKey = (e: KeyboardEvent) => {
    if (e.defaultPrevented) return;
    const all = items();
    const at = all.indexOf(e.target as HTMLElement);
    if (at < 0) return;
    const { axis } = options;
    const flip = axis !== 'vertical' && isRtl(root);
    let step = 0;
    if (e.key === 'Home') step = -Infinity;
    else if (e.key === 'End') step = Infinity;
    else if (axis !== 'horizontal' && (e.key === 'ArrowDown' || e.key === 'ArrowUp'))
      step = e.key === 'ArrowDown' ? 1 : -1;
    else if (axis !== 'vertical' && (e.key === 'ArrowRight' || e.key === 'ArrowLeft'))
      step = (e.key === 'ArrowRight') !== flip ? 1 : -1;
    else return;
    const n = all.length;
    const next = step === -Infinity ? all[0] : step === Infinity ? all[n - 1] : all[(at + step + n) % n];
    e.preventDefault();
    next.focus();
  };
  // Scoped (dom §6): detached while the scope is paused, gone with it (or with the returned function).
  const offKey = on(root, 'keydown', onKey);
  const offFocus = on(root, 'focusin', settle);
  // An observer is not a listener `on` can pause: connected while the scope is active.
  const observer = new MutationObserver(settle);
  const attributeFilter = ['hidden', ...(options.attributes ?? [])];
  let ended = false;
  onActive((active) => {
    if (active && !ended) {
      observer.observe(root, { childList: true, subtree: true, attributes: true, attributeFilter });
      settle();
    } else observer.disconnect();
  });
  const end = () => {
    ended = true;
    offKey();
    offFocus();
    observer.disconnect();
  };
  onCleanup(end);
  return end;
}

/**
 * A tab list as one Tab stop (`rovingFocus`: the focused tab, else the selected one, else the first) whose tabs the
 * arrows move between (both axes; Home, End); Enter or Space activates the focused tab (it is a button). Bind as
 * `ref={(el) => onCleanup(rovingTabs(el))}`, in a scope; it follows tabs added, removed and selected.
 */
export function rovingTabs(list: HTMLElement): () => void {
  return rovingFocus(list, {
    items: '[role="tab"]',
    axis: 'both',
    preferred: (all) => all.find((t) => t.getAttribute('aria-selected') === 'true'),
    attributes: ['aria-selected'],
  });
}

let ids = 0;

/** A new element id, unique on the page (ARIA's references: `aria-labelledby`, `aria-activedescendant`). */
export function uniqueId(prefix = 'id'): string {
  return `${prefix}-${++ids}`;
}
