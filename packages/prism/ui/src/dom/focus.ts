/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyboard's focus around what opens and closes (README §3.7): a dialog keeps Tab inside it and gives the focus back
// to what opened it; a menu moves between its items by the arrows (WAI-ARIA's menu pattern); a tab list is one Tab stop
// whose tabs the arrows move between (the tabs pattern, manual activation). Across shadow roots: a viewer's shell lives
// in one.

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

/** A menu's items that can be used now (not disabled, not hidden), in order. */
export function menuItems(menu: Element): HTMLElement[] {
  return [...menu.querySelectorAll<HTMLElement>(MenuItem)].filter(
    (i) => !i.hasAttribute('disabled') && i.getAttribute('aria-disabled') !== 'true' && !i.closest('[hidden]'),
  );
}

/** Focuses a menu's first item (or its checked one, `checked`); true when there was one. */
export function focusFirstItem(menu: Element, checked = false): boolean {
  const items = menuItems(menu);
  const item = (checked ? items.find((i) => i.getAttribute('aria-checked') === 'true') : undefined) ?? items[0];
  item?.focus();
  return !!item;
}

/**
 * The arrows on a `role="menu"`: ↓ and ↑ move between its items (wrapping), Home and End jump. Bind as
 * `onKeyDown={menuKeys}` on the menu.
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

/**
 * A tab list as one Tab stop (roving `tabindex`: the focused tab, else the selected one, else the first) whose tabs the
 * arrows move between (both axes; Home, End); Enter or Space activates the focused tab (it is a button). Bind as
 * `ref={(el) => onCleanup(rovingTabs(el))}`; it follows tabs added, removed and selected.
 */
export function rovingTabs(list: HTMLElement): () => void {
  const tabs = () => [...list.querySelectorAll<HTMLElement>('[role="tab"]')].filter((t) => !t.closest('[hidden]'));
  const settle = () => {
    const all = tabs();
    const active = deepActiveElement(list.ownerDocument);
    const current =
      all.find((t) => t === active) ?? all.find((t) => t.getAttribute('aria-selected') === 'true') ?? all[0];
    for (const t of all) {
      const v = t === current ? 0 : -1;
      if (t.tabIndex !== v) t.tabIndex = v;
    }
  };
  const onKey = (e: KeyboardEvent) => {
    if (e.defaultPrevented) return;
    const all = tabs();
    const at = all.indexOf(e.target as HTMLElement);
    if (at < 0) return;
    let next: HTMLElement | undefined;
    switch (e.key) {
      case 'ArrowDown':
      case 'ArrowRight':
        next = all[(at + 1) % all.length];
        break;
      case 'ArrowUp':
      case 'ArrowLeft':
        next = all[(at - 1 + all.length) % all.length];
        break;
      case 'Home':
        next = all[0];
        break;
      case 'End':
        next = all[all.length - 1];
        break;
      default:
        return;
    }
    e.preventDefault();
    next.focus();
  };
  list.addEventListener('keydown', onKey);
  list.addEventListener('focusin', settle);
  const observer = new MutationObserver(settle);
  observer.observe(list, { childList: true, subtree: true, attributes: true, attributeFilter: ['aria-selected'] });
  settle();
  return () => {
    list.removeEventListener('keydown', onKey);
    list.removeEventListener('focusin', settle);
    observer.disconnect();
  };
}

let ids = 0;

/** A new element id, unique on the page (ARIA's references: `aria-labelledby`, `aria-activedescendant`). */
export function uniqueId(prefix = 'msp'): string {
  return `${prefix}-${++ids}`;
}
