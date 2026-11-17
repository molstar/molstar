/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Events across shadow roots (dom §7): a listener on the document sees an event from inside a shadow root retargeted to
// its host; the composed path has the real target.

/** The innermost target of an event, across shadow roots. */
export function eventTarget(e: Event): EventTarget | null {
  const path = typeof e.composedPath === 'function' ? e.composedPath() : [];
  return path[0] ?? e.target;
}

/** Whether an event happened inside `el` (across shadow roots). */
export function isInside(e: Event, el: Node | null | undefined): boolean {
  if (!el) return false;
  const path = typeof e.composedPath === 'function' ? e.composedPath() : [];
  if (path.length) return path.includes(el);
  return el.contains(e.target as Node);
}

/** The focused element as seen from `node`: through its shadow root, and into nested ones. */
export function activeElement(node: Node): Element | null {
  const root = node.getRootNode() as Document | ShadowRoot;
  let active = root.activeElement ?? null;
  while (active?.shadowRoot?.activeElement) active = active.shadowRoot.activeElement;
  return active;
}

/** Whether keys typed at `target` belong to it (a text field, an editable). */
export function isEditable(target: EventTarget | null): boolean {
  if (typeof HTMLElement === 'undefined' || !(target instanceof HTMLElement)) return false;
  if (target.isContentEditable) return true;
  if (target instanceof HTMLTextAreaElement || target instanceof HTMLSelectElement) return true;
  return target instanceof HTMLInputElement && !['checkbox', 'radio', 'range', 'button', 'color'].includes(target.type);
}
