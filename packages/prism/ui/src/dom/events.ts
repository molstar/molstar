/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Events across a shadow root, and keys scoped to one viewer (spec library.md §3, README decision 29). A listener on
// the document sees an event from inside a `<prism-viewer>` retargeted to the viewer's element: the composed path has
// the real target. Keys go to the viewer that has the focus, or the one under the pointer while the focus is on the
// page's body, never to every viewer on the page.

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

/**
 * The element a view's keys are scoped to: the shadow root's host (a `<prism-viewer>`), else the prism-ui root it is
 * in; undefined outside both (a story, a test): no scoping.
 */
export function viewerHost(node: Node): Element | undefined {
  const root = node.getRootNode();
  if (typeof ShadowRoot !== 'undefined' && root instanceof ShadowRoot) return root.host;
  return (node as Element).closest?.('.msp-root') ?? node.parentElement?.closest('.msp-root') ?? undefined;
}

/**
 * Whether a key event is this view's viewer's: the focus is in it, or nothing on the page has the focus and the
 * pointer is over it. Without a viewer (`viewerHost` undefined) every key is.
 */
export function ownsKey(node: Node, e: KeyboardEvent): boolean {
  const host = viewerHost(node);
  if (!host) return true;
  if (isInside(e, host)) return true;
  const doc = host.ownerDocument;
  const focused = doc.activeElement;
  if (focused && focused !== doc.body && focused !== doc.documentElement) return false;
  return host.matches(':hover');
}

/** Whether keys typed at `target` belong to it (a text field, an editable). */
export function isEditable(target: EventTarget | null): boolean {
  if (typeof HTMLElement === 'undefined' || !(target instanceof HTMLElement)) return false;
  if (target.isContentEditable) return true;
  if (target instanceof HTMLTextAreaElement || target instanceof HTMLSelectElement) return true;
  return target instanceof HTMLInputElement && !['checkbox', 'radio', 'range', 'button', 'color'].includes(target.type);
}
