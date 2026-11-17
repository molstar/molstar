/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { deepActiveElement, eventTarget, isInside, on, onActive, onCleanup } from '../../dom/index.js';

// How an overlay (a popover, a select's options, a menu, a toolbar's menu, a context menu) closes: one rule for all,
// so layers opened from one another (a select in a popover, a submenu) agree on whose a press or a key is.

/** Why a layer closed: a press outside it, or a key (Escape, Tab: from the keyboard). */
export type DismissReason = 'outside' | 'escape' | 'tab';

export interface DismissOptions {
  /** The layer's element: a press or a key in it is its own. A getter: a ref sets it after the call. */
  layer(): Element | undefined;
  /** What opened it (its button): a press on it is the button's (it toggles); Escape on it is the layer's. */
  anchor?: Element | (() => Element | undefined);
  /** Closes it. */
  close(reason: DismissReason): void;
  /**
   * Where the focus goes after Escape or Tab closed it while the focus was in it (its button), so it never drops to
   * the page; without it the caller places it (told by the reason). Tab then moves on from there.
   */
  refocus?(reason: DismissReason): HTMLElement | undefined;
  /** Tab from inside it closes it (a menu, a list: its rows are no Tab stops); default false (a dialog keeps it). */
  tab?: boolean;
  /** Escape closes it (default true). */
  escape?: boolean;
  /** A press outside closes it (default true). */
  outside?: boolean;
  /**
   * Whether Escape pressed outside it and its anchor is its own. Default: only with the focus lost (on the body); a key
   * with the focus on another surface (a dialog, a field) is that surface's.
   */
  ownsKey?(e: KeyboardEvent): boolean;
}

interface Layer {
  readonly options: DismissOptions;
  active: boolean;
}

/** The open layers, in the order they opened: the last is the innermost. */
const layers: Layer[] = [];
/** Keys already taken by a layer: every layer listens, one acts. */
const taken = new WeakSet<Event>();

function anchorOf(layer: Layer): Element | undefined {
  const a = layer.options.anchor;
  return typeof a === 'function' ? a() : a;
}

/**
 * Whether the event happened in `el` and not in a floating layer inside it (a top-layer `[popover]` between the target
 * and `el`: a select's options opened from a panel, registered as a layer or not).
 */
function inOwn(e: Event, el: Element | undefined): boolean {
  if (!el) return false;
  const path = typeof e.composedPath === 'function' ? e.composedPath() : [];
  if (!path.length) return el.contains(e.target as Node);
  for (const n of path) {
    if (n === el) return true;
    if (n instanceof Element && n.hasAttribute('popover') && !n.contains(el)) return false;
  }
  return false;
}

function focusLost(e: Event): boolean {
  const t = eventTarget(e);
  if (!(t instanceof Element)) return true;
  const doc = t.ownerDocument;
  return t === doc.body || t === doc.documentElement;
}

/** Whether the key is the layer's: in it, on its anchor (Escape), or (Escape) one it owns from outside. */
function claims(layer: Layer, e: KeyboardEvent): boolean {
  if (inOwn(e, layer.options.layer())) return true;
  if (e.key === 'Tab') return false;
  if (isInside(e, anchorOf(layer))) return true;
  return layer.options.ownsKey ? layer.options.ownsKey(e) : focusLost(e);
}

/** The layer a key closes: the innermost open one that takes that key and claims it. */
function ownerOf(e: KeyboardEvent): Layer | undefined {
  const tab = e.key === 'Tab';
  for (let i = layers.length - 1; i >= 0; i--) {
    const layer = layers[i];
    if (!layer.active) continue;
    if (tab ? !layer.options.tab : layer.options.escape === false) continue;
    if (claims(layer, e)) return layer;
  }
  return undefined;
}

/**
 * Whether a press outside `self` is in a layer opened after it (its select's options, its submenu, a context menu over
 * it), or in a floating layer (`[popover]`) that is no layer of this registry (a list a control floats on its own):
 * that layer's, not a press outside. A press in a layer opened before it is outside.
 */
function inLaterLayer(e: Event, self: Layer): boolean {
  const at = layers.indexOf(self);
  for (let i = at + 1; i > 0 && i < layers.length; i++) {
    if (isInside(e, layers[i].options.layer())) return true;
  }
  const own = self.options.layer();
  const path = typeof e.composedPath === 'function' ? e.composedPath() : [];
  for (const n of path) {
    if (!(n instanceof Element) || !n.hasAttribute('popover') || (own && n.contains(own))) continue;
    if (!layers.some((l) => n.contains(l.options.layer() ?? null))) return true;
  }
  return false;
}

/**
 * Closes an overlay on a press outside it and its anchor (a press in a layer opened after it, nested in it, is that
 * layer's), on Escape and, with `tab`, on Tab from inside it. Layers nest: a key goes to the innermost open layer that
 * claims it, so Escape in a select's options closes them and not the popover they are in. A keyboard close moves the
 * focus to `refocus()` when the focus was in the layer (or lost). In a scope: its listeners are `on`'s (paused with
 * it), the layer leaves the stack with it.
 */
export function dismissable(options: DismissOptions): void {
  const self: Layer = { options, active: true };
  layers.push(self);
  onCleanup(() => {
    const i = layers.indexOf(self);
    if (i >= 0) layers.splice(i, 1);
  });
  onActive((active) => {
    self.active = active;
  });

  const onKey = (e: KeyboardEvent) => {
    if ((e.key !== 'Escape' && e.key !== 'Tab') || e.defaultPrevented || taken.has(e)) return;
    if (ownerOf(e) !== self) return;
    taken.add(e);
    const el = options.layer();
    const reason: DismissReason = e.key === 'Escape' ? 'escape' : 'tab';
    // Escape is this layer's alone (not also the dialog it is in); Tab goes on, from where the focus is put back.
    if (reason === 'escape') {
      e.preventDefault();
      e.stopPropagation();
    }
    const before = deepActiveElement();
    const hadFocus = !before || !!el?.contains(before);
    options.close(reason);
    if (!hadFocus || !options.refocus) return;
    const now = deepActiveElement();
    if (now && !el?.contains(now)) return;
    options.refocus(reason)?.focus({ preventScroll: true });
  };

  const onDown = (e: PointerEvent) => {
    if (options.outside === false || !self.active) return;
    if (isInside(e, options.layer()) || isInside(e, anchorOf(self))) return;
    if (inLaterLayer(e, self)) return;
    options.close('outside');
  };

  on(document, 'keydown', onKey, { capture: true });
  on(document, 'pointerdown', onDown, { capture: true });
}

/** Moves the focus to `target` when it is in `layer` or lost (a pick closed the layer that had it). */
export function refocusFrom(layer: Element | undefined, target: HTMLElement | undefined): void {
  const now = deepActiveElement();
  if (now && !layer?.contains(now)) return;
  target?.focus({ preventScroll: true });
}
