/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observeResize, on } from '../../dom/scoped.js';
import { classes } from './_generated/tooltip.module.css.js';

/** Shown after the pointer rests this long; moving to a neighbour while one is shown switches at once. */
const ShowDelayMs = 350;

/**
 * Tooltips for icon buttons: any element in `root` with `data-tip` (the text) gets one, placed by `data-tip-side`
 * (`left` for the rail and the canvas stack, `top` for the bottom of the canvas, `bottom` by default), on hover or keyboard focus (never for an element in an `inert` part: hidden UI takes no pointer). One element is reused,
 * positioned in the viewport; buttons keep their `aria-label` for assistive technology (no native `title` as well).
 * `data-tip-overflow` (no text) shows the element's own text, only while it is cut off (a long label, faded at its end).
 * `data-tip-row` (a tree row) shows its parts' text (`data-tip-part`: the label, then the description, muted) while
 * one is cut off, from a pointer over a part: a card beside the row, outside its panel (`role="tabpanel"`, or
 * `data-tip-edge`) and clear of its edge, at the row's height and its parts' sizes and weights, its text level with the
 * row's (on its left, or its right where there is no room).
 * Called in a scope (its listeners are `on`'s: paused and released with it); returns the uninstall function.
 */
export function installTooltips(root: HTMLElement): () => void {
  const tip = document.createElement('div');
  tip.className = classes.tooltip;
  tip.setAttribute('role', 'tooltip');
  tip.hidden = true;
  root.append(tip);

  let target: HTMLElement | undefined;
  // Whether the pointer (not the keyboard focus) brought the tip: it then ends when the pointer is not over the anchor.
  let byPointer = false;
  let timer: ReturnType<typeof setTimeout> | undefined;
  // While one is shown: the anchor's size. A view closing or hiding takes its anchor out without a pointerout, and the
  // tip would stay; an anchor that is gone or has no size ends it (no listener per event, one observer on one element).
  const watch = observeResize(undefined, () => !usable(target) && hide());

  const place = (el: HTMLElement) => {
    const row = el.dataset.tipRow !== undefined;
    if (row) fillRow(tip, el);
    else tip.textContent = el.dataset.tip ?? el.textContent ?? '';
    tip.dataset.kind = row ? 'row' : '';
    tip.style.paddingBlock = '';
    tip.hidden = false;
    // A row's card on one line is the row's height, its text level with the row's; wrapped (a long text), it keeps its
    // own padding.
    if (row) {
      const line = Number.parseFloat(getComputedStyle(tip).lineHeight);
      const height = el.getBoundingClientRect().height;
      if (line > 0 && height > line) {
        tip.style.paddingBlock = `${(height - line) / 2}px`;
        if (tip.getBoundingClientRect().height > height + 0.5) tip.style.paddingBlock = '';
      }
    }
    const side = row ? 'beside' : (el.dataset.tipSide ?? 'bottom');
    const r = el.getBoundingClientRect();
    const t = tip.getBoundingClientRect();
    const gap = 8;
    let x: number;
    let y: number;
    let panel: DOMRect | undefined;
    if (side === 'beside') {
      // Outside the panel the row is in (its tab panel), clear of its edge; on its left, or its right without room.
      const edge = el.closest('[data-tip-edge], [role="tabpanel"]')?.getBoundingClientRect() ?? r;
      panel = edge;
      const clear = 8;
      x = edge.left - clear - t.width;
      if (x < 4) x = edge.right + clear;
      y = r.top + (r.height - t.height) / 2;
    } else if (side === 'left') {
      x = r.left - gap - t.width;
      y = r.top + (r.height - t.height) / 2;
    } else if (side === 'top') {
      x = r.left + (r.width - t.width) / 2;
      y = r.top - gap - t.height;
    } else {
      x = r.left + (r.width - t.width) / 2;
      y = r.bottom + gap;
    }
    x = Math.max(4, Math.min(x, window.innerWidth - t.width - 4));
    y = Math.max(4, Math.min(y, window.innerHeight - t.height - 4));
    // A card beside its panel never lies over it: with no room outside either edge there is no card (it would cover the
    // next rows).
    if (panel && x < panel.right && x + t.width > panel.left && panel !== r) return hide();
    tip.style.transform = `translate(${Math.round(x)}px, ${Math.round(y)}px)`;
    tip.dataset.side = side;
  };

  const hide = () => {
    if (timer !== undefined) clearTimeout(timer);
    timer = undefined;
    watch?.disconnect();
    target = undefined;
    tip.hidden = true;
  };

  const show = (el: HTMLElement, immediate: boolean, pointer: boolean) => {
    if (timer !== undefined) clearTimeout(timer);
    byPointer = pointer;
    // Its own text, only when cut off (by more than rounding: scrollWidth rounds a fractional width up, clientWidth
    // down, so text that fits can read 1 px over).
    if (el.dataset.tip === undefined && !cutOff(el)) {
      target = undefined;
      tip.hidden = true;
      return;
    }
    target = el;
    const showNow = () => {
      // The anchor may have gone while the pointer rested (a panel closed): no tip for what is not there.
      if (target !== el || !usable(el)) return hide();
      place(el);
      if (tip.hidden) return;
      watch?.disconnect();
      watch?.observe(el);
    };
    if (immediate) showNow();
    else timer = setTimeout(showNow, ShowDelayMs);
  };

  const tipOf = (e: Event) => {
    const at = e.target as Element | null;
    const el = at?.closest?.<HTMLElement>('[data-tip], [data-tip-overflow], [data-tip-row]') ?? undefined;
    // Hidden UI is inert (a panel behind a card, a folded section): no tip for it, whatever hit it.
    if (el?.closest('[inert]')) return undefined;
    // A row answers only over its text (not its chevron or its buttons).
    if (el?.dataset.tipRow !== undefined && !at?.closest('[data-tip-part]')) return undefined;
    return el;
  };

  const onOver = (e: PointerEvent) => {
    // A pointer event anywhere: a tip whose anchor was removed ends now (a connected check, nothing more).
    if (target && !target.isConnected) hide();
    const el = tipOf(e);
    if (!el || el === target) return;
    show(el, !tip.hidden, true);
  };
  const onOut = (e: PointerEvent) => {
    const el = tipOf(e);
    if (el && !el.contains(e.relatedTarget as Node)) hide();
  };
  const onFocus = (e: FocusEvent) => {
    const el = tipOf(e);
    if (el?.matches(':focus-visible')) show(el, true, false);
  };

  const onMove = (e: PointerEvent) => {
    // The anchor went, or is no longer what is under the pointer (it was hidden, made inert, covered): the tip ends. A
    // contains check on the one anchor, only while a tip is pending or shown.
    if (target && (!target.isConnected || (byPointer && !target.contains(e.target as Node)))) hide();
  };
  // Scoped (rendering.md §2.4): detached while the scope is paused, gone with it (or with the returned function).
  const offs = [
    on(root, 'pointerover', onOver),
    on(root, 'pointermove', onMove, { passive: true }),
    on(root, 'pointerout', onOut),
    on(root, 'focusin', onFocus),
    on(root, 'focusout', hide),
    on(root, 'pointerdown', hide),
    on(window, 'keydown', hide),
    on(window, 'scroll', hide, { capture: true }),
  ];
  return () => {
    hide();
    tip.remove();
    watch?.disconnect();
    for (const off of offs) off();
  };
}

/**
 * Whether the anchor can still carry a tip: in the document, rendered (not under a `display: none` or `hidden`), and not
 * in an inert part (a hidden panel).
 */
function usable(el: HTMLElement | undefined): el is HTMLElement {
  return !!el && el.isConnected && (el.checkVisibility?.() ?? true) && !el.closest('[inert]');
}

/** Cut off by more than rounding (scrollWidth rounds a fractional width up, clientWidth down: text that fits can read
 * 1 px over); a row, when one of its parts is. */
function cutOff(el: HTMLElement): boolean {
  if (el.dataset.tipRow === undefined) return el.scrollWidth - el.clientWidth > 1;
  for (const part of el.querySelectorAll<HTMLElement>('[data-tip-part]')) {
    if (part.scrollWidth - part.clientWidth > 1) return true;
  }
  return false;
}

/** A row's tip: its parts' text, the first as it is, the rest muted, each in its part's size and weight. */
function fillRow(tip: HTMLElement, row: HTMLElement) {
  tip.textContent = '';
  let first = true;
  for (const part of row.querySelectorAll<HTMLElement>('[data-tip-part]')) {
    const text = part.textContent?.trim();
    if (!text) continue;
    const span = document.createElement('span');
    span.textContent = text;
    const font = getComputedStyle(part);
    span.style.fontSize = font.fontSize;
    span.style.fontWeight = font.fontWeight;
    if (!first) span.className = classes.muted;
    if (!first) tip.append(' ');
    tip.append(span);
    first = false;
  }
}
