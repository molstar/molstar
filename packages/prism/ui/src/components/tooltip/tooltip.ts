/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from './_generated/tooltip.module.css.js';

/** Shown after the pointer rests this long; moving to a neighbour while one is shown switches at once. */
const ShowDelayMs = 350;

/**
 * Tooltips for icon buttons: any element in `root` with `data-tip` (the text) gets one, placed by `data-tip-side`
 * (`left` for the rail and the canvas stack, `top` for the bottom of the canvas, `bottom` by default), on hover or keyboard focus. One element is reused,
 * positioned in the viewport; buttons keep their `aria-label` for assistive technology (no native `title` as well).
 * `data-tip-overflow` (no text) shows the element's own text, only while it is cut off (a long label, faded at its end).
 * `data-tip-row` (a tree row) shows its parts' text (`data-tip-part`: the label, then the description, muted) while
 * one is cut off, from a pointer over a part: a card beside the row, outside its panel (`role="tabpanel"`, or
 * `data-tip-edge`) and clear of its edge, at the row's text size (on its left, or its right where there is no room).
 * Returns the uninstall function.
 */
export function installTooltips(root: HTMLElement): () => void {
  const tip = document.createElement('div');
  tip.className = classes.tooltip;
  tip.setAttribute('role', 'tooltip');
  tip.hidden = true;
  root.append(tip);

  let target: HTMLElement | undefined;
  let timer: ReturnType<typeof setTimeout> | undefined;

  const place = (el: HTMLElement) => {
    const row = el.dataset.tipRow !== undefined;
    if (row) fillRow(tip, el);
    else tip.textContent = el.dataset.tip ?? el.textContent ?? '';
    tip.dataset.kind = row ? 'row' : '';
    tip.hidden = false;
    const side = row ? 'beside' : (el.dataset.tipSide ?? 'bottom');
    const r = el.getBoundingClientRect();
    const t = tip.getBoundingClientRect();
    const gap = 8;
    let x: number;
    let y: number;
    if (side === 'beside') {
      // Outside the panel the row is in (its tab panel), clear of its edge; on its left, or its right without room.
      const edge = el.closest('[data-tip-edge], [role="tabpanel"]')?.getBoundingClientRect() ?? r;
      const clear = 12;
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
    tip.style.transform = `translate(${Math.round(x)}px, ${Math.round(y)}px)`;
    tip.dataset.side = side;
  };

  const hide = () => {
    if (timer !== undefined) clearTimeout(timer);
    timer = undefined;
    target = undefined;
    tip.hidden = true;
  };

  const show = (el: HTMLElement, immediate: boolean) => {
    if (timer !== undefined) clearTimeout(timer);
    // Its own text, only when cut off (by more than rounding: scrollWidth rounds a fractional width up, clientWidth
    // down, so text that fits can read 1 px over).
    if (el.dataset.tip === undefined && !cutOff(el)) {
      target = undefined;
      tip.hidden = true;
      return;
    }
    target = el;
    if (immediate) place(el);
    else timer = setTimeout(() => target === el && place(el), ShowDelayMs);
  };

  const tipOf = (e: Event) => {
    const at = e.target as Element | null;
    const el = at?.closest?.<HTMLElement>('[data-tip], [data-tip-overflow], [data-tip-row]') ?? undefined;
    // A row answers only over its text (not its chevron or its buttons).
    if (el?.dataset.tipRow !== undefined && !at?.closest('[data-tip-part]')) return undefined;
    return el;
  };

  const onOver = (e: PointerEvent) => {
    const el = tipOf(e);
    if (!el || el === target) return;
    show(el, !tip.hidden);
  };
  const onOut = (e: PointerEvent) => {
    const el = tipOf(e);
    if (el && !el.contains(e.relatedTarget as Node)) hide();
  };
  const onFocus = (e: FocusEvent) => {
    const el = tipOf(e);
    if (el?.matches(':focus-visible')) show(el, true);
  };

  root.addEventListener('pointerover', onOver);
  root.addEventListener('pointerout', onOut);
  root.addEventListener('focusin', onFocus);
  root.addEventListener('focusout', hide);
  root.addEventListener('pointerdown', hide);
  window.addEventListener('keydown', hide);
  window.addEventListener('scroll', hide, true);
  return () => {
    hide();
    tip.remove();
    root.removeEventListener('pointerover', onOver);
    root.removeEventListener('pointerout', onOut);
    root.removeEventListener('focusin', onFocus);
    root.removeEventListener('focusout', hide);
    root.removeEventListener('pointerdown', hide);
    window.removeEventListener('keydown', hide);
    window.removeEventListener('scroll', hide, true);
  };
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

/** A row's tip: its parts' text, the first as it is, the rest muted. */
function fillRow(tip: HTMLElement, row: HTMLElement) {
  tip.textContent = '';
  let first = true;
  for (const part of row.querySelectorAll<HTMLElement>('[data-tip-part]')) {
    const text = part.textContent?.trim();
    if (!text) continue;
    const span = document.createElement('span');
    span.textContent = text;
    if (!first) span.className = classes.muted;
    if (!first) tip.append(' ');
    tip.append(span);
    first = false;
  }
}
