/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MarkdownInline } from '@molstar/core/util/markdown';
import { cardHasLink, type HoverCardPortFactory, type HoverCardShowOptions, placeCard } from './hover-card.js';

// The hover card's default adapter (_spec/README.md "Ports", Hover card): `@molstar/ui`'s primitives over the canvas: a
// surface (`msui-surface`: its tokens, its squircle where the browser draws one, its shadow as what floats), a header
// with the title and a form of rows (`msui-form`, `msui-label`, `msui-value`). DOM over the canvas: never in the draw
// lists, never in an export. Hidden from assistive technology: the live region says the same words (Announcements).
// Placed by `placeCard` at its measured size (the place's suggestion is for an estimated size). Its width is capped
// (`maxWidth`) and bounds the value column: a long value wraps there (at `,` `;` `|` `:` first, anywhere in a token
// with none) to `hover.tooltipLines`, then an ellipsis. A row's Markdown (a tooltip's) is drawn as elements, never
// HTML; a card with a link takes the pointer (entered across the gap, its link run through the viewer's commands), one
// without never does.

/** The card's widest, CSS px: the value column wraps within it. */
export const CardMaxPx = 320;

/** Where a long value may break besides its spaces: after these, a zero-width space (shown, never read: aria-hidden). */
const Breaks = /([,;|:])(?=\S)/g;

/** A value's text with its break opportunities. */
export const breakable = (text: string) => text.replace(Breaks, '$1​');

/** A row's Markdown as elements (strong, emphasis, code, a link a button-like anchor, a break), into `into`. */
function drawRich(doc: Document, into: HTMLElement, nodes: readonly MarkdownInline[], link: (href: string) => void) {
  for (const n of nodes) {
    if (n.kind === 'text') into.append(doc.createTextNode(breakable(n.text)));
    else if (n.kind === 'code') {
      const c = doc.createElement('code');
      c.textContent = n.text;
      into.append(c);
    } else if (n.kind === 'break') into.append(doc.createElement('br'));
    else if (n.kind === 'link') {
      const a = doc.createElement('a');
      a.href = '#';
      a.dataset.seqLink = n.href;
      a.addEventListener('click', (e) => {
        e.preventDefault();
        link(n.href);
      });
      drawRich(doc, a, n.children, link);
      into.append(a);
    } else {
      const el = doc.createElement(n.kind === 'strong' ? 'strong' : 'em');
      drawRich(doc, el, n.children, link);
      into.append(el);
    }
  }
}

/** The value column's wrap: within the card's width, to `lines` lines, then an ellipsis. */
function clamp(el: HTMLElement, lines: number) {
  Object.assign(el.style, {
    display: '-webkit-box',
    webkitBoxOrient: 'vertical',
    webkitLineClamp: String(lines),
    overflow: 'hidden',
    overflowWrap: 'anywhere',
    whiteSpace: 'normal',
    minWidth: '0',
  });
  el.style.setProperty('-webkit-line-clamp', String(lines));
  el.style.setProperty('-webkit-box-orient', 'vertical');
}

/** `@molstar/ui`'s card, mounted in the canvas's element: one, reused. */
export function hoverCardOverlay(options: { readonly maxWidth?: number } = {}): HoverCardPortFactory {
  return (container) => {
    if (!container) return { show() {}, hide() {} };
    const doc = container.ownerDocument;
    const card = doc.createElement('div');
    card.className = 'msui-surface seqx-hover-card';
    card.setAttribute('aria-hidden', 'true');
    Object.assign(card.style, {
      position: 'absolute',
      display: 'none',
      left: '0',
      top: '0',
      zIndex: 'var(--msui-z-tooltip)',
      pointerEvents: 'none',
      maxWidth: `${options.maxWidth ?? CardMaxPx}px`,
      boxSizing: 'border-box',
      boxShadow: '0 0 0 1px var(--msui-line), var(--msui-shadow-float)',
    });
    // compact: a card is read at a glance (its title a line, its rows a text line each, the value close to its label)
    const title = doc.createElement('div');
    title.className = 'msui-header';
    Object.assign(title.style, { height: 'auto', paddingBlock: 'var(--msui-space-1)', minWidth: '0' });
    const rows = doc.createElement('div');
    rows.className = 'msui-form';
    // the label column as wide as its longest label, the value column the rest of the capped width (never past it)
    rows.style.setProperty('--msui-label', 'max-content');
    Object.assign(rows.style, {
      gridTemplateColumns: 'max-content minmax(0, 1fr)',
      gridAutoRows: 'auto',
      columnGap: 'var(--msui-space-2)',
      paddingBlockEnd: 'var(--msui-space-1)',
    });
    card.append(title, rows);
    const pool: [HTMLElement, HTMLElement][] = [];
    container.append(card);
    let shown: HoverCardShowOptions | undefined;
    const enter = () => shown?.enter?.();
    const leave = () => shown?.leave?.();
    card.addEventListener('pointerenter', enter);
    card.addEventListener('pointerleave', leave);
    return {
      show(c, place, opts) {
        shown = opts;
        const lines = opts?.lines ?? 3;
        title.textContent = breakable(c.title);
        clamp(title, 2);
        const link = (href: string) => opts?.link?.(href);
        // the rows' elements are kept across shows: made when a card has more rows than any before, hidden past its own
        for (let i = 0; i < c.rows.length; i++) {
          let pair = pool[i];
          if (!pair) {
            const label = doc.createElement('span');
            label.className = 'msui-label';
            const value = doc.createElement('span');
            value.className = 'msui-value msui-num';
            Object.assign(value.style, { textAlign: 'start', paddingInlineStart: '0' });
            rows.append(label, value);
            pool.push((pair = [label, value]));
          }
          const r = c.rows[i];
          pair[0].textContent = '';
          pair[1].textContent = '';
          // a row with no value is its text alone (a tooltip with no label): it spans the row and wraps as a value
          // does (the label column is as wide as its longest label: a long text there would never wrap)
          const alone = !r.value;
          if (alone && r.rich) drawRich(doc, pair[0], r.rich, link);
          else pair[0].textContent = alone ? breakable(r.label) : r.label;
          if (!alone) {
            if (r.rich) drawRich(doc, pair[1], r.rich, link);
            else pair[1].textContent = breakable(r.value);
          }
          pair[0].style.gridColumn = alone ? '1 / -1' : '';
          // shown as clamped boxes (`clamp`'s display: a plain display would drop the line clamp)
          clamp(pair[0], lines);
          clamp(pair[1], lines);
          if (alone) pair[1].style.display = 'none';
        }
        for (let i = c.rows.length; i < pool.length; i++) pool[i][0].style.display = pool[i][1].style.display = 'none';
        // `display`, not `hidden`: the surface's own `display` would win over the attribute's
        rows.style.display = c.rows.length ? '' : 'none';
        // a card with a link takes the pointer (entered across the gap); one without never steals a hover from the data
        card.style.pointerEvents = cardHasLink(c) ? 'auto' : 'none';
        card.style.display = '';
        // placed at its own size (jsdom measures none: the estimate's place then)
        const width = card.offsetWidth;
        const height = card.offsetHeight;
        const at = width && height ? placeCard(place.anchor, place.pointer, place.view, { width, height }) : place;
        card.style.transform = `translate(${Math.round(at.x)}px, ${Math.round(at.y)}px)`;
        card.dataset.side = at.side;
      },
      hide() {
        shown = undefined;
        card.style.display = 'none';
        card.style.pointerEvents = 'none';
      },
      dispose() {
        card.removeEventListener('pointerenter', enter);
        card.removeEventListener('pointerleave', leave);
        card.remove();
      },
    };
  };
}
