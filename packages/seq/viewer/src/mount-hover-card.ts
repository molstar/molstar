/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type HoverCardPortFactory, placeCard } from './hover-card.js';

// The hover card's default adapter (_spec/README.md "Ports", Hover card): `@molstar/ui`'s primitives over the canvas: a
// surface (`msui-surface`: its tokens, its squircle where the browser draws one, its shadow as what floats), a header
// with the title and a form of rows (`msui-form`, `msui-label`, `msui-value`). DOM over the canvas: never in the draw
// lists, never in an export. Hidden from assistive technology: the live region says the same words (Announcements).
// Placed by `placeCard` at its measured size (the place's suggestion is for an estimated size).

/** `@molstar/ui`'s card, mounted in the canvas's element: one, reused. */
export function hoverCardOverlay(): HoverCardPortFactory {
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
      maxWidth: '320px',
      boxShadow: '0 0 0 1px var(--msui-line), var(--msui-shadow-float)',
    });
    const title = doc.createElement('div');
    title.className = 'msui-header';
    const rows = doc.createElement('div');
    rows.className = 'msui-form';
    // the label column as wide as its longest label: a card is narrow
    rows.style.setProperty('--msui-label', 'max-content');
    card.append(title, rows);
    const pool: [HTMLElement, HTMLElement][] = [];
    container.append(card);
    return {
      show(c, place) {
        title.textContent = c.title;
        // the rows' elements are kept across shows: made when a card has more rows than any before, hidden past its own
        for (let i = 0; i < c.rows.length; i++) {
          let pair = pool[i];
          if (!pair) {
            const label = doc.createElement('span');
            label.className = 'msui-label';
            const value = doc.createElement('span');
            value.className = 'msui-value msui-num';
            value.style.textAlign = 'start';
            rows.append(label, value);
            pool.push((pair = [label, value]));
          }
          pair[0].textContent = c.rows[i].label;
          pair[1].textContent = c.rows[i].value;
          pair[0].style.display = pair[1].style.display = '';
        }
        for (let i = c.rows.length; i < pool.length; i++) pool[i][0].style.display = pool[i][1].style.display = 'none';
        // `display`, not `hidden`: the surface's own `display` would win over the attribute's
        rows.style.display = c.rows.length ? '' : 'none';
        card.style.display = '';
        // placed at its own size (jsdom measures none: the estimate's place then)
        const width = card.offsetWidth;
        const height = card.offsetHeight;
        const at = width && height ? placeCard(place.anchor, place.pointer, place.view, { width, height }) : place;
        card.style.transform = `translate(${Math.round(at.x)}px, ${Math.round(at.y)}px)`;
        card.dataset.side = at.side;
      },
      hide() {
        card.style.display = 'none';
      },
      dispose() {
        card.remove();
      },
    };
  };
}
