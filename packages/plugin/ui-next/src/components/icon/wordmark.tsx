/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from './_generated/wordmark.module.css.js';

// The Mol* wordmark as outlines, after the logo, drawn on the pixel grid at its size (20 px, with a 1 px margin for the
// line, which a heading clips): every straight edge sits mid-pixel, so the line is the same on all of them, every stroke
// of a letter is 3 px across (the M's diagonals too, measured square to them: their edges 6.18 apart vertically). The M
// symmetric, its diagonals parallel; the o a rounded rectangle; the l running into an asterisk whose arms are exactly 45° (built as the Mol* mark is, molecular.ts).
const M =
  'M0.5 1.5 3.5 1.5 8.5 10.5 13.5 1.5 16.5 1.5 16.5 19.5 13.5 19.5 13.5 7.68 8.5 16.68 3.5 7.68 3.5 19.5 0.5 19.5Z';
const LStar =
  'M46.5 3.5 49.5 3.5 49.5 6.38 51.54 4.34 53.66 6.46 51.62 8.5 54.5 8.5 54.5 11.5 51.62 11.5 53.66 13.54 51.54 15.66 49.5 13.62 49.5 16.5 46.5 16.5 46.5 13.62 44.46 15.66 42.34 13.54 44.38 11.5 36.5 11.5 36.5 19.5 33.5 19.5 33.5 0.5 36.5 0.5 36.5 8.5 44.38 8.5 42.34 6.46 44.46 4.34 46.5 6.38Z';

let gradients = 0;

/** The outline's shapes. */
function Shapes() {
  return (
    <>
      <path d={M} />
      <rect x="19.5" y="6.5" width="11" height="13" rx="5" />
      <rect x="22.5" y="9.5" width="5" height="7" rx="2.5" />
      <path d={LStar} />
    </>
  );
}

/**
 * The Mol* wordmark, outlined in the text color: Home's heading. A flash of the logo's colors runs along the outline
 * now and then (often on the empty canvas) and on hover (a copy stroked with the logo's gradient, shown through a
 * moving band: CSS only).
 */
export function MolstarWordmark(props: {
  class?: string;
  /** A whole multiple of its drawn size (the empty canvas shows it at 3): lines stay on the pixel grid. */
  scale?: number;
  /** Its colors: a flash now and then (a heading), often, or wiped in once and kept (the empty canvas). */
  flash?: 'rare' | 'often' | 'reveal';
  /**
   * Alive under the pointer (the empty canvas): around the cursor the logo's colors in rings centred on it, moving with
   * it, fading in on enter and back out on leave.
   */
  live?: boolean;
}) {
  const id = `msp-wordmark-${++gradients}`;
  const scale = Math.max(1, Math.round(props.scale ?? 1));
  const svg = (extra: { class: string; stroke: string; hidden?: boolean }) => (
    <svg
      class={extra.class}
      viewBox="-1 -1 57 22"
      width={57 * scale}
      height={22 * scale}
      fill="none"
      stroke={extra.stroke}
      stroke-linejoin="miter"
      aria-hidden={extra.hidden ? 'true' : undefined}
    >
      {extra.hidden ? (
        <defs>
          {/* The logo's colors, left to right. */}
          <linearGradient id={id} gradientUnits="userSpaceOnUse" x1="0" y1="0" x2="55" y2="0">
            {/* The theme's brand colors (theme.css): lighter on dark backgrounds, from the nearest theme. */}
            <stop offset="0" style={{ 'stop-color': 'var(--msp-color-brand-1)' }} />
            <stop offset="0.35" style={{ 'stop-color': 'var(--msp-color-brand-2)' }} />
            <stop offset="0.6" style={{ 'stop-color': 'var(--msp-color-brand-3)' }} />
            <stop offset="1" style={{ 'stop-color': 'var(--msp-color-brand-4)' }} />
          </linearGradient>
        </defs>
      ) : null}
      <Shapes />
    </svg>
  );
  const glow = props.live ? live(id, scale) : undefined;
  const root = (
    <span
      class={`${classes.wordmark} ${props.class ?? ''}`}
      style={{ '--msp-wordmark-scale': scale }}
      data-flash={props.flash ?? 'rare'}
      data-live={props.live ? '' : undefined}
      role="img"
      aria-label="Mol*"
    >
      {svg({ class: classes.outline, stroke: 'currentColor' })}
      {svg({ class: classes.flash, stroke: `url(#${id})`, hidden: true })}
      {glow?.svg}
    </span>
  ) as HTMLElement;
  if (glow) {
    // Follow the pointer: the rings' centre (drawing units) and the shown circle (px); nothing allocated per move.
    root.addEventListener('pointermove', (e) => {
      const r = root.getBoundingClientRect();
      const x = e.clientX - r.left;
      const y = e.clientY - r.top;
      glow.gradient.setAttribute('cx', String(x / scale - 1));
      glow.gradient.setAttribute('cy', String(y / scale - 1));
      root.style.setProperty('--msp-wordmark-x', `${x}px`);
      root.style.setProperty('--msp-wordmark-y', `${y}px`);
    });
  }
  return root;
}

/** The live layer: the outline stroked with the logo's colors in rings around a centre the pointer moves. */
function live(id: string, scale: number) {
  const gradient = (
    <radialGradient id={`${id}-live`} gradientUnits="userSpaceOnUse" cx="27" cy="10" r="16">
      <stop offset="0" style={{ 'stop-color': 'var(--msp-color-brand-2)' }} />
      <stop offset="0.35" style={{ 'stop-color': 'var(--msp-color-brand-3)' }} />
      <stop offset="0.7" style={{ 'stop-color': 'var(--msp-color-brand-4)' }} />
      <stop offset="1" style={{ 'stop-color': 'var(--msp-color-brand-1)' }} />
    </radialGradient>
  ) as SVGRadialGradientElement;
  const svg = (
    <svg
      class={classes.live}
      viewBox="-1 -1 57 22"
      width={57 * scale}
      height={22 * scale}
      fill="none"
      stroke={`url(#${id}-live)`}
      stroke-linejoin="miter"
      aria-hidden="true"
    >
      <defs>{gradient}</defs>
      <Shapes />
    </svg>
  );
  return { svg, gradient };
}
