/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from './_generated/wordmark.module.css.js';

// The Mol* wordmark as outlines, after the logo, drawn on the pixel grid at its size (20 px, with a 1 px margin for the
// line, which a heading clips): every straight edge sits mid-pixel, so the line is the same on all of them, every stroke
// of a letter is 3 px across (the M's diagonals too, measured square to them: their edges 6.18 apart vertically). The M
// symmetric, its diagonals parallel; the o a rounded rectangle; the l running into an asterisk whose arms are exactly 45° (built as the Mol* mark is, molecular.ts).
// Then Prism's mark: a triangle as tall as the asterisk, pointing away from it (a prism, and a play button: stories,
// animation), its left face on the pixel grid; inside it a smaller one in a thicker line, a core (a tetrahedron seen
// from above).
const M =
  'M0.5 1.5 3.5 1.5 8.5 10.5 13.5 1.5 16.5 1.5 16.5 19.5 13.5 19.5 13.5 7.68 8.5 16.68 3.5 7.68 3.5 19.5 0.5 19.5Z';
const LStar =
  'M46.5 3.5 49.5 3.5 49.5 6.38 51.54 4.34 53.66 6.46 51.62 8.5 54.5 8.5 54.5 11.5 51.62 11.5 53.66 13.54 51.54 15.66 49.5 13.62 49.5 16.5 46.5 16.5 46.5 13.62 44.46 15.66 42.34 13.54 44.38 11.5 36.5 11.5 36.5 19.5 33.5 19.5 33.5 0.5 36.5 0.5 36.5 8.5 44.38 8.5 42.34 6.46 44.46 4.34 46.5 6.38Z';

// The triangle: 3.5 to 16.5 like the asterisk, equilateral (its depth 13·√3/2); the core inset 2.5 square to each side.
const Prism = 'M57.5 3.5 68.76 10 57.5 16.5Z';
const PrismCore = 'M60 7.83 63.76 10 60 12.17Z';

/** Drawn width with its margin (the viewBox runs from -1). */
const Width = 72;
/** Without Prism's mark: the Mol* wordmark alone, to the asterisk's last arm (54.5) and the margin. */
const WordWidth = 57;
/** The asterisk's middle: where the empty canvas's splash splits the light (`data-part="star"`). */
const StarX = 48;
const StarY = 10;

let gradients = 0;

/** The outline's shapes. */
function Shapes(props: { prism: boolean }) {
  return (
    <>
      <path d={M} />
      <rect x="19.5" y="6.5" width="11" height="13" rx="5" />
      <rect x="22.5" y="9.5" width="5" height="7" rx="2.5" />
      <path d={LStar} />
      {props.prism ? <path d={Prism} data-part="prism" /> : null}
      {props.prism ? <path d={PrismCore} class={classes.core} /> : null}
    </>
  );
}

/**
 * The Mol* wordmark with Prism's mark after it, outlined in the text color: Home's heading and the empty canvas. A flash of the logo's colors runs along the outline
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
   * Alive under the pointer (the empty canvas): a vertical strip at the cursor, the logo's colors running across it
   * centred on the cursor, moving with it, fading in on enter and back out on leave.
   */
  live?: boolean;
  /** Prism's mark after the word (on by default); the empty canvas leaves it out (its splash makes the star the prism). */
  prism?: boolean;
}) {
  const id = `msp-wordmark-${++gradients}`;
  const scale = Math.max(1, Math.round(props.scale ?? 1));
  const prism = props.prism ?? true;
  const width = prism ? Width : WordWidth;
  const svg = (extra: { class: string; stroke: string; hidden?: boolean }) => (
    <svg
      class={extra.class}
      viewBox={`-1 -1 ${width} 22`}
      width={width * scale}
      height={22 * scale}
      fill="none"
      stroke={extra.stroke}
      stroke-linejoin="miter"
      aria-hidden={extra.hidden ? 'true' : undefined}
    >
      {extra.hidden ? (
        <defs>
          {/* The logo's colors, left to right. */}
          <linearGradient id={id} gradientUnits="userSpaceOnUse" x1="0" y1="0" x2={String(width - 2)} y2="0">
            {/* The theme's brand colors (theme.css): lighter on dark backgrounds, from the nearest theme. */}
            <stop offset="0" style={{ 'stop-color': 'var(--msp-color-brand-1)' }} />
            <stop offset="0.35" style={{ 'stop-color': 'var(--msp-color-brand-2)' }} />
            <stop offset="0.6" style={{ 'stop-color': 'var(--msp-color-brand-3)' }} />
            <stop offset="1" style={{ 'stop-color': 'var(--msp-color-brand-4)' }} />
          </linearGradient>
        </defs>
      ) : null}
      <Shapes prism={prism} />
      {/* The star's middle, unpainted: what the splash measures. */}
      {extra.hidden ? null : (
        <rect data-part="star" x={StarX - 0.5} y={StarY - 0.5} width="1" height="1" stroke="none" visibility="hidden" />
      )}
    </svg>
  );
  const glow = props.live ? live(id, scale, width, prism) : undefined;
  const root = (
    <span
      class={`${classes.wordmark} ${props.class ?? ''}`}
      style={{ '--msp-wordmark-scale': scale, '--msp-wordmark-width': width }}
      data-flash={props.flash ?? 'rare'}
      data-live={props.live ? '' : undefined}
      role="img"
      aria-label="Mol* Prism"
    >
      {svg({ class: classes.outline, stroke: 'currentColor' })}
      {svg({ class: classes.flash, stroke: `url(#${id})`, hidden: true })}
      {glow?.svg}
    </span>
  ) as HTMLElement;
  if (glow) {
    // Follow the pointer across: the colors' centre (drawing units) and the shown strip (px); nothing allocated per move.
    root.addEventListener('pointermove', (e) => {
      const r = root.getBoundingClientRect();
      const x = e.clientX - r.left;
      const at = x / scale - 1;
      glow.gradient.setAttribute('x1', String(at - 14));
      glow.gradient.setAttribute('x2', String(at + 14));
      root.style.setProperty('--msp-wordmark-x', `${x}px`);
    });
  }
  return root;
}

/** The live layer: the outline stroked with the logo's colors running across, centred where the pointer is. */
function live(id: string, scale: number, width: number, prism: boolean) {
  const gradient = (
    <linearGradient id={`${id}-live`} gradientUnits="userSpaceOnUse" x1="13" y1="0" x2="41" y2="0">
      <stop offset="0" style={{ 'stop-color': 'var(--msp-color-brand-1)' }} />
      <stop offset="0.3" style={{ 'stop-color': 'var(--msp-color-brand-3)' }} />
      <stop offset="0.5" style={{ 'stop-color': 'var(--msp-color-brand-2)' }} />
      <stop offset="0.7" style={{ 'stop-color': 'var(--msp-color-brand-3)' }} />
      <stop offset="1" style={{ 'stop-color': 'var(--msp-color-brand-4)' }} />
    </linearGradient>
  ) as SVGLinearGradientElement;
  const svg = (
    <svg
      class={classes.live}
      viewBox={`-1 -1 ${width} 22`}
      width={width * scale}
      height={22 * scale}
      fill="none"
      stroke={`url(#${id}-live)`}
      stroke-linejoin="miter"
      aria-hidden="true"
    >
      <defs>{gradient}</defs>
      <Shapes prism={prism} />
    </svg>
  );
  return { svg, gradient };
}
