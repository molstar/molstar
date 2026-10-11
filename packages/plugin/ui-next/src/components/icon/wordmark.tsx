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
 * when it appears and on hover (a copy stroked with the logo's gradient, shown through a moving band: CSS only).
 */
export function MolstarWordmark(props: {
  class?: string;
  /** A whole multiple of its drawn size (the empty canvas shows it at 4): lines stay on the pixel grid. */
  scale?: number;
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
            <stop offset="0" stop-color="#7a4f2f" />
            <stop offset="0.35" stop-color="#e8b93a" />
            <stop offset="0.6" stop-color="#7fa38a" />
            <stop offset="1" stop-color="#3f86ad" />
          </linearGradient>
        </defs>
      ) : null}
      <Shapes />
    </svg>
  );
  return (
    <span
      class={`${classes.wordmark} ${props.class ?? ''}`}
      style={{ '--msp-wordmark-scale': scale }}
      role="img"
      aria-label="Mol*"
    >
      {svg({ class: classes.outline, stroke: 'currentColor' })}
      {svg({ class: classes.flash, stroke: `url(#${id})`, hidden: true })}
    </span>
  );
}
