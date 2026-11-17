/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { on } from '../../dom/index.js';
import { classes } from './_generated/wordmark.module.css.js';

// The Mol* wordmark as outlines, after the logo, drawn on the pixel grid at its size (20 px, with a 1 px margin for the
// line, which a heading clips): every straight edge sits mid-pixel, so the line is the same on all of them, every stroke
// of a letter is 3 px across (the M's diagonals too, measured square to them: their edges 6.18 apart vertically). The M
// symmetric, its diagonals parallel; the o a rounded rectangle; the l running into an asterisk whose arms are exactly
// 45° (built as the Mol* mark is, icon/molecular.ts).
const M =
  'M0.5 1.5 3.5 1.5 8.5 10.5 13.5 1.5 16.5 1.5 16.5 19.5 13.5 19.5 13.5 7.68 8.5 16.68 3.5 7.68 3.5 19.5 0.5 19.5Z';
const LStar =
  'M46.5 3.5 49.5 3.5 49.5 6.38 51.54 4.34 53.66 6.46 51.62 8.5 54.5 8.5 54.5 11.5 51.62 11.5 53.66 13.54 51.54 15.66 49.5 13.62 49.5 16.5 46.5 16.5 46.5 13.62 44.46 15.66 42.34 13.54 44.38 11.5 36.5 11.5 36.5 19.5 33.5 19.5 33.5 0.5 36.5 0.5 36.5 8.5 44.38 8.5 42.34 6.46 44.46 4.34 46.5 6.38Z';

/** Drawn width with its margin (the viewBox runs from -1): to the asterisk's last arm (54.5) and the margin. */
const Width = 57;
const Height = 22;
/** The asterisk's middle (`data-part="star"`): where an effect drawn from the word starts. */
const StarX = 48;
const StarY = 10;

/** The o's two rounded rectangles as path data (the outline draws them as `rect`s), for drawing the word elsewhere. */
const OOuter =
  'M24.5 6.5H25.5A5 5 0 0 1 30.5 11.5V14.5A5 5 0 0 1 25.5 19.5H24.5A5 5 0 0 1 19.5 14.5V11.5A5 5 0 0 1 24.5 6.5Z';
const OInner = 'M25 9.5A2.5 2.5 0 0 1 27.5 12V14A2.5 2.5 0 0 1 25 16.5A2.5 2.5 0 0 1 22.5 14V12A2.5 2.5 0 0 1 25 9.5Z';

/**
 * The word as geometry, in drawing units (the viewBox's: x and y from -1, 1 unit a pixel at 1×), for drawing it
 * elsewhere as the outline draws it (a canvas effect). The line is `stroke(scale)` units wide, as the outline's CSS has
 * it.
 */
export const WordmarkGeometry = {
  paths: [M, OOuter, OInner, LStar] as readonly string[],
  width: Width,
  height: Height,
  star: { x: StarX, y: StarY },
  stroke: (scale: number) => 1.5 / Math.sqrt(scale),
} as const;

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

/** A gradient stop in one of the logo's colors (`--msui-logo-1` … `-4`, wordmark.module.css). */
function Stop(props: { offset: string; color: 1 | 2 | 3 | 4 }) {
  return <stop offset={props.offset} style={{ 'stop-color': `var(--msui-logo-${props.color})` }} />;
}

/** The wordmarks' own text (`MolstarWordmark`, the glass's; `DefaultWordmarkLabels`). */
export interface WordmarkLabels {
  /** What the image says (its accessible name): the word it draws. */
  mark: string;
}

/** The word the wordmark draws, as its accessible name (the image is that word, so its name is too). */
export const DefaultWordmarkLabels: WordmarkLabels = { mark: 'Mol*' };

/**
 * The Mol* wordmark, outlined in the text color (a heading, an empty canvas). A flash of the logo's colors runs along
 * the outline now and then (often, or wiped in once and kept) and on hover: a copy stroked with the logo's gradient,
 * shown through a moving band (CSS only).
 */
export function MolstarWordmark(props: {
  class?: string;
  /** A whole multiple of its drawn size (an empty canvas shows it at 3): lines stay on the pixel grid. */
  scale?: number;
  /** Its colors: a flash now and then (a heading), often, or wiped in once and kept (an empty canvas). */
  flash?: 'rare' | 'often' | 'reveal';
  /**
   * Alive under the pointer (an empty canvas): a vertical strip at the cursor, the logo's colors running across it
   * centred on the cursor, moving with it, fading in on enter and back out on leave.
   */
  live?: boolean;
  labels?: Partial<WordmarkLabels>;
}) {
  const id = `msui-wordmark-${++gradients}`;
  const scale = Math.max(1, Math.round(props.scale ?? 1));
  const svg = (extra: { class: string; stroke: string; hidden?: boolean }) => (
    <svg
      class={extra.class}
      viewBox={`-1 -1 ${Width} ${Height}`}
      width={Width * scale}
      height={Height * scale}
      fill="none"
      stroke={extra.stroke}
      stroke-linejoin="miter"
      aria-hidden={extra.hidden ? 'true' : undefined}
    >
      {extra.hidden ? (
        <defs>
          {/* The logo's colors, left to right. */}
          <linearGradient id={id} gradientUnits="userSpaceOnUse" x1="0" y1="0" x2={String(Width - 2)} y2="0">
            <Stop offset="0" color={1} />
            <Stop offset="0.35" color={2} />
            <Stop offset="0.6" color={3} />
            <Stop offset="1" color={4} />
          </linearGradient>
        </defs>
      ) : null}
      <Shapes />
      {/* The star's middle, unpainted: what an effect drawn from the word measures. */}
      {extra.hidden ? null : (
        <rect data-part="star" x={StarX - 0.5} y={StarY - 0.5} width="1" height="1" stroke="none" visibility="hidden" />
      )}
    </svg>
  );
  const glow = props.live ? live(id, scale) : undefined;
  const root = (
    <span
      class={`${classes.wordmark} ${props.class ?? ''}`}
      style={{ '--msui-wordmark-scale': scale }}
      data-flash={props.flash ?? 'rare'}
      data-live={props.live ? '' : undefined}
      role="img"
      aria-label={props.labels?.mark ?? DefaultWordmarkLabels.mark}
    >
      {svg({ class: classes.outline, stroke: 'currentColor' })}
      {svg({ class: classes.flash, stroke: `url(#${id})`, hidden: true })}
      {glow?.svg}
    </span>
  ) as HTMLElement;
  if (glow) {
    // Follow the pointer across: the colors' centre (drawing units) and the shown strip (px); nothing allocated per move.
    on<PointerEvent>(root, 'pointermove', (e) => {
      const r = root.getBoundingClientRect();
      const x = e.clientX - r.left;
      const at = x / scale - 1;
      glow.gradient.setAttribute('x1', String(at - 14));
      glow.gradient.setAttribute('x2', String(at + 14));
      root.style.setProperty('--msui-wordmark-x', `${x}px`);
    });
  }
  return root;
}

/** The live layer: the outline stroked with the logo's colors running across, centred where the pointer is. */
function live(id: string, scale: number) {
  const gradient = (
    <linearGradient id={`${id}-live`} gradientUnits="userSpaceOnUse" x1="13" y1="0" x2="41" y2="0">
      <Stop offset="0" color={1} />
      <Stop offset="0.3" color={3} />
      <Stop offset="0.5" color={2} />
      <Stop offset="0.7" color={3} />
      <Stop offset="1" color={4} />
    </linearGradient>
  ) as SVGLinearGradientElement;
  const svg = (
    <svg
      class={classes.live}
      viewBox={`-1 -1 ${Width} ${Height}`}
      width={Width * scale}
      height={Height * scale}
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
