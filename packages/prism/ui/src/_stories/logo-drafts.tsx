/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Drafts for the Prism logo (Overview › Foundations › Logo): the Mol* wordmark with a delta after it, a beam from the
// star's arm into it, split into the logo's colors on the way out; and the Home icon as the delta. Not used by the app
// until one is picked.

// The wordmark's letters, as in components/icon/wordmark.tsx (outlines on the pixel grid, 20 units tall).
const M =
  'M0.5 1.5 3.5 1.5 8.5 10.5 13.5 1.5 16.5 1.5 16.5 19.5 13.5 19.5 13.5 7.68 8.5 16.68 3.5 7.68 3.5 19.5 0.5 19.5Z';
const LStar =
  'M46.5 3.5 49.5 3.5 49.5 6.38 51.54 4.34 53.66 6.46 51.62 8.5 54.5 8.5 54.5 11.5 51.62 11.5 53.66 13.54 51.54 15.66 49.5 13.62 49.5 16.5 46.5 16.5 46.5 13.62 44.46 15.66 42.34 13.54 44.38 11.5 36.5 11.5 36.5 19.5 33.5 19.5 33.5 0.5 36.5 0.5 36.5 8.5 44.38 8.5 42.34 6.46 44.46 4.34 46.5 6.38Z';

const Brand = ['--msp-color-brand-1', '--msp-color-brand-2', '--msp-color-brand-3', '--msp-color-brand-4'];

/** The delta as an outlined band (3 units across, like the letters' strokes): an outer and an inner triangle. */
function delta(x0: number) {
  const top = 0.5;
  const base = 19.5;
  const half = (base - top) / Math.sqrt(3);
  const cx = x0 + half;
  // Inset by 3 square to each side: the apex moves down twice that, the base up by it.
  const iTop = top + 6;
  const iBase = base - 3;
  const iHalf = (iBase - iTop) / Math.sqrt(3);
  const f = (n: number) => n.toFixed(2);
  return {
    outer: `M${f(cx)} ${top} ${f(x0 + 2 * half)} ${base} ${f(x0)} ${base}Z`,
    inner: `M${f(cx)} ${iTop} ${f(cx + iHalf)} ${iBase} ${f(cx - iHalf)} ${iBase}Z`,
    /** Where the left and right outer faces are at a height. */
    left: (y: number) => x0 + ((base - y) / (base - top)) * half,
    right: (y: number) => cx + ((y - top) / (base - top)) * half,
    innerLeft: (y: number) => cx - ((y - iTop) / (iBase - iTop)) * iHalf,
    innerRight: (y: number) => cx + ((y - iTop) / (iBase - iTop)) * iHalf,
  };
}

export type WordmarkDraft =
  | 'delta'
  | 'beam'
  | 'spectrum'
  | 'refracted'
  | 'play'
  | 'play-spectrum'
  | 'play-nested'
  | 'play-levels'
  | 'play-tetra'
  | 'play-core';

/**
 * The triangle turned to point away from the star, the same height: a prism, and a play button (animations, stories).
 * Its left face is vertical on the pixel grid; outlined as the letters are (3 across).
 */
function play(x0: number) {
  const top = 0.5;
  const base = 19.5;
  const mid = (top + base) / 2;
  const w = ((base - top) * Math.sqrt(3)) / 2;
  // Inset by 3 square to each side: the left face moves right by 3, the apex left by 6, the corners in along the faces.
  const iLeft = x0 + 3;
  const iApex = x0 + w - 6;
  const iHalf = (iApex - iLeft) / Math.sqrt(3);
  const f = (n: number) => n.toFixed(2);
  return {
    apex: x0 + w,
    outer: `M${x0} ${top} ${f(x0 + w)} ${mid} ${x0} ${base}Z`,
    inner: `M${iLeft} ${f(mid - iHalf)} ${f(iApex)} ${mid} ${iLeft} ${f(mid + iHalf)}Z`,
  };
}

/** The wordmark with the delta: alone, with a beam and rays in the text color, or the rays in the logo's colors. */
export function PrismWordmarkDraft(props: { variant: WordmarkDraft; scale?: number }) {
  const scale = props.scale ?? 1;
  const x0 = 57.5;
  if (props.variant === 'play' || props.variant === 'play-spectrum') return playWordmark(props.variant, scale, x0);
  if (
    props.variant === 'play-nested' ||
    props.variant === 'play-levels' ||
    props.variant === 'play-tetra' ||
    props.variant === 'play-core'
  )
    return tetraWordmark(props.variant, scale, x0);
  const d = delta(x0);
  const beamY = 10;
  const exitY = 11;
  const width = props.variant === 'delta' ? x0 + 23 : x0 + 39;
  const rays = [10.5, 13.5, 16.5, 19.5];
  const colored = props.variant === 'spectrum' || props.variant === 'refracted';
  return (
    <svg
      viewBox={`-1 -1 ${width + 1} 22`}
      width={(width + 1) * scale}
      height={22 * scale}
      fill="none"
      stroke="currentColor"
      stroke-linejoin="miter"
      style={{ 'stroke-width': String(1.5 / Math.sqrt(scale)), overflow: 'visible' }}
    >
      <path d={M} />
      <rect x="19.5" y="6.5" width="11" height="13" rx="5" />
      <rect x="22.5" y="9.5" width="5" height="7" rx="2.5" />
      <path d={LStar} />
      <path d={d.outer} />
      <path d={d.inner} />
      {props.variant === 'delta' ? null : (
        <>
          {/* The beam: from the star's right arm to the delta's left face. */}
          <path d={`M54.5 ${beamY}H${d.left(beamY).toFixed(2)}`} />
          {props.variant === 'refracted' ? (
            <path
              d={`M${d.innerLeft(beamY + 0.6).toFixed(2)} ${beamY + 0.6}L${d.innerRight(exitY - 0.4).toFixed(2)} ${exitY - 0.4}`}
              style={{ opacity: '0.6' }}
            />
          ) : null}
          {rays.map((y, i) => (
            <path
              d={`M${d.right(exitY).toFixed(2)} ${exitY}L${(x0 + 38).toFixed(2)} ${y}`}
              style={colored ? { stroke: `var(${Brand[i]})` } : undefined}
            />
          ))}
        </>
      )}
    </svg>
  );
}

/** The pointing-away triangle inset by `d` square to each side (it shrinks toward its centroid). */
function playInset(x0: number, d: number) {
  const top = 0.5;
  const base = 19.5;
  const mid = (top + base) / 2;
  const w = ((base - top) * Math.sqrt(3)) / 2;
  const left = x0 + d;
  const apex = x0 + w - 2 * d;
  const half = (apex - left) / Math.sqrt(3);
  const f = (n: number) => n.toFixed(2);
  return { path: `M${f(left)} ${f(mid - half)} ${f(apex)} ${mid} ${f(left)} ${f(mid + half)}Z`, left, apex, mid, half };
}

/**
 * Seen from above, a tetrahedron: triangles inside the triangle, toward its centre, each line thinner (`nested`: the
 * letters' band and two more inside it; `levels`: evenly spaced, no band; `tetra`: the levels and the three edges up
 * to the apex; `core`: the levels thickest at the centre, thinning out).
 */
function tetraWordmark(variant: 'play-nested' | 'play-levels' | 'play-tetra' | 'play-core', scale: number, x0: number) {
  const base = 1.5 / Math.sqrt(scale);
  const levels: [number, number][] =
    variant === 'play-core'
      ? // Thickest at the centre, thinning out to the outline.
        [
          [0, 0.45],
          [1.6, 0.6],
          [3, 0.8],
          [4.15, 1],
          [5, 1.25],
        ]
      : variant === 'play-nested'
        ? [
            [0, 1],
            [3, 1],
            [4.1, 0.6],
            [4.9, 0.35],
          ]
        : [
            [0, 1],
            [1.6, 0.75],
            [3, 0.55],
            [4.15, 0.4],
            [5, 0.25],
          ];
  const outer = playInset(x0, 0);
  const width = outer.apex + 1;
  const cx = x0 + (outer.apex - x0) / 3;
  return (
    <svg
      viewBox={`-1 -1 ${width + 1} 22`}
      width={(width + 1) * scale}
      height={22 * scale}
      fill="none"
      stroke="currentColor"
      stroke-linejoin="miter"
      style={{ 'stroke-width': String(base), overflow: 'visible' }}
    >
      <path d={M} />
      <rect x="19.5" y="6.5" width="11" height="13" rx="5" />
      <rect x="22.5" y="9.5" width="5" height="7" rx="2.5" />
      <path d={LStar} />
      {levels.map(([d, k]) => (
        <path d={playInset(x0, d).path} style={{ 'stroke-width': String(base * k) }} />
      ))}
      {variant === 'play-tetra' ? (
        <path
          d={`M${x0} ${outer.mid - outer.half}L${cx.toFixed(2)} ${outer.mid}M${outer.apex.toFixed(2)} ${outer.mid}L${cx.toFixed(2)} ${outer.mid}M${x0} ${outer.mid + outer.half}L${cx.toFixed(2)} ${outer.mid}`}
          style={{ 'stroke-width': String(base * 0.4) }}
        />
      ) : null}
    </svg>
  );
}

function playWordmark(variant: 'play' | 'play-spectrum', scale: number, x0: number) {
  const p = play(x0);
  const rays = variant === 'play-spectrum';
  const end = p.apex + 13;
  const width = rays ? end + 1 : p.apex + 1;
  const ends = [3.5, 8, 12, 16.5];
  return (
    <svg
      viewBox={`-1 -1 ${width + 1} 22`}
      width={(width + 1) * scale}
      height={22 * scale}
      fill="none"
      stroke="currentColor"
      stroke-linejoin="miter"
      style={{ 'stroke-width': String(1.5 / Math.sqrt(scale)), overflow: 'visible' }}
    >
      <path d={M} />
      <rect x="19.5" y="6.5" width="11" height="13" rx="5" />
      <rect x="22.5" y="9.5" width="5" height="7" rx="2.5" />
      <path d={LStar} />
      <path d={p.outer} />
      <path d={p.inner} />
      {rays
        ? ends.map((y, i) => (
            <path
              d={`M${(p.apex + 2).toFixed(2)} ${10 + (y - 10) * 0.12}L${end.toFixed(2)} ${y}`}
              style={{ stroke: `var(${Brand[i]})` }}
            />
          ))
        : null}
    </svg>
  );
}

export type IconDraft = 'delta' | 'rays' | 'spectrum' | 'beam' | 'play' | 'play-spectrum';

/**
 * The Home icon on the 24 grid, in Lucide's stroke (2, round): the delta alone, with three rays fanning out of its right
 * face (in the text color or the logo's), or with a beam coming in as well.
 */
export function PrismIconDraft(props: { variant: IconDraft; size: number }) {
  const rays = ['M14 9.5 21.5 7', 'M14.8 12 22 12.5', 'M15.6 14.5 21 18'];
  const colors = ['--msp-color-brand-1', '--msp-color-brand-3', '--msp-color-brand-4'];
  return (
    <svg
      viewBox="0 0 24 24"
      width={props.size}
      height={props.size}
      fill="none"
      stroke="currentColor"
      stroke-width="2"
      stroke-linecap="round"
      stroke-linejoin="round"
    >
      {props.variant === 'play' ? (
        <path d="M5 3.5 19.5 12 5 20.5z" />
      ) : props.variant === 'play-spectrum' ? (
        <>
          <path d="M3 4.5 15.5 12 3 19.5z" />
          <path d="M18 10 22 7.5" style={{ stroke: `var(${colors[0]})` }} />
          <path d="M18.5 12H22.5" style={{ stroke: `var(${colors[1]})` }} />
          <path d="M18 14 22 16.5" style={{ stroke: `var(${colors[2]})` }} />
        </>
      ) : props.variant === 'delta' ? (
        <path d="M12 3.5 21 19.5H3z" />
      ) : props.variant === 'beam' ? (
        <>
          <path d="M12 5 18.5 18.5H5.5z" />
          <path d="M1.5 12.5H8.6" />
          <path d="M15.4 12.5 22.5 10M15.4 12.5 22.5 15" />
        </>
      ) : (
        <>
          <path d="M9.5 4 16.5 19H2.5z" />
          {rays.map((d, i) => (
            <path d={d} style={props.variant === 'spectrum' ? { stroke: `var(${colors[i]})` } : undefined} />
          ))}
        </>
      )}
    </svg>
  );
}
