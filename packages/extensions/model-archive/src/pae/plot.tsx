/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, onCleanup, type ReadonlyState, Show } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { bind } from '@molstar/prism-ui/dom/scope';
import { drag, frames, on } from '@molstar/prism-ui/dom/scoped';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import {
  type PaeMatrix,
  PaeColorStyles,
  paeBoxMean,
  paeCellAt,
  paeCellLabel,
  paeChains,
  paeImage,
  paeRangeLabel,
  paeTicks,
} from './data.js';
import { PlddtBands, plddtColor } from '../quality-assessment/plddt-bands.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// The PAE plot (parity.md §6; AlphaFold DB's, in Prism's look): the heat map, scored residue across and aligned residue
// down, one pixel per cell scaled crisply; over it an SVG in the matrix's coordinates for the chains' boundaries, the
// crosshair of the cell under the pointer, and a box dragged over it with its two ranges on the diagonal (scored in
// orange, aligned in yellow, as the 3D view and the sequence mark them) and its edges projected, dashed, to them.
// Residue ticks on both axes, chain bands for a complex, a pLDDT strip along both (its toggle hides it), the scale in
// Å. Pointer-rate input is handled once per frame (styling.md §5).

/** A box of cells: aligned rows `i0`–`i1`, scored columns `j0`–`j1` (all in, ordered). */
export interface PaeBox {
  readonly i0: number;
  readonly i1: number;
  readonly j0: number;
  readonly j1: number;
}

/**
 * The marks' colors (marks.md tones, AlphaFold DB's roles): the aligned range in the highlight tone (yellow), the scored
 * range in the accent; the viewer's, so the plot, the 3D view and the sequence agree. These until it says.
 */
export interface PaeMarkColors {
  readonly highlight: number;
  readonly accent: number;
}
const Fallback: PaeMarkColors = { highlight: 0xf2b400, accent: 0xae5d04 };
const hex = (c: number) => `#${c.toString(16).padStart(6, '0')}`;

export interface PaePlotProps {
  readonly matrix: PaeMatrix;
  /** The cell under the pointer (aligned `i`, scored `j`); undefined: none. */
  onHover?(cell: { i: number; j: number } | undefined): void;
  /** A click on a cell (aligned `i`, scored `j`), while no box is shown. */
  onPick?(i: number, j: number): void;
  /** A box dragged over the plot; undefined: cleared. */
  onBox?(box: PaeBox | undefined): void;
  /** The marks' colors (the viewer's tones). */
  colors?: ReadonlyState<PaeMarkColors>;
}

const Quiet = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs)';
const Tick =
  'position: absolute; color: var(--msp-color-text-subtle); font-size: 10px; line-height: 1; white-space: nowrap';
const Svg = 'http://www.w3.org/2000/svg';
/** A dark rim under the box's white lines: they read on the plot's light cells (a complex's empty blocks) too. */
const Halo = 'filter: drop-shadow(0 0 1px rgba(0, 0, 0, 0.55))';

function svg<K extends keyof SVGElementTagNameMap>(
  tag: K,
  attrs: Record<string, string | number>,
): SVGElementTagNameMap[K] {
  const el = document.createElementNS(Svg, tag);
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, String(v));
  return el;
}

export function PaePlot(props: PaePlotProps) {
  const m = props.matrix;
  const n = m.size;
  const canvas = (
    <canvas
      width={n}
      height={n}
      style="position: absolute; inset: 0; width: 100%; height: 100%; image-rendering: pixelated"
    />
  ) as HTMLCanvasElement;
  const g = canvas.getContext('2d');
  if (g) {
    const image = g.createImageData(n, n);
    paeImage(m, image.data);
    g.putImageData(image, 0, 0);
  }

  // The overlay, in the matrix's coordinates; strokes keep their width however the plot is scaled.
  const overlay = svg('svg', {
    viewBox: `0 0 ${n} ${n}`,
    preserveAspectRatio: 'none',
    style: 'position: absolute; inset: 0; width: 100%; height: 100%; overflow: visible; pointer-events: none',
  });
  const line = (attrs: Record<string, string | number>) =>
    overlay.appendChild(svg('line', { 'vector-effect': 'non-scaling-stroke', ...attrs }));
  const chains = paeChains(m);
  for (const c of chains.slice(1)) {
    const edge = { stroke: 'rgba(255, 255, 255, 0.85)', 'stroke-width': 1 };
    line({ x1: c.from, y1: 0, x2: c.from, y2: n, ...edge });
    line({ x1: 0, y1: c.from, x2: n, y2: c.from, ...edge });
  }
  // The cell under the pointer: its row and column, and its square.
  const hair = { stroke: 'rgba(255, 255, 255, 0.7)', 'stroke-width': 1, visibility: 'hidden' };
  const across = line({ ...hair });
  const down = line({ ...hair });
  // The box, and its two ranges on the diagonal: the scored a narrower line over the aligned, so where they overlap
  // both show (AlphaFold DB gives the overlap a third color; the marks have two tones).
  const box = overlay.appendChild(
    svg('rect', {
      fill: 'rgba(255, 255, 255, 0.16)',
      stroke: 'rgba(255, 255, 255, 0.95)',
      'stroke-width': 1.5,
      'vector-effect': 'non-scaling-stroke',
      style: Halo,
      visibility: 'hidden',
    }),
  );
  // Its projections (AlphaFold DB's): each edge of the box carried, dashed, to the diagonal where it lies outside the
  // box, so each range's run on the diagonal reads from the box.
  const projections = [0, 1, 2, 3].map(() =>
    line({
      stroke: 'rgba(255, 255, 255, 0.95)',
      'stroke-width': 1,
      'stroke-dasharray': '4 3',
      style: Halo,
      visibility: 'hidden',
    }),
  );
  const seg = (color: number, width: number) =>
    line({ stroke: hex(color), 'stroke-width': width, 'stroke-linecap': 'round', visibility: 'hidden' });
  const alignedSeg = seg(Fallback.highlight, 6);
  const scoredSeg = seg(Fallback.accent, 2.5);
  const tones = observableState<PaeMarkColors>(Fallback);
  // Bound, not subscribed: it pauses with the card's scope while hidden (styling.md §5) and goes with it.
  if (props.colors)
    bind(props.colors, (c) => {
      tones.set(c);
      scoredSeg.setAttribute('stroke', hex(c.accent));
      alignedSeg.setAttribute('stroke', hex(c.highlight));
    });

  const plot = (
    <div
      role="img"
      aria-label={`${m.name}: ${n} by ${n} residues; drag a box to compare two ranges`}
      style="position: relative; aspect-ratio: 1; border-radius: var(--msp-radius-sm); overflow: hidden; background: var(--msp-color-well); box-shadow: 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 78%); cursor: crosshair; touch-action: none"
    >
      {canvas}
      {overlay}
    </div>
  ) as HTMLElement;

  // ✕ on the box's top-right corner (AlphaFold DB's), kept inside the plot; beside the plot, not in it (an image).
  const close = (
    <button
      type="button"
      aria-label="Clear the box"
      data-tip="Clear the box"
      style="position: absolute; display: none; place-items: center; width: 20px; height: 20px; padding: 0; border: 1px solid color-mix(in oklab, var(--msp-color-text), transparent 80%); border-radius: 50%; background: var(--msp-color-panel); color: var(--msp-color-text); cursor: pointer; transform: translate(-50%, -50%)"
      onClick={() => clear()}
    >
      <Icon name="x" size={12} />
    </button>
  ) as HTMLButtonElement;
  const placeClose = (b: PaeBox | undefined) => {
    close.style.display = b ? 'grid' : 'none';
    if (!b) return;
    close.style.left = `clamp(10px, ${(100 * (b.j1 + 1)) / n}%, calc(100% - 10px))`;
    close.style.top = `clamp(10px, ${(100 * b.i0) / n}%, calc(100% - 10px))`;
  };
  const plotArea = (
    <div style="position: relative; min-width: 0">
      {plot}
      {close}
    </div>
  );

  // State: the box (committed or being dragged), and the line under the plot.
  const shown = observableState<PaeBox | undefined>(undefined);
  const hoverText = observableState<string | undefined>(undefined);
  const setVisible = (el: SVGElement, on: boolean) => el.setAttribute('visibility', on ? 'visible' : 'hidden');
  const drawBox = (b: PaeBox | undefined) => {
    for (const el of [box, scoredSeg, alignedSeg]) setVisible(el, !!b);
    if (!b) {
      for (const el of projections) setVisible(el, false);
      return;
    }
    box.setAttribute('x', String(b.j0));
    box.setAttribute('y', String(b.i0));
    box.setAttribute('width', String(b.j1 - b.j0 + 1));
    box.setAttribute('height', String(b.i1 - b.i0 + 1));
    const diag = (el: SVGLineElement, a: number, z: number) => {
      el.setAttribute('x1', String(a));
      el.setAttribute('y1', String(a));
      el.setAttribute('x2', String(z + 1));
      el.setAttribute('y2', String(z + 1));
    };
    diag(scoredSeg, b.j0, b.j1);
    diag(alignedSeg, b.i0, b.i1);
    // Row edges run across to the diagonal (x = y), column edges down or up to it (y = x).
    const project = (el: SVGLineElement, from: [number, number], to: [number, number] | undefined) => {
      setVisible(el, !!to);
      if (!to) return;
      el.setAttribute('x1', String(from[0]));
      el.setAttribute('y1', String(from[1]));
      el.setAttribute('x2', String(to[0]));
      el.setAttribute('y2', String(to[1]));
    };
    const [top, bottom, left, right] = projections;
    const across = (el: SVGLineElement, y: number) =>
      project(el, [y, y], y < b.j0 ? [b.j0, y] : y > b.j1 + 1 ? [b.j1 + 1, y] : undefined);
    const upDown = (el: SVGLineElement, x: number) =>
      project(el, [x, x], x < b.i0 ? [x, b.i0] : x > b.i1 + 1 ? [x, b.i1 + 1] : undefined);
    across(top, b.i0);
    across(bottom, b.i1 + 1);
    upDown(left, b.j0);
    upDown(right, b.j1 + 1);
  };
  const drawHair = (cell: { i: number; j: number } | undefined) => {
    setVisible(across, !!cell);
    setVisible(down, !!cell);
    if (!cell) return;
    for (const [el, x1, y1, x2, y2] of [
      [across, 0, cell.i + 0.5, n, cell.i + 0.5],
      [down, cell.j + 0.5, 0, cell.j + 0.5, n],
    ] as const) {
      el.setAttribute('x1', String(x1));
      el.setAttribute('y1', String(y1));
      el.setAttribute('x2', String(x2));
      el.setAttribute('y2', String(y2));
    }
  };

  // Pointer-rate input, once per frame: a hover names and lights the cell; a drag draws a box, or on the box's edge or
  // corner resizes it, or inside it moves it; a press without a drag picks the cell. The box is set when the drag ends.
  let px = 0;
  let py = 0;
  let inside = false;
  let cellKey = -1;
  // A drag hid the crosshair: the first frame after it redraws it (and sends the hover) even on the same cell.
  let hairHidden = false;
  // The cell sent on to the 3D view and the sequence (`onHover`) at 24 a second at most, the last one always: they
  // redraw on each, the plot's own crosshair and line keep the frame rate.
  const HoverEvery = 1000 / 24;
  let hoverSent = 0;
  let hoverTimer: ReturnType<typeof setTimeout> | undefined;
  let hoverNext: { i: number; j: number } | undefined;
  const sendHover = () => {
    hoverTimer = undefined;
    hoverSent = performance.now();
    props.onHover?.(hoverNext);
  };
  const hover = (at: { i: number; j: number } | undefined) => {
    hoverNext = at;
    if (hoverTimer !== undefined) return;
    const wait = HoverEvery - (performance.now() - hoverSent);
    if (wait <= 0) sendHover();
    else hoverTimer = setTimeout(sendHover, wait);
  };
  /** An edge within this many pixels of the pointer is grabbed. */
  const Grab = 6;
  type Edges = { top: boolean; bottom: boolean; left: boolean; right: boolean };
  type Drag = {
    readonly kind: 'new' | 'resize' | 'move';
    readonly i: number;
    readonly j: number;
    readonly x: number;
    readonly y: number;
    readonly from?: PaeBox;
    readonly edges?: Edges;
    moved: boolean;
  };
  let pressed: Drag | undefined;
  // Gone: a drag ended by the scope's end asks for no frame.
  let disposed = false;
  const cellAt = (x: number, y: number) => {
    const r = plot.getBoundingClientRect();
    return r.width > 0 ? paeCellAt(m, (x - r.left) / r.width, (y - r.top) / r.height) : undefined;
  };
  const clampCell = (x: number, y: number) => {
    const r = plot.getBoundingClientRect();
    const f = (v: number) => Math.min(0.999999, Math.max(0, v));
    return paeCellAt(m, f((x - r.left) / r.width), f((y - r.top) / r.height))!;
  };
  /** What the pointer is on of the box: its edges (a corner is two), inside it, or neither. */
  const onBox = (x: number, y: number): { edges?: Edges; inside: boolean } | undefined => {
    const b = shown.value;
    const r = plot.getBoundingClientRect();
    if (!b || !(r.width > 0)) return undefined;
    const cell = r.width / n;
    const lx = x - r.left;
    const ly = y - r.top;
    const [x0, x1, y0, y1] = [b.j0 * cell, (b.j1 + 1) * cell, b.i0 * cell, (b.i1 + 1) * cell];
    const alongX = lx > x0 - Grab && lx < x1 + Grab;
    const alongY = ly > y0 - Grab && ly < y1 + Grab;
    const edges = {
      top: alongX && Math.abs(ly - y0) <= Grab,
      bottom: alongX && Math.abs(ly - y1) <= Grab,
      left: alongY && Math.abs(lx - x0) <= Grab,
      right: alongY && Math.abs(lx - x1) <= Grab,
    };
    if (edges.top || edges.bottom || edges.left || edges.right) return { edges, inside: false };
    return lx > x0 && lx < x1 && ly > y0 && ly < y1 ? { inside: true } : undefined;
  };
  const cursorFor = (hit: ReturnType<typeof onBox>) => {
    if (!hit) return 'crosshair';
    if (!hit.edges) return 'move';
    const { top, bottom, left, right } = hit.edges;
    if ((top && left) || (bottom && right)) return 'nwse-resize';
    if ((top && right) || (bottom && left)) return 'nesw-resize';
    return top || bottom ? 'ns-resize' : 'ew-resize';
  };
  /** The box a drag makes with the pointer on cell `at`. */
  const dragged = (d: Drag, at: { i: number; j: number }): PaeBox => {
    if (d.kind === 'move' && d.from) {
      const b = d.from;
      const di = Math.min(n - 1 - b.i1, Math.max(-b.i0, at.i - d.i));
      const dj = Math.min(n - 1 - b.j1, Math.max(-b.j0, at.j - d.j));
      return { i0: b.i0 + di, i1: b.i1 + di, j0: b.j0 + dj, j1: b.j1 + dj };
    }
    if (d.kind === 'resize' && d.from && d.edges) {
      const b = d.from;
      const e = d.edges;
      // The grabbed edges follow the pointer; past the opposite edge the box turns over.
      const ia = e.top ? at.i : b.i0;
      const ib = e.bottom ? at.i : b.i1;
      const ja = e.left ? at.j : b.j0;
      const jb = e.right ? at.j : b.j1;
      return { i0: Math.min(ia, ib), i1: Math.max(ia, ib), j0: Math.min(ja, jb), j1: Math.max(ja, jb) };
    }
    return { i0: Math.min(d.i, at.i), i1: Math.max(d.i, at.i), j0: Math.min(d.j, at.j), j1: Math.max(d.j, at.j) };
  };
  const update = () => {
    if (pressed) {
      const b = dragged(pressed, clampCell(px, py));
      drawBox(b);
      placeClose(b);
      drawHair(undefined);
      hairHidden = true;
      return;
    }
    plot.style.cursor = inside ? cursorFor(onBox(px, py)) : 'crosshair';
    const at = inside ? cellAt(px, py) : undefined;
    const key = at ? at.i * n + at.j : -1;
    if (key === cellKey && !hairHidden) return;
    cellKey = key;
    hairHidden = false;
    drawHair(at);
    hoverText.set(at ? paeCellLabel(m, at.i, at.j) : undefined);
    hover(at);
  };
  // Once a frame on the page's one tick, asked for from the pointer's listeners (`frames`: cancelled with the scope).
  const updating = frames(update);
  const schedule = () => {
    if (!disposed) updating.request();
  };
  // The box's drag (new, move, resize by its edges; `drag`: attached once, every way it can end ends it). A release
  // without a move picks the cell; a cancel puts the box back as it was.
  drag<Drag>(plot, {
    start: (e) => {
      if (e.button !== 0) return undefined;
      const at = cellAt(e.clientX, e.clientY);
      if (!at) return undefined;
      const hit = onBox(e.clientX, e.clientY);
      const kind = hit?.edges ? 'resize' : hit?.inside ? 'move' : 'new';
      pressed = { kind, ...at, x: e.clientX, y: e.clientY, from: shown.value, edges: hit?.edges, moved: false };
      return pressed;
    },
    move: (e, d) => {
      if (!d.moved && Math.hypot(e.clientX - d.x, e.clientY - d.y) > 4) d.moved = true;
      if (d.moved) schedule();
    },
    end: (d, how, e) => {
      pressed = undefined;
      if (how === 'cancel' || !e) {
        drawBox(shown.value);
        placeClose(shown.value);
        schedule();
        return;
      }
      if (!d.moved) {
        // A pair or a box, one at a time: while a box is shown, a click selects nothing.
        if (!shown.value) props.onPick?.(d.i, d.j);
        return;
      }
      const b = dragged(d, clampCell(e.clientX, e.clientY));
      shown.set(b);
      drawBox(b);
      placeClose(b);
      props.onBox?.(b);
      schedule();
    },
  });
  // The hover (the drag's moves are its own, above).
  on<PointerEvent>(plot, 'pointermove', (e) => {
    px = e.clientX;
    py = e.clientY;
    inside = true;
    if (!pressed) schedule();
  });
  on(plot, 'pointerleave', () => {
    inside = false;
    if (!pressed) schedule();
  });
  const clear = () => {
    shown.set(undefined);
    drawBox(undefined);
    placeClose(undefined);
    props.onBox?.(undefined);
  };
  onCleanup(() => {
    disposed = true;
    if (hoverTimer !== undefined) clearTimeout(hoverTimer);
    // Leaving: what the plot lit goes with it.
    if (cellKey >= 0) props.onHover?.(undefined);
    if (shown.value) props.onBox?.(undefined);
  });

  // The axes: residue ticks (scored across, under the plot; aligned down, left of it), and for a complex its chains as
  // bands along the top and the left.
  const ticks = paeTicks(m, 6);
  const pct = (k: number) => `${(100 * k) / n}%`;
  const complex = chains.length > 1;
  const band = (vertical: boolean) => (
    <div
      aria-hidden="true"
      style={`position: relative; ${vertical ? 'width: 14px; height: 100%' : 'height: 14px; width: 100%'}`}
    >
      {chains.map((c, k) => (
        <span
          style={`position: absolute; ${
            vertical
              ? `top: ${pct(c.from)}; height: ${pct(c.to - c.from)}; left: 0; right: 0`
              : `left: ${pct(c.from)}; width: ${pct(c.to - c.from)}; top: 0; bottom: 0`
          }; display: grid; place-items: center; overflow: hidden; font-size: 10px; font-weight: 600; color: var(--msp-color-text-muted); background: color-mix(in oklab, var(--msp-color-text), transparent ${k % 2 ? 90 : 84}%); border-radius: 2px`}
        >
          {c.to - c.from >= n / 24 ? c.label : ''}
        </span>
      ))}
    </div>
  );
  const yTicks = (
    <div aria-hidden="true" style="position: relative; width: 28px">
      {ticks.map((t) => (
        <span style={`${Tick}; right: 4px; top: ${pct(t.at + 0.5)}; transform: translateY(-50%)`}>{t.label}</span>
      ))}
    </div>
  );
  const xTicks = (
    <div aria-hidden="true" style="position: relative; height: 14px">
      {ticks.map((t) => (
        <span style={`${Tick}; top: 3px; left: ${pct(t.at + 0.5)}; transform: translateX(-50%)`}>{t.label}</span>
      ))}
    </div>
  );

  // The scale: the ramp with a tick each 5 Å (each 10 Å past 40).
  // pLDDT along both axes (AlphaFold DB's confidence bands), a pixel a residue as the plot's; shown until turned off.
  const plddt = m.plddt;
  const showPlddt = observableState(true);
  const stripCanvas = (vertical: boolean) => {
    const c = (
      <canvas
        width={vertical ? 1 : n}
        height={vertical ? n : 1}
        aria-hidden="true"
        style="position: absolute; inset: 0; width: 100%; height: 100%; image-rendering: pixelated; border-radius: 1px"
      />
    ) as HTMLCanvasElement;
    const g2 = plddt && c.getContext('2d');
    if (g2) {
      const img = g2.createImageData(c.width, c.height);
      for (let k = 0; k < n; k++) {
        const color = plddtColor(plddt[k]);
        img.data[k * 4] = (color >> 16) & 0xff;
        img.data[k * 4 + 1] = (color >> 8) & 0xff;
        img.data[k * 4 + 2] = color & 0xff;
        img.data[k * 4 + 3] = Number.isNaN(plddt[k]) ? 0 : 255;
      }
      g2.putImageData(img, 0, 0);
    }
    return c;
  };
  /** A grid item shown with the strips only. */
  const withStrip = (child: unknown) => (
    <div
      style={derive(showPlddt.asReadonly, (on) =>
        on ? 'display: block; position: relative; min-width: 0; min-height: 0' : 'display: none',
      )}
    >
      {child as never}
    </div>
  );
  const plddtLegend = (
    <div
      style={derive(
        showPlddt.asReadonly,
        (on) => `display: ${on ? 'flex' : 'none'}; flex-wrap: wrap; gap: 4px 10px; ${Quiet}`,
      )}
      aria-label="pLDDT bands"
    >
      {PlddtBands.map((b) => (
        <span>
          <span
            style={`display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 4px; vertical-align: -1px; background: ${hex(b.color)}`}
          />
          {b.label} {b.range}
        </span>
      ))}
    </div>
  );

  const step = m.max > 40 ? 10 : 5;
  const scaleTicks = Array.from({ length: Math.floor(m.max / step) + 1 }, (_, k) => k * step);
  const scale = (
    <div style="display: grid; gap: 4px" aria-label={`Expected position error, 0 to ${m.max.toFixed(0)} Å`}>
      <div
        style={`height: 8px; border-radius: 4px; background: linear-gradient(to right, ${PaeColorStyles.low}, ${PaeColorStyles.high}); box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 85%)`}
      />
      <div aria-hidden="true" style="position: relative; height: 12px">
        {scaleTicks.map((v) => (
          <span
            style={`${Tick}; top: 0; left: ${(100 * v) / (m.max || 1)}%; transform: translateX(${v === 0 ? '0' : v >= m.max ? '-100%' : '-50%'})`}
          >
            {String(v)}
          </span>
        ))}
      </div>
      <div style="display: flex; align-items: center; justify-content: space-between; gap: var(--msp-space-2)">
        <span style={Quiet}>Expected position error (Å)</span>
        {plddt ? (
          // Inside the card's edge, as everything else: pressed, its background ends where the plot does.
          <Button
            variant="ghost"
            size="xs"
            pressed={showPlddt.asReadonly}
            title="Show each residue's pLDDT along the axes"
            onClick={() => showPlddt.set(!showPlddt.value)}
          >
            pLDDT
          </Button>
        ) : null}
      </div>
      {plddt ? plddtLegend : null}
    </div>
  );

  const idle = 'Point at a cell; drag a box to compare two ranges.';
  return (
    <div style="display: grid; gap: var(--msp-space-3); min-width: 0">
      <div
        style={derive(showPlddt.asReadonly, (on) => {
          const strip = on && plddt ? '6px ' : '';
          return `display: grid; grid-template-columns: 28px ${complex ? '14px ' : ''}${strip}minmax(0, 1fr); grid-template-rows: ${complex ? '14px ' : ''}${strip}auto auto; column-gap: 2px; row-gap: 2px`;
        })}
      >
        {complex ? (
          <>
            <span />
            <span />
            {plddt ? withStrip(null) : null}
            {band(false)}
          </>
        ) : null}
        {plddt ? (
          <>
            {withStrip(null)}
            {complex ? withStrip(null) : null}
            {withStrip(null)}
            {withStrip(stripCanvas(false))}
          </>
        ) : null}
        {yTicks}
        {complex ? band(true) : null}
        {plddt ? withStrip(stripCanvas(true)) : null}
        {plotArea}
        <span />
        {complex ? <span /> : null}
        {plddt ? withStrip(null) : null}
        {xTicks}
      </div>
      <div style={`display: flex; justify-content: space-between; ${Quiet}`} aria-hidden="true">
        <span>↓ Aligned residue</span>
        <span>Scored residue →</span>
      </div>
      {/* What the pointer is on: one line, always there (a box's summary goes at the end, so nothing above moves). */}
      <p style={`margin: 0; min-height: 1.5em; ${Quiet}`} aria-live="polite">
        {derive(hoverText.asReadonly, (t) => t ?? idle)}
      </p>
      {scale}
      {/* The box's two ranges in their colors and its mean error, last. */}
      <Show when={shown.asReadonly}>
        {(b) => <BoxSummary matrix={m} box={b} colors={tones.asReadonly} onClear={clear} />}
      </Show>
    </div>
  );
}

/** The box's ranges, each with its color's swatch, and its mean error; ✕ clears it (and its marks). */
function BoxSummary(props: { matrix: PaeMatrix; box: PaeBox; colors: ReadonlyState<PaeMarkColors>; onClear(): void }) {
  const { matrix: m, box: b } = props;
  const mean = paeBoxMean(m, b.i0, b.i1, b.j0, b.j1);
  const swatch = (tone: keyof PaeMarkColors) =>
    derive(
      props.colors,
      (c) =>
        `display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; background: ${hex(c[tone])}`,
    );
  return (
    <div
      aria-live="polite"
      style="display: grid; grid-template-columns: minmax(0, 1fr) auto; align-items: center; gap: var(--msp-space-2); padding: var(--msp-space-2) var(--msp-space-3); border-radius: var(--msp-radius-md); background: var(--msp-color-well); font-size: var(--msp-font-size-xs)"
    >
      <div style="display: grid; gap: 2px; min-width: 0">
        <span>
          <span style={swatch('highlight')} />
          Aligned {paeRangeLabel(m, b.i0, b.i1)}
        </span>
        <span>
          <span style={swatch('accent')} />
          Scored {paeRangeLabel(m, b.j0, b.j1)}
        </span>
        <span style="color: var(--msp-color-text-muted)">
          Mean error {Number.isNaN(mean) ? '–' : `${mean.toFixed(1)} Å`}
        </span>
      </div>
      <Button
        variant="ghost"
        size="xs"
        icon={<Icon name="x" size={14} />}
        title="Clear the box"
        onClick={props.onClear}
      />
    </div>
  );
}
