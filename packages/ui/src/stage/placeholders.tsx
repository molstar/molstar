/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { frames, on, uniqueId } from '../dom/index.js';
import { classes } from './_generated/placeholders.module.css.js';
import type { HoverLine } from './stage-model.js';

// Stand-ins for stories and tests (stage §6): a canvas and a band of the right shape, with no WebGL and no sequence
// viewer. They feed the stage's hover lines as the real ones would.

/** The stand-in structure's shape: chains, residues a chain, atoms a residue. */
const Chains = 3;
const ChainResidues = 70;
const ResidueAtoms = 4;
/** A chain's color (Mol*'s chain-id look: a blue, a gold, a rose). */
const ChainColors = ['#5b8fd6', '#e3b04b', '#d4667f'];
/** Depth bands: the farther, the darker (a cheap depth cue, one gradient a chain and band). */
const Depths = 5;

interface Atom {
  readonly x: number;
  readonly y: number;
  readonly z: number;
  readonly r: number;
  readonly chain: number;
  /** The residue, over every chain (what hover tells). */
  readonly residue: number;
}

/** A small seeded generator (mulberry32): the same structure every time. */
function generator(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

let structure: { atoms: readonly Atom[]; box: readonly [number, number, number, number] } | undefined;

/**
 * The stand-in structure, made once: each chain a self-avoiding-ish walk of residues (3.8 Å apart, turning a little
 * each step) kept inside a globule around its own center, the centers touching, so the whole reads as a three-chain
 * assembly; a residue's atoms around its trace point. Seen from a fixed angle, sorted back to front.
 */
function makeStructure() {
  if (structure) return structure;
  const random = generator(1);
  const atoms: Atom[] = [];
  const centers = [
    [-9, 3, 0],
    [9, 5, -2],
    [1, -11, 3],
  ];
  const radius = 9;
  for (let c = 0; c < Chains; c++) {
    const [cx, cy, cz] = centers[c];
    let x = cx + (random() - 0.5) * 6;
    let y = cy + (random() - 0.5) * 6;
    let z = cz + (random() - 0.5) * 6;
    let dx = random() - 0.5;
    let dy = random() - 0.5;
    let dz = random() - 0.5;
    for (let i = 0; i < ChainResidues; i++) {
      // Turn a little, and back toward the center once outside the globule.
      dx += (random() - 0.5) * 1.2;
      dy += (random() - 0.5) * 1.2;
      dz += (random() - 0.5) * 1.2;
      const ox = x - cx;
      const oy = y - cy;
      const oz = z - cz;
      const out = Math.sqrt(ox * ox + oy * oy + oz * oz);
      if (out > radius) {
        dx -= ox / out;
        dy -= oy / out;
        dz -= oz / out;
      }
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1;
      dx /= d;
      dy /= d;
      dz /= d;
      x += dx * 3.8;
      y += dy * 3.8;
      z += dz * 3.8;
      const residue = c * ChainResidues + i;
      for (let a = 0; a < ResidueAtoms; a++) {
        atoms.push({
          x: x + (random() - 0.5) * 3.2,
          y: y + (random() - 0.5) * 3.2,
          z: z + (random() - 0.5) * 3.2,
          r: 2 + random() * 0.4,
          chain: c,
          residue,
        });
      }
    }
  }
  // A fixed view: turned about y, then tilted about x.
  const ay = 0.6;
  const ax = -0.35;
  const seen = atoms.map((p): Atom => {
    const x1 = p.x * Math.cos(ay) + p.z * Math.sin(ay);
    const z1 = -p.x * Math.sin(ay) + p.z * Math.cos(ay);
    const y2 = p.y * Math.cos(ax) - z1 * Math.sin(ax);
    const z2 = p.y * Math.sin(ax) + z1 * Math.cos(ax);
    return { ...p, x: x1, y: y2, z: z2 };
  });
  seen.sort((a, b) => a.z - b.z);
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  for (const a of seen) {
    minX = Math.min(minX, a.x - a.r);
    maxX = Math.max(maxX, a.x + a.r);
    minY = Math.min(minY, a.y - a.r);
    maxY = Math.max(maxY, a.y + a.r);
  }
  structure = { atoms: seen, box: [minX, minY, maxX - minX, maxY - minY] };
  return structure;
}

/** A color darkened by a factor (0: as is, 1: black), as hex. */
function darken(hex: string, by: number): string {
  const n = Number.parseInt(hex.slice(1), 16);
  const f = 1 - by;
  const c = (v: number) =>
    Math.round(v * f)
      .toString(16)
      .padStart(2, '0');
  return `#${c((n >> 16) & 255)}${c((n >> 8) & 255)}${c(n & 255)}`;
}

/** A residue's hover lines (its index over every chain), as a structure's label would read. */
export function structureLines(index: number): readonly HoverLine[] {
  const chain = Math.floor(index / ChainResidues);
  const c = Residues[(index * 7) % Residues.length];
  return [
    { runs: [{ text: `${c} ${(index % ChainResidues) + 1}`, bold: true }] },
    {
      runs: [{ text: 'Chain ' }, { text: String.fromCharCode(65 + chain), code: true }, { text: ` · ${kind(c)}` }],
    },
  ];
}

export interface StagePlaceholderProps {
  /** Nothing drawn (an empty canvas). */
  empty?: boolean;
  /** A backdrop of its own (a CSS background); default none: the stage's color shows through. */
  background?: string;
  /** The residue under the pointer (its index over every chain), or none. */
  onHover?: (index: number | undefined) => void;
  /** The chains' colors (three; default a blue, a gold, a rose): another pane, another look. */
  colors?: readonly string[];
}

/**
 * A stand-in for a 3D canvas: a three-chain assembly drawn space-filling (shaded atoms, the farther darker, colored by
 * chain) on the stage's color, scaled to fit, so layouts and the chrome over them can be judged against what a loaded
 * structure looks like. `onHover` gets the residue under the pointer.
 */
export function StagePlaceholder(props: StagePlaceholderProps) {
  const style = props.background ? { background: props.background } : undefined;
  if (props.empty) return <div class={classes.canvas} style={style} />;
  const { atoms, box } = makeStructure();
  const id = uniqueId('msui-atom');
  const gradients: SVGElement[] = [];
  for (let c = 0; c < Chains; c++) {
    for (let d = 0; d < Depths; d++) {
      const base = darken((props.colors ?? ChainColors)[c % Chains], 0.45 * (1 - d / (Depths - 1)));
      gradients.push(
        (
          <radialGradient id={`${id}-${c}-${d}`} cx="38%" cy="35%" r="70%">
            <stop offset="0%" stop-color="#fff" stop-opacity="0.6" />
            <stop offset="35%" stop-color={base} />
            <stop offset="100%" stop-color={darken(base, 0.45)} />
          </radialGradient>
        ) as SVGElement,
      );
    }
  }
  const minZ = atoms[0].z;
  const spanZ = atoms[atoms.length - 1].z - minZ || 1;
  const circles = atoms.map((a) => {
    const d = Math.min(Depths - 1, Math.floor(((a.z - minZ) / spanZ) * Depths));
    return (
      <circle
        cx={a.x.toFixed(2)}
        cy={a.y.toFixed(2)}
        r={a.r.toFixed(2)}
        fill={`url(#${id}-${a.chain}-${d})`}
        data-index={a.residue}
      />
    ) as SVGElement;
  });
  // Room around the structure: a tenth of its size each way.
  const [bx, by, bw, bh] = box;
  const viewBox = `${(bx - bw * 0.1).toFixed(1)} ${(by - bh * 0.1).toFixed(1)} ${(bw * 1.2).toFixed(1)} ${(bh * 1.2).toFixed(1)}`;
  let hovered: number | undefined;
  const hover = (index: number | undefined) => {
    if (index === hovered) return;
    hovered = index;
    props.onHover?.(index);
  };
  const root = (
    <div class={classes.canvas} style={style}>
      <svg class={classes.structure} viewBox={viewBox} preserveAspectRatio="xMidYMid meet" aria-hidden="true">
        <defs>{gradients}</defs>
        <g>{circles}</g>
      </svg>
    </div>
  ) as HTMLElement;
  if (props.onHover) {
    // Over and out bubble: two listeners for the structure, not one per atom.
    on<PointerEvent>(root, 'pointerover', (e) => {
      const index = (e.target as Element).getAttribute?.('data-index');
      hover(index === null || index === undefined ? undefined : Number(index));
    });
    on(root, 'pointerleave', () => hover(undefined));
  }
  return root;
}

/** One-letter residues the band's rows repeat (a protein-like mix). */
const Residues =
  'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAACRESEHAKAMLCLLLAGTDVNLRDLRTGLTPLHLAARVGHRDAAEALLSAGA';
/** A cell (px): a residue's box. */
const CellW = 12;
const CellH = 16;
/** The ruler's height (px). */
const RulerH = 16;

/** The residue kinds, by the letter (a box's color). */
function kind(c: string): string {
  if ('AVLIMFWC'.includes(c)) return 'hydrophobic';
  if ('KRH'.includes(c)) return 'positive';
  if ('DE'.includes(c)) return 'negative';
  if ('STNQ'.includes(c)) return 'polar';
  return 'other';
}

export interface BandPlaceholderProps {
  /** Rows (chains), default 2. */
  rows?: number;
  /** Residues a row, default 120. */
  length?: number;
  /** A selected range on the first row: [first, last] residue indices. */
  selection?: readonly [number, number];
  /** The residue under the pointer (its row and index), or none. */
  onHover?: (hit: { row: number; index: number } | undefined) => void;
}

/** A residue's hover lines, as a sequence's label would read. */
export function residueLines(row: number, index: number): readonly HoverLine[] {
  const c = Residues[(index + row * 37) % Residues.length];
  return [
    { runs: [{ text: `${c} ${index + 1}`, bold: true }] },
    { runs: [{ text: `Chain ${String.fromCharCode(65 + row)} · ${kind(c)}` }] },
  ];
}

/**
 * A stand-in for a sequence band (stage §6): rows of residue-like boxes under a ruler, a selection on the first row and
 * the hovered column marked, scrolling across when wider than the stage. The column under the pointer is written once
 * a frame (one transform), and `onHover` gets the residue. Marked `data-msui-hover-source`: the stage's hover label
 * shows over it.
 */
export function BandPlaceholder(props: BandPlaceholderProps) {
  const rows = Math.max(1, props.rows ?? 2);
  const length = Math.max(1, props.length ?? 120);
  const width = length * CellW;
  const height = RulerH + rows * CellH;
  const xs = Array.from({ length }, (_, i) => String(i * CellW + CellW / 2)).join(' ');
  const ticks: SVGElement[] = [];
  for (let i = 9; i < length; i += 10) {
    const x = i * CellW + CellW / 2;
    ticks.push((<line x1={x} x2={x} y1={RulerH - 4} y2={RulerH} />) as SVGElement);
    ticks.push(
      (
        <text x={x} y={RulerH - 6}>
          {String(i + 1)}
        </text>
      ) as SVGElement,
    );
  }
  const boxes: SVGElement[] = [];
  const letters: SVGElement[] = [];
  for (let r = 0; r < rows; r++) {
    let text = '';
    for (let i = 0; i < length; i++) {
      const c = Residues[(i + r * 37) % Residues.length];
      text += c;
      boxes.push(
        (
          <rect
            x={i * CellW + 1}
            y={RulerH + r * CellH + 2}
            width={CellW - 2}
            height={CellH - 4}
            rx="2"
            data-kind={kind(c)}
          />
        ) as SVGElement,
      );
    }
    letters.push(
      (
        <text x={xs} y={RulerH + r * CellH + CellH / 2}>
          {text}
        </text>
      ) as SVGElement,
    );
  }
  const selection = props.selection;
  const column = (<rect class={classes.column} x="0" y="0" width={CellW} height={height} />) as SVGElement;
  const svg = (
    <svg class={classes.residues} width={width} height={height} viewBox={`0 0 ${width} ${height}`} aria-hidden="true">
      <g class={classes.ruler}>{ticks}</g>
      {selection ? (
        <rect
          class={classes.selection}
          x={selection[0] * CellW}
          y={RulerH}
          width={(selection[1] - selection[0] + 1) * CellW}
          height={CellH}
          rx="3"
        />
      ) : null}
      <g class={classes.boxes}>{boxes}</g>
      <g class={classes.letters}>{letters}</g>
      {column}
    </svg>
  ) as SVGSVGElement;
  const root = (
    <div
      class={classes.band}
      data-msui-hover-source=""
      role="img"
      aria-label={`A sequence band: ${rows} rows of ${length} residues`}
    >
      {svg}
    </div>
  ) as HTMLElement;

  // The pointer's column, applied once a frame: one attribute write, and the hover told when the residue changes.
  let x = 0;
  let y = 0;
  let over = false;
  let shownIndex = -1;
  let shownRow = -1;
  const apply = frames(() => {
    const index = over ? Math.floor(x / CellW) : -1;
    const row = over ? Math.floor((y - RulerH) / CellH) : -1;
    const hit = index >= 0 && index < length && row >= 0 && row < rows;
    const i = hit ? index : -1;
    const r = hit ? row : -1;
    if (i === shownIndex && r === shownRow) return;
    if (i !== shownIndex) column.setAttribute('transform', `translate(${Math.max(0, i) * CellW} 0)`);
    shownIndex = i;
    shownRow = r;
    column.toggleAttribute('data-shown', hit);
    props.onHover?.(hit ? { row: r, index: i } : undefined);
  });
  // The svg's offset: read when the pointer comes in and after a scroll of the band, not per move.
  let left = 0;
  let top = 0;
  let stale = true;
  on<PointerEvent>(svg, 'pointermove', (e) => {
    if (stale) {
      stale = false;
      const r = svg.getBoundingClientRect();
      left = r.left;
      top = r.top;
    }
    x = e.clientX - left;
    y = e.clientY - top;
    over = true;
    apply.request();
  });
  on(svg, 'pointerenter', () => {
    stale = true;
  });
  on(root, 'scroll', () => {
    stale = true;
  });
  on(svg, 'pointerleave', () => {
    over = false;
    apply.request();
  });
  return root;
}
