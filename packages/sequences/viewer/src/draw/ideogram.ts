/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ideogram navigator (ST5; seqspec:viewer#extension-points): each block's axis drawn as a chromosome's karyotype
// strip, its cytobands as Giemsa stains (grays from light to dark through the skin's ink, so it follows the theme), the
// centromere's two `acen` bands as triangles meeting in the middle, a stalk or a variable band lighter; the view's
// window over it, washed outside and outlined in two tones, as over any map. The bands are named data at app level
// (`CytoBands`, UCSC's cytoBandIdeo rows), found by the block's name; a block with none gets a hairline. The strip is
// drawn on the nav pass (a layout, size or skin change, never a pan), so a frame draws nothing here but the window's
// own, and every path is allocation-free once warm: module scratch, grown only. Its gestures are the `navigator`
// behavior's, as any navigator's.

import type { Axis } from "../axis.ts";
import { type NavigatorEntry, navigatorFrom } from "../navigators.ts";
import { PathBatch, PathVerb, RectBatch, SkinToken } from "./backend.ts";
import type { NavView } from "./nav.ts";

/** A genome's bands by chromosome name: [start, stop) 0-based, the band's name and its Giemsa stain (cytoBandIdeo). */
export type CytoBands = Readonly<Record<string, readonly (readonly [number, number, string, string])[]>>;

// A band's look by stain code: 0 negative (the ground shows), 1 to 4 positive by quarters, 5 and 6 the centromere's
// two halves (triangles), 7 variable, 8 stalk. Ink at an alpha; the centromere a fixed red.
const STAINS = ["gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar", "stalk"];
const [ACEN_LEFT, ACEN_RIGHT] = [5, 6];
const ALPHA = new Uint8Array([0, 56, 112, 176, 236, 0, 0, 80, 40]);
const CENTROMERE = 0xb5483f;
// Tokens and verbs read once: loaded from their exported objects per frame they allocated about 41 bytes in the heap
// guard (ideogram.heap.ts), under the test runner's module wrapper.
const [SURFACE, INK, QUIET, FAINT] = [SkinToken.surface, SkinToken.ink, SkinToken.inkQuiet, SkinToken.inkFaint];
const [MOVE, LINE, CLOSE] = [PathVerb.move, PathVerb.line, PathVerb.close];

interface Chromosome {
  readonly start: Int32Array;
  readonly stop: Int32Array;
  readonly stain: Uint8Array;
}
/** A chromosome's name as the blocks name it: "chr17" and "17" are one. */
const keyOf = (name: string) => name.replace(/^chr/i, "").toLowerCase();
/** The chromosome a block's name stands for: its last segment that names one ("hg38/chr17", "chr17/unit"). */
function chromosomeOf(all: ReadonlyMap<string, Chromosome>, name: string): Chromosome | undefined {
  const parts = name.split("/");
  for (let k = parts.length - 1; k >= 0; k--) {
    const c = all.get(keyOf(parts[k]));
    if (c) return c;
  }
  return undefined;
}

function chromosomes(bands: CytoBands): Map<string, Chromosome> {
  const out = new Map<string, Chromosome>();
  for (const [name, rows] of Object.entries(bands)) {
    const stain = new Uint8Array(rows.length);
    let acen = 0;
    rows.forEach(([, , , s], k) => {
      // The first `acen` band is the short arm's, the second the long arm's.
      const code = STAINS.indexOf(s);
      stain[k] = code < 0 ? 0 : code === 5 ? (acen++ === 0 ? ACEN_LEFT : ACEN_RIGHT) : code < 5 ? code : code + 1;
    });
    out.set(keyOf(name), {
      start: Int32Array.from(rows, (r) => r[0]),
      stop: Int32Array.from(rows, (r) => r[1]),
      stain,
    });
  }
  return out;
}

const ground = new RectBatch(32);
const strips = new RectBatch(512);
const outline = new RectBatch(32);
const hair = new RectBatch(32);
const pinch = new PathBatch(32);
let shown = new Int32Array(8);
/** The device pixel ratio and a block's CSS px a shown index, across the runs' calls (a double kept over a call is boxed). */
const scale = new Float64Array(2);

/** The blocks' chromosomes, found again when the axis's names are another array. */
class Blocks {
  names: readonly string[] | null = null;
  chromosomes: (Chromosome | undefined)[] = [];
}

function drawIdeogram(n: NavView, all: ReadonlyMap<string, Chromosome>, blocks: Blocks): void {
  const ax = n.axis as Axis | null;
  n.overMap = false;
  if (!ax || n.backend === null) return;
  if (blocks.names !== n.names) {
    blocks.names = n.names;
    blocks.chromosomes = n.names.map((name) => chromosomeOf(all, name));
  }
  n.overMap = true;
  ground.reset();
  strips.reset();
  outline.reset();
  hair.reset();
  pinch.reset();
  scale[0] = n.dpr;
  const th = Math.min(n.parts.thickness, n.height);
  const y0 = Math.round((n.top + (n.height - th) / 2) * scale[0]) | 0;
  const h = Math.max(1, (Math.round((n.top + (n.height + th) / 2) * scale[0]) | 0) - y0);
  const lw = Math.max(1, Math.round(scale[0])) | 0;
  ground.begin(SURFACE);
  outline.begin(QUIET, 255, lw);
  hair.begin(FAINT);
  pinch.begin(CENTROMERE);
  for (let b = 0; b < ax.count; b++) {
    const x0 = Math.round(n.x0s[b] * scale[0]) | 0;
    const w = Math.max(1, (Math.round(n.x1s[b] * scale[0]) | 0) - x0);
    const c = b < blocks.chromosomes.length ? blocks.chromosomes[b] : undefined;
    if (c === undefined || ax.shown[b] === 0) {
      hair.push(x0, y0 + ((h - lw) >> 1), w, lw);
      continue;
    }
    ground.push(x0, y0, w, h);
    scale[1] = (n.x1s[b] - n.x0s[b]) / ax.shown[b];
    for (let q = 0; q < c.stain.length; q++) {
      const s = c.stain[q];
      if (s === 0) continue;
      let runs = ax.runToShown(b, c.start[q], c.stop[q], shown);
      if (2 * runs > shown.length) {
        shown = new Int32Array(4 * runs);
        runs = ax.runToShown(b, c.start[q], c.stop[q], shown);
      }
      for (let r = 0; r < runs; r++) {
        const xa = Math.round((n.x0s[b] + shown[2 * r] * scale[1]) * scale[0]) | 0;
        const xb = Math.max(xa + 1, Math.round((n.x0s[b] + shown[2 * r + 1] * scale[1]) * scale[0]) | 0);
        if (s === ACEN_LEFT || s === ACEN_RIGHT) {
          // Two triangles meeting at the centromere: the short arm's points toward the long arm, which points back; on
          // a reversed piece the arms swap sides, and so do the apexes (GM-06).
          const left = (s === ACEN_LEFT) === ax.dir[b] > 0;
          const o = pinch.verb(MOVE);
          pinch.coords[o] = left ? xa : xb;
          pinch.coords[o + 1] = y0;
          const p = pinch.verb(LINE);
          pinch.coords[p] = left ? xb : xa;
          pinch.coords[p + 1] = y0 + (h >> 1);
          const e = pinch.verb(LINE);
          pinch.coords[e] = left ? xa : xb;
          pinch.coords[e + 1] = y0 + h;
          pinch.verb(CLOSE);
        } else {
          strips.begin(INK, ALPHA[s]);
          strips.push(xa, y0, xb - xa, h);
        }
      }
    }
    outline.push(x0, y0, w, h);
  }
  if (ground.n > 0) n.backend.rects(ground);
  if (strips.n > 0) n.backend.rects(strips);
  if (pinch.nv > 0) n.backend.paths(pinch);
  if (outline.n > 0) n.backend.rects(outline);
  if (hair.n > 0) n.backend.rects(hair);
}

/**
 * A navigator drawing each block's chromosome as its ideogram from `bands` (24 px, the strip 14): register it
 * (`{ navigators: [ideogramNavigator(bands)] }`) and name it in a preset's `navigator`. A block is found by the last
 * segment of its name that names a chromosome, with or without "chr" (`hg38/chr17`, `17/unit`).
 */
export function ideogramNavigator(bands: CytoBands, name = "ideogram"): NavigatorEntry {
  const all = chromosomes(bands);
  const blocks = new Blocks();
  const base = navigatorFrom({ height: 24, thickness: 14, axis: [], ground: false, ticks: false }, name);
  return { ...base, draw: (n) => drawIdeogram(n, all, blocks) };
}
