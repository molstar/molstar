/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  DrawLists,
  drawListsToSvg,
  type Frame,
  type GlyphList,
  type ImageList,
  type PassList,
  type PathList,
  PathCommandPoints,
  type RectList,
  type SectionLayout,
  walkPass,
  WarpKind,
} from '@molstar/seqx-core/render';
import { CanvasBackend, type CanvasLike, type Context2D } from '@molstar/seqx-viewer/canvas';
import { describerFor } from '@molstar/seqx-viewer/describe';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The hover contract's sweep (_spec/README.md "Render", Emphasis, The hover contract): points on the marks each layer
// instance drew, read from its own data pass (drawn again alone into scratch lists, so a mark is a layer's by
// construction, never guessed from the frame's merged lists), carried to the page through its strip's warp; the
// pointer moved to each, and what the viewer answers checked.

/** A point on a drawn mark of a layer: its page point, the section it is in and the layer's key and type. */
export interface Probe {
  readonly x: number;
  readonly y: number;
  /** Its section's address, or a bridge's id. */
  readonly section: string;
  readonly layer: string;
  readonly type: string;
  /** What the mark is: `rect`, `path`, `letter`, `image`. */
  readonly mark: string;
}

type Point = { x: number; y: number; mark: string };
type Bounds = { x0: number; y0: number; x1: number; y1: number };

const scratch = new DrawLists();

/** Points on a pass's marks, in its frame: a rect's centre, a path's middle point on its curve, a letter's, an image's. */
function marksOf(p: PassList, out: Point[]) {
  walkPass(p, {
    rect(l: RectList, i: number) {
      if (l.w[i] < 0.5 || l.h[i] < 0.5) return;
      out.push({ x: l.x[i] + l.w[i] / 2, y: l.y[i] + l.h[i] / 2, mark: 'rect' });
    },
    path(l: PathList, i: number) {
      // the end point of its middle command: on the curve (a cubic's control points are not)
      const c0 = l.cmdStart[i];
      const n = l.cmdCount[i];
      if (!n) return;
      let pt = l.pointStart[i];
      const mid = c0 + Math.floor(n / 2);
      for (let c = c0; c < mid; c++) pt += PathCommandPoints[l.commands[c]];
      const k = l.commands[mid] === 3 ? pt - 1 : pt + PathCommandPoints[l.commands[mid]] - 1;
      if (k < 0) return;
      // a pixel in from the extent drawn so far (an area's baseline point lies on its band's foot, outside the band)
      let x0 = Infinity;
      let y0 = Infinity;
      let x1 = -Infinity;
      let y1 = -Infinity;
      for (let q = l.pointStart[i]; q <= k; q++) {
        x0 = Math.min(x0, l.points[2 * q]);
        x1 = Math.max(x1, l.points[2 * q]);
        y0 = Math.min(y0, l.points[2 * q + 1]);
        y1 = Math.max(y1, l.points[2 * q + 1]);
      }
      const x = l.points[2 * k];
      const y = l.points[2 * k + 1];
      out.push({
        x: Math.min(Math.max(x, x0 + 1), x1 - 1) || x,
        y: Math.min(Math.max(y, y0 + 1), y1 - 1) || y,
        mark: 'path',
      });
    },
    glyphs(l: GlyphList, i: number) {
      const j = l.start[i] + Math.floor(l.count[i] / 2);
      if (!l.count[i]) return;
      // a letter's centre, a little above its baseline
      const y = l.placed[i] ? l.ly[j] : l.y[i];
      out.push({ x: l.x[j], y: y - 3, mark: 'letter' });
    },
    text() {},
    image(l: ImageList, i: number) {
      // the pixels it paints (a density strip's empty bins are clear): a few spread over them, at their centres
      const im = l.images[i];
      const lit: number[] = [];
      for (let k = 0; k < im.width * im.height; k++) if (im.data[4 * k + 3] > 0) lit.push(k);
      for (const k of spread(lit, 3)) {
        const px = k % im.width;
        const py = Math.floor(k / im.width);
        out.push({
          x: l.x[i] + ((px + 0.5) * l.w[i]) / im.width,
          y: l.y[i] + ((py + 0.5) * l.h[i]) / im.height,
          mark: 'image',
        });
      }
    },
    clip() {},
    unclip() {},
  });
}

/** Up to `per` of a list, spread over it (its first, its last and evenly between). */
function spread<T>(xs: readonly T[], per: number): T[] {
  if (xs.length <= per) return [...xs];
  const out: T[] = [];
  for (let k = 0; k < per; k++) out.push(xs[Math.round((k * (xs.length - 1)) / (per - 1))]);
  return out;
}

/** Whether a page point lies in a section's visible box. */
function inBox(s: SectionLayout, x: number, y: number): boolean {
  const b = s.placement.box;
  return x >= b.x && x < b.x + b.width && y >= b.y && y < b.y + b.height;
}

/**
 * The probes of a frame: for every layer instance (a strip's band's layers, a section's own layers, a bridge's
 * ribbons), up to `per` points on the marks its data pass drew, inside its section's visible box (in its band).
 */
export function probesOf(frame: Frame, per = 4): Probe[] {
  const out: Probe[] = [];
  const pass = scratch.passes.data;
  const w = scratch.writer('data');
  const pts: Point[] = [];
  const take = (
    section: string,
    type: string,
    layer: string,
    map: (x: number, y: number) => [number, number],
    inside: (x: number, y: number) => boolean,
  ) => {
    pts.length = 0;
    marksOf(pass, pts);
    const placed: Probe[] = [];
    for (const p of pts) {
      const [x, y] = map(p.x, p.y);
      if (!Number.isFinite(x) || !Number.isFinite(y) || !inside(x, y)) continue;
      placed.push({ x, y, section, layer, type, mark: p.mark });
    }
    out.push(...spread(placed, per));
  };
  for (const s of frame.layout.sections) {
    const st = s.placement.strips;
    const nb = st.bands;
    const inSection = (x: number, y: number) => inBox(s, x, y);
    for (let i = 0; i < st.n; i++) {
      const polar = st.warp[i] === WarpKind.polar ? s.placement.warpOf?.(i) : undefined;
      const dx = st.dx[i];
      const dy = st.dy[i];
      const map = polar
        ? (x: number, y: number): [number, number] => {
            polar.map(x, y);
            return [polar.px, polar.py];
          }
        : (x: number, y: number): [number, number] => [x + dx, y + dy];
      for (let j = 0; j < nb; j++) {
        // a mark is seen inside its band alone (the data is clipped to its well, a block's to its strip's box; hit
        // tests read the band)
        const by = st.bandY[i * nb + j];
        const bh = st.bandH[i * nb + j];
        const box = s.boxes[i * nb + j];
        const inBand = (x: number, y: number): [number, number] =>
          y >= by && y < by + bh && x >= box.x && x < box.x + box.width ? map(x, y) : [Number.NaN, Number.NaN];
        for (const l of s.lays[i * nb + j] ?? []) {
          if (l.zoomIn) continue;
          scratch.reset(['data']);
          l.type.draw(l.lay, 'data', w);
          take(s.section.address, l.type.name, l.key, l.layer.params.across ? map : inBand, inSection);
        }
      }
    }
    for (const l of s.sectionLays) {
      scratch.reset(['data']);
      l.type.draw(l.lay, 'data', w);
      take(s.section.address, l.type.name, l.key, (x, y) => [x, y], inSection);
    }
  }
  const bridges = frame.layout.bridges;
  if (bridges) {
    scratch.reset(['data']);
    bridges.draw('data', w, undefined);
    for (const b of bridges.list)
      take(
        b.id,
        'bridge',
        b.id,
        (x, y) => [x, y],
        (x, y) => bridges.at(x, y) === b.index,
      );
  }
  return out;
}

/** The boxes of a pass's primitives on the page (rects, paths' extents, images): where its emphasis was drawn. */
export function passBoxes(p: PassList): Bounds[] {
  const out: Bounds[] = [];
  walkPass(p, {
    rect(l: RectList, i: number) {
      out.push({ x0: l.x[i], y0: l.y[i], x1: l.x[i] + l.w[i], y1: l.y[i] + l.h[i] });
    },
    path(l: PathList, i: number) {
      let pt = l.pointStart[i];
      let n = 0;
      for (let c = l.cmdStart[i]; c < l.cmdStart[i] + l.cmdCount[i]; c++) n += PathCommandPoints[l.commands[c]];
      let x0 = Infinity;
      let y0 = Infinity;
      let x1 = -Infinity;
      let y1 = -Infinity;
      for (const end = pt + n; pt < end; pt++) {
        x0 = Math.min(x0, l.points[2 * pt]);
        x1 = Math.max(x1, l.points[2 * pt]);
        y0 = Math.min(y0, l.points[2 * pt + 1]);
        y1 = Math.max(y1, l.points[2 * pt + 1]);
      }
      if (n) out.push({ x0, y0, x1, y1 });
    },
    glyphs() {},
    text() {},
    image(l: ImageList, i: number) {
      out.push({ x0: l.x[i], y0: l.y[i], x1: l.x[i] + l.w[i], y1: l.y[i] + l.h[i] });
    },
    clip() {},
    unclip() {},
  });
  return out;
}

/** Whether a pass drew anything over a section's visible box. */
export function drawnIn(boxes: readonly Bounds[], s: SectionLayout): boolean {
  const b = s.placement.box;
  return boxes.some((r) => r.x1 >= b.x && r.x0 <= b.x + b.width && r.y1 >= b.y && r.y0 <= b.y + b.height);
}

/**
 * The sections other than the hovered one where the projected hover (`projected`, the viewer's emphasis) has positions
 * in a strip that draws a band: each must draw a hover mark (a linked view).
 */
export function linkedSections(
  frame: Frame,
  projected: { readonly parts: readonly { readonly system: string; readonly ranges: ArrayLike<number> }[] } | undefined,
  hovered: string | undefined,
): SectionLayout[] {
  if (!projected) return [];
  const out: SectionLayout[] = [];
  for (const s of frame.layout.sections) {
    if (s.section.address === hovered) continue;
    const st = s.placement.strips;
    let shown = false;
    for (let i = 0; i < st.n && !shown; i++) {
      // a strip with no band drawn (its tracks hidden) shows nothing to mark
      if (!(s.columns[i]?.height ?? 0)) continue;
      const a = st.axes[i];
      const part = projected.parts.find((p) => p.system === a.system);
      if (!part) continue;
      const vis = a.visible([]);
      for (let k = 0; k + 1 < part.ranges.length && !shown; k += 2)
        for (let m = 0; m + 1 < vis.length; m += 2)
          if (part.ranges[k] < vis[m + 1] && vis[m] < part.ranges[k + 1]) shown = true;
    }
    if (shown) out.push(s);
  }
  return out;
}

/** A 2D context that counts the calls that put ink down (the canvas backend's). */
function countingContext(): { ctx: Context2D; painted: () => number; reset: () => void } {
  let n = 0;
  const ink = new Set(['fill', 'stroke', 'fillRect', 'strokeRect', 'fillText', 'drawImage', 'putImageData']);
  const ctx = new Proxy({} as Record<string | symbol, unknown>, {
    get(t, k) {
      if (k in t) return t[k];
      if (k === 'createImageData') return (w: number, h: number) => ({ data: new Uint8ClampedArray(4 * w * h) });
      if (k === 'createLinearGradient') return () => ({ addColorStop() {} });
      return typeof k === 'string' && ink.has(k) ? () => void n++ : () => {};
    },
    set(t, k, v) {
      t[k] = v;
      return true;
    },
  }) as unknown as Context2D;
  return { ctx, painted: () => n, reset: () => void (n = 0) };
}

const counting = countingContext();
const canvas = (): CanvasLike => ({ width: 1, height: 1, getContext: () => counting.ctx });
const backend = new CanvasBackend(canvas(), canvas(), canvas());
const NoLists = new DrawLists();

/** The frame's hover pass alone, as draw lists both backends read (the other passes empty). */
function hoverOnly(lists: DrawLists): DrawLists {
  return Object.assign(Object.create(lists) as DrawLists, { passes: { ...NoLists.passes, hover: lists.passes.hover } });
}

/** What one view's sweep found: what is wrong, the probes, and the layer types and layouts it reached. */
export interface Swept {
  readonly wrong: string[];
  readonly probes: number;
  /** Probes in a rail's grab zone (the rail's hover is its own). */
  readonly rails: number;
  readonly types: string[];
  readonly layouts: string[];
}

/**
 * The hover contract over one settled view (_spec/README.md "Render", Emphasis, The hover contract): each probe's
 * pointer move answers with (1) a hit, (2) a hover mark over its own section (a bridge's: anywhere) and every linked
 * section with positions in view, drawn by the canvas's backend and in the SVG's hover group, (3) a card with words.
 */
export function sweepHover(v: SeqViewer, id: string, per = 4): Swept {
  const wrong: string[] = [];
  const types = new Set<string>();
  const layouts = new Set<string>();
  let rails = 0;
  const probes = probesOf(v.frame!, per);
  const words = describerFor({ ...v, layout: () => v.frame?.layout });
  backend.resize(v.frame!.layout.width, v.frame!.layout.height, 1);
  for (const p of probes) {
    const where = `${id} ${p.type} ${p.layer} ${p.mark}@${p.x.toFixed(0)},${p.y.toFixed(0)}`;
    // a rail's grab zone over a ruler takes the pointer: its hover is its own (no card, no data hover)
    if (v.hit(p.x, p.y)?.target === 'window') {
      rails++;
      continue;
    }
    v.input.move(p.x, p.y);
    const f = v.draw() ?? v.frame!;
    const h = v.hover.current;
    types.add(p.type);
    const s = f.layout.sections.find((x) => x.section.address === p.section);
    if (s) {
      layouts.add(s.layout.name);
      if (s.section.blocks) layouts.add('blocks');
      if (p.type === 'rows') layouts.add('rows');
    } else layouts.add('bridge');
    if (h.owner === undefined || h.layer === undefined) {
      wrong.push(`${where}: no hit`);
      v.input.leave();
      continue;
    }
    const boxes = passBoxes(f.lists.passes.hover);
    const own = f.layout.sections.find((x) => x.section.address === h.section);
    if (own ? !drawnIn(boxes, own) : !boxes.length) wrong.push(`${where}: no hover mark (hit ${h.layer})`);
    for (const l of linkedSections(f, v.emphasis.projected?.hover, h.section))
      if (!drawnIn(boxes, l)) wrong.push(`${where}: no hover mark in linked ${l.section.address} (hit ${h.layer})`);
    // both backends draw it: the canvas's paints, the SVG's hover group holds a mark
    if (boxes.length) {
      counting.reset();
      backend.drawOverlay(hoverOnly(f.lists));
      if (!counting.painted()) wrong.push(`${where}: the canvas painted no hover`);
      if (!/<g data-pass="hover">\s*<(rect|path|image|g)/.test(drawListsToSvg(hoverOnly(f.lists))))
        wrong.push(`${where}: the SVG drew no hover`);
    }
    const card = words.describe(h);
    if (!card || (!card.title && !card.rows.length)) wrong.push(`${where}: no card (hit ${h.layer})`);
    v.input.leave();
  }
  return { wrong, probes: probes.length, rails, types: [...types], layouts: [...layouts] };
}

/** What the sweep finds on this commit, a view at a time (each the probes' finding and their count): the next commit fixes them. */
export const Known: ReadonlyMap<string, readonly string[]> = new Map([
  ['across/across-bins', ['bridge peach-grape path: no card ×4']],
  ['across/across-chain', ['bridge grape-plum path: no card ×4', 'bridge peach-grape path: no card ×4']],
  ['across/across-fit', ['bridge peach-grape path: no card ×4']],
  ['across/across-pending', ['bridge peach-grape path: no card ×4']],
  ['across/across-self', ['bridge self path: no card ×4']],
  ['across/across-stub', ['bridge peach-grape path: no card ×4']],
  ['across/across-zoomed', ['bridge peach-grape path: no card ×4']],
  [
    'complex/text',
    [
      'ranges /views/main/complex/features/cdk2 image: no hover mark ×9',
      'ranges /views/main/complex/features/cyclin rect: no hit ×2',
      'ranges /views/main/complex/interface/cyclin rect: no hit ×1',
    ],
  ],
  ['complex/units-atoms', ['pairs /views/main/complex/contacts/contacts rect: no hover mark ×1']],
  [
    'genome/bed-over-rows',
    [
      'ranges /views/main/genome/features/features image: no hover mark ×15',
      'ranges /views/main/genome/features/features rect: no hit ×1',
    ],
  ],
  ['genome/inline-matrix', ['matrix /views/main/genome/levels/levels rect: no hover mark ×8']],
  ['genomes/vcf-genotypes', ['matrix /views/main/genotypes/gt/gt rect: no hover mark ×4']],
  ['p53/chromosome', ['ranges /views/main/chr17/genes/ranges image: no hover mark ×3']],
  [
    'p53/text',
    [
      'ranges /views/main/protein/domains/ranges image: no hit ×7',
      'ranges /views/main/protein/domains/ranges image: no hover mark ×23',
    ],
  ],
  ['p53/text-centered', ['ranges /views/main/protein/domains/ranges image: no hover mark ×42']],
  [
    'p53/text-under-tracks',
    [
      'ranges /views/main/text/domains/ranges image: no hit ×7',
      'ranges /views/main/text/domains/ranges image: no hover mark ×23',
    ],
  ],
  ['p53/variants', ['variants /views/main/gene/variants/variants image: no hover mark ×3']],
  [
    'plasmid/scale',
    [
      'ranges /views/main/detail/minus/features image: no hover mark ×2',
      'ranges /views/main/detail/plus/features image: no hover mark ×2',
      'ranges /views/main/ring/minus/features image: no hover mark ×1',
      'ranges /views/main/ring/plus/features image: no hover mark ×1',
    ],
  ],
  ['regions/matrix-methylation', ['matrix /views/main/chr/methylation/methylation image: no hover mark ×3']],
  ['regions/matrix-ring', ['matrix /views/main/chr/gt/gt image: no hover mark ×3']],
  ['regions/matrix-samples', ['matrix /views/main/chr/gt/gt image: no hover mark ×3']],
]);

/** A sweep's findings summed: each probe's finding with its place dropped, and how many probes found it. */
export function summed(wrong: readonly string[]): string[] {
  const n = new Map<string, number>();
  for (const w of wrong) {
    const k = w
      .replace(/^\S+ /, '')
      .replace(/@-?\d+,-?\d+/, '')
      .replace(/ \(hit [^)]*\)/, '');
    n.set(k, (n.get(k) ?? 0) + 1);
  }
  return [...n].map(([k, c]) => `${k} ×${c}`).sort();
}
