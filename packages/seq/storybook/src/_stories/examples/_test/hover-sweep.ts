/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Location, rangesLocation } from '@molstar/seqx-core/query';
import { DrawLists, drawListsToSvg, type Frame, type SectionLayout, WarpKind } from '@molstar/seqx-core/render';
import { projectOntoView } from '@molstar/seqx-core/view';
import { CanvasBackend, type CanvasLike, type Context2D } from '@molstar/seqx-viewer/canvas';
import { describerFor } from '@molstar/seqx-viewer/describe';
import { legendRows } from '@molstar/seqx-viewer/legends';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { type Bounds, passBoxes, probesOf } from './hover-probes.js';

// The hover contract's sweep (_spec/README.md "Render", Emphasis, The hover contract): each probe on a drawn mark
// (`probesOf`) is hovered by the pointer, and the viewer's answer checked: the hit is the probe's layer (or an upper
// layer over it, recorded), a hover mark over the probe's point in its own section and over the hover's positions in
// each section they are carried to (computed here from the hover's extent, independently of the viewer's emphasis),
// in the draw lists, painted by the canvas's backend and drawn in the SVG's hover group alike, and a card with words.

/** A page box the hover must reach: the probe's point (its own section), a linked strip's span of the hover. */
type Region = Bounds;

/** Whether one of the boxes reaches a region (a point: within the outline's offset and width, 3 px). */
const reaches = (boxes: readonly Bounds[], r: Region, slop = 3) =>
  boxes.some((b) => b.x1 >= r.x0 - slop && b.x0 <= r.x1 + slop && b.y1 >= r.y0 - slop && b.y0 <= r.y1 + slop);

/**
 * The regions of a section where a hover over `loc` must draw: each translated strip's band column over the
 * positions of `loc` in view (a polar strip's: the section's box), none when no position is in view.
 */
function linkedRegions(s: SectionLayout, loc: Location): Region[] {
  const out: Region[] = [];
  const st = s.placement.strips;
  for (let i = 0; i < st.n; i++) {
    const col = s.columns[i];
    if (!col?.height) continue;
    const a = st.axes[i];
    const part = loc.parts.find((p) => p.system === a.system);
    if (!part) continue;
    const spans: number[] = [];
    for (let k = 0; k + 1 < part.ranges.length; k += 2) a.spans(part.ranges[k], part.ranges[k + 1], spans);
    if (!spans.length) continue;
    if (st.warp[i] === WarpKind.polar) {
      const b = s.placement.box;
      out.push({ x0: b.x, y0: b.y, x1: b.x + b.width, y1: b.y + b.height });
      continue;
    }
    for (let k = 0; k + 2 < spans.length; k += 3)
      out.push({
        x0: spans[k] + st.dx[i],
        x1: spans[k + 1] + st.dx[i],
        y0: col.y + st.dy[i],
        y1: col.y + col.height + st.dy[i],
      });
  }
  return out;
}

// ---- the two backends, read back as boxes ----

/** A 2D context that keeps the box of each thing painted (a filled or stroked path, a rect), in CSS px. */
function recordingContext() {
  const boxes: Bounds[] = [];
  let path: Bounds | undefined;
  const add = (x: number, y: number) => {
    path ??= { x0: x, y0: y, x1: x, y1: y };
    path.x0 = Math.min(path.x0, x);
    path.y0 = Math.min(path.y0, y);
    path.x1 = Math.max(path.x1, x);
    path.y1 = Math.max(path.y1, y);
  };
  const paint = () => {
    if (path) boxes.push({ ...path });
  };
  const own: Record<string, unknown> = {
    beginPath: () => void (path = undefined),
    moveTo: add,
    lineTo: add,
    // the curve's hull (its controls and its end), as the lists' and the SVG's boxes: a chord bulges between its ends
    bezierCurveTo: (a: number, b: number, c: number, d: number, x: number, y: number) => {
      add(a, b);
      add(c, d);
      add(x, y);
    },
    rect: (x: number, y: number, w: number, h: number) => {
      add(x, y);
      add(x + w, y + h);
    },
    fill: paint,
    stroke: paint,
    fillRect: (x: number, y: number, w: number, h: number) => boxes.push({ x0: x, y0: y, x1: x + w, y1: y + h }),
    strokeRect: (x: number, y: number, w: number, h: number) => boxes.push({ x0: x, y0: y, x1: x + w, y1: y + h }),
    drawImage: (_s: unknown, x: number, y: number, w: number, h: number) =>
      boxes.push({ x0: x, y0: y, x1: x + w, y1: y + h }),
    createImageData: (w: number, h: number) => ({ data: new Uint8ClampedArray(4 * w * h) }),
    createLinearGradient: () => ({ addColorStop() {} }),
  };
  const ctx = new Proxy({} as Record<string | symbol, unknown>, {
    get(t, k) {
      if (k in t) return t[k];
      if (typeof k === 'string' && k in own) return own[k];
      return () => {};
    },
    set(t, k, v) {
      t[k] = v;
      return true;
    },
  }) as unknown as Context2D;
  return { ctx, boxes };
}

const recording = recordingContext();
const canvas = (): CanvasLike => ({ width: 1, height: 1, getContext: () => recording.ctx });
const backend = new CanvasBackend(canvas(), canvas(), canvas());
const NoLists = new DrawLists();

/** The frame's hover pass alone, as draw lists both backends read (the other passes empty). */
function hoverOnly(lists: DrawLists): DrawLists {
  return Object.assign(Object.create(lists) as DrawLists, { passes: { ...NoLists.passes, hover: lists.passes.hover } });
}

/** What the canvas's backend painted of the hover pass, as boxes. */
function canvasBoxes(lists: DrawLists): Bounds[] {
  recording.boxes.length = 0;
  backend.drawOverlay(hoverOnly(lists));
  return [...recording.boxes];
}

/** What the SVG's hover group draws, as boxes: its rects and its paths' points (a clip's group draws nothing). */
function svgBoxes(lists: DrawLists): Bounds[] {
  const svg = drawListsToSvg(hoverOnly(lists));
  const at = svg.indexOf('<g data-pass="hover">');
  if (at < 0) return [];
  const body = svg.slice(at);
  const out: Bounds[] = [];
  const num = (tag: string, a: string) => Number(new RegExp(`\\b${a}="(-?[\\d.e-]+)"`).exec(tag)?.[1] ?? Number.NaN);
  for (const m of body.matchAll(/<rect\b[^>]*>/g)) {
    const [x, y, w, h] = ['x', 'y', 'width', 'height'].map((a) => num(m[0], a));
    if ([x, y, w, h].every(Number.isFinite)) out.push({ x0: x, y0: y, x1: x + w, y1: y + h });
  }
  for (const m of body.matchAll(/<path\b[^>]*\bd="([^"]*)"/g)) {
    const n = (m[1].match(/-?[\d.]+(?:e-?\d+)?/g) ?? []).map(Number);
    if (n.length < 2) continue;
    const b = { x0: Infinity, y0: Infinity, x1: -Infinity, y1: -Infinity };
    for (let k = 0; k + 1 < n.length; k += 2) {
      b.x0 = Math.min(b.x0, n[k]);
      b.x1 = Math.max(b.x1, n[k]);
      b.y0 = Math.min(b.y0, n[k + 1]);
      b.y1 = Math.max(b.y1, n[k + 1]);
    }
    out.push(b);
  }
  return out;
}

// ---- the sweep ----

/** What one view's sweep found: what is wrong, the probes, and the layer types and layouts it reached. */
export interface Swept {
  readonly wrong: string[];
  readonly probes: number;
  /** Probes in a rail's grab zone (the rail's hover is its own). */
  readonly rails: number;
  /** The layer types that answered a probe on their own marks (by the hit, never the probe). */
  readonly types: string[];
  readonly layouts: string[];
  /** Linked sections where the hover was drawn over its positions: a projection dropped leaves none. */
  readonly linked: number;
  /** Whether a hover of this view has positions in view in another section (computed from the model and view). */
  readonly links: boolean;
  /** The legend entries hovered items read as, each with a probe point of one (`key`, NUL, `value` → point). */
  readonly entries: ReadonlyMap<string, Region>;
  /** Every card the probes showed, in probe order, a line of JSON each: what the card providers' parity reads. */
  readonly cards: string;
}

/** The type of the layer a key names in a frame ('bridge' for a bridge's). */
function typeOf(f: Frame, key: string): string | undefined {
  for (const s of f.layout.sections)
    for (const l of [...s.lays.flat(), ...s.sectionLays]) if (l?.key === key) return l.type.name;
  return f.layout.bridges?.list.some((b) => key.endsWith(`/${b.id}`) || key === b.id) ? 'bridge' : undefined;
}

/**
 * The hover contract over one settled view: each probe's pointer move answers with (1) the probe's layer's hit (or an
 * upper layer's over the same point, recorded: the probe's type must answer somewhere in the view), (2) a hover mark
 * over the probe's point and over the hover's positions in each linked section, in the lists, the canvas and the SVG,
 * (3) a card with words.
 */
export function sweepHover(v: SeqViewer, id: string, per = 4): Swept {
  const wrong: string[] = [];
  const answered = new Set<string>();
  const covered = new Map<string, string>();
  const layouts = new Set<string>();
  const entries = new Map<string, Region>();
  let rails = 0;
  let linked = 0;
  let links = false;
  let cards = '';
  const probes = probesOf(v.frame!, per);
  const words = describerFor({ ...v, layout: () => v.frame?.layout });
  backend.resize(v.frame!.layout.width, v.frame!.layout.height, 1);
  const options = { own: true, close: v.config.hover.closeCodons, linked: v.config.hover.projectLinked };
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
    // (1) the probe's own layer, else an upper one over its mark in the same section
    const hitType = typeOf(f, h.layer) ?? (s ? undefined : 'bridge');
    if (h.section !== p.section) wrong.push(`${where}: hit in ${h.section}, not the probe's section`);
    if (h.layer === p.layer) answered.add(p.type);
    else if (!covered.has(p.type)) covered.set(p.type, `${where}: answered by ${h.layer}`);
    if (hitType) answered.add(hitType);
    // (2) over the probe's point, and over the hover's positions in every other section they reach
    const point: Region = { x0: p.x, y0: p.y, x1: p.x, y1: p.y };
    const loc =
      h.extent ?? (Number.isNaN(h.position) ? undefined : rangesLocation(h.system!, [h.position, h.position + 1]));
    const shown = loc && projectOntoView(v.view.model.value, v.view.resolved.value!, loc, options);
    // its own section: over the hovered extent's columns in the probe's band (an item's, a position's), else the point
    const ownRegions =
      s && p.band && shown
        ? linkedRegions(s, shown).map((r) => ({ x0: r.x0, x1: r.x1, y0: p.band!.y0, y1: p.band!.y1 }))
        : [];
    const needs: { what: string; regions: Region[] }[] = [
      { what: 'the hovered extent', regions: ownRegions.length ? ownRegions : [point] },
    ];
    for (const o of f.layout.sections) {
      if (o.section.address === h.section || !shown) continue;
      const rs = linkedRegions(o, shown);
      if (!rs.length) continue;
      links = true;
      needs.push({ what: `linked ${o.section.address}`, regions: rs });
    }
    const backends = {
      lists: passBoxes(f.lists.passes.hover),
      canvas: canvasBoxes(f.lists),
      svg: svgBoxes(f.lists),
    };
    for (const [name, boxes] of Object.entries(backends))
      for (const n of needs) {
        const ok = n.regions.some((r) => reaches(boxes, r));
        if (!ok) wrong.push(`${where}: ${name} drew no hover over ${n.what} (hit ${h.layer})`);
        else if (name === 'lists' && n.what.startsWith('linked')) linked++;
      }
    // (3) a card with words
    const card = words.describe(h);
    if (!card || (!card.title && !card.rows.length)) wrong.push(`${where}: no card (hit ${h.layer})`);
    cards += `${JSON.stringify(card ?? null)}\n`;
    // the legend entry the item reads as (outlined): the way back is swept by `sweepLegends`
    const entry = v.legends.hover.value;
    if (entry?.value !== undefined && !entries.has(`${entry.key}\u0000${entry.value}`))
      entries.set(`${entry.key}\u0000${entry.value}`, point);
    v.input.leave();
  }
  // a type covered by an upper layer everywhere answers no hover of its own
  for (const [type, first] of covered) if (!answered.has(type)) wrong.push(`${first}; ${type} answered nowhere`);
  return {
    wrong,
    probes: probes.length,
    rails,
    types: [...answered],
    layouts: [...layouts],
    linked,
    links,
    entries,
    cards,
  };
}

/**
 * Legends both ways (_spec/README.md "Legends"; the hover contract): each entry a hovered item reads as has an action,
 * and hovering it lights its items (a highlight over that item's point and over every item it names in view); an entry without an action (the default
 * cursor) lights nothing when hovered.
 */
export function sweepLegends(v: SeqViewer, id: string, entries: ReadonlyMap<string, Region>): string[] {
  const wrong: string[] = [];
  const rows = legendRows(v);
  const lit = (layer: string, value: string) => {
    void v.execute('legend.hover', { address: layer, entry: value });
    const boxes = passBoxes((v.draw() ?? v.frame!).lists.passes.highlight);
    void v.execute('legend.hover', {});
    v.draw();
    return boxes;
  };
  for (const [e, point] of entries) {
    const [key, value] = e.split('\u0000');
    const row = rows.find((r) => r.key === key);
    const entry = row?.entries.find((x) => x.value === value);
    if (!row || !entry?.action) {
      wrong.push(`${id} legend ${key} ${value}: an item reads as it, but the entry has no action`);
      continue;
    }
    const boxes = lit(row.layer, value);
    if (!reaches(boxes, point)) wrong.push(`${id} legend ${key} ${value}: its item not lit`);
    // and every item it names in view, over its positions (the items drawn: `entryItems`)
    const items = v.legends.entryItems(row.layer, value);
    let unlit = 0;
    for (const part of items?.location.parts ?? [])
      for (let k = 0; k + 1 < part.ranges.length; k += 2) {
        const one = rangesLocation(part.system, [part.ranges[k], part.ranges[k + 1]]);
        const regions = v.frame!.layout.sections.flatMap((s) => linkedRegions(s, one));
        if (regions.length && !regions.some((r) => reaches(boxes, r))) unlit++;
      }
    if (unlit) wrong.push(`${id} legend ${key} ${value}: ${unlit} of its items in view not lit`);
  }
  for (const row of rows)
    for (const entry of row.entries)
      if (!entry.action && entry.value !== undefined && lit(row.layer, entry.value).length)
        wrong.push(`${id} legend ${row.key} ${entry.value}: no action (the default cursor), yet it lights`);
  return wrong;
}
