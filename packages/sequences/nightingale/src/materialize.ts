/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `materialize` (seqspec:libraries#nightingale): a resolved view, resolved with `NightingaleSupported`, as
// Nightingale's elements and their data. Platform-free: it imports no Nightingale. What Nightingale cannot draw is
// reported under the format's codes, the message naming the limit; `path` is the node's ref (a resolved view keeps no
// paths).

import {
  Arrangement,
  Color,
  endsOf,
  KernelColorSchemes,
  Layout,
  layerPalette,
  letterColors,
  letterMapTable,
  makePaletteFunction,
  type Palette,
  paletteValue,
  type ResolvedView,
  type StreamedData,
  StreamProp,
  type SvsRegistries,
  schemeTheme,
  splitRows,
  type ThemeRef,
} from "@molstar/seq-core/format";
import {
  type Annotation,
  Annotations,
  type ComposedPath,
  CoordinateSystemKind,
  Coords,
  type Issue,
  type SeqSet,
  threeLetter,
} from "@molstar/seq-core/model";
import { type AxisBlock, type AxisMeta, AxisRuns, AxisSink, axesOf, firstFrom, indexAt, sectionRuns } from "./axis.ts";
import { type ElementDesc, el } from "./element.ts";
import { nightingaleIssue } from "./issues.ts";
import type { StreamedPart, StreamedValues } from "./stream.ts";

const [FALLBACK, THEME, GLYPH] = ["svs.render.fallback", "svs.render.theme", "svs.render.glyph"];
/** The letter dictionaries nightingale-msa has a scheme for, and the scheme. */
const SCHEMES: Readonly<Record<string, string>> = Object.fromEntries(
  "Zappo:zappo Taylor:taylor Hydrophobicity:hydro HelixPropensity:helix StrandPropensity:strand TurnPropensity:turn BuriedIndex:buried Nucleotide:nucleotide PurinePyrimidine:purine"
    .split(" ")
    .map((x) => x.split(":")),
);
/** The kernel's color schemes nightingale-msa has its own for: drawn as those (`SCHEMES`, Clustal X's `clustal2`). */
const MSA_NAMED: Readonly<Record<string, ThemeRef>> = {
  clustal_x: { theme: "rules", params: { rules: [] }, table: "clustal_x" },
  ...Object.fromEntries(
    "zappo:Zappo taylor:Taylor hydrophobicity:Hydrophobicity nucleotide:Nucleotide purine_pyrimidine:PurinePyrimidine"
      .split(" ")
      .map((x) => x.split(":"))
      .map(([s, d]) => [s, { theme: "letter", params: { colors: d } }]),
  ),
};
const kernelSchemes = { get: (_: string, n: string) => KernelColorSchemes.find((e) => e.name === n) } as SvsRegistries;
/**
 * A `scheme` theme as the theme it stands for, the light variant (seqspec:interaction#color-schemes-at-runtime); others
 * as they are.
 */
function unscheme(r: ThemeRef | null): ThemeRef | null {
  if (r?.theme !== "scheme") return r;
  const p = r.params ?? {};
  const own = typeof p.name === "string" && Object.keys(p.overrides ?? {}).length === 0 ? MSA_NAMED[p.name] : undefined;
  return own ?? schemeTheme(kernelSchemes, p, { dark: false }) ?? r;
}
/** Element heights in pixels (hints Nightingale takes as its defaults). */
const HEIGHT: Readonly<Record<string, number>> = {
  ...{ letters: 20, colors: 20, ranges: 24, points: 18, variants: 200, secondary: 18, values: 60, color: 14 },
  ...{ ruler: 40, pairs: 60, logo: 80, blocks: 18, tile: 18 },
};
const SECONDARY: Readonly<Record<string, string>> = { H: "helix", G: "helix", I: "helix", E: "strand", B: "strand" };
const COMMON = { "highlight-event": "onmouseover", "use-ctrl-to-zoom": "" };
/** The data side's gutters in pixels (each element's `margin-left`, `margin-right`): Nightingale's default is 10. */
const GUTTER = { side: 14, variants: 30 } as const;
/** The custom property a pairs layer's one color (`uniform`) is given in, for `mount`'s style over Nightingale's. */
export const PAIR_COLOR = "--svs-pair-color";
/** The class of our element stacking a row's elements in one place (`mount`'s style): letters over their bar. */
export const STACK = "svs-stack";
/** The class of a nightingale-sequence whose own axis `mount` hides (only a ruler draws positions). */
export const AXISLESS = "svs-axisless";
/** A theme as the document named it: a named rule table arrives expanded as `rules`, its name in `table` (VW4). */
const themeName = (r: ThemeRef) => r.table ?? r.theme;
/** Letters' codes as text: latin1 is the ASCII the codes are, with no per-character spread. */
const latin1 = new TextDecoder("latin1");
/** A color number that is no color (an `Int32Array` of colors starts filled with it). */
const NO = -1;

/** nightingale-colored-sequence's ground, as nightingale-msa's own placeholders: a letter left uncolored. */
const GROUND = 0xffffff;
/** The distinct letters of `seq` (upper case, as the component reads them) a dictionary does not color. */
function lettersOutside(seq: string, dict: Readonly<Record<string, unknown>>): string[] {
  const seen = new Uint8Array(65536);
  const out: string[] = [];
  for (let i = 0; i < seq.length; i++) {
    const c = seq.charCodeAt(i);
    if (seen[c]) continue;
    seen[c] = 1;
    const k = String.fromCharCode(c).toUpperCase();
    if (!(k in dict) && !out.includes(k) && k !== "," && k !== ":") out.push(k);
  }
  return out;
}

/** One end of a pair as placed: its manager, axis position, side and placement. */
interface PairEnd {
  readonly manager: number;
  readonly at: number;
  readonly side: number;
  readonly o: On;
  readonly row: number;
}

/**
 * Nightingale's limit, in positions, of a section (seqspec:libraries#nightingale): its tracks take an element or a
 * datum per position (a sequence's letters, a line graph's points), so a longer section (a streamed chromosome is 83
 * to 249 Mb) is left out, reported (`svs.render.omitted`), while the view's other sections render.
 */
export const NIGHTINGALE_MAX_POSITIONS = 1_000_000;

/** A section's positions: its blocks' pieces summed. */
function sectionLength(view: ResolvedView, s: number): number {
  const [S, B] = [view.sections, view.blocks];
  let n = 0;
  for (let b = S.blockStart[s]; b < S.blockStart[s] + S.blockCount[s]; b++) n += B.stop[b] - B.start[b];
  return n;
}

/** The attribute of a track some of whose features' labels go inside them (their `svs` record's `inside`): `mount` lays them over it. */
export const LABELED = "data-svs-labels";

/** Whether two features are drawn alike: everything but their accession and `svs`, field by field. */
function alike(a: Feature, b: Feature): boolean {
  const [x, y] = [a.locations, b.locations];
  if (x?.length !== y?.length) return false;
  for (let i = 0; x && y && i < x.length; i++) {
    const [u, v] = [x[i].fragments, y[i].fragments];
    if (u.length !== v.length) return false;
    for (let k = 0; k < u.length; k++) if (u[k].start !== v[k].start || u[k].end !== v[k].end) return false;
  }
  return (
    a.start === b.start &&
    a.end === b.end &&
    a.color === b.color &&
    a.fill === b.fill &&
    a.shape === b.shape &&
    a.opacity === b.opacity &&
    a.tooltipContent === b.tooltipContent
  );
}

/**
 * Features as a track: canvas from 1,000 features, the given attributes. Features drawn alike (merged copies' items at
 * the same place, one item in two tracks of a collapsed group) are drawn once, so they do not stack into lanes.
 */
function featureTrack(features: readonly Feature[], attributes: Record<string, string>): ElementDesc {
  const byStart = new Map<number, Feature[]>();
  const data: Feature[] = [];
  for (const f of features) {
    const same = byStart.get(f.start) ?? byStart.set(f.start, []).get(f.start);
    if (same?.some((g) => alike(f, g))) continue;
    same?.push(f);
    data.push(f);
  }
  const labeled: Record<string, string> = data.some((f) => f.svs.inside) ? { [LABELED]: "" } : {};
  // Lanes keep their features readable: Nightingale divides the height among them (a 2 px gap each), so the track
  // grows with its lanes, up to LANES of them, rather than squeezing them.
  const lanes = attributes.layout === "non-overlapping" ? lanesOf(data) : 1;
  const height = String(Math.max(Number(attributes.height), lanes * LANE));
  const tag = data.length >= 1000 ? "nightingale-track-canvas" : "nightingale-track";
  return el(tag, { ...attributes, height, ...labeled }, { data });
}

/** A lane's height in a track of several (its feature 12 px and Nightingale's gap), the most lanes given room. */
const [LANE, LANES] = [14, 10];
/** A run of a lane's disjoint intervals, sorted by start; a lane is its chunks in order, each at most `2 * CHUNK`. */
interface Chunk {
  readonly starts: number[];
  readonly ends: number[];
}
const CHUNK = 256;
/** What `lanesOf` did: binary search steps, and intervals moved by inserts (tests bound both). */
export interface LaneStats {
  probes: number;
  moves: number;
}
/**
 * The lanes Nightingale's non-overlapping layout puts features in (first fit in data order, ends inclusive), at most
 * `LANES` (a feature that fits none returns at once). Per lane its disjoint intervals sorted by start, in chunks: a
 * feature beyond a lane's hull fits in O(1), else a binary search finds its two neighbors, and an insert moves at most
 * a chunk, so any order is O(n log n) and sorted input O(n).
 */
export function lanesOf(features: readonly { start: number; end: number }[], stats?: LaneStats): number {
  const lanes: Chunk[][] = [];
  for (const f of features) {
    let k = 0;
    while (k < lanes.length && !fits(lanes[k], f.start, f.end, stats)) k++;
    if (k < lanes.length) continue;
    if (k >= LANES) return LANES;
    lanes.push([{ starts: [f.start], ends: [f.end] }]);
  }
  return lanes.length;
}
/** Puts [start, end] in the lane when it overlaps none of its intervals; a chunk past `2 * CHUNK` splits in two. */
function fits(lane: Chunk[], start: number, end: number, stats?: LaneStats): boolean {
  const last = lane[lane.length - 1];
  let [i, at] = [0, 0];
  if (start > last.ends[last.ends.length - 1]) [i, at] = [lane.length - 1, last.starts.length];
  else if (end >= lane[0].starts[0]) {
    let [lo, hi] = [0, lane.length - 1];
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (stats) stats.probes++;
      if (lane[mid].starts[0] <= start) lo = mid;
      else hi = mid - 1;
    }
    const { starts, ends } = lane[lo];
    let [a, b] = [0, starts.length];
    while (a < b) {
      const mid = (a + b) >> 1;
      if (stats) stats.probes++;
      if (starts[mid] <= start) a = mid + 1;
      else b = mid;
    }
    const next = a < starts.length ? starts[a] : (lane[lo + 1]?.starts[0] ?? Number.POSITIVE_INFINITY);
    if ((a > 0 && ends[a - 1] >= start) || next <= end) return false;
    [i, at] = [lo, a];
  }
  const c = lane[i];
  if (stats) stats.moves += c.starts.length - at;
  c.starts.splice(at, 0, start);
  c.ends.splice(at, 0, end);
  if (c.starts.length > 2 * CHUNK)
    lane.splice(i + 1, 0, { starts: c.starts.splice(CHUNK), ends: c.ends.splice(CHUNK) });
  return true;
}

/** A feature's outline: its fill darkened, so features of one color that touch (lanes, neighbors) stay apart. */
const outline = (fill: string | undefined) => {
  const c = Color.decode(fill);
  return c === undefined ? undefined : Color.toHex(Color.interpolate(c, 0x000000 as Color, 0.35));
};

/** The Nightingale elements under the given descriptions, in order. */
function nightingaleElements(ds: readonly ElementDesc[], out: ElementDesc[] = []): ElementDesc[] {
  for (const d of ds) {
    if (d.tag.startsWith("nightingale-")) out.push(d);
    else nightingaleElements(d.children, out);
  }
  return out;
}

interface Feature {
  accession: string;
  start: number;
  end: number;
  locations?: { fragments: { start: number; end: number }[] }[];
  color?: string;
  fill?: string;
  shape?: string;
  opacity?: number;
  tooltipContent?: string;
  /**
   * For mount: the view's layer, the annotation and its item (event mapping); its label (a block's name, layer -1; a
   * range's `label`; a residue's letter over its codon), the hover hint's unless `hint` says more, and laid inside the
   * feature where it fits when `inside`.
   */
  svs: { layer: number; annotation: number; item: number; label?: string; inside?: boolean; hint?: string };
}
/** Where a placement puts its indices: the axis block and the path onto it (`read` 1: its map from → to; -1: back). */
interface At {
  b: AxisBlock;
  path: ComposedPath | null;
  read: number;
}
/** One placement's rows on one axis block: the rows the filter keeps. */
interface On extends At {
  ann: Annotation;
  side: 0 | 1;
  rows: number[];
  /** The end's slot (dense values: its array). */
  slot: number;
}
/** A layer's colors: a color number per row or value (`NO`: none). */
interface Colorer {
  row(a: Annotation, r: number, b: AxisBlock): number;
  value(x: number): number;
}
/** Visits a path map's runs onto the sink (`Coords`' visitor shape). */
function toSink(a: number, z: number, sink: AxisSink): void {
  sink.add(a, z);
}
/** The path's `from` element whose image is the block index `x` (the first copy); -1: none. */
function sourceAt(p: ComposedPath, read: number, x: number): number {
  return read < 0 ? Coords.toAt(p.map, x, 0) : firstFrom(p.map, x);
}

export function materialize(view: ResolvedView, set: SeqSet): { root: ElementDesc; issues: Issue[] } {
  const issues: Issue[] = [];
  const seen = new Set<string>();
  const name = (id: number) => (id >= 0 ? set.strings.get(id) : "");
  const report = (code: string, ref: number, message: string) => {
    const key = `${code}|${ref}|${message}`;
    if (!seen.has(key) && seen.add(key)) issues.push(nightingaleIssue(code, name(ref), `Nightingale: ${message}`));
  };
  const [S, G, T, L, P] = [view.sections, view.groups, view.tracks, view.layers, view.placements];
  /** What a streamed layer takes of the fallback, once a layer (`svs.stream.fallback`). */
  const streamNote = (l: number, message: string) => report("svs.stream.fallback", L.ref[l], message);
  const typeOf = (l: number) => name(L.type[l]);
  const lref = (l: number) => L.ref[l];
  const params = (l: number) => L.params[l] as Record<string, unknown>;
  const mapOf = (l: number) => letterMapTable(params(l).letter_map as never);
  /** A layer's theme of a channel (`color` first, then `glyph`), as the resolution left it. */
  const themeOf = (l: number, channel: number) => unscheme((L.themes[l]?.[channel] ?? null) as ThemeRef | null);

  if (view.layout === Layout.text)
    report(FALLBACK, view.ref, "no wrapped component: the text layout is drawn on one unwrapped axis");

  /** The one converter onto the axis: block-system runs into `runs` (axis.ts), the placements' and the highlight's. */
  const sink = new AxisSink();
  const runs = sink.runs;
  /**
   * The indices [s, t) of a placement's system onto the axis, as sorted joined runs in `runs`; a boundary (s === t) as
   * its two neighbors. Whether anything is shown.
   */
  const into = (o: At, a: number, z: number) => {
    if (!o.path) sink.add(a, z);
    else if (o.read > 0) Coords.forEachTo(o.path.map, a, z, toSink, sink);
    else Coords.forEachFrom(o.path.map, a, z, toSink, sink);
  };
  const place = (o: At, s: number, t: number): boolean => {
    runs.count = 0;
    runs.sorted = true;
    sink.block = o.b;
    if (s !== t) into(o, s, t);
    else {
      if (s > 0) into(o, s - 1, s);
      into(o, s, s + 1);
    }
    AxisRuns.finish(runs);
    return runs.count > 0;
  };
  /** The first axis position (0-based) of the index `i`; -1: not shown. */
  const firstAt = (o: At, i: number) => (place(o, i, i + 1) ? runs.start[0] : -1);
  /**
   * The runs a feature draws as fragments: each run on its own block (a block's origin seam or reversed piece parts
   * them, NG-04) or through a codon step; through another path their hull (a gap does not split).
   */
  const fragmentsOf = (o: At): { start: number; end: number }[] => {
    const out: { start: number; end: number }[] = [];
    if (!o.path || o.path.heldSpans)
      for (let k = 0; k < runs.count; k++) out.push({ start: runs.start[k] + 1, end: runs.stop[k] });
    else out.push({ start: runs.start[0] + 1, end: runs.stop[runs.count - 1] });
    return out;
  };
  const placementOf = (p: number, axis: AxisMeta): At | null => {
    const b = axis.blocks.find((x) => x.block === P.block[p]);
    const path = P.path[p] >= 0 ? view.paths[P.path[p]] : null;
    return b ? { b, path, read: path && Coords.sameSystem(path.map.from, path.from) ? 1 : -1 } : null;
  };

  /** Every placement of a layer on an axis, with its rows. */
  const placements = (l: number, axis: AxisMeta): On[] => {
    const out: On[] = [];
    for (let p = L.placementStart[l]; p < L.placementStart[l] + L.placementCount[l]; p++) {
      const at = placementOf(p, axis);
      const ann = set.annotations[P.annotation[p]];
      const end = ann && endsOf(ann).find((e) => e.slot === P.slot[p]);
      if (!at || !ann || !end) continue;
      // A values layer's bins are fetched by `mount` (stream.ts); other types draw none.
      if (ann.streamed && typeOf(l) !== "values")
        streamNote(l, `streamed ${typeOf(l)} is not drawn: mount fetches values only`);
      let k = L.dataStart[l];
      while (k < L.dataStart[l] + L.dataCount[l] && view.data[k] !== ann.id) k++;
      const kept = view.kept[k];
      const rows = end.rows.filter((r) => !kept || kept[r]);
      out.push({ ...at, ann, side: end.side, rows, slot: end.slot });
    }
    return out;
  };

  /** Hex strings of colors, made once each (a categorical palette's are made up front, a continuous one's as it goes). */
  const hexes = new Map<number, string>();
  const hexOf = (c: number | undefined): string | undefined => {
    if (c === undefined || c < 0) return undefined;
    const h = hexes.get(c) ?? Color.toHex(c);
    hexes.set(c, h);
    return h;
  };
  /**
   * The layer's color per row (its theme), one per layer: its palette is built once over all the layer's data
   * (`layerPalette`), so a value has the viewer's color in every block, section and filter.
   */
  const colorers = new Map<number, Colorer>();
  const colorer = (l: number): Colorer => {
    const c = colorers.get(l) ?? colorerOf(l);
    colorers.set(l, c);
    return c;
  };
  const colorerOf = (l: number): Colorer => {
    const ref = themeOf(l, 0);
    const p = (ref?.params ?? {}) as Record<string, unknown>;
    const none: Colorer = { row: () => NO, value: () => NO };
    if (!ref) return none;
    if (ref.theme === "uniform") {
      const c = Color.decodePrimary(String(p.color));
      return { row: () => c, value: () => c };
    }
    if (ref.theme === "field_color") {
      const field = String(p.field ?? "color");
      const of = (x: unknown) => Color.decodeFirst(String(x)) ?? NO;
      return { row: (a, r) => of(Annotations.rowValue(a, set, r, field)), value: of };
    }
    if (ref.theme === "field_palette") {
      // The field's default (no `category` column: the item) is the kernel's, shared by every viewer.
      const field = typeof p.field === "string" ? p.field : null;
      const data = Array.from(
        { length: L.dataCount[l] },
        (_, k) => set.annotations[view.data[L.dataStart[l] + k]] as Annotation,
      );
      const f = layerPalette(set, data, (p.palette ?? { kind: "categorical" }) as Palette, field, [
        L.domainMin[l],
        L.domainMax[l],
      ]);
      for (const key of f.keys ?? []) hexOf(f(key) ?? NO);
      return { row: (a, r) => f(paletteValue(a, set, r, field)) ?? NO, value: (x) => f(x) ?? NO };
    }
    if (ref.theme === "by_unit" || ref.theme === "by_copy_group") {
      const units = Array.from({ length: set.units.count }, (_, u) => u).filter((u) => set.units.alive[u]);
      const key = (u: number) =>
        name(ref.theme === "by_unit" ? set.units.key[u] : set.units.key[set.groups.first[set.units.group[u]] ?? u]);
      const f = makePaletteFunction((p.palette ?? { kind: "categorical" }) as Palette, units.map(key));
      return { row: (_a, _r, b) => (b && b.unit >= 0 ? (f(key(b.unit)) ?? NO) : NO), value: () => NO };
    }
    report(
      THEME,
      lref(l),
      `no slot for the "${themeName(ref)}" theme on ${typeOf(l)}: drawn in the component's colors`,
    );
    return none;
  };

  const element = (tag: string, l: number, kind: string, attributes: Record<string, string> = {}, props = {}) =>
    el(tag, { ...COMMON, height: String(HEIGHT[kind]), "data-svs-layer": String(l), ...attributes }, props);
  /** The letters of every block on the axis, as text: one byte each, decoded once. */
  const lettersOf = (axis: AxisMeta, row: number | null, map: Uint8Array | null = null) => {
    const out = new Uint8Array(axis.blocks.reduce((n, b) => n + b.count, 0));
    for (const b of axis.blocks) {
      const aln = set.alignments[b.alignment];
      const mol = set.molecules[set.units.molecule[b.unit]];
      const slot = aln ? aln.rows.slotOf((row ?? aln.referenceRow) as never) : -1;
      for (let k = 0; k < b.count; k++) {
        const i = indexAt(b, k);
        out[b.offset + k] = aln ? (slot < 0 ? 45 : aln.rows.codes[slot * aln.columnCount + i]) : (mol?.code[i] ?? 88);
      }
    }
    // `letter_map` (seqspec:format#layer-types): the letters shown.
    for (let i = 0; map && i < out.length; i++) out[i] = map[out[i] & 127];
    return latin1.decode(out);
  };

  /**
   * Features as a track: canvas from 1,000 features, non-overlapping lanes. Features drawn alike (merged copies' items
   * at the same place) are drawn once, so copies merged into one block do not stack into lanes.
   */
  const track = (l: number, kind: string, features: Feature[], attributes = { layout: "non-overlapping" }) =>
    featureTrack(features, { ...COMMON, height: String(HEIGHT[kind]), "data-svs-layer": String(l), ...attributes });

  /** Unique feature ids: an item's id, else the layer's ref and the row; a repeat takes the base's next suffix. */
  const ids = new Set<string>();
  const suffix = new Map<string, number>();
  const accession = (l: number, a: Annotation, r: number) => {
    const key = a.itemKey[a.item[r]];
    const base = key >= 0 ? name(key) : `${name(lref(l))}:${a.id}:${r}`;
    let id = base;
    if (ids.has(id)) {
      let k = suffix.get(base) ?? 2;
      while (ids.has(`${base}~${k}`)) k++;
      suffix.set(base, k + 1);
      id = `${base}~${k}`;
    }
    ids.add(id);
    return id;
  };
  const tooltip = (a: Annotation, r: number) =>
    ["label", "tooltip", "category", "alt"]
      .map((f) => Annotations.rowValue(a, set, r, f))
      .filter((v) => v !== undefined)
      .join(" · ");
  const shapeOf = (l: number) => {
    const g = themeOf(l, 1);
    const p = (g?.params ?? {}) as Record<string, unknown>;
    const shape = (x: unknown): string => {
      const s = typeof x === "string" ? x : ((x as { shape?: string; fallback?: string })?.shape ?? "circle");
      if (s === "circle" || s === "diamond" || s === "triangle") return s;
      report(GLYPH, lref(l), `no "${s}" shape: drawn as its fallback`);
      return shape((x as { fallback?: string })?.fallback ?? "circle");
    };
    if (g?.theme === "uniform") return () => shape(p.glyph);
    const field = String(p.field ?? (g?.theme === "field_glyph" ? "glyph" : "category"));
    const map = (p.map ?? {}) as Record<string, unknown>;
    return (a: Annotation, r: number) => {
      const v = String(Annotations.rowValue(a, set, r, field));
      return shape(g?.theme === "field_map" ? (map[v] ?? p.default ?? "circle") : v);
    };
  };

  const ranges = (l: number, axis: AxisMeta, kind = "ranges"): ElementDesc => {
    const color = colorer(l);
    // `underline` is a text-layout display: on Nightingale's plain axis, without letters to sit under, a bare line, so
    // drawn as a rectangle; the text layout's loss covers it.
    const shape = kind === "secondary" ? null : params(l).display === "brackets" ? "bridge" : "rectangle";
    const items = new Map<string, Feature>();
    let boundaries = false;
    const labels = kind === "ranges" && params(l).labels === true;
    for (const o of placements(l, axis)) {
      for (const r of o.rows) {
        if (!place(o, o.ann.start[r], o.ann.stop[r])) continue;
        const [first, last] = [runs.start[0], runs.stop[runs.count - 1] - 1];
        boundaries ||= o.ann.start[r] === o.ann.stop[r];
        const cat = String(Annotations.rowValue(o.ann, set, r, "category") ?? "").toLowerCase();
        const sec = SECONDARY[cat.toUpperCase()] ?? cat;
        const sh = shape ?? (sec === "helix" || sec === "strand" ? sec : "rectangle");
        const key = `${o.ann.id}:${o.ann.item[r]}:${o.b.block}`;
        const c = hexOf(color.row(o.ann, r, o.b));
        const f = items.get(key) ?? {
          ...{ accession: accession(l, o.ann, r), start: first + 1, end: last + 1, locations: [{ fragments: [] }] },
          ...{ shape: sh, opacity: Number(params(l).opacity ?? 1) },
          ...(kind === "ranges" ? filled(c) : { color: c }),
          ...{ tooltipContent: tooltip(o.ann, r), svs: { layer: l, annotation: o.ann.id, item: o.ann.item[r] } },
        };
        for (const x of fragmentsOf(o)) f.locations?.[0].fragments.push(x);
        // A range's label: the hover hint's, and laid inside it by `mount` where it fits when the layer asks (`labels`).
        if (kind === "ranges" && !f.svs.label && o.ann.rowLabel && o.ann.rowLabel[r] >= 0) {
          f.svs.label = name(o.ann.rowLabel[r]);
          if (labels) f.svs.inside = true;
        }
        [f.start, f.end] = [Math.min(f.start, first + 1), Math.max(f.end, last + 1)];
        if (o.ann.start[r] === o.ann.stop[r]) f.shape = "bridge";
        // Nightingale's helix is a row of loops: filled, they read as a row of glyphs; outlined, as one coil.
        if (f.shape === "helix") f.fill = "transparent";
        items.set(key, f);
      }
    }
    if (boundaries)
      report(FALLBACK, lref(l), "nightingale-track does not draw boundaries (a bridge over k..k+1 instead)");
    // An arrowhead at a feature's end is the format's, not a hint (NG-09).
    if (params(l).direction_field != null)
      report("svs.render.direction", lref(l), "nightingale-track draws no direction: every feature reads as forward");
    return track(l, kind, [...items.values()]);
  };

  /** A feature's fill and its outline (Nightingale strokes in `color`, fills in `fill`). */
  const filled = (fill: string | undefined) => (fill === undefined ? {} : { fill, color: outline(fill) });
  /** A layer's `opacity` where it is not the default (1). */
  const opacityOf = (l: number) => {
    const x = params(l).opacity;
    return typeof x === "number" && x !== 1 ? x : undefined;
  };
  /** A layer whose element draws no opacity says so when it asks for one (NG-08). */
  const opaque = (l: number, tag: string) => {
    if (opacityOf(l) !== undefined) report(FALLBACK, lref(l), `${tag} has no opacity: drawn opaque`);
  };
  /** A point feature at an axis index. */
  const point = (l: number, o: On, r: number, at: number, color?: string, shape = "circle"): Feature => ({
    ...{ accession: accession(l, o.ann, r), start: at + 1, end: at + 1, ...filled(color), shape },
    ...(opacityOf(l) === undefined ? {} : { opacity: opacityOf(l) }),
    ...{ tooltipContent: tooltip(o.ann, r), svs: { layer: l, annotation: o.ann.id, item: o.ann.item[r] } },
  });
  const points = (l: number, axis: AxisMeta) => {
    const [color, shape, out] = [colorer(l), shapeOf(l), [] as Feature[]];
    for (const o of placements(l, axis)) {
      for (const r of o.rows) {
        for (let i = o.ann.start[r]; i < Math.max(o.ann.stop[r], o.ann.start[r] + 1); i++) {
          // Through a codon step, a point across its codon: one feature over its fragments.
          if (!place(o, i, i + 1)) continue;
          const pt = point(l, o, r, runs.start[0], hexOf(color.row(o.ann, r, o.b)), shape(o.ann, r));
          if (o.path?.heldSpans)
            Object.assign(pt, { end: runs.stop[runs.count - 1], locations: [{ fragments: fragmentsOf(o) }] });
          out.push(pt);
        }
      }
    }
    return track(l, "points", out);
  };

  const values = (l: number, axis: AxisMeta): ElementDesc => {
    const [lo, hi] = [L.domainMin[l], L.domainMax[l]];
    const v: (number | null)[] = new Array(axis.length).fill(null);
    /** The value `x` on every axis position the indices [s, t) reach. */
    const put = (o: On, s: number, t: number, x: number) => {
      if (Number.isNaN(x) || !place(o, s, t)) return;
      const y = Number.isNaN(lo) ? x : Math.min(hi, Math.max(lo, x));
      for (let k = 0; k < runs.count; k++) v.fill(y, runs.start[k], runs.stop[k]);
    };
    // Dense values are one array per target (seqspec:model#annotations); sparse ones rows, the last covering row
    // winning. A streamed annotation has neither: `mount` fetches its bins (stream.ts) into the nulls left here.
    const parts: StreamedPart[] = [];
    for (const o of placements(l, axis)) {
      const d = o.ann.dense?.[o.slot];
      for (let i = 0; d && i < d.length; i++) put(o, i, i + 1, d[i]);
      for (const r of o.rows)
        if (o.ann.start[r] < o.ann.stop[r]) put(o, o.ann.start[r], o.ann.stop[r], o.ann.sparse?.[r] ?? NaN);
      const lazy = o.ann.streamed ? (o.ann.custom.props.get(StreamProp) as StreamedData | undefined) : undefined;
      // Through a path its regions would be asked back through it (ST7): not drawn here.
      if (lazy?.source.data_kind === "values" && o.path) streamNote(l, "streamed values through a path are not drawn");
      else if (lazy?.source.data_kind === "values") {
        parts.push({ stream: lazy, system: [...lazy.systems.values()][o.slot] ?? 0, block: o.b });
        streamNote(l, "streamed values are fetched at one level, the finest whose bins are a pixel wide");
      }
    }
    const ref = themeOf(l, 0);
    if (ref && ref.theme !== "uniform")
      report(THEME, lref(l), `nightingale-linegraph-track has one color: "${themeName(ref)}" is not drawn`);
    const enc = params(l).encoding;
    const color = ref?.theme === "uniform" ? hexOf(Color.decodePrimary(String(ref.params?.color))) : "#4a6fa5";
    const values = v.map((value, k) => ({ position: k + 1, value }));
    const range = Number.isNaN(lo) ? [0, 1] : [lo, hi];
    const series = { name: name(lref(l)), range, color, fill: enc === "area" ? color : undefined, values };
    const domain = Number.isNaN(lo) ? null : ([lo, hi] as const);
    const svs: StreamedValues | null = parts.length ? { ref: name(lref(l)), domain, parts } : null;
    // Its own hover highlight (`highlight-on-click` unset) is its line's point, positions off the pointer's, and
    // replaced the host's: off, so the pointer's hover by position lights it as every row.
    const own = { "highlight-on-click": "" };
    opaque(l, "nightingale-linegraph-track");
    return element("nightingale-linegraph-track", l, "values", own, svs ? { data: [series], svs } : { data: [series] });
  };

  /**
   * A `color` track: a feature per run of one color. Rows color their positions (one palette call per row, the last
   * placement winning); dense values, which have no rows, color every position through the layer's palette, walked in
   * runs of one value (G1-2).
   */
  const color = (l: number, axis: AxisMeta): ElementDesc => {
    const c = colorer(l);
    const [per, owner, from] = [0, 1, 2].map(() => new Int32Array(axis.length).fill(NO));
    const paint = (o: On, color: number, row: number) => {
      for (let k = 0; k < runs.count; k++) {
        per.fill(color, runs.start[k], runs.stop[k]);
        owner.fill(o.ann.id, runs.start[k], runs.stop[k]);
        from.fill(row, runs.start[k], runs.stop[k]);
      }
    };
    for (const o of placements(l, axis)) {
      const d = o.ann.dense?.[o.slot];
      for (let i = 0; d && i < d.length; ) {
        let j = i + 1;
        while (j < d.length && (d[j] === d[i] || (Number.isNaN(d[i]) && Number.isNaN(d[j])))) j++;
        const color = Number.isNaN(d[i]) ? NO : c.value(d[i]);
        if (color >= 0 && place(o, i, j)) paint(o, color, -1);
        i = j;
      }
      for (const r of o.rows) {
        const color = c.row(o.ann, r, o.b);
        if (place(o, o.ann.start[r], Math.max(o.ann.stop[r], o.ann.start[r] + 1))) paint(o, color, r);
      }
    }
    const out: Feature[] = [];
    const op = opacityOf(l);
    for (let k = 0; k < axis.length; k++) {
      if (per[k] < 0) continue;
      let e = k;
      while (e + 1 < axis.length && per[e + 1] === per[k]) e++;
      const item = from[k] < 0 ? -1 : (set.annotations[owner[k]]?.item[from[k]] ?? -1);
      const [accession, svs] = [`${name(lref(l))}:${k}`, { layer: l, annotation: owner[k], item }];
      const h = hexOf(per[k]);
      out.push({
        accession,
        start: k + 1,
        end: e + 1,
        color: h,
        fill: h,
        svs,
        ...(op === undefined ? {} : { opacity: op }),
      });
      k = e;
    }
    return track(l, "color", out);
  };

  /**
   * `letters` with `of`: through index steps that polymer's letters where they land, `-` elsewhere; through a codon
   * step a feature per residue over its codon (fragments where an intron splits it), its one-letter code laid inside it
   * by `mount` where it fits.
   */
  const lettersOfOther = (l: number, axis: AxisMeta): ElementDesc[] => {
    const seq = new Uint8Array(axis.length).fill(45);
    const map = mapOf(l);
    const features: Feature[] = [];
    const ref = themeOf(l, 0);
    const dict = ref?.theme === "letter" ? (letterColors((ref.params?.colors ?? null) as never, "protein") ?? {}) : {};
    for (let p = L.placementStart[l]; p < L.placementStart[l] + L.placementCount[l]; p++) {
      const at = placementOf(p, axis);
      const path = at?.path ?? null;
      const b = at?.b;
      const from = path?.from ?? { kind: CoordinateSystemKind.unit, id: b?.unit ?? -1 };
      const mol = set.molecules[from.kind === CoordinateSystemKind.unit ? set.units.molecule[from.id] : from.id];
      if (!at || !b || !mol) continue;
      for (let k = 0; !path?.heldSpans && k < b.count; k++) {
        const j = path ? sourceAt(path, at.read, indexAt(b, k)) : indexAt(b, k);
        if (j >= 0) seq[b.offset + k] = map ? map[mol.code[j] & 127] : mol.code[j];
      }
      for (let r = 0; path?.heldSpans && r < mol.length; r++) {
        if (!place(at, r, r + 1)) continue;
        // The letter shown (`letter_map`): its color, and one letter where the map changes it (as the viewer draws).
        const stored = mol.code[r];
        const code = map ? map[stored & 127] : stored;
        const here = {
          start: runs.start[0] + 1,
          end: runs.stop[runs.count - 1],
          locations: [{ fragments: fragmentsOf(at) }],
        };
        const letter = String.fromCharCode(code);
        const tooltipContent = `${code === stored ? threeLetter(code) : letter} ${r + 1}`;
        const look = { ...filled(hexOf(dict[letter])), tooltipContent };
        // Its letter drawn inside it (mount's label layer); where it does not fit, the hover's hint names the residue.
        const svs = { layer: l, annotation: -1, item: r, label: letter, inside: true, hint: tooltipContent };
        // Several polymers on a block (`of` a list, `mapped`): each residue keyed by its polymer too.
        const by = L.placementCount[l] > 1 ? `${from.kind}.${from.id}:` : "";
        features.push({ accession: `${name(lref(l))}:${b.block}:${by}${r}`, ...here, ...look, svs });
      }
    }
    if (!features.length) return [element("nightingale-sequence", l, "letters", {}, { data: latin1.decode(seq) })];
    report(
      FALLBACK,
      lref(l),
      "no three-letter residues: a residue over its codon is a feature, its one-letter code inside it where it fits",
    );
    return [track(l, "ranges", features)];
  };

  const letters = (l: number, axis: AxisMeta): ElementDesc[] => {
    opaque(l, "nightingale-sequence");
    if (params(l).of != null) return lettersOfOther(l, axis);
    const seq = lettersOf(axis, (params(l).row as number | null) ?? null, mapOf(l));
    const out = [element("nightingale-sequence", l, "letters", {}, { data: seq })];
    const ref = themeOf(l, 0);
    const b0 = axis.blocks[0];
    const alphabet =
      b0 && (set.alignments[b0.alignment]?.alphabet ?? set.molecules[set.units.molecule[b0.unit]]?.alphabet);
    if (!ref || !b0) return out;
    if (ref.theme !== "letter") {
      report(THEME, lref(l), `no slot for the "${themeName(ref)}" theme on letters: drawn uncolored`);
      return out;
    }
    const dict = letterColors((ref.params?.colors ?? null) as never, alphabet ?? null) ?? {};
    const colors = [...new Set(Object.values(dict))];
    const scale = Object.entries(dict)
      .filter(([k]) => k.length === 1)
      .map(([k, c]) => `${k}:${colors.indexOf(c)}`);
    // nightingale-colored-sequence colors a letter not in its scale as value 0, the first color (NG-06): each such letter
    // (`X`, `-`, `*`, …) takes the ground, as the viewer leaves it uncolored.
    const uncolored = lettersOutside(seq, dict);
    if (uncolored.length) {
      colors.push(GROUND);
      for (const k of uncolored) scale.push(`${k}:${colors.length - 1}`);
    }
    const range = (colors.length > 1 ? colors : [colors[0] ?? 0, colors[0] ?? 0]).map((c, k) => `${hexOf(c)}:${k}`);
    const attributes = { scale: scale.join(","), "color-range": range.join(",") };
    // The letters over their colored bar in one row of one height at any zoom (zoomed out, only the bar shows), not
    // stacked under it: an empty band above the bar while the letters do not fit.
    const bar = element("nightingale-colored-sequence", l, "colors", attributes, { data: seq });
    return [el("div", { class: STACK }, {}, [bar, ...out])];
  };

  /** Tracks labeled by what they draw rather than their own label (a lazy alignment's rows as its reference row). */
  const relabeled = new Map<number, { label: string; hint: string }>();
  const clustal = (r: ThemeRef | null) => r?.table === "clustal_x";
  const rows = (l: number, axis: AxisMeta, t: number): ElementDesc[] => {
    const b = axis.blocks[0];
    const aln = b && set.alignments[b.alignment];
    if (!aln || aln.lazy) {
      report(FALLBACK, lref(l), "nightingale-msa needs every row: a lazy alignment is its reference row's letters");
      // The row is labeled by what it draws: the reference row (its label once known), the rest said in the hint.
      if (aln) {
        const ref = aln.labelOf(aln.referenceRow as never);
        const others = (aln.rowCount - 1).toLocaleString("en-US");
        const label = `${ref >= 0 ? name(ref) : "Reference row"} (reference row; ${others} rows not drawn)`;
        const track = name(T.label[t]);
        const hint = `Track ${track ? `"${track}"` : "of rows"}: the alignment's ${aln.rowCount.toLocaleString("en-US")} rows load as they are viewed, and nightingale-msa needs every row, so only the reference row's letters are drawn.`;
        relabeled.set(t, { label, hint });
      }
      return letters(l, axis);
    }
    const p = params(l);
    const lanes = Array.from({ length: T.laneCount[t] }, (_, k) => T.laneStart[t] + k);
    const data = lanes.map((k) => ({
      name: name(aln.labelOf(view.lanes.row[k] as never)),
      sequence: lettersOf(axis, view.lanes.row[k], mapOf(l)),
    }));
    const ref = themeOf(l, 0);
    // A `rules` theme (a named table but Clustal X's, as `identity_to_row`) draws its first rule's nested theme.
    const base = (r: ThemeRef | null): ThemeRef | null => {
      const inner =
        r?.theme === "rules" && !clustal(r)
          ? ((r.params?.rules ?? []) as { theme?: ThemeRef }[]).find((x) => x.theme)
          : null;
      return inner?.theme ? base(inner.theme) : r;
    };
    const scheme = (r: ThemeRef | null) => {
      if (clustal(r)) return "clustal2";
      const c = r?.theme === "letter" ? (r.params?.colors ?? null) : null;
      if (typeof c === "string" && SCHEMES[c]) return SCHEMES[c];
      return aln.alphabet === "protein" ? "zappo" : "nucleotide";
    };
    const shown = base(ref ?? { theme: "letter" });
    const lost = [
      p.show_tree !== false && aln.tree && "the tree",
      Array.isArray(p.data_columns) && p.data_columns.length > 0 && "data columns",
      lanes.some((k) => view.lanes.depth[k] > 0) && "variant indentation",
      ref && ref !== shown && `the "${themeName(ref)}" theme (its base is drawn)`,
      clustal(shown) && "Clustal X (its static clustal2 scheme instead)",
      shown?.theme === "letter" &&
        shown.params?.colors != null &&
        !SCHEMES[String(shown.params.colors)] &&
        "these letter colors",
    ].filter(Boolean);
    if (lost.length)
      report(ref && ref !== shown ? THEME : FALLBACK, lref(l), `nightingale-msa does not draw ${lost.join(", ")}`);
    // A theme with no scheme of nightingale-msa's (`percent_identity`, `uniform`, `by_unit`, …) draws a default (NG-03).
    if (shown && shown.theme !== "letter" && !clustal(shown))
      report(THEME, lref(l), `nightingale-msa has no "${themeName(shown)}" scheme: drawn in ${scheme(shown)}`);
    const height = String(Math.min(lanes.length, 16) * HEIGHT.tile);
    const attrs = { "color-scheme": scheme(shown), "tile-height": String(HEIGHT.tile), height };
    opaque(l, "nightingale-msa");
    return [element("nightingale-msa", l, "tile", attrs, { data })];
  };

  const variants = (l: number, axis: AxisMeta): ElementDesc => {
    const [c, out] = [colorer(l), [] as Record<string, unknown>[]];
    let lost = false;
    for (const o of placements(l, axis)) {
      for (const r of o.rows) {
        const alt = String(Annotations.rowValue(o.ann, set, r, "alt") ?? "");
        const at = firstAt(o, o.ann.start[r]);
        if (o.ann.stop[r] - o.ann.start[r] !== 1 || alt.length !== 1 || at < 0) {
          lost ||= at >= 0;
          continue;
        }
        const { accession: id, tooltipContent, svs } = point(l, o, r, at);
        const nightingale = { xrefNames: [], hasPredictions: false, consequenceType: "" };
        out.push({
          accession: id,
          variant: alt,
          start: at + 1,
          color: hexOf(c.row(o.ann, r, o.b)),
          tooltipContent,
          svs,
          ...nightingale,
        });
      }
    }
    if (lost) report(FALLBACK, lref(l), "nightingale-variation-canvas draws no multi-residue changes or indels");
    opaque(l, "nightingale-variation-canvas");
    const data = { sequence: lettersOf(axis, null), variants: out };
    // `split.empty: "hide"`: the grid's rows with a variant only, nightingale-variation's own.
    const hide = (params(l).split as { empty?: string } | null)?.empty === "hide";
    return element("nightingale-variation-canvas", l, "variants", hide ? { "condensed-view": "" } : {}, { data });
  };

  const pairs = (l: number, axis: AxisMeta, all: AxisMeta[]): ElementDesc[] => {
    // Each pair's ends, by annotation and row (NG-19: no packed key): where each lies, on which side and block.
    const ends = new Map<number, Map<number, PairEnd[]>>();
    all.forEach((m, k) => {
      for (const o of placements(l, m)) {
        const rows = ends.get(o.ann.id) ?? ends.set(o.ann.id, new Map()).get(o.ann.id);
        for (const r of o.rows) {
          const [s, t] =
            o.side === 1 && o.ann.pairs
              ? [o.ann.pairs.bStart[r], o.ann.pairs.bStop[r]]
              : [o.ann.start[r], o.ann.stop[r]];
          const at = s < 0 ? o.b.offset : place(o, s, t) ? runs.start[0] : -1;
          if (at >= 0 && rows)
            (rows.get(r) ?? rows.set(r, []).get(r))?.push({ manager: k, at, side: o.side, o, row: r });
        }
      }
    });
    const me = all.indexOf(axis);
    const contacts: number[][] = [];
    const apart: [On, number, number][] = [];
    // A side-0 end joins the side-1 end on its own block first (each copy's pair on its copy, NG-02). Only the ends
    // left then join across blocks, and only where a pair can lie across them: two blocks of one unit (or alignment),
    // or a pair whose ends are on different targets (LA-05: a copy whose partner its range cut took the next copy's,
    // before that copy's own end). A pair is between managers (points) only when its two ends are.
    const B = view.blocks;
    const sameThing = (x: number, y: number) =>
      (B.unit[x] >= 0 && B.unit[x] === B.unit[y]) || (B.alignment[x] >= 0 && B.alignment[x] === B.alignment[y]);
    const systemOf = (x: PairEnd) => endsOf(x.o.ann).find((t) => t.slot === x.o.slot)?.system;
    const across = (a: PairEnd, b: PairEnd) => {
      if (sameThing(a.o.b.block, b.o.b.block)) return true;
      const [s, t] = [systemOf(a), systemOf(b)];
      return s !== undefined && t !== undefined && !Coords.sameSystem(s, t);
    };
    for (const rows of ends.values())
      for (const e of rows.values()) {
        const taken = new Set<PairEnd>();
        const join = (a: PairEnd, b: PairEnd) => {
          taken.add(b);
          if (a.manager === me && b.manager === me) contacts.push([a.at + 1, b.at + 1, 1]);
          else for (const x of [a, b]) if (x.manager === me) apart.push([x.o, x.row, x.at]);
        };
        const left: PairEnd[] = [];
        for (const a of e) {
          if (a.side === 1) continue;
          const b = e.find((x) => x.side === 1 && !taken.has(x) && x.o.b === a.o.b);
          if (b) join(a, b);
          else left.push(a);
        }
        for (const a of left) {
          const b = e.find((x) => x.side === 1 && !taken.has(x) && across(a, x));
          if (b) join(a, b);
        }
      }
    const p = params(l);
    const a = placements(l, axis)[0]?.ann;
    const lost = [
      a?.category && "categories",
      p.color != null && themeOf(l, 0) && themeOf(l, 0)?.theme !== "uniform" && "per-pair colors",
      p.side === "below" && "side",
    ].filter(Boolean);
    if (lost.length) report(FALLBACK, lref(l), `nightingale-links does not draw ${lost.join(", ")}`);
    // A lost direction is the format's own code (LA-30), as a track's.
    if (a?.pairs?.direction)
      report("svs.render.direction", lref(l), "nightingale-links draws no direction: every pair reads as undirected");
    opaque(l, "nightingale-links");
    // nightingale-links colors its points by its own ramp; a uniform color goes in by `mount`'s style (PAIR_COLOR).
    const ref = themeOf(l, 0);
    const one = ref?.theme === "uniform" ? hexOf(Color.decodePrimary(String(ref.params?.color))) : undefined;
    const style: Record<string, string> = one ? { style: `${PAIR_COLOR}: ${one}` } : {};
    const out = [element("nightingale-links", l, "pairs", { "min-distance": "0", ...style }, { contacts })];
    if (apart.length === 0) return out;
    report(FALLBACK, lref(l), "nightingale-links joins no stacked managers: pairs between blocks are points");
    return [
      ...out,
      track(
        l,
        "points",
        apart.map(([o, r, at]) => point(l, o, r, at)),
      ),
    ];
  };

  /** Matrix rows (a position each) on the axis: the states and, per axis position, its cells. */
  const cells = (l: number, axis: AxisMeta) => {
    const out = new Map<number, Float64Array>();
    let states: string[] = [];
    for (const o of placements(l, axis)) {
      const m = o.ann.matrix;
      if (!m) continue;
      states = Array.from(m.states, (s) => name(s));
      for (const r of o.rows) {
        const at = firstAt(o, o.ann.start[r]);
        if (at >= 0) out.set(at, m.values.subarray(r * m.states.length, (r + 1) * m.states.length));
      }
    }
    return { states, out };
  };

  const logo = (l: number, axis: AxisMeta): ElementDesc => {
    const { states, out } = cells(l, axis);
    const index = [...out.keys()].sort((x, y) => x - y);
    const p = (at: number, k: number) => out.get(at)?.[k] || 0; // NaN: 0
    const probabilities = Object.fromEntries(states.map((s, k) => [s, index.map((at) => p(at, k))]));
    const data = { index: index.map((x) => x + 1), probabilities };
    opaque(l, "nightingale-conservation-track");
    return element("nightingale-conservation-track", l, "logo", {}, { data });
  };

  const matrix = (l: number, axis: AxisMeta, id: string): ElementDesc => {
    const { states, out } = cells(l, axis);
    const seenStates = new Map<string, number>();
    const unique = states.map((s) => {
      const k = (seenStates.get(s) ?? 0) + 1;
      seenStates.set(s, k);
      return k > 1 ? `${s} (${k})` : s;
    });
    // The states drawn, in their keyed order (seqspec:format#keyed-sets): a state left out has no row.
    const members = L.members[l];
    const shown = members ? new Set(members) : null;
    const data: { xValue: number; yValue: string; score: number }[] = [];
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    // The layer's domain clamps the scores (NG-07), as the viewer draws them.
    const [min, max] = [L.domainMin[l], L.domainMax[l]];
    const clamp = (x: number) => (Number.isNaN(min) ? x : Math.min(max, Math.max(min, x)));
    for (const [at, v] of [...out].sort((x, y) => x[0] - y[0]))
      v.forEach((raw, k) => {
        if (Number.isNaN(raw) || (shown && !shown.has(k))) return;
        const score = clamp(raw);
        data.push({ xValue: at + 1, yValue: unique[k], score });
        lo = Math.min(lo, score);
        hi = Math.max(hi, score);
      });
    // The heatmap's color scale spans the cells' scores: anchors under the first cell, which overwrites them in the
    // grid, spread it to the domain's ends where the cells do not reach them (NG-07); a single score throws ("not
    // numeric and monotonous"), drawing nothing (a fingerprint's counts all 1), so it gets one at 0 (1 for 0), the
    // domain from zero as the viewer's counts ([0, hi]).
    const anchors = data.length === 0 ? [] : [min < lo && min, max > hi && max].filter((x) => x !== false);
    if (data.length > 0 && anchors.length === 0 && lo === hi) anchors.push(lo === 0 ? 1 : 0);
    data.unshift(...anchors.map((score) => ({ ...data[0], score: score as number })));
    if (params(l).color != null)
      report(THEME, lref(l), "nightingale-sequence-heatmap has its own colors: the color theme is not drawn");
    const x = Array.from({ length: axis.length }, (_, k) => k + 1);
    const rows = members ? Array.from(members, (k) => unique[k]) : unique;
    const height = String(Math.min(160, Math.max(40, rows.length * 8)));
    opaque(l, "nightingale-sequence-heatmap");
    const e = element("nightingale-sequence-heatmap", l, "matrix", { "heatmap-id": id, height });
    e.calls.push({ method: "setHeatmapData", args: [x, rows, data] });
    return e;
  };

  const ruler = (l: number, axis: AxisMeta): ElementDesc => {
    const s = axis.section;
    const own = params(l).numbering;
    const differs =
      S.blockCount[s] > 1 ||
      axis.blocks.some((b) => indexAt(b, 0) !== 0) ||
      (own != null && own !== "ordinal") ||
      axis.blocks.some(
        (b) => view.blocks.numbering[b.block] && name(view.blocks.numbering[b.block]?.name ?? -1) !== "ordinal",
      );
    // A ruler `of` another system loses that system's positions (LA-31); else the section's labels.
    if (params(l).of != null)
      report(
        FALLBACK,
        lref(l),
        "nightingale-navigation: the ruler shows axis positions, not those of the system it is of",
      );
    else if (differs)
      report(FALLBACK, lref(l), "nightingale-navigation: the ruler shows axis positions, not the section's labels");
    opaque(l, "nightingale-navigation");
    return element("nightingale-navigation", l, "ruler");
  };

  /** One layer's elements on one manager, by the type drawn. */
  type Draw = (l: number, axis: AxisMeta, all: AxisMeta[], t: number, id: string) => ElementDesc | ElementDesc[];
  const draw: Readonly<Record<string, Draw>> = {
    ...{ letters, points, variants, values, color, ruler, logo, pairs, ranges: (l, a) => ranges(l, a) },
    ...{ rows: (l, a, _, t) => rows(l, a, t), secondary: (l, a) => ranges(l, a, "secondary") },
    matrix: (l, a, _, __, id) => matrix(l, a, id),
  };

  const hidden = (t: number) => {
    if (T.hidden[t]) return true;
    for (let g = T.group[t]; g >= 0; g = G.parent[g]) if (G.hidden[g]) return true;
    return false;
  };
  /** A track's hint: its data's provenance (the format has no track description). */
  const hintOf = (ls: number[]) => {
    const said = new Set<string>();
    for (const l of ls) {
      for (let k = L.dataStart[l]; k < L.dataStart[l] + L.dataCount[l]; k++) {
        const p = set.annotations[view.data[k]]?.provenance;
        const text =
          p &&
          [p.source, p.method, p.note]
            .filter((x) => x >= 0)
            .map(name)
            .join(" · ");
        if (text) said.add(text);
      }
    }
    return [...said].join("; ");
  };
  /** A label cell: its text and an optional hint (`data-svs-hint`, for the host's tooltip), or a group's toggle. */
  const text = (label: string) => el("span", { class: "svs-label-text" }, { textContent: label });
  const hintMark = (hint: string) =>
    el("span", { class: "svs-hint", tabindex: "0", "aria-label": hint, "data-svs-hint": hint }, { textContent: "i" });
  const row = (label: ElementDesc[], children: ElementDesc[], depth: number, cls = "svs-track") =>
    el("div", { class: cls, style: `--svs-depth: ${depth}` }, {}, [
      el("div", { class: "svs-label" }, {}, label),
      el("div", { class: "svs-layers" }, {}, children),
    ]);
  /**
   * A group: a header row (a disclosure toggle and its label; collapsed, its tracks' features merged into one row
   * beside it) over a body of its tracks and groups, indented. `mount` toggles `data-collapsed`.
   */
  interface Open {
    g: number;
    body: ElementDesc[];
    header: ElementDesc;
    features: Feature[];
  }
  const group = (label: ElementDesc[], ref: number, collapsed: boolean, depth: number, into: ElementDesc[]) => {
    const toggle = el("button", { type: "button", class: "svs-toggle", "aria-expanded": String(!collapsed) }, {}, [
      el("span", { class: "svs-disclosure", "aria-hidden": "true" }, { textContent: collapsed ? "▸" : "▾" }),
      ...label,
    ]);
    const header = row([toggle], [], depth, "svs-track svs-group-header");
    const body = el("div", { class: "svs-group-body" });
    const attributes = {
      class: "svs-group",
      "data-svs-group": name(ref),
      ...(collapsed ? { "data-collapsed": "" } : {}),
    };
    into.push(el("div", attributes, {}, [header, body]));
    return { body: body.children, header };
  };
  const openGroup = (g: number, depth: number, into: ElementDesc[]): Open => ({
    g,
    ...group([text(name(G.label[g]))], G.ref[g], G.collapsed[g] === 1, depth, into),
    features: [],
  });
  /**
   * A split track: a group's rows, its header the track's label with the layer drawn whole beside it (collapsed), its
   * body a track per value of the layer's features (`splitRows`, as every viewer reads them).
   */
  const split = (
    t: number,
    l: number,
    label: ElementDesc[],
    els: ElementDesc[],
    depth: number,
    into: ElementDesc[],
  ) => {
    const s = group(label, T.ref[t], T.expanded[t] === 0, depth, into);
    s.header.children[1].children.push(...els);
    const range = (xs: readonly unknown[]) => xs.slice(L.dataStart[l], L.dataStart[l] + L.dataCount[l]);
    const anns = range(Array.from(view.data)).map((id) => set.annotations[id as number] as Annotation);
    const rows = splitRows(set, anns, range(view.kept) as (Uint8Array | null)[], params(l).split as never);
    /** Per annotation, an item's value (its first row's): rows from the last to the first, so the first wins. */
    const itemValue = new Map<number, Int32Array>();
    for (const [d, a] of anns.entries()) {
      const items = new Int32Array(a.itemCount).fill(NO);
      for (let r = a.count - 1; r >= 0; r--) items[a.item[r]] = rows.rowOf[d][r];
      itemValue.set(a.id, items);
    }
    const mine = els.filter((x) => x.attributes["data-svs-layer"] === String(l));
    const at = { ...COMMON, height: String(HEIGHT[typeOf(l)]), layout: "non-overlapping", "data-svs-layer": String(l) };
    /** An element's features bucketed into the split values, in one pass. */
    const buckets = mine.map((x) => {
      const out: Feature[][] = rows.values.map(() => []);
      for (const f of x.properties.data as Feature[]) {
        const v = itemValue.get(f.svs.annotation)?.[f.svs.item] ?? NO;
        if (v >= 0) out[v].push(f);
      }
      return out;
    });
    rows.values.forEach((value, v) => {
      s.body.push(
        row(
          [text(value)],
          buckets.map((b) => featureTrack(structuredClone(b[v]), at)),
          depth + 1,
        ),
      );
    });
  };
  /** The features an element draws, for its groups' merged rows: a track's, a variant as a point. */
  const featuresOf = (x: ElementDesc): Feature[] => {
    if (x.tag === "nightingale-track" || x.tag === "nightingale-track-canvas") return x.properties.data as Feature[];
    if (x.tag !== "nightingale-variation-canvas") return [];
    const { variants } = x.properties.data as { variants: (Feature & { variant: string })[] };
    return variants.map(({ accession, start, color, tooltipContent, svs }) => {
      return { accession, start, end: start, color, shape: "circle", tooltipContent, svs };
    });
  };

  const sections: ElementDesc[] = [];
  let shownSections = 0;
  for (let s = 0; s < S.count; s++) if (!S.hidden[s]) shownSections++;
  for (let s = 0; s < S.count; s++) {
    if (S.hidden[s]) continue;
    const sref = S.ref[s];
    const length = sectionLength(view, s);
    if (length > NIGHTINGALE_MAX_POSITIONS) {
      const [n, max] = [length.toLocaleString("en-US"), NIGHTINGALE_MAX_POSITIONS.toLocaleString("en-US")];
      report(
        "svs.render.omitted",
        sref,
        `a section of ${n} positions, past its limit of ${max} (an element per position)`,
      );
      continue;
    }
    const all = axesOf(view, set, s);
    const n = S.blockCount[s];
    if (S.arrangement[s] === Arrangement.side && n > 1)
      report(FALLBACK, sref, "one axis per manager: no gaps between blocks; block labels are a first track");
    if (S.arrangement[s] === Arrangement.stacked && n > 1)
      report(FALLBACK, sref, "one manager per stacked block: zoom is not linked");
    for (let b = S.blockStart[s]; b < S.blockStart[s] + n; b++)
      if (view.blocks.varying[b]) report(FALLBACK, sref, "merged copies: positions where they differ are not marked");
    // The tracks of more than one layer, said once for the section with the first few (GV-07's rule).
    const composited = new Set<number>();
    const managers = all.map((axis, m) => {
      const children: ElementDesc[] = [];
      if (S.arrangement[s] === Arrangement.side && n > 1) {
        // Block names are laid over their bars by `mount` (Nightingale's tracks draw no text).
        const data: Feature[] = axis.blocks.map((b) => {
          const label = name(view.blocks.label[b.block]);
          return {
            ...{ accession: `block-${b.block}`, start: b.offset + 1, end: b.offset + b.count, fill: "#e8e8e8" },
            ...{ tooltipContent: label, svs: { layer: -1, annotation: -1, item: b.block, label, inside: true } },
          };
        });
        const attributes = { ...COMMON, height: String(HEIGHT.blocks), [LABELED]: "" };
        children.push(row([text("")], [el("nightingale-track", attributes, { data })], 0));
      }
      // Groups open as their first drawn track comes (tracks come in document order, a group's together), so a group
      // that draws nothing has no header.
      const [open, opened]: Open[][] = [[], []];
      for (let t = 0; t < T.count; t++) {
        if (T.section[t] !== s || hidden(t)) continue;
        const ls: number[] = [];
        // A layer placed nowhere (its data unavailable or reaching no block; the resolution reports it) draws nothing.
        for (let l = 0; l < L.count; l++) if (L.track[l] === t && L.type[l] >= 0 && L.placementCount[l] > 0) ls.push(l);
        if (ls.length > 1) composited.add(t);
        const els = ls.flatMap((l) => draw[typeOf(l)]?.(l, axis, all, t, `svs-heatmap-${s}-${m}-${l}`) ?? []);
        if (!els.length) continue;
        const chain: number[] = [];
        for (let g = T.group[t]; g >= 0; g = G.parent[g]) chain.unshift(g);
        let k = 0;
        while (k < open.length && k < chain.length && open[k].g === chain[k]) k++;
        open.length = k;
        for (; k < chain.length; k++) {
          const o = openGroup(chain[k], k, k ? open[k - 1].body : children);
          open.push(o);
          opened.push(o);
        }
        const own = relabeled.get(t);
        const hint = [own?.hint, hintOf(ls)].filter(Boolean).join(" ");
        const label = [text(own?.label ?? name(T.label[t])), ...(hint ? [hintMark(hint)] : [])];
        const into = open.length ? open[open.length - 1].body : children;
        // A split layer's rows; variants keep nightingale-variation's own grid, a row per residue.
        const sl = T.split[t];
        if (sl >= 0 && ls.includes(sl) && typeOf(sl) !== "variants") split(t, sl, label, els, open.length, into);
        else into.push(row(label, els, open.length));
        if (sl >= 0 && typeOf(sl) === "variants" && (params(sl).split as { by?: string } | null)?.by !== "alt")
          report(
            FALLBACK,
            T.ref[t],
            "nightingale-variation's rows are residues: a split by another field is not drawn",
          );
        for (const o of open) for (const x of els) for (const f of featuresOf(x)) o.features.push(f);
      }
      // A group's merged row: its tracks' features (copies: Nightingale keeps what it is given), drawn while collapsed.
      for (const o of opened) {
        if (!o.features.length) continue;
        const attributes = { ...COMMON, height: String(HEIGHT.ranges), layout: "non-overlapping" };
        o.header.children[1].children.push(featureTrack(structuredClone(o.features), attributes));
      }
      const focus = S.focus[s];
      const fb = focus && axis.blocks.find((b) => b.block === focus.block);
      sink.block = fb as AxisBlock;
      runs.count = 0;
      if (fb && focus) sink.add(focus.start, focus.stop);
      // The section's `highlights`, which `mount` keeps under a host's and under Nightingale's own hover.
      const lit = sectionRuns(view, set, axis);
      const attributes: Record<string, string> = { "reflected-attributes": "display-start,display-end,highlight" };
      // The focus's range and the highlights open on each element, never on the manager: a manager re-applies its own
      // `display-*` and `highlight` after every child's change (NG-01), undoing each zoom and hover. `length` too: a
      // manager re-sets it at each change, and nightingale-track redraws whole on any `length` set, the same value too.
      // nightingale-sequence draws an axis of its own while its letters do not fit: a scale the document did not ask
      // for (a `ruler` layer is how it asks), so always hidden.
      const start: Record<string, string> = {
        length: String(axis.length),
        ...(runs.count > 0
          ? { "display-start": String(runs.start[0] + 1), "display-end": String(runs.stop[runs.count - 1]) }
          : {}),
        ...(lit.count > 0 ? { highlight: AxisRuns.text(lit) } : {}),
      };
      for (const x of nightingaleElements(children)) {
        Object.assign(x.attributes, start);
        if (x.tag === "nightingale-sequence") x.attributes.class = AXISLESS;
      }
      return el(
        "nightingale-manager",
        attributes,
        { svs: lit.count > 0 ? { ...axis, highlights: lit } : axis },
        children,
      );
    });
    if (composited.size > 0) {
      const refs = [...composited].map((t) => `"${name(T.ref[t])}"`);
      const more = refs.length > 3 ? ` and ${refs.length - 3} more` : "";
      const which = `${refs.length} track(s): ${refs.slice(0, 3).join(", ")}${more}`;
      report(FALLBACK, sref, `no layer compositing: each layer is its own element, in order (${which})`);
    }
    // One gutter on every element of the section, so its tracks (and stacked managers, and mount's labels) keep one
    // axis; wider on the left where a variation track's residue axis is drawn in it.
    const elements = managers.flatMap((m) => nightingaleElements(m.children));
    const variants = elements.some((x) => x.tag === "nightingale-variation-canvas");
    for (const x of elements) {
      x.attributes["margin-left"] = String(variants ? GUTTER.variants : GUTTER.side);
      x.attributes["margin-right"] = String(GUTTER.side);
    }
    // The section's label, else its ref; shown when more than one section is.
    const given = S.params[s].label;
    const label =
      shownSections > 1
        ? [el("div", { class: "svs-section-label" }, { textContent: typeof given === "string" ? given : name(sref) })]
        : [];
    sections.push(el("div", { class: "svs-section", "data-svs-section": name(sref) }, {}, [...label, ...managers]));
  }
  return { root: el("div", { class: "svs-nightingale", "data-svs-view": name(view.ref) }, {}, sections), issues };
}
