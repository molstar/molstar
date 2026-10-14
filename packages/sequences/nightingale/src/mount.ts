/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `mount` (seqspec:libraries#nightingale): the one module that imports Nightingale. It builds the described elements,
// sets their data once they have rendered (a streamed linegraph's after fetching its bins), maps Nightingale's `change`
// events back to block, index, layer, item and the item's extent, copies display ranges within a scroll group, lays
// labels over labeled tracks (block names, range labels; Nightingale's tracks draw no text), and lights a host's loci
// through each manager's `highlight`.

import "@nightingale-elements/nightingale-colored-sequence";
import "@nightingale-elements/nightingale-conservation-track";
import "@nightingale-elements/nightingale-linegraph-track";
import "@nightingale-elements/nightingale-links";
import "@nightingale-elements/nightingale-manager";
import "@nightingale-elements/nightingale-msa";
import "@nightingale-elements/nightingale-navigation";
import "@nightingale-elements/nightingale-sequence";
import "@nightingale-elements/nightingale-sequence-heatmap";
import "@nightingale-elements/nightingale-track";
import "@nightingale-elements/nightingale-track-canvas";
import "@nightingale-elements/nightingale-variation-canvas";
import { type Closure, projectionOf, type ResolvedView } from "@molstar/seq-core/format";
import { Extents, type Issue, type Loci, type SeqSet } from "@molstar/seq-core/model";
import {
  type AxisMeta,
  AxisRuns,
  AxisSink,
  type BlockRun,
  BlockRuns,
  extentOf,
  highlightRuns,
  keepLit,
  type Located,
  locate,
  prepareLoci,
} from "./axis.ts";
import type { ElementDesc } from "./element.ts";
import {
  arbitrate,
  type Detail,
  HoverResult,
  type HoverTarget,
  hoverAt,
  type Light,
  mergeSeries,
  newLight,
  nextHighlight,
  noteLit,
  positionOf,
  rangesOf,
} from "./events.ts";
import {
  axisView,
  INK,
  inkOn,
  type LabelBox,
  type Laid,
  type LayoutLike,
  layLabels,
  type NightingaleAxis,
  watchDimensions,
} from "./labels.ts";
import { AXISLESS, LABELED, PAIR_COLOR, STACK } from "./materialize.ts";
import { fetchStreamed, type StreamedValues } from "./stream.ts";

/** A Nightingale event in the view's terms; -1: not known (a position only, or a whole element). */
export interface NightingaleEvent {
  readonly kind: "hover" | "click" | "leave";
  readonly section: number;
  readonly block: number;
  readonly index: number;
  readonly layer: number;
  readonly annotation: number;
  readonly item: number;
  /**
   * What the event covers, as runs of block indices: a hovered item's every fragment (a range, a secondary element,
   * a color run, a contact's ends), a residue's one position; empty on leave. `extentLoci` makes it the host's loci.
   */
  readonly extent: readonly BlockRun[];
  /**
   * On hover, a label not shown in full (cut short, hidden where not even its first characters fit, or not laid) and
   * the drawn element it belongs to, for the host's tooltip; null otherwise.
   */
  readonly hint: { readonly anchor: Element; readonly text: string } | null;
}
export interface Mounted {
  dispose(): void;
  /**
   * Lights the loci in every manager (null: clears) through the shared projection of the mount's `source`, `element` as
   * the closure unless `opts` says `containing` (a pointer hover); selection uses the same channel. Loci of another
   * set light nothing.
   */
  highlight(loci: Loci | null, opts?: { closure?: Closure }): void;
  /** Settles once every element has its data, streamed bins fetched; the issues of tiles that failed. */
  readonly ready: Promise<Issue[]>;
}

/** A linegraph track, whose first series names the one key of its events' features. */
interface LineGraph {
  data?: { name?: string }[];
}
const KIND: Readonly<Record<string, NightingaleEvent["kind"]>> = {
  mouseover: "hover",
  click: "click",
  mouseout: "leave",
  onResidueMouseEnter: "hover",
  onResidueClick: "click",
};

/** A feature as `materialize` describes it and Nightingale draws it (the same object: its data is ours). */
interface LabelFeature {
  start: number;
  end: number;
  locations?: { fragments: { start: number; end: number }[] }[];
  color?: string;
  fill?: string;
  shape?: string;
  opacity?: number;
  svs: { layer: number; label?: string; inside?: boolean; hint?: string };
}
/** The state of a feature track `mount` reads to place labels: what its `xScale` and layout are made of. */
interface TrackState {
  "display-start"?: number;
  "display-end"?: number;
  length?: number;
  width?: number;
  "margin-left"?: number;
  "margin-right"?: number;
  layoutObj?: LayoutLike;
  onDimensionsChange?(): void;
}
/** The view a track shows, from its own state (its xScale's inputs), into `out`. */
const trackView = (out: NightingaleAxis, t: TrackState): NightingaleAxis =>
  axisView(
    out,
    t["display-start"],
    t["display-end"],
    t.length ?? 1,
    t.width ?? 0,
    t["margin-left"] ?? 0,
    t["margin-right"] ?? 0,
  );

/** One label laid over a feature, with its element. */
interface FeatureLabel extends Laid {
  readonly span: HTMLElement;
}

/** Text widths for labels, measured once each in the layer's font. */
let measure: CanvasRenderingContext2D | null = null;
function textWidth(text: string, font: string): number {
  measure ??= document.createElement("canvas").getContext("2d");
  if (!measure) return text.length * 6;
  measure.font = font;
  return measure.measureText(text).width;
}

/**
 * A label layer over a labeled track: our own elements, not Nightingale's, pointer events passing through. Each label
 * sits in its feature's widest visible fragment at the range the track shows, read from the track's own state when its
 * display range, length or height changes (attributes the manager sets) or its size does, then placed in the next
 * frame, with Nightingale's redraw; only positions change, nothing is rebuilt.
 */
function labelLayer(
  node: Element,
  layer: HTMLElement,
  features: readonly LabelFeature[],
  all: Map<unknown, FeatureLabel>,
) {
  const labels: FeatureLabel[] = [];
  for (const feature of features) {
    const text = feature.svs.label;
    if (!text || !feature.svs.inside) continue;
    const span = document.createElement("span");
    span.className = feature.svs.layer < 0 ? "svs-feature-label svs-block-label" : "svs-feature-label";
    span.textContent = text;
    const rect = (feature.shape ?? "rectangle") === "rectangle";
    span.style.color = rect ? inkOn(feature.fill ?? feature.color, feature.opacity || 0.9) : INK.dark;
    span.style.display = "none";
    layer.appendChild(span);
    const fragments = feature.locations?.flatMap((l) => l.fragments) ?? [{ start: feature.start, end: feature.end }];
    const box: LabelBox = { mode: "none", left: 0, width: 0 };
    const label = {
      feature,
      text,
      span,
      fragments,
      style: span.style,
      width: -1,
      cut: -1,
      box,
      top: 0,
      height: 0,
      visible: false,
    };
    labels.push(label);
    all.set(feature, label);
  }
  const view: NightingaleAxis = { start: 1, end: 1, left: 0, right: 0 };
  const track = node as unknown as TrackState;
  let frame = 0;
  // The computed font is read only while a label is unmeasured; a label is written only where its box changed.
  let font = "";
  const measureText = (t: string) => {
    font ||= getComputedStyle(labels[0].span).font;
    return textWidth(t, font);
  };
  const place = () => {
    frame = 0;
    if (!layer.isConnected) return;
    layLabels(labels, trackView(view, track), track.layoutObj, measureText);
  };
  const schedule = () => {
    if (!frame) frame = requestAnimationFrame(place);
  };
  const attributes = new MutationObserver(schedule);
  attributes.observe(node, { attributes: true, attributeFilter: ["display-start", "display-end", "length", "height"] });
  const size = new ResizeObserver(schedule);
  size.observe(layer);
  // The track's `width` is the one its xScale draws with: set from its measured size a frame or more after the layer
  // has that size (and after `mount` gave it the container's), with no attribute or event; place again when it is set.
  const unwatch = watchDimensions(track, schedule);
  // A text measured in a fallback font before the page's font loaded is measured again once fonts have loaded.
  const refont = () => {
    font = "";
    for (const l of labels) l.width = -1;
    schedule();
  };
  document.fonts?.addEventListener("loadingdone", refont);
  schedule();
  return {
    schedule,
    dispose() {
      unwatch();
      document.fonts?.removeEventListener("loadingdone", refont);
      attributes.disconnect();
      size.disconnect();
      cancelAnimationFrame(frame);
    },
  };
}

/** Mounts so far: heatmap ids made unique on the page (the heatmap draws into the element of that id). */
let mounts = 0;
/**
 * Our rules over Nightingale's light DOM: a sequence under a ruler without its own axis; contact points one size; and
 * groups (materialize's): a collapsed group's body and an expanded group's merged row are hidden in place, at their
 * width (`display: none` would measure them at their `min-width` and redraw them on showing), so a toggle redraws
 * nothing.
 */
function addStyle() {
  if (document.getElementById("svs-nightingale-style")) return;
  const style = document.createElement("style");
  style.id = "svs-nightingale-style";
  // nightingale-links sizes its points by the base width (dots when zoomed out, discs on short sequences): fixed.
  style.textContent = [
    `nightingale-sequence.${AXISLESS} g.axis { display: none; }`,
    "nightingale-links circle.contact-point { r: 4px; }",
    `nightingale-links[style*="${PAIR_COLOR}"] circle.contact-point { fill: var(${PAIR_COLOR}); }`,
    // Nightingale draws a variation track's right-hand residue axis inside the plot (by twice the right margin); the
    // left one, in its gutter, says the same.
    "nightingale-variation-canvas g.variation-y-right, nightingale-variation g.variation-y-right { display: none; }",
    // A feature track's gutters are rects over its plot, flush with the axis ends: they hid the outer half of the
    // 1 px outline of a feature at either end (a block's bar, a range from position 1). Moved out by half that, so the
    // outline shows whole; the features stay where the axis puts them.
    ":is(nightingale-track, nightingale-track-canvas) g.margin rect.margin-left { transform: translateX(-0.5px); }",
    ":is(nightingale-track, nightingale-track-canvas) g.margin rect.margin-right { transform: translateX(0.5px); }",
    ".svs-group[data-collapsed] > .svs-group-body { height: 0; overflow: hidden; visibility: hidden; }",
    ".svs-group:not([data-collapsed]) > .svs-group-header > .svs-layers > * { visibility: hidden; }",
    // Elements stacked in one place: letters over their colored bar.
    // One column of the cell's width (not its widest element's: they start at the container's width until measured).
    `.${STACK} { display: grid; grid-template-columns: minmax(0, 1fr); } .${STACK} > * { grid-area: 1 / 1; }`,
    // Over a bar, the letters' own striped ground (drawn where they fit) would hide it; hover is by position.
    `.${STACK} > nightingale-sequence rect.base_bg { display: none; }`,
    // A labeled track and its label layer, over it; pointer events pass through to Nightingale's features.
    ".svs-labeled { position: relative; }",
    ".svs-label-layer { position: absolute; inset: 0; overflow: hidden; pointer-events: none; }",
    ".svs-feature-label { position: absolute; box-sizing: border-box; padding: 0 3px; overflow: hidden;" +
      " text-align: center; text-overflow: ellipsis; white-space: nowrap; font-size: 10px; }",
  ].join("\n");
  document.head.appendChild(style);
}

/**
 * Mounts `root` in `container`. `source` is the resolved view and set `materialize` drew: `highlight` goes through its
 * projection. Required: without it every host highlight lit nothing, silently (the Nightingale stories, from SW1 on).
 */
export function mount(
  root: ElementDesc,
  container: HTMLElement,
  source: { view: ResolvedView; set: SeqSet },
  onEvent?: (e: NightingaleEvent) => void,
): Mounted {
  const managers = new Map<Element, AxisMeta>();
  /** Every laid label by its feature (hover hints), and each layer's placement. */
  const labels = new Map<unknown, FeatureLabel>();
  const layers: { schedule(): void; dispose(): void }[] = [];
  const pending: Promise<Issue[]>[] = [];
  const abort = new AbortController();
  const page = ++mounts;
  addStyle();
  // Nightingale's elements draw first at their `min-width` (below their margins: negative widths, and a navigation
  // brush that dispatches a one-position range) until a ResizeObserver measures them; start them at the container's,
  // then at their own cell's once attached (below), before their first draw.
  const width = container.clientWidth;
  const build = (d: ElementDesc): Element => {
    const node = document.createElement(d.tag);
    if (width > 0 && d.tag.startsWith("nightingale-") && d.tag !== "nightingale-manager")
      (node as unknown as { width: number }).width = width;
    for (const [k, v] of Object.entries(d.attributes)) node.setAttribute(k, k === "heatmap-id" ? `${v}-m${page}` : v);
    for (const c of d.children) node.appendChild(build(c));
    if (d.tag === "nightingale-manager" && d.properties.svs) managers.set(node, d.properties.svs as AxisMeta);
    const props = Object.entries(d.properties).filter(([k]) => k !== "svs");
    const lazy = d.tag === "nightingale-linegraph-track" ? (d.properties.svs as StreamedValues | undefined) : null;
    // Our own elements (labels, toggles) take their text at once; Nightingale's their data once rendered.
    let drawn: Promise<unknown> | null = null;
    if (!d.tag.includes("-")) for (const [k, v] of props) (node as unknown as Record<string, unknown>)[k] = v;
    else if (props.length || d.calls.length) {
      const data = apply(node, d, props).then(() => (lazy ? fill(node, d, lazy) : []));
      pending.push(data);
      drawn = data;
    }
    if (!(LABELED in d.attributes)) return node;
    // A labeled track: wrapped with its label layer, laid out again once its data has drawn (its lanes known).
    const wrap = document.createElement("div");
    wrap.className = "svs-labeled";
    const layer = document.createElement("div");
    layer.className = "svs-label-layer";
    layer.setAttribute("aria-hidden", "true");
    wrap.append(node, layer);
    const placed = labelLayer(node, layer, d.properties.data as LabelFeature[], labels);
    layers.push(placed);
    drawn?.then(() => placed.schedule());
    return wrap;
  };
  let disposed = false;
  const apply = async (node: Element, d: ElementDesc, props: [string, unknown][]) => {
    if (d.tag.includes("-")) await customElements.whenDefined(d.tag);
    await (node as { updateComplete?: Promise<unknown> }).updateComplete;
    if (disposed) return;
    const target = node as unknown as Record<string, unknown>;
    for (const [k, v] of props) target[k] = v;
    for (const c of d.calls) (target[c.method] as (...a: unknown[]) => void)?.apply(node, c.args);
  };
  /** A streamed linegraph: its bins fetched for the element's width, then its one series' values set. */
  const fill = async (node: Element, d: ElementDesc, p: StreamedValues): Promise<Issue[]> => {
    const length = Number(d.attributes.length);
    const r = await fetchStreamed(p, length, width, abort.signal);
    if (disposed) return [];
    const [series] = d.properties.data as { values: { position: number; value: number | null }[] }[];
    (node as unknown as { data: unknown }).data = [{ ...series, values: mergeSeries(series.values, r.values) }];
    return r.issues;
  };
  const top = build(root);
  /** The arbiter between a component's item hover and the pointer's, the pointer's hover, and scratch reused per event. */
  const arbiter = { onItem: false };
  const hover: HoverTarget = { manager: null, position: -1, layer: null };
  const where: Located = { block: -1, index: -1 };
  const [axisRuns, blockRuns] = [AxisRuns.create(), BlockRuns.create()];
  /** What each manager shows (`events.ts`), made on its first use. */
  const lights = new Map<Element, Light>();
  const light = (m: Element) => lights.get(m) ?? lights.set(m, newLight()).get(m) ?? newLight();

  /** Scratch for a hover's highlight and the section's kept under it. */
  const keptRuns = AxisRuns.create();
  /**
   * What `manager` shows after Nightingale's own `highlight`: with the section's `highlights` added (sent again only
   * where they were not in it, once per hover change).
   */
  const keptUnder = (manager: Element, meta: AxisMeta, highlight: string | undefined): string | null => {
    const kept = meta.highlights;
    if (!kept || kept.count === 0) return highlight || null;
    keptRuns.count = 0;
    keptRuns.sorted = true;
    AxisRuns.parse(highlight ?? "", keptRuns);
    const text = AxisRuns.text(keepLit(keptRuns, kept));
    if (text !== highlight) manager.dispatchEvent(new CustomEvent("change", { detail: { highlight: text } }));
    return text;
  };
  const listener = (e: Event) => {
    const detail = ((e as CustomEvent).detail ?? {}) as Detail;
    const target = e.target as Element;
    const manager = target.closest("nightingale-manager");
    const meta = manager && managers.get(manager);
    if (!manager || !meta || target === manager) return;
    // Nightingale's own highlight (a hover's) reached the manager, which showed it on every element: the section's
    // `highlights` are put back under it (its leave clears all), and what it now shows is noted for `highlight`.
    if ("highlight" in detail)
      noteLit(light(manager), keptUnder(manager, meta, detail.highlight as string | undefined));
    if ("display-start" in detail && meta.scrollGroup !== null) {
      // A change on the other managers, as their own children send one: an attribute set on a manager would be
      // re-applied over every later zoom of its own (NG-01). Its target is the manager: this listener skips it.
      for (const [m, other] of managers) {
        if (m === manager || other.scrollGroup !== meta.scrollGroup) continue;
        const range = { "display-start": detail["display-start"], "display-end": detail["display-end"] };
        m.dispatchEvent(new CustomEvent("change", { detail: range }));
      }
      return;
    }
    const kind = KIND[String(detail.eventType ?? detail.eventtype ?? e.type)] ?? (detail.highlight ? "hover" : null);
    if (!onEvent || !kind) return;
    const svs = detail.feature?.svs as { layer: number; annotation: number; item: number } | undefined;
    // Hover by position is the pointer's (`move`); a component's own hover wins only for an item (a feature, a
    // variant, a contact), so one move is one hover. Its leave hands hover back to the pointer.
    const item = !!svs || target.tagName === "NIGHTINGALE-LINKS";
    if (!arbitrate(arbiter, kind, item)) return;
    // An item's hover holds the pointer's hover (its manager, no position): leaving the viewer from the item, whose
    // own leave comes after the pointer's (mouse events follow pointer events), still says leave once.
    if (kind !== "click") {
      hover.manager = kind === "hover" ? manager : null;
      hover.position = -1;
    }
    // A group's merged row has no one layer: its features name theirs.
    const own = Number(target.closest("[data-svs-layer]")?.getAttribute("data-svs-layer") ?? -1);
    const layer = svs && svs.layer >= 0 ? svs.layer : own;
    const msa = e.type.startsWith("onResidue") ? Number(detail.position) + 1 : -1;
    const series =
      target.tagName === "NIGHTINGALE-LINEGRAPH-TRACK" ? (target as unknown as LineGraph).data?.[0]?.name : undefined;
    const position = msa > 0 ? msa : positionOf(detail, series);
    const at = kind !== "leave" && locate(meta, position - 1, where) ? where : null;
    const [annotation, index] = svs && svs.layer >= 0 ? [svs.annotation, svs.item] : [-1, -1];
    const block = at?.block ?? (svs && svs.layer < 0 ? svs.item : -1);
    // A block label (layer -1) names its block, not a range of it.
    const extent =
      kind === "leave" || (svs && svs.layer < 0)
        ? []
        : BlockRuns.toArray(extentOf(meta, rangesOf(detail, position, axisRuns), blockRuns));
    // A label not shown in full (cut, no room, or not laid inside: `labels` off) is the hover's hint, at the feature.
    const said = kind === "hover" ? (detail.feature?.svs as LabelFeature["svs"] | undefined) : undefined;
    const text = said?.hint ?? said?.label;
    const full = labels.get(detail.feature)?.box.mode === "full";
    const anchor = (detail.target as Element | undefined) ?? target;
    const hint = text && !full ? { anchor, text } : null;
    onEvent({
      kind,
      section: meta.section,
      block,
      index: at?.index ?? -1,
      layer,
      annotation,
      item: index,
      extent,
      hint,
    });
  };
  const view: NightingaleAxis = { start: 1, end: 1, left: 0, right: 0 };
  /** The cell under the pointer: its manager, an element showing its range, and that element's left edge (read once). */
  interface Cell {
    readonly cell: Element;
    readonly manager: Element;
    readonly element: Element & TrackState;
    readonly left: number;
  }
  let cached: Cell | null = null;
  const cellOf = (target: Element): Cell | null => {
    const cell = target.closest?.(".svs-layers");
    if (!cell) return null;
    if (cached?.cell === cell) return cached;
    const manager = cell.closest("nightingale-manager");
    const element = cell.querySelector("[length]") as (Element & TrackState) | null;
    cached = manager && element ? { cell, manager, element, left: element.getBoundingClientRect().left } : null;
    return cached;
  };
  /** Scrolling and resizing move the cell: its left edge is read again (and when the pointer comes back, `leave`). */
  const moved = () => {
    cached = null;
  };
  const resized = new ResizeObserver(moved);
  resized.observe(container);
  document.addEventListener("scroll", moved, { capture: true, passive: true });
  /**
   * Hover by position on every row at any zoom: the pointer's x in a row's data cell through the range its elements
   * show (their display range and gutters, one per manager) to an axis position, then to block
   * and index; emitted when the position changes, not while a component's item is hovered. The position, manager and
   * layer are compared before anything is located (`hoverAt`), so a move within one position allocates nothing.
   */
  const move = (e: Event) => {
    if (!onEvent || arbiter.onItem) return;
    const cell = cellOf(e.target as Element);
    const meta = cell && managers.get(cell.manager);
    if (!cell || !meta) return away();
    const layerEl = (e.target as Element).closest("[data-svs-layer]");
    const x = (e as PointerEvent).clientX - cell.left;
    const result = hoverAt(hover, cell.manager, meta, trackView(view, cell.element), x, layerEl, where);
    if (result === HoverResult.off) return away();
    if (result === HoverResult.same) return;
    const layer = Number(layerEl?.getAttribute("data-svs-layer") ?? -1);
    const extent = [{ block: where.block, start: where.index, stop: where.index + 1 }];
    const { block, index } = where;
    onEvent({
      kind: "hover",
      section: meta.section,
      block,
      index,
      layer,
      annotation: -1,
      item: -1,
      extent,
      hint: null,
    });
  };
  /**
   * The pointer left the viewer: nothing is hovered (a component's item whose leave has not come yet included), and the
   * cell's edge is read again on its return (the page may move the viewer meanwhile, with no scroll or resize).
   */
  const leave = () => {
    arbiter.onItem = false;
    cached = null;
    away();
  };
  /** The pointer left the rows' data (or a position off the axis): the pointer's hover ends. */
  const away = () => {
    if (!hover.manager || !onEvent) return;
    const section = managers.get(hover.manager as Element)?.section ?? -1;
    hover.manager = null;
    hover.position = -1;
    onEvent({
      kind: "leave",
      section,
      block: -1,
      index: -1,
      layer: -1,
      annotation: -1,
      item: -1,
      extent: [],
      hint: null,
    });
  };
  /** A group's header toggles it in place: only attributes change, every element keeps its data and drawing. */
  const toggle = (e: Event) => {
    const button = (e.target as Element).closest?.(".svs-toggle");
    const group = button?.closest(".svs-group");
    if (!button || !group) return;
    const collapsed = group.toggleAttribute("data-collapsed");
    button.setAttribute("aria-expanded", String(!collapsed));
    const mark = button.querySelector(".svs-disclosure");
    if (mark) mark.textContent = collapsed ? "▸" : "▾";
  };
  const types = ["change", "onResidueMouseEnter", "onResidueClick"];
  // The shared projection of the mount's view: its routes are cached on it; its extents, sink and runs are reused.
  const projection = projectionOf(source.view, source.set);
  const [ext, sink, shownRuns] = [Extents.create(16), new AxisSink(), AxisRuns.create()];
  for (const t of types) top.addEventListener(t, listener);
  top.addEventListener("click", toggle);
  top.addEventListener("pointermove", move);
  top.addEventListener("pointerleave", leave);
  container.appendChild(top);
  // The container is wider than a track's cell (the row labels, the card's padding): an element first drawn at the
  // container's width is drawn too wide until its own ResizeObserver sets the measured one, a frame or more later (its
  // features and our labels past the ruler's end). Each takes its cell's width now, before its first update.
  for (const el of top.querySelectorAll("*")) {
    if (!el.localName.startsWith("nightingale-") || el.localName === "nightingale-manager") continue;
    const target = el as unknown as { width?: number; "min-width"?: number };
    const cell = el.getBoundingClientRect().width;
    if (cell > 0 && target.width !== undefined && el.getAttribute("width") === null)
      target.width = Math.max(cell, target["min-width"] ?? 0);
  }

  return {
    ready: Promise.all(pending).then((all) => all.flat()),
    highlight(loci, opts) {
      const lit = prepareLoci(projection, ext, loci, opts?.closure ?? "element");
      for (const [m, meta] of managers) {
        shownRuns.count = 0;
        if (lit) highlightRuns(projection, ext, meta, sink, shownRuns);
        // The section's `highlights` stay lit under the host's.
        keepLit(shownRuns, meta.highlights);
        // Only a changed highlight goes to a manager (and its string is built then): each change re-applies its
        // attributes to every element.
        const highlight = nextHighlight(light(m), shownRuns);
        if (highlight !== undefined)
          m.dispatchEvent(new CustomEvent("change", { detail: { highlight }, bubbles: true }));
      }
    },
    dispose() {
      disposed = true;
      abort.abort();
      for (const l of layers) l.dispose();
      for (const t of types) top.removeEventListener(t, listener);
      top.removeEventListener("click", toggle);
      top.removeEventListener("pointermove", move);
      top.removeEventListener("pointerleave", leave);
      document.removeEventListener("scroll", moved, { capture: true });
      resized.disconnect();
      for (const msa of top.querySelectorAll("nightingale-msa")) (msa as { worker?: Worker }).worker?.terminate();
      top.remove();
    },
  };
}
