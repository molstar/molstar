/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the layout, hit testing and text tests (LY1): the format's acceptance cases through its fixture loaders,
// stub renderers with heights by type (a split layer prepares its values' rows, as `Features.split`), the fixture text
// measurer, and a manual clock; helpers that find a point on a section from the layout (`t.xOf` and `t.yOf`).

import { type Loaders, type SeqViewSpecDoc, splitRows } from "@molstar/seq-core/format";
import type { Annotation } from "@molstar/seq-core/model";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { fixtureLoaders, MolstarFixture } from "../../core/_test/format/load-fixtures.js";
import type { ShownIndex } from "../src/axis.js";
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  type LayerInput,
  type LayerRenderer,
  manualScheduler,
  type PresetEntry,
  type RenderProps,
  type SeqViewerContext,
  type SeqViewerRegistryEntry,
} from "../src/index.js";
import type { Layouts } from "../src/layout.js";
import { PreparedParts } from "../src/prepare.js";
import { type LaneSegments, laneLineHeights } from "../src/runs.js";
import { fixtureMeasurer } from "../src/testing/index.js";
import type { Viewport } from "../src/viewport.js";
/** A type's own height (free sizing, and every type under `natural`). */
export const HEIGHTS: Record<string, number> = {
  letters: 26,
  ruler: 30,
  ranges: 14,
  points: 18,
  values: 40,
  pairs: 48,
  color: 12,
  rows: 22,
  variants: 18,
  secondary: 14,
};
const TYPES = Object.keys(HEIGHTS);

/**
 * A split layer's values and each data entry's value per row (seqspec:format#layer-types), as `Features.split` holds
 * them.
 */
function split(input: LayerInput) {
  const sp = input.params.split as { by: string } | null | undefined;
  if (!sp) return null;
  const anns = input.layer.data.map((id) => input.set.annotations[id] as Annotation);
  const { values, rowOf } = splitRows(input.set, anns, input.layer.kept, sp);
  return { labels: values, rows: values.map(() => ({})), rowOf };
}

export function renderers(): LayerRenderer[] {
  return TYPES.map((name) => ({
    name,
    geometries: ["linear", "wrapped"],
    prepare: (input) => ({ ref: input.layer.ref, split: split(input) }),
    height: () => HEIGHTS[name],
    splitHeight: name === "variants" ? 8 : undefined,
  }));
}

/** The kernel types' prepares (`PreparedParts`, LY2) with the fixture heights: prepared state with no draws. */
export function preparedRenderers(): LayerRenderer[] {
  return Object.entries(PreparedParts).map(([name, part]) => ({
    name,
    geometries: ["linear", "wrapped"],
    ...part,
    height: () => HEIGHTS[name] ?? 20,
    splitHeight: name === "variants" ? 8 : undefined,
  }));
}

/** A layer's index, prepared state and track, by its track's label (the first layer of that track). */
export function layer(ctx: SeqViewerContext, label: string) {
  const m = ctx.graph.assembled.value?.view;
  if (!m) throw new Error("no view");
  const v = m.resolved;
  for (let k = 0; k < v.layers.count; k++) {
    const t = v.layers.track[k];
    if (m.set.strings.get(v.tracks.label[t]) === label) return { k, t, state: m.prepared[k].state, m };
  }
  throw new Error(`no track "${label}"`);
}
/**
 * Ranges prepared as lanes of features over block coordinate system indices (`LaneSegments`, LY2's `Features` in
 * small): packed greedily into lanes per block, so a wrapped line's band counts the lanes with an item on it
 * (`laneLineHeights`).
 */
export function laneRanges(): LayerRenderer {
  return {
    name: "ranges",
    geometries: ["linear", "wrapped"],
    prepare: (input): LaneSegments => {
      const blocks = input.section.blocks.length;
      const ann = input.set.annotations[input.layer.data[0]] as Annotation;
      // Data on the block's own coordinate system (the tests' is), packed by start.
      const items = input.layer.placements.flatMap((p) =>
        Array.from({ length: ann.count }, (_, r) => ({
          block: p.block,
          start: ann.start[r],
          stop: ann.stop[r],
          lane: 0,
        })),
      );
      items.sort((x, y) => x.block - y.block || x.start - y.start);
      const ends: number[][] = Array.from({ length: blocks }, () => []);
      for (const f of items) {
        const e = ends[f.block];
        f.lane = e.findIndex((x) => x <= f.start);
        if (f.lane < 0) f.lane = e.push(0) - 1;
        e[f.lane] = f.stop;
      }
      const lanes = Math.max(1, ...ends.map((e) => e.length));
      items.sort((x, y) => x.block - y.block || x.lane - y.lane || x.start - y.start);
      const seg = new Int32Array(blocks * lanes + 1);
      for (let g = 0, k = 0; g <= blocks * lanes; g++) {
        while (k < items.length && items[k].block * lanes + items[k].lane < g) k++;
        seg[g] = k;
      }
      const reach = new Int32Array(items.length);
      for (let g = 0; g < blocks * lanes; g++)
        for (let k = seg[g], r = -(2 ** 31); k < seg[g + 1]; k++) reach[k] = r = Math.max(r, items[k].stop);
      return {
        lanes,
        seg,
        start: Int32Array.from(items, (f) => f.start),
        stop: Int32Array.from(items, (f) => f.stop),
        reach,
      };
    },
    lineHeights: laneLineHeights,
  };
}

/** Presets with geometry only: msa with a navigator and a vertical overview. */
export const PRESETS: PresetEntry[] = [
  { name: "tracks", geometry: "linear", navigator: { height: 20 } },
  {
    name: "msa",
    geometry: "linear",
    render: { rowHeight: 18, verticalOverview: 24, verticalNav: "both", verticalZoom: true },
    navigator: { height: 20 },
  },
  {
    name: "msa-fixed",
    geometry: "linear",
    render: { rowHeight: 18, verticalOverview: 24, verticalNav: "both", verticalZoom: false },
    navigator: { height: 20 },
  },
  { name: "wrapped", geometry: "wrapped", render: { headerWidth: 0 } },
  { name: "compact", geometry: "wrapped", render: { headerWidth: 0, wrappedFontPx: 13, wrappedGutter: 8 } },
  { name: "bare", geometry: "linear" },
  { name: "railed", geometry: "linear", navigator: "rulerRail" },
];

export interface Options {
  size?: Record<string, readonly [number, number]>;
  view?: string;
  render?: Partial<RenderProps>;
  loaders?: Loaders;
  registry?: SeqViewerRegistryEntry[];
  /** In place of the stubs. */
  renderers?: LayerRenderer[];
  presets?: boolean;
}

/** A context showing `doc`, every viewport sized (1000 × 400 unless given) before the graph settles. */
export async function opened(doc: SeqViewSpecDoc, o: Options = {}) {
  const scheduler = manualScheduler();
  const measurer = fixtureMeasurer();
  const spec = DefaultSeqViewerSpec();
  const entries: SeqViewerRegistryEntry[] = [{ svs: MolstarFixture }, { renderers: o.renderers ?? renderers() }];
  if (o.presets !== false) entries.push({ presets: PRESETS });
  const ctx = createSeqViewer(
    // The fixture presets (geometry and a navigator row) in place of the built-ins.
    {
      ...spec,
      registry: [{ ...DefaultRegistry, presets: undefined, renderers: undefined }, ...entries, ...(o.registry ?? [])],
      render: o.render,
    },
    { scheduler, loaders: o.loaders ?? fixtureLoaders(), text: measurer },
  );
  const size = () => {
    for (const ref of ctx.viewports.sections.value) {
      const [w, h] = o.size?.[ref] ?? [1000, 400];
      ctx.viewports.get(ref)?.setSize(w, h, 1);
    }
  };
  const done = ctx.setDocument(doc);
  size();
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  if (o.view) {
    ctx.views.setCurrent(o.view);
    await ctx.graph.update();
  }
  return { ctx, scheduler, measurer, layouts: ctx.layout as Layouts, ...helpers(ctx) };
}
export type Opened = Awaited<ReturnType<typeof opened>>;

export function helpers(ctx: SeqViewerContext) {
  const vp = (section: string) => ctx.viewports.get(section) as Viewport;
  const layout = (section: string) => {
    const l = ctx.layout.get(section);
    if (!l) throw new Error(`no layout for ${section}`);
    return l.current;
  };
  const trackIndex = (ref: string) => {
    const a = ctx.graph.assembled.value;
    const v = a?.view?.resolved;
    for (let t = 0; v && a && t < v.tracks.count; t++) if (a.set.strings.get(v.tracks.ref[t]) === ref) return t;
    throw new Error(`no track ${ref}`);
  };
  /** The center x of a block index of a block, from the axis view. */
  const xOf = (section: string, block: number, index: number) => {
    const v = ctx.layout.axis(section);
    if (!v) throw new Error(`no axis for ${section}`);
    const k = v.shownOf(block, index);
    return (v.x(block, k) + v.x(block, (k + 1) as ShownIndex)) / 2;
  };
  /** The center y of a track's entry (its `row`th display row on a rows track), in surface px. */
  const yOf = (section: string, track: string, row = -1) => {
    const l = layout(section);
    const t = trackIndex(track);
    const s = vp(section).live.current;
    for (let e = 0; e < l.entryCount; e++) {
      if (l.entryTrack[e] !== t) continue;
      const lazy = l.entryRun[e] > 0;
      if (row >= 0 && !lazy && l.entryRow[e] !== row) continue;
      const h = lazy ? l.rows.height : l.entryHeight[e];
      const top = l.entryTop[e] + (lazy && row >= 0 ? row * h : 0);
      return (e < l.pinnedCount ? top : l.scrolledTop + top - s.scrollY) + h / 2;
    }
    throw new Error(`no entry for ${track}`);
  };
  /** Entries as [kind, track ref or group or block, row, top, height], for readable checks. */
  const entries = (section: string) => {
    const l = layout(section);
    const a = ctx.graph.assembled.value;
    const v = a?.view?.resolved;
    const out: [number, string, number, number, number][] = [];
    for (let e = 0; v && a && e < l.entryCount; e++) {
      const t = l.entryTrack[e];
      const name =
        t >= 0
          ? a.set.strings.get(v.tracks.ref[t])
          : l.entryGroup[e] >= 0
            ? a.set.strings.get(v.groups.ref[l.entryGroup[e]])
            : `block ${l.entryBlock[e]}`;
      out.push([l.entryKind[e], name, l.entryRow[e], l.entryTop[e], l.entryHeight[e]]);
    }
    return out;
  };
  return { vp, layout, trackIndex, xOf, yOf, entries };
}

/** An annotation's item rows by key, for readable checks. */
export const itemKey = (ctx: SeqViewerContext, id: number) => {
  const set = ctx.graph.assembled.value?.set;
  const a = set?.annotations[id];
  return set && a ? set.strings.get(a.key) : "";
};

/** Loaders serving one procedural alignment of `rows` × `columns` by view (`seqstar-tiles/1`). */
export function lazyLoaders(uri: string, rows: number, columns: number): Loaders {
  const server = tileServer({ rowCount: rows, columnCount: columns });
  return {
    load: (r) => Promise.reject(new Error(`no fixture for ${r.uri ?? r.hostFile}`)),
    stream: tilesStream((u) => (u === uri ? server.transport : null)),
  };
}
