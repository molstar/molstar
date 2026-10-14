/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// XR1b (seqspec:viewer#extension-points; seqspec:format#extensions): the conformance kit over every kernel entry's
// headless half (hits sampled, legends, fallbacks without the entry and in a geometry it does not declare) and over
// the `org.example:stars` example drawn on the batch API (the same batches on a tee, hits back to their items, the
// legend); field providers read wherever a field is (a `ranges` layer colored by `org.example:affinity`, a `rules`
// measure); and a host's own `ranges` replacing the built-in on case 13.

import { numberColumn, type SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it, vi } from "vitest";
import { case3, case13 } from "../../core/_test/format/cases.ts";
import type { ShownIndex } from "../src/axis.ts";
import { band, type DrawView, type Features, hitItem, PreparedParts, RectBatch } from "../src/draw/api.ts";
import { Hit, HoverKind, KernelRenderers, type LayerRenderer, type SeqViewerContext } from "../src/index.ts";
import { DefaultPresets } from "../src/presets.ts";
import { RangesRenderer } from "../src/renderers/ranges.ts";
import {
  AffinityField,
  type ConformanceFixtures,
  EdgeKinds,
  edgeDocument,
  rendererConformance,
  StarsFormat,
  StarsRenderer,
} from "../src/testing/index.ts";
import { pinsDoc } from "./boundary-pins-fixtures.ts";
import {
  alignedDoc,
  asText,
  directed,
  edgeTypes,
  mapped,
  matrixDoc,
  miniatureDoc,
  ringed,
  starsDoc,
} from "./conformance-fixtures.ts";
import { dotPlotDoc } from "./dot-plot-fixtures.ts";
import { opened } from "./layout-fixtures.ts";
import { allTypes } from "./prepare-fixtures.ts";

const without = (name: string) => KernelRenderers.filter((r) => r.name !== name);
/** The built-in `circular` preset beside the fixture presets: sections on rings (CG1b). */
const RINGS = { presets: DefaultPresets.filter((p) => p.name === "circular") };
async function open(doc: SeqViewSpecDoc, renderers: LayerRenderer[], extra: LayerRenderer[] = []) {
  const t = await opened(doc, { renderers: [...renderers, ...extra], registry: [{ svs: StarsFormat }, RINGS] });
  return t.ctx;
}
/** The kit's fixtures for a kernel entry: its document in tracks or in the text layout, with or without the entry. */
function kernel(entry: LayerRenderer, doc: () => SeqViewSpecDoc, reversed = false, edges = false): ConformanceFixtures {
  return {
    reversed,
    edges,
    open: async ({ geometry, entry: has, direction, edge }) => {
      const d = edge ? edgeDocument(doc(), edge) : direction ? directed(doc(), direction) : doc();
      const as = geometry === "wrapped" ? asText(d) : geometry === "circular" ? ringed(d) : d;
      const ctx = await open(as, has ? [...KernelRenderers] : without(entry.name));
      // Wrapped lines run past 400 px: tall enough that every line drawn can be hit (RN2: underlines at a line's foot).
      if (geometry === "wrapped") {
        for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1000, 1400, 1);
        await ctx.graph.update();
      }
      return ctx;
    },
    // Letters hit residues (`of`, `row`), not items; secondary's own `hit` is its wrapped underlines' (each an item), its
    // coil line linear and on rings the position's.
    items:
      entry.name === "letters"
        ? false
        : entry.name === "secondary"
          ? { wrapped: true, linear: false, circular: false }
          : undefined,
  };
}
const DOCS: Record<string, () => SeqViewSpecDoc> = {
  rows: case3,
  logo: case3,
  matrix: matrixDoc,
  dot_plot: () => dotPlotDoc(),
  "seqv:map": miniatureDoc,
};
/** The edge geometries' documents, where the kit's own differs: an alignment given by value for rows and logo. */
const EDGE_DOCS: Record<string, () => SeqViewSpecDoc> = { ...DOCS, rows: alignedDoc, logo: alignedDoc };
/** What an alignment's columns cannot take. */
const ALIGNED_SKIPS = {
  origin: "an alignment's columns have no cyclic origin",
  stacked: "an axis over an alignment is its one block of columns",
  copies: "an axis over an alignment is its one block of columns",
};
/** Reversed, a dot plot mirrors down the track too (its y is the axis): the kit's mirror is across only (dot-plot.test.ts). */
const TWO_AXES = new Set(["dot_plot"]);
const failures = (checks: Awaited<ReturnType<typeof rendererConformance>>) =>
  checks.filter((c) => !c.ok).map((c) => `${c.check} ${c.geometry}: ${c.detail}`);

describe("the conformance kit over the kernel entries (headless halves; RN1's draws on batches)", () => {
  for (const entry of KernelRenderers)
    it(`${entry.name}: hits, legend and fallbacks`, async () => {
      const checks = await rendererConformance(entry, kernel(entry, DOCS[entry.name] ?? (() => allTypes())));
      expect(failures(checks)).toEqual([]);
      // Every geometry checked its fallback; the declared ones their hits and legend.
      expect(new Set(checks.map((c) => `${c.check} ${c.geometry}`))).toEqual(
        new Set([
          "fallback linear",
          "fallback wrapped",
          ...(entry.geometries.includes("circular") ? [] : ["fallback circular"]),
          ...entry.geometries.flatMap((g) => [`hits ${g}`, `legend ${g}`, ...(entry.draw ? [`batches ${g}`] : [])]),
        ]),
      );
    });
});

describe("the conformance kit over reversed and ranged pieces (AX1)", () => {
  for (const entry of KernelRenderers)
    it(`${entry.name}: every check on reversed pieces, and the reversed draw the forward one mirrored`, async () => {
      // Reversed and ranged pieces through every check; then forward against reversed, mirrored.
      const doc = DOCS[entry.name] ?? (() => allTypes());
      const checks = [
        ...(await rendererConformance(
          entry,
          kernel(entry, () => directed(doc(), -1, true)),
        )),
        // At the fit: a short document, so letters and labels are drawn where every position is in view.
        ...(
          await rendererConformance(entry, kernel(entry, DOCS[entry.name] ?? (() => allTypes({ length: 50 })), true))
        ).filter((c) => c.check === "mirror" && !TWO_AXES.has(entry.name)),
      ];
      expect(failures(checks)).toEqual([]);
      if (entry.draw && entry.geometries.includes("linear") && !TWO_AXES.has(entry.name))
        expect(checks.filter((c) => c.check === "mirror").length).toBeGreaterThan(0);
    });
});

describe("the conformance kit over the shared edge geometries (reversed, a cyclic origin, stacked, filtered, copies)", () => {
  for (const entry of KernelRenderers)
    it(`${entry.name}: drawn and hit back in each`, async () => {
      const aligned = EDGE_DOCS[entry.name] === alignedDoc;
      const fx = { ...kernel(entry, EDGE_DOCS[entry.name] ?? edgeTypes, false, true) };
      const checks = await rendererConformance(entry, aligned ? { ...fx, edgesSkipped: ALIGNED_SKIPS } : fx);
      const edges = checks.filter((c) => c.check === "edges");
      expect(failures(edges)).toEqual([]);
      // Reversed, the forward draw mirrored (a dot plot mirrors down the track too).
      const mirror = checks.filter((c) => c.check === "mirror");
      if (!TWO_AXES.has(entry.name)) expect(failures(mirror)).toEqual([]);
      expect(mirror.length > 0).toBe(entry.draw !== undefined);
      // Each geometry checked, or skipped by name.
      for (const edge of EdgeKinds)
        expect(edges.filter((c) => c.detail.startsWith(`${edge}:`)).length, edge).toBeGreaterThan(
          aligned && edge in ALIGNED_SKIPS ? 0 : 1,
        );
    });
});

describe("the edge kit fails what it does not check (LC-14 to LC-16)", () => {
  const rows = KernelRenderers.find((r) => r.name === "rows") as LayerRenderer;
  it("a geometry its document keeps out fails as not applied, unless the fixture skips it by name", async () => {
    // Case 3 reads its alignment from a file: no rewrite applies but `reversed`.
    const checks = (await rendererConformance(rows, kernel(rows, case3, false, true))).filter(
      (c) => c.check === "edges",
    );
    const kept = checks.filter((c) => !c.ok && c.detail.includes("the geometry did not apply"));
    expect(kept.map((c) => c.detail.split(":")[0]).sort()).toEqual(["copies", "filtered", "origin", "stacked"]);
  });
  it("copies show two blocks of one unit", async () => {
    const ctx = await open(edgeDocument(edgeTypes(), "copies"), [...KernelRenderers]);
    // allTypes' two sequences, the first shown twice.
    expect(ctx.layout.axis("v.main")?.axis?.count).toBe(3);
    ctx.dispose();
  });
  it("a renderer that draws nothing once the view is stacked, filtered or past an origin fails `edges`", async () => {
    const points = KernelRenderers.find((r) => r.name === "points") as LayerRenderer;
    const draw = points.draw as NonNullable<LayerRenderer["draw"]>;
    const edgeOnly = (v: DrawView) =>
      v.onlyBlock >= 0 ||
      v.positions.some((p) => p !== null) ||
      Array.from(v.first).some((f, b) => f > 0 && v.dir[b] === 1);
    const mutant = {
      ...points,
      draw: (v: DrawView, s: unknown, p: unknown) => (edgeOnly(v) ? undefined : draw(v, s, p)),
    };
    const fx: ConformanceFixtures = {
      edges: true,
      open: async ({ entry: has, edge }) =>
        open(
          edge ? edgeDocument(edgeTypes(), edge) : edgeTypes(),
          without("points"),
          has ? [mutant as LayerRenderer] : [],
        ),
    };
    const checks = await rendererConformance(mutant as LayerRenderer, fx);
    const failed = new Set(checks.filter((c) => c.check === "edges" && !c.ok).map((c) => c.detail.split(":")[0]));
    expect(failed).toEqual(new Set(["origin", "stacked", "filtered"]));
  });
});

describe("the kit's hit round trip and the reversed mirror catch a hit or a draw gone wrong (LC-17)", () => {
  /** The kernel fixture with `entry` registered in place of its built-in, in the edge geometries. */
  const fixture = (entry: LayerRenderer): ConformanceFixtures => ({
    edges: true,
    open: async ({ entry: has, edge }) =>
      open(edge ? edgeDocument(edgeTypes(), edge) : edgeTypes(), without(entry.name), has ? [entry] : []),
  });
  it("a hit naming another item of its data fails `edges`", async () => {
    const ranges = RangesRenderer as LayerRenderer;
    const hit = ranges.hit as NonNullable<LayerRenderer["hit"]>;
    const first: LayerRenderer = {
      ...ranges,
      hit: (v, s, p, out) => {
        if (!hit(v, s, p, out)) return false;
        // The data's first item, wherever the point is.
        out.item = out.itemLoc.item = 0 as never;
        return true;
      },
    };
    const failed = (await rendererConformance(first, fixture(first))).filter((c) => c.check === "edges" && !c.ok);
    expect(failed.length).toBeGreaterThan(0);
  });
  it("a draw that ignores a reversed piece fails the reversed geometry's mirror", async () => {
    const values = KernelRenderers.find((r) => r.name === "values") as LayerRenderer;
    const draw = values.draw as NonNullable<LayerRenderer["draw"]>;
    const forward: LayerRenderer = {
      ...values,
      draw: (v, s, p) => {
        const rev = v.dir.subarray(0, v.axis?.count ?? 0).includes(-1);
        if (!rev) return draw(v, s, p);
        // Drawn as forward: the reversed draw mirrored back within the view.
        const be = v.backend;
        const [l, r] = [Math.round(v.left * v.dpr), Math.round(v.right * v.dpr)];
        v.backend = new Proxy(be, {
          get: (t, key) =>
            key === "rects"
              ? (b: RectBatch) => {
                  for (let i = 0; i < b.n; i++) b.x[i] = l + r - b.x[i] - b.w[i];
                  t.rects(b);
                }
              : (t as never)[key],
        });
        try {
          draw(v, s, p);
        } finally {
          v.backend = be;
        }
      },
    };
    const checks = await rendererConformance(forward, fixture(forward));
    expect(checks.filter((c) => c.check === "mirror" && !c.ok).length).toBeGreaterThan(0);
  });
});

describe("the conformance kit over boundaries pinned over letters", () => {
  it("ranges: every check, forward and reversed, and the reversed draw mirrored", async () => {
    const entry = RangesRenderer as LayerRenderer;
    const checks = [
      ...(await rendererConformance(
        entry,
        kernel(entry, () => pinsDoc()),
      )),
      ...(await rendererConformance(
        entry,
        kernel(entry, () => directed(pinsDoc(), -1, true)),
      )),
      ...(
        await rendererConformance(
          entry,
          kernel(entry, () => pinsDoc(), true),
        )
      ).filter((c) => c.check === "mirror"),
    ];
    expect(failures(checks)).toEqual([]);
    expect(checks.filter((c) => c.check === "mirror").length).toBeGreaterThan(0);
  });
});

describe("the conformance kit over letter maps (LM1)", () => {
  // `letters` and `rows` with a map, linear and wrapped, forward and reversed (and the reversed draw mirrored at the fit).
  for (const name of ["letters", "rows"]) {
    const entry = KernelRenderers.find((r) => r.name === name) as LayerRenderer;
    it(`${name}: every check through a complement map, forward and reversed, and the mirror`, async () => {
      const doc = () => mapped((DOCS[name] ?? (() => allTypes()))());
      const short = () => mapped((DOCS[name] ?? (() => allTypes({ length: 50 })))());
      const checks = [
        ...(await rendererConformance(entry, kernel(entry, doc))),
        ...(await rendererConformance(
          entry,
          kernel(entry, () => directed(doc(), -1, true)),
        )),
        ...(await rendererConformance(entry, kernel(entry, short, true))).filter((c) => c.check === "mirror"),
      ];
      expect(failures(checks)).toEqual([]);
      expect(checks.filter((c) => c.check === "mirror").length).toBeGreaterThan(0);
      expect(checks.some((c) => c.check === "batches" && c.geometry === "wrapped")).toBe(
        entry.geometries.includes("wrapped"),
      );
    });
  }
});

describe("org.example:stars, drawn on the batch API", () => {
  const fixtures = (entry: LayerRenderer): ConformanceFixtures => ({
    open: ({ geometry, entry: has }) =>
      open(
        geometry === "wrapped" ? asText(starsDoc()) : geometry === "circular" ? ringed(starsDoc()) : starsDoc(),
        [...KernelRenderers],
        has ? [entry] : [],
      ),
    legend: ["(ramp)", "affinity (pKd)"],
  });

  it("passes the kit: batches on a tee, hits back to their items, its legend, `points` without it", async () => {
    const checks = await rendererConformance(StarsRenderer as LayerRenderer, fixtures(StarsRenderer as LayerRenderer));
    expect(failures(checks)).toEqual([]);
    // Linear at three zooms, wrapped, and on rings.
    expect(checks.filter((c) => c.check === "batches")).toHaveLength(5);
  });

  it("passes the kit on a long document, below a px a position at the fit (LC-18)", async () => {
    const entry = StarsRenderer as LayerRenderer;
    const long: ConformanceFixtures = {
      ...fixtures(entry),
      open: ({ geometry, entry: has }) =>
        open(
          geometry === "wrapped"
            ? asText(starsDoc({ long: true }))
            : geometry === "circular"
              ? ringed(starsDoc({ long: true }))
              : starsDoc({ long: true }),
          [...KernelRenderers],
          has ? [entry] : [],
        ),
    };
    const checks = await rendererConformance(entry, long);
    expect(failures(checks)).toEqual([]);
    expect(checks.some((c) => c.check === "batches" && / at 0\.\d+ px/.test(c.detail))).toBe(true);
  });

  it("a renderer that draws nothing fails `batches` at every zoom (VB-04)", async () => {
    const blank = { ...StarsRenderer, draw: () => {} } as LayerRenderer;
    const batches = (await rendererConformance(blank, fixtures(blank))).filter((c) => c.check === "batches");
    // Each zoom, and the layer over its zooms in each geometry.
    expect(batches.filter((c) => !c.ok && c.detail.includes(" 0 items")).length).toBe(5);
    expect(batches.filter((c) => !c.ok && c.detail.endsWith("no item drawn at any zoom")).length).toBe(3);
    expect(batches.length).toBe(8);
  });

  it("an entry declaring only linear falls back in wrapped, with svs.render.fallback once", async () => {
    const linear = { ...StarsRenderer, geometries: ["linear"] } as LayerRenderer;
    const checks = await rendererConformance(linear, fixtures(linear));
    expect(failures(checks)).toEqual([]);
    expect(checks.find((c) => c.check === "fallback" && c.geometry === "wrapped")?.detail).toBe(
      "1 layer(s), 1 issue(s) [svs.render.fallback]",
    );
  });

  it("sizes by the field provider and reads custom in prepare; hits name the item drawn", async () => {
    const ctx = await open(starsDoc(), [...KernelRenderers], [StarsRenderer as LayerRenderer]);
    const m = ctx.graph.assembled.value?.view;
    const k = m?.prepared.findIndex((p) => p.renderer === "org.example:stars") ?? -1;
    const s = m?.prepared[k].state as { radius: Float64Array; unit: string };
    // pKd 9, 8, 7, 6 → radius pKd − 3, clamped to [3, 10].
    expect([...s.radius, s.unit]).toEqual([6, 5, 4, 3, "pKd"]);
    const hit = Hit.create();
    const ax = ctx.layout.axis("v.main");
    const y = (() => {
      const l = ctx.layout.get("v.main")?.current;
      const e = l ? Array.from(l.entryLayer.subarray(0, l.entryCount)).indexOf(k) : -1;
      return l ? l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2 : 0;
    })();
    const k12 = ax?.shownOf(0, 12) ?? (0 as ShownIndex);
    const x = ax ? (ax.x(0, k12) + ax.x(0, (k12 + 1) as never)) / 2 : 0;
    ctx.layout.hitTest("v.main", x, y, hit);
    expect([hit.kind, hit.item, hit.layer]).toEqual([HoverKind.item, 1, k]);
    ctx.dispose();
  });
});

describe("field providers", () => {
  it("a `ranges` layer colors by org.example:affinity through `field_palette`, read once per annotation", async () => {
    const read = vi.spyOn(AffinityField, "read");
    const ctx = await open(starsDoc({ ranges: true }), [...KernelRenderers], [StarsRenderer as LayerRenderer]);
    const m = ctx.graph.assembled.value?.view;
    const k = m?.prepared.findIndex((p) => p.renderer === "ranges") ?? -1;
    const theme = m?.prepared[k].themes[0];
    const loc = { annotation: m?.resolved.data[m.resolved.layers.dataStart[k]] ?? -1, row: 0 };
    const at = (row: number) => theme?.value({ ...loc, row, item: row, value: Number.NaN } as never);
    // pKd 9 is the domain's top (white), 6 its bottom (black).
    expect([at(0), at(3)]).toEqual([0xffffff, 0x000000]);
    expect(theme?.granularity).toBeDefined();
    expect(read).toHaveBeenCalledTimes(1);
    const set = ctx.graph.assembled.value?.set;
    const a = set?.annotations[loc.annotation];
    expect(a && set ? Array.from(numberColumn(a, set, "org.example:affinity", ctx.registry.svs) ?? []) : []).toEqual([
      9, 8, 7, 6,
    ]);
    read.mockRestore();
    ctx.dispose();
  });

  it("a `rules` measure reads a provider's field: pKd 8 and over in red", async () => {
    const ctx = await open(starsDoc({ rules: true }), [...KernelRenderers], [StarsRenderer as LayerRenderer]);
    const m = ctx.graph.assembled.value?.view;
    const k = m?.prepared.findIndex((p) => p.renderer === "points") ?? -1;
    const theme = m?.prepared[k].themes[0];
    const annotation = m?.resolved.data[m.resolved.layers.dataStart[k]] ?? -1;
    const at = (row: number) => theme?.value({ annotation, row, item: row, value: Number.NaN } as never);
    expect([0, 1, 2, 3].map(at)).toEqual([0xff0000, 0xff0000, 0x808080, 0x808080]);
    ctx.dispose();
  });
});

describe("replacing a built-in renderer", () => {
  it("a host's own `ranges` on case 13 draws, hits, takes its themes and lists its legend; the built-in is unused", async () => {
    const spies = (["prepare", "hit", "legend", "height", "lineHeights"] as const).flatMap((m) =>
      RangesRenderer[m] ? [vi.spyOn(RangesRenderer as Required<LayerRenderer<Features>>, m)] : [],
    );
    const rects = new RectBatch();
    const legend: string[] = [];
    const drawn = new Set<number>();
    // The same type semantics (the built-in prepare, from the API), another look: one flat bar per item, no lanes.
    const host: LayerRenderer<Features> = {
      name: "ranges",
      geometries: ["linear"],
      stable: true,
      minHeight: 4,
      content: false,
      prepare: (input, rc) => PreparedParts.ranges.prepare(input, rc) as Features,
      height: () => 12,
      draw(v: DrawView, f: Features) {
        rects.reset();
        for (let i = 0; v.axis && i < f.count; i++) {
          const loc = v.loc;
          loc.annotation = f.annotation[i];
          loc.row = f.row[i];
          const color = v.color ? v.color.value(loc) : 0;
          drawn.add(color);
          rects.begin(color);
          band(v, f.block[i], f.start[i], Math.max(f.stop[i], f.start[i] + 1), v.top + 2, v.top + 10, rects);
        }
        v.backend.rects(rects);
      },
      hit(v, f, _p, out) {
        for (let i = f.count - 1; i >= 0; i--)
          if (f.block[i] === out.block && f.start[i] <= out.index && out.index < Math.max(f.stop[i], f.start[i] + 1))
            return hitItem(out, (v.axis as NonNullable<typeof v.axis>).set, f.annotation[i], f.row[i]);
        return false;
      },
      legend: (_p, themes, out) => {
        for (const t of themes) t?.legend?.(out);
        legend.push("host");
      },
    };
    const checks = await rendererConformance(host as LayerRenderer, {
      // At 16 px a position from the start, some of case 13's layers have no item in view.
      empty: ["middle", "close"],
      open: async ({ geometry, entry }) => {
        const doc = geometry === "wrapped" ? asText(case13()) : geometry === "circular" ? ringed(case13()) : case13();
        const renderers = [...without("ranges"), ...(entry ? [host as LayerRenderer] : [])];
        return (await opened(doc, { renderers, registry: [RINGS] })).ctx;
      },
    });
    // Without any `ranges`, case 13's ranges layers have no fallback: omitted, one issue each.
    expect(failures(checks)).toEqual([]);
    const frames = checks.filter((c) => c.check === "batches");
    expect(frames.length).toBeGreaterThan(0);
    expect(frames.every((c) => /: [1-9]\d* entries/.test(c.detail))).toBe(true);
    // The open groups' ranges layers (Chain, Domain, Region, Motif) drew in their uniform colors, through the host.
    expect(drawn).toEqual(new Set([0xcc9933, 0x7a8ff0, 0xd17b3a, 0x3aa58f]));
    expect(legend.length).toBeGreaterThan(0);
    const ctx: SeqViewerContext = (await opened(case13(), { renderers: [...without("ranges"), host as LayerRenderer] }))
      .ctx;
    expect(ctx.registry.get("renderers", "ranges")).toBe(host);
    for (const s of spies) expect(s).not.toHaveBeenCalled();
    for (const s of spies) s.mockRestore();
    ctx.dispose();
  });
});
