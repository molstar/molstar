/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// KS1: keyed sets (seqspec:format#layer-types) and the kernel providers `row_stats` and `redundancy`
// (seqspec:format#computed-nodes): a rows layer's `order` and `filter` over key tables, the permutation and kept mask
// resolved into the lanes, and the providers' key tables.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Issue, SeqSet } from "@molstar/seq-core/model";
import { beforeEach, describe, expect, it } from "vitest";
import { compile, resolveView, type SeqViewSpecDoc, supportedAll } from "../../src/format/index.js";
import { clearNodeProviderCache } from "../../src/format/load/node-providers.js";
import { fixtureLoaders, registries } from "./load-fixtures.js";
import { scripted } from "./stream-provider.js";

// Row 0 is the reference. Row 1 differs at one of ten, 2 at two, 3 is row 1 with one more change, 4 is far, 5 is row 4 with
// a gap run (so only 6 of its columns have a letter).
const ROWS = ["ACDEFGHIKL", "ACDEFGHIKM", "ACDEFGHIAM", "ACDEFGHIKV", "WYWYWYWYWY", "WYWYWY----"];
const reg = registries();
const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const errors = (issues: Issue[]) => issues.filter((i) => i.severity === "error");

interface Opts {
  rows?: Record<string, unknown>;
  stats?: Record<string, unknown> | null;
  nr?: Record<string, unknown> | null;
  mid?: (b: ReturnType<typeof createSeqViewSpec>) => void;
}

function doc(o: Opts = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "Keyed" });
  const aln = b.data.alignment({ ref: "aln", value: ROWS, labels: ROWS.map((_, i) => `r${i}`) });
  if (o.stats !== null)
    b.data.computed({ ref: "stats", provider: "row_stats", alignment: aln, ...(o.stats ?? {}) } as never);
  if (o.nr !== null)
    b.data.computed({ ref: "nr", provider: "redundancy", alignment: aln, threshold: 0.8, ...(o.nr ?? {}) } as never);
  o.mid?.(b);
  const v = b.view({ ref: "v", label: "V" }).block({ alignment: aln });
  v.track({ label: "Rows" }).rows((o.rows ?? {}) as never);
  return b.getDocument();
}
const load = async (d: SeqViewSpecDoc, signal?: AbortSignal) => compile(d, reg, fixtureLoaders(), signal);
async function lanes(d: SeqViewSpecDoc) {
  const { set, issues } = await load(d);
  const r = resolveView(d, set, "v", reg, supportedAll(reg));
  return {
    set,
    issues: [...issues, ...r.issues],
    rows: Array.from(r.view.lanes.row),
    depth: Array.from(r.view.lanes.depth),
  };
}
const table = (set: SeqSet, key: string) => {
  const t = { kind: -1, id: -1, row: -1 };
  const ann = set.byName(key, t) ? set.annotations[t.id] : undefined;
  const col = (name: string) => {
    const id = set.strings.find(name);
    const f = ann?.fields.find((x) => x.name === id);
    if (f?.kind === "string") return Array.from(f.data, (c) => (c < 0 ? null : set.strings.get(c)));
    return f ? Array.from(f.data as ArrayLike<unknown>) : undefined;
  };
  return { ann, col, rowOf: ann ? Array.from(ann.target.row) : [] };
};

beforeEach(clearNodeProviderCache);

describe("row_stats", () => {
  it("is a rows table of identity to the reference row, coverage and length", async () => {
    const { set, issues } = await load(doc({ nr: null }));
    expect(errors(issues)).toEqual([]);
    const t = table(set, "stats");
    expect(t.ann?.geometry).toBe("rows");
    expect(t.rowOf).toEqual([0, 1, 2, 3, 4, 5]);
    expect(t.col("identity")).toEqual([1, 0.9, 0.8, 0.9, 0, 0]);
    expect(t.col("coverage")).toEqual([1, 1, 1, 1, 1, 0.6]);
    expect(t.col("length")).toEqual([10, 10, 10, 10, 10, 6]);
  });

  it("takes a reference row, and reports one that is not there", async () => {
    const r = await load(doc({ nr: null, stats: { reference_row: 4 } }));
    expect(table(r.set, "stats").col("identity")?.[4]).toBe(1);
    const bad = await load(doc({ nr: null, stats: { reference_row: 9 } }));
    expect(codes(bad.issues)).toContain("svs.ref.derived_unknown");
  });

  it("counts an annotation's items on each row, of a category and within columns, with the largest score", async () => {
    const mid = (b: ReturnType<typeof createSeqViewSpec>) => {
      b.data.annotation({
        ref: "hits",
        geometry: "ranges",
        target: "aln",
        value: {
          target: ["aln/1", "aln/1", "aln/2", "aln/5"],
          start: [0, 7, 3, 1],
          stop: [2, 9, 5, 3],
          category: ["a", "b", "a", "a"],
          score: [0.5, 2, 1, 9],
        },
      } as never);
    };
    const stats = { items: "hits", score: "score" };
    const all = await load(doc({ nr: null, stats, mid }));
    expect(errors(all.issues)).toEqual([]);
    expect(table(all.set, "stats").col("items")).toEqual([0, 2, 1, 0, 0, 1]);
    expect(
      table(all.set, "stats")
        .col("score")
        ?.map((x) => (Number.isNaN(x) ? null : x)),
    ).toEqual([null, 2, 1, null, null, 9]);
    const cat = await load(doc({ nr: null, stats: { ...stats, category: "a" }, mid }));
    expect(table(cat.set, "stats").col("items")).toEqual([0, 1, 1, 0, 0, 1]);
    const cols = await load(doc({ nr: null, stats: { ...stats, columns: [6, 10] }, mid }));
    expect(table(cols.set, "stats").col("items")).toEqual([0, 1, 0, 0, 0, 0]);
  });

  it("counts the items of an annotation on a sequence on the rows its mappings bind it to, through the path", async () => {
    const mid = (b: ReturnType<typeof createSeqViewSpec>) => {
      const seq = b.data.sequence({ ref: "s", value: "ACDEFGHIKL" });
      b.data.mapping({ ref: "m0", from: "aln/0", to: seq, method: "identity" } as never);
      b.data.mapping({ ref: "m2", from: "aln/2", to: seq, method: "identity" } as never);
      b.data.annotation({
        ref: "site",
        geometry: "ranges",
        target: seq,
        value: { start: [2, 8], stop: [4, 9], score: [3, 5] },
      } as never);
    };
    const { set, issues } = await load(doc({ nr: null, stats: { items: "site", score: "score" }, mid }));
    expect(errors(issues)).toEqual([]);
    expect(table(set, "stats").col("items")).toEqual([2, 0, 2, 0, 0, 0]);
    const some = await load(doc({ nr: null, stats: { items: "site", score: "score", columns: [0, 5] }, mid }));
    expect(table(some.set, "stats").col("items")).toEqual([1, 0, 1, 0, 0, 0]);
    expect(
      table(some.set, "stats")
        .col("score")
        ?.map((x) => (Number.isNaN(x) ? null : x)),
    ).toEqual([3, null, 3, null, null, null]);
  });

  it("is svs.provider.lazy over a streamed alignment's rows (a warning), and an error for an alignment that is not there", async () => {
    const manifest = {
      systems: [{ key: "aln", extent: [0, 8] as [number, number], numbering: null, label: null }],
      levels: [[1, 1] as [number, number]],
      rows: 50,
      columns: 8,
    };
    const s = scripted(manifest, (r) => ({ ...r, kind: "rows", labels: [`row${r.row_start}`], rows: ["ACDEFGHI"] }));
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: "msa", stream: "msa", data_kind: "alignment" });
    const aln = b.data.alignment({ ref: "aln", source_ref: src });
    b.data.computed({ ref: "stats", provider: "row_stats", alignment: aln } as never);
    b.view({ ref: "v", label: "V" }).block({ alignment: "aln" }).track({ label: "R" }).ruler();
    const lazy = await compile(b.getDocument(), reg, {
      load: (r) => Promise.reject(new Error(`no fixture ${r.uri}`)),
      stream: () => s.provider,
    });
    expect(lazy.issues.map((i) => [i.code, i.severity])).toEqual([["svs.provider.lazy", "warning"]]);
    expect(table(lazy.set, "stats").ann).toBeUndefined();
    const gone = await load(doc({ nr: null, stats: { alignment: "nope" } }));
    expect(gone.issues.some((i) => i.severity === "error")).toBe(true);
  });
});

describe("redundancy", () => {
  it("clusters greedily in data order: a row joins the first representative it is at least threshold identical to", async () => {
    const { set, issues } = await load(doc({ stats: null }));
    expect(errors(issues)).toEqual([]);
    const t = table(set, "nr");
    // Row 1 (0.9) and 2 (0.8) join row 0; row 3 (0.9) too; rows 4 and 5 are far, 5 shares 6 columns with 4 (all equal).
    expect(t.col("cluster")).toEqual([0, 0, 0, 0, 4, 4]);
    expect(t.col("representative")).toEqual([true, false, false, false, true, false]);
  });

  it("a higher threshold keeps more representatives, and the order decides who represents", async () => {
    const hi = await load(doc({ stats: null, nr: { threshold: 0.9 } }));
    expect(table(hi.set, "nr").col("cluster")).toEqual([0, 0, 2, 0, 4, 4]);
    const byId = await load(
      doc({ nr: { threshold: 0.9, order: { data: "stats", field: "identity", direction: "asc" } } }),
    );
    // Taken from the least identical to row 0: 4, 5 (joins 4), 2, 1 (0.9 to 2), 3 (0.8 to 2), 0 (0.9 to 3).
    expect(table(byId.set, "nr").col("cluster")).toEqual([3, 2, 2, 3, 4, 4]);
  });

  it("is unavailable past Redundancy.MaxCells (svs.provider.size, a warning); the rest of the document loads", async () => {
    // 1,300 rows x 1,300 rows x 1,300 columns is over 2^31 cells; no work is done before the check.
    const wide = Array.from({ length: 1300 }, (_, r) => (r % 2 ? "A" : "C").repeat(1300));
    const b = createSeqViewSpec({ title: "Big" });
    const aln = b.data.alignment({ ref: "aln", value: wide });
    b.data.computed({ ref: "nr", provider: "redundancy", alignment: aln, threshold: 0.9 } as never);
    b.data.computed({ ref: "stats", provider: "row_stats", alignment: aln } as never);
    b.view({ ref: "v", label: "V" }).block({ alignment: aln });
    const r = await load(b.getDocument());
    expect(r.issues.filter((i) => i.code === "svs.provider.size").map((i) => i.severity)).toEqual(["warning"]);
    expect(errors(r.issues)).toEqual([]);
    expect(table(r.set, "nr").ann).toBeUndefined();
    expect(table(r.set, "stats").ann?.count).toBe(1300);
  });

  it("stops at the stage's abort signal", async () => {
    const ac = new AbortController();
    ac.abort();
    await expect(load(doc(), ac.signal)).rejects.toThrow();
  });

  it("caches by the rows' content and its options", async () => {
    const a = await load(doc({ stats: null }));
    const b = await load(doc({ stats: null }));
    expect(table(b.set, "nr").col("cluster")).toEqual(table(a.set, "nr").col("cluster"));
  });
});

describe("a rows layer's order and filter", () => {
  it("default: tree order when there is a tree, else data order", async () => {
    const r = await lanes(doc({ stats: null, nr: null }));
    expect(r.rows).toEqual([0, 1, 2, 3, 4, 5]);
  });

  it("orders by a key table's field, descending, missing members last, ties in data order", async () => {
    const r = await lanes(doc({ nr: null, rows: { order: { data: "stats", field: "identity", direction: "desc" } } }));
    expect(errors(r.issues)).toEqual([]);
    expect(r.rows).toEqual([0, 1, 3, 2, 4, 5]);
    const asc = await lanes(doc({ nr: null, rows: { order: { data: "stats", field: "identity" } } }));
    expect(asc.rows).toEqual([4, 5, 2, 1, 3, 0]);
  });

  it("members the table has no row for come last, whichever the direction", async () => {
    const mid = (b: ReturnType<typeof createSeqViewSpec>) => {
      b.data.annotation({
        ref: "pick",
        geometry: "rows",
        target: "aln",
        value: { row: [3, 1], rank: [2, 7] },
      } as never);
    };
    const asc = await lanes(doc({ stats: null, nr: null, mid, rows: { order: { data: "pick", field: "rank" } } }));
    expect(asc.rows).toEqual([3, 1, 0, 2, 4, 5]);
    const desc = await lanes(
      doc({ stats: null, nr: null, mid, rows: { order: { data: "pick", field: "rank", direction: "desc" } } }),
    );
    expect(desc.rows).toEqual([1, 3, 0, 2, 4, 5]);
    // A filter keeps only members the table has a value for.
    const kept = await lanes(
      doc({ stats: null, nr: null, mid, rows: { filter: [{ data: "pick", field: "rank", min: 0 }] } }),
    );
    expect(kept.rows).toEqual([1, 3]);
  });

  it("clusters, then by identity within a cluster (a list is lexicographic)", async () => {
    const order = [
      { data: "nr", field: "cluster" },
      { data: "stats", field: "identity", direction: "desc" },
    ];
    const r = await lanes(doc({ rows: { order } }));
    expect(errors(r.issues)).toEqual([]);
    expect(r.rows).toEqual([0, 1, 3, 2, 4, 5]);
    const hi = await lanes(doc({ nr: { threshold: 0.9 }, rows: { order } }));
    // clusters 0 (rows 0, 1, 3), 2 (row 2), 4 (rows 4, 5).
    expect(hi.rows).toEqual([0, 1, 3, 2, 4, 5]);
    const by2 = await lanes(doc({ nr: { threshold: 0.95 }, rows: { order } }));
    expect(by2.rows).toEqual([0, 1, 2, 3, 4, 5].sort((x, y) => x - y));
  });

  it("filters to the representatives; row_subset applies first and every filter must keep a row", async () => {
    const nr = { data: "nr", field: "representative", values: [true] };
    const r = await lanes(doc({ rows: { filter: [nr] } }));
    expect(r.rows).toEqual([0, 4]);
    const sub = await lanes(doc({ rows: { filter: [nr], row_subset: [1, 2, 4, 5] } }));
    expect(sub.rows).toEqual([4]);
    const both = await lanes(doc({ rows: { filter: [nr, { data: "stats", field: "identity", max: 0.5 }] } }));
    expect(both.rows).toEqual([4]);
    const range = await lanes(
      doc({ nr: null, rows: { filter: [{ data: "stats", field: "identity", min: 0.85, max: 0.95 }] } }),
    );
    expect(range.rows).toEqual([1, 3]);
  });

  it("order and filter together: the filter keeps, the order sorts what is kept", async () => {
    const r = await lanes(
      doc({
        rows: {
          order: { data: "stats", field: "identity", direction: "desc" },
          filter: [{ data: "stats", field: "coverage", min: 1 }],
        },
      }),
    );
    expect(r.rows).toEqual([0, 1, 3, 2, 4]);
  });

  it("a table, a field or a member list that does not fit is svs.layer.data and the key is skipped", async () => {
    const nope = await lanes(doc({ rows: { order: { data: "nothing", field: "x" } } }));
    expect(codes(nope.issues).filter((c) => c.startsWith("svs.layer"))).toEqual([]);
    expect(nope.issues.some((i) => i.code === "svs.ref.unresolved")).toBe(true);
    const field = await lanes(doc({ rows: { order: { data: "stats", field: "idnetity" } } }));
    const i = field.issues.find((x) => x.code === "svs.layer.data");
    expect(i?.message).toContain('no field "idnetity"');
    expect(i?.message).toContain("identity");
    expect(field.rows).toEqual([0, 1, 2, 3, 4, 5]);
  });

  it("key tables keyed by state have no target and read as string fields", async () => {
    const d = createSeqViewSpec({ title: "States" });
    const aln = d.data.alignment({ ref: "aln", value: ROWS });
    d.data.annotation({
      ref: "samples",
      geometry: "rows",
      value: { state: ["A", "B", "C"], population: ["x", "y", "x"], age: [30, 20, 10] },
    } as never);
    d.view({ ref: "v", label: "V" }).block({ alignment: aln });
    const { set, issues } = await load(d.getDocument());
    expect(errors(issues)).toEqual([]);
    const t = table(set, "samples");
    expect(t.col("state")).toEqual(["A", "B", "C"]);
    expect(t.rowOf).toEqual([-1, -1, -1]);
    expect(t.col("age")).toEqual([30, 20, 10]);
  });
});

describe("a template's order and filter (for_each)", () => {
  it("orders and filters instances by a key table keyed by node; key order stays the default", async () => {
    const { SeqSets } = await import("@molstar/seq-core/model");
    const { indexDoc, templateNodes } = await import("../../src/format/index.js");
    const P = "MKVLAAGGHHKK";
    const ann = (key: string) => ({
      key,
      geometry: "ranges" as const,
      target: "p",
      columns: { start: [1], stop: [5] },
    });
    const set = SeqSets.build({
      molecules: [{ key: "p", residues: P, ownUnit: true }],
      annotations: [
        ann("genes/a"),
        ann("genes/b"),
        ann("genes/c"),
        {
          key: "info",
          geometry: "rows" as const,
          columns: { node: ["genes/c", "genes/a", "genes/b"], size: [30, 20, 10], big: [true, false, false] },
        },
      ],
    }).set;
    const track = (fe: Record<string, unknown>) => {
      const b = createSeqViewSpec({ title: "t" });
      const p = b.data.sequence({ ref: "p", value: P });
      b.data.source({ ref: "genes", stream: "genes", data_kind: "items" });
      const v = b.view({ ref: "v", label: "v" });
      v.block({ sequence: p });
      v.track({ ref: "tpl", for_each: { kind: "annotation", source: "genes", ...fe } as never }).ranges({
        data: "$node",
      });
      return indexDoc(b.getDocument(), reg);
    };
    const pick = (fe: Record<string, unknown>) => templateNodes("tpl", track(fe), set, undefined, reg, () => {});
    expect(pick({})).toEqual(["genes/a", "genes/b", "genes/c"]);
    expect(pick({ order: { data: "info", field: "size", direction: "desc" } })).toEqual([
      "genes/c",
      "genes/a",
      "genes/b",
    ]);
    expect(pick({ order: { data: "info", field: "size" } })).toEqual(["genes/b", "genes/a", "genes/c"]);
    expect(pick({ filter: [{ data: "info", field: "big", values: [false] }] })).toEqual(["genes/a", "genes/b"]);
    expect(pick({ order: "data", filter: [{ data: "info", field: "size", min: 20 }] })).toEqual(["genes/a", "genes/c"]);
  });
});

describe("the committed PF00042 (the story's data)", () => {
  it("32 globins: HBA_HUMAN leads by identity; 30 representatives at 90%, 22 at 50%", async () => {
    const { pf00042 } = await import("../model/alignment-fixtures.js");
    const family = pf00042();
    const by = { data: "stats", field: "identity", direction: "desc" };
    const cluster = (ref: string, threshold: number) => (b: ReturnType<typeof createSeqViewSpec>) => {
      b.data.computed({ ref, provider: "redundancy", alignment: "aln", threshold, order: by } as never);
    };
    const b = createSeqViewSpec({ title: "PF00042" });
    const aln = b.data.alignment({
      ref: "aln",
      value: family.map((r) => r.gapped),
      labels: family.map((r) => r.label),
    });
    b.data.computed({ ref: "stats", provider: "row_stats", alignment: aln } as never);
    cluster("nr90", 0.9)(b);
    cluster("nr50", 0.5)(b);
    const v = b.view({ ref: "v", label: "V" }).block({ alignment: aln });
    v.track({ label: "Rows" }).rows({
      order: by,
      filter: [{ data: "nr90", field: "representative", values: [true] }],
    } as never);
    const r = await lanes(b.getDocument());
    expect(errors(r.issues)).toEqual([]);
    expect(r.rows).toHaveLength(30);
    expect(r.rows[0]).toBe(0);
    const ids = table(r.set, "stats").col("identity") as number[];
    expect(r.rows.map((x) => ids[x])).toEqual([...r.rows.map((x) => ids[x])].sort((x, y) => y - x));
    expect((table(r.set, "nr50").col("representative") as boolean[]).filter(Boolean)).toHaveLength(22);
  });
});

// The lazy-alignment path (the Opus review of KS1): a streamed alignment is ordered by tree or data only.
describe("a streamed alignment and key orders", () => {
  const manifest = {
    systems: [{ key: "aln", extent: [0, 8] as [number, number], numbering: null, label: null }],
    levels: [[1, 1] as [number, number]],
    rows: 6,
    columns: 8,
  };
  async function lazyView(
    layer: Record<string, unknown>,
    tree?: string,
    more?: (b: ReturnType<typeof createSeqViewSpec>) => void,
  ) {
    const { indexDoc } = await import("../../src/format/index.js");
    const { resolveSection, resolveLayer } = await import("../../src/format/index.js");
    const s = scripted(manifest, (r) => ({ ...r, kind: "rows", labels: [`row${r.row_start}`], rows: ["ACDEFGHI"] }));
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: "msa", stream: "msa", data_kind: "alignment" });
    const aln = b.data.alignment({ ref: "aln", source_ref: src });
    if (tree) aln.tree({ value: tree, leaves: "order" } as never);
    b.data.annotation({
      ref: "A",
      geometry: "rows",
      target: "aln",
      value: { row: [0, 1, 2, 3, 4, 5], k: [5, 4, 3, 2, 1, 0] },
    } as never);
    more?.(b);
    const v = b.view({ ref: "v", label: "V" }).block({ alignment: aln });
    const t = v.track({ label: "R" });
    t.rows(layer as never);
    t.ranges({ ref: "L", data: "site" } as never);
    const d = b.getDocument();
    const c = await compile(d, reg, { load: () => Promise.reject(new Error("x")), stream: () => s.provider });
    const idx = indexDoc(d, reg);
    const sec = [...idx.nodes.values()].find((n) => n.kind === "section")?.ref as string;
    const r = resolveSection(idx, c.set, sec, reg);
    return { c, idx, set: c.set, section: r.section, issues: [...c.issues, ...r.issues], resolveLayer };
  }
  const SITE = (b: ReturnType<typeof createSeqViewSpec>) =>
    b.data.annotation({
      ref: "site",
      geometry: "ranges",
      target: "aln",
      value: { target: [0, 1, 2, 3, 4, 5].map((k) => `aln/${k}`), start: [0, 0, 0, 0, 0, 0], stop: [1, 1, 1, 1, 1, 1] },
    } as never);

  it("leaves a key order and a filter out, in tree or data order, and says why", async () => {
    const r = await lazyView(
      { order: { data: "A", field: "k" }, filter: [{ data: "A", field: "k", min: 3 }] },
      undefined,
      SITE,
    );
    const lazy = r.issues.filter((i) => i.code === "svs.provider.lazy");
    expect(lazy.map((i) => i.severity)).toEqual(["warning", "warning"]);
    const said = lazy.map((i) => i.message);
    expect(said).toHaveLength(2);
    expect(said.every((m) => m.includes("a streamed alignment is ordered by tree or data only"))).toBe(true);
    const run = r.section.lanes.lazy[0];
    expect([run.length, run.order]).toEqual([6, null]);
    // plain tree and data orders are not reported
    const plain = await lazyView({ order: "data" }, undefined, SITE);
    expect(plain.issues.filter((i) => i.code === "svs.provider.lazy")).toEqual([]);
  });

  it("a run's order is in its axis key: another permutation is another axis", async () => {
    const { orderHash } = await import("../../src/format/resolve/section.js");
    const r = await lazyView({}, "((row0,row1),((row2,row3),(row4,row5)));", SITE);
    const run = r.section.lanes.lazy[0];
    expect(run.order).not.toBeNull();
    expect(r.section.axisKey).toContain(`${run.start}+6t${orderHash(run.order as Int32Array)}`);
    const h = (o: number[]) => orderHash(Int32Array.from(o));
    expect(h([0, 1, 2, 3, 4, 5])).toBe(h([0, 1, 2, 3, 4, 5]));
    expect(h([0, 1, 2, 3, 4, 5])).not.toBe(h([1, 0, 2, 3, 4, 5]));
    expect(h([0, 1, 2, 3, 4, 5])).not.toBe(h([0, 1, 2, 3, 5, 4]));
    // a data-order run has none
    const data = await lazyView({ order: "data" }, undefined, SITE);
    expect(data.section.axisKey).toContain("+6d");
  });

  it("a run with an order places only the rows it holds: the others are left out, not at lane -1", async () => {
    const r = await lazyView({}, "((row0,row1),((row2,row3),(row4,row5)));", SITE);
    const run = r.section.lanes.lazy[0];
    const cut = {
      ...r.section,
      lanes: { ...r.section.lanes, lazy: [{ ...run, length: 3, order: Int32Array.of(3, 4, 5) }] },
    };
    const { supportedAll } = await import("../../src/format/index.js");
    const l = r.resolveLayer(cut as never, "L", r.idx, r.set, reg, supportedAll(reg)).layer;
    expect(l.placements.map((p) => p.lane - run.start).sort()).toEqual([0, 1, 2]);
    expect(l.placements.every((p) => p.lane >= run.start)).toBe(true);
  });
});

describe("keyed details", () => {
  it("a field the table has but no member has a value for is 'no values', and a filter over it keeps none", async () => {
    const mid = (b: ReturnType<typeof createSeqViewSpec>) => {
      b.data.annotation({
        ref: "empty",
        geometry: "rows",
        target: "aln",
        value: { row: [0, 1, 2], k: [null, null, null], j: [1, 2, 3] },
      } as never);
    };
    const r = await lanes(
      doc({ stats: null, nr: null, mid, rows: { filter: [{ data: "empty", field: "k", min: 0 }] } }),
    );
    expect(r.rows).toEqual([]);
    const i = r.issues.find((x) => x.code === "svs.layer.data");
    expect(i?.message).toContain('has no values for "k"');
    // an unknown field is still "no field", with the fields it has, and the filter is skipped
    const bad = await lanes(
      doc({ stats: null, nr: null, mid, rows: { filter: [{ data: "empty", field: "zzz", min: 0 }] } }),
    );
    expect(bad.rows).toEqual([0, 1, 2, 3, 4, 5]);
    expect(bad.issues.find((x) => x.code === "svs.layer.data")?.message).toContain('has no field "zzz"');
  });

  it("rows with a parent are indented under it in tree or data order, and not under a key order", async () => {
    const make = (rows: Record<string, unknown>) => {
      const b = createSeqViewSpec({ title: "parents" });
      const aln = b.data.alignment({ ref: "aln", value: ROWS, parent_rows: [null, 0, null, 2, null, null] } as never);
      b.data.annotation({
        ref: "T",
        geometry: "rows",
        target: "aln",
        value: { row: [0, 1, 2, 3, 4, 5], k: [5, 4, 3, 2, 1, 0] },
      } as never);
      b.view({ ref: "v", label: "V" })
        .block({ alignment: aln })
        .track({ label: "R" })
        .rows(rows as never);
      return b.getDocument();
    };
    expect((await lanes(make({ order: "data" }))).depth).toEqual([0, 1, 0, 1, 0, 0]);
    expect((await lanes(make({}))).depth).toEqual([0, 1, 0, 1, 0, 0]);
    const keyedOrder = await lanes(make({ order: { data: "T", field: "k" } }));
    expect(keyedOrder.rows).toEqual([5, 4, 3, 2, 1, 0]);
    expect(keyedOrder.depth).toEqual([0, 0, 0, 0, 0, 0]);
  });
});

describe("a matrix's states (KS2)", () => {
  /** Five samples as states over a short sequence, and a sample table keyed by state (a state it does not name: none). */
  const SAMPLES = {
    state: ["n0", "n1", "n2", "n3", "nX"],
    population: ["x", "y", "x", "x", "y"],
    age: [30, 20, 10, 40, 5],
  };
  function matrixDoc(
    layer: Record<string, unknown>,
    states = ["n0", "n1", "n2", "n3", "n4"],
    samples: object = SAMPLES,
  ) {
    const b = createSeqViewSpec({ title: "Genotypes" });
    const s = b.data.sequence({ ref: "s", value: "ACGTACGTAC" });
    // A table keyed by row: an alignment's, naming no state.
    const aln = b.data.alignment({ ref: "aln", value: ["ACGT", "ACGA"] });
    b.data.annotation({ ref: "byRow", geometry: "rows", target: aln, value: { row: [0, 1], k: [1, 2] } } as never);
    const index = Array.from({ length: 10 }, (_, i) => i);
    const value = { index, value: index.map((i) => states.map((_, q) => (i + q) % 3)) };
    const m = b.data.annotation({ ref: "m", geometry: "matrix", target: s, states, value } as never);
    b.data.annotation({
      ref: "samples",
      geometry: "rows",
      value: samples,
    } as never);
    b.view({ ref: "v", label: "V" })
      .block({ sequence: s })
      .track({ label: "G" })
      .matrix({ data: m, ...layer } as never);
    return b.getDocument();
  }
  async function members(layer: Record<string, unknown>, states?: string[], samples?: object) {
    const d = matrixDoc(layer, states, samples);
    const { set, issues } = await load(d);
    const r = resolveView(d, set, "v", reg, supportedAll(reg));
    const at = r.view.layers.members[0];
    return { members: at ? Array.from(at) : null, issues: [...issues, ...r.issues] };
  }

  it("orders by a key table's field, missing values last, and filters; data order by default", async () => {
    expect((await members({})).members).toBeNull();
    expect((await members({ order: "data" })).members).toBeNull();
    const byAge = await members({ order: { data: "samples", field: "age", direction: "desc" } });
    expect(errors(byAge.issues)).toEqual([]);
    expect(byAge.members).toEqual([3, 0, 1, 2, 4]);
    const kept = await members({
      order: [
        { data: "samples", field: "population" },
        { data: "samples", field: "age" },
      ],
      filter: [{ data: "samples", field: "age", min: 15 }],
    });
    expect(kept.members).toEqual([0, 3, 1]);
    const only = await members({ filter: [{ data: "samples", field: "population", values: ["y"] }] });
    expect(only.members).toEqual([1]);
  });

  it("tree order is reported and is data order; a filter that keeps nothing is noted; a row-keyed table names none", async () => {
    const tree = await members({ order: [{ data: "samples", field: "age" }, "tree"] });
    const said = tree.issues.find((x) => x.code === "svs.layer.data");
    expect(said?.path).toContain("params/order/1");
    expect(said?.message).toContain("no tree");
    expect(tree.members).toEqual([2, 1, 0, 3, 4]);
    const alone = await members({ order: "tree" });
    expect(alone.members).toBeNull();
    expect(alone.issues.some((x) => x.code === "svs.layer.data" && x.message.includes("no tree"))).toBe(true);
    const none = await members({ filter: [{ data: "samples", field: "age", min: 100 }] });
    expect(none.members).toEqual([]);
    // Valid input read as written: a note, not an error (GV-10).
    const empty = none.issues.find((x) => x.code === "svs.layer.empty");
    expect([empty?.severity, empty?.message]).toEqual(["info", expect.stringContaining("keeps no state")]);
    const byRow = await members({ order: { data: "byRow", field: "k" } });
    expect(byRow.issues.find((x) => x.code === "svs.layer.data")?.message).toContain("no row for any member");
    expect(byRow.members).toEqual([0, 1, 2, 3, 4]);
  });

  it('a key table\'s row with no state names no state, not one named "undefined" (GV-11)', async () => {
    const samples = { state: ["n0", "n1", "n2", "n3", null], age: [30, 20, 10, 40, 5] };
    const r = await members(
      { order: { data: "samples", field: "age" } },
      ["n0", "n1", "n2", "n3", "undefined"],
      samples,
    );
    // The state "undefined" has no row: last, not first with the stateless row's age 5.
    expect(r.members).toEqual([2, 1, 0, 3, 4]);
  });

  it("a state named twice: a key table's row names the first of that name, the other is missing", async () => {
    const r = await members({ order: { data: "samples", field: "age" } }, ["n0", "n1", "n0", "n3", "n4"]);
    // Ages by name: n0 30, n1 20, n3 40; the second n0 and n4 have none (last, in data order).
    expect(r.members).toEqual([1, 0, 3, 2, 4]);
  });

  it("a field the table lacks is svs.layer.data at the layer's order, naming the fields it has; the key is skipped", async () => {
    const r = await members({ order: { data: "samples", field: "agee" } });
    const i = r.issues.find((x) => x.code === "svs.layer.data");
    expect(i?.path).toContain("params/order");
    expect(i?.message).toContain("population");
    expect(r.members).toEqual([0, 1, 2, 3, 4]);
  });
});
