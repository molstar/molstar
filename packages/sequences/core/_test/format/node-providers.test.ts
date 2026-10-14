/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// NP1a: the tree node providers `nj` and `upgma` (seqspec:model#alignments-and-trees): computed at load from the parent
// alignment's rows, checked like a given tree, cached by content and options, aborted by the stage's signal.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Alignment, Issue, SeqSet, Tree } from "@molstar/seq-core/model";
import { beforeEach, describe, expect, it } from "vitest";
import {
  compile,
  createRegistries,
  type NodeProviderEntry,
  register,
  type SeqViewSpecDoc,
  validate,
} from "../../src/format/index.ts";
import { KernelRegistry } from "../../src/format/kernel-registry.ts";
import { clearNodeProviderCache } from "../../src/format/load/node-providers.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";
import { scripted } from "./stream-provider.ts";

// Two close pairs, (a, b) and (c, d), and a distant e: a clade each, whatever the method.
const ROWS = [
  "ACDEFGHIKLMNPQRSTVWY",
  "ACDEFGHIKLMNPQRSTVWF",
  "WYVTSRQPNMLKIHGFEDCA",
  "WYVTSRQPNMLKIHGFEDCG",
  "AAAAAAAAAAAAAAAAAAAA",
];
const LABELS = ["a", "b", "c", "d", "e"];

function doc(
  provider: string | null,
  options: Record<string, unknown> = {},
  rows = ROWS,
  flavor?: "a3m",
): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "Guide tree" });
  const aln = b.data.alignment({
    ref: "aln",
    value: rows,
    labels: rows.map((_, i) => LABELS[i] ?? `r${i}`),
    ...(flavor ? { flavor } : {}),
  });
  if (provider) aln.tree({ provider, options } as never);
  const v = b.view({ ref: "v", label: "Rows", style: "msa" });
  v.block({ alignment: aln });
  v.track({ label: "Rows" }).rows({});
  return b.getDocument();
}
const first = (set: SeqSet) => set.alignments[0] as Alignment;
const errors = (issues: Issue[]) => issues.filter((i) => i.severity === "error");
const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const load = (d: SeqViewSpecDoc, signal?: AbortSignal) => compile(d, registries(), fixtureLoaders(), signal);

beforeEach(clearNodeProviderCache);

describe.each(["nj", "upgma"])("%s", (provider) => {
  it("computes a tree whose leaves are the rows, in a leaf order that keeps the close pairs together", async () => {
    const { set, issues } = await load(doc(provider));
    expect(errors(issues)).toEqual([]);
    const tree = first(set).tree as Tree;
    expect(tree?.leafOrder.length).toBe(5);
    expect(Array.from(tree?.leafOrder ?? []).sort()).toEqual([0, 1, 2, 3, 4]);
    const order = Array.from(tree?.leafOrder ?? []);
    const next = (a: number, b: number) => Math.abs(order.indexOf(a) - order.indexOf(b)) === 1;
    expect(next(0, 1)).toBe(true);
    expect(next(2, 3)).toBe(true);
    // The root has no branch; every other node does.
    const root = tree?.parent.indexOf(-1) ?? -1;
    expect(Number.isNaN(tree?.branchLength[root])).toBe(true);
  });
  it("takes the one distance, p, and says so for identity", () => {
    const reg = registries();
    expect(validate(doc(provider, { distance: "p" }), reg)).toEqual([]);
    expect(validate(doc(provider), reg)).toEqual([]);
    const issues = validate(doc(provider, { distance: "identity" }), reg);
    expect(issues.map((i) => i.code)).toEqual(["svs.params.type"]);
    expect(issues[0].message).toContain('write "p"');
    expect(codes(validate(doc(provider, { distance: "kimura" }), reg))).toEqual(["svs.params.type"]);
  });
});

describe("the distances", () => {
  it("compare the columns where both rows have a letter, gaps and an A3M's inserts out", async () => {
    const rows = ["AC-DE", "ACGDE", "AC-DF", "TT-TT"];
    const t = first((await load(doc("upgma", {}, rows))).set).tree as Tree;
    // a and b are equal where both have a letter: distance 0, joined first, at a height of 0.
    const [a, b] = [0, 1].map((r) => t?.leafRow.indexOf(r) ?? -1);
    expect(t?.parent[a]).toBe(t?.parent[b]);
    expect(t?.branchLength[a]).toBe(0);
    // The same rows as an A3M: lower case is an insert, not a column.
    const a3m = ["AC-DE", "AC-gDE", "AC-DF", "TT-TT"];
    const t3 = first((await load(doc("upgma", {}, a3m, "a3m"))).set).tree as Tree;
    expect(t3?.parent[t3.leafRow.indexOf(0)]).toBe(t3?.parent[t3.leafRow.indexOf(1)]);
  });
});

describe("a missing piece or a bad request", () => {
  it("leaves a provider that is not registered as an unavailable tree, quietly for the rest", async () => {
    const { set, issues } = await load(doc("nope"));
    expect(codes(issues)).toEqual(["svs.provider.unavailable"]);
    expect(first(set).tree).toBeNull();
    expect(first(set).rowCount).toBe(5);
  });
  it("is validated: an unknown distance, and a tree given twice", () => {
    const reg = registries();
    expect(codes(validate(doc("nj", { distance: "kimura" }), reg))).toContain("svs.params.type");
    const twice = doc("nj");
    const tree = twice.root.children?.[0].children?.[0].children?.[0];
    expect(tree?.kind).toBe("tree");
    (tree as NonNullable<typeof tree>).params = { provider: "nj", value: "(a,b);" };
    expect(codes(validate(twice, reg))).toContain("svs.params.exclusive");
  });
  it("is not computed over a streamed alignment: the alignment stays, the tree is reported", async () => {
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
    aln.tree({ provider: "nj", leaves: "order" } as never);
    b.view({ ref: "v", label: "V" }).block({ alignment: "aln" }).track({ label: "R" }).ruler();
    const { set, issues } = await compile(b.getDocument(), registries(), {
      load: (r) => Promise.reject(new Error(`no fixture ${r.uri}`)),
      stream: () => s.provider,
    });
    expect(codes(issues)).toEqual(["svs.provider.lazy"]);
    expect([first(set).lazy, first(set).tree]).toEqual([true, null]);
  });
  it("is unavailable, with its size code, past the row limit", async () => {
    const rows = Array.from({ length: 2049 }, () => "A");
    const { set, issues } = await compile(doc("nj", {}, rows), registries(), fixtureLoaders());
    expect(codes(issues)).toEqual(["svs.provider.size"]);
    expect(first(set).tree).toBeNull();
  });
});

describe("registered node providers", () => {
  it("run for a tree, and one that computes another kind is unavailable", async () => {
    const reg = createRegistries();
    register(reg, KernelRegistry);
    register(reg, {
      nodeProviders: [
        {
          name: "x:star",
          kind: "tree",
          options: KernelRegistry.nodeProviders?.[0].options as never,
          async compute(_ctx, inputs) {
            const n = (inputs[0] as { gapped: string[] }).gapped.length;
            return {
              parent: Int32Array.from({ length: n + 1 }, (_, i) => (i === n ? -1 : n)),
              branchLength: Float64Array.from({ length: n + 1 }, (_, i) => (i === n ? Number.NaN : 1)),
              label: Array.from({ length: n + 1 }, () => null),
              leafRow: Int32Array.from({ length: n + 1 }, (_, i) => (i === n ? -1 : i)),
            };
          },
        },
        {
          name: "x:wrong",
          kind: "mapping",
          options: KernelRegistry.nodeProviders?.[0].options as never,
          compute: async () => null,
        },
      ],
    });
    const star = await compile(doc("x:star"), reg, fixtureLoaders());
    expect(errors(star.issues)).toEqual([]);
    expect(first(star.set).tree?.count).toBe(6);
    const wrong = await compile(doc("x:wrong"), reg, fixtureLoaders());
    expect(codes(wrong.issues)).toEqual(["svs.provider.unavailable"]);
  });
  it("must be namespaced unless the kernel's", () => {
    const reg = registries();
    expect(() =>
      register(reg, {
        nodeProviders: [
          {
            name: "plain",
            kind: "tree",
            options: (reg.get("nodeProviders", "nj") as NodeProviderEntry).options,
            compute: async () => null,
          },
        ],
      }),
    ).toThrow(/not namespaced/);
  });
});

describe("caching and aborting", () => {
  it("gives the very same tree for the same rows and options, however often it is asked", async () => {
    const reg = registries();
    const nj = reg.get("nodeProviders", "nj") as NonNullable<ReturnType<typeof reg.get<"nodeProviders">>>;
    const input = { key: "aln", gapped: ROWS, flavor: "aligned" as const, labels: LABELS };
    const ask = (rows = ROWS, o = { distance: "p" }) =>
      nj.compute(
        { signal: new AbortController().signal } as never,
        [{ ...input, gapped: [...rows] }],
        o,
        new AbortController().signal,
      );
    const first = await ask();
    expect(await ask()).toBe(first);
    expect(await ask(ROWS, { distance: "q" })).not.toBe(first);
    expect(await ask([...ROWS.slice(0, 4), "CCCCCCCCCCCCCCCCCCCC"])).not.toBe(first);
  });
  it("keys content by each string's bounds and the options at every depth (GV-01)", async () => {
    const { contentKey } = await import("../../src/format/load/node-providers.ts");
    const rows = ["ACGT", "ACGA", "ACGG"];
    // Labels that read alike joined by newlines: another document's labels, before.
    expect(contentKey("nj", {}, rows, ["aligned", "a\nb", "c", "d"])).not.toBe(
      contentKey("nj", {}, rows, ["aligned", "a", "b\nc", "d"]),
    );
    expect(contentKey("nj", {}, ["AB", "C"], [])).not.toBe(contentKey("nj", {}, ["A", "BC"], []));
    // Nested option keys in another order are the same options (a hit, not a miss).
    expect(contentKey("realign", { a: { x: 1, y: 2 } }, rows, [])).toBe(
      contentKey("realign", { a: { y: 2, x: 1 } }, rows, []),
    );
  });

  it("keeps results by content: the same input gives the very same object, a changed row a new one", async () => {
    const { cached, inputKey } = await import("../../src/format/load/node-providers.ts");
    const input = (rows: string[]) => ({ key: "aln", gapped: rows, flavor: "aligned" as const, labels: LABELS });
    const k = inputKey("nj", { distance: "p" }, input(ROWS));
    expect(inputKey("nj", { distance: "p" }, input([...ROWS]))).toBe(k);
    expect(inputKey("nj", { distance: "q" }, input(ROWS))).not.toBe(k);
    expect(inputKey("upgma", { distance: "p" }, input(ROWS))).not.toBe(k);
    expect(inputKey("nj", { distance: "p" }, input([...ROWS.slice(0, 4), "AAAAAAAAAAAAAAAAAAAC"]))).not.toBe(k);
    let runs = 0;
    const tree = {
      parent: new Int32Array(1).fill(-1),
      branchLength: new Float64Array(1),
      label: [null],
      leafRow: new Int32Array(1),
    };
    const make = async () => {
      runs++;
      return tree;
    };
    expect(await cached(k, make)).toBe(tree);
    expect(await cached(k, make)).toBe(tree);
    expect(runs).toBe(1);
    // An aborted run is not kept.
    await expect(
      cached("other", async () => {
        throw new Error("aborted");
      }),
    ).rejects.toThrow();
    expect(await cached("other", make)).toBe(tree);
    expect(runs).toBe(2);
  });
  it("rejects a compile with the signal's reason when a newer request aborts a long build", async () => {
    const rows = Array.from({ length: 900 }, (_, i) =>
      Array.from({ length: 40 }, (_, c) => "ACDEFGHIKLMNPQRSTVWY"[(i * 7 + c * (i % 13) + (i >> 3)) % 20]).join(""),
    );
    const c = new AbortController();
    const running = compile(doc("nj", {}, rows), registries(), fixtureLoaders(), c.signal);
    setTimeout(() => c.abort(new Error("newer request")), 1);
    await expect(running).rejects.toThrow("newer request");
    // Nothing aborted is kept: the same request runs again to the end.
    const again = await compile(doc("nj", {}, rows.slice(0, 30)), registries(), fixtureLoaders());
    expect(first(again.set).tree?.leafOrder.length).toBe(30);
  });
});
