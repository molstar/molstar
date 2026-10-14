/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqSets } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  applyPatch,
  canonical,
  createRegistries,
  type DocNode,
  diff,
  int,
  normalize,
  nullable,
  type Op,
  OptionalField,
  type Patch,
  rebase,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  snapshotDoc,
  str,
  TreeSchema,
  UnionParamsSchema,
  validate,
} from "../../src/format/index.ts";
import { normalizeNodeParams } from "../../src/format/normalize.ts";
import { diffTree, patchTree } from "../../src/format/patch/tree-patch.ts";
import { itWithin } from "../model/timing.ts";
import { CASES, case2, case10 } from "./cases.ts";
import { registries } from "./load-fixtures.ts";

function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** A small schema: `n` simple, `m` a union whose discriminator switches null the old case's params. */
const kids = ["root", "n", "m"];
const schema = TreeSchema({
  rootKind: "root",
  nodes: {
    root: { description: "", parent: [], ref: "implied", params: SimpleParamsSchema({}) },
    n: {
      description: "",
      parent: kids,
      params: SimpleParamsSchema({ v: OptionalField(int, 0, ""), w: OptionalField(nullable(int), null, "") }),
    },
    m: {
      description: "",
      parent: kids,
      params: UnionParamsSchema("type", "", {
        x: SimpleParamsSchema({ p: OptionalField(int, 1, "") }),
        y: SimpleParamsSchema({ q: OptionalField(str, "", "") }),
      }),
    },
  },
});
const opts = { schema, reg: createRegistries() };

type Rnd = () => number;
const pick = <T>(rnd: Rnd, xs: T[]) => xs[Math.floor(rnd() * xs.length)] as T;
const POOL = Array.from({ length: 60 }, (_, i) => `r${i}`);

function node(rnd: Rnd, ref: string): DocNode {
  const custom = rnd() < 0.3 ? { [pick(rnd, ["ns1", "ns2"])]: { a: Math.floor(rnd() * 3) } } : undefined;
  const n: DocNode =
    rnd() < 0.5
      ? { kind: "n", ref, params: { v: Math.floor(rnd() * 5) } }
      : { kind: "m", ref, params: rnd() < 0.5 ? { type: "x", p: Math.floor(rnd() * 3) } : { type: "y", q: "s" } };
  if (custom) n.custom = custom;
  n.children = [];
  return n;
}

/** A random tree with unique refs from the pool. */
function tree(rnd: Rnd, size: number): DocNode {
  let next = 0;
  const make = (depth: number): DocNode => {
    const n = depth === 0 ? { kind: "root", children: [] as DocNode[] } : node(rnd, POOL[next++]);
    const count = depth < 3 ? Math.floor(rnd() * 4) : 0;
    for (let k = 0; k < count && next < size; k++) n.children?.push(make(depth + 1));
    return n;
  };
  return make(0);
}

function refs(n: DocNode, out: DocNode[] = []): DocNode[] {
  out.push(n);
  for (const c of n.children ?? []) refs(c, out);
  return out;
}
const find = (n: DocNode, ref: string) => refs(n).find((x) => x.ref === ref);
const anchorOf = (rnd: Rnd, parent: DocNode) =>
  pick(rnd, ["first", "last", ...(parent.children ?? []).map((c) => c.ref as string)]);

/** A random valid op: params and custom set (discriminators switched), removes, inserts reusing refs, moves. */
function randomOp(rnd: Rnd, doc: DocNode): Op | null {
  const all = refs(doc).slice(1);
  const r = rnd();
  if (r < 0.3 && all.length) {
    const n = pick(rnd, all);
    const params =
      n.kind === "n" ? { v: Math.floor(rnd() * 5), w: rnd() < 0.3 ? null : 1 } : { type: pick(rnd, ["x", "y"]) };
    const custom = rnd() < 0.5 ? { [pick(rnd, ["ns1", "ns2"])]: rnd() < 0.3 ? null : { a: 9 } } : undefined;
    return { op: "set", ref: n.ref as string, params, ...(custom && { custom }) };
  }
  if (r < 0.45 && all.length) return { op: "remove", ref: pick(rnd, all).ref as string };
  if (r < 0.75) {
    const free = POOL.filter((x) => !find(doc, x));
    const parent = pick(rnd, refs(doc));
    return {
      op: "insert",
      parent: parent.ref ?? "root",
      after: anchorOf(rnd, parent),
      node: node(rnd, pick(rnd, free)),
    };
  }
  if (!all.length) return null;
  const moving = pick(rnd, all);
  const targets = refs(doc).filter((c) => !find(moving, c.ref ?? "root"));
  const parent = pick(rnd, targets);
  const after = pick(rnd, [
    "first",
    "last",
    ...(parent.children ?? []).filter((c) => c !== moving).map((c) => c.ref as string),
  ]);
  return { op: "move", ref: moving.ref as string, parent: parent.ref ?? "root", after };
}

function edit(rnd: Rnd, root: DocNode, count: number): DocNode {
  let doc = root;
  for (let k = 0; k < count; k++) {
    const op = randomOp(rnd, doc);
    if (!op) continue;
    const r = patchTree(doc, [op], opts);
    expect(r.issues).toEqual([]);
    doc = r.root;
  }
  return doc;
}

/** The tree with params compared as normalized values (seqspec:format#what-a-change-affects). */
const norm = (n: DocNode): unknown => {
  const { params, custom } = normalizeNodeParams(n, schema, opts.reg);
  return { kind: n.kind, ref: n.ref, params, custom, children: (n.children ?? []).map(norm) };
};
const same = (a: DocNode, b: DocNode) => expect(canonical(norm(a))).toBe(canonical(norm(b)));

describe("tree patches (the spike's tests, with anchors)", () => {
  it("diff then apply gives the other version, both ways, for random edits with custom, kinds and refs reused", () => {
    for (let seed = 1; seed <= 300; seed++) {
      const rnd = random(seed);
      const a = tree(rnd, 40);
      const b = edit(rnd, a, 1 + Math.floor(rnd() * 8));
      for (const [x, y] of [
        [a, b],
        [b, a],
      ]) {
        const r = patchTree(x, diffTree(x, y, opts) as Op[], opts);
        expect(r.issues).toEqual([]);
        same(r.root, y);
      }
    }
  });

  it("the inverse undoes, exactly", () => {
    let ran = 0;
    for (let seed = 1; seed <= 200; seed++) {
      const rnd = random(seed);
      const a = tree(rnd, 30);
      const ops = Array.from({ length: 5 }, () => randomOp(rnd, a)).filter((o): o is Op => !!o);
      const r = patchTree(a, ops, opts);
      if (r.issues.length > 0) continue; // ops drawn over `a` may not compose
      ran++;
      const back = patchTree(r.root, r.inverse ?? [], opts);
      expect(back.issues).toEqual([]);
      expect(canonical(back.root)).toBe(canonical(a));
    }
    expect(ran).toBeGreaterThan(50);
  });

  it("diff removes a ref reused under another kind first, and inserts it with its descendants", () => {
    const y: DocNode = { kind: "n", ref: "y", params: { v: 2 } };
    const a: DocNode = { kind: "root", children: [{ kind: "n", ref: "x", params: { v: 1 }, children: [y] }] };
    const b: DocNode = { kind: "root", children: [{ kind: "m", ref: "x", params: { type: "x" }, children: [y] }] };
    const ops = diffTree(a, b, opts) as Op[];
    expect(ops.map((o) => o.op)).toEqual(["remove", "insert"]);
    expect(canonical(patchTree(a, ops, opts).root)).toBe(canonical(b));
  });

  it("shares every subtree, params and custom object a patch does not touch", () => {
    let seed = 3;
    while (refs(tree(random(seed), 60)).length < 12) seed++;
    const a = tree(random(seed), 60);
    const target = refs(a)[5].ref as string;
    const b = patchTree(a, [{ op: "set", ref: target, params: { v: 42 } }], opts).root;
    for (const n of refs(b)) {
      const old = find(a, n.ref ?? "") ?? a;
      if (find(n, target)) {
        expect(n.params === old.params).toBe(n.ref !== target);
        expect(n.custom).toBe(old.custom);
      } else expect(n).toBe(old);
    }
  });

  it("an empty diff for equal documents; nothing changes for no ops", () => {
    const a = tree(random(6), 30);
    expect(diffTree(a, JSON.parse(JSON.stringify(a)), opts)).toEqual([]);
    expect(patchTree(a, [], opts).root).toBe(a);
  });

  it("a discriminator switch nulls the old case's params; the inverse restores them", () => {
    const a: DocNode = { kind: "root", children: [{ kind: "m", ref: "m1", params: { type: "x", p: 2 } }] };
    const r = patchTree(a, [{ op: "set", ref: "m1", params: { type: "y", q: "z" } }], opts);
    expect(r.root.children?.[0].params).toEqual({ type: "y", q: "z" });
    expect(canonical(patchTree(r.root, r.inverse ?? [], opts).root)).toBe(canonical(a));
    expect(diffTree(a, r.root, opts)).toEqual([{ op: "set", ref: "m1", params: { type: "y", p: null, q: "z" } }]);
  });

  it("anchors, refs and cycles: all or nothing, pathed /ops/<k>", () => {
    const a: DocNode = {
      kind: "root",
      children: [
        { kind: "n", ref: "p", children: [{ kind: "n", ref: "c" }] },
        { kind: "n", ref: "q" },
      ],
    };
    const fails = (ops: Op[]) => {
      const r = patchTree(a, ops, opts);
      expect(r.root).toBe(a);
      return r.issues.map((i) => [i.code, i.path]);
    };
    const ok: Op = { op: "set", ref: "q", params: { v: 1 } };
    expect(fails([ok, { op: "insert", parent: "p", after: "q", node: { kind: "n", ref: "z" } }])).toEqual([
      ["svs.patch.anchor", "/ops/1"],
    ]);
    expect(
      fails([ok, { op: "insert", parent: "p", node: { kind: "n", ref: "z", children: [{ kind: "n", ref: "c" }] } }]),
    ).toEqual([["svs.patch.ref", "/ops/1"]]);
    expect(fails([{ op: "move", ref: "p", parent: "c" }])).toEqual([["svs.patch.cycle", "/ops/0"]]);
    expect(fails([{ op: "remove", ref: "root" }])).toEqual([["svs.patch.ref", "/ops/0"]]);
    expect(fails([{ op: "set", ref: "p/x", params: {} }])).toEqual([["svs.patch.ref", "/ops/0"]]);
    const moved = patchTree(a, [{ op: "move", ref: "q", parent: "p", after: "first" }], opts).root;
    expect(moved.children?.[0].children?.map((c) => c.ref)).toEqual(["q", "c"]);
  });
});

describe("SeqViewSpec patches", () => {
  const reg = registries();
  const revised = (doc: SeqViewSpecDoc): SeqViewSpecDoc => ({ ...doc, revision: "r1" });

  it("applies with the revision; the inverse undoes; untouched params are shared", () => {
    const a = revised(case10());
    const p = {
      format: "seq_view_spec_patch" as const,
      base: "r1",
      revision: "r2",
      ops: [
        { op: "set" as const, ref: "chains.main", params: { positions: "pocket", copies: "each" } },
        { op: "move" as const, ref: "strip", parent: "root", after: "first" },
      ],
    };
    const r = applyPatch(a, p, reg);
    expect(r.issues).toEqual([]);
    expect(r.doc.revision).toBe("r2");
    expect(validate(r.doc, reg)).toEqual([]);
    const data = (d: SeqViewSpecDoc) => d.root.children?.find((c) => c.kind === "data");
    expect(data(r.doc)).toBe(data(a));
    const view = (d: SeqViewSpecDoc, ref: string) => d.root.children?.find((c) => c.ref === ref);
    expect(view(r.doc, "chains")?.params).toBe(view(a, "chains")?.params);
    expect(view(r.doc, "alignment")).toBe(view(a, "alignment"));
    expect(r.inverse).toMatchObject({ base: "r2", revision: "r1" });
    const back = applyPatch(r.doc, r.inverse as never, reg);
    expect(back.doc).toEqual(a);
  });

  it("rejects a stale base, returning the document unchanged", () => {
    const a = revised(case2());
    const r = applyPatch(a, { format: "seq_view_spec_patch", base: "r0", ops: [] }, reg);
    expect([r.doc === a, r.inverse, r.issues.map((i) => i.code)]).toEqual([true, null, ["svs.patch.base"]]);
  });

  it("CF-09: a patch whose format is given and is not seq_view_spec_patch is svs.format.invalid", () => {
    const a = case2();
    const r = applyPatch(a, { format: "not_a_patch", ops: [] } as never, reg);
    expect([r.doc === a, r.inverse, r.issues.map((i) => [i.code, i.path])]).toEqual([
      true,
      null,
      [["svs.format.invalid", "/format"]],
    ]);
    expect(applyPatch(a, { ops: [] } as never, reg).issues).toEqual([]);
  });

  it("validates against a set's keys when given one (views over data given without a document)", () => {
    const { set } = SeqSets.build({ molecules: [{ key: "ext", residues: "MKV" }] });
    const block = { kind: "block", ref: "v.b", params: { sequence: "ext" } };
    const a: SeqViewSpecDoc = {
      format: "seq_view_spec",
      version: "1.0.0",
      root: {
        kind: "root",
        children: [
          { kind: "data" },
          { kind: "view", ref: "v", params: { label: "V" }, children: [{ kind: "section", ref: "v.main" }] },
        ],
      },
      snapshots: [{ ref: "s", ops: [{ op: "insert", parent: "v.main", node: block }] }],
    };
    const p = {
      format: "seq_view_spec_patch" as const,
      ops: [{ op: "insert" as const, parent: "v.main", node: { ...block, ref: "v.c" } }],
    };
    expect(applyPatch(a, p, reg).issues.map((i) => i.code)).toEqual(["svs.ref.unresolved"]);
    const r = applyPatch(a, p, reg, { set });
    expect([r.issues, r.inverse === null]).toEqual([[], false]);
    expect(snapshotDoc(a, "s", reg).issues.map((i) => i.code)).toEqual(["svs.ref.unresolved", "svs.snapshot.invalid"]);
    expect(snapshotDoc(a, "s", reg, { set }).issues).toEqual([]);
    expect(validate(a, reg, { set })).toEqual([]);
  });

  it("fails as a whole on an invalid result; warns per snapshot that no longer applies", () => {
    const a = case10();
    const bad = applyPatch(
      a,
      { format: "seq_view_spec_patch", ops: [{ op: "set", ref: "chains.main", params: { copies: "all" } }] },
      reg,
    );
    expect([bad.doc === a, bad.issues.map((i) => i.code)]).toEqual([true, ["svs.params.type"]]);
    const r = applyPatch(a, { format: "seq_view_spec_patch", ops: [{ op: "remove", ref: "chains.contacts" }] }, reg);
    expect(r.issues.map((i) => [i.code, i.severity, i.path])).toEqual([
      ["svs.snapshot.invalid", "warning", "/snapshots/2"],
    ]);
  });

  it("no diff, and no hang, when a ref is met twice (a root under the root, its `data` a second time)", () => {
    const reg = registries();
    const a = CASES.case10();
    const b = structuredClone(a);
    (b.root.children as DocNode[])[0].kind = "root";
    expect(diff(b, a, reg)).toBeNull();
    expect(diff(a, b, reg)).toBeNull();
  });

  it("a snapshot entry that is not one (null) is no snapshot: no throw, `svs.ref.unresolved`", () => {
    const doc = { ...CASES.case10(), snapshots: [null] } as unknown as SeqViewSpecDoc;
    expect(snapshotDoc(doc, undefined as never, registries()).issues.map((i) => i.code)).toEqual([
      "svs.ref.unresolved",
    ]);
  });

  it("each case's snapshots apply, and diff moves between them", () => {
    for (const make of Object.values(CASES)) {
      const doc = make();
      for (const s of doc.snapshots ?? []) {
        const r = snapshotDoc(doc, s.ref, reg);
        expect(r.issues).toEqual([]);
        const p = diff(doc, r.doc, reg) as Patch;
        expect(p.ops.length).toBe((s.ops ?? []).length);
        const back = applyPatch(doc, p, reg);
        expect(normalize(back.doc, reg).root).toEqual(normalize(r.doc, reg).root);
        expect(applyPatch(r.doc, diff(r.doc, doc) as Patch, reg).doc.root).toEqual(doc.root);
      }
    }
  });

  it("validates snapshot ops statically", () => {
    const doc = case10();
    const snaps = (ops: Op[][]) => validate({ ...doc, snapshots: ops.map((o, k) => ({ ref: `x${k}`, ops: o })) }, reg);
    const focus = (parent: string): Op => ({
      op: "insert",
      parent,
      node: { kind: "focus", ref: "f", params: { target: "1iep/A", start: 0, stop: 5 } },
    });
    const at = (issues: ReturnType<typeof validate>) => issues.map((i) => [i.code, i.path]);
    expect(at(snaps([[focus("strip.main")], [focus("chains.main")]]))).toEqual([
      ["svs.snapshot.ref", "/snapshots/1/ops/0"],
    ]);
    expect(at(snaps([[{ op: "insert", parent: "strip.main", node: { kind: "track", ref: "data" } }]]))).toEqual([
      ["svs.patch.ref", "/snapshots/0/ops/0"],
    ]);
    expect(at(snaps([[{ op: "set", ref: "strip.main", params: { height: "tall" } }]]))).toEqual([
      ["svs.snapshot.invalid", "/snapshots/0"],
    ]);
  });
});

describe("rebase at scale (CF-04)", () => {
  itWithin(
    "rebasing a document on itself, or on a copy of new objects, compares in place: no text of its data",
    (within) => {
      const reg = registries();
      const a = CASES.case1();
      // A 3,000² grid of floats, held where no schema reads it: canonical text of it took 5 s.
      a.root.custom = { acme: { grid: new Float32Array(3000 * 3000).map((_, i) => i % 97) } };
      expect(within(50, () => rebase(a, a, reg)).root).toBe(a.root);
      const copy = structuredClone(a);
      const r = within(1500, () => rebase(a, copy, reg));
      expect(r.root.custom).toBe(a.root.custom);
    },
  );
});
