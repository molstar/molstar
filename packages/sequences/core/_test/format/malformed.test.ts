/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// V1's tests: malformed input always becomes issues, never an exception; the static checks match the
// spec. Each case here fails on the old stack.

import { SeqSets } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  applyPatch,
  Color,
  createRegistries,
  diff,
  KernelRegistry,
  type LayerTypeEntry,
  MvsRegistry,
  normalize,
  type Op,
  type Patch,
  rebase,
  register,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  snapshotDoc,
  validate,
} from "../../src/format/index.js";
import { namespaceOfCode, registeredNamespaces } from "../../src/format/registry.js";
import { CASES } from "./cases.js";
import { dataOf, fullDoc, type J, kernel, sectionOf } from "./fixtures.js";
import { registries } from "./load-fixtures.js";

const issues = (d: unknown, reg = kernel(), o = {}) => validate(d, reg, o).map((i) => [i.code, i.path]);
const errors = (d: unknown, reg = kernel(), o = {}) =>
  validate(d, reg, o)
    .filter((i) => i.severity === "error")
    .map((i) => [i.code, i.path]);

const D = "/root/children/0/children";
const S = "/root/children/1/children/0/children";
/** `v.ranges.l`, the ranges layer; its params. */
const RANGES = `${S}/3/children/0/children/0`;
const rangesLayer = (d: J) => sectionOf(d)[3].children[0].children[0];

describe("children that are absent or not a list", () => {
  it("a root, a view or a section without children reports the grammar counts", () => {
    expect(issues({ format: "seq_view_spec", version: "1.0.0", root: { kind: "root" } })).toEqual([
      ["svs.node.count", "/root"],
      ["svs.node.count", "/root"],
    ]);
    const d = fullDoc();
    delete d.root.children[1].children;
    expect(issues(d)).toEqual([["svs.node.count", "/root/children/1"]]);
    const e = fullDoc();
    delete e.root.children[1].children[0].children;
    expect(issues(e)).toEqual([]); // a section needs no children
  });

  /** Every node of the full document that holds children, by path. */
  const parents = (n: J, path: string, out: [J, string][] = []) => {
    if (n.children) {
      out.push([n, path]);
      for (let i = 0; i < n.children.length; i++) parents(n.children[i], `${path}/children/${i}`, out);
    }
    return out;
  };

  it("children that are an object, a string, a number or null: one type error, and no exception", () => {
    const paths = parents(fullDoc().root, "/root").map(([, p]) => p);
    expect(paths.length).toBeGreaterThan(10);
    for (const path of paths) {
      for (const bad of [{}, "x", 5, null, true]) {
        const d = fullDoc();
        const at = parents(d.root, "/root").find(([, p]) => p === path) as [J, string];
        at[0].children = bad;
        const own = errors(d).filter(([code, p]) => code === "svs.params.type" && p === `${path}/children`);
        expect(own, `${path} = ${JSON.stringify(bad)}`).toHaveLength(1);
      }
    }
  });

  it("a child that is not a node is a type error at its index", () => {
    const d = fullDoc();
    d.root.children[1].children.push(null, 5, "x", {});
    expect(errors(d).map(([, p]) => p)).toEqual([
      "/root/children/1/children/1",
      "/root/children/1/children/2",
      "/root/children/1/children/3",
      "/root/children/1/children/4",
    ]);
  });
});

describe("static checks", () => {
  /** [what, expected [code, path] list, the mutation]. */
  const cases: [string, string[][], (d: J) => void][] = [
    [
      "a letter dictionary value that is not a string",
      [["svs.params.type", `${RANGES}/params/color/params/colors/A`]],
      (d) => (rangesLayer(d).params.color = { theme: "letter", params: { colors: { A: null } } }),
    ],
    [
      "a letter dictionary value that is not a color",
      [["svs.color.invalid", `${RANGES}/params/color/params/colors/A`]],
      (d) => (rangesLayer(d).params.color = { theme: "letter", params: { colors: { A: "notacolor" } } }),
    ],
    [
      "a categorical palette's dictionary value",
      [["svs.params.type", `${RANGES}/params/color/params/palette/colors/x`]],
      (d) => {
        rangesLayer(d).params.color.params.palette = { kind: "categorical", colors: { x: 5 } };
      },
    ],
    [
      "a categorical palette's dictionary color",
      [["svs.color.invalid", `${RANGES}/params/color/params/palette/colors/x`]],
      (d) => {
        rangesLayer(d).params.color.params.palette = { kind: "categorical", colors: { x: "#ggg" } };
      },
    ],
    [
      "an extension type's data against its fallback, unregistered",
      [["svs.layer.data", `${RANGES}/params/data`]],
      (d) => Object.assign(rangesLayer(d).params, { type: "org.example:heat", fallback: "values", data: "regions" }),
    ],
    [
      "a protocol source without a uri",
      [["svs.params.missing", `${D}/12/params/uri`]],
      (d) =>
        dataOf(d).push({
          kind: "source",
          ref: "tiles",
          params: { service: "uniprot", id: "P1", protocol: "seqstar-tiles/1" },
        }),
    ],
    ["max_bytes of 0", [["svs.params.type", `${D}/0/params/max_bytes`]], (d) => (dataOf(d)[0].params.max_bytes = 0)],
    [
      "a negative max_bytes",
      [["svs.params.type", `${D}/0/params/max_bytes`]],
      (d) => (dataOf(d)[0].params.max_bytes = -1),
    ],
    [
      "a name with whitespace in a names numbering",
      [["svs.position.syntax", `${D}/1/children/1/params/labels/1`]],
      (d) =>
        dataOf(d)[1].children.push({
          kind: "numbering",
          ref: "p53.chain",
          params: { name: "chain", method: "names", labels: ["N1", "C 1"] },
        }),
    ],
    [
      "a glyph box that is not two positive numbers",
      [["svs.glyph.invalid", `${RANGES}/params/glyph/params/glyph/box`]],
      (d) => {
        rangesLayer(d).params.glyph.params.glyph = { path: "M0 0L1 1", box: ["a", -1] };
      },
    ],
    [
      "a glyph box of 0",
      [["svs.glyph.invalid", `${RANGES}/params/glyph/params/glyph/box`]],
      (d) => {
        rangesLayer(d).params.glyph.params.glyph = { path: "M0 0L1 1", box: [0, 5] };
      },
    ],
    [
      "a glyph box without a path",
      [["svs.glyph.invalid", `${RANGES}/params/glyph/params/glyph/box`]],
      (d) => {
        rangesLayer(d).params.glyph.params.glyph = { shape: "circle", box: [10, 10] };
      },
    ],
    ["snapshot ref next", [["svs.ref.syntax", "/snapshots/1/ref"]], (d) => d.snapshots.push({ ref: "next" })],
    ["snapshot ref prev", [["svs.ref.syntax", "/snapshots/1/ref"]], (d) => d.snapshots.push({ ref: "prev" })],
    [
      "snapshot ref with a space",
      [["svs.ref.syntax", "/snapshots/1/ref"]],
      (d) => d.snapshots.push({ ref: "after the fix" }),
    ],
    [
      "a type in a namespace that no extension can register",
      [["svs.params.type", `${RANGES}/params/type`]],
      (d) => Object.assign(rangesLayer(d).params, { type: "Org.Example:heat", fallback: "ranges" }),
    ],
  ];
  for (const [what, expected, mutate] of cases) {
    it(what, () => {
      const d = fullDoc();
      mutate(d);
      expect(errors(d)).toEqual(expected);
    });
  }

  it("a color that is not a string decodes to the fallback, never throws", () => {
    for (const v of [null, undefined, 5, {}, ["#fff"]]) expect(Color.decodePrimary(v)).toBe(0x000000);
    expect(Color.decodePrimary("#fff/red")).toBe(0xffffff);
  });

  it("valid dictionaries and a path glyph with a good box still validate", () => {
    const d = fullDoc();
    rangesLayer(d).params.color = { theme: "letter", params: { colors: { A: "#fff", C: "red" } } };
    rangesLayer(d).params.glyph.params.glyph = { path: "M0 0L1 1", box: [10, 12] };
    expect(issues(d)).toEqual([]);
  });

  it("a registered extension type's data is checked against its fallback and its own accepts", () => {
    const entry = (accepts: LayerTypeEntry["accepts"]): LayerTypeEntry => ({
      name: "org.example:heat",
      options: SimpleParamsSchema({}),
      accepts,
      layouts: ["tracks"],
      axes: ["blocks"],
      sizing: "row",
      fallback: "values",
      channels: [],
    });
    const doc = (fallback: string) => {
      const d = fullDoc();
      Object.assign(rangesLayer(d).params, { type: "org.example:heat", fallback, data: "regions" });
      return d;
    };
    const complaints = (d: J, reg = kernel()) =>
      validate(d, reg)
        .filter((i) => i.code === "svs.layer.data")
        .map((i) => i.message.split(" takes ")[0]);
    // Ranges data against the fallback `values`: refused whether or not the type is registered, and by the type too
    // when it takes only values.
    for (const [accepts, expected] of [
      [["ranges"], ["values"]],
      [["values"], ["values", "org.example:heat"]],
    ] as const) {
      const reg = kernel();
      register(reg, { layerTypes: [entry(accepts)] });
      expect(complaints(doc("values"), reg)).toEqual(expected);
    }
    const reg = kernel();
    register(reg, { layerTypes: [entry(["values"])] });
    expect(complaints(doc("ranges"), reg)).toEqual(["org.example:heat"]);
    expect(complaints(doc("ranges"))).toEqual([]); // unregistered, the fallback's `ranges` takes it
  });

  it("a key of the set that is not the kind the param names", () => {
    const { set } = SeqSets.build({
      molecules: [
        { key: "ext", residues: "MKV" },
        { key: "het", kind: "ligand", elements: ["C", "N"] },
      ],
    });
    const doc = (params: J): SeqViewSpecDoc => ({
      format: "seq_view_spec",
      version: "1.0.0",
      root: {
        kind: "root",
        children: [
          { kind: "data" },
          {
            kind: "view",
            ref: "v",
            params: { label: "V" },
            children: [{ kind: "section", ref: "v.main", children: [{ kind: "block", ref: "v.b", params }] }],
          },
        ],
      },
    });
    const at = "/root/children/1/children/0/children/0/params";
    expect(errors(doc({ sequence: "ext" }), kernel(), { set })).toEqual([]);
    expect(errors(doc({ ligand: "het" }), kernel(), { set })).toEqual([]);
    expect(errors(doc({ alignment: "ext" }), kernel(), { set })).toEqual([["svs.ref.wrong_kind", `${at}/alignment`]]);
    expect(errors(doc({ sequence: "het" }), kernel(), { set })).toEqual([["svs.ref.wrong_kind", `${at}/sequence`]]);
    expect(errors(doc({ unit: "ext" }), kernel(), { set })).toEqual([["svs.ref.wrong_kind", `${at}/unit`]]);
  });

  it("normalize keeps the unknown keys of a ThemeRef", () => {
    const d = fullDoc();
    Object.assign(rangesLayer(d).params.color, { shiny: 1 });
    const out: J = normalize(d, kernel());
    expect(out.root.children[1].children[0].children[3].children[0].children[0].params.color.shiny).toBe(1);
  });

  it("unknown glyph keys are warnings in a newer minor, errors otherwise", () => {
    const d = fullDoc();
    rangesLayer(d).params.glyph.params.glyph = { char: "x", zoom: 2 };
    const at = `${RANGES}/params/glyph/params/glyph/zoom`;
    expect(issues(d)).toEqual([["svs.glyph.invalid", at]]);
    d.version = "1.1.0";
    expect(validate(d, kernel()).map((i) => [i.code, i.severity, i.path])).toEqual([
      ["svs.version.newer_minor", "warning", "/version"],
      ["svs.version.newer_minor", "warning", at],
    ]);
  });
});

describe("namespaces", () => {
  const type = (name: string): LayerTypeEntry => ({
    name,
    options: SimpleParamsSchema({}),
    accepts: ["ranges"],
    layouts: ["tracks"],
    axes: ["blocks"],
    sizing: "row",
    fallback: "ranges",
    channels: [],
  });
  const conflicts = (entry: Parameters<typeof register>[1]) => {
    const reg = createRegistries();
    register(reg, KernelRegistry);
    try {
      register(reg, entry);
    } catch (e) {
      return (e as { conflicts: string[] }).conflicts;
    }
    return [];
  };

  it("mvs is reserved: no extension registers a custom table, a name or a code in it; the kernel's own table is exempt", () => {
    const schema = SimpleParamsSchema({});
    expect(conflicts({ custom: [{ namespace: "mvs", kind: "other", schema }] })).toEqual(['custom: "mvs" is reserved']);
    expect(conflicts({ layerTypes: [type("mvs:x")] })).toEqual(['layerTypes: "mvs" is reserved']);
    expect(conflicts({ issueCodes: [{ code: "mvs.bad", severity: "error" }] })).toEqual([
      'issueCodes: "mvs.bad" is reserved',
    ]);
    expect(conflicts(MvsRegistry)).toEqual([]);
  });

  it("one syntax: a name in a namespace with capitals cannot be registered, and a document cannot use it", () => {
    expect(conflicts({ layerTypes: [type("Org.Example:heat")] })).toEqual([
      'layerTypes: "Org.Example:heat" is not namespaced (<ns>:<name>)',
    ]);
    expect(conflicts({ layerTypes: [type("org.example:Heat")] })).toEqual([]);
  });

  it("an issue code belongs to the longest registered namespace it starts with, and owns no namespace by itself", () => {
    expect(namespaceOfCode(new Set(["org", "org.example"]), "org.example.bad")).toBe("org.example");
    expect(namespaceOfCode(new Set(["org"]), "org.example.bad")).toBe("org");
    expect(namespaceOfCode(new Set(["lab"]), "org.example.bad")).toBeNull();
    expect(namespaceOfCode(new Set(), "svs.params.type")).toBe("svs");
    const reg = createRegistries();
    register(reg, [
      KernelRegistry,
      { layerTypes: [type("org.example:heat")], issueCodes: [{ code: "org.example.bad", severity: "warning" }] },
      { issueCodes: [{ code: "lab.stray", severity: "warning" }] },
    ]);
    expect([...registeredNamespaces(reg)].sort()).toEqual(["org.example", "svs"]);
    const doc = (requires: string[]) => ({ ...fullDoc(), requires });
    expect(errors(doc(["org"]), reg)).toEqual([["svs.requires.unknown", "/requires/0"]]);
    expect(errors(doc(["lab"]), reg)).toEqual([["svs.requires.unknown", "/requires/0"]]);
    expect(errors(doc(["org.example"]), reg)).toEqual([]);
  });
});

describe("malformed patches", () => {
  const reg = registries();
  const doc = (): SeqViewSpecDoc => ({ ...fullDoc(), revision: "r1" });
  const ok: Op = { op: "set", ref: "v.main", params: { numbering: "kabat" } };
  const track = (more: J): J => ({ kind: "track", ref: "t.new", ...more });

  /** [what, a malformed op, the code]: each alone is the one issue of its patch, at /ops/1. */
  const bad: [string, unknown, string?][] = [
    ["set: params a number", { op: "set", ref: "v.main", params: 5 }],
    ["set: params a list", { op: "set", ref: "v.main", params: [1] }],
    ["set: params null", { op: "set", ref: "v.main", params: null }],
    ["set: custom a string", { op: "set", ref: "v.main", custom: "x" }],
    ["insert: children an object", { op: "insert", parent: "v.main", node: track({ children: {} }) }],
    ["insert: children a string", { op: "insert", parent: "v.main", node: track({ children: "x" }) }],
    ["insert: children a list with null", { op: "insert", parent: "v.main", node: track({ children: [null] }) }],
    ["insert: children a list with a number", { op: "insert", parent: "v.main", node: track({ children: [5] }) }],
    [
      "insert: a child's params a number",
      { op: "insert", parent: "v.main", node: track({ children: [{ kind: "layer", ref: "l", params: 5 }] }) },
    ],
    [
      "insert: a grandchild's children an object",
      { op: "insert", parent: "v.main", node: track({ children: [{ kind: "layer", ref: "l", children: {} }] }) },
    ],
    ["insert: params a number", { op: "insert", parent: "v.main", node: track({ params: 5 }) }],
    ["insert: node a number", { op: "insert", parent: "v.main", node: 5 }],
    ["insert: node null", { op: "insert", parent: "v.main", node: null }],
    ["insert: node a list", { op: "insert", parent: "v.main", node: [] }],
    ["insert: node without a kind", { op: "insert", parent: "v.main", node: { ref: "x" } }],
    ["insert: node without a ref", { op: "insert", parent: "v.main", node: { kind: "track" } }, "svs.ref.missing"],
    [
      "insert: a descendant without a ref",
      { op: "insert", parent: "v.main", node: track({ children: [{ kind: "layer" }] }) },
      "svs.ref.missing",
    ],
    ["an op that is a number", 5],
    ["an op that is null", null],
    ["an op that is a string", "set"],
    ["an op that is a list", []],
    ["an op without op", { ref: "v.main" }],
    ["an unknown op", { op: "swap", ref: "v.main" }],
    ["set: ref a number", { op: "set", ref: 5, params: {} }],
    ["remove: no ref", { op: "remove" }],
    ["insert: no parent", { op: "insert", node: track({}) }],
    ["move: parent a list", { op: "move", ref: "v.main", parent: ["root"] }],
    ["move: after a number", { op: "move", ref: "v.main", parent: "root", after: 0 }],
  ];
  for (const [what, op, code] of bad) {
    it(`${what}: one issue at /ops/1, nothing applied`, () => {
      const a = doc();
      const r = applyPatch(a, { format: "seq_view_spec_patch", ops: [ok, op as Op] }, reg);
      expect(r.issues.map((i) => [i.code, i.path])).toEqual([[code ?? "svs.params.type", "/ops/1"]]);
      expect([r.doc === a, r.inverse]).toEqual([true, null]);
    });
  }

  it("every malformed op is reported, each once", () => {
    const r = applyPatch(
      doc(),
      { format: "seq_view_spec_patch", ops: [5, ok, { op: "set", ref: "v.main", params: 1 }] as Op[] },
      reg,
    );
    expect(r.issues.map((i) => [i.code, i.path])).toEqual([
      ["svs.params.type", "/ops/0"],
      ["svs.params.type", "/ops/2"],
    ]);
  });

  it("ops that are not a list, absent, or a patch that is not an object: one issue", () => {
    for (const patch of [
      { format: "seq_view_spec_patch", ops: {} },
      { format: "seq_view_spec_patch", ops: "x" },
      { format: "seq_view_spec_patch", ops: null },
      { format: "seq_view_spec_patch" },
    ]) {
      const r = applyPatch(doc(), patch as never, reg);
      expect(r.issues.map((i) => [i.code, i.path])).toEqual([["svs.params.type", "/ops"]]);
    }
    for (const patch of [null, 5, "x", []]) {
      const r = applyPatch(doc(), patch as never, reg);
      expect(r.issues.map((i) => [i.code, i.path])).toEqual([["svs.params.type", ""]]);
    }
  });

  it("base: null applies to any revision, as no base does; a wrong base is rejected", () => {
    const a = doc();
    const p = (base: unknown) => ({ format: "seq_view_spec_patch" as const, base, ops: [ok] }) as unknown as Patch;
    for (const base of [null, undefined]) {
      const r = applyPatch(a, p(base), reg);
      expect([r.issues, r.inverse === null]).toEqual([[], false]);
    }
    expect(applyPatch(a, p("r0"), reg).issues.map((i) => i.code)).toEqual(["svs.patch.base"]);
  });

  it("snapshot ops are checked the same", () => {
    const d = fullDoc();
    d.snapshots = [{ ref: "s1", ops: [{ op: "set", ref: "v.main", params: 5 }] }];
    expect(issues(d)).toEqual([["svs.params.type", "/snapshots/0/ops/0"]]);
    d.snapshots = [{ ref: "s1", ops: {} }];
    expect(issues(d)).toEqual([["svs.params.type", "/snapshots/0/ops"]]);
    expect(snapshotDoc(d, "s1", kernel()).issues.map((i) => i.code)).toContain("svs.params.type");
  });

  it("a document with a node without a ref cannot be patched: svs.ref.missing, not an exception", () => {
    const a = doc();
    delete (dataOf(a)[2] as J).ref;
    const r = applyPatch(a, { format: "seq_view_spec_patch", ops: [ok] }, reg);
    expect(r.issues.map((i) => [i.code, i.path])).toEqual([["svs.ref.missing", "/ops"]]);
    expect(r.issues[0].message).toContain("/root/children/0/children/2");
  });

  it("a snapshot inserting a node whose children are not a list is an issue, not an exception", () => {
    for (const children of [{}, 5, "x"]) {
      const d = fullDoc();
      d.snapshots = [{ ref: "s1", ops: [{ op: "insert", parent: "v.main", node: track({ children }) }] }];
      expect(issues(d)).toEqual([["svs.params.type", "/snapshots/0/ops/0"]]);
    }
  });
});

describe("issue paths of a patch's result", () => {
  const reg = registries();
  const patch = (...ops: Op[]): Patch => ({ format: "seq_view_spec_patch", ops });

  it("an issue of the result is pathed at the last op that touched its node, naming the node and its path", () => {
    const a = fullDoc();
    const r = applyPatch(
      a,
      patch(
        { op: "set", ref: "v.b2", params: { label: "two" } },
        { op: "set", ref: "v.b1", params: { sequence: "nope" } },
        { op: "set", ref: "v.b1", params: { label: "one" } },
      ),
      reg,
    );
    expect(r.issues.map((i) => [i.code, i.path])).toEqual([["svs.ref.unresolved", "/ops/2"]]);
    expect(r.issues[0].message).toContain(`(at v.b1: ${S}/0/params/sequence)`);
  });

  it("an issue at an inserted node, or at a parent whose children changed, is pathed at that op", () => {
    const a = fullDoc();
    const added = applyPatch(
      a,
      patch({ op: "insert", parent: "v.main", node: { kind: "block", ref: "v.b3", params: { sequence: "nope" } } }),
      reg,
    );
    expect(added.issues.map((i) => [i.code, i.path])).toEqual([["svs.ref.unresolved", "/ops/0"]]);
    const removed = applyPatch(a, patch({ op: "remove", ref: "v.b1" }, { op: "remove", ref: "v.b2" }), reg);
    expect(removed.issues).toEqual([]); // a section needs no blocks
  });

  it("an issue at a node no op touched is pathed /ops, still naming the node", () => {
    const r = applyPatch(fullDoc(), patch({ op: "remove", ref: "regions" }), reg);
    expect(r.issues.length).toBeGreaterThan(0);
    for (const i of r.issues) expect(i.path).toBe("/ops");
    expect(r.issues.map((i) => i.message).join("\n")).toMatch(
      /\(at v\.main: \/root\/children\/1\/children\/0\/params\/positions\)/,
    );
  });

  it("the same for a snapshot: /snapshots/<k>/ops/<j>", () => {
    const d = fullDoc();
    d.snapshots = [
      {
        ref: "bad",
        ops: [
          { op: "set", ref: "v.b1", params: { label: "x" } },
          { op: "set", ref: "v.b1", params: { sequence: "nope" } },
        ],
      },
    ];
    const r = snapshotDoc(d, "bad", reg);
    expect(r.issues.map((i) => [i.code, i.path])).toEqual([
      ["svs.ref.unresolved", "/snapshots/0/ops/1"],
      ["svs.snapshot.invalid", "/snapshots/0"],
    ]);
    expect(r.issues[0].message).toContain("(at v.b1: ");
  });
});

describe("diff and rebase", () => {
  const reg = registries();
  const docs = Object.entries(CASES).map(([name, make]) => [name, make] as const);

  it("diff applied to a structured clone of a gives b, for every case and snapshot", () => {
    let checked = 0;
    for (const [, make] of docs) {
      const a = make();
      for (const s of a.snapshots ?? []) {
        const b = snapshotDoc(a, s.ref, reg).doc;
        const p = diff(a, b);
        const r = applyPatch(structuredClone(a), p as Patch, reg);
        expect(r.issues.filter((i) => i.severity === "error")).toEqual([]);
        expect(normalize(r.doc, reg).root).toEqual(normalize(b, reg).root);
        checked++;
      }
      expect(diff(a, structuredClone(a))?.ops).toEqual([]);
    }
    expect(checked).toBeGreaterThan(2);
  });

  it("diff is null when either document holds a node without a ref, and the implied ones need none", () => {
    const a = CASES.case10();
    const bare = structuredClone(a);
    const nodes = (n: J, out: J[] = []): J[] => {
      out.push(n);
      for (const c of n.children ?? []) nodes(c, out);
      return out;
    };
    expect(
      nodes(bare.root)
        .filter((n) => n.ref === undefined)
        .map((n) => n.kind),
    ).toEqual(["root", "data"]);
    expect(diff(a, bare)?.ops).toEqual([]);
    const victim = nodes(bare.root).find((n) => n.kind === "view") as J;
    delete victim.ref;
    expect(diff(a, bare)).toBeNull();
    expect(diff(bare, a)).toBeNull();
    for (const root of [
      { kind: "root", children: {} },
      { kind: "root", children: [5] },
      { kind: "root", children: [null] },
    ]) {
      expect(diff(a, { ...a, root } as never)).toBeNull();
    }
  });

  it("rebase gives next sharing what equals base's, and next itself when there is no diff", () => {
    const a = CASES.case10();
    const next = structuredClone(a);
    const same = rebase(a, next, reg);
    expect(same.root).toBe(a.root);
    const view = (d: SeqViewSpecDoc) => d.root.children?.find((c) => c.ref === "chains") as J;
    const changed = structuredClone(a);
    view(changed).children[0].params = { ...view(changed).children[0].params, positions: "pocket" };
    const r = rebase(a, changed, reg);
    expect(view(r).children[0].params.positions).toBe("pocket");
    const data = (d: SeqViewSpecDoc) => d.root.children?.find((c) => c.kind === "data");
    expect(data(r)).toBe(data(a));
    expect(view(r).children[1]).toBe(view(a).children[1]);
    const bare = structuredClone(a);
    delete (view(bare) as J).ref;
    expect(rebase(a, bare, reg)).toBe(bare);
  });
});
