/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The graph's assembly stages what `compile` does (seqspec:viewer#updates-and-status; VS-25): a node that reads the set
// sees the computed alignments before it, a complex by value waits for its units' availability, a unit naming the wrong
// kind is left out (core's `checkUnitKinds`), and the assembly's issues name document paths (core's `remapIssues`).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { compile, type SeqViewSpecDoc, SimpleParamsSchema } from "@molstar/seq-core/format";
import { KeyKind, type KeyTarget, type SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { fixtureLoaders } from "../../core/_test/format/load-fixtures.js";
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler } from "../src/index.js";
import { fixtureMeasurer } from "../src/testing/index.js";

const ligands = {
  svs: {
    complexProviders: [
      {
        name: "org.example:ligands",
        options: SimpleParamsSchema({}),
        formats: ["json"],
        build: async () => ({ molecules: [{ key: "LIG", kind: "ligand" as const, elements: ["C"] }], units: [] }),
      },
    ],
  },
};

/** The graph's set and issues, and `compile`'s, for one document over the same registries. */
async function both(make: () => SeqViewSpecDoc) {
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, ligands] },
    { scheduler: manualScheduler(), loaders: fixtureLoaders(), text: fixtureMeasurer() },
  );
  const compiled = await compile(make(), ctx.registry.svs, fixtureLoaders());
  const r = await ctx.setDocument(make());
  await ctx.graph.update();
  const set = ctx.graph.assembled.value?.set as SeqSet;
  ctx.dispose();
  const issues = (is: readonly { code: string; path: string }[]) => is.map((i) => [i.code, i.path]).sort();
  return { set, issues: issues(r.issues), compiled: compiled.set, compiledIssues: issues(compiled.issues) };
}
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const has = (set: SeqSet, key: string) => set.byName(key, t);
const complexUnits = (set: SeqSet, key: string) =>
  set.byName(key, t) && t.kind === KeyKind.complex ? (set.complexes[t.id]?.units.length ?? -1) : -1;

describe("the assembly stages what `compile` does (VS-25)", () => {
  it("a computed alignment over an earlier one, and an `align` mapping from its row, are built", async () => {
    const make = () => {
      const b = createSeqViewSpec();
      for (const [k, v] of ["MKVLAGWTQRS", "MKVLAGWTRRS", "KVLAGWTQRS"].entries())
        b.data.sequence({ ref: `s${k}`, value: v });
      b.data.alignment({ ref: "A", provider: "realign", inputs: ["s0", "s1", "s2"] } as never);
      b.data.alignment({ ref: "B", provider: "realign", inputs: "A" } as never);
      b.data.sequence({ ref: "q", value: "MKVLAGWTQRS" });
      b.data.mapping({ ref: "m", from: "A/0", to: "q", method: "align" } as never);
      b.view({ ref: "v", label: "V", style: "msa" }).block({ alignment: "B" } as never);
      return b.getDocument();
    };
    const r = await both(make);
    expect([has(r.compiled, "B"), has(r.compiled, "m")]).toEqual([true, true]);
    expect([has(r.set, "B"), has(r.set, "m")]).toEqual([true, true]);
    expect(r.issues).toEqual(r.compiledIssues);
  });

  it("a complex by value with an unavailable unit reports only that unit, as `compile` does", async () => {
    const make = () => {
      const b = createSeqViewSpec();
      const complex = b.data.complex({ ref: "c" });
      complex.unit({ ref: "ua", sequence: "a" } as never);
      complex.unit({ ref: "ub", sequence: "b" } as never);
      const ext = b.data.source({ ref: "ext", value: ">x\nMKV\n", format: "genbank" });
      b.data.sequence({ ref: "a", source_ref: ext });
      b.data.sequence({ ref: "b", value: "MKVLA" });
      b.view({ ref: "v", label: "V" }).block({ complex }).track({ label: "L" }).letters();
      return b.getDocument();
    };
    const r = await both(make);
    expect(r.compiledIssues.map((i) => i[0])).toEqual(["svs.provider.unavailable"]);
    expect(r.issues).toEqual(r.compiledIssues);
    expect(complexUnits(r.set, "c")).toBe(1);
  });

  it("a unit naming a ligand as a sequence is left out, of its complex too, with `compile`'s issue", async () => {
    const make = () => {
      const b = createSeqViewSpec();
      b.data.complex({
        ref: "prov",
        provider: "org.example:ligands",
        source: { value: "{}", format: "json" },
      } as never);
      const c = b.data.complex({ ref: "c" });
      c.unit({ ref: "bad", sequence: "prov/LIG" } as never);
      c.unit({ ref: "good", sequence: b.data.sequence({ ref: "s", value: "MKV" }) });
      b.view({ ref: "v", label: "V" }).block({ complex: "c" } as never);
      return b.getDocument();
    };
    const r = await both(make);
    expect(r.compiledIssues.map((i) => i[0])).toEqual(["svs.ref.wrong_kind"]);
    expect(r.issues).toEqual(r.compiledIssues);
    expect([has(r.set, "bad"), complexUnits(r.set, "c")]).toEqual([false, 1]);
  });

  it("the assembly's issues name the document node, not the model's input", async () => {
    const make = () => {
      const b = createSeqViewSpec();
      b.data.sequence({ ref: "s", value: "MKV" });
      b.data.annotation({ ref: "a", geometry: "ranges", target: "s", value: { start: [0], stop: [9] } });
      b.view({ ref: "v", label: "V" }).track({ label: "L" }).letters();
      return b.getDocument();
    };
    const r = await both(make);
    expect(r.compiledIssues.length).toBeGreaterThan(0);
    expect(r.issues).toEqual(r.compiledIssues);
  });
});
