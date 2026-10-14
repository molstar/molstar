/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer reads the shown set's systems through core's `Systems` (seqspec:viewer#host-coordinate-systems,
// seqspec:model#coordinate-systems; VS-27): an own system's descriptor lists a numbering a unit and its molecule share
// once, a ruler names a numbering by its namespace alone as a document does, and an alignment's columns keep the key
// `seqstar:<alignment>` (GS-13).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler } from "../src/index.js";
import { fixtureMeasurer } from "../src/testing/index.js";

async function open() {
  const b = createSeqViewSpec();
  b.data
    .sequence({ ref: "p", value: "MKVLATGHEK" })
    .numbering({ ref: "p.auth", name: "auth", method: "offset", first: 101 });
  const c = b.data.complex({ ref: "c" });
  c.unit({ ref: "u", sequence: "p" })
    .numbering({ ref: "u.auth", name: "auth", method: "offset", first: 201 })
    .numbering({ ref: "u.up", name: "uniprot:P04637", method: "offset", first: 1 });
  b.data.alignment({ ref: "msa", value: ["MKV-A", "MKVLA"], alphabet: "protein" });
  const v = b.view({ ref: "v", label: "V" });
  v.block({ unit: "u" } as never);
  v.track({ ref: "t" }).ruler({ ref: "r", numbering: "uniprot" } as never);
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry] },
    { scheduler: manualScheduler(), text: fixtureMeasurer() },
  );
  const r = await ctx.setDocument(b.getDocument());
  await ctx.graph.update();
  return { ctx, issues: r.issues };
}

describe("own systems through core's `Systems` (VS-27)", () => {
  it("a descriptor lists a numbering shared by a unit and its molecule once", async () => {
    const { ctx } = await open();
    expect(ctx.systems.descriptor("seqstar:u")?.numberings).toEqual(["auth", "uniprot:P04637", "ordinal"]);
  });

  it("a ruler names a numbering by its namespace alone", async () => {
    const { ctx, issues } = await open();
    expect(issues.filter((i) => i.severity === "error")).toEqual([]);
    const m = ctx.graph.assembled.value?.view as unknown as {
      resolved: { layers: { ref: Int32Array } };
      set: { strings: { get(i: number): string } };
      prepared: { state: { sources: ({ numberings: { name: number }[] }[] | null)[] } }[];
    };
    const k = Array.from(m.resolved.layers.ref).findIndex((x) => m.set.strings.get(x) === "r");
    const sources = m.prepared[k].state.sources[0];
    expect(sources?.map((s) => s.numberings.map((n) => m.set.strings.get(n.name)))).toEqual([["uniprot:P04637"]]);
  });

  it("an alignment's columns keep the key `seqstar:<alignment>`", async () => {
    const { ctx } = await open();
    const loci = ctx.refs.toLoci("msa", { start: 0, stop: 2 });
    expect(loci && ctx.systems.coordinatesOf(loci as never).systems).toEqual(["seqstar:msa"]);
    expect(ctx.systems.descriptor("seqstar:msa")?.kind).toBe("columns");
  });
});
