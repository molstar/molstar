/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// GI-14: onto's row maps live for one source's routes; a throw while building them leaves nothing for the next run.
import { bindQuery, readQuery } from "@molstar/seq-core/query";
import { describe, expect, it, vi } from "vitest";
import { SeqSets } from "../../src/model/index.js";
import { Mappings } from "../../src/model/mapping.js";
import { compileQuery, createQueryEnv, SeqSelection } from "../../src/model/query/index.js";
import { ontoFixture } from "./query-fixtures.js";

const ENV = createQueryEnv([]);
const run = (set: unknown, text: string) => {
  const b = bindQuery(readQuery(text).node as never, ENV.language);
  const r = compileQuery(b.expr as never, ENV, b.source).run({ set, inputs: {} } as never);
  return SeqSelection.summary(r.selection, { maxParts: 8 });
};

describe("onto's row maps", () => {
  it("a run that throws mid-route leaves no row maps for another set's run", () => {
    const A = SeqSets.build(ontoFixture()).set;
    const rotated = ontoFixture();
    const aln = rotated.alignments?.[0];
    if (!aln) throw new Error("no alignment");
    rotated.alignments = [{ ...aln, gapped: aln.gapped.map((g: string) => g.slice(7) + g.slice(0, 7)) }];
    const B = SeqSets.build(rotated).set;
    const Q = "(onto fam/0 fam)";
    const fresh = run(B, Q);
    const spy = vi.spyOn(Mappings, "compose").mockImplementationOnce(() => {
      throw new Error("injected");
    });
    expect(() => run(A, "(onto fam gene)")).toThrow("injected");
    spy.mockRestore();
    expect(run(B, Q)).toBe(fresh);
  });
});
