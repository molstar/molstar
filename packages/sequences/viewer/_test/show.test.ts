/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// "Show me" across views (QR2a): `showExtents` focuses in the current view when it shows the extents, else switches to
// the first view a block of which shows one of their systems; case 10's views hold the kinase chain (Chains, Strip) and
// the gatekeeper alignment (Alignment, Annotations and MSA) apart.

import { Extents, CoordinateSystemKind as K, type SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { case10 } from "../../core/_test/format/cases.ts";
import { showExtents } from "../src/index.ts";
import { viewShowing } from "../src/show.ts";
import { shown } from "./projection-fixtures.ts";

const extents = (set: SeqSet, kind: number, ref: string, start: number, stop: number) => {
  const e = Extents.create(2);
  const key = set.strings.find(ref) as number;
  const id =
    kind === K.columns ? (set.alignments.findIndex((a) => a?.key === key) as number) : set.units.key.indexOf(key);
  Extents.add(e, { kind, id, row: -1 }, start, stop);
  return e;
};

describe("showExtents", () => {
  it("switches to the view that shows the extents, focuses there, and back", async () => {
    const t = await shown(case10(), "chains");
    const set = t.ctx.graph.assembled.value?.set as SeqSet;
    const columns = extents(set, K.columns, "gatekeeper-msa", 3, 9);
    expect(viewShowing(t.ctx, columns)).toBe("alignment");
    expect(await showExtents(t.ctx, columns, { transitionMs: 0 })).toBe(true);
    expect(t.ctx.views.current.value).toBe("alignment");
    const [entry] = [...t.ctx.focus.current.value.values()];
    expect([entry.start, entry.stop]).toEqual([3, 9]);
    // the kinase chain is in Chains, the first view that shows it
    const chain = extents(set, K.unit, "1iep/A", 80, 90);
    expect(await showExtents(t.ctx, chain, { transitionMs: 0 })).toBe(true);
    expect(t.ctx.views.current.value).toBe("chains");
  });

  it("prefers a view that holds them to one the projection reaches them in", async () => {
    // Text shows the kinase domain, a unit of ABL1 that 1iep/A reaches through a mapping: Chains holds the chain itself
    const t = await shown(case10(), "text");
    const set = t.ctx.graph.assembled.value?.set as SeqSet;
    const chain = extents(set, K.unit, "1iep/A", 80, 90);
    expect(t.ctx.focus.fromExtents(chain, { transitionMs: 0 })).toBe(true);
    expect(await showExtents(t.ctx, chain, { transitionMs: 0 })).toBe(true);
    expect(t.ctx.views.current.value).toBe("chains");
  });

  it("stays in the current view when it shows them, and reports false when no view does", async () => {
    const t = await shown(case10(), "strip");
    const set = t.ctx.graph.assembled.value?.set as SeqSet;
    expect(await showExtents(t.ctx, extents(set, K.unit, "1iep/A", 80, 90), { transitionMs: 0 })).toBe(true);
    expect(t.ctx.views.current.value).toBe("strip");
    // a docked pose's ligand: no block shows it, no path reaches it
    const pose = extents(set, K.unit, "poses/0/ligand", 0, 1);
    expect([viewShowing(t.ctx, pose), await showExtents(t.ctx, pose)]).toEqual([null, false]);
    expect(t.ctx.views.current.value).toBe("strip");
  });
});
