/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyed-sets story's five modes (`Seq*/Views and layouts/Alignments › Order and filter rows`): each mode's patch, applied to the story's document, sets the rows
// layer's order and filter, and the lanes, the layout's rows and the overview's extent follow. (The story's first
// version patched the track's ref, which no layer has: nothing changed, and the buttons looked dead.)

import { describe, expect, it } from "vitest";
import { pf00042 } from "../../core/_test/model/alignment-fixtures.js";
import { identityOf } from "../../core/src/algorithms/index.js";
import { keyedFamilyDoc, MEMBERS_LAYER, MODES, modePatch } from "./keyed-family.js";
import { opened } from "./layout-fixtures.js";

const FAMILY = pf00042();
const IDENTITY = FAMILY.map((r) => {
  const a = FAMILY[0].gapped;
  let both = 0;
  let same = 0;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== "-" && r.gapped[i] !== "-") {
      both++;
      if (a[i] === r.gapped[i]) same++;
    }
  }
  return same / both;
});
const mode = (key: string) => MODES.find((m) => m.key === key) ?? MODES[0];

async function view() {
  const t = await opened(keyedFamilyDoc(FAMILY), { size: { "v.main": [900, 500] } });
  const lanes = () => {
    const v = t.ctx.graph.assembled.value?.view?.resolved;
    return v ? Array.from(v.lanes.row) : [];
  };
  const displayed = () => {
    const r = t.ctx.layout.visibleRows("v.main")?.current;
    return r ? Array.from(r.dataRow.subarray(0, r.count)) : [];
  };
  const apply = async (key: string) => {
    const before = t.ctx.issues.current.value.length;
    await t.ctx.applyPatch(modePatch(mode(key)) as never);
    expect(t.ctx.issues.current.value.filter((i) => i.severity === "error")).toEqual([]);
    void before;
  };
  return { t, lanes, displayed, apply };
}

describe("the keyed-sets story's modes", () => {
  it("the layer the patches name exists, and the document opens in data order: 32 rows", async () => {
    const { t, lanes } = await view();
    expect(t.ctx.graph.assembled.value?.set.strings.find(MEMBERS_LAYER)).toBeGreaterThanOrEqual(0);
    expect(lanes()).toEqual(Array.from({ length: 32 }, (_, i) => i));
  });

  it("each mode's patch sets the lanes and the layout's rows; the overview's extent follows the kept rows", async () => {
    const { t, lanes, apply } = await view();
    const extent = () => {
      const l = t.layout("v.main");
      return [l.rows.count, l.scrolledHeight];
    };
    // identity to the reference row, descending (ties in data order)
    await apply("identity");
    const byId = lanes();
    expect(byId).toHaveLength(32);
    expect(byId[0]).toBe(0);
    expect(byId.map((r) => IDENTITY[r])).toEqual([...byId.map((r) => IDENTITY[r])].sort((a, b) => b - a));
    expect(extent()).toEqual([32, 32 * 18]);
    // clustered at 90%, then by identity: every row, the cluster's rows together (here only two pairs merge)
    await apply("clusters");
    expect(lanes()).toHaveLength(32);
    expect(new Set(lanes()).size).toBe(32);
    // representatives only
    await apply("reps90");
    expect(lanes()).toHaveLength(30);
    expect(extent()).toEqual([30, 30 * 18]);
    expect(lanes().map((r) => IDENTITY[r])).toEqual([...lanes().map((r) => IDENTITY[r])].sort((a, b) => b - a));
    await apply("reps50");
    expect(lanes()).toHaveLength(22);
    expect(extent()).toEqual([22, 22 * 18]);
    // back to data order
    await apply("data");
    expect(lanes()).toEqual(Array.from({ length: 32 }, (_, i) => i));
    expect(extent()).toEqual([32, 32 * 18]);
  });

  it("the displayed rows are the lanes' (the layout's dataRow), and the clusters' rows are adjacent", async () => {
    const { t, displayed, apply } = await view();
    await apply("reps90");
    const rows = displayed();
    expect(rows.length).toBeGreaterThan(10);
    expect(rows.every((r) => IDENTITY[r] !== undefined)).toBe(true);
    // A non-representative at 90% (identity 0.94 with HBA_HUMAN, row 12) is not shown.
    expect(identityOf).toBeDefined();
    expect(t.layout("v.main").rows.count).toBe(30);
    expect(rows).not.toContain(12);
  });

  it("a patch that names the track instead of its layer changes nothing (what the first version of the story did)", async () => {
    const { t, lanes } = await view();
    const m = mode("reps90");
    const r = await t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: "members", params: { order: m.order, filter: m.filter } }],
    } as never);
    expect(lanes()).toHaveLength(32);
    expect(JSON.stringify(r)).toContain("unknown param");
  });
});
