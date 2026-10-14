/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tree beside the rows (LY1), headless: laid out by branch length for the layout's display rows, its nodes hit with
// a forgiving hitbox through the rows renderer's `header.hit`, a clade's span, a collapsed clade as one wedge row that
// hits its node. A stub `rows` renderer prepares what the tree reads (`TreeRows`); the hover, label, cursor and clicks
// of the old test are IN2's and IN3's.

import { describe, expect, it } from "vitest";
import { Hit, HitRegion, HoverKind, type LayerRenderer } from "../src/index.js";
import {
  hitTree,
  layoutTree,
  TREE_INSET,
  TREE_LEFT,
  TREE_RIGHT,
  TreeNode,
  type TreeRows,
  treeOf,
  treeSpan,
  treeX,
} from "../src/tree.js";
import { FAMILY_SECTION as F, familyWithTree } from "./families.js";
import { opened, renderers } from "./layout-fixtures.js";

/** A rows renderer that prepares the tree and its rows' order, with the tree's header. */
function rowsWithTree(): LayerRenderer<TreeRows> {
  return {
    name: "rows",
    geometries: ["linear"],
    prepare: (input) => {
      const a = input.section.blocks[0]?.alignment ?? -1;
      const aln = input.set.alignments[a];
      const track = input.section.tracks[input.layer.track];
      const rowOf = Int32Array.from(
        { length: track.laneCount },
        (_, k) => input.section.lanes.row[track.laneStart + k],
      );
      return {
        set: input.set,
        alignment: a,
        rowOf,
        order: null,
        source: null,
        tree: aln?.tree ? treeOf(aln.tree) : null,
      };
    },
    header: { inset: () => TREE_INSET, hit: (v, s, _p, out) => hitTree(v, s, out), span: treeSpan },
  };
}

async function family() {
  const t = await opened(familyWithTree(40, 10), {
    renderers: [...renderers().filter((r) => r.name !== "rows"), rowsWithTree() as LayerRenderer],
    size: { [F]: [1000, 300] },
  });
  const m = t.ctx.graph.assembled.value?.view;
  const k = m ? Array.from(m.resolved.layers.type).findIndex((x) => x >= 0 && m.set.strings.get(x) === "rows") : -1;
  const s = m?.prepared[k]?.state as TreeRows;
  if (!s?.tree) throw new Error("no tree");
  return { ...t, s, g: s.tree };
}

describe("the tree beside the rows", () => {
  it("lays out by branch length, hits its nodes, collapses a clade into a wedge row that spans it", async () => {
    const { s, g, ...t } = await family();
    const l = t.layout(F);
    layoutTree(s, l.rows);
    const hit = Hit.create();
    // A clade of seed 0 (the root's first child's first child): its branch on its display row.
    const node = 2;
    const top = l.scrolledTop + l.entryTop[l.rows.entry];
    const width = TREE_INSET - TREE_LEFT - TREE_RIGHT;
    t.ctx.layout.hitTest(F, treeX(g, node, TREE_LEFT, width) - 1, top + g.y[node] * l.rows.height, hit);
    expect([hit.region, hit.kind, hit.treeNode]).toEqual([HitRegion.header, HoverKind.header, node]);
    const span = new Int32Array(2);
    treeSpan(s, node, span);
    expect(span[1] - span[0] + 1).toBe(10);
    // Folded: one wedge row for the clade's 10 rows; the wedge spans its own row.
    t.ctx.ui.toggleClade("family", node);
    const l2 = t.layout(F);
    expect(l2.rows.count).toBe(31);
    t.ctx.layout.hitTest(F, 300, top + (span[0] + 0.5) * l2.rows.height, hit);
    expect(hit.treeNode).toBe(node);
    layoutTree(s, l2.rows);
    treeSpan(s, node, span);
    expect([span[0], span[1]]).toEqual([span[0], span[0]]);
    expect([g.state[node], g.wedgeText.get(node)]).toEqual([TreeNode.wedge, expect.stringMatching(/ \+9$/)]);
    t.ctx.ui.toggleClade("family", node);
    expect(t.layout(F).rows.count).toBe(40);
  });

  it("a forgiving hitbox: the nearest node within 6 px; a leaf's branch is its row; none past the inset", async () => {
    const { s, g, ...t } = await family();
    const l = t.layout(F);
    layoutTree(s, l.rows);
    const hit = Hit.create();
    const h = l.rows.height;
    const top = l.scrolledTop + l.entryTop[l.rows.entry];
    const w = TREE_INSET - TREE_LEFT - TREE_RIGHT;
    const X = (k: number) => treeX(g, k, TREE_LEFT, w);
    const Y = (rows: number) => top + rows * h;
    const at = (x: number, y: number) => {
      t.ctx.layout.hitTest(F, x, y, hit);
      return hit.treeNode;
    };
    // Node 2 (seed 0's clade, 10 rows): its connector, 4 px to its right, midway down a row away from its children.
    const node = 2;
    const mid = (g.lo[node] + g.hi[node]) / 2;
    expect(at(X(node) + 4, Y(mid))).toBe(node);
    expect(at(X(node) + 7, Y(mid))).toBe(-1);
    // An inner node on a row boundary, reached from the row above, 4 px off its branch.
    const inner = g.byY.find((k) => g.t.leafRow[k] < 0 && g.y[k] === Math.round(g.y[k]) && g.t.parent[k] >= 0) ?? -1;
    expect(inner).toBeGreaterThan(0);
    const xb = (X(g.t.parent[inner]) + X(inner)) / 2;
    expect(at(xb, Y(g.y[inner]) - 4)).toBe(inner);
    expect(hit.displayRow).toBe(g.y[inner] - 1);
    // A leaf's branch is its row: no node; past the tree's inset + 8 px, none either.
    const leaf = g.byY.find((k) => g.t.leafRow[k] >= 0) ?? -1;
    expect(at(X(leaf) - 1, Y(g.y[leaf]))).toBe(-1);
    expect(at(TREE_INSET + 9, Y(mid))).toBe(-1);
    expect(hit.kind).toBe(HoverKind.header);
    const span = new Int32Array(2);
    treeSpan(s, node, span);
    expect([span[0], span[1]]).toEqual([g.first[node], g.last[node]]);
  });

  it("a resize and a row zoom lay the tree out once; a fold lays it out again", async () => {
    const { s, g, ...t } = await family();
    layoutTree(s, t.layout(F).rows);
    const at = g.at;
    t.vp(F).setSize(900, 260, 1);
    t.vp(F).setRowScale(0.6, 0);
    layoutTree(s, t.layout(F).rows);
    expect(g.at).toBe(at);
    t.ctx.ui.toggleClade("family", 2);
    layoutTree(s, t.layout(F).rows);
    expect(g.at).toBe(at + 1);
  });
});
