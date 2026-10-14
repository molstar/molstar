/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The findings IN2 closes, through the dispatcher and the default behaviors (VI1's list): a hover per fingerprint band,
// stripe and split row; the wheel with Shift; a double click on a group header; each behavior's teardown; a select drag
// closed, then widened, on release; a long drag closed by runs.

import { OrderedSet, type SeqLoci } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { case3, case7, case13, case14 } from "../../core/_test/format/cases.ts";
import { EntryKind, type LayerRenderer } from "../src/index.ts";
import type { Viewports } from "../src/viewport.ts";
import { driven } from "./input-fixtures.ts";
import { preparedRenderers, renderers } from "./layout-fixtures.ts";

const S = "family.main";

/** The stub renderers, the `points` one with a hit of its own: bands of 5 px down, stripes of 3 px across. */
function withBands(): LayerRenderer[] {
  return renderers().map((r) =>
    r.name === "points"
      ? {
          ...r,
          hit: (_v, _s, _p, out) => {
            out.state = Math.floor(out.point[1] / 5) % 4;
            out.part = Math.floor(out.point[0] / 3) % 3;
            return true;
          },
        }
      : r,
  );
}

describe("the hover element is every field of the hit", () => {
  it("a move between a fingerprint's ligand bands (state) or stripes (part) is a new element; one inside a band is not", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] }, renderers: withBands() });
    // Band edges every 5 px down and stripe edges every 3 px across, inside the entry and the fourth column's cell.
    const y0 = 5 * Math.ceil((t.yOf(S, "family.flvr-arginine") - 9) / 5);
    const x0 = 3 * Math.round(t.xOf(S, 0, 4) / 3);
    let hovers = 0;
    t.ctx.interaction.hover.subscribe(() => hovers++);
    t.pointer("move", S, x0 + 0.25, y0 + 0.25);
    const first = hovers;
    t.pointer("move", S, x0 + 0.5, y0 + 1.25); // the same band and stripe
    expect(hovers).toBe(first);
    t.pointer("move", S, x0 + 0.5, y0 + 5.25); // the next band
    expect(hovers).toBe(first + 1);
    t.pointer("move", S, x0 + 3.5, y0 + 5.25); // the next stripe
    expect(hovers).toBe(first + 2);
    const h = t.ctx.hover.live.current;
    expect([h.state, h.part]).toEqual([Math.floor((y0 + 5.25) / 5) % 4, Math.floor((x0 + 3.5) / 3) % 3]);
  });

  it("a move between split rows of one cell column is a new element", async () => {
    const F = "features.main";
    const t = await driven(case13(), { size: { [F]: [1000, 900] } });
    const NV = "features.natural-variants";
    const x = t.xOf(F, 0, 174);
    const seen: number[] = [];
    t.ctx.interaction.hover.subscribe((e) => seen.push(e.hit.split));
    for (const row of [3, 4, 4, 5]) t.pointer("move", F, x, t.yOf(F, NV, row));
    expect(seen.slice(-3)).toEqual([3, 4, 5]);
  });
});

describe("Shift turns a vertical wheel into a pan, by whichever delta is set", () => {
  it("deltaY alone, deltaX alone (macOS), and a horizontal wheel without Shift", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const vp = t.vp(S);
    vp.zoomTo(0, 5, 15);
    const s = vp.live.current;
    const y = t.layout(S).scrolledTop + 9;
    const at = s.scrollX;
    expect(t.wheel(S, 500, y, 0, 20, { shift: true })).toBe(true);
    expect(s.scrollX - at).toBe(20);
    expect(t.wheel(S, 500, y, 25, 0, { shift: true })).toBe(true);
    expect(s.scrollX - at).toBe(45);
    expect(t.wheel(S, 500, y, 10, 2)).toBe(true);
    expect(s.scrollX - at).toBe(55);
    // Vertical without Shift scrolls rows, not columns.
    const sx = s.scrollX;
    t.wheel(S, 500, y, 0, 20);
    expect(s.scrollX).toBe(sx);
  });
});

describe("a group's header toggles on a single click only", () => {
  it("a double click toggles once", async () => {
    const t = await driven(case7());
    const A = "complex-one.main";
    const l = t.layout(A);
    const e = t.entries(A).findIndex(([kind]) => kind === EntryKind.group);
    const y = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2 - t.vp(A).live.current.scrollY;
    const ref = "complex-one.structure-1brs";
    const was = t.ctx.ui.expanded.value.has(ref);
    t.click(A, 30, y);
    t.advance(50);
    expect(t.ctx.ui.expanded.value.has(ref)).toBe(!was);
    t.click(A, 30, y); // the second click of a double click
    expect(t.ctx.ui.expanded.value.has(ref)).toBe(!was);
    t.advance(1000);
    t.click(A, 30, y); // a new single click
    expect(t.ctx.ui.expanded.value.has(ref)).toBe(was);
  });
});

describe("a split track's chevron after its label toggles its rows", () => {
  it("a click on the chevron toggles; one on the label's text or well past the chevron does not", async () => {
    const A = "features.main";
    const NV = "features.natural-variants";
    const t = await driven(case13(), { size: { [A]: [1000, 900] }, renderers: preparedRenderers() });
    const m = t.ctx.graph.assembled.value?.view;
    const l = t.layout(A);
    if (!m) throw new Error("no view");
    const e = t.entries(A).findIndex(([kind, name, row]) => kind === EntryKind.track && name === NV && row === -2);
    const y = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2 - t.vp(A).live.current.scrollY;
    let depth = 0;
    for (let g = l.entryGroup[e]; g >= 0; g = m.resolved.groups.parent[g]) depth++;
    const font = t.ctx.text.font("ui", 12, 400);
    const label = t.ctx.text.label("Natural variants", font).width;
    const chevron = 10 + depth * 14 + label + 8;
    const open = () => t.ctx.ui.expanded.value.has(NV);
    const was = open();
    t.click(A, 12 + depth * 14, y);
    t.advance(1000);
    t.click(A, chevron + 12, y);
    t.advance(1000);
    expect(open()).toBe(was);
    t.click(A, chevron, y);
    expect(open()).toBe(!was);
  });
});

describe("removing a behavior mid-interaction leaves nothing lit or stuck", () => {
  it("hover clears its owner; select its pending range; cursor its cursors and lit window; hoverLabel its label", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const y = t.layout(S).scrolledTop + 9;
    t.vp(S).zoomTo(0, 5, 15);
    t.pointer("move", S, 500, y);
    expect([t.ctx.hover.live.current.owner, t.ctx.ui.cursor(S).value, t.ctx.labels.current.value === null]).toEqual([
      "pointer",
      "grab",
      false,
    ]);
    t.ctx.behaviors.get("hoverLabel")?.remove();
    expect(t.ctx.labels.current.value).toBeNull();
    t.ctx.behaviors.get("cursor")?.remove();
    expect(t.ctx.ui.cursor(S).value).toBe("default");
    t.ctx.behaviors.get("hover")?.remove();
    expect(t.ctx.hover.live.current.owner).toBe("");

    // A select drag in progress, and a lit window.
    const u = await driven(case3(), { size: { [S]: [1000, 150] } });
    const pending = u.ctx.selection.pending;
    u.pointer("down", S, 500, u.layout(S).scrolledTop + 9, { alt: true });
    u.pointer("move", S, 480, u.layout(S).scrolledTop + 9, { alt: true });
    expect(pending.current.active).toBe(true);
    const v = pending.version.value;
    u.ctx.behaviors.get("select")?.remove();
    expect([pending.current.active, pending.version.value > v]).toEqual([false, true]);
    u.pointer("up", S, 480, u.layout(S).scrolledTop + 9, { alt: true });

    const w = await driven(case3(), { size: { [S]: [1000, 150] } });
    const g = w.frame(S, 0);
    w.pointer("move", S, (g[8] + g[9]) / 2, 10);
    expect(w.ctx.ui.window.current.kind).toBe(0);
    w.ctx.behaviors.get("cursor")?.remove();
    expect([w.ctx.ui.window.current.kind, w.ctx.ui.window.current.part]).toEqual([-1, 0]);
  });
});

describe("a select drag is closed, then widened, on release", () => {
  it("a drag over nucleotides snaps to whole codons on release when the view holds their translation", async () => {
    const G = "gene-view.main";
    const t = await driven(case14());
    const dna = t.yOf(G, "gene-view.dna");
    // Alt-drag on data selects a range; the nucleotides 13 to 15 (inside the first codon, 12 to 14).
    t.vp(G).zoomTo(0, 8, 30);
    t.drag(G, t.xOf(G, 0, 13), dna, t.xOf(G, 0, 15), dna, { alt: true });
    const seq = t.ctx.selection.current.value.seq as SeqLoci;
    const set = t.ctx.graph.assembled.value?.set;
    const byUnit = new Map<string, number[]>();
    for (const el of seq.elements) {
      const xs: number[] = [];
      OrderedSet.forEach(el.indices, (i: number) => void xs.push(i), null);
      byUnit.set(set?.strings.get(set.units.label[el.unit]) ?? "", xs);
    }
    // The three nucleotides dragged and the two residues that contain them: the codons light whole (marks project a
    // residue onto its nucleotides), with the drag's nucleotides as dragged.
    expect([...byUnit]).toEqual([
      ["Gene fragment", [13, 14, 15]],
      ["Protein", [0, 1]],
    ]);
  });

  it("a drag over a million positions closes by runs: one range, no array of positions", async () => {
    const { createSeqViewSpec } = await import("@molstar/seq-core/builder");
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "chr", alphabet: "dna", value: "A".repeat(1_000_000) });
    b.view({ ref: "v", label: "V" }).track({ label: "DNA" }).letters();
    const t = await driven(b.getDocument());
    const set = t.ctx.graph.assembled.value?.set;
    const projector = (t.ctx.viewports as Viewports).projector;
    const loci = {
      kind: "seq-loci",
      set,
      elements: [{ unit: 0, indices: OrderedSet.ofRange(10, 999_990) }],
    } as unknown as SeqLoci;
    const closed = projector.closeLoci(loci, "containing") as SeqLoci;
    expect(closed.elements).toHaveLength(1);
    const indices = closed.elements[0].indices;
    expect([ArrayBuffer.isView(indices), OrderedSet.size(indices), OrderedSet.min(indices)]).toEqual([
      false,
      999_980,
      10,
    ]);
  });
});
