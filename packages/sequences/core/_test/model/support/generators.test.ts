/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
import { alignment, gene, table } from "./generators.ts";

describe("generators", () => {
  it("builds a 20k-row alignment with two bound rows, the same on every call", () => {
    const a = alignment();
    expect(a.rows).toHaveLength(20_000);
    expect(a.rows.every((r) => r.length === a.length)).toBe(true);
    expect(a.bound.map((b) => b.row)).toEqual([0, 10_000]);
    expect(a.bound[0].sequence).toBe(a.rows[0].replaceAll("-", ""));
    expect(alignment().rows[777]).toBe(a.rows[777]);
  });

  it("builds a DMD-sized gene whose exons join to a CDS that codes the protein", () => {
    const g = gene();
    expect(g.exons).toHaveLength(79);
    expect(g.dna.length).toBeGreaterThan(2_200_000);
    expect(g.protein).toHaveLength(3685);
    expect(g.exons.map(([s, e]) => g.dna.slice(s, e)).join("")).toBe(g.cds);
    expect(g.cds.length).toBe(3685 * 3 + 3);
    expect(g.cds.startsWith("ATG")).toBe(true);
    expect(g.dna.slice(g.exons[1][0] - 2, g.exons[1][0])).toBe("AG");
  });

  it("builds a 3M-row table of typed columns, ordered by start", () => {
    const t = table();
    expect(t.rows).toBe(3_000_000);
    expect(t.start.length).toBe(t.rows);
    expect(t.start[t.rows - 1]).toBeGreaterThan(t.start[0]);
    expect(Math.max(...t.category.subarray(0, 1000))).toBeLessThan(t.categories.length);
  });
});
