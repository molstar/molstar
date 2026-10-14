/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The layout and draw review's nits (2026-10-14), where a test can hold them: a path batch made with no room grows
// (VL-20); a view drawn without marks reads empty runs, not null (VL-22).

import { describe, expect, it } from "vitest";
import { DrawView, PathBatch, PathVerb } from "../src/draw/api.js";

describe("draw nits", () => {
  it("VL-20: `new PathBatch(0)` takes its verbs", () => {
    const p = new PathBatch(0);
    p.begin(0x112233);
    p.moveTo(1, 2);
    p.lineTo(3, 4);
    expect([p.n, p.nv, Array.from(p.coords.subarray(0, p.nc))]).toEqual([1, 2, [1, 2, 3, 4]]);
    expect(p.verbs[1]).toBe(PathVerb.line);
  });

  it("VL-22: a view's runs without marks are empty", () => {
    const r = new DrawView().marks.runs(1);
    expect([r.count, r.from(0), r.to(0)]).toEqual([0, 0, 0]);
  });
});
