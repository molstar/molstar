/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guard for the label layer (NG-12), on RB0b's stepped `heapScaling`: `layLabels` over a few hundred labels at an
// unchanged view (a manager's attribute storm on each hover lays them again and writes nothing). A file of its own
// (LA-28: guards sharing a file share V8's feedback). The positive control is the template's: the same step with a
// leak.

import { describe, expect, it } from "vitest";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "../../core/_test/model/heap.js";
import { layLabels } from "../src/labels.js";
import { labels, layout, measure, N, view } from "./labels-heap-fixture.js";

describe("the label layer allocates nothing per frame", () => {
  it("layLabels at an unchanged view: laid, compared, nothing written", () => {
    layLabels(labels, view, layout, measure);
    let writes = 0;
    const step = () => {
      writes += layLabels(labels, view, layout, measure);
    };
    const r = heapScaling(N / 10, step);
    expect(writes).toBe(0);
    expect(isFlat(r), describeScaling(r, N / 10)).toBe(true);
    expect(controlCaught(() => heapScaling(N / 10, leaking(step)))).toBe(true);
  });
});
