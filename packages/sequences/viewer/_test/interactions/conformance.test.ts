/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// XR1's conformance kit over the fingerprint (seqspec:viewer#extension-points): case 6 drawn at three zooms and below
// the minimum on the recording backend, every drawn rect hitting back to its layer (`items: false`: the hit gives a row
// and a kind, not an item), the legend listing the kinds in the vocabulary's order, and a context with only the format
// half drawing the `matrix` fallback with one `svs.render.fallback` for the layer; wrapped is not declared (the type's
// layouts are `tracks`), nor is circular (on rings it is drawn as bars, with one warning: CG1b). The kit's tee has its
// SVG side (SV1a): the recording and the SVG backend see the same calls.

import { describe, expect, it } from "vitest";
import { case6 } from "../../../core/_test/format/cases.js";
import { KernelRenderers, type LayerRenderer } from "../../src/index.js";
import { InteractionsRegistry } from "../../src/interactions/index.js";
import { FingerprintRenderer, InteractionsViewerRegistry } from "../../src/interactions/viewer/index.js";
import { DefaultPresets } from "../../src/presets.js";
import { type ConformanceFixtures, rendererConformance } from "../../src/testing/index.js";
import { asText, directed, ringed } from "../conformance-fixtures.js";
import { opened } from "../layout-fixtures.js";
import { ShortFingerprint } from "./scenes.js";

/** The built-in `circular` preset beside the fixture presets: sections on rings. */
const RINGS = { presets: DefaultPresets.filter((p) => p.name === "circular") };
const fixtures = (
  doc: () => ReturnType<typeof case6>,
  view: string,
  legend?: string[],
  reversed = false,
): ConformanceFixtures => ({
  reversed,
  open: async ({ geometry, entry, direction }) => {
    const d = direction ? directed(doc(), direction) : doc();
    const as = geometry === "wrapped" ? asText(doc()) : geometry === "circular" ? ringed(d) : d;
    const registry = entry ? [InteractionsViewerRegistry, ShortFingerprint] : [{ svs: InteractionsRegistry }];
    return (await opened(as, { renderers: [...KernelRenderers], view, registry: [...registry, RINGS] })).ctx;
  },
  legend,
  items: false,
});
const failures = (checks: Awaited<ReturnType<typeof rendererConformance>>) =>
  checks.filter((c) => !c.ok).map((c) => `${c.check} ${c.geometry}: ${c.detail}`);

describe("the conformance kit", () => {
  it("case 6: batches and hits at three zooms (tiny too), the legend, the fallback without the entry", async () => {
    const checks = await rendererConformance(
      FingerprintRenderer as LayerRenderer,
      fixtures(case6, "poses-view", ["hydrogen bond", "π stacking", "cation–π", "hydrophobic contact"]),
    );
    expect(failures(checks)).toEqual([]);
    expect(checks.filter((c) => c.check === "batches")).toHaveLength(6);
    expect(checks.filter((c) => c.check === "hits" && c.geometry === "linear")).toHaveLength(3);
    // Declared linear only: wrapped is checked for its fallback alone, one issue for the layer; on rings, bars.
    expect(checks.filter((c) => c.geometry === "wrapped").map((c) => [c.check, c.detail])).toEqual([
      ["fallback", "1 layer(s), 1 issue(s) [svs.render.fallback]"],
    ]);
    expect(checks.find((c) => c.geometry === "circular")?.detail).toMatch(/^1 layer\(s\), 1 issue\(s\), \d+ sectors$/);
    expect(checks.find((c) => c.check === "fallback" && c.geometry === "linear")?.detail).toBe(
      "1 layer(s), 1 issue(s) [svs.render.fallback]",
    );
  });

  it("case 6 reversed (AX1): every check on reversed pieces, and the reversed draw the forward one mirrored", async () => {
    const legend = ["hydrogen bond", "π stacking", "cation–π", "hydrophobic contact"];
    const entry = FingerprintRenderer as LayerRenderer;
    const checks = [
      ...(await rendererConformance(
        entry,
        fixtures(() => directed(case6(), -1), "poses-view", legend),
      )),
      ...(await rendererConformance(entry, fixtures(case6, "poses-view", legend, true))),
    ];
    expect(failures(checks)).toEqual([]);
    expect(checks.filter((c) => c.check === "mirror").length).toBeGreaterThan(0);
  });
});
