/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A renderer that throws in any of its parts is omitted with one `seqv.renderer.error` naming its layer, and the frame
// draws the rest (seqspec:viewer#extension-points): its draw, its rows band (dense rows), its header, its gutter and
// its overlay, each caught where the composer calls it (LC-41: the proof's corpus reached none of these catches).

import { describe, expect, it } from "vitest";
import { Dirty } from "../src/draw/index.js";
import { KernelRenderers, type LayerRenderer } from "../src/index.js";
import { composed, recording } from "./composer-fixtures.js";
import { FAMILY_SECTION, familyWithTree } from "./families.js";

const ROWS = KernelRenderers.find((r) => r.name === "rows") as LayerRenderer;
const boom = () => {
  throw new Error("boom");
};
const PARTS: Record<string, (r: LayerRenderer) => LayerRenderer> = {
  draw: (r) => ({ ...r, draw: boom }),
  rows: (r) => ({ ...r, rows: boom }),
  header: (r) => ({ ...r, header: r.header && { ...r.header, draw: boom } }) as LayerRenderer,
  overlay: (r) => ({ ...r, overlay: boom }),
};

describe("a renderer that throws", () => {
  for (const [part, mutant] of Object.entries(PARTS))
    it(`in its ${part}: omitted with one issue, the frame drawn`, async () => {
      // Rows under `denseRowHeight` are one band (`rows`); else a row each (`draw`).
      const dense = part === "rows";
      const t = await composed(familyWithTree(40, 10), {
        section: FAMILY_SECTION,
        replace: { rows: mutant(ROWS) },
        render: dense ? { denseRowHeight: 100 } : undefined,
        width: 1000,
        height: 400,
      });
      const said: string[] = [];
      const sub = t.ctx.events.issue.subscribe((i) => {
        if (i.code === "seqv.renderer.error") said.push(i.message);
      });
      const be = t.frame(recording(), Dirty.all);
      t.frame(recording(), Dirty.all);
      sub.unsubscribe();
      expect(said.length, part).toBe(1);
      expect(said[0]).toMatch(/threw while drawing and is omitted: boom/);
      expect(be.opCount).toBeGreaterThan(0);
      t.ctx.dispose();
    });
});
