/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyboard cursor's ring and when it shows (seqspec:viewer#behaviors): the ring sits evenly around the cell as
// drawn (its gap on the right and the lane's padding left out), Escape puts the cursor away with the selection, and a
// press puts it away and gives a focus that shows none until a key.

import { describe, expect, it } from "vitest";
import { Dirty, Pass, SectionComposer, SkinToken } from "../src/draw/index.js";
import { KernelRenderers, type LayerRenderer } from "../src/index.js";
import { recording } from "./composer-fixtures.js";
import { driven } from "./input-fixtures.js";
import { passRects, twoBlocks } from "./overlay-fixtures.js";

async function opened() {
  const t = await driven(twoBlocks(), { renderers: KernelRenderers as LayerRenderer[], size: { v: [800, 300] } });
  const section = t.ctx.viewports.sections.value[0];
  t.ctx.viewports.get(section)?.setView(30, 0, "api");
  return { t, section };
}

describe("the keyboard cursor", () => {
  it("its ring's inner edge is the cell as drawn, on all four sides", async () => {
    const { t, section } = await opened();
    const c = new SectionComposer(t.ctx, section);
    t.focus(section);
    t.key(section, "ArrowRight");
    const rec = recording();
    c.frame(rec, Dirty.all);
    const [ring] = passRects(rec, Pass.overlay).filter((r) => r[4] === SkinToken.focusRing);
    const [x, y, w, h, , , stroke] = ring;
    const inner = [x + stroke, y + stroke, w - 2 * stroke, h - 2 * stroke];
    const cell = passRects(rec, Pass.main).find((r) => r[0] === inner[0] && r[2] > 5 && r[3] > 5);
    expect(cell?.slice(0, 4)).toEqual(inner);
  });

  it("Escape puts it away with the selection; the next key brings it back where it was", async () => {
    const { t, section } = await opened();
    t.focus(section);
    t.key(section, "ArrowRight");
    t.key(section, "ArrowRight");
    expect(t.ctx.ui.keyCursor.value?.index).toBe(2);
    t.key(section, "Escape");
    expect(t.ctx.ui.keyCursor.value).toBeNull();
    t.key(section, "ArrowRight");
    expect(t.ctx.ui.keyCursor.value?.index).toBe(3);
  });

  it("a press puts it away, and the focus the press gives shows none until a key", async () => {
    const { t, section } = await opened();
    t.focus(section);
    t.key(section, "ArrowRight");
    expect(t.ctx.ui.keyCursor.value).not.toBeNull();
    t.click(section, 300, 40);
    expect(t.ctx.ui.keyCursor.value).toBeNull();
    // Focus lost, then a press that focuses: no cursor; a key shows it; a Tab back in shows it at once.
    t.focus(section, false);
    t.click(section, 300, 40);
    t.focus(section);
    expect(t.ctx.ui.keyCursor.value).toBeNull();
    t.key(section, "ArrowRight");
    expect(t.ctx.ui.keyCursor.value).not.toBeNull();
    t.focus(section, false);
    t.focus(section);
    expect(t.ctx.ui.keyCursor.value).not.toBeNull();
  });
});
