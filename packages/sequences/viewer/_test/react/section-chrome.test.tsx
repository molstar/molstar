/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's chrome in the React UI: the header, its switches and its hidden chips read one value, computed once per
// change of its inputs (the DOM's twin is _test/dom/subscribe.test.ts).

import { act } from "react";
import { describe, expect, it, vi } from "vitest";
import { shownSections } from "../../src/index.ts";
import { harness, mount, openedContext, settle, stubSurface } from "./fixtures.tsx";

const { cleanups } = harness();

describe("a section's chrome", () => {
  it("is computed once per change for the header, the switches and the chips together (React)", async () => {
    const { ctx } = await openedContext();
    const sections = shownSections(ctx).length;
    expect(sections).toBeGreaterThan(0);
    // `sectionInfo` asks the registry for the presets once per section it computes.
    const list = vi.spyOn(ctx.registry, "list");
    const computes = () => list.mock.calls.filter((c) => c[0] === "presets").length;
    const m = await mount(ctx, { surface: stubSurface });
    cleanups.push(() => ctx.dispose(), m.unmount);
    expect(m.container.querySelector(".seqv-section-header")).not.toBeNull();
    expect(m.container.querySelector(".seqv-switches")).not.toBeNull();
    expect(computes()).toBe(sections);
    await act(async () => {
      expect(ctx.actions.run("section.style", { section: "chains.main", value: "railed" })).toBe(true);
    });
    await settle();
    expect(computes()).toBe(2 * sections);
  });
});
