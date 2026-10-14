/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Nightingale has no two-axis tracks (seqspec:libraries#nightingale): case 8's dot plot (its predicted aligned error and
// the contacts over it) is left out and reported, its fallback `omit`; everything else of the view is materialized.

import { compile, resolveView } from "@molstar/seq-core/format";
import { expect, it } from "vitest";
import { case8 } from "../../core/_test/format/cases.js";
import { fixtureLoaders, registries } from "../../core/_test/format/load-fixtures.js";
import { materialize, NightingaleSupported } from "../src/index.js";

it("leaves a dot plot out and says so", async () => {
  const reg = registries();
  const doc = case8();
  const { set } = await compile(doc, reg, fixtureLoaders());
  const { view, issues } = resolveView(doc, set, "chains", reg, NightingaleSupported);
  const omitted = issues.filter((i) => i.code === "svs.render.omitted");
  expect(omitted.length).toBe(2);
  expect(omitted.every((i) => i.message.includes("dot_plot"))).toBe(true);
  expect(materialize(view, set).root).toBeTruthy();
});
