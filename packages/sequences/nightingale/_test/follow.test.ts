/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `follow` in Nightingale (seqspec:libraries#nightingale, LK1): managers have no shared viewport, so a section that
// follows another materializes as it would without: the same elements, the same issues, nothing reported for it.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import {
  compile,
  createRegistries,
  KernelRegistry,
  MvsRegistry,
  register,
  resolveView,
} from "@molstar/seq-core/format";
import { expect, it } from "vitest";
import { materialize, NightingaleSupported } from "../src/index.js";

const kernel = createRegistries();
register(kernel, [KernelRegistry, MvsRegistry]);
const none = { load: async () => Promise.reject(new Error("no files")) };

async function draw(follow: boolean) {
  const b = createSeqViewSpec();
  const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: "GAATTCGCCACCATGGCTTCCAAGGAACTGACC" });
  const protein = translate(b.data, { ref: "protein", from: gene, cds: [[12, 33]] });
  const v = b.view({ ref: "v", label: "V" });
  const dna = v.section({ ref: "v.dna" }).block({ sequence: gene });
  dna.track({ label: "DNA" }).letters();
  v.section({ ref: "v.protein", ...(follow ? { follow: dna.handle, follow_anchor: null } : {}) })
    .block({ sequence: protein })
    .track({ label: "Protein" })
    .letters();
  const doc = b.getDocument();
  const { set, issues } = await compile(doc, kernel, none);
  expect(issues.filter((i) => i.severity === "error")).toEqual([]);
  const r = resolveView(doc, set, "v", kernel, NightingaleSupported);
  return { follows: r.view.follows.count, resolved: r.issues, ...materialize(r.view, set) };
}

it("a section that follows another materializes as one that does not", async () => {
  const [plain, linked] = [await draw(false), await draw(true)];
  expect([plain.follows, linked.follows]).toEqual([0, 1]);
  expect(linked.resolved).toEqual(plain.resolved);
  expect(linked.issues).toEqual(plain.issues);
  expect(linked.root).toEqual(plain.root);
});
