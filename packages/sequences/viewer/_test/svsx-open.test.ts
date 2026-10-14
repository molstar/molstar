/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A viewer opens an `.svsx` archive (seqspec:sources#archives): the document's relative URIs are valid against the
// archive's base, which its loaders (`archiveLoaders`) carry, so `setDocument` and the patches over it validate as
// `compile` does. Before, `setDocument` validated without the base and refused every relative URI (`svs.source.uri`).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { bundleUris, writeSvsx } from "@molstar/seq-core/builder/svsx";
import { archiveLoaders, type Loaders, readSvsx } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { FASTA_URI, fixtureFetch, SITES_URI } from "../../core/_test/format/svsx-fixture.ts";
import { createSeqViewer, DefaultRegistry, type SeqIssue } from "../src/index.ts";
import { renderSvg } from "../src/svg/index.ts";
import { patch } from "./graph-fixtures.ts";
import { all, parseXml, rectsOf } from "./svg/xml.ts";

const nothing: Loaders = { load: (r) => Promise.reject(new Error(`only the archive's files (${r.uri})`)) };

/** The fixture's FASTA and sites by absolute URI: a letters track and a ranges track. */
function absolute() {
  const b = createSeqViewSpec();
  const src = b.data.source({ ref: "fa", uri: FASTA_URI, format: "fasta" });
  const seq = b.data.sequence({ ref: "abl", source_ref: src, label: "ABL kinase excerpt" });
  const sites = b.data.annotation({ ref: "sites", geometry: "ranges", target: seq, uri: SITES_URI, format: "tsv" });
  const block = b.view({ ref: "v", label: "ABL excerpt" }).block({ sequence: seq });
  block.track({ ref: "v.seq", label: "Sequence" }).letters();
  block.track({ ref: "v.sites", label: "Sites" }).ranges({ data: sites });
  return b.getDocument();
}

/** The document bundled and written by the builder, and read back with its loaders. */
async function archive() {
  const bundled = await bundleUris(absolute(), fixtureFetch().fetch);
  const svsx = await readSvsx(await writeSvsx(bundled.doc, bundled.files));
  if (!svsx.doc) throw new Error(svsx.issues.map((i) => i.message).join("; "));
  return { doc: svsx.doc, loaders: archiveLoaders(svsx.files, nothing, svsx.base) };
}

/** Over a fetch of the fixture files, the loaders of the document with absolute URIs. */
const direct: Loaders = {
  load: async (r) => new Uint8Array(await (await fixtureFetch().fetch(r.uri as string)).arrayBuffer()),
};

/** The filled rectangles of each track's data, by its clip (the sections' data areas in order). */
function cells(svg: string): number[] {
  const root = parseXml(svg);
  return all(root, "g")
    .filter((g) => g.attrs["clip-path"] && all(g, "path").some((p) => p.attrs.fill))
    .map((g) => all(g, "path").reduce((n, p) => n + rectsOf(p.attrs.d ?? "").length, 0));
}

describe("an .svsx archive in a viewer", () => {
  it("opens with no issue, its relative URIs answered from the archive, its data drawn as from the URIs", async () => {
    const { doc, loaders } = await archive();
    let issues: readonly SeqIssue[] = [];
    const svg = await renderSvg(doc, { width: 900, loaders, issues: (i) => (issues = i) });
    expect(issues).toEqual([]);
    expect(svg).toBe(await renderSvg(absolute(), { width: 900, loaders: direct }));
    // A cell per residue of the FASTA's 106, and the TSV's three sites.
    const [letters, sites] = cells(svg);
    expect([letters, sites]).toEqual([106, 3]);
  });

  it("a patch over the archive's document validates against its base too", async () => {
    const { doc, loaders } = await archive();
    const ctx = createSeqViewer({ registry: [DefaultRegistry], behaviors: [] }, { loaders });
    const opened = await ctx.setDocument(doc);
    expect([opened.ok, opened.issues]).toEqual([true, []]);
    const r = await ctx.applyPatch(patch([{ op: "set", ref: "v.sites", params: { label: "Kinase sites" } }]));
    expect([r.ok, r.issues]).toEqual([true, []]);
    ctx.dispose();
  });
});
