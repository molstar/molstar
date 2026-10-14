/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ideogram navigator's story (`Seq*/Views and layouts/Navigators`) and the genome browser demo's navigator: hg38's committed
// cytobands (lib/seqstar/viewer/_test/fixtures/hg38-cytobands.json, UCSC's cytoBandIdeo, recorded by
// scripts/record-cytobands.ts) registered as the `ideogram` navigator, over the offline streamed genome of the format's
// stream-provider.ts (24 chromosomes, synthetic letters, genes and variants). Kept out of the story file so the smoke
// test can open the document. No React here.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Loaders, SeqViewSpecDoc } from "@molstar/seq-core/format";
import { type CytoBands, ideogramNavigator, type SeqViewerRegistryEntry } from "@molstar/seq-viewer";
import { hg38Streams } from "../../../../core/_test/format/stream-provider.js";
import cytobands from "../../../../viewer/_test/fixtures/hg38-cytobands.json";

/** The registry entry that adds the `ideogram` navigator (name it with `ViewerOptions.navigator`). Keep it stable. */
export const IDEOGRAM_REGISTRY: readonly SeqViewerRegistryEntry[] = [
  { navigators: [ideogramNavigator(cytobands.chromosomes as unknown as CytoBands)] },
];

export const IDEOGRAM_SECTION = "browser.main";
/** The genome as blocks, one a chromosome: the reference's letters, genes and variants (all synthetic). */
export function ideogramDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "hg38" });
  const genome = b.data.source({ ref: "genome", stream: "genome", data_kind: "letters" });
  const hg38 = b.data.complex({ ref: "hg38", source_ref: genome, label: "GRCh38" });
  for (const k of ["genes", "variants"]) {
    b.data.source({ ref: `${k}_src`, stream: k, data_kind: "items", coordinate_system: "hg38" });
    b.data.annotation({ ref: k, geometry: "ranges", source_ref: `${k}_src` });
  }
  const v = b.view({ ref: "browser", label: "hg38" });
  v.block({ complex: hg38 });
  v.track({ label: "Ruler" }).ruler();
  v.track({ label: "Sequence" }).letters({ color: b.colors.scheme("nucleotide") });
  v.track({ label: "Genes", height: 4 }).ranges({ data: "genes", labels: true });
  v.track({ label: "Variants" }).variants({ data: "variants" });
  return b.getDocument();
}
export function ideogramLoaders(): Loaders {
  const s = hg38Streams({ letters: true });
  return { load: () => Promise.reject(new Error("no files")), stream: s.stream };
}

/** Places to go: label, block (chr1 is 0, chr17 16, chrX 22) and its [start, stop). */
export const IDEOGRAM_PLACES: readonly (readonly [string, number, number, number])[] = [
  ["chr1", 0, 0, 248_956_422],
  ["chr1 centromere", 0, 121_700_000, 125_100_000],
  ["chr17 p53 region", 16, 7_000_000, 8_500_000],
  ["chr17 q arm", 16, 25_100_000, 83_257_441],
  ["chrX", 22, 0, 156_040_895],
];
