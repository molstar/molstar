/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The document of `Seq*/Cases/Genome browser › hg38, live` (lib/seqstar/storybook, ST4b): four streamed sources (`genome` letters,
// `gc` values, `genes` and `clinvar` items) over one complex, shared by the story (which serves them live through
// @molstar/seq-core/stream/gmod) and jbrowse-story.test.ts (which opens it offline, over the fixture providers), so the story's
// document cannot silently fail to open.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";

export const JBROWSE_SECTION = "jbrowse.main";

/**
 * The demo on one chromosome (its unit, so the navigator spans it; the complex's other members stay in the set),
 * opening on [start, stop) of it: by default TP53, about 26 kb (ST4b: never a whole genome at first).
 */
export function jbrowseDoc(
  chr = "chr17",
  start = 7_661_779,
  stop = 7_687_538,
  o: { genome?: boolean } = {},
): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "hg38, live" });
  const genome = b.data.source({ ref: "genome", stream: "genome", data_kind: "letters" });
  const hg38 = b.data.complex({ ref: "hg38", source_ref: genome, label: "GRCh38" });
  b.data.source({ ref: "gc_src", stream: "gc", data_kind: "values", coordinate_system: "hg38" });
  const gc = b.data.annotation({ ref: "gc", geometry: "values", source_ref: "gc_src", label: "GC %" });
  for (const k of ["genes", "clinvar"]) {
    b.data.source({ ref: `${k}_src`, stream: k, data_kind: "items", coordinate_system: "hg38" });
    b.data.annotation({ ref: k, geometry: "ranges", source_ref: `${k}_src` });
  }
  const v = b.view({ ref: "jbrowse", label: "hg38" });
  // `genome`: every chromosome a block (the demo before David's review; its navigation test's).
  if (o.genome) v.block({ complex: hg38 });
  else v.block({ unit: `hg38/${chr}` });
  v.focus({ target: `hg38/${chr}`, start, stop });
  v.track({ label: "" }).ruler();
  v.track({ label: "Reference" }).letters({ color: b.colors.scheme("nucleotide") });
  v.track({ label: "GENCODE v50", height: 12 }).ranges({
    data: "genes",
    labels: true,
    // Where the lanes do not fit (TP53 has about 30 transcripts), one lane a gene: its transcripts folded, as JBrowse 2's
    // "collapse"; the track menu's "Expand rows" shows them all.
    collapse: "label",
    color: { theme: "field_palette", params: { field: "category" } },
  });
  v.track({ label: "ClinVar" }).variants({
    data: "clinvar",
    color: { theme: "field_palette", params: { field: "category" } },
    glyph: {
      theme: "field_map",
      params: {
        field: "category",
        map: { SNV: "circle", deletion: "triangle", insertion: "diamond" },
        default: "square",
      },
    },
  } as never);
  // GC % by a sequential ramp over the file's range (a values track draws in the ink without a theme).
  v.track({ label: "GC content", height: 2 }).values({
    data: gc,
    color: { theme: "field_palette", params: { field: "value", palette: { kind: "continuous", colors: "Viridis" } } },
  } as never);
  return b.getDocument();
}
