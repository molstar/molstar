/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Small hand-written texts for the parsers and alignments (no network). The SH2 excerpt is written in the
// shape of Pfam's PF00017 seed (Stockholm, interleaved blocks, `#=GC RF` with `.` over insert columns,
// lower-case residues and `.` gaps in them, `NAME/start-end` headers); its rows are short made-up stretches
// of SH2 domains, not the seed's. PF00042 is the viewer spike's fixture (Pfam, CC0).

import { readFileSync } from "node:fs";

/** 33 columns: 30 match columns and inserts at 10, 11 and 22. */
export const SH2_STOCKHOLM = `# STOCKHOLM 1.0
#=GF ID   SH2
#=GF AC   PF00017.30
#=GF CC   An excerpt in the seed's shape,
#=GF CC   written by hand.
#=GS SRC_HUMAN/151-180     AC P12931.3
#=GS PTN11_HUMAN/6-36      AC Q06124.2

SRC_HUMAN/151-180          WYFGKITRRE..SERLLLNA
PTN11_HUMAN/6-36           WFHPNITGVEkdAENLLLTR
#=GR PTN11_HUMAN/6-36 SS   -EEEE---HHhhHHHHH---
GRB2_HUMAN/60-90           WFFGKIPRAKe.AEEMLSKQ
P85A_HUMAN/333-360         WYWGDISREE..VNEKLRDT
#=GC RF                    xxxxxxxxxx..xxxxxxxx

SRC_HUMAN/151-180          EN.PRGTFLVRES
PTN11_HUMAN/6-36           GV.DGSFLARPS-
#=GR PTN11_HUMAN/6-36 SS   ---EEEEE-----
GRB2_HUMAN/60-90           RHgDGAFLIRES-
P85A_HUMAN/333-360         A-.DGTFLVRDA-
#=GC RF                    xx.xxxxxxxxxx
//
`;

/** The same rows without `RF`: insert columns come from case. */
export const SH2_NO_RF = SH2_STOCKHOLM.split("\n")
  .filter((l) => !l.startsWith("#=GC RF"))
  .join("\n");

/** The same rows as A3M: lower case an insert, `.` dropped. */
export const SH2_A3M = `#A3M written by hand
>SRC_HUMAN/151-180 Proto-oncogene tyrosine-protein kinase Src
WYFGKITRRESERLLLNAEN
PRGTFLVRES
>PTN11_HUMAN/6-36
WFHPNITGVEkdAENLLLTRGVDGSFLARPS-
>GRB2_HUMAN/60-90
WFFGKIPRAKeAEEMLSKQRHgDGAFLIRES-
>P85A_HUMAN/333-360
WYWGDISREEVNEKLRDTA-DGTFLVRDA-
`;

/** A tRNA-like RNA alignment with WUSS structure and a pseudoknot (Aa). */
export const RNA_STOCKHOLM = `# STOCKHOLM 1.0
seq1  GCGGAUUUAGCUC
seq2  GCGGA-UUAGCUC
#=GC SS_cons <<<A.[..]a>>>
//
# STOCKHOLM 1.0
seq3  ACGU
//
`;

/** Rows: [label, first residue number, gapped]. */
export function pf00042(): { label: string; start: number; gapped: string }[] {
  const url = new URL("../fixtures/pf00042.json", import.meta.url);
  const json = JSON.parse(readFileSync(url, "utf8")) as { rows: [string, number, string][] };
  return json.rows.map(([label, start, gapped]) => ({ label, start, gapped }));
}
