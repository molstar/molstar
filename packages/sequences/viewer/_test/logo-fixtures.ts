/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The logo's document for RN2's tests (browser-safe: the builder and the synthetic family only).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { familyRows } from "./families.ts";

/** The synthetic family's alignment (30 columns, `rows` members; or `columns` random ones) with its profile as a logo and a matrix. */
export function logoDoc(rows = 12, columns = 30) {
  const { labels, rows: family } = familyRows(rows);
  let seed = 5;
  const rnd = () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return seed / 2 ** 32;
  };
  const seqs =
    columns === 30
      ? family
      : family.map(() =>
          Array.from({ length: columns }, () => "ACDEFGHIKLMNPQRSTVWY"[Math.floor(rnd() * 20)]).join(""),
        );
  const sto = ["# STOCKHOLM 1.0", ...seqs.map((r, k) => `${labels[k].padEnd(14)}${r}`), "//", ""].join("\n");
  const b = createSeqViewSpec({ title: "Logo" });
  const aln = b.data.alignment({ ref: "family", source: { value: sto, format: "stockholm" } });
  const profile = b.data.computed({ ref: "profile", provider: "profile", alignment: aln });
  const v = b.view({ ref: "fam", label: "Family", style: "msa", height: "box", height_rows: 10 });
  v.block({ alignment: aln });
  v.track({ label: "Logo" }).logo({ data: profile });
  v.track({ label: "Matrix" }).matrix({ data: profile });
  return b.getDocument();
}
