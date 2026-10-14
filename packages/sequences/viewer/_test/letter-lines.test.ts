/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Letters centered on their cells (seqspec:viewer#the-render-backend, text placement): every letter the `letters`,
// `rows`, `of` and `emphasis` draws place sits on the alphabetic baseline that centers a capital's ink on its cell, a
// lower-case letter's x-height (a gap's too), from the font's ink heights at the DPR (`TextManager.ink`). The pixels
// are _test/canvas/letter-centering.test.ts's; here the recordings, at DPR 1 and 2.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { LETTER_LINE, TextBaseline } from "../src/draw/api.ts";
import type { SeqViewerContext } from "../src/index.ts";
import { LayerFrame, Op, RecordingBackend } from "../src/testing/index.ts";
import { drawn, entryOf, layerOf, record, rects, texts } from "./draw-fixtures.ts";

const PROTEIN = "MSEEKLKELAEKAIRLAKEGNPEE";

/** A protein's letters over their cells, `letter_map` as given, with an emphasis over residues 3 to 6. */
function protein(map: string | null = null) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: "p", label: "P", alphabet: "protein", value: PROTEIN });
  const f = b.data.annotation({ ref: "f", geometry: "ranges", target: a, value: { start: [3], stop: [6] } });
  const s = b.view({ ref: "v", label: "V" }).section({ ref: "v.s", height: "box" }).block({ sequence: a });
  s.track({ label: "P" })
    .letters({ letter_map: map } as never)
    .emphasis({ data: f });
  return b.getDocument();
}
/** An alignment with lower case and gaps. */
function alignment() {
  const b = createSeqViewSpec();
  const a = b.data.alignment({
    ref: "aln",
    value: ["MSEEKLKELAEKAIRL", "msee-lkelaEK..RL", "MSDEKLKE-AEKAIRL"],
    flavor: "aligned",
  } as never);
  b.view({ ref: "v", label: "V", height: "box" }).block({ alignment: a }).track({ label: "Rows" }).rows();
  return b.getDocument();
}

/** A gene's protein through its codons (`of`), M shown lower case. */
function gene() {
  const b = createSeqViewSpec();
  const g = b.data.sequence({ ref: "gene", alphabet: "dna", value: "ATGGCTTCCAAGCACTGA" });
  const p = translate(b.data, { ref: "p", from: g, cds: [[0, 18]] });
  const v = b.view({ ref: "v", label: "V", height: "box" }).block({ sequence: g });
  v.track({ label: "Protein" }).letters({ of: p, letter_map: { M: "m" } } as never);
  return b.getDocument();
}

/** The batches' baselines in a recording. */
function baselines(rec: RecordingBackend): number[] {
  const out: number[] = [];
  for (let k = 0; k < rec.opCount; k++) if (rec.ops[4 * k] === Op.text) out.push(rec.args[rec.ops[4 * k + 3] + 2]);
  return out;
}
/** The baseline that centers `ch`'s line (a capital's or the x-height) on the band [y0, y1) of device px. */
function line(ctx: SeqViewerContext, font: number, dpr: number, ch: string, y0: number, y1: number): number {
  const ink = ctx.text.ink(font, Math.round(dpr * 64));
  return Math.ceil((y0 + y1 + ink[LETTER_LINE[ch.charCodeAt(0) & 127]]) / 2 - 0.5);
}

describe("letters on the line that centers them on their cell", () => {
  for (const dpr of [1, 2])
    for (const map of [null, "lower"])
      it(`letters and their emphasis, DPR ${dpr}, letter_map ${map}`, async () => {
        const t = await drawn(protein(map), { width: 900, height: 200, dpr });
        const { rec } = record(t.ctx, "P");
        const ts = texts(rec);
        expect(ts.map((d) => d.text).join("")).toBe(map ? PROTEIN.toLowerCase() : PROTEIN);
        expect(baselines(rec).every((b) => b === TextBaseline.alphabetic)).toBe(true);
        // Each letter on the line that centers it on the cell under it.
        const cells = rects(rec).filter((r) => r[4] >= 0 && r[3] > 4);
        for (const d of ts) {
          const c = cells.find(([x, , w]) => d.x >= x && d.x < x + w);
          expect(c, d.text).toBeDefined();
          const [, y, , h] = c as number[];
          expect(d.y, d.text).toBe(line(t.ctx, d.font, dpr, d.text, y, y + h));
        }
        // The emphasis redraws its letters where the letters are.
        const { k, section } = layerOf(t.ctx, "P");
        const over = new RecordingBackend();
        const entry = t.ctx.registry.get("renderers", "emphasis");
        if (entry) new LayerFrame(t.ctx, entry, section, k + 1).draw(over);
        const bold = texts(over);
        expect(bold.length).toBe(3);
        for (const d of bold) expect(d.y, d.text).toBe(ts.find((x) => x.x === d.x)?.y);
        t.ctx.dispose();
      });

  for (const dpr of [1, 2])
    it(`another polymer's letters (\`of\`), three letters and one, DPR ${dpr}`, async () => {
      const t = await drawn(gene(), { width: 1200, height: 200, dpr });
      const { rec } = record(t.ctx, "Protein");
      const ts = texts(rec).sort((x, y) => x.x - y.x);
      expect(ts.map((d) => d.text)).toEqual(["m", "Ala", "Ser", "Lys", "His"]);
      const cells = rects(rec).filter((r) => r[4] >= 0 && r[3] > 4);
      const [, y, , h] = cells[0];
      // "Met" shown as "m": the x-height's line; the three-letter forms a capital's.
      for (const d of ts) expect(d.y, d.text).toBe(line(t.ctx, d.font, dpr, d.text, y, y + h));
      expect(ts[0].y).toBeLessThan(ts[1].y);
      t.ctx.dispose();
    });

  for (const dpr of [1, 2])
    it(`an alignment's rows, capitals, lower case and gaps, DPR ${dpr}`, async () => {
      const t = await drawn(alignment(), { width: 600, height: 300, dpr });
      const { rec } = record(t.ctx, "Rows");
      const e = entryOf(t.ctx, "Rows");
      const ts = texts(rec);
      expect(ts.length).toBeGreaterThan(40);
      expect(baselines(rec).every((b) => b === TextBaseline.alphabetic)).toBe(true);
      // A row each entry height from the entry's top, its cells inset a px above and below.
      const lines = new Set<string>();
      for (const d of ts) {
        const row = Math.floor((d.y / dpr - e.top) / e.height);
        const top = e.top + row * e.height;
        expect(d.y, d.text).toBe(
          line(t.ctx, d.font, dpr, d.text, Math.round((top + 1) * dpr), Math.round((top + e.height - 1) * dpr)),
        );
        lines.add(`${row} ${LETTER_LINE[d.text.charCodeAt(0)]} ${d.y}`);
      }
      // A capital's line and the x-height's are a device px or more apart at these sizes.
      const rows = [...lines].map((l) => l.split(" ").map(Number));
      const row1 = rows.filter(([r]) => r === 1);
      expect(new Set(row1.map(([, l]) => l))).toEqual(new Set([0, 1]));
      expect(new Set(row1.map(([, , y]) => y)).size).toBe(2);
      t.ctx.dispose();
    });
});
