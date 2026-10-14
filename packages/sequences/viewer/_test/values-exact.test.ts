/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Values are JSON numbers (seqspec:format#values), the prepare half: a decimal such as 0.15 written in a document comes
// back exactly through load, compile, resolve and prepare, dense and per row, in memory, as Float32 and after JSON; the
// pyramid's base bins hold its exact minimum and maximum. The colors and the legend are the theme providers' (LY3).

import { readFileSync } from "node:fs";
import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { fromJSON, type Palette, type SeqViewSpecDoc, toJSON } from "@molstar/seq-core/format";
import { Annotations } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import type { PerPosition } from "../src/prepare.js";
import { layer, opened, preparedRenderers } from "./layout-fixtures.js";

const recorded = JSON.parse(
  readFileSync(new URL("../../core/_test/format/palettes/mvs-recorded.json", import.meta.url), "utf8"),
) as { cases: { name: string; palette: Palette; values: number[] }[] };
const DECIMALS = recorded.cases.filter((c) => c.name.startsWith("decimals:"));

type Form = "JSON" | "Float32" | "Float32, then JSON";
/** Dense (no position column) and per-row values, a `color` track over each and a `values` track. */
function doc(values: number[], form: Form): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "Decimals" });
  const s = b.data.sequence({ ref: "s", value: "MKVLAGHEKL".slice(0, values.length) });
  const value = form === "JSON" ? values : Float32Array.from(values);
  const dense = b.data.annotation({ ref: "dense", geometry: "values", target: s, value: { value } });
  const index = values.map((_, k) => k);
  const rows = b.data.annotation({ ref: "rows", geometry: "values", target: s, value: { index, value } });
  const v = b.view({ ref: "v", label: "V" });
  v.block({ sequence: s });
  v.track({ label: "Dense" }).color({ data: dense });
  v.track({ label: "Rows" }).color({ data: rows });
  v.track({ label: "Values" }).values({ data: dense });
  const built = b.getDocument();
  if (form === "Float32") return built;
  expect(JSON.parse(toJSON(built)).root.children[0].children[1].params.value.value).toEqual(values);
  const { doc: d, issues } = fromJSON(toJSON(built));
  expect(issues).toEqual([]);
  return d as SeqViewSpecDoc;
}

const FORMS: Form[] = ["JSON", "Float32", "Float32, then JSON"];

describe("values are JSON numbers, exactly, through prepare", () => {
  for (const [form, c] of FORMS.flatMap((f) => DECIMALS.map((c) => [f, c] as const))) {
    it(`${c.name} (${form})`, async () => {
      const t = await opened(doc(c.values, form), { renderers: preparedRenderers() });
      const dense = layer(t.ctx, "Dense");
      const a = dense.m.set.annotations[dense.m.resolved.data[dense.m.resolved.layers.dataStart[dense.k]]];
      if (!a?.dense) throw new Error("no dense values");
      expect(c.values.map((_, i) => Annotations.valueAt(a, 0, i))).toEqual(c.values);
      expect(Array.from((dense.state as PerPosition).value[0] ?? [])).toEqual(c.values);
      // Per row: each position its row's value, and its row.
      const rows = layer(t.ctx, "Rows").state as PerPosition;
      expect([Array.from(rows.value[0] ?? []), Array.from(rows.row[0] ?? [])]).toEqual([
        c.values,
        c.values.map((_, k) => k),
      ]);
      const values = layer(t.ctx, "Values").state as PerPosition;
      expect(Array.from(values.value[0] ?? [])).toEqual(c.values);
      expect([values.pyramid[0]?.min[0][0], values.pyramid[0]?.max[0][0]]).toEqual([
        Math.min(...c.values.slice(0, 4)),
        Math.max(...c.values.slice(0, 4)),
      ]);
    });
  }
});
