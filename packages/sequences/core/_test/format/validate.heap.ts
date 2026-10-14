/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for validation (VW2): a table cell is checked with nothing allocated, in every column form, on RB0b's
// stepped `heapScaling` (a step is one cell); and validating the 3M-row coded genome table grows the heap by
// the same small amount as validating a quarter of it (the old check kept about 52 MB of objects, a path and a ref per
// cell).

import { describe, expect, it } from "vitest";
import { validate } from "../../src/format/index.js";
import { CellKind, CellScan, type Col, ColumnForm, scanCells } from "../../src/format/tables/cells.js";
import {
  controlCaught,
  describeScaling,
  heapGrowth,
  heapScaling,
  isFlat,
  leak,
  leaking,
  SCALED_LIMIT,
} from "../model/heap.js";
import { kernel } from "./fixtures.js";
import { CHROMOSOMES, codedDoc, freshParams } from "./table-fixtures.js";

/** Steps in the shorter run; the longer runs four times as many. */
const N = 250_000;
const LEN = 4 * N;

const col = (form: number, data: ArrayLike<unknown>, dict: readonly unknown[] | null = null, name = "c"): Col => ({
  name,
  form,
  data,
  dict,
  length: data.length,
});
const ints = Int32Array.from({ length: LEN }, (_, i) => i);
const doubles = Array.from({ length: LEN }, (_, i) => (i % 100) / 100 + 0.001);
const columns: [string, Col, number][] = [
  ["typed indices", col(ColumnForm.typed, ints), CellKind.index],
  ["typed values", col(ColumnForm.typed, Float64Array.from(doubles)), CellKind.value],
  [
    "coded refs",
    col(
      ColumnForm.dict,
      Uint8Array.from(ints, (i) => i % 24),
      CHROMOSOMES,
    ),
    CellKind.ref,
  ],
  [
    "string refs",
    col(
      ColumnForm.array,
      Array.from(ints, (i) => CHROMOSOMES[i % 24]),
    ),
    CellKind.ref,
  ],
  [
    "label strings",
    col(
      ColumnForm.array,
      Array.from(ints, (i) => `${i % 900}A`),
    ),
    CellKind.label,
  ],
  ["an array of doubles", col(ColumnForm.array, doubles), CellKind.value],
  [
    "rows",
    col(
      ColumnForm.rows,
      Array.from(ints, (i) => ({ start: i, label: "x" })),
      null,
      "start",
    ),
    CellKind.index,
  ],
];

describe("table cells allocate nothing", () => {
  it.each(columns)("%s", (_, c, kind) => {
    const s = new CellScan();
    s.reset(c);
    const r = heapScaling(N, (i) => scanCells(c, kind, i % LEN, (i % LEN) + 1, s));
    expect(s.bad, "every cell is valid").toBe(0);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("the control: a cell's check with a leak is caught (the template, CF-05)", () => {
    const [, c, kind] = columns[0];
    const s = new CellScan();
    s.reset(c);
    expect(
      controlCaught(() =>
        heapScaling(
          N,
          leaking((i) => scanCells(c, kind, i % LEN, (i % LEN) + 1, s)),
        ),
      ),
    ).toBe(true);
  });
});

describe("validating the coded genome table", () => {
  it("3M rows grow the heap as 750k do: no allocation per row", () => {
    const reg = kernel();
    /** `perRow` after each validation, once a row: the control's leak, else nothing. */
    const grow = (n: number, perRow: (k: number) => void = () => {}) => {
      const d = codedDoc(n);
      validate(d, reg);
      // The least of three: the first runs still compile (a few hundred KB, whatever the size).
      let least = Number.POSITIVE_INFINITY;
      for (let k = 0; k < 3; k++)
        least = Math.min(
          least,
          heapGrowth(() => {
            validate(freshParams(d), reg);
            for (let r = 0; r < n; r++) perRow(r);
          }),
        );
      return least;
    };
    const small = grow(750_000);
    const large = grow(3_000_000);
    expect(large - small, `${(small / 1024).toFixed(1)} KB, then ${(large / 1024).toFixed(1)} KB`).toBeLessThan(
      SCALED_LIMIT,
    );
    // The template's control (CF-05): an object a row, the old check's fault, is caught by the same comparison.
    expect(
      controlCaught(() => {
        const [short, long] = [grow(750_000, leak), grow(3_000_000, leak)];
        return { short, long, scaled: long - short };
      }),
    ).toBe(true);
  });
});
