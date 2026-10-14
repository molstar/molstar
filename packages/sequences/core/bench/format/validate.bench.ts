/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Validation cost at genome scale (VW2; V5): the 3M-row coded genome table, and 1M rows with 10 snapshots,
// each checked cold (a new params object, so the node's own checks run) and again (the memo). The old stack, measured
// once on the same documents: 1,783 ms for 3M rows, 532 ms for 1M, 6,933 ms for 1M with 10 snapshots (VW2).

import { bench, describe } from "vitest";
import { kernel } from "../../_test/format/fixtures.ts";
import { codedDoc, freshParams } from "../../_test/format/table-fixtures.ts";
import { validate } from "../../src/format/index.ts";

const reg = kernel();
const genome = codedDoc(3_000_000);
const snapshots = codedDoc(1_000_000, { snapshots: 10 });

describe("validate", () => {
  bench("3M coded rows, cold", () => {
    validate(freshParams(genome), reg);
  });
  bench("3M coded rows, again", () => {
    validate(genome, reg);
  });
  bench("1M rows with 10 snapshots, cold", () => {
    validate(freshParams(snapshots), reg);
  });
});
