/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The package's entries name what they export, the host's contract (NG-15): no `export *`, internals reached by path.

import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import * as entry from "../src/index.js";

describe("the entries", () => {
  it("name their exports: `.` is materialize, its limit, the supported set and extentLoci", () => {
    for (const f of ["../src/index.ts", "../src/mount.ts"])
      expect(readFileSync(new URL(f, import.meta.url), "utf8")).not.toMatch(/^export \* from/m);
    expect(Object.keys(entry).sort()).toEqual([
      "NIGHTINGALE_MAX_POSITIONS",
      "NightingaleSupported",
      "extentLoci",
      "materialize",
    ]);
  });
});
