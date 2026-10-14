/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The format's coded columns (the format's resolve/codes.ts) are compared with their named constants, never with a
// number, in every package that reads or writes them. A sweep over the sources: a coded column or `heightMode` against
// a numeric literal, either way round, fails.

import { readdirSync, readFileSync, statSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";

const root = join(import.meta.dirname, "../..");
const CODED =
  /(?:\bheightMode|sections\.(?:axis|arrangement|ligandForm|height|rowSizing)\[[^\]]*\]|layers\.sizing\[[^\]]*\]|\blayout|section\.layout|row\.axis)\s*[!=]==?\s*-?\d|[!=]==?\s*-?\d+\s*\?\s*"(?:box|columns|text|wrapped)"/;

function sources(dir: string, out: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    const f = join(dir, name);
    if (statSync(f).isDirectory()) {
      if (name !== "node_modules") sources(f, out);
    } else if (/\.tsx?$/.test(name) && !name.endsWith(".d.ts")) out.push(f);
  }
  return out;
}

describe("the format's codes are named", () => {
  for (const dir of ["viewer/src", "core/src/format", "nightingale/src"])
    it(`${dir} compares no coded column with a number`, () => {
      const found = sources(join(root, dir)).flatMap((f) =>
        readFileSync(f, "utf8")
          .split("\n")
          .map((line, i) => ({ at: `${f.slice(root.length + 1)}:${i + 1}`, line }))
          .filter(({ line }) => CODED.test(line))
          .map(({ at, line }) => `${at}: ${line.trim()}`),
      );
      expect(found).toEqual([]);
    });
});
