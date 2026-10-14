/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The exports guard (lib/seqstar/scripts/exports.ts): no name is exported by two entries of the packages as different
// declarations, in the type or the value namespace; and the guard itself, on a folder of packages written here.

import { mkdirSync, mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { collisions, type Entry, entries, exported, unused } from "../../scripts/exports.js";

const pkg = (dir: string, name: string, exports: Record<string, string>, files: Record<string, string>) => {
  mkdirSync(join(dir, name, "src"), { recursive: true });
  writeFileSync(join(dir, name, "package.json"), JSON.stringify({ name, exports }));
  for (const [f, text] of Object.entries(files)) writeFileSync(join(dir, name, "src", f), text);
};

describe("export guard", () => {
  // A TypeScript program over every entry: about 1.5 s alone, past the default 5 s under a full parallel run. The slow
  // tier's: `mise run seqstar:exports` runs the same guard before every commit.
  it("finds no collision across the entries of the packages", { tags: ["slow"], timeout: 30_000 }, () => {
    expect(entries().length).toBeGreaterThan(10);
    expect(collisions(exported(entries()))).toEqual([]);
  });

  it("names a clash in either namespace, and lets a re-export and the two namespaces be", () => {
    const dir = mkdtempSync(join(tmpdir(), "seqstar-exports-"));
    pkg(
      dir,
      "a",
      { ".": "./src/index.ts" },
      {
        "index.ts":
          "export interface Field { a: 1 }\nexport const Same = 1;\nexport interface Geometry { g: 1 }\nexport const Clash = 1;",
      },
    );
    pkg(
      dir,
      "b",
      { ".": "./src/index.ts", "./x/*": "./src/x/*.ts" },
      {
        "index.ts":
          "export type Field = { b: 2 };\nexport { Same } from '../../a/src/index.js';\nexport const Geometry = 1;\nexport const Clash = 2;",
      },
    );
    mkdirSync(join(dir, "b", "src", "x"));
    writeFileSync(join(dir, "b", "src", "x", "ruler.ts"), "export const Ruler = 1;");
    const list = entries(dir);
    expect(list.map((e: Entry) => e.subpath)).toEqual([".", ".", "./x/ruler"]);
    expect(collisions(exported(list))).toEqual([
      "type Field: @molstar/seq-a, @molstar/seq-b",
      "value Clash: @molstar/seq-a, @molstar/seq-b",
    ]);
  });
});

describe("the two UIs", () => {
  it("may share a name (one slot contract), and are checked against everything else", () => {
    const dir = mkdtempSync(join(tmpdir(), "seqstar-exports-"));
    const ui = { ".": "./src/index.ts", "./react": "./src/react.ts", "./dom": "./src/dom.ts" };
    pkg(dir, "a", ui, {
      "index.ts": "export const Root = 1;\nexport const Other = 1;",
      "react.ts": "export const ViewTabs = 1;\nexport interface PartProps { r: 1 }\nexport const Other = 2;",
      "dom.ts": "export const ViewTabs = 2;\nexport interface PartProps { d: 1 }",
    });
    // Another package's UI entries are not alternatives of this one's.
    pkg(dir, "b", { "./dom": "./src/dom.ts" }, { "dom.ts": "export const ViewTabs = 3;" });
    expect(collisions(exported(entries(dir)))).toEqual([
      "value Other: @molstar/seq-a, @molstar/seq-a/react",
      "value ViewTabs: @molstar/seq-a/react, @molstar/seq-a/dom, @molstar/seq-b/dom",
    ]);
  });
});

describe("unused-exports listing", () => {
  // A TypeScript program over a written folder of packages, for a listing run by hand: the slow tier's.
  it("lists a name no consumer outside the package mentions, and not one that a script does", {
    tags: ["slow"],
  }, () => {
    const entry = { pkg: "viewer", subpath: ".", file: "" };
    const list = [
      { name: "zzNobodyMentionsThis", space: "value" as const, decl: "x", entry },
      { name: "collisions", space: "value" as const, decl: "y", entry },
    ];
    expect([...unused(list)]).toEqual([["@molstar/seq-viewer", ["zzNobodyMentionsThis"]]]);
  });
});
