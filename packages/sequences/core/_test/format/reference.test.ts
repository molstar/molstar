/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// V6: the generated reference and the CLIs.

import { spawnSync } from "node:child_process";
import { mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { describe, expect, it } from "vitest";
import { printReference, toJSON } from "../../src/format/index.js";
import { CASES } from "./cases.js";
import { fullDoc, kernel } from "./fixtures.js";

describe("the generated reference", () => {
  it("matches lib/seqstar/_spec/reference/format.md (a change to a table, a type or a code shows here)", async () => {
    await expect(printReference(kernel())).toMatchFileSnapshot("../../../_spec/reference/format.md");
  });

  it("defines the types the tables name", () => {
    const text = printReference(kernel());
    const types = text.slice(text.indexOf("## Types"), text.indexOf("## Layer types"));
    for (const heading of ["Table<geometry>", "Palette", "Glyph", "Color", "Label", "ThemeRef"]) {
      expect(types).toContain(`### \`${heading}\``);
    }
    // The columns of every geometry, and the three palette kinds with their params.
    for (const geometry of ["ranges", "values", "pairs", "matrix", "grid", "rows"])
      expect(types).toContain(`| \`${geometry}\` |`);
    for (const kind of ["categorical", "discrete", "continuous"]) expect(types).toContain(`\`kind: "${kind}"\``);
    expect(types).toContain("`value_domain`");
    expect(types).toContain("`b_target`");
  });

  it("lists the codes the spec defines", () => {
    const text = printReference(kernel());
    for (const code of ["svs.annotation.rows", "svs.unit.orphaned", "svs.provider.unmapped"]) {
      expect(text).toContain(`| \`${code}\` | warning |`);
    }
    for (const code of ["svs.alignment.inserts", "svs.ref.pending", "svs.query.note"])
      expect(text).toContain(`| \`${code}\` | info |`);
  });
});

describe("svs-validate and svs-print-schema", { tags: ["slow"] }, () => {
  const dir = mkdtempSync(join(tmpdir(), "svs-cli-"));
  const path = (name: string) => join(dir, name);
  const cli = (name: string, ...args: string[]) =>
    spawnSync(process.execPath, [fileURLToPath(new URL(`../../src/format/cli/${name}.ts`, import.meta.url)), ...args], {
      encoding: "utf8",
    });
  writeFileSync(path("good.svsj"), toJSON(CASES.case1()));
  const broken = fullDoc();
  broken.root.children[1].params.label = 5;
  writeFileSync(path("bad.svsj"), JSON.stringify(broken));
  const other = fullDoc();
  other.root.children[1].children[0].children[0].params.sequence = "nope";
  writeFileSync(path("bad2.svsj"), JSON.stringify(other));

  it("exits 0 for a valid document, 1 for errors, 2 for usage or a file that cannot be read", () => {
    const good = cli("svs-validate", path("good.svsj"));
    expect([good.status, good.stdout]).toEqual([0, ""]);
    const bad = cli("svs-validate", path("bad.svsj"));
    expect(bad.status).toBe(1);
    expect(bad.stdout).toContain("svs.params.type");
    expect(cli("svs-validate").status).toBe(2);
    expect(cli("svs-validate", "--nope", path("good.svsj")).status).toBe(2);
    expect(cli("svs-validate", path("missing.svsj")).status).toBe(2);
  });

  it("validates every file: the exit code and the issues are those of all of them, not the first", () => {
    const run = cli("svs-validate", path("good.svsj"), path("bad.svsj"), path("bad2.svsj"));
    expect(run.status).toBe(1);
    expect(run.stdout).toContain(`${path("bad.svsj")}#`);
    expect(run.stdout).toContain(`${path("bad2.svsj")}#`);
    expect(run.stdout).not.toContain(`${path("good.svsj")}#`);
    expect(run.stderr).toContain(`${path("good.svsj")}: 0 error(s)`);
    const last = cli("svs-validate", path("good.svsj"), path("bad2.svsj"));
    expect(last.status).toBe(1);
    const unreadable = cli("svs-validate", path("missing.svsj"), path("bad.svsj"));
    expect([unreadable.status, unreadable.stdout.includes("svs.params.type")]).toEqual([2, true]);
  });

  it("svs-print-schema writes the reference", () => {
    const run = cli("svs-print-schema");
    expect(run.status).toBe(0);
    expect(run.stdout).toBe(printReference(kernel()));
  });
});
