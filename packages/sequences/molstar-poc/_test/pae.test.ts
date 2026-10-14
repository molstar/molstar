/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// AlphaFold DB's predicted aligned error (seqspec:sources#formats-and-services, seqspec:libraries#the-mol-adapter): the
// `afdb` service's `pae` resource, the `afdb_pae` format in both of AlphaFold DB's forms, and `molstar:pae`'s grid over
// a complex's units in order (scored across, aligned on down), over the shared synthetic fixture of the fixture model.

import { readFileSync } from "node:fs";
import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { compile, createRegistries, KernelRegistry, register } from "@molstar/seq-core/format";
import { KeyKind } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { MolstarServices, PaeFormat, PaeMatrix, PaeProvider } from "../src/index.js";

const FIXTURE = readFileSync(new URL("../../core/_test/fixtures/af-barnase-barstar-pae.json", import.meta.url), "utf8");
const ROWS = JSON.parse(FIXTURE)[0].predicted_aligned_error as number[][];
const N = ROWS.length;
const bytes = (s: string) => new TextEncoder().encode(s);
const issues: { code: string; message: string }[] = [];
const kernel = createRegistries();
register(kernel, KernelRegistry);
const ctx = { registries: kernel, report: (i: { code: string; message: string }) => issues.push(i) } as never;

describe("AlphaFold DB's PAE", () => {
  it("is the afdb service's pae resource, by entry id or version", () => {
    const afdb = MolstarServices.find((s) => s.name === "afdb");
    expect(afdb?.resolve("AF-P04637-F1", "pae")).toEqual({
      uri: "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-predicted_aligned_error_v6.json",
      format: "afdb_pae",
    });
    expect(afdb?.resolve("AF-P04637-F1-model_v4", "pae")?.uri).toBe(
      "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-predicted_aligned_error_v4.json",
    );
  });

  it("reads both forms into one matrix, by scored token (the JSON transposed); anything else is reported", () => {
    const m = PaeFormat.read(bytes(FIXTURE), {} as never, ctx) as { value: PaeMatrix };
    expect(PaeMatrix.is(m.value)).toBe(true);
    // Scored token j's row: its error aligned on each token i, the JSON's [i][j].
    expect([m.value.size, m.value.values[1], m.value.values[N + 5]]).toEqual([
      N,
      Math.fround(ROWS[1][0]),
      Math.fround(ROWS[5][1]),
    ]);
    // The current layout as a plain object, not a list of one.
    const plain = {
      predicted_aligned_error: [
        [0, 1.5],
        [2.5, 0],
      ],
      max_predicted_aligned_error: 31.75,
    };
    const p = PaeFormat.read(bytes(JSON.stringify(plain)), {} as never, ctx) as { value: PaeMatrix };
    expect(Array.from(p.value.values)).toEqual([0, 2.5, 1.5, 0]);
    const old = {
      residue1: [1, 1, 2, 2],
      residue2: [1, 2, 1, 2],
      distance: [0, 3, 4, 0],
      max_predicted_aligned_error: 31,
    };
    const o = PaeFormat.read(bytes(JSON.stringify([old])), {} as never, ctx) as { value: PaeMatrix };
    expect(Array.from(o.value.values)).toEqual([0, 4, 3, 0]);
    issues.length = 0;
    expect((PaeFormat.read(bytes("{}"), {} as never, ctx) as { value: unknown }).value).toBe(null);
    // The older form's residues are 1..N: a 0 or one past N says which pair.
    const bad = { residue1: [1, 1, 2, 3], residue2: [1, 2, 1, 2], distance: [0, 3, 4, 0] };
    expect((PaeFormat.read(bytes(JSON.stringify([bad])), {} as never, ctx) as { value: unknown }).value).toBe(null);
    expect(issues.map((i) => [i.code, i.message])).toEqual([
      ["svs.source.format", expect.stringContaining("no predicted_aligned_error")],
      ["svs.source.format", "AlphaFold DB's PAE: pair 3 (residue1 3, residue2 2) is not two residues 1..2"],
    ]);
  });

  it("molstar:pae: a grid row per pair of the complex's units, scored across, aligned on down", async () => {
    const reg = createRegistries();
    register(reg, [KernelRegistry, { formats: [PaeFormat], providers: [PaeProvider] }]);
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "A".repeat(157) });
    const s = b.data.sequence({ ref: "s", value: "S".repeat(N - 157) });
    const m = b.data.complex({ ref: "model" });
    m.unit({ ref: "model-a", sequence: a });
    m.unit({ ref: "model-s", sequence: s });
    // A ligand in the complex is no token of the PAE: the polymer units alone count.
    m.unit({ ref: "model-lig", ligand: b.data.ligand({ ref: "lig", value: ["C", "N", "O"] }) });
    const src = b.data.source({ ref: "pae-src", value: FIXTURE, format: "afdb_pae" });
    const pae = b.data.computed({ ref: "pae", provider: "molstar:pae", options: { of: m, source: src } });
    b.view({ ref: "v", label: "V" }).track({ label: "PAE" }).dot_plot({ data: pae });
    const { set, issues: found } = await compile(b.getDocument(), reg, {
      load: () => Promise.reject(new Error("none")),
    });
    expect(found.filter((i) => i.severity === "error")).toEqual([]);
    const t = { kind: -1, id: -1, row: -1 };
    expect(set.byName("pae", t) && t.kind === KeyKind.annotation).toBe(true);
    const g = set.annotations[t.id];
    const w = N - 157;
    expect([g?.geometry, g?.count, g?.grid?.rows.map((r) => r.length)]).toEqual([
      "grid",
      4,
      [157 * 157, 157 * w, w * 157, w * w],
    ]);
    // Float32 end to end: the provider's rows are the grid's.
    expect(g?.grid?.rows.every((r) => r instanceof Float32Array)).toBe(true);
    // A × S at (scored a 3, aligned on s 7): the matrix's row 157 + 7, column 3.
    expect(g?.grid?.rows[1][3 * w + 7]).toBe(Math.fround(ROWS[157 + 7][3]));
  });

  it("molstar:pae on a single chain: its one row is the matrix itself, no copy", async () => {
    // The matrix as read, to compare the grid's row with.
    let read: PaeMatrix | null = null;
    const format = {
      ...PaeFormat,
      read: (...args: Parameters<typeof PaeFormat.read>) => {
        const r = PaeFormat.read(...args) as { value: PaeMatrix };
        read = r.value;
        return r;
      },
    } as typeof PaeFormat;
    const reg = createRegistries();
    register(reg, [KernelRegistry, { formats: [format], providers: [PaeProvider] }]);
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "A".repeat(N) });
    const m = b.data.complex({ ref: "model" });
    m.unit({ ref: "model-a", sequence: a });
    const src = b.data.source({ ref: "pae-src", value: FIXTURE, format: "afdb_pae" });
    const pae = b.data.computed({ ref: "pae", provider: "molstar:pae", options: { of: m, source: src } });
    b.view({ ref: "v", label: "V" }).track({ label: "PAE" }).dot_plot({ data: pae });
    const { set } = await compile(b.getDocument(), reg, { load: () => Promise.reject(new Error("none")) });
    const t = { kind: -1, id: -1, row: -1 };
    expect(set.byName("pae", t)).toBe(true);
    const row = set.annotations[t.id]?.grid?.rows[0];
    // (scored 3, aligned on 7): the JSON's row 7, column 3.
    expect(row?.[3 * N + 7]).toBe(Math.fround(ROWS[7][3]));
    expect(row).toBe((read as PaeMatrix | null)?.values);
  });
});
