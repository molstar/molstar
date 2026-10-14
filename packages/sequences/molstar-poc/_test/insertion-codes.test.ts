/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
// Insertion codes end to end on the Mol* side (seqspec:libraries#the-mol-adapter): a hand-written mmCIF of nine
// residues numbered 50, 51, 52, 52A, 52B, 53, 54, 55, 55A (the last unobserved, its code from pdbx_poly_seq_scheme),
// SIFTS onto a made-up accession's 70-78. The unit's `auth` labels carry the codes; a label range selects them; a
// selection reaches Mol* as the residue with the code (by label_seq_id, never by number) and comes back the same;
// UniProt numbers run across the inserted residues; atom_site's rows by auth_seq_id and pdbx_PDB_ins_code land on
// their own residues.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { compile } from "@molstar/seq-core/format";
import type { SeqSet } from "@molstar/seq-core/model";
import { compileQuery, createQueryEnv, SeqSelection } from "@molstar/seq-core/model/query";
import { bindQuery, readQuery, renderIssue } from "@molstar/seq-core/query";
import { Translator } from "@molstar/seq-core/translate";
import { StructureProperties as P, StructureElement } from "molstar/lib/mol-model/structure";
import { beforeAll, describe, expect, it } from "vitest";
import {
  complexInputs,
  MolstarSvsEntry,
  readStructure,
  StructureLink,
  selectionLoci,
  structureInput,
  structureOf,
} from "../src/index.ts";
import { registries } from "./fixtures.ts";

const AUTH: [string, number, string][] = [
  ["TRP", 50, ""],
  ["VAL", 51, ""],
  ["ALA", 52, ""],
  ["ARG", 52, "A"],
  ["ILE", 52, "B"],
  ["TYR", 53, ""],
  ["PRO", 54, ""],
  ["THR", 55, ""],
  ["ASN", 55, "A"],
];
const OBSERVED = 8;
const q = (s: string) => (s === "" ? "?" : s);
const atoms = AUTH.slice(0, OBSERVED)
  .map(
    ([comp, num, ins], i) =>
      `ATOM ${i + 1} C CA . ${comp} A 1 ${i + 1} ${q(ins)} ${(3.8 * i).toFixed(3)} 0.000 0.000 1.00 ${(10 * (i + 1)).toFixed(2)} ${num} ${comp} H CA 1`,
  )
  .join("\n");
const scheme = AUTH.map(
  ([comp, num, ins], i) => `A 1 ${i + 1} ${comp} ${num} ${num} ${q(ins)} H ${comp} ${i < OBSERVED ? num : "?"}`,
).join("\n");
/** The mmCIF text: one entity, one chain H (label_asym_id A), its SIFTS segment. */
const INSERTIONS_CIF = `data_INS1
#
_entry.id INS1
#
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'Heavy chain fragment'
#
_entity_poly.entity_id 1
_entity_poly.type 'polypeptide(L)'
_entity_poly.pdbx_seq_one_letter_code_can WVARIYPTN
#
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
${AUTH.map(([comp], i) => `1 ${i + 1} ${comp} n`).join("\n")}
#
_struct_asym.id A
_struct_asym.entity_id 1
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.ndb_seq_num
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_mon_id
_pdbx_poly_seq_scheme.auth_seq_num
${scheme}
#
loop_
_pdbx_sifts_unp_segments.entity_id
_pdbx_sifts_unp_segments.asym_id
_pdbx_sifts_unp_segments.unp_acc
_pdbx_sifts_unp_segments.segment_id
_pdbx_sifts_unp_segments.instance_id
_pdbx_sifts_unp_segments.unp_start
_pdbx_sifts_unp_segments.unp_end
_pdbx_sifts_unp_segments.seq_id_start
_pdbx_sifts_unp_segments.seq_id_end
_pdbx_sifts_unp_segments.best_mapping
_pdbx_sifts_unp_segments.identity
1 A Q0INS1 1 1 70 78 1 9 y 1.0
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
${atoms}
#
`;
const URI = "https://example.org/ins1.cif";
const bytes = () => new TextEncoder().encode(INSERTIONS_CIF);

describe("insertion codes on the Mol* side", () => {
  it("the unit's auth labels carry the codes, observed or not; label and SIFTS are by label_seq_id", async () => {
    const src = await readStructure(bytes(), "mmcif");
    const structure = await structureOf(src, 0, null);
    const issues: string[] = [];
    const r = complexInputs([{ structure, model: null }], src.frame, { unit_name: "auth_asym_id" }, (code) =>
      issues.push(code),
    );
    expect(issues).toEqual([]);
    expect(r.molecules[0].numberings).toEqual([
      { name: "label", method: "offset", first: 1 },
      { name: "uniprot:Q0INS1", method: "segments", segments: [[0, 70, 9]] },
    ]);
    const auth = r.units[0].numberings?.[0];
    expect(auth?.method === "labels" && auth.labels).toEqual([50, 51, 52, "52A", "52B", 53, 54, 55, "55A"]);
  });

  it("a code that is not letters is read as none, once reported, naming the residue and the code", async () => {
    // 52B's code written as "1" (atom_site and pdbx_poly_seq_scheme): it would read as 521.
    const odd = INSERTIONS_CIF.replace(" ILE A 1 5 B ", " ILE A 1 5 1 ").replace(
      "A 1 5 ILE 52 52 B ",
      "A 1 5 ILE 52 52 1 ",
    );
    expect(odd).not.toBe(INSERTIONS_CIF);
    const src = await readStructure(new TextEncoder().encode(odd), "mmcif");
    const structure = await structureOf(src, 0, null);
    const issues: [string, string, string][] = [];
    const r = complexInputs(
      [{ structure, model: null }],
      src.frame,
      { unit_name: "auth_asym_id" },
      (code, _, message) =>
        issues.push([code, MolstarSvsEntry.issueCodes?.find((c) => c.code === code)?.severity ?? "", message]),
    );
    const auth = r.units[0].numberings?.[0];
    expect(auth?.method === "labels" && auth.labels).toEqual([50, 51, 52, "52A", 52, 53, 54, 55, "55A"]);
    expect(issues).toEqual([
      [
        "molstar.ins_code.invalid",
        "warning",
        'chain H residue 52 has insertion code "1", not letters: read with no code, so its auth label is the number alone and may repeat its neighbor\'s; select it by label_seq_id (numbering "label") to tell them apart',
      ],
    ]);
  });

  let set: SeqSet;
  let link: StructureLink;
  const env = createQueryEnv([]);
  beforeAll(async () => {
    const b = createSeqViewSpec();
    b.data.complex({
      ref: "ab",
      provider: "molstar:structure",
      source: { uri: URI, format: "mmcif" },
      options: { unit_name: "auth_asym_id" },
    } as never);
    const uni = b.data
      .sequence({ ref: "uni", value: `${"M".repeat(69)}WVARIYPTN${"K".repeat(10)}` })
      .numbering({ name: "uniprot:Q0INS1", method: "offset", first: 1 });
    b.data.mapping({
      ref: "sifts",
      from: uni,
      to: "ab/entity:1",
      method: "by_numbering",
      numbering: "uniprot:Q0INS1",
    } as never);
    b.view({ ref: "v", label: "V" }).track({ label: "" }).ruler();
    const loaders = { load: async () => bytes() };
    const c = await compile(b.getDocument(), registries(), loaders as never);
    expect(c.issues.filter((i) => i.severity === "error")).toEqual([]);
    set = c.set;
    const structure = await structureOf(await readStructure(bytes(), "mmcif"), 0, null);
    link = StructureLink.attach(new Translator(), "3d", structure, set, "ab");
  });
  const run = (text: string, structure?: SeqSelection) => {
    const bound = bindQuery(readQuery(text).node as never, env.language);
    if (!bound.expr) throw new Error(bound.issues.map((i) => renderIssue(i)).join("\n"));
    const r = compileQuery(bound.expr, env, bound.source).run({ set, inputs: structure ? { structure } : {} });
    expect(r.issues.map((i) => renderIssue(i).split("\n")[0])).toEqual([]);
    return r.selection as SeqSelection;
  };
  /** The residues a loci holds, as auth number and code. */
  const residues = (loci: StructureElement.Loci) => {
    const out: string[] = [];
    StructureElement.Loci.forEachLocation(loci, (l) => {
      out.push(`${P.residue.auth_seq_id(l)}${P.residue.pdbx_PDB_ins_code(l)}`);
    });
    return out;
  };

  it("a label with a code selects that residue in 3D, and Mol*'s residue comes back as the same label", () => {
    const sel = run('(sel "ab/H[52A]")');
    expect(SeqSelection.summary(sel, { maxParts: 4 })).toBe("1 group, 1 position: ab/H 52A");
    const loci = selectionLoci(link, sel);
    expect(residues(loci)).toEqual(["52A"]);
    expect(SeqSelection.summary(structureInput(link, loci), { maxParts: 4 })).toBe("1 group, 1 position: ab/H 52A");
    // A range by plain numbers holds the inserted residues between them, in index order.
    expect(residues(selectionLoci(link, run('(sel "ab/H[52-53]")')))).toEqual(["52", "52A", "52B", "53"]);
  });

  /** B-factors (10 per residue in file order) from atom_site as a values annotation, its columns renamed by `fields`. */
  async function bfactors(numbering: string, fields: Record<string, string>) {
    const b = createSeqViewSpec();
    b.data.complex({
      ref: "ab",
      provider: "molstar:structure",
      source: { uri: URI, format: "mmcif" },
      options: { unit_name: "auth_asym_id" },
    } as never);
    b.data.annotation({
      ref: "bf",
      geometry: "values",
      from_source: "ab",
      category_name: "atom_site",
      numbering,
      fields: { target: "auth_asym_id", value: "B_iso_or_equiv", ...fields },
    } as never);
    b.view({ ref: "v", label: "V" }).track({ label: "" }).ruler();
    const c = await compile(b.getDocument(), registries(), { load: async () => bytes() } as never);
    const a = c.set.annotations.find((x) => x && c.set.strings.get(x.key).endsWith("bf"));
    return { issues: c.issues.map((i) => i.code), start: Array.from(a?.start ?? []) };
  }
  it("a CIF category by label_seq_id puts each row on its own residue", async () => {
    expect(await bfactors("label", { pos: "label_seq_id" })).toEqual({ issues: [], start: [0, 1, 2, 3, 4, 5, 6, 7] });
  });
  it("by auth numbers, pdbx_PDB_ins_code as ins_code puts 52A's row on 52A ('?' none)", async () => {
    expect(await bfactors("auth", { pos: "auth_seq_id", ins_code: "pdbx_PDB_ins_code" })).toEqual({
      issues: [],
      start: [0, 1, 2, 3, 4, 5, 6, 7],
    });
    // By the number alone, 52A's and 52B's rows are 52's: what the code is for.
    expect((await bfactors("auth", { pos: "auth_seq_id" })).start).toEqual([0, 1, 2, 2, 2, 5, 6, 7]);
  });

  it("onto UniProt through SIFTS: numbers run across the inserted residues", () => {
    expect(SeqSelection.summary(run('(onto (sel "ab/H[52A-53]") uni)'), { maxParts: 4 })).toBe(
      "1 group, 3 positions: uni/unit 73-75",
    );
  });
});
