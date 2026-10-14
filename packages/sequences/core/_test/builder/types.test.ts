/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder's types (seqspec:format#builders): a handle stands for its ref only in a `ref<kind>` param of its kind,
// each layer type's method takes that type's params, and a mapping's params follow its `method`. These are compile-time
// checks, run by `tsc` (`mise run typecheck`): a `@ts-expect-error` that stops erroring is an error.
// The functions are never called; the one test only keeps vitest from reporting an empty file.

import { expect, it } from "vitest";
import { basePairs, createSeqViewSpec, translate } from "../../src/builder/index.ts";

function types() {
  const b = createSeqViewSpec();
  const seq = b.data.sequence({ ref: "s", value: "MKV" });
  const lig = b.data.ligand({ ref: "l", value: "ATP" } as never);
  const ann = b.data.annotation({ geometry: "ranges", target: seq, value: [{ start: 0, stop: 1 }] });
  const com = b.data.computed({ provider: "coverage" } as never);
  const map = b.data.mapping({ from: seq, to: seq, method: "identity" });
  const aln = b.data.alignment({ ref: "a", value: ["MKV", "MK-"] });
  const cpx = b.data.complex({ ref: "c" } as never);
  const v = b.view({ label: "V" });
  const t = v.track({});

  // A ref<kind> param takes a handle of one of its kinds, or a string.
  b.data.annotation({ geometry: "ranges", target: seq, value: [] });
  b.data.annotation({ geometry: "ranges", target: aln.row(0), value: [] });
  b.data.annotation({ geometry: "ranges", target: "s", value: [] });
  // @ts-expect-error an annotation handle is no target
  b.data.annotation({ geometry: "ranges", target: ann, value: [] });
  // @ts-expect-error a mapping handle is no target
  b.data.annotation({ geometry: "ranges", target: map, value: [] });
  // @ts-expect-error a source handle is no sequence
  v.block({ sequence: b.data.source({ host_file: "x", format: "fasta" }) });
  v.block({ sequence: seq });
  v.block({ unit: cpx.unit("A") });
  v.block({ alignment: aln });
  // @ts-expect-error an alignment is no unit
  v.block({ unit: aln });

  // A param that is not a ref takes no handle.
  // @ts-expect-error a label is a string
  b.data.sequence({ label: seq, value: "M" });
  // @ts-expect-error the same for a layer's opacity
  t.ranges({ opacity: seq });
  // @ts-expect-error a ref string is no handle for a number
  t.values({ domain: "0-1" });
  // @ts-expect-error a tuple keeps its length
  t.values({ domain: [0, 1, 2] });
  // @ts-expect-error the same in a mapping's ranges
  b.data.mapping({ from: seq, to: seq, method: "codon", cds: [[0, 3, 6]] });

  // Lists of refs.
  t.ranges({ data: ann });
  t.ranges({ data: [ann, com, "x"] });
  // @ts-expect-error a sequence is no data
  t.ranges({ data: [ann, seq] });
  t.letters({ of: [seq, "p"] });
  t.letters({ of: "mapped" });
  t.letters({ of: seq });
  // @ts-expect-error letters `of` takes polymers of sequence or unit kind: a source is none
  t.letters({ of: b.data.source({ host_file: "x", format: "fasta" }) });
  // @ts-expect-error nor is a ligand
  t.letters({ of: lig });

  // One method per kernel layer type, with that type's params.
  t.ranges({ display: "brackets", split: "category", labels: true });
  t.variants({ split: "alt" });
  t.points({ split: { by: "x", order: ["b"], empty: "hide" } });
  t.values({ encoding: "line", domain: [0, 1] });
  t.letters({ row: 1, residue_form: "three" });
  t.letters({ letter_map: "complement" });
  t.rows({ letter_map: { A: "T", T: "A" } });
  // @ts-expect-error a letter map is a name or a dict of letters, not a list
  t.letters({ letter_map: ["A"] });
  // @ts-expect-error ranges have none
  t.ranges({ letter_map: "complement" });
  t.rows({ order: "tree", row_subset: [0, 1] });
  t.pairs({ side: "above" });
  // @ts-expect-error a misspelled method
  t.rangez({});
  // @ts-expect-error `display` is ranges' and pairs', not letters'
  t.letters({ display: "lanes" });
  // @ts-expect-error letters' `row` is a number
  t.letters({ row: "1" });
  // @ts-expect-error a value `display` does not have
  t.ranges({ display: "bogus" });
  // @ts-expect-error only ranges, points and variants split
  t.letters({ split: "category" });
  // @ts-expect-error values has no `split`
  t.values({ split: "category" });
  // @ts-expect-error an unknown param
  t.ranges({ colour: "red" });

  // `layer`: a kernel type with its params, or an extension's.
  t.layer({ type: "ranges", split: "category" });
  t.layer({ type: "org.example:dots", fallback: "points", options: { size: 3 } });
  // @ts-expect-error a misspelled kernel type is not an extension either
  t.layer({ type: "rangez" });
  // @ts-expect-error the type's params follow `type`
  t.layer({ type: "letters", display: "lanes" });

  // Mappings: `method` chooses the case, `bind` takes the row and `to` apart from it.
  b.data.mapping({ from: seq, to: seq, method: "identity" });
  b.data.mapping({ from: seq, to: lig, method: "offset", offset: 2 });
  b.data.mapping({ from: seq, to: seq, method: "pairwise", gapped: ["M-", "MK"] });
  b.data.mapping({ from: seq, to: seq, method: "codon", cds: [[0, 3]], strand: "+" });
  // @ts-expect-error a codon mapping has no `offset`
  b.data.mapping({ from: seq, to: seq, method: "codon", offset: 1 });
  // @ts-expect-error an annotation is no side
  b.data.mapping({ from: ann, to: seq, method: "identity" });
  aln.bind(0, seq, { method: "align" });
  aln.bind(0, seq, { method: "align", aligner: "biomsa", options: { gap: 2 } });
  aln.bind(1, seq, { method: "identity", ref: "m" });
  aln.bind(1, "s", { method: "offset", offset: 1 });
  // @ts-expect-error `bind` takes `from` as the row
  aln.bind(0, seq, { method: "identity", from: seq });
  // @ts-expect-error not a method
  aln.bind(0, seq, { method: "bogus" });
  // @ts-expect-error align's aligner is a name
  aln.bind(0, seq, { method: "align", aligner: 3 });
  // @ts-expect-error `bind` needs a method
  aln.bind(0, seq, {});

  // A complex's unit by value takes `Of<"unit">`: no `source`.
  cpx.unit({ ref: "u", sequence: seq });
  // @ts-expect-error a unit has no inline source
  cpx.unit({ ref: "u", source: { host_file: "x", format: "cif" } });
  // `translate` reads a sequence by value, `basePairs` writes over a target.
  // @ts-expect-error translate's `from` is a sequence handle
  translate(b.data, { from: aln.row(0) });
  // @ts-expect-error an annotation is no target
  basePairs(b.data, { target: ann, structure: "()" });
  const named: string = cpx.unit("A").ref;

  // A section's style: a preset's name, or an inline preset of the curated keys (seqspec:format#sections).
  b.view({ label: "S", style: "msa" });
  b.view({ label: "S", style: { preset: "msa", column_heads: "own_row", data_floor: 0, row_borders: true } });
  v.section({ style: { vertical_nav: "both", wrapped_group_labels: "inline_start", row_height: 18, header_width: 0 } });
  // @ts-expect-error the format's values are snake_case: own_row
  v.section({ style: { column_heads: "ownRow" } });
  // @ts-expect-error only the curated keys
  v.section({ style: { dataFloor: 0 } });
  // @ts-expect-error px are numbers
  v.section({ style: { data_floor: "0" } });
  return [named, com];
}

it("compiles (checked by tsc)", () => {
  expect(types).toBeTypeOf("function");
});
