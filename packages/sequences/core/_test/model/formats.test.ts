/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import { Clustal, Fasta, type Issue, Newick, Stockholm } from "../../src/model/index.js";
import { RNA_STOCKHOLM, SH2_A3M, SH2_STOCKHOLM } from "./alignment-fixtures.js";

const codes = (issues: Issue[]) => issues.map((i) => `${i.severity} ${i.code}`);

describe("FASTA", () => {
  it("takes the record id from the header's first word and joins sequence lines", () => {
    const { records, issues } = Fasta.parse(">sp|P00648|RNBR_BACAM Barnase  precursor\nAQVIN\n TFDG \r\n>b\n\nKKA\n");
    expect(issues).toEqual([]);
    expect(records).toEqual([
      { id: "sp|P00648|RNBR_BACAM", description: "Barnase  precursor", sequence: "AQVINTFDG" },
      { id: "b", description: "", sequence: "KKA" },
    ]);
  });

  it("reads A3M with its comment line, keeping case for the alignment", () => {
    const { records, issues } = Fasta.parse(SH2_A3M);
    expect(issues).toEqual([]);
    expect(records.map((r) => r.id)).toEqual([
      "SRC_HUMAN/151-180",
      "PTN11_HUMAN/6-36",
      "GRB2_HUMAN/60-90",
      "P85A_HUMAN/333-360",
    ]);
    expect(records[0].sequence).toBe("WYFGKITRRESERLLLNAENPRGTFLVRES");
    const { input, issues: skipped } = Fasta.alignment(records, "sh2", "a3m");
    expect(skipped).toEqual([]);
    expect([input.flavor, input.labels?.[2], input.gapped[2]]).toEqual([
      "a3m",
      "GRB2_HUMAN/60-90",
      "WFFGKIPRAKeAEEMLSKQRHgDGAFLIRES-",
    ]);
  });

  it("reports a sequence before a header, duplicate ids, empty input", () => {
    // An unreadable file is model.source.format, an error; a missing or repeated id model.source.id, a warning: one
    // code, one severity (GV-05).
    expect(codes(Fasta.parse("ACGT\n>a\nAC\n>a\nGT").issues)).toEqual([
      "error model.source.format",
      "warning model.source.id",
    ]);
    expect(codes(Fasta.parse("").issues)).toEqual(["error model.source.format"]);
    expect(codes(Fasta.parse("> \nAC").issues)).toEqual(["warning model.source.id"]);
  });

  it("reads LF, CRLF and CR line ends alike", () => {
    const lf = ">a one\nAC\nGT\n>b\nTT\n";
    const want = Fasta.parse(lf);
    expect(want.records.map((r) => [r.id, r.description, r.sequence])).toEqual([
      ["a", "one", "ACGT"],
      ["b", "", "TT"],
    ]);
    expect(Fasta.parse(lf.replaceAll("\n", "\r\n"))).toEqual(want);
    expect(Fasta.parse(lf.replaceAll("\n", "\r"))).toEqual(want);
    // A classic-Mac file used to be one header, its description swallowing the rest.
    expect(Fasta.parse(">a\r>b\rACGT").records.map((r) => [r.id, r.sequence])).toEqual([
      ["a", ""],
      ["b", "ACGT"],
    ]);
  });

  it("keeps a header with no sequence, so record indices hold, and notes it with its line", () => {
    const { records, issues } = Fasta.parse(">a\n>b\nACGT");
    expect(records.map((r) => [r.id, r.sequence])).toEqual([
      ["a", ""],
      ["b", "ACGT"],
    ]);
    expect(issues).toEqual([
      {
        code: "model.source.empty_record",
        severity: "info",
        path: "",
        message: 'FASTA: 1 record(s) with no sequence: "a" (line 1).',
      },
    ]);
    // At the end, with blank lines, and with a path.
    const last = Fasta.parse("\n>a\nAC\n\n>b\n\n", { path: "/nodes/1/params/value" });
    expect(last.records.map((r) => r.id)).toEqual(["a", "b"]);
    expect(last.issues.map((i) => [i.code, i.path, i.message])).toEqual([
      ["model.source.empty_record", "/nodes/1/params/value", 'FASTA: 1 record(s) with no sequence: "b" (line 5).'],
    ]);
    expect(codes(Fasta.parse(">a\n;only a comment\n>b\nAC").issues)).toEqual(["info model.source.empty_record"]);
  });

  it("says each repeated cause once for the file, with its count and first examples (GV-07)", () => {
    const hits = Array.from({ length: 300 }, (_, k) => `>id${k % 30}\nAC`).join("\n");
    expect(Fasta.parse(hits).issues.map((i) => i.message)).toEqual([
      'FASTA: 30 id(s) name several records: "id0", "id1", "id2", "id3", "id4" and 25 more.',
    ]);
    const empties = Array.from({ length: 8 }, (_, k) => `>e${k}`).join("\n");
    expect(Fasta.parse(`${empties}\n>x\nAC`).issues.map((i) => i.message)).toEqual([
      'FASTA: 8 record(s) with no sequence: "e0" (line 1), "e1" (line 2), "e2" (line 3), "e3" (line 4), "e4" (line 5) and 3 more.',
    ]);
    expect(Fasta.parse("AC\nGT\n>a\nAC").issues.map((i) => i.message)).toEqual([
      "FASTA: a sequence before the first header (line 1, 2).",
    ]);
  });

  it("skips `;` and `#` lines anywhere, not only before the first header", () => {
    expect(Fasta.parse(">a\nAC\n;note\n#x\nGT").records[0].sequence).toBe("ACGT");
    const text = ";legacy\n#also\n>a one\n;comment\nAC\n  ;indented\nGT\n>b\n#c\nTT\n";
    const { records, issues } = Fasta.parse(text);
    expect(issues).toEqual([]);
    expect(records.map((r) => r.sequence)).toEqual(["ACGT", "TT"]);
    expect(codes(Fasta.parse("ACGT\n>a\nAC").issues)).toEqual(["error model.source.format"]);
  });

  it("leaves HH-suite's ss_*, sa_* and Consensus records out of an A3M, with one note", () => {
    // addss.pl writes ss_dssp and sa_dssp (DSSP's solvent accessibility, letters A to F and -) besides ss_pred.
    const text =
      ">ss_pred PSIPRED predicted secondary structure\nCCHHH\n>ss_conf\n98765\n>sa_dssp\nABC-F\n>Consensus\nWYFGK\n" +
      ">q\nWYFGK\n>r\nWY-GK\n";
    const { records, issues } = Fasta.parse(text);
    expect(issues).toEqual([]);
    expect(records).toHaveLength(6);
    const a3m = Fasta.alignment(records, "hh", "a3m", { path: "/nodes/3/params/value" });
    expect(a3m.input.labels).toEqual(["q", "r"]);
    expect(a3m.input.gapped).toEqual(["WYFGK", "WY-GK"]);
    expect(a3m.issues.map((i) => [i.severity, i.code, i.path])).toEqual([
      ["info", "model.source.skipped", "/nodes/3/params/value"],
    ]);
    expect(a3m.issues[0].message).toContain("4 HH-suite record");
    // Only an A3M is HH-suite's; other records are rows, and a clean file has no warning.
    expect(Fasta.alignment(records, "hh", "aligned").input.labels).toEqual([
      "ss_pred",
      "ss_conf",
      "sa_dssp",
      "Consensus",
      "q",
      "r",
    ]);
    expect(Fasta.alignment(records.slice(4), "hh", "a3m").issues).toEqual([]);
  });

  it("bounds the size in UTF-8 bytes", () => {
    const text = ">a é\nACGT\n"; // 10 characters, 11 bytes
    expect(Fasta.parse(text, { maxBytes: 11 }).issues).toEqual([]);
    const over = Fasta.parse(text, { maxBytes: 10, path: "/nodes/2/params/value" });
    expect(over.records).toEqual([]);
    expect(over.issues.map((i) => [i.code, i.path])).toEqual([["model.source.size_bound", "/nodes/2/params/value"]]);
    expect(codes(Fasta.parse(">😀\nA", { maxBytes: 6 }).issues)).toEqual(["error model.source.size_bound"]);
    expect(Fasta.parse(">😀\nA", { maxBytes: 7 }).issues).toEqual([]);
  });
});

describe("Stockholm", () => {
  it("joins interleaved blocks and keeps GF, GS, GR and GC tags", () => {
    const { records, issues } = Stockholm.parse(SH2_STOCKHOLM);
    expect(issues).toEqual([]);
    const [r] = records;
    expect(r.names).toEqual(["SRC_HUMAN/151-180", "PTN11_HUMAN/6-36", "GRB2_HUMAN/60-90", "P85A_HUMAN/333-360"]);
    expect(r.rows[1]).toBe("WFHPNITGVEkdAENLLLTRGV.DGSFLARPS-");
    expect(r.rows.every((row) => row.length === 33)).toBe(true);
    expect(r.gf.get("CC")).toBe("An excerpt in the seed's shape,\nwritten by hand.");
    expect(r.gs.get("SRC_HUMAN/151-180")?.get("AC")).toBe("P12931.3");
    expect(r.gr.get("PTN11_HUMAN/6-36")?.get("SS")).toBe("-EEEE---HHhhHHHHH------EEEEE-----");
    expect(r.gc.get("RF")).toBe("xxxxxxxxxx..xxxxxxxxxx.xxxxxxxxxx");
    const input = Stockholm.alignment(r, "sh2");
    expect(Array.from(input.insertColumns ?? []).flatMap((x, c) => (x ? [c] : []))).toEqual([10, 11, 22]);
  });

  it("reads several records; WUSS pairs with brackets and letter pseudoknots", () => {
    const { records, issues } = Stockholm.parse(RNA_STOCKHOLM);
    expect(issues).toEqual([]);
    expect(records.map((r) => r.names)).toEqual([["seq1", "seq2"], ["seq3"]]);
    const ss = records[0].gc.get("SS_cons") as string;
    const pairs = Stockholm.pairs(ss, [], "") as Int32Array;
    expect(Array.from(pairs)).toEqual([12, 11, 10, 9, -1, 8, -1, -1, 5, 3, 2, 1, 0]);
    const bad: Issue[] = [];
    expect(Stockholm.pairs("<<.>", bad, "/x")).toBeNull();
    expect(Stockholm.pairs("<.>>", bad, "/x")).toBeNull();
    expect(Stockholm.pairs("A..b", bad, "/x")).toBeNull();
    expect(codes(bad)).toEqual(Array(3).fill("error model.source.format"));
  });

  it("reports tags of the wrong length, a GR without a row, a missing `//`, junk lines", () => {
    const text = "# STOCKHOLM 1.0\na ACGU\n#=GC SS_cons <<>\n#=GR b SS ....\nnot a row line\na AC";
    const { records, issues } = Stockholm.parse(text);
    expect(records[0].rows).toEqual(["ACGUAC"]);
    expect(records[0].gc.has("SS_cons")).toBe(false);
    expect(issues.map((i) => i.message)).toEqual([
      "Stockholm line 5: expected a row name and its residues.",
      "Stockholm: the last record does not end with `//`.",
      "Stockholm record 0: #=GC SS_cons has 3 columns, not 6.",
      'Stockholm record 0: #=GR for "b", which has no row.',
      "Stockholm record 0: #=GR b SS has 4 columns, not 6.",
    ]);
    expect(codes(Stockholm.parse("# STOCKHOLM 1.0\na AC").issues)).toEqual(["info model.source.unterminated"]);
    expect(codes(Stockholm.parse("//\n").issues)).toEqual(["error model.source.format"]);
    expect(codes(Stockholm.parse("# STOCKHOLM 1.0\n//").issues)).toEqual(["error model.source.format"]);
  });
});

describe("Clustal", () => {
  // PF00042's seed, its first eight rows, in CLUSTAL W 2.1's layout (fixtures README): two blocks of 60 and 58
  // columns, each line's residue count, a conservation line under each block.
  const seed = readFileSync(new URL("../fixtures/pf00042-seed-8.aln", import.meta.url), "utf8");
  const globins = JSON.parse(readFileSync(new URL("../fixtures/pf00042.json", import.meta.url), "utf8")).rows as [
    string,
    number,
    string,
  ][];

  it("reads the header, the blocks joined per name in order, counts and conservation lines left out", () => {
    const { record, issues } = Clustal.parse(seed);
    expect(issues).toEqual([]);
    expect(record?.header).toBe("CLUSTAL W (2.1) multiple sequence alignment");
    expect(record?.names).toEqual(
      globins.slice(0, 8).map(([label, start, g]) => {
        const n = g.replace(/[-.]/g, "").length;
        return `${label}/${start}-${start + n - 1}`;
      }),
    );
    expect(record?.rows).toEqual(globins.slice(0, 8).map(([, , g]) => g.replace(/\./g, "-")));
    const input = record ? Clustal.alignment(record, "globins") : null;
    expect([input?.flavor, input?.labels?.[0], input?.gapped.length]).toEqual(["aligned", "HBA_HUMAN/27-137", 8]);
  });

  it("reads Clustal Omega's and MUSCLE's headers, no counts, CRLF, and residues split by spaces", () => {
    const text =
      "CLUSTAL O(1.2.4) multiple sequence alignment\r\n\r\na  MK-V\r\nb  MKLV\r\n   ** *\r\n\r\na  LA AG\r\nb  LA-G\r\n";
    expect(Clustal.parse(text)).toEqual({
      record: {
        header: "CLUSTAL O(1.2.4) multiple sequence alignment",
        names: ["a", "b"],
        rows: ["MK-VLAAG", "MKLVLA-G"],
      },
      issues: [],
    });
    expect(Clustal.parse("MUSCLE (3.8) multiple sequence alignment\n\na MK\nb MR\n").record?.rows).toEqual([
      "MK",
      "MR",
    ]);
  });

  it("an error saying what to write: no header, a block naming other rows, rows of other lengths", () => {
    const err = (text: string) => {
      const { record, issues } = Clustal.parse(text);
      expect(record).toBeNull();
      return issues.map((i) => `${i.severity} ${i.code}: ${i.message}`);
    };
    expect(err(">a\nMK\n")).toEqual([
      'error model.source.format: Clustal: line 1 is not a Clustal header; the file starts with "CLUSTAL" (or "MUSCLE", "PROBCONS"), e.g. "CLUSTAL W (2.1) multiple sequence alignment".',
    ]);
    expect(err("CLUSTAL W\n\na MK\nb MR\n\nb LA\na LA\n")).toEqual([
      'error model.source.format: Clustal: line 6 names "b" where the first block has "a": every block lists the rows in the first block\'s order.',
    ]);
    expect(err("CLUSTAL W\n\na MK\nb MR\n\na LA\n")).toEqual([
      "error model.source.format: Clustal: the block at line 6 has 1 row(s); the first block has 2: every block lists every row.",
    ]);
    expect(err("CLUSTAL W\n\na MKV\nb MR\n")).toEqual([
      'error model.source.format: Clustal: the block at line 3 is ragged: "a" has 3 columns there, "b" 2; pad rows with "-".',
    ]);
    expect(err("CLUSTAL W\n\n")).toEqual(["error model.source.format: Clustal: no rows after the header."]);
  });

  it("bounds its size as the other text readers", () => {
    const { record, issues } = Clustal.parse(seed, { maxBytes: 100 });
    expect([record, issues.map((i) => i.code)]).toEqual([null, ["model.source.size_bound"]]);
  });
});

describe("line ends of Stockholm", () => {
  it("reads CRLF and CR as LF", () => {
    const want = Stockholm.parse(SH2_STOCKHOLM);
    expect(want.issues).toEqual([]);
    expect(Stockholm.parse(SH2_STOCKHOLM.replaceAll("\n", "\r\n"))).toEqual(want);
    expect(Stockholm.parse(SH2_STOCKHOLM.replaceAll("\n", "\r"))).toEqual(want);
  });
});

describe("Newick", () => {
  it("parses nodes in preorder with labels and branch lengths", () => {
    const { tree, issues } = Newick.parse("((A:0.1,B:2e-1)90:0.05,(C,'D x''s':.4)[a comment]:1)root;\n");
    expect(issues).toEqual([]);
    expect(tree?.label).toEqual(["root", "90", "A", "B", null, "C", "D x's"]);
    expect(Array.from(tree?.parent ?? [])).toEqual([-1, 0, 1, 1, 0, 4, 4]);
    expect(Array.from(tree?.branchLength ?? []).map((x) => (Number.isNaN(x) ? null : +x.toFixed(3)))).toEqual([
      null,
      0.05,
      0.1,
      0.2,
      1,
      null,
      0.4,
    ]);
    // `90` is the inner node's label (a support value), kept as written.
    expect(Newick.parse("((A,B)90);").tree?.label[1]).toBe("90");
  });

  it("reports unbalanced parentheses, two labels, bad lengths, trailing text, a missing `;`", () => {
    const first = (t: string) => Newick.parse(t).issues.map((i) => `${i.severity}: ${i.message}`)[0];
    expect(first("((A,B);")).toMatch(/^error: .*not balanced/);
    expect(first("(A,B));")).toMatch(/^error: .*outside parentheses/);
    expect(first("(A B,C);")).toMatch(/two labels/);
    expect(first("(A:x,B);")).toMatch(/not a number/);
    expect(first("(A,B);(C,D);")).toMatch(/after the tree/);
    expect(first("(A,'B);")).toMatch(/not closed/);
    expect(first("(A,B)")).toMatch(/^info: .*does not end with `;`/);
    expect(Newick.parse("(A,B)").issues[0].code).toBe("model.source.unterminated");
    expect(Newick.parse("(A,B)", { maxBytes: 2 }).tree).toBeNull();
  });

  it("reports empty or comment-only input, a group after a group, and two lengths on one node", () => {
    // Each was a tree (one node, or merged children, or the last length) and at most a warning.
    for (const [text, message] of [
      ["", /no tree/],
      ["  \n", /no tree/],
      ["[c]", /no tree/],
      ["[a][b] ", /no tree/],
      ["(a,b)(c,d);", /group follows a group/],
      ["(a,b)x(c,d);", /group follows a group/],
      ["(a,b):1(c,d);", /group follows a group/],
      ["a(b,c);", /group follows a group/],
      ["(a:1:2);", /two lengths/],
      ["(a,b):1:2;", /two lengths/],
    ] as [string, RegExp][]) {
      const { tree, issues } = Newick.parse(text, { path: "/nodes/2/params/value" });
      expect(tree, text).toBeNull();
      expect(
        issues.map((i) => [i.severity, i.code, i.path]),
        text,
      ).toEqual([["error", "model.source.format", "/nodes/2/params/value"]]);
      expect(issues[0].message, text).toMatch(message);
    }
  });

  it("still reads the nested and one-node trees that are trees", () => {
    expect(Newick.parse("((a,b),(c,d));").issues).toEqual([]);
    expect(Newick.parse("((a,b)x:1,c:2)root:3;").tree?.label).toEqual(["root", "x", "a", "b", "c"]);
    expect(Array.from(Newick.parse(";").tree?.parent ?? [])).toEqual([-1]);
    expect(Newick.parse("A;").tree?.label).toEqual(["A"]);
    expect(Newick.parse("(A,B)").issues.map((i) => i.severity)).toEqual(["info"]);
  });
});
