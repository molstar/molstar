/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1: streamed sources compile to their manifests (seqspec:sources#streamed-sources): only the manifest is asked;
// a letters stream's complex is lazy units over lazy sequences (3.1 G positions, nothing per position), a sequence takes
// its record's system, a values annotation is streamed over the complex's members with the manifest's domain, an
// alignment's shape is the manifest's (declared counts checked) and its named rows come as `rows` patches; a source
// with no provider is unavailable, a failed or malformed manifest is reported. And N alignments over N records
// of one Stockholm file parse it once.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { type Issue, KeyKind, type KeyTarget, type SeqSet, Stockholm } from "@molstar/seq-core/model";
import { describe, expect, it, vi } from "vitest";
import {
  compile,
  type DataPatch,
  type Loaders,
  resolveView,
  type SeqViewSpecDoc,
  type StreamedData,
  StreamProp,
  supportedAll,
  validate,
} from "../../src/format/index.js";
import { registries } from "./load-fixtures.js";
import { GRCH38, hg38Streams, lettersManifest, scripted } from "./stream-provider.js";

const reg = registries();
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
/** An entry of a set's list that the test knows is there. */
const got = <T>(xs: readonly (T | null)[], i: number) => xs[i] as T;
const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const find = (set: SeqSet, key: string) => (set.byKey(set.strings.find(key), t) ? { ...t } : null);
/** A document needs a view: one over `node`'s block. */
const viewed = (b: ReturnType<typeof createSeqViewSpec>, block: Record<string, string>) => {
  b.view({ ref: "v", label: "V" }).block(block).track({ label: "R" }).ruler();
  return b.getDocument();
};
const loaders = (stream: Loaders["stream"]): Loaders => ({
  load: (r) => Promise.reject(new Error(`no fixture ${r.uri}`)),
  stream,
});

function genomeDoc(o: { record?: string } = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "hg38" });
  const genome = b.data.source({ ref: "genome", stream: "genome", data_kind: "letters" });
  const hg38 = b.data.complex({ ref: "hg38", source_ref: genome, label: "GRCh38" });
  b.data.source({ ref: "phylop_src", stream: "phylop", data_kind: "values", coordinate_system: "hg38" });
  const phylop = b.data.annotation({ ref: "phylop", geometry: "values", source_ref: "phylop_src" });
  if (o.record) b.data.sequence({ ref: "chr21", source_ref: genome, record: o.record, alphabet: "dna" });
  const v = b.view({ ref: "v", label: "hg38" });
  v.block({ complex: hg38 });
  v.track({ label: "phyloP" }).values({ data: phylop });
  return b.getDocument();
}

describe("a streamed source compiles to its manifest", () => {
  it("validates: stream needs data_kind, a letters source takes no coordinate_system", () => {
    expect(validate(genomeDoc(), reg)).toEqual([]);
    const b = createSeqViewSpec();
    b.data.source({ ref: "s", stream: "x" });
    b.data.source({ ref: "c", stream: "y", data_kind: "letters", coordinate_system: "s" as never });
    const want = ["svs.params.exclusive", "svs.params.missing", "svs.ref.wrong_kind"];
    expect(codes(validate(viewed(b, { sequence: "x" }), reg))).toEqual(expect.arrayContaining(want));
  });

  it("a letters stream's complex: lazy units and sequences of the manifest's lengths, only the manifest asked", async () => {
    const s = hg38Streams();
    const { set, issues } = await compile(genomeDoc({ record: "chr21" }), reg, loaders(s.stream));
    expect(issues).toEqual([]);
    expect(s.genome.requests).toEqual([]);
    expect(s.phylop.requests).toEqual([]);
    expect(set.units.count).toBe(25);
    let total = 0;
    for (const [key, n] of GRCH38) {
      const u = find(set, `hg38/${key}`);
      expect(u?.kind).toBe(KeyKind.unit);
      const m = got(set.molecules, set.units.molecule[u?.id ?? -1]);
      expect([m.lazy, m.length, m.code.length, m.alphabet]).toEqual([true, n, 0, "dna"]);
      total += m.length;
    }
    expect(total).toBe(3_088_269_832);
    const chr21 = got(set.molecules, (find(set, "chr21") as KeyTarget).id);
    expect([chr21.lazy, chr21.length]).toEqual([true, 46_709_983]);
  });

  it("a values annotation: streamed over the complex's members, the manifest's domain the layer's", async () => {
    const s = hg38Streams();
    const doc = genomeDoc();
    const { set } = await compile(doc, reg, loaders(s.stream));
    const a = got(set.annotations, (find(set, "phylop") as KeyTarget).id);
    expect([a.streamed, a.count]).toEqual([true, 0]);
    const data = a.custom.props.get(StreamProp) as StreamedData;
    expect(data.systems.size).toBe(24);
    expect(data.systems.get("hg38/chrX")).toBe(22);
    const { view, issues } = resolveView(doc, set, "v", reg, supportedAll(reg));
    expect(issues.filter((i) => i.severity === "error")).toEqual([]);
    const k = Array.from(view.layers.type).findIndex((x) => set.strings.get(x) === "values");
    expect([view.layers.domainMin[k], view.layers.domainMax[k]]).toEqual([0, 100]);
    const blocks = new Set<number>();
    for (let i = 0; i < view.placements.count; i++)
      if (view.placements.layer[i] === k) blocks.add(view.placements.block[i]);
    expect(blocks.size).toBe(24);
  });

  it("a manifest key naming no member is left out (svs.stream.ref); no provider: unavailable; a failed manifest", async () => {
    const extra = [...GRCH38, ["chrM", 16_569] as const];
    const s = hg38Streams({ chromosomes: extra });
    const genome = scripted(lettersManifest(), () => new Error("none"));
    const stream = (src: { stream: string | null }) => (src.stream === "genome" ? genome.provider : s.stream(src));
    const r = await compile(genomeDoc(), reg, loaders(stream));
    expect(codes(r.issues)).toEqual(["svs.stream.ref"]);
    expect(r.issues[0].message).toContain("chrM");
    const none = await compile(
      genomeDoc(),
      reg,
      loaders(() => null),
    );
    expect(codes(none.issues)).toEqual(["svs.provider.unavailable"]);
    const bad = await compile(
      genomeDoc(),
      reg,
      loaders((src) =>
        src.stream === "phylop"
          ? { fetch: () => Promise.reject(new Error("offline")) }
          : { fetch: () => Promise.resolve({ kind: "manifest", manifest: { levels: [] } } as unknown as DataPatch) },
      ),
    );
    expect(codes(bad.issues)).toEqual(["svs.source.format", "svs.source.load"]);
  });

  it("a streamed alignment: the manifest's shape, the declared counts checked, its reference row a rows patch", async () => {
    const manifest = {
      systems: [{ key: "aln", extent: [0, 300] as [number, number], numbering: null, label: null }],
      levels: [[1, 1] as [number, number]],
      rows: 5000,
      columns: 300,
    };
    const s = scripted(manifest, (r) => ({
      ...r,
      kind: "rows",
      labels: [`row${r.row_start}`],
      rows: ["A".repeat(300)],
    }));
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: "msa", stream: "msa", data_kind: "alignment" });
    b.data.alignment({ ref: "aln", source_ref: src, row_count: 4999, reference_row: 7 });
    const { set, issues } = await compile(
      viewed(b, { alignment: "aln" }),
      reg,
      loaders(() => s.provider),
    );
    expect(codes(issues)).toEqual(["svs.alignment.shape"]);
    const a = got(set.alignments, 0);
    expect([a.lazy, a.rowCount, a.columnCount]).toEqual([true, 5000, 300]);
    expect(s.requests.map((r) => [r.row_start, r.row_stop, r.level, r.labels])).toEqual([[7, 8, -1, false]]);
    expect((a.custom.props.get(StreamProp) as StreamedData).labels.get(7)).toBe("row7");
  });

  it("a reader of another data kind is svs.source.format", async () => {
    const s = hg38Streams();
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: "phylop_src", stream: "phylop", data_kind: "values" });
    b.data.sequence({ ref: "seq", source_ref: src, alphabet: "dna" });
    const r = await compile(viewed(b, { sequence: "seq" }), reg, loaders(s.stream));
    expect(codes(r.issues)).toEqual(["svs.source.format"]);
  });
});

describe("one parse per source for the alignments of its records", () => {
  it("three alignment nodes over three records of one Stockholm file parse it once", async () => {
    const one = (k: number) => `# STOCKHOLM 1.0\nA${k}/1-4 ACDE\nB${k}/1-4 AC-E\n//\n`;
    const file = [0, 1, 2].map(one).join("");
    const parse = vi.spyOn(Stockholm, "parse");
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: "f", value: file, format: "stockholm" });
    for (const k of [0, 1, 2]) b.data.alignment({ ref: `a${k}`, source_ref: src, record: k });
    const { set, issues } = await compile(viewed(b, { alignment: "a0" }), reg, loaders(undefined));
    expect(issues).toEqual([]);
    expect(set.alignments.map((a) => a?.labelText(0 as never))).toEqual(["A0/1-4", "A1/1-4", "A2/1-4"]);
    expect(parse).toHaveBeenCalledTimes(1);
    parse.mockRestore();
  });
});
