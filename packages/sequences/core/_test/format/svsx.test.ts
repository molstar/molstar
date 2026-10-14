/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SX1's tests (seqspec:sources#archives, seqspec:sources#data-by-reference): Mol*'s zip
// through its shims, the archive's main file, relative URIs against its directory, and the loaders over an archive's
// files. Node only.

import { readFileSync } from "node:fs";
import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { writeSvsx } from "@molstar/seq-core/builder/svsx";
import type { Issue, KeyTarget, SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  archiveLoaders,
  compile,
  type Loaders,
  type LoadRequest,
  loadStreamNode,
  readSvsx,
  type StageContext,
  type StreamSource,
  snapshotDoc,
  validate,
} from "../../src/format/index.ts";
import { isRelativeUri, resolveArchivePath } from "../../src/format/uri-path.ts";
import { RuntimeContext } from "../../src/format/zip/shims.ts";
import { deflate, inflateRaw, unzip, zip } from "../../src/format/zip/zip.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";

const enc = (s: string) => new TextEncoder().encode(s);
const dec = (b: Uint8Array) => new TextDecoder().decode(b);
const reg = registries(false);
const none: Loaders = {
  async load(r: LoadRequest) {
    throw new Error(`fallback: ${r.uri}`);
  },
};
const fixture = () => new Uint8Array(readFileSync(new URL("fixtures/kinase.svsx", import.meta.url)));

/** A document with one FASTA source at `uri` read by a sequence in a view. */
function docAt(uri: string, label = "x") {
  const b = createSeqViewSpec();
  const src = b.data.source({ ref: "fa", uri, format: "fasta" });
  const seq = b.data.sequence({ ref: "s", source_ref: src, label });
  b.view({ ref: "v", label: "V" }).block({ sequence: seq }).track({ label: "L" }).letters();
  return b.getDocument();
}

describe("Mol*'s zip, through its shims", () => {
  const sample = (n: number) => Uint8Array.from({ length: n }, (_, i) => (i * 7 + (i >> 5)) % 251);
  const files = { "index.svsm": sample(5000), "data/a.txt": enc("ABC".repeat(4000)), "data/empty": new Uint8Array(0) };

  it("zips and unzips every file as it was, compressing the repetitive one", async () => {
    const bytes = new Uint8Array(await zip(RuntimeContext, files as never));
    expect(bytes.length).toBeLessThan(5000 + 12000);
    const out = await unzip(RuntimeContext, bytes.slice().buffer);
    expect(Object.keys(out).sort()).toEqual(Object.keys(files).sort());
    for (const [k, v] of Object.entries(files)) expect(out[k]).toEqual(v);
  });

  it("deflates and inflates (zlib wrapper) back to the input", async () => {
    const data = enc("the same words, the same words, the same words ".repeat(200));
    const packed = await deflate(RuntimeContext, data as never);
    expect(packed.length).toBeLessThan(data.length / 4);
    // `deflate` writes a zlib stream: a 2-byte header and a 4-byte checksum around the raw deflate data.
    expect(await inflateRaw(RuntimeContext, packed.subarray(2, packed.length - 4) as never)).toEqual(data);
  });
});

describe("relative URIs", () => {
  it("resolve against the main file's directory; `.` and empty segments drop", () => {
    expect(resolveArchivePath("data/abl1.fa", "")).toBe("data/abl1.fa");
    expect(resolveArchivePath("./tiles/aln/", "")).toBe("tiles/aln");
    expect(resolveArchivePath("../data/x.fa", "docs")).toBe("data/x.fa");
    expect(resolveArchivePath("x.fa", "docs/sub")).toBe("docs/sub/x.fa");
  });

  it("`..` past the archive's root is refused, and so is anything not a relative path", () => {
    for (const [uri, base] of [
      ["../x.fa", ""],
      ["../../x.fa", "docs"],
      ["a/../../x", ""],
      ["/etc/x", ""],
      ["a\\b", ""],
      ["", ""],
      ["./", ""],
    ]) {
      expect(resolveArchivePath(uri, base)).toBeNull();
    }
    expect(["https://a/b", "http://a/b", "file:///x", "//host/x", "/x", ""].map(isRelativeUri)).toEqual(
      Array(6).fill(false),
    );
  });

  it("validate takes them only with a base: outside an archive a relative uri is svs.source.uri", () => {
    const doc = docAt("data/a.fa");
    expect(validate(doc, reg).map((i) => i.code)).toEqual(["svs.source.uri"]);
    expect(validate(doc, reg, { base: "" })).toEqual([]);
    expect(validate(docAt("../a.fa"), reg, { base: "" }).map((i) => i.code)).toEqual(["svs.source.uri"]);
    expect(validate(docAt("../a.fa"), reg, { base: "docs" })).toEqual([]);
    expect(validate(docAt("http://example.org/a.fa"), reg, { base: "" }).map((i) => i.code)).toEqual([
      "svs.source.uri",
    ]);
  });
});

describe("readSvsx", () => {
  it("opens the committed fixture: its document validated, its files by path", async () => {
    const r = await readSvsx(fixture());
    expect(r.issues).toEqual([]);
    expect([r.main, r.base]).toEqual(["index.svsm", ""]);
    expect([...r.files.keys()].sort()).toEqual([
      "data/example.org/kinase/abl.fa",
      "data/example.org/kinase/sites.tsv",
      "index.svsm",
    ]);
    expect(dec(r.files.get("data/example.org/kinase/sites.tsv") as Uint8Array)).toContain("gatekeeper");
  });

  it("compiles over archiveLoaders with no network: a failing fallback is never asked", async () => {
    const r = await readSvsx(fixture());
    const { set, issues } = await compile(r.doc as never, registries(), archiveLoaders(r.files, none, r.base));
    expect(issues).toEqual([]);
    const t: KeyTarget = { kind: -1, id: -1, row: -1 };
    const find = (s: SeqSet, key: string) => s.byKey(s.strings.find(key), t);
    expect(find(set, "abl")).toBe(true);
    expect(set.molecules[t.id]?.length).toBe(106);
    expect(find(set, "sites")).toBe(true);
    expect(set.annotations[t.id]?.count).toBe(3);
  });

  it("finds index.svsm before index.svsj, and index.svsj when it is alone", async () => {
    const fa = { "a.fa": enc(">s\nMKV\n") };
    const [m, j] = [docAt("a.fa", "from msgpack"), docAt("a.fa", "from json")];
    const msgpack = (await readSvsx(await writeSvsx(m, fa))).files.get("index.svsm") as Uint8Array;
    // The JSON file is listed first in the zip: the reader looks by name, not by order.
    const both = new Uint8Array(
      await zip(RuntimeContext, { "index.svsj": enc(JSON.stringify(j)), "index.svsm": msgpack, ...fa } as never),
    );
    const r = await readSvsx(both);
    expect([r.main, r.issues]).toEqual(["index.svsm", []]);
    expect(JSON.stringify(r.doc)).toContain("from msgpack");
    expect(JSON.stringify((await readSvsx(both, { main: "index.svsj" })).doc)).toContain("from json");
    const alone = await readSvsx(await writeSvsx(j, fa, { format: "json" }));
    expect([alone.main, alone.issues]).toEqual(["index.svsj", []]);
    expect(JSON.stringify(alone.doc)).toContain("from json");
    expect([...(await readSvsx(await writeSvsx(m, fa))).files.keys()]).toEqual(["index.svsm", "a.fa"]);
  });

  it("reports an archive with no main file, a missing named one, and bytes that are no zip", async () => {
    const empty = new Uint8Array(await zip(RuntimeContext, { "a.fa": enc(">s\nM\n") as never }));
    expect((await readSvsx(empty)).issues.map((i) => [i.code, i.severity])).toEqual([["svs.format.invalid", "error"]]);
    expect((await readSvsx(fixture(), { main: "other.svsj" })).issues[0].message).toContain("other.svsj");
    const r = await readSvsx(enc("not a zip"));
    expect([r.doc, r.issues[0].message]).toEqual([null, expect.stringContaining("the archive cannot be read")]);
  });

  it("a main file below the root resolves its relative uris from there; `..` reaches the data beside it", async () => {
    const doc = docAt("../data/a.fa");
    const bytes = await writeSvsx(doc, { "data/a.fa": enc(">s\nMKV\n") });
    // The writer's main file is at the root; the same files with the document under `docs/` read with `main`.
    const files = (await readSvsx(bytes)).files;
    const nested = new Uint8Array(
      await zip(RuntimeContext, {
        "docs/index.svsm": files.get("index.svsm") as never,
        "data/a.fa": files.get("data/a.fa") as never,
      }),
    );
    const r = await readSvsx(nested, { main: "docs/index.svsm" });
    expect([r.base, r.issues]).toEqual(["docs", []]);
    expect(
      dec(
        await archiveLoaders(r.files, none, r.base).load(
          { uri: "../data/a.fa", hostFile: null, maxBytes: null },
          new AbortController().signal,
        ),
      ),
    ).toBe(">s\nMKV\n");
    expect((await readSvsx(bytes, { main: "index.svsm" })).issues.map((i) => i.code)).toEqual(["svs.source.uri"]);
  });
});

describe("archiveLoaders", () => {
  const files = new Map([["data/a.fa", enc(">s\nMKV\n")]]);
  const signal = new AbortController().signal;
  const req = (uri: string | null, hostFile: string | null = null): LoadRequest => ({ uri, hostFile, maxBytes: null });

  it("answers a relative uri from the archive and passes everything else to the fallback", async () => {
    const asked: (string | null)[] = [];
    const fallback: Loaders = {
      async load(r) {
        asked.push(r.uri ?? r.hostFile);
        return enc("remote");
      },
    };
    const l = archiveLoaders(files, fallback);
    expect(dec(await l.load(req("data/a.fa"), signal))).toBe(">s\nMKV\n");
    expect(dec(await l.load(req("./data/a.fa"), signal))).toBe(">s\nMKV\n");
    expect(dec(await l.load(req("https://example.org/a.fa"), signal))).toBe("remote");
    expect(dec(await l.load(req(null, "down.fasta"), signal))).toBe("remote");
    expect(asked).toEqual(["https://example.org/a.fa", "down.fasta"]);
  });

  it("fails a relative uri the archive lacks, or that climbs out of it", async () => {
    const l = archiveLoaders(files, none);
    await expect(l.load(req("data/missing.fa"), signal)).rejects.toThrow("not in the archive");
    await expect(l.load(req("../a.fa"), signal)).rejects.toThrow("outside the archive");
  });

  it("shares a parse only for absolute uris, and gives a streamed source's relative uri as its archive path", async () => {
    const seen: (StreamSource | string | null)[] = [];
    const fallback: Loaders = {
      load: none.load,
      async parsed(r) {
        seen.push(r.uri);
        return null;
      },
      stream(source) {
        seen.push(source);
        return null;
      },
    };
    const l = archiveLoaders(files, fallback, "docs");
    await l.parsed?.(req("a.fa"), "fasta", signal);
    await l.parsed?.(req("https://example.org/a.fa"), "fasta", signal);
    const src = {
      ref: "s",
      stream: null,
      protocol: "seqstar-tiles/1",
      uri: "./tiles/aln/",
      data_kind: "alignment",
      parent: null,
      declared: {},
    } as const;
    l.stream?.(src);
    l.stream?.({ ...src, uri: "https://example.org/aln" });
    expect(seen).toEqual([
      "https://example.org/a.fa",
      { ...src, uri: "docs/tiles/aln" },
      { ...src, uri: "https://example.org/aln" },
    ]);
  });

  it("carries the base for validation: compile reads it", () => {
    expect(archiveLoaders(files, none, "docs").base).toBe("docs");
  });
});

describe("nodes added after the document read from an archive", () => {
  it("a snapshot's inserted source may name a path inside the archive, in validate and in snapshotDoc", () => {
    const doc = docAt("https://example.org/a.fa");
    const insert = (uri: string) => ({
      ref: "s1",
      ops: [
        {
          op: "insert" as const,
          parent: "data",
          node: { kind: "source", ref: "more", params: { uri, format: "fasta" } },
        },
      ],
    });
    doc.snapshots = [insert("data/b.fa")];
    expect(validate(doc, reg, { base: "" })).toEqual([]);
    expect(snapshotDoc(doc, "s1", reg, { base: "" }).issues).toEqual([]);
    const out = validate(doc, reg).map((i) => i.code);
    expect(out).toContain("svs.snapshot.invalid");
    expect(snapshotDoc(doc, "s1", reg).issues.map((i) => i.code)).toContain("svs.source.uri");
    doc.snapshots = [insert("../b.fa")];
    expect(validate(doc, reg, { base: "" }).map((i) => i.code)).toContain("svs.snapshot.invalid");
  });

  it("a streamed manifest's dynamic node may name one, when the loaders are the archive's", async () => {
    const node = {
      kind: "annotation" as const,
      key: "sites",
      params: { geometry: "ranges", target: "src/p1/unit", uri: "sites.tsv", format: "tsv" },
    };
    const run = async (loaders: Loaders) => {
      const issues: Issue[] = [];
      const ctx: StageContext = {
        signal: new AbortController().signal,
        registries: reg,
        set: null,
        sources: new Map(),
        unavailable: new Set(),
        source: () => null,
        report: (i) => issues.push(i),
      };
      const part = await loadStreamNode("src", node, (r) => r === "src/p1/unit", ctx, loaders);
      return { part, codes: issues.map((i) => i.code) };
    };
    const files = new Map([["sites.tsv", enc("start\tstop\n0\t2\n")]]);
    const inArchive = await run(archiveLoaders(files, fixtureLoaders()));
    expect(inArchive.codes).toEqual([]);
    expect(inArchive.part?.annotations?.[0]).toMatchObject({
      key: "src/sites",
      geometry: "ranges",
      target: "src/p1/unit",
    });
    const outside = await run(fixtureLoaders());
    expect([outside.part, outside.codes]).toEqual([null, ["svs.source.uri"]]);
  });
});
