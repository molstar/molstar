/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SX1's tests of the writer and the bundle step (seqspec:sources#archives). Node only; the archive read back with the reader.

import { readSvsx, type SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { createSeqViewSpec } from "../../src/builder/index.ts";
import { bundleUris, writeSvsx } from "../../src/builder/svsx.ts";

const enc = (s: string) => new TextEncoder().encode(s);
const dec = (b: Uint8Array | undefined) => new TextDecoder().decode(b);

/** Two nodes read the same URI; a third reads another; a streamed source and a hosted file stay as they are. */
function doc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const a = b.data.source({ ref: "a", uri: "https://example.org/x/a.fa", format: "fasta" });
  const a2 = b.data.source({ ref: "a2", uri: "https://example.org/x/a.fa", format: "fasta" });
  const p = b.data.sequence({ ref: "p", source_ref: a });
  b.data.sequence({ ref: "p2", source_ref: a2 });
  b.data.annotation({
    ref: "t",
    geometry: "ranges",
    target: p,
    uri: "https://example.org/x/t.tsv?v=1",
    format: "tsv",
  } as never);
  b.data.alignment({
    ref: "aln",
    row_count: 10,
    column_count: 5,
    source: { uri: "https://example.org/aln", protocol: "seqstar-tiles/1" },
  });
  b.view({ ref: "v", label: "V" }).block({ sequence: p }).track({ label: "L" }).letters();
  return b.getDocument();
}
const FILES: Record<string, string> = {
  "https://example.org/x/a.fa": ">p\nMKV\n",
  "https://example.org/x/t.tsv?v=1": "start\tstop\n0\t2\n",
};
function fakeFetch(status = 200) {
  const requests: string[] = [];
  return {
    requests,
    fetch: async (uri: string) => {
      requests.push(uri);
      return new Response(FILES[uri] ?? "", { status: FILES[uri] === undefined ? 404 : status });
    },
  };
}
const uris = (d: SeqViewSpecDoc) =>
  (d.root.children?.[0].children ?? []).map((n) => [
    n.ref,
    n.params?.uri ?? (n.params?.source as never as { uri: string })?.uri,
  ]);

describe("writeSvsx", () => {
  it("writes index.svsm first by default and index.svsj on request, each read back by the reader", async () => {
    const d = doc();
    const files = { "data/a.fa": enc(">p\nMKV\n") };
    const [m, j] = [
      await readSvsx(await writeSvsx(d, files)),
      await readSvsx(await writeSvsx(d, files, { format: "json" })),
    ];
    expect([m.main, j.main]).toEqual(["index.svsm", "index.svsj"]);
    expect([...m.files.keys()]).toEqual(["index.svsm", "data/a.fa"]);
    expect(dec(m.files.get("data/a.fa"))).toBe(">p\nMKV\n");
    expect(JSON.parse(JSON.stringify(m.doc))).toEqual(JSON.parse(JSON.stringify(d)));
    expect(JSON.parse(JSON.stringify(j.doc))).toEqual(JSON.parse(JSON.stringify(d)));
  });

  it("takes files as a map too, and refuses a path outside the archive or a main file's name", async () => {
    const d = doc();
    expect([...(await readSvsx(await writeSvsx(d, new Map([["a.fa", enc("x")]])))).files.keys()]).toContain("a.fa");
    for (const bad of ["../a.fa", "/a.fa", "a//b", "./a", "dir/", "", "index.svsm", "index.svsj"]) {
      await expect(writeSvsx(d, { [bad]: enc("x") })).rejects.toThrow();
    }
  });
});

describe("bundleUris", () => {
  it("fetches each absolute URI once, in document order, and rewrites it to its stable path", async () => {
    const f = fakeFetch();
    const d = doc();
    const { doc: out, files } = await bundleUris(d, f.fetch);
    expect(f.requests).toEqual(["https://example.org/x/a.fa", "https://example.org/x/t.tsv?v=1"]);
    expect([...files.keys()]).toEqual(["data/example.org/x/a.fa", "data/example.org/x/t.tsv_v_1"]);
    expect(dec(files.get("data/example.org/x/t.tsv_v_1"))).toBe(FILES["https://example.org/x/t.tsv?v=1"]);
    expect(uris(out).slice(0, 2)).toEqual([
      ["a", "data/example.org/x/a.fa"],
      ["a2", "data/example.org/x/a.fa"],
    ]);
    expect(uris(out)[4][1]).toBe("data/example.org/x/t.tsv_v_1");
    // The streamed source is a tile server: kept as it was. The input document is not changed.
    expect(uris(out)[5]).toEqual(["aln.src", "https://example.org/aln"]);
    expect(uris(d)[0]).toEqual(["a", "https://example.org/x/a.fa"]);
  });

  it("also bundles the nodes a snapshot inserts, after the tree's, sharing a file with an equal URI", async () => {
    const f = fakeFetch();
    const d = doc();
    const insert = (ref: string, uri: string) => ({
      op: "insert" as const,
      parent: "data",
      node: { kind: "source", ref, params: { uri, format: "tsv" } },
    });
    d.snapshots = [
      {
        ref: "s1",
        ops: [insert("late", "https://example.org/x/t.tsv?v=1"), insert("again", "https://example.org/x/a.fa")],
      },
      { ref: "empty" },
    ];
    const { doc: out, files } = await bundleUris(d, f.fetch);
    expect(f.requests).toEqual(["https://example.org/x/a.fa", "https://example.org/x/t.tsv?v=1"]);
    expect(files.size).toBe(2);
    const ops = out.snapshots?.[0].ops as unknown as { node: { params: { uri: string } } }[];
    expect(ops.map((o) => o.node.params.uri)).toEqual(["data/example.org/x/t.tsv_v_1", "data/example.org/x/a.fa"]);
    expect(out.snapshots?.[1]).toEqual({ ref: "empty" });
    expect((d.snapshots[0].ops as never as typeof ops)[0].node.params.uri).toBe("https://example.org/x/t.tsv?v=1");
    const r = await readSvsx(await writeSvsx(out, files));
    expect(r.issues.map((i) => i.code)).toEqual([]);
  });

  it("a clash of stable paths gets a number; a failed fetch rejects naming the URI", async () => {
    const b = createSeqViewSpec();
    b.data.source({ ref: "a", uri: "https://example.org/a/b.fa", format: "fasta" });
    b.data.source({ ref: "c", uri: "https://example.org/a_b.fa", format: "fasta" });
    b.data.source({ ref: "d", uri: "https://example.org/a/b.fa?", format: "fasta" });
    const f = {
      "https://example.org/a/b.fa": "1",
      "https://example.org/a_b.fa": "2",
      "https://example.org/a/b.fa?": "3",
    };
    const { files } = await bundleUris(b.getDocument(), async (u) => new Response(f[u as keyof typeof f]));
    expect([...files.keys()]).toEqual([
      "data/example.org/a/b.fa",
      "data/example.org/a_b.fa",
      "data/example.org/a/b.fa-2",
    ]);
    await expect(bundleUris(doc(), fakeFetch(404).fetch)).rejects.toThrow("https://example.org/x/a.fa: HTTP 404");
  });

  it("the bundle written and read back validates with no network and compiles to the same data", async () => {
    const { doc: out, files } = await bundleUris(doc(), fakeFetch().fetch);
    const r = await readSvsx(await writeSvsx(out, files));
    expect(r.issues.map((i) => i.code)).toEqual([]);
    expect(r.files.size).toBe(3);
  });
});
