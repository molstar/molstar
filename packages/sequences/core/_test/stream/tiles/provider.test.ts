/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `seqstar-tiles/1` on the stream port (seqspec:sources#the-tiles-protocol, "The provider"): the manifest each data
// kind gets from what its readers declare, the patches that answer its requests, what is not served.

import { StreamNotServed, type StreamSource, streamRequest, TILE } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { PROTOCOL, tileClient, tilesStream } from "../../../src/stream/tiles/index.js";
import { proceduralRow, tileServer } from "../../../src/stream/tiles/server.js";

const signal = new AbortController().signal;
/** Two typed arrays of one kind holding the same bytes (a deep match of a tile's 65,536 cells took over a second). */
const same = (a: unknown, b: ArrayBufferView) => {
  if (!ArrayBuffer.isView(a) || a.constructor !== b.constructor) return false;
  const bytes = (v: ArrayBufferView) => Buffer.from(v.buffer, v.byteOffset, v.byteLength);
  return bytes(a).equals(bytes(b));
};
const source = (over: Partial<StreamSource> = {}): StreamSource => ({
  ref: "s",
  stream: null,
  protocol: PROTOCOL,
  uri: "https://example.org/t",
  data_kind: "alignment",
  parent: null,
  declared: { row_count: 300, column_count: 600 },
  ...over,
});
const server = tileServer({
  rowCount: 300,
  columnCount: 600,
  length: 1000,
  values: (i) => i % 7,
  counts: () => 1,
  stats: { conservation: (c) => (c % 10) / 10 },
});
const provider = (over: Partial<StreamSource> = {}) => {
  const p = tilesStream(() => server.transport)(source(over));
  if (!p) throw new Error("no provider");
  return p;
};
const manifest = async (over: Partial<StreamSource> = {}) => {
  const p = await provider(over).fetch(streamRequest("s", null), signal);
  if (p.kind !== "manifest") throw new Error(p.kind);
  return p.manifest;
};

describe("which sources it serves", () => {
  it("only seqstar-tiles/1 by uri, where the host opens one", () => {
    const stream = tilesStream(() => server.transport);
    expect(stream(source({ protocol: "other/1" }))).toBeNull();
    expect(stream(source({ protocol: null, uri: null, stream: "genome" }))).toBeNull();
    expect(tilesStream(() => null)(source())).toBeNull();
    expect(stream(source())).not.toBeNull();
  });
});

describe("an alignment", () => {
  it("is a manifest of the declared shape: every [2^ly, 2^lx], labels, one system keyed by the source", async () => {
    const m = await manifest();
    expect([m.rows, m.columns, m.labels, m.levels.length, m.levels[16 * 3 + 5]]).toEqual([
      300,
      600,
      true,
      256,
      [8, 32],
    ]);
    expect(m.systems).toEqual([{ key: "s", extent: [0, 600], numbering: null, label: null }]);
    await expect(manifest({ declared: {} })).rejects.toThrow(/declares row_count/);
  });

  it("answers a level with the tile's class and share, labels with labels, a row whole with its label", async () => {
    const p = provider();
    const client = tileClient("x", server.transport);
    const at = { start: 256, stop: 512, row_start: 0, row_stop: 256 };
    const bins = await p.fetch(streamRequest("s", "s", { ...at, level: 0 }), signal);
    const tile = await client.tile(0, 0, 0, 1);
    expect(bins.kind === "bins" && [same(bins.bins.class, tile.cls), same(bins.bins.share, tile.share)]).toEqual([
      true,
      true,
    ]);
    const coarse = await p.fetch(streamRequest("s", "s", { ...at, start: 0, row_start: 0, level: 2 * 16 + 1 }), signal);
    expect(coarse.kind === "bins" && same(coarse.bins.class, (await client.tile(2, 1, 0, 0)).cls)).toBe(true);
    const labels = await p.fetch(
      streamRequest("s", "s", { stop: 600, row_start: 10, row_stop: 12, labels: true }),
      signal,
    );
    expect(labels).toMatchObject({
      kind: "rows",
      rows: null,
      labels: [expect.stringMatching(/^SEQ10\//), expect.stringMatching(/^SEQ11\//)],
    });
    const rows = await p.fetch(streamRequest("s", "s", { stop: 600, row_start: 3, row_stop: 5 }), signal);
    expect(rows).toMatchObject({ kind: "rows", rows: [proceduralRow(3, 600), proceduralRow(4, 600)] });
  });
});

describe("values and items", () => {
  const values = { data_kind: "values" as const, declared: { length: 1000, target: "aln" } };

  it("a values manifest: one system keyed by the target, widths 1 to the one-tile level, the domain of the coarsest tile", async () => {
    const m = await manifest(values);
    expect([m.systems, m.levels, m.domain]).toEqual([
      [{ key: "aln", extent: [0, 1000], numbering: null, label: null }],
      [1, 2, 4],
      [0, 6],
    ]);
    const items = await manifest({ data_kind: "items", declared: { length: 1000 } });
    expect([items.systems[0].key, items.levels, items.vocabulary, items.domain]).toEqual([
      "s",
      [1, 2, 4],
      {},
      undefined,
    ]);
  });

  it("a level is one tile's bins with their mean; the data is level 0's bins over the region, tiles stitched", async () => {
    const p = provider(values);
    const level = await p.fetch(streamRequest("s", "aln", { start: 0, stop: 512, level: 1 }), signal);
    expect(level).toMatchObject({ kind: "bins", bins: { min: expect.any(Float64Array) } });
    const bins = level.kind === "bins" ? level.bins : {};
    expect([bins.min[0], bins.max[0], bins.mean[0], bins.min[249]]).toEqual([0, 1, 0.5, 1]);
    expect(bins.mean.length).toBe(256);
    // Positions 250 to 262 lie in two tiles of level 0.
    const data = await p.fetch(streamRequest("s", "aln", { start: 250, stop: 262 }), signal);
    const d = data.kind === "bins" ? data.bins : {};
    expect(Array.from(d.mean as ArrayLike<number>)).toEqual(Array.from({ length: 12 }, (_, j) => (250 + j) % 7));
    // A region ends with the source.
    const last = await p.fetch(streamRequest("s", "aln", { start: 768, stop: 1000, level: 2 }), signal);
    expect(last.kind === "bins" ? last.bins.mean.length : 0).toBe(58);
  });

  it("items are counts; a statistic is the server's stats/<measure>, its domain the coarsest tile's, else not served", async () => {
    const items = await provider({ data_kind: "items", declared: { length: 1000 } }).fetch(
      streamRequest("s", "s", { start: 0, stop: 256, level: 0 }),
      signal,
    );
    expect(items.kind === "bins" ? Array.from(items.bins.count).slice(0, 3) : []).toEqual([1, 1, 1]);
    const stat = { data_kind: "values" as const, declared: { measure: "conservation", length: 600, target: "aln" } };
    expect((await manifest(stat)).domain).toEqual([0, 0.9]);
    const e = await manifest({ ...stat, declared: { ...stat.declared, measure: "identity" } }).catch((x) => x);
    expect([e instanceof StreamNotServed, e.message]).toEqual([true, "https://example.org/t serves no stats/identity"]);
  });

  it("serves no letters", async () => {
    await expect(manifest({ data_kind: "letters" })).rejects.toBeInstanceOf(StreamNotServed);
    expect(TILE).toBe(256);
  });
});
