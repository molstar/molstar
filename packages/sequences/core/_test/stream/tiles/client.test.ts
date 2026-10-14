/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `seqstar-tiles/1` client over the reference server (moved from the format's load and protocol tests, ST2): the
// tiles' bytes, rows and labels, the 1-D tiles, and the transport's failures.

import { StreamNotServed, TILE } from "@molstar/seq-core/format";
import { describe, expect, it, vi } from "vitest";
import { fetchTransport, NotServed, tileClient } from "../../../src/stream/tiles/index.ts";
import { proceduralRow, tileServer } from "../../../src/stream/tiles/server.ts";

describe("the client on the reference server", () => {
  const [R, C] = [300, 600];
  const server = tileServer({ rowCount: R, columnCount: C, length: 1000, values: (i) => i % 7, counts: () => 1 });
  const client = tileClient("https://example.org/t", server.transport);
  const row = (k: number) => proceduralRow(k, C);

  it("serves level [0, 0] as the cells, and aggregates above it", async () => {
    const t0 = await client.tile(0, 0, 1, 2);
    expect(String.fromCharCode(...t0.cls.subarray(5 * TILE, 5 * TILE + 10))).toBe(row(256 + 5).slice(512, 522));
    expect(t0.share[5 * TILE]).toBe(255);
    expect(t0.cls[50 * TILE]).toBe(0); // row 306: past the rows, no data
    const t1 = await client.tile(1, 2, 0, 0);
    const cells = [0, 1].flatMap((r) => [...row(2 * 7 + r).slice(4 * 3, 4 * 3 + 4)]);
    const n = (x: string) => cells.filter((c) => c === x).length;
    const best = [...new Set(cells)]
      .filter((c) => c !== "-")
      .sort()
      .reduce((b, c) => (n(c) > n(b) ? c : b));
    expect(String.fromCharCode(t1.cls[7 * TILE + 3])).toBe(best);
    expect(t1.share[7 * TILE + 3]).toBe(Math.round((n(best) / 8) * 255));
  });

  it("a coarse sample takes its most frequent residue at its share of all cells; only gaps stay a gap", async () => {
    // Row 0 mostly gaps, row 1 gaps only: samples of 2 rows × 4 columns.
    const rows = ["--A---CA", "--------"];
    const gappy = tileClient(
      "https://example.org/g",
      tileServer({ rowCount: 2, columnCount: 8, row: (k) => rows[k] }).transport,
    );
    const t = await gappy.tile(1, 2, 0, 0);
    expect([String.fromCharCode(t.cls[0]), t.share[0]]).toEqual(["A", Math.round((1 / 8) * 255)]);
    expect([String.fromCharCode(t.cls[1]), t.share[1]]).toEqual(["A", Math.round((1 / 8) * 255)]);
    const g = await gappy.tile(1, 1, 0, 0);
    expect([String.fromCharCode(g.cls[0]), g.share[0]]).toEqual(["-", 255]);
  });

  it("fetches rows whole and labels; 1-D tiles; statistics not served are null; a short answer is an error", async () => {
    const whole = await client.row(3);
    expect([whole.row, whole.label]).toEqual([row(3), expect.stringMatching(/^SEQ3\/1-/)]);
    expect((await client.labels(10, 12)).map((l) => l.split("/")[0])).toEqual(["SEQ10", "SEQ11"]);
    const v = await client.values(2, 0);
    expect([v[0], v[1], v[2 * 249], v[2 * 249 + 1], v[2 * 250]]).toEqual([0, 3, 2, 5, Number.NaN]);
    expect((await client.counts(0, 3)).slice(230, 234)).toEqual(Uint32Array.from([1, 1, 0, 0]));
    expect(await client.stats("conservation", 0, 0)).toBeNull();
    const calls = server.calls.length;
    await expect(tileClient("x", async () => new Uint8Array(3)).tile(0, 0, 0, 0)).rejects.toThrow(/3 bytes/);
    expect(server.calls.length).toBe(calls);
  });
});

describe("the transport", () => {
  const outage = new Error("503 Service Unavailable");
  const calls = (client: ReturnType<typeof tileClient>) => ({
    tile: () => client.tile(0, 0, 0, 0),
    values: () => client.values(0, 0),
    counts: () => client.counts(0, 0),
    labels: () => client.labels(0, 2),
    row: () => client.row(0),
    stats: () => client.stats("conservation", 0, 0),
  });

  it("a JSON answer of another shape is an error naming the path, not data (GI-13)", async () => {
    const answering = (body: unknown) => tileClient("x", async () => new TextEncoder().encode(JSON.stringify(body)));
    await expect(answering({ length: 1_000_000_000 }).labels(0, 100)).rejects.toThrow(
      "x/labels/0/100: the answer is not a list of 100 labels",
    );
    await expect(answering(Array.from({ length: 2001 }, () => "a")).labels(0, 100)).rejects.toThrow(/100 labels/);
    await expect(answering(["a", 7]).labels(0, 2)).rejects.toThrow(/2 labels/);
    expect(await answering(["a", "b"]).labels(0, 2)).toEqual(["a", "b"]);
    await expect(answering({ label: "a" }).row(3)).rejects.toThrow("x/rows/3: the answer is not a row");
    expect(await answering({ label: "a", row: "AC-" }).row(3)).toEqual({ label: "a", row: "AC-" });
  });

  it("not served: every call rejects with NotServed (a StreamNotServed) but stats, which is null", async () => {
    const client = tileClient("x", async (path) => {
      throw new NotServed(`404 ${path}`);
    });
    const { stats, ...rest } = calls(client);
    for (const [name, call] of Object.entries(rest)) await expect(call(), name).rejects.toBeInstanceOf(StreamNotServed);
    expect(await stats()).toBeNull();
  });

  it("an outage and an abort are failures for every call, stats included", async () => {
    const control = new AbortController();
    for (const failure of [outage, new DOMException("aborted", "AbortError")]) {
      const client = tileClient("x", async () => {
        throw failure;
      });
      for (const [name, call] of Object.entries(calls(client))) {
        const e = await call().then(
          () => null,
          (x) => x,
        );
        expect([name, e, e instanceof NotServed]).toEqual([name, failure, false]);
      }
    }
    control.abort();
    const server = tileServer({ rowCount: 4, columnCount: 4, stats: { conservation: () => 1 } });
    const client = tileClient("x", (p, signal) => {
      signal?.throwIfAborted();
      return server.transport(p, signal);
    });
    await expect(client.stats("conservation", 0, 0, control.signal)).rejects.toThrow();
  });

  it("fetchTransport: 404 and 410 are not served, a 503 is a failure", async () => {
    const answer = (status: number) =>
      vi.stubGlobal("fetch", async () => new Response("", { status, statusText: "answer" }));
    try {
      for (const status of [404, 410]) {
        answer(status);
        await expect(fetchTransport("https://example.org/t")("rows/0")).rejects.toBeInstanceOf(NotServed);
      }
      answer(503);
      const e = await fetchTransport("https://example.org/t")("rows/0").catch((x) => x);
      expect([e instanceof Error, e instanceof NotServed]).toEqual([true, false]);
    } finally {
      vi.unstubAllGlobals();
    }
  });

  it("binary payloads are little-endian", async () => {
    const bytes = new Uint8Array(16 * 256);
    new DataView(bytes.buffer).setFloat64(0, 1.5, true);
    new DataView(bytes.buffer).setFloat64(8, -2, true);
    const counts = new Uint8Array(4 * 256);
    new DataView(counts.buffer).setUint32(4, 0x01020304, true);
    const client = tileClient("x", async (p) => (p.startsWith("counts") ? counts : bytes));
    expect(Array.from((await client.values(0, 0)).subarray(0, 2))).toEqual([1.5, -2]);
    expect((await client.counts(0, 0))[1]).toBe(0x01020304);
  });
});
