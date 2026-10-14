/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The transports: HTTP ranges against a local server, and the page cache.

import { createServer, type Server } from "node:http";
import type { AddressInfo } from "node:net";
import { afterAll, beforeAll, describe, expect, it } from "vitest";
import { fetchRanges, gmodProvider, gmodStream, httpFile } from "../../../src/stream/gmod/index.ts";
import { once, paged } from "../../../src/stream/gmod/ranges.ts";
import { req } from "./cases.ts";

const bytes = Uint8Array.from({ length: 1000 }, (_, i) => i % 251);
let server: Server;
let tbiAsked = 0;
let url: string;
beforeAll(async () => {
  server = createServer((req, res) => {
    const m = /bytes=(\d+)-(\d+)/.exec(req.headers.range ?? "");
    if (req.url === "/ignores") return void res.end(Buffer.from(bytes));
    if (req.url?.endsWith(".tbi")) {
      tbiAsked++;
      return void res.writeHead(404).end("not here");
    }
    if (!m || Number(m[1]) >= bytes.length) return void res.writeHead(416).end();
    const [a, b] = [Number(m[1]), Math.min(Number(m[2]), bytes.length - 1)];
    res.writeHead(206).end(Buffer.from(bytes.subarray(a, b + 1)));
  });
  await new Promise<void>((ok) => server.listen(0, "127.0.0.1", ok));
  url = `http://127.0.0.1:${(server.address() as AddressInfo).port}`;
});
afterAll(() => new Promise<void>((ok) => server.close(() => ok())));

describe("fetchRanges", () => {
  it("asks for the range, gets fewer bytes at the end of the file and none past it", async () => {
    const read = fetchRanges(`${url}/f`);
    expect([...(await read(10, 4))]).toEqual([10, 11, 12, 13]);
    expect((await read(990, 100)).length).toBe(10);
    expect((await read(2000, 10)).length).toBe(0);
  });
  it("is an error, not a download, when the server ignores Range", async () => {
    await expect(fetchRanges(`${url}/ignores`)(0, 10)).rejects.toThrow("answered 200");
  });
});

describe("paged", () => {
  it("reads whole pages once and cuts the reads from them", async () => {
    const asked: [number, number][] = [];
    const read = paged(async (o, l) => {
      asked.push([o, l]);
      return bytes.subarray(o, o + l);
    }, 100);
    expect([...(await read(95, 10))]).toEqual([...bytes.subarray(95, 105)]);
    expect([...(await read(150, 5))]).toEqual([...bytes.subarray(150, 155)]);
    expect((await read(990, 100)).length).toBe(10);
    expect(asked).toEqual([
      [0, 100],
      [100, 100],
      [900, 100],
      [1000, 100],
    ]);
  });

  it("keeps a bounded number of pages, the oldest going first", async () => {
    const seen: number[] = [];
    const read = paged(
      async (o, l) => {
        seen.push(o);
        return bytes.subarray(o, o + l);
      },
      100,
      2,
    );
    for (const at of [0, 100, 200, 0]) await read(at, 1);
    expect(seen).toEqual([0, 100, 200, 0]); // page 0 was the oldest of three when two are kept
  });
});

describe("shared reads (CB-01)", () => {
  /** A reader that answers after `release()`, honoring its signal; `signals` the ones it was given. */
  const slow = () => {
    const signals: AbortSignal[] = [];
    let release = () => {};
    const gate = new Promise<void>((ok) => {
      release = ok;
    });
    const read = async (o: number, l: number, signal?: AbortSignal) => {
      if (signal) signals.push(signal);
      await Promise.race([
        gate,
        new Promise((_, no) => signal?.addEventListener("abort", () => no(signal.reason), { once: true })),
      ]);
      return bytes.subarray(o, o + l);
    };
    return { read, signals, release: () => release() };
  };

  it("one caller's abort fails that caller only; the page reads on for the other", async () => {
    const r = slow();
    const read = paged(r.read, 100);
    const [a, b] = [new AbortController(), new AbortController()];
    const first = read(0, 10, a.signal);
    const second = read(5, 10, b.signal);
    a.abort();
    await expect(first).rejects.toThrow(/aborted/);
    r.release();
    expect([...(await second)]).toEqual([...bytes.subarray(5, 15)]);
    expect(r.signals.length).toBe(1);
    expect(r.signals[0]?.aborted).toBe(false);
  });

  it("the last caller to leave aborts the read, and the next asks again", async () => {
    const r = slow();
    const read = paged(r.read, 100);
    const [a, b] = [new AbortController(), new AbortController()];
    const first = read(0, 10, a.signal);
    const second = read(5, 10, b.signal);
    a.abort();
    b.abort();
    await expect(Promise.all([first, second])).rejects.toThrow(/aborted/);
    expect(r.signals[0]?.aborted).toBe(true);
    const third = read(0, 10);
    r.release();
    expect([...(await third)]).toEqual([...bytes.subarray(0, 10)]);
    expect(r.signals.length).toBe(2);
  });

  it("a call whose signal is already aborted starts no read (LA-37: an orphan nothing aborted)", async () => {
    const r = slow();
    const read = paged(r.read, 100);
    const opened = once((s) => r.read(0, 10, s));
    const gone = AbortSignal.abort();
    await expect(read(0, 10, gone)).rejects.toThrow(/aborted/);
    await expect(opened(gone)).rejects.toThrow(/aborted/);
    expect(r.signals.length).toBe(0);
  });
});

describe("gmodStream", () => {
  it("maps a source to a provider, or to none", () => {
    const source = {
      ref: "r",
      stream: null,
      protocol: "gmod",
      uri: "x",
      data_kind: "values",
      parent: null,
      declared: {},
    } as const;
    expect(gmodStream(() => null)(source)).toBeNull();
    expect(gmodStream(() => ({ format: "bigwig", read: async () => new Uint8Array(0) }))(source)?.fetch).toBeTypeOf(
      "function",
    );
  });

  it("a VCF's .tbi is fetched when first needed, its HTTP status is the failure, and the next fetch asks again (CB-06)", async () => {
    const file = httpFile("vcf", `${url}/a.vcf.gz`);
    expect(tbiAsked).toBe(0);
    const p = gmodProvider(file);
    const manifest = req(null);
    const signal = new AbortController().signal;
    await expect(p.fetch(manifest, signal)).rejects.toThrow(/a\.vcf\.gz\.tbi: HTTP 404/);
    await expect(p.fetch(manifest, signal)).rejects.toThrow(/HTTP 404/);
    expect(tbiAsked).toBe(2);
  });

  it("passes a 2bit's `only` to its reader: the manifest reads the sizes it keeps (CB-07)", async () => {
    const asked: number[] = [];
    const file = { format: "2bit" as const, read: twoBitIndex(asked), only: (n: string) => n === "b" };
    const p = gmodProvider(file);
    const manifest = req(null);
    const m = await p.fetch(manifest, new AbortController().signal);
    expect(m.kind === "manifest" && m.manifest.systems.map((s) => s.key)).toEqual(["b"]);
    expect(asked.filter((o) => o >= 100)).toEqual([116]); // b's record only
  });
});

/** A 2bit of sequences "a" and "b" (records at 100 and 116, each a size word of 1), recording read offsets. */
function twoBitIndex(asked: number[]) {
  const b = new Uint8Array(132);
  const dv = new DataView(b.buffer);
  dv.setUint32(0, 0x1a412743, true);
  dv.setUint32(8, 2, true);
  let at = 16;
  for (const [name, offset] of [
    ["a", 100],
    ["b", 116],
  ] as const) {
    b[at] = 1;
    b[at + 1] = name.charCodeAt(0);
    dv.setUint32(at + 2, offset, true);
    at += 6;
  }
  dv.setUint32(100, 1, true);
  dv.setUint32(116, 1, true);
  return async (o: number, n: number) => {
    asked.push(o);
    return b.subarray(o, o + n);
  };
}
