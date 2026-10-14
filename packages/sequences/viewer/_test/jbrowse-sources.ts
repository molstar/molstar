/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The JBrowse 2 demo's four sources (as the Seq* Storybook's `Seq*/Cases/Genome browser` makes them, @molstar/seq-core/stream/gmod over
// UCSC's and NCBI's files) over recorded slices of the files (core/_test/stream/gmod/fixtures/<set>/; the .tbi is the main
// fixtures'), replayed; `SEQSTAR_RECORD=1` records a set again, by hand, over curl (Node's fetch is refused in the
// sandbox this is recorded in), never in a test run. Every read is logged, every request too.

import { execFile } from "node:child_process";
import { existsSync, mkdirSync, readFileSync, writeFileSync } from "node:fs";
import type { DataPatch, Loaders, StreamProvider, StreamRequest } from "@molstar/seq-core/format";
import { GRCH38 } from "../../core/_test/format/stream-provider.ts";
import { URLS } from "../../core/_test/stream/gmod/cases.ts";
import { whole } from "../../core/_test/stream/gmod/slices.ts";
import {
  bigBedProvider,
  bigWigProvider,
  chrNames,
  type RangeReader,
  renamed,
  twoBitProvider,
  vcfProvider,
} from "../../core/src/stream/gmod/index.ts";

export const RECORD = process.env.SEQSTAR_RECORD === "1";
export type Key = "twobit" | "bigwig" | "bigbed" | "vcf";
interface Recorded {
  url: string;
  chunks: [number, number, number][];
}

/** Curl reads four at a time (more fail), each in a lane. */
const lanes: Promise<unknown>[] = [Promise.resolve(), Promise.resolve(), Promise.resolve(), Promise.resolve()];
let lane = 0;
const curl = (url: string, offset: number, length: number) => {
  const k = lane++ % lanes.length;
  const p = lanes[k].then(() => curlOne(url, offset, length));
  lanes[k] = p.catch(() => {});
  return p;
};
const curlOne = (url: string, offset: number, length: number) =>
  new Promise<Uint8Array>((ok, no) =>
    execFile(
      "curl",
      ["-s", "-f", "--retry", "3", "-H", `Range: bytes=${offset}-${offset + length - 1}`, url],
      { encoding: "buffer", maxBuffer: 1 << 28 },
      (e, out) => (e ? no(e) : ok(new Uint8Array(out))),
    ),
  );

/** The four providers of the demo over set `name`'s slices: their transports' reads, the requests, and `save`. */
export function demoSources(name: string, only?: readonly string[]) {
  const dir = new URL(`../../core/_test/stream/gmod/fixtures/${name}/`, import.meta.url);
  const saved: Record<string, Recorded> = existsSync(new URL("slices.json", dir))
    ? JSON.parse(readFileSync(new URL("slices.json", dir), "utf8"))
    : {};
  const got = new Map<Key, Map<string, Uint8Array>>();
  const reader = (key: Key) => {
    const reads: [number, number][] = [];
    const bin = RECORD ? null : readFileSync(new URL(`${key}.bin`, dir));
    let at = 0;
    const chunks = RECORD
      ? []
      : saved[key].chunks.map(([offset, asked, n]) => {
          const c = { offset, asked, at, n };
          at += n;
          return c;
        });
    const map = new Map<string, Uint8Array>();
    got.set(key, map);
    const read: RangeReader = async (offset, length) => {
      reads.push([offset, length]);
      if (RECORD) {
        const k = `${offset}:${length}`;
        if (!map.has(k)) map.set(k, await curl(URLS[key], offset, length));
        return map.get(k) as Uint8Array;
      }
      const c = chunks.find((x) => offset >= x.offset && offset + length <= x.offset + x.asked);
      if (!c || !bin) throw new Error(`${key}: bytes ${offset}..${offset + length} are not in the recorded slices`);
      const from = Math.min(offset - c.offset, c.n);
      return new Uint8Array(bin.subarray(c.at + from, c.at + Math.min(offset - c.offset + length, c.n)));
    };
    return { read, reads };
  };
  const io = { twobit: reader("twobit"), bigwig: reader("bigwig"), bigbed: reader("bigbed"), vcf: reader("vcf") };
  const asked: StreamRequest[] = [];
  const logged = (p: StreamProvider): StreamProvider => ({
    fetch: (r, s): Promise<DataPatch> => {
      asked.push(r);
      return p.fetch(r, s);
    },
  });
  const lengths = new Map<string, number>(
    [...GRCH38.map(([k, n]) => [k, n] as [string, number]), ["chrM", 16_569] as [string, number]].filter(
      ([k]) => !only || only.includes(k),
    ),
  );
  const primary = (n: string) => (lengths.has(n) ? n : null);
  const providers: Record<string, StreamProvider> = {
    genome: logged(twoBitProvider(io.twobit.read, { only: (n) => lengths.has(n) })),
    gc: logged(renamed(bigWigProvider(io.bigwig.read), primary)),
    genes: logged(
      renamed(
        bigBedProvider(io.bigbed.read, {
          label: "geneName",
          category: "transcriptClass",
          vocabulary: ["coding", "nonCoding", "pseudo"],
        }),
        primary,
      ),
    ),
    clinvar: logged(
      renamed(
        vcfProvider(io.vcf.read, whole("tbi"), {
          info: ["CLNSIG"],
          lengths: Object.fromEntries([...lengths].map(([k, n]) => [k === "chrM" ? "MT" : k.slice(3), n])),
        }),
        (n) => primary(chrNames(n) ?? ""),
      ),
    ),
  };
  const loaders: Loaders = {
    load: () => Promise.reject(new Error("none")),
    stream: (s) => (s.stream ? (providers[s.stream] ?? null) : null),
  };
  /** In record mode, the reads made written as the set's slices. */
  const save = () => {
    // Nothing read (its test skipped): the set's slices are kept as they are.
    if (!RECORD || [...got.values()].every((m) => m.size === 0)) return;
    mkdirSync(dir, { recursive: true });
    for (const [key, map] of got) {
      const entries = [...map].map(([k, b]) => [...k.split(":").map(Number), b] as [number, number, Uint8Array]);
      entries.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
      writeFileSync(new URL(`${key}.bin`, dir), Buffer.concat(entries.map((e) => e[2])));
      saved[key] = { url: URLS[key], chunks: entries.map(([o, l, b]) => [o, l, b.length]) };
    }
    writeFileSync(new URL("slices.json", dir), `${JSON.stringify(saved, null, 1)}\n`);
  };
  return { io, asked, loaders, save };
}
