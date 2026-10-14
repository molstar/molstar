/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The host's `loaders.stream` for indexed genome files (README, "Mapping a source"): a file's format picks the reader.

import {
  type DataPatch,
  type Manifest,
  StreamNotServed,
  type StreamProvider,
  type StreamSource,
} from "@molstar/seq-core/format";
import { bigBedProvider } from "./bigbed.ts";
import { bigWigProvider } from "./bigwig.ts";
import { fetchRanges, type RangeReader } from "./ranges.ts";
import { type VcfIndex, vcfProvider } from "./tabix.ts";
import { twoBitProvider } from "./twobit.ts";

/**
 * A file a host serves for a source: its format and its reads. A VCF's index file is read whole, by `index`: its
 * bytes, or a function that fetches them, called again after a failure (CB-06). A 2bit's `only` keeps the sequences
 * a source needs, so its manifest reads only their sizes (CB-07).
 */
export type GmodFile =
  | { format: "2bit"; read: RangeReader; only?: (name: string) => boolean }
  | { format: "bigwig" | "bigbed"; read: RangeReader }
  | { format: "vcf"; read: RangeReader; index: VcfIndex; info?: string[] };

/** The reader of a file; BigWig and BigBed share their container, so a mismatched file is an error at its manifest. */
export function gmodProvider(file: GmodFile): StreamProvider {
  if (file.format === "vcf") return vcfProvider(file.read, file.index, { info: file.info });
  if (file.format === "2bit") return twoBitProvider(file.read, { only: file.only });
  return file.format === "bigwig" ? bigWigProvider(file.read) : bigBedProvider(file.read);
}

/** A `loaders.stream` over `open`, which maps a source to its file (null: not served here), over HTTP ranges by default. */
export function gmodStream(open: (source: StreamSource) => GmodFile | null) {
  return (source: StreamSource): StreamProvider | null => {
    const file = open(source);
    return file && gmodProvider(file);
  };
}

/**
 * A file at `url` over HTTP ranges; a VCF's `.tbi` is fetched whole from `${url}.tbi`, when first needed and again
 * after a failure, and an answer other than 2xx is that failure (CB-06: a 404 body was parsed as an index, and one
 * failure failed the source for good).
 */
export function httpFile(format: GmodFile["format"], url: string): GmodFile {
  const read = fetchRanges(url);
  if (format !== "vcf") return { format, read } as GmodFile;
  // Under the open's signal: an open every request left stops its download (LA-36).
  const index = async (signal?: AbortSignal) => {
    const r = await fetch(`${url}.tbi`, { signal });
    if (!r.ok) throw new Error(`${url}.tbi: HTTP ${r.status} ${r.statusText}`);
    return new Uint8Array(await r.arrayBuffer());
  };
  return { format, read, index };
}

/**
 * A provider whose coordinate systems carry the document's names (ST4b): `name` maps a file's name (a VCF's `1`, `MT`)
 * to the document's (`chr1`, `chrM`), null leaving it out; requests are mapped back to the file's, and a name the file
 * does not hold is `StreamNotServed`.
 */
export function renamed(p: StreamProvider, name: (fileName: string) => string | null): StreamProvider {
  const back = new Map<string, string>();
  return {
    async fetch(r, signal) {
      if (r.system === null) {
        const patch = await p.fetch(r, signal);
        if (patch.kind !== "manifest") return patch;
        const systems: Manifest["systems"] = [];
        for (const x of patch.manifest.systems) {
          const key = name(x.key);
          if (key === null) continue;
          back.set(key, x.key);
          systems.push({ ...x, key });
        }
        return { ...patch, manifest: { ...patch.manifest, systems } };
      }
      const file = back.get(r.system);
      if (file === undefined) throw new StreamNotServed(`no sequence "${r.system}" in this file`);
      const patch = await p.fetch({ ...r, system: file }, signal);
      return { ...patch, system: r.system } as DataPatch;
    },
  };
}

/** UCSC's names for a file's (Ensembl and NCBI's): `1` → `chr1`, `MT` → `chrM`; a `chr` name kept; others left out. */
export function chrNames(fileName: string): string | null {
  if (fileName.startsWith("chr")) return fileName;
  if (fileName === "MT") return "chrM";
  return /^([1-9]|1[0-9]|2[0-2]|X|Y)$/.test(fileName) ? `chr${fileName}` : null;
}
