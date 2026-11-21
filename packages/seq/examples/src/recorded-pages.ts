/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssetRangeReader } from '@molstar/core/util/assets';
import type { AssetReader, FormatRuntime } from '@molstar/seqx-core/data';

// Real files recorded as pages (examples.md, "Real slices"): the 64 KiB pages a runtime read of a published file (core's
// page size), kept byte for byte in one JSON file beside where the file would be (`<uri>.pages.json`), so a test reads
// the container the real tools wrote. A read of a page the recording did not make fails, naming the page and the file:
// a change to a runtime's reads is a re-recording (`scripts/record-pages.ts`). Browser-safe: no `node:` import.
//
//   recordingRanges(bytes)        a range reader over a downloaded file that notes the pages read (to record them)
//   encodePages(file, bytes, pages)   the recording of those pages (or 'all')
//   recordedAssets(inner)         an asset port reading `<uri>` from `<uri>.pages.json` through `inner`
//   withRecordedPages(runtime)    a runtime whose sources read their files so (a view's registry entry)

/** Core's page size (`openRanges`' default): the unit a recording keeps. */
export const PageSize = 64 * 1024;

/** The file suffix of a recording: `fixtures/genomes/cytoBandIdeo.txt.gz.pages.json` records `…/cytoBandIdeo.txt.gz`. */
export const PagesSuffix = '.pages.json';

/** A recording: the file's name and size, and its pages read, each in base64 by its number. */
export interface RecordedPages {
  readonly format: 'seqx-pages/1';
  readonly file: string;
  readonly size: number;
  readonly pageSize: number;
  readonly pages: Readonly<Record<string, string>>;
}

/** A read the recording did not make: the page and the file named. Matched by `name`. */
export class PageNotRecorded extends Error {
  readonly file: string;
  readonly page: number;
  // no parameter properties: the recording script imports this module with Node's type stripping
  constructor(file: string, page: number) {
    super(`Page ${page} of ${file} was not recorded: record the file again (scripts/record-pages.ts)`);
    this.name = 'PageNotRecorded';
    this.file = file;
    this.page = page;
  }
}

/** The pages `[offset, offset + length)` ranges touch, within a file of `size` bytes. */
export function pagesOf(
  ranges: readonly { readonly offset: number; readonly length: number }[],
  size: number,
  pageSize = PageSize,
): number[] {
  const out = new Set<number>();
  for (const r of ranges) {
    const end = Math.min(size, r.offset + Math.max(0, r.length));
    for (let p = Math.floor(r.offset / pageSize); p * pageSize < end; p++) out.add(p);
  }
  return [...out].sort((a, b) => a - b);
}

function toBase64(bytes: Uint8Array): string {
  let s = '';
  for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(s);
}

function fromBase64(text: string): Uint8Array {
  const s = atob(text);
  const out = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
  return out;
}

/** The recording of `pages` of a file's bytes (`'all'`: every page, the file whole). */
export function encodePages(file: string, bytes: Uint8Array, pages: Iterable<number> | 'all'): RecordedPages {
  const count = Math.ceil(bytes.length / PageSize);
  const list = pages === 'all' ? Array.from({ length: count }, (_, i) => i) : [...new Set(pages)].sort((a, b) => a - b);
  const out: Record<string, string> = {};
  for (const p of list) {
    if (!(p >= 0 && p < count)) throw new Error(`No page ${p} in ${file} (${count} pages)`);
    out[p] = toBase64(bytes.subarray(p * PageSize, (p + 1) * PageSize));
  }
  return { format: 'seqx-pages/1', file, size: bytes.length, pageSize: PageSize, pages: out };
}

/** A recording checked: its format, page size and pages' lengths. */
export function readRecording(text: string): RecordedPages {
  const r = JSON.parse(text) as RecordedPages;
  if (
    r?.format !== 'seqx-pages/1' ||
    r.pageSize !== PageSize ||
    typeof r.size !== 'number' ||
    typeof r.pages !== 'object'
  )
    throw new Error('Not a recording of pages (seqx-pages/1)');
  return r;
}

/** A recording's bytes over `[offset, offset + length)`, every page they touch recorded or `PageNotRecorded`. */
function recordedRead(r: RecordedPages, decoded: Map<number, Uint8Array>, offset: number, length: number): Uint8Array {
  const end = Math.min(r.size, offset + Math.max(0, length));
  if (end <= offset) return new Uint8Array(0);
  const out = new Uint8Array(end - offset);
  for (const p of pagesOf([{ offset, length: end - offset }], r.size)) {
    let page = decoded.get(p);
    if (!page) {
      const b64 = r.pages[p];
      if (b64 === undefined) throw new PageNotRecorded(r.file, p);
      page = fromBase64(b64);
      decoded.set(p, page);
    }
    const from = Math.max(offset, p * PageSize);
    const to = Math.min(end, (p + 1) * PageSize);
    out.set(page.subarray(from - p * PageSize, to - p * PageSize), from - offset);
  }
  return out;
}

/** A range reader over a recording. */
export function recordedRanges(r: RecordedPages): AssetRangeReader {
  const decoded = new Map<number, Uint8Array>();
  return {
    size: r.size,
    read: async (ranges, signal) => {
      signal?.throwIfAborted();
      return ranges.map((g) => recordedRead(r, decoded, g.offset, g.length));
    },
    dispose() {},
  };
}

/** A range reader over a downloaded file that notes every page read (`pages`): what a recording keeps. */
export function recordingRanges(bytes: Uint8Array): AssetRangeReader & { readonly pages: ReadonlySet<number> } {
  const pages = new Set<number>();
  return {
    size: bytes.length,
    pages,
    read: async (ranges) => {
      for (const p of pagesOf(ranges, bytes.length)) pages.add(p);
      return ranges.map((g) => bytes.subarray(g.offset, Math.min(bytes.length, g.offset + g.length)));
    },
    dispose() {},
  };
}

/**
 * An asset port reading each file from its recording (`<uri>.pages.json`, through `inner`): whole (every page
 * recorded) by `readBytes`, by range by `openRange`. A location without a recording is `inner`'s own.
 */
export function recordedAssets(inner: AssetReader): AssetReader {
  const recordings = new Map<string, Promise<RecordedPages | undefined>>();
  const recording = (uri: string, signal: AbortSignal) => {
    let p = recordings.get(uri);
    if (!p) {
      p = inner.readText(`${uri}${PagesSuffix}`, signal).then(readRecording, () => undefined);
      recordings.set(uri, p);
    }
    return p;
  };
  return {
    readText: async (uri, signal) => {
      const r = await recording(uri, signal);
      return r ? new TextDecoder().decode(recordedRead(r, new Map(), 0, r.size)) : inner.readText(uri, signal);
    },
    readBytes: async (uri, signal) => {
      const r = await recording(uri, signal);
      return r ? recordedRead(r, new Map(), 0, r.size) : inner.readBytes(uri, signal);
    },
    openRange: async (uri, signal) => {
      const r = await recording(uri, signal);
      if (r) return recordedRanges(r);
      if (!inner.openRange) throw new Error(`No range reads of ${uri}`);
      return inner.openRange(uri, signal);
    },
  };
}

/** A runtime whose sources read their files from recordings (`recordedAssets` over the context's port). */
export function withRecordedPages(runtime: FormatRuntime): FormatRuntime {
  return {
    name: `${runtime.name}-recorded`,
    formats: runtime.formats,
    open: (source, ctx) => runtime.open(source, { ...ctx, assets: recordedAssets(ctx.assets) }),
  };
}
