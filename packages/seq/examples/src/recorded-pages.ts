/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssetRangeReader } from '@molstar/core/util/assets';
import type { AssetReader, FormatRuntime } from '@molstar/seqx-core/data';

// Real files recorded as pages (_spec/README.md "Fixtures"): the 64 KiB pages a runtime read of a published file
// (core's page size), kept byte for byte in one JSON file beside where the file would be (`<uri>.pages.json`), so a
// test reads the container the real tools wrote. A read of a page the recording did not make fails, naming the page and
// the file: a change to a runtime's reads is a re-recording (`scripts/record-pages.ts`). Browser-safe: no `node:`
// import.
//
//   recordingRanges(bytes)        a range reader over a downloaded file that notes the pages read (to record them)
//   encodePages(file, bytes, pages)   the recording of those pages (or 'all')
//   recordedAssets(inner, uris)   an asset port reading those `uris` from `<uri>.pages.json` through `inner`
//   withRecordedPages(runtime, uris)  a runtime whose sources read those files so (a view's registry entry)

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

/** The bytes a base64 text decodes to, told from its length and padding. */
const base64Bytes = (t: string) => (t.length / 4) * 3 - (t.endsWith('==') ? 2 : t.endsWith('=') ? 1 : 0);

/**
 * A recording checked: its format and page size, and each page's length (a page a full 64 KiB, the file's last its
 * rest), so a page cut short is refused here, never read as zeros; the message names the page and the file.
 */
export function readRecording(text: string): RecordedPages {
  const r = JSON.parse(text) as RecordedPages;
  if (
    r?.format !== 'seqx-pages/1' ||
    r.pageSize !== PageSize ||
    !Number.isInteger(r.size) ||
    r.size < 0 ||
    typeof r.pages !== 'object' ||
    r.pages === null
  )
    throw new Error('Not a recording of pages (seqx-pages/1)');
  const count = Math.ceil(r.size / PageSize);
  for (const [k, b64] of Object.entries(r.pages)) {
    const p = Number(k);
    if (!Number.isInteger(p) || p < 0 || p >= count)
      throw new Error(`Page ${k} of ${r.file} is past its ${count} pages`);
    const want = Math.min(PageSize, r.size - p * PageSize);
    if (typeof b64 !== 'string' || b64.length % 4 !== 0 || base64Bytes(b64) !== want)
      throw new Error(`Page ${p} of ${r.file} is not its ${want} bytes: record the file again`);
  }
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
 * An asset port reading the files `recorded` lists from their recordings (`<uri>.pages.json`, through `inner`): whole
 * (every page recorded) by `readBytes`, by range by `openRange`; any other location is `inner`'s own, never asked for a
 * recording. A recording that cannot be read (missing, malformed, a page cut short) or an abort fails the read, naming
 * the cause: never a fallback to the file.
 */
export function recordedAssets(inner: AssetReader, recorded: readonly string[]): AssetReader {
  const recordings = new Map<string, Promise<RecordedPages>>();
  const recording = (uri: string, signal: AbortSignal): Promise<RecordedPages> | undefined => {
    if (!recorded.includes(uri)) return undefined;
    let p = recordings.get(uri);
    if (!p) {
      p = inner.readText(`${uri}${PagesSuffix}`, signal).then(readRecording);
      // a failure is not kept: an aborted read is asked again by the next
      p.catch(() => recordings.delete(uri));
      recordings.set(uri, p);
    }
    return p;
  };
  return {
    readText: async (uri, signal) => {
      const p = recording(uri, signal);
      if (!p) return inner.readText(uri, signal);
      const r = await p;
      return new TextDecoder().decode(recordedRead(r, new Map(), 0, r.size));
    },
    readBytes: async (uri, signal) => {
      const p = recording(uri, signal);
      if (!p) return inner.readBytes(uri, signal);
      const r = await p;
      return recordedRead(r, new Map(), 0, r.size);
    },
    openRange: async (uri, signal) => {
      const p = recording(uri, signal);
      if (p) return recordedRanges(await p);
      if (!inner.openRange) throw new Error(`No range reads of ${uri}`);
      return inner.openRange(uri, signal);
    },
  };
}

/** A runtime whose sources read the files `recorded` lists from their recordings (`recordedAssets`). */
export function withRecordedPages(runtime: FormatRuntime, recorded: readonly string[]): FormatRuntime {
  return {
    name: `${runtime.name}-recorded`,
    formats: runtime.formats,
    open: (source, ctx) => runtime.open(source, { ...ctx, assets: recordedAssets(ctx.assets, recorded) }),
  };
}
