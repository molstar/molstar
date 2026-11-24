/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { inflateStream } from '@molstar/core/util/zip/inflate';
import { checkCompressionStreamSupport } from '@molstar/core/util/zip/util';
import { gzipDeclaredSize, isGzip, ungzip } from '@molstar/core/util/zip/zip';
import type { BgzfError } from '@molstar/io/common/bgzf';

// Bytes to text for the genome runtimes (_spec/README.md "Data", Genome runtimes' internals): a file is plain, gzip or
// BGZF, told by its magic, never its extension; bounded by `maxBytes` before and after inflating. Shared by the BED,
// GFF and VCF runtimes, not published (`./runtimes/internal/*`). io's BGZF reader is loaded on a gzip file's magic.

export type Compression = 'none' | 'gzip' | 'bgzf';

export type TextBytes =
  | { readonly ok: true; readonly bytes: Uint8Array; readonly compression: Compression }
  | { readonly ok: false; readonly kind: 'too-large' | 'malformed'; readonly message: string };

/**
 * A file's text bytes: BGZF (`isBgzf`) inflated block by block (`inflateBgzf`, the ISIZEs summed against `maxBytes`
 * first), plain gzip (one member, the platform's decoder; core's `ungzip` where it has none) into a buffer of its
 * declared size (refused past `maxBytes` before it is allocated), else the bytes as they are, within `maxBytes`.
 * Never a throw for a file's content.
 */
export async function textBytes(bytes: Uint8Array, maxBytes: number, signal?: AbortSignal): Promise<TextBytes> {
  if (bytes.length > maxBytes) return tooLarge(`The file is ${bytes.length} bytes, over the bound of ${maxBytes}`);
  // io's BGZF reader loaded only for a gzip file (BGZF is one): a plain file's read never loads it
  const bgzf = isGzip(bytes) ? await import('@molstar/io/common/bgzf') : undefined;
  if (bgzf?.isBgzf(bytes)) {
    try {
      return {
        ok: true,
        bytes: await bgzf.inflateBgzf(bytes, { maxInflatedBytes: maxBytes }, signal),
        compression: 'bgzf',
      };
    } catch (e) {
      if ((e as Error)?.name === 'BgzfError') {
        const b = e as BgzfError;
        return b.kind === 'too-large' ? tooLarge(b.message) : malformed(b.message);
      }
      throw e;
    }
  }
  if (isGzip(bytes)) {
    const declared = gzipDeclaredSize(bytes);
    if (declared > maxBytes) return tooLarge(`The file inflates to ${declared} bytes, over the bound of ${maxBytes}`);
    signal?.throwIfAborted();
    const file = bytes as Uint8Array<ArrayBuffer>;
    const into = new Uint8Array(declared);
    try {
      // the platform's decoder checks the member's CRC-32 and size (a damaged file is refused, not logged); core's
      // `ungzip` where it has none
      const out = checkCompressionStreamSupport('gzip')
        ? await inflateStream(SyncRuntimeContext, file, 'gzip', into)
        : await ungzip(SyncRuntimeContext, file, into);
      if (out.length !== declared) return malformed('The gzip file is shorter than it declares');
      return { ok: true, bytes: out, compression: 'gzip' };
    } catch (e) {
      if ((e as Error)?.name === 'AbortError') throw e;
      // a damaged file, or several plain members (bgzip the file, or gunzip it)
      return malformed(`The gzip file does not inflate as one member: ${(e as Error)?.message ?? e}`);
    }
  }
  return { ok: true, bytes, compression: 'none' };
}

function tooLarge(message: string): TextBytes {
  return { ok: false, kind: 'too-large', message };
}

function malformed(message: string): TextBytes {
  return { ok: false, kind: 'malformed', message };
}

const decoder = new TextDecoder('utf-8');

/**
 * UTF-8 text of `bytes[start, end)` (the records' rule), and whether it holds U+FFFD: a runtime says it once a file
 * (malformed sequences decode to it; a file that writes it is told the same).
 */
export function decodeText(bytes: Uint8Array, start = 0, end = bytes.length): { text: string; replaced: boolean } {
  const text = decoder.decode(bytes.subarray(start, end));
  return { text, replaced: text.includes('�') };
}
