/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StringLike } from './string-like.js';
import { asciiIndexOf, asciiStartsWith, utf8Encode } from './utf8.js';

/**
 * Data a format is told by (a format provider's `isApplicable`): its first bytes, kept whatever form the rest is held
 * in, and searches over all of it that work on that form, bytes or text, so asking neither decodes nor encodes it.
 */
export interface Sniffable {
  /** The first bytes, up to `Sniffable.HeaderBytes`: magic numbers, first lines, a binary header. */
  readonly header: Uint8Array;
  /** Whether the header holds `needle` (bytes, or an ASCII string) at byte `offset`. */
  sniffBytes(needle: ArrayLike<number> | string, offset?: number): boolean;
  /** Whether the text holds ASCII `needle` anywhere, or with `'start'` at its beginning. */
  sniffText(needle: string, where?: 'anywhere' | 'start'): boolean;
  /**
   * All the bytes when the data is held as bytes (a binary format read whole, such as a BinaryCIF header); `undefined`
   * when it is held as text, which this does not encode.
   */
  readonly heldBytes: Uint8Array | undefined;
}

export namespace Sniffable {
  /** The bytes a header keeps: 64 KiB, cheap to hold per data and enough for every header a format is told by. */
  export const HeaderBytes = 1 << 16;

  export function is(x: unknown): x is Sniffable {
    // not by `header`, which a text's makes on first use
    const o = x as Sniffable | undefined;
    return typeof o?.sniffText === 'function' && typeof o.sniffBytes === 'function';
  }

  /** `header` at `offset` is `needle` */
  export function headerHas(header: Uint8Array, needle: ArrayLike<number> | string, offset = 0): boolean {
    if (offset + needle.length > header.length) return false;
    for (let i = 0; i < needle.length; i++) {
      const c = typeof needle === 'string' ? needle.charCodeAt(i) : needle[i];
      if (header[offset + i] !== c) return false;
    }
    return true;
  }

  /** Searches over UTF-8 bytes, without decoding them */
  export function searchBytes(data: Uint8Array, needle: string, where: 'anywhere' | 'start' = 'anywhere') {
    return where === 'start' ? asciiStartsWith(data, needle) : asciiIndexOf(data, needle) >= 0;
  }

  /** Searches over a text */
  export function searchText(data: StringLike, needle: string, where: 'anywhere' | 'start' = 'anywhere') {
    return where === 'start' ? data.startsWith(needle) : data.includes(needle);
  }

  /** Raw data as a `Sniffable`: bytes are their own header (no copy), a string's header is encoded on first use. */
  export function of(data: Sniffable | StringLike | Uint8Array): Sniffable {
    if (is(data)) return data;
    if (data instanceof Uint8Array) {
      const header = data.length > HeaderBytes ? data.subarray(0, HeaderBytes) : data;
      return {
        header,
        heldBytes: data,
        sniffBytes: (needle, offset) => headerHas(header, needle, offset),
        sniffText: (needle, where) => searchBytes(data, needle, where),
      };
    }
    const text = data as StringLike;
    let header: Uint8Array | undefined;
    const headerOf = () => (header ??= utf8Encode(text.substring(0, HeaderBytes)).subarray(0, HeaderBytes));
    return {
      get header() {
        return headerOf();
      },
      heldBytes: undefined,
      sniffBytes: (needle, offset) => headerHas(headerOf(), needle, offset),
      sniffText: (needle, where) => searchText(text, needle, where),
    };
  }
}
