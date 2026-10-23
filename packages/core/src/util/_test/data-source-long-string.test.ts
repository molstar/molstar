/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ajaxGet } from '../data-source.js';

// The string limit and the chunk size made small, so a download past them needs no 512 MB. This file mocks a module,
// so it runs in a registry of its own.
vi.mock('../string-like.js', async (importOriginal) => {
  const actual = await importOriginal<typeof import('../string-like.js')>();
  class SmallChunks extends actual.ChunkedBigString {
    static override fromUtf8Data(data: Uint8Array, start?: number, end?: number) {
      return actual.ChunkedBigString.fromUtf8Data(data, start, end, 4);
    }
  }
  return { ...actual, MAX_STRING_LENGTH: 100, ChunkedBigString: SmallChunks };
});

describe('ajaxGet past the string limit', () => {
  it('decodes a string download into a ChunkedBigString', async () => {
    // multi-byte characters straddle the 16-byte chunk boundaries
    const text = 'data_big\n' + 'Aé€😀'.repeat(40);
    const bytes = new TextEncoder().encode(text);
    const chunks = [bytes.subarray(0, 77), bytes.subarray(77)];
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => {
        let i = 0;
        const body = new ReadableStream<Uint8Array>({
          pull(c) {
            if (i < chunks.length) c.enqueue(chunks[i++]);
            else c.close();
          },
        });
        return new Response(body);
      }),
    );

    const result = await ajaxGet({ url: 'https://example.org/big.cif', type: 'string' }).run();
    const { ChunkedBigString } = await vi.importActual<typeof import('../string-like.js')>('../string-like.js');
    expect(result).toBeInstanceOf(ChunkedBigString);
    expect(result.length).toBe(text.length);
    expect(result.substring(5, 40)).toBe(text.substring(5, 40));
    expect(result.toString()).toBe(text);
  });
});
