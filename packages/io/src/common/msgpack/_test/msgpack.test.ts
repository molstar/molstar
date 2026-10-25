/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createMsgPackState, decodeMsgPack, skipMsgPackValue, type MsgPackExtDecoder } from '../decode.js';
import { encodeMsgPack, type MsgPackExtEncoder } from '../encode.js';

const roundTrip = (v: unknown) => decodeMsgPack(encodeMsgPack(v));

describe('msgpack', () => {
  it('roundtrips strings, containers and the number forms', () => {
    const v = { s: 'Abé☃🧬', n: [0, -1, 127, 128, -33, 65536, -70000, 0.5], m: { t: true, f: false, z: null } };
    expect(roundTrip(v)).toEqual(v);
  });

  it('decodes bin as an exact-size copy, also from a Node Buffer at an offset (perf-audit-3 13/P1)', () => {
    const bytes = Uint8Array.from({ length: 300 }, (_, i) => (i * 7) & 255);
    const encoded = encodeMsgPack({ a: bytes, b: bytes.subarray(10, 13) });
    // the document inside a larger buffer, as a Buffer view and as a Uint8Array view
    const host = new Uint8Array(encoded.length + 9);
    host.set(encoded, 5);
    for (const input of [Buffer.from(host.buffer, 5, encoded.length), host.subarray(5, 5 + encoded.length)]) {
      const v = decodeMsgPack(input) as { a: Uint8Array; b: Uint8Array };
      for (const [out, expected] of [
        [v.a, bytes],
        [v.b, bytes.subarray(10, 13)],
      ] as const) {
        expect(out.constructor).toBe(Uint8Array);
        expect(out.byteOffset).toBe(0);
        expect(out.buffer.byteLength).toBe(expected.length);
        expect(Array.from(out)).toEqual(Array.from(expected));
      }
      // a copy: changing the input leaves the decoded bytes as they were
      input.fill(0);
      expect(v.a[1]).toBe(7);
      host.set(encoded, 5);
    }
  });

  it('writes a character outside the BMP as one 4-byte sequence', () => {
    expect(Array.from(encodeMsgPack('🧬'))).toEqual([0xa4, 0xf0, 0x9f, 0xa7, 0xac]);
  });

  describe('JSON rules', () => {
    it('writes a number that is not finite as nil', () => {
      expect(roundTrip([Number.NaN, Infinity, -Infinity])).toEqual([null, null, null]);
    });

    it('writes an integer outside the 32-bit ranges as a float 64', () => {
      const v = [2 ** 32, 2 ** 40, -(2 ** 31) - 1, -(2 ** 40), Number.MAX_SAFE_INTEGER, 2 ** 32 - 1, -(2 ** 31)];
      expect(roundTrip(v)).toEqual(v);
      expect(encodeMsgPack(2 ** 32)[0]).toBe(0xcb);
      expect(encodeMsgPack(2 ** 32 - 1)[0]).toBe(0xce);
      expect(encodeMsgPack(-(2 ** 31))[0]).toBe(0xd2);
    });

    it('leaves out the keys whose value is undefined, keeps undefined items as nil', () => {
      const r = roundTrip({ a: 1, gone: undefined, b: [1, undefined] }) as object;
      expect(Object.keys(r)).toEqual(['a', 'b']);
      expect(r).toEqual({ a: 1, b: [1, null] });
    });
  });

  describe('extensions', () => {
    const encoder: MsgPackExtEncoder = (v) =>
      v instanceof Float64Array
        ? { type: 0x18, data: new Uint8Array(v.buffer, v.byteOffset, v.byteLength) }
        : undefined;
    const decoder: MsgPackExtDecoder = (type, data) =>
      type === 0x18 ? new Float64Array(data.slice().buffer) : undefined;

    it('asks the extension encoder before the Uint8Array branch, and writes bin without one', () => {
      const bytes = Uint8Array.of(1, 2, 3, 4);
      const claim: MsgPackExtEncoder = (v) => (v instanceof Uint8Array ? { type: 0x12, data: v } : undefined);
      expect(Array.from(encodeMsgPack(bytes, claim))).toEqual([0xd6, 0x12, 1, 2, 3, 4]);
      expect(Array.from(encodeMsgPack(bytes))).toEqual([0xc4, 4, 1, 2, 3, 4]);
      expect(Array.from(encodeMsgPack(bytes, () => undefined))).toEqual([0xc4, 4, 1, 2, 3, 4]);
    });

    it('writes fixext for payloads of 1, 2, 4, 8 and 16 bytes, ext 8, 16 and 32 otherwise', () => {
      const head = (n: number) => encodeMsgPack({ type: 5 }, () => ({ type: 5, data: new Uint8Array(n) }))[0];
      expect([1, 2, 4, 8, 16].map(head)).toEqual([0xd4, 0xd5, 0xd6, 0xd7, 0xd8]);
      expect([0, 3, 17, 255].map(head)).toEqual([0xc7, 0xc7, 0xc7, 0xc7]);
      expect([256, 65535].map(head)).toEqual([0xc8, 0xc8]);
      expect(head(65536)).toBe(0xc9);
    });

    it('roundtrips every payload length through the extension codec', () => {
      for (const n of [0, 1, 2, 3, 4, 5, 31, 32, 8191, 8192, 10_000]) {
        const a = Float64Array.from({ length: n }, (_, i) => i / 10);
        // A string first puts each payload at an odd offset
        const back = decodeMsgPack(encodeMsgPack(['x', a, 'yy'], encoder), decoder) as [string, Float64Array, string];
        expect(back[0]).toBe('x');
        expect(back[2]).toBe('yy');
        expect(Array.from(back[1]), `n = ${n}`).toEqual(Array.from(a));
      }
    });

    it('gives { type, data } with a copy of the payload without a decoder', () => {
      const bytes = encodeMsgPack(Float64Array.of(1), encoder);
      const r = decodeMsgPack(bytes) as { type: number; data: Uint8Array };
      expect(r.type).toBe(0x18);
      expect(r.data.byteLength).toBe(8);
      r.data[0] = 9;
      expect(bytes[2]).toBe(0);
    });

    it('skips extension values', () => {
      for (const n of [1, 5, 16, 300, 70000]) {
        const bytes = encodeMsgPack([new Float64Array(n), 7], encoder);
        const state = createMsgPackState(bytes);
        state.offset = 1; // past the array header
        skipMsgPackValue(state);
        expect(state.offset, `n = ${n}`).toBe(bytes.length - 1);
        expect(bytes[state.offset]).toBe(7);
      }
    });
  });

  describe('untrusted input', () => {
    const header = (type: number) => Uint8Array.of(0x81, 0xa1, 0x61, type, 0xff, 0xff, 0xff, 0xff);

    it('refuses a length header with no payload, for every form', () => {
      // array 32, map 32, bin 32, str 32, ext 32
      for (const type of [0xdd, 0xdf, 0xc6, 0xdb, 0xc9]) {
        expect(() => decodeMsgPack(header(type)), type.toString(16)).toThrow(/truncated/);
      }
    });

    it('refuses a cut document', () => {
      const whole = encodeMsgPack({ a: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], s: 'hello world, hello world, hello' });
      expect(decodeMsgPack(whole)).toBeTruthy();
      for (let at = 1; at < whole.length; at++) {
        expect(() => decodeMsgPack(whole.slice(0, at)), `cut at ${at}`).toThrow(/truncated/);
      }
      expect(() => decodeMsgPack(new Uint8Array(0))).toThrow(/truncated/);
    });

    it('does not size arrays by their header (nested headers each claiming every byte left)', () => {
      const depth = 2000;
      const claim = 1_000_000;
      const b = new Uint8Array(5 * depth + claim);
      const v = new DataView(b.buffer);
      for (let k = 0; k < depth; k++) {
        b[5 * k] = 0xdd;
        v.setUint32(5 * k + 1, claim);
      }
      expect(() => decodeMsgPack(b)).toThrow(/truncated/);
    });

    it('reads a __proto__ key as an own key, as JSON.parse does', () => {
      const v = decodeMsgPack(encodeMsgPack(JSON.parse('{"__proto__": {"hidden": 1}, "a": 2}'))) as Record<
        string,
        unknown
      >;
      expect(Object.keys(v)).toEqual(['__proto__', 'a']);
      expect(Object.getPrototypeOf(v)).toBe(Object.prototype);
      expect('hidden' in v).toBe(false);
    });
  });
});
