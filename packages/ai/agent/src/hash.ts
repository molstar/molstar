/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message } from './messages.js';
import type { Json } from './values.js';

/**
 * Hashes (harness.md §2, §14). `hashMessage` walks a message's fields into bytes, never through JSON text:
 * each value is a type tag plus its content, strings and object keys length-prefixed UTF-8, object
 * keys in sorted order, `undefined` properties skipped. `digest` is RFC 8785 canonical JSON. Both use
 * SHA-256 through WebCrypto, so both are async, in browsers and Node alike.
 */

const TAG = { null: 0, false: 1, true: 2, number: 3, string: 4, array: 5, object: 6 } as const;
const text = new TextEncoder();

class Bytes {
  buffer: Uint8Array<ArrayBuffer> = new Uint8Array(256);
  /** A view of `buffer`, replaced only when `reserve` grows it. */
  private view = new DataView(this.buffer.buffer);
  length = 0;
  private reserve(more: number): void {
    if (this.length + more <= this.buffer.length) return;
    const next = new Uint8Array(Math.max(this.buffer.length * 2, this.length + more));
    next.set(this.buffer.subarray(0, this.length));
    this.buffer = next;
    this.view = new DataView(next.buffer);
  }
  byte(value: number): void {
    this.reserve(1);
    this.buffer[this.length++] = value;
  }
  uint32(value: number): void {
    this.reserve(4);
    this.view.setUint32(this.length, value);
    this.length += 4;
  }
  float64(value: number): void {
    this.reserve(8);
    this.view.setFloat64(this.length, value);
    this.length += 8;
  }
  /** Length-prefixed UTF-8, encoded in place (at most 3 bytes per UTF-16 unit). */
  string(value: string): void {
    this.reserve(4 + value.length * 3);
    const { written } = text.encodeInto(value, this.buffer.subarray(this.length + 4));
    this.view.setUint32(this.length, written);
    this.length += 4 + written;
  }
}

function walk(out: Bytes, value: unknown): void {
  if (value === null) {
    out.byte(TAG.null);
  } else if (typeof value === 'boolean') {
    out.byte(value ? TAG.true : TAG.false);
  } else if (typeof value === 'number') {
    if (!Number.isFinite(value)) throw new TypeError('Cannot hash a non-finite number.');
    out.byte(TAG.number);
    out.float64(value === 0 ? 0 : value);
  } else if (typeof value === 'string') {
    out.byte(TAG.string);
    out.string(value);
  } else if (Array.isArray(value)) {
    out.byte(TAG.array);
    out.uint32(value.length);
    for (const item of value) walk(out, item);
  } else if (typeof value === 'object') {
    const record = value as Record<string, unknown>;
    const keys = Object.keys(record)
      .filter((key) => record[key] !== undefined)
      .sort();
    out.byte(TAG.object);
    out.uint32(keys.length);
    for (const key of keys) {
      out.string(key);
      walk(out, record[key]);
    }
  } else {
    throw new TypeError(`Cannot hash a ${typeof value}.`);
  }
}

/** The bytes `hashMessage` hashes. */
export function encodeMessage(message: Message): Uint8Array<ArrayBuffer> {
  const out = new Bytes();
  walk(out, message);
  return out.buffer.slice(0, out.length);
}

async function sha256(bytes: Uint8Array<ArrayBuffer>): Promise<string> {
  const hash = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes));
  return Array.from(hash, (byte) => byte.toString(16).padStart(2, '0')).join('');
}

/** An entry's id: SHA-256 of the message's structural encoding, lowercase hex. */
export async function hashMessage(message: Message): Promise<string> {
  return await sha256(encodeMessage(message));
}

/** RFC 8785 (JCS): keys sorted by UTF-16 code units, ECMAScript number and string serialization. */
export function canonicalJson(value: Json): string {
  if (value === null || typeof value !== 'object') {
    if (typeof value === 'number' && !Number.isFinite(value)) throw new TypeError('Non-finite number.');
    return JSON.stringify(value);
  }
  if (Array.isArray(value)) return `[${value.map((item) => canonicalJson(item)).join(',')}]`;
  const keys = Object.keys(value).filter((key) => (value as Record<string, Json | undefined>)[key] !== undefined);
  keys.sort();
  return `{${keys.map((key) => `${JSON.stringify(key)}:${canonicalJson((value as Record<string, Json>)[key] as Json)}`).join(',')}}`;
}

/** SHA-256 of `canonicalJson(value)`, lowercase hex. */
export async function digest(value: Json): Promise<string> {
  return await sha256(text.encode(canonicalJson(value)));
}
