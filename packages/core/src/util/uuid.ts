/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { now } from './now.js';

/** A UUID, either standard 36 characters or 22 characters base64 encoded. */
type UUID = string & { '@type': 'uuid' };

namespace UUID {
  const Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_';
  /** Random characters of the prefix: 14 × 6 = 84 bits, renewed when the counter wraps. */
  const PrefixLength = 14;
  /** Counter characters after the prefix: 7 slow ones, recomputed every 64 ids, and the fast last one. */
  const CounterHighLength = 7;
  const CounterHighMax = 64 ** CounterHighLength;

  let prefix = '';
  let counterHigh = 0;
  let counterLow = 0;
  let head = '';

  function randomChars(n: number) {
    const bytes = new Uint8Array(n);
    if (typeof crypto !== 'undefined' && typeof crypto.getRandomValues === 'function') {
      crypto.getRandomValues(bytes);
    } else {
      for (let i = 0; i < n; i++) bytes[i] = (Math.random() * 256) | 0;
    }
    let s = '';
    for (let i = 0; i < n; i++) s += Alphabet[bytes[i] & 63];
    return s;
  }

  function nextHead() {
    if (counterHigh >= CounterHighMax || !prefix) {
      prefix = randomChars(PrefixLength);
      counterHigh = 0;
    }
    let s = '';
    for (let i = 0, c = counterHigh; i < CounterHighLength; i++, c = Math.floor(c / 64)) s = Alphabet[c % 64] + s;
    head = prefix + s;
    counterHigh++;
  }

  /**
   * Creates a 22 characters id in the URL-safe base64 alphabet: a random prefix (per process, renewed every 2^48 ids)
   * and a counter. Unique within the process by the counter and across processes by the prefix; not an RFC 4122 UUID,
   * not random per call, not a secret. Creating one is a string concatenation (ids are made for every state transform
   * and its version, many per frame of an animation).
   */
  export function create22(): UUID {
    if (counterLow === 0) nextHead();
    const id = head + Alphabet[counterLow];
    counterLow = (counterLow + 1) & 63;
    return id as UUID;
  }

  export function createv4(): UUID {
    let d = +new Date() + now();
    const uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function (c) {
      const r = ((d + Math.random() * 16) % 16) | 0;
      d = Math.floor(d / 16);
      return (c === 'x' ? r : (r & 0x3) | 0x8).toString(16);
    });
    return uuid as any;
  }

  export function is(x: any): x is UUID {
    return typeof x === 'string';
  }
}

export { UUID };
