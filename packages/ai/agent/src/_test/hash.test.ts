/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { canonicalJson, digest, hashMessage, type Json, type Message } from '../index.js';

const call: Message = {
  role: 'assistant',
  parts: [
    { kind: 'text', text: 'Loading 1BL8.' },
    { kind: 'tool-call', requestId: 'r1', tool: { id: 'load_structure' }, arguments: { pdb_id: '1BL8', assembly: 1 } },
  ],
};

describe('hashMessage', () => {
  it('is a 64-character lowercase hex SHA-256', async () => {
    expect(await hashMessage(call)).toMatch(/^[0-9a-f]{64}$/);
  });

  it('does not depend on key order', async () => {
    const reordered: Message = {
      parts: [
        { text: 'Loading 1BL8.', kind: 'text' },
        {
          arguments: { assembly: 1, pdb_id: '1BL8' },
          tool: { id: 'load_structure' },
          requestId: 'r1',
          kind: 'tool-call',
        },
      ],
      role: 'assistant',
    };
    expect(await hashMessage(reordered)).toBe(await hashMessage(call));
  });

  it('is unchanged by an absent label, and changed by a label (it is part of the entry)', async () => {
    const hash = await hashMessage(call);
    expect(await hashMessage({ ...call, label: undefined })).toBe(hash);
    const labelled: Message = { ...call, label: { kind: 'host:step' } };
    expect(await hashMessage(labelled)).not.toBe(hash);
    expect(await hashMessage({ ...call, label: { kind: 'host:step', key: 'k' } })).not.toBe(
      await hashMessage(labelled),
    );
  });

  it('treats an undefined property as absent', async () => {
    const withUndefined = {
      role: 'user',
      parts: [{ kind: 'text', text: 'hi', extra: undefined }],
    } as unknown as Message;
    expect(await hashMessage(withUndefined)).toBe(
      await hashMessage({ role: 'user', parts: [{ kind: 'text', text: 'hi' }] }),
    );
  });

  it('changes with any part', async () => {
    const base = await hashMessage(call);
    const variants: Message[] = [
      { ...call, role: 'assistant', parts: call.parts.slice(0, 1) },
      { role: 'assistant', parts: [{ kind: 'text', text: 'Loading 1BL9.' }, ...call.parts.slice(1)] },
      {
        role: 'assistant',
        parts: [
          call.parts[0] as never,
          {
            kind: 'tool-call',
            requestId: 'r1',
            tool: { id: 'load_structure' },
            arguments: { pdb_id: '1BL8', assembly: 2 },
          },
        ],
      },
      { role: 'user', parts: [{ kind: 'text', text: 'Loading 1BL8.' }] },
    ];
    for (const variant of variants) expect(await hashMessage(variant)).not.toBe(base);
  });

  it('keeps types apart: a string is not a number, an array is not an object', async () => {
    const as = (value: Json): Message => ({
      role: 'assistant',
      parts: [{ kind: 'tool-call', requestId: 'r', tool: { id: 't' }, arguments: value }],
    });
    const hashes = await Promise.all([as('1'), as(1), as([1]), as({ 0: 1 }), as(null), as(false)].map(hashMessage));
    expect(new Set(hashes).size).toBe(6);
  });

  it('orders arrays: the same parts in another order is another message', async () => {
    const swapped: Message = { role: 'assistant', parts: [call.parts[1] as never, call.parts[0] as never] };
    expect(await hashMessage(swapped)).not.toBe(await hashMessage(call));
  });

  it('treats -0 as 0 and refuses non-finite numbers', async () => {
    const as = (n: number): Message => ({
      role: 'assistant',
      parts: [{ kind: 'tool-call', requestId: 'r', tool: { id: 't' }, arguments: n }],
    });
    expect(await hashMessage(as(-0))).toBe(await hashMessage(as(0)));
    await expect(hashMessage(as(Number.NaN))).rejects.toThrow(TypeError);
  });
});

describe('canonicalJson and digest (RFC 8785)', () => {
  it('serializes numbers, strings and literals as RFC 8785 section 3.2.2 shows', () => {
    const input = JSON.parse(
      '{"numbers":[333333333.33333329,1E30,4.50,2e-3,0.000000000000000000000000001],"string":"\\u20ac$\\u000F\\u000aA\'\\u0042\\u0022\\u005c\\\\\\"\\/","literals":[null,true,false]}',
    ) as Json;
    expect(canonicalJson(input)).toBe(
      '{"literals":[null,true,false],"numbers":[333333333.3333333,1e+30,4.5,0.002,1e-27],"string":"€$\\u000f\\nA\'B\\"\\\\\\\\\\"/"}',
    );
  });

  it('sorts keys by UTF-16 code units, as RFC 8785 section 3.2.3 shows', () => {
    const input = JSON.parse(
      '{"\\u20ac":"Euro Sign","\\r":"Carriage Return","\\ufb33":"Hebrew Letter Dalet With Dagesh","1":"One","\\ud83d\\ude00":"Emoji: Grinning Face","\\u0080":"Control","\\u00f6":"Latin Small Letter O With Diaeresis"}',
    ) as Record<string, string>;
    const canonical = canonicalJson(input);
    const names = [
      'Carriage Return',
      'One',
      'Control',
      'Latin Small Letter O With Diaeresis',
      'Euro Sign',
      'Emoji: Grinning Face',
      'Hebrew Letter Dalet With Dagesh',
    ];
    const positions = names.map((name) => canonical.indexOf(`"${name}"`));
    expect(positions.every((p) => p >= 0)).toBe(true);
    expect([...positions].sort((x, y) => x - y)).toEqual(positions);
  });

  it('keeps integer-like keys in code-unit order, not property order (defect 32)', () => {
    expect(canonicalJson({ 10: 'a', 2: 'b', b: 'c', a: 'd' })).toBe('{"10":"a","2":"b","a":"d","b":"c"}');
  });

  it('digests to SHA-256 hex, equal for equal values', async () => {
    expect(await digest({ b: 1, a: [true, null] })).toBe(await digest({ a: [true, null], b: 1 }));
    expect(await digest({ a: 1 })).toMatch(/^[0-9a-f]{64}$/);
  });
});
