/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { hashNode, type Node } from '@molstar/core/repr/tree/index';
import { describe, expect, it } from 'vitest';
import { MsgPackCodec } from '../msgpack.js';

// A large spec of columns: ~25 MB of float32, int32 and uint8 columns beside a few hundred small nodes.
function large(): Node {
  const n = 1_000_000;
  const f = new Float32Array(n).map((_, i) => Math.sin(i));
  const j = new Int32Array(n).map((_, i) => i * 3);
  const u = new Uint8Array(n).map((_, i) => i & 255);
  const children: Node[] = [];
  for (let i = 0; i < 300; i++)
    children.push({ kind: 'row', id: `r${i}`, params: { values: f.slice(i * 1000, i * 1000 + 1000) } });
  return { kind: 'big', params: { f, j, u, f2: f.slice(), j2: j.slice() }, children };
}

describe('io/repr/msgpack throughput', () => {
  it('decodes at least 200 MB/s (performance.md), in place and copying', { tags: ['heavy'], timeout: 60_000 }, () => {
    const node = large();
    const bytes = MsgPackCodec.encode(node);
    const mb = bytes.length / 1e6;
    // the document one byte further in a host: the offsets of its payloads shift, so the aligned ones are copied
    const host = new Uint8Array(bytes.length + 1);
    host.set(bytes, 1);
    const shifted = host.subarray(1);
    const rate = (input: Uint8Array) => {
      let best = 0;
      for (let k = 0; k < 5; k++) {
        const t0 = performance.now();
        const decoded = MsgPackCodec.decode(input);
        const seconds = (performance.now() - t0) / 1000;
        expect(decoded.node).toBeDefined();
        best = Math.max(best, mb / seconds);
      }
      return best;
    };
    const asGiven = rate(bytes);
    const other = rate(shifted);
    console.log(
      `msgpack decode of ${mb.toFixed(1)} MB: ${asGiven.toFixed(0)} MB/s, shifted by one byte ${other.toFixed(0)} MB/s`,
    );
    expect(asGiven).toBeGreaterThan(200);
    expect(other).toBeGreaterThan(200);
    expect(hashNode(MsgPackCodec.decode(bytes).node!)).toBe(hashNode(node));
  });
});
