/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, hashNode, type Node } from '../../../tree/index.js';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { MsgPackCodec } from '../msgpack.js';

// Decode throughput (at least 200 MB/s, over the bytes of the document: every byte a copy or a parse touches), on the
// two paths that read bytes. A column decoded in place reads none, so a document of such columns says nothing: the
// node-heavy document is the parse, the unaligned column the copy. The quick lane runs the same bodies with no bound.

const MinMBps = 200;
const Repeats = 5;
/**
 * A node-heavy document is bound by the parse and the checks of its many small values (strings, keys, extensions of a
 * few bytes), not by bytes: about 45 MB/s measured, 130 MB/s for the bare MessagePack parse of it. The bound guards a
 * regression there; the 200 MB/s of the performance targets is that of documents of columns (below).
 */
const NodeHeavyMBps = 20;

/** Milliseconds `repeats` decodes of `bytes` may take at `MinMBps`. */
const budget = (bytes: number, mbps = MinMBps, repeats = Repeats) => (bytes / (mbps * 1e6)) * 1000 * repeats;

/** Many small nodes: strings, numbers, a short column and a category each. */
function nodeHeavy(count: number): Node {
  const children: Node[] = [];
  for (let i = 0; i < count; i++)
    children.push({
      kind: 'row',
      id: `r${i}`,
      params: {
        name: `residue ${i} ä`,
        index: i,
        score: i / 7,
        flags: [true, false, null, i],
        start: new Int32Array([i, i + 1, i + 2]),
        label: categoryColumn(['a', 'b', null, `c${i % 5}`]),
      },
    });
  return { kind: 'doc', children };
}

describe('repr/codecs/msgpack throughput', () => {
  itWithin('parses a node-heavy document at NodeHeavyMBps', (within) => {
    const node = nodeHeavy(20_000);
    const bytes = MsgPackCodec.encode(node);
    let last: ReturnType<typeof MsgPackCodec.decode> | undefined;
    within(budget(bytes.length, NodeHeavyMBps), () => {
      for (let k = 0; k < Repeats; k++) last = MsgPackCodec.decode(bytes);
    });
    expect(last?.node).toBeDefined();
    expect(hashNode(last!.node!)).toBe(hashNode(node));
  });

  itWithin('copies unaligned columns at MinMBps', (within) => {
    const f64 = new Float64Array(1_000_000).map((_, i) => i * 0.5);
    const i32 = new Int32Array(1_000_000).map((_, i) => i * 3);
    const node: Node = { kind: 'doc', params: { f64, i32 } };
    const bytes = MsgPackCodec.encode(node);
    // a document placed so that the first column's payload is aligned to 8 plus one byte: neither column is in place
    // (the second payload follows the first by its 6-byte ext32 header)
    const payload = bytes.length - f64.byteLength - 6 - i32.byteLength;
    const shift = (8 - (payload % 8) + 1) % 8;
    const host = new Uint8Array(bytes.length + 16);
    host.set(bytes, shift);
    const input = host.subarray(shift, shift + bytes.length);
    let last: ReturnType<typeof MsgPackCodec.decode> | undefined;
    within(budget(bytes.length), () => {
      for (let k = 0; k < Repeats; k++) last = MsgPackCodec.decode(input);
    });
    const out = last!.node!.params!;
    expect((out.f64 as Float64Array).buffer).not.toBe(host.buffer);
    expect((out.i32 as Int32Array).buffer).not.toBe(host.buffer);
    expect(hashNode(last!.node!)).toBe(hashNode(node));
  });
});
