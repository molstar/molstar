/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeMsgPackTree, encodeMsgPackTree } from '@molstar/core/repr/codecs/msgpack';
import { hashNode, nodesEqual } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { ExampleViews } from '../views.js';

// The views' specs through MessagePack (_spec/README.md; format.md, Encodings): each example view's spec encodes to
// bytes and decodes to an equal tree (the same hash), whichever columns it holds; a quick representative and the sweep
// over every view.

const Representatives = ['p53/tracks', 'globins/family', 'antibody/schemes'];

function roundTrip(view: (typeof ExampleViews)[number]) {
  const spec = view.build();
  const bytes = encodeMsgPackTree(spec);
  const read = decodeMsgPackTree(bytes);
  expect(read.issues).toEqual([]);
  expect(read.node).toBeDefined();
  expect(nodesEqual(read.node!, spec)).toBe(true);
  expect(hashNode(read.node!)).toBe(hashNode(spec));
  // encoding the decoded tree gives the same bytes
  expect(encodeMsgPackTree(read.node!)).toEqual(bytes);
}

describe('the example views through MessagePack', () => {
  it('there are views to check', () => {
    expect(ExampleViews.length).toBeGreaterThan(50);
    for (const r of Representatives) expect(ExampleViews.some((v) => `${v.dataset}/${v.id}` === r)).toBe(true);
  });

  for (const view of ExampleViews) {
    const name = `${view.dataset}/${view.id}`;
    it(`${name}: its spec encodes to MessagePack and decodes equal`, () => roundTrip(view));
  }
});
