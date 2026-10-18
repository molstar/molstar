/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { added, common, productionLines } from '../../../../scripts/ported-lines.js';

const source = `import type { Space } from "@molstar/seq-core/model";

/** Docs. */
export function f(s: Space) {
  /* a
     block */
  const a = 1;
  return a;
}
`;

describe('seqstar-ported-lines', () => {
  it('counts production lines: not blank, not comments', () => {
    expect(productionLines(source)).toHaveLength(5);
  });

  it('takes the longest common subsequence of lines', () => {
    expect(common(['a', 'b', 'c', 'd'], ['b', 'x', 'd'])).toBe(2);
    expect(common([], ['a'])).toBe(0);
  });

  it('counts the lines a base lacks as new, and a file with no base as all new', () => {
    const next = source.replace('return a;', 'return a + 1;\n  // note\n  log(a);');
    expect(added(next, source)).toBe(2);
    expect(added(source, '')).toBe(5);
  });
});
