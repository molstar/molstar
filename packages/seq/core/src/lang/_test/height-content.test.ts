/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { checkSeqSpec } from '../check.js';

/** A tracks section built with the given section and track `height`s. */
function spec(section: unknown, track: unknown, version?: string) {
  const b = createSeqSpec(version ? { version } : {});
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
  b.views.view('v', (v) =>
    v
      .section('sec', { of: s, height: section } as never)
      .track('t', { height: track } as never, (t) => t.layer('l', { type: 'letters' })),
  );
  return b.build();
}

const codes = (n: ReturnType<typeof spec>) => checkSeqSpec(n).map((i) => i.code);

describe('height: content (language 0.11.0)', () => {
  it('is accepted on a track and a section, beside a px number and an absent height', () => {
    expect(codes(spec('content', 'content'))).toEqual([]);
    expect(codes(spec(200, 40))).toEqual([]);
    expect(codes(spec(undefined, undefined))).toEqual([]);
  });

  it('refuses any other string', () => {
    expect(codes(spec('content', 'auto')).length).toBeGreaterThan(0);
    expect(codes(spec('auto', 'content')).length).toBeGreaterThan(0);
  });

  it('leaves a 0.10.0 spec with a number height unchanged', () => {
    expect(codes(spec(200, 40, '0.10.0'))).toEqual([]);
  });
});
