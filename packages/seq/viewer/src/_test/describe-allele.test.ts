/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A variant allele's card words (_spec/README.md "Behaviors", hoverCard, A variant allele's card): the record's ID,
// QUAL and FILTER, INFO fields bounded in number ("+k more") and in length (an ellipsis), core's fixed fields never
// listed as INFO, and a hover held over one reading it once.

import { categoryColumn } from '@molstar/core/repr/tree';
import type { RangesList } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { AlleleInfoChars, AlleleInfoRows, alleleOf } from '../describe-allele.js';

function list(infoCount: number, long = ''): RangesList {
  const info = Array.from({ length: infoCount }, (_, k) => ({
    name: `K${k}`,
    column: categoryColumn([k === 0 && long ? long : `v${k}`]),
  }));
  return {
    form: 'list',
    start: Int32Array.of(10),
    end: Int32Array.of(11),
    fields: [
      { name: 'label', column: categoryColumn(['rs1']) },
      { name: 'name', column: categoryColumn(['rs1']) },
      { name: 'ref', column: categoryColumn(['A']) },
      { name: 'alt', column: categoryColumn(['G']) },
      { name: 'qual', column: Float32Array.of(50) },
      { name: 'filter', column: categoryColumn(['PASS']) },
      { name: 'category', column: categoryColumn(['SNV']) },
      ...info,
    ],
  };
}

describe('a variant allele’s card words', () => {
  it('lists the ID, QUAL and FILTER, the first INFO fields and "+k more"; the fixed fields are no INFO', () => {
    const a = alleleOf(list(AlleleInfoRows + 3), 0)!;
    expect(a.rows.slice(0, 3)).toEqual([
      ['ID', 'rs1'],
      ['QUAL', '50'],
      ['FILTER', 'PASS'],
    ]);
    expect(a.rows.filter(([n]) => n.startsWith('K')).length).toBe(AlleleInfoRows);
    expect(a.rows.at(-1)).toEqual(['INFO', '+3 more']);
    expect(a.info.map(([n]) => n)).not.toContain('name');
  });

  it('no "+k more" at the bound', () => {
    const a = alleleOf(list(AlleleInfoRows), 0)!;
    expect(a.rows.some(([n]) => n === 'INFO')).toBe(false);
  });

  it('bounds an INFO value’s length with an ellipsis', () => {
    const a = alleleOf(list(1, 'x'.repeat(200)), 0)!;
    const [, v] = a.info[0];
    expect(v.length).toBe(AlleleInfoChars);
    expect(v.endsWith('…')).toBe(true);
  });

  it('a hover held over one reads it once', () => {
    const t = list(2);
    expect(alleleOf(t, 0)).toBe(alleleOf(t, 0));
  });
});
