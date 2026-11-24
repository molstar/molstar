/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A VCF allele's card words read once a hover (the viewer's _spec "Behaviors", hoverCard, A variant allele's card): a
// hover held over one allele, its words asked again (each card's describe), allocates nothing. By hand, never in CI:
// `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { categoryColumn } from '@molstar/core/repr/tree';
import type { RangesList } from '@molstar/seqx-core/data';
import { alleleOf } from '@molstar/seqx-viewer/describe-allele';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

/** A ClinVar-like row: the fixed fields and a dozen INFO fields (more than the card lists). */
const Alleles: RangesList = {
  form: 'list',
  start: Int32Array.of(7_676_000, 7_676_010),
  end: Int32Array.of(7_676_001, 7_676_011),
  fields: [
    { name: 'label', column: categoryColumn(['12345', '12346']) },
    { name: 'name', column: categoryColumn(['12345', '12346']) },
    { name: 'ref', column: categoryColumn(['C', 'G']) },
    { name: 'alt', column: categoryColumn(['T', 'A']) },
    { name: 'qual', column: Float32Array.of(Number.NaN, Number.NaN) },
    { name: 'filter', column: categoryColumn(['.', '.']) },
    { name: 'category', column: categoryColumn(['SNV', 'SNV']) },
    ...Array.from({ length: 12 }, (_, k) => ({ name: `INFO${k}`, column: categoryColumn([`a${k}`, `b${k}`]) })),
  ],
};

describe('a VCF allele’s card words', () => {
  it('a hover held over one allele: its words again', ({ task }) => {
    alleleOf(Alleles, 1);
    const { r } = quickCase(task, 'vcf allele card words, a held hover', () => void alleleOf(Alleles, 1), {});
    console.info(`vcf allele card words: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
