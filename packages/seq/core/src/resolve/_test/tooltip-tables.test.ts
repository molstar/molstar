/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tooltips' tables at resolve (_spec/README.md "Resolve", Tooltip tables): what loads behind the tooltips' door
// (a spec without a tooltip loads nothing; one with waits for it), and what the check cannot say: an MVS table off its
// schema (EC-252), an `address` or MVS table from a source (not read yet: said, never silent), a table read whole past
// `tooltipRows`.

import { createSeqSpec, DefaultSpecBounds, TooltipSchemas } from '@molstar/seqx-core/lang';
import { beforeAll, describe, expect, it } from 'vitest';
import { resolveSpec } from '../model.js';
import { hasTooltip, loadTooltips, tooltipsLoading } from '../tooltip-door.js';
import { MvsSchemaColumns } from '../tooltip-tables.js';

beforeAll(() => loadTooltips());

const tipCodes = (spec: ReturnType<ReturnType<typeof createSeqSpec>['build']>) =>
  resolveSpec(spec)
    .issues.filter((i) => i.code.startsWith('seq.resolve.tooltip'))
    .map((i) => [i.code, i.address]);

describe('tooltip tables at resolve', () => {
  it('a spec with a tooltip is seen (under views, under a view); one without loads nothing', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    expect(hasTooltip(b.build())).toBe(false);
    expect(tooltipsLoading(b.build())).toBeUndefined();
    b.views.view('v', (v) => {
      v.section('a', { of: s });
      v.tooltip('t', { on: s, text: 'x' });
    });
    expect(hasTooltip(b.build())).toBe(true);
  });

  it('every MVS schema has its columns', () => {
    expect(TooltipSchemas.length).toBe(Object.keys(MvsSchemaColumns).length + 2);
  });

  it('an MVS table off its schema is said; MVS’s schemas need one of their columns, others ignored (EC-252)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' });
    b.views.tooltip('b', {
      target: s,
      schema: 'chain',
      rows: { label_asym_id: ['A'], beg_label_seq_id: [3], tooltip: ['x'] } as never,
    });
    b.views.tooltip('c', { target: s, schema: 'residue', rows: { start: [0], end: [1], tooltip: ['x'] } });
    expect(tipCodes(b.build())).toEqual([['seq.resolve.tooltip.schema', '/views/c']]);
  });

  it('an address or MVS table from a source is said, not silent; a ranges one is read', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' });
    const src = b.data.source('notes', { format: 'bed', uri: 'notes.bed' });
    b.views.tooltip('a', { target: s, source: src, schema: 'address' });
    b.views.tooltip('m', { target: s, source: src, schema: 'residue' });
    b.views.tooltip('r', { target: s, source: src });
    expect(tipCodes(b.build())).toEqual([
      ['seq.resolve.tooltip.source', '/views/a'],
      ['seq.resolve.tooltip.source', '/views/m'],
    ]);
  });

  it('a table read whole past `tooltipRows` is said (read to the bound)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' });
    const n = DefaultSpecBounds.tooltipRows + 1;
    b.views.tooltip('a', {
      target: s,
      schema: 'address',
      rows: { address: Array.from({ length: n }, () => '(range :from 1 :to 2)'), tooltip: ['x'] } as never,
    });
    expect(tipCodes(b.build())).toEqual([['seq.resolve.tooltip.rows', '/views/a']]);
  });
});
