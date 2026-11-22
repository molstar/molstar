/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `tooltip` kind (_spec/README.md, Language, Tooltips): the design's builder lines check clean; by data or by a
// table, never both; a template field the table lacks (EC-251), a table off its schema (EC-252), an `on` naming nothing
// (EC-250), a template past the bound (EC-256); MVS's `{}` reads `field`; a view's `export.titles`.

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { specBounds } from '../bounds.js';
import { checkSeqSpec } from '../check.js';
import { SeqLanguageVersion } from '../language.js';
import { templateFields } from '../tooltip-language.js';

/** p53's protein, its notes and the design's four tooltips. */
function designSpec() {
  const b = createSeqSpec({ title: 'tooltips' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: 'M'.repeat(393) });
  const notes = b.data.ranges('notes', {
    target: p53,
    rows: { start: [93, 324], end: [292, 355], tooltip: ['DNA-binding', 'tetramerisation'] },
  });
  b.refs.ref('hotspots', { to: `${p53}/(range :from 175 :to 175)` });
  b.views.tooltip('notes', { on: notes, label: 'Note' });
  b.views.tooltip('dbd', { on: b.refs.ref('dbd', { to: `${p53}/(range :from 94 :to 292)` }), text: 'binds **DNA**' });
  b.views.tooltip('curated', {
    target: p53,
    schema: 'address',
    rows: { address: ['(at :label 175)'], tooltip: ['R175H, [a hotspot](!seq-focus=/refs/hotspots)'] },
  });
  b.views.tooltip('mvs', {
    target: p53,
    source: b.data.source('ann', { format: 'json', uri: 'ann.json' }),
    schema: 'auth_residue',
    field: 'note',
  });
  b.views.view('v', (v) => {
    v.section('s', { of: p53 }).track('t').layer('l', { type: 'letters' });
    v.tooltip('own', { target: p53, rows: { start: [0], end: [10], tooltip: ['the start'] } });
  });
  return b;
}

const codes = (b: ReturnType<typeof createSeqSpec>) => b.check().map((i) => i.code);

describe('the tooltip kind', () => {
  it('is the language 0.9.0’s, under views and a view; the design’s lines check clean', () => {
    expect(SeqLanguageVersion).toBe('0.9.0');
    const b = designSpec();
    expect(b.check()).toEqual([]);
    const spec = b.build();
    const views = spec.children!.find((c) => c.kind === 'views')!;
    expect(views.children!.filter((c) => c.kind === 'tooltip').map((c) => c.id)).toEqual([
      'notes',
      'dbd',
      'curated',
      'mvs',
    ]);
    const view = views.children!.find((c) => c.kind === 'view')!;
    expect(view.children!.find((c) => c.kind === 'tooltip')?.id).toBe('own');
  });

  it('is by data or by a table, never both nor neither', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.tooltip('a', { target: s });
    b.views.tooltip('b', {
      target: s,
      rows: { start: [0], end: [1], tooltip: ['x'] },
      source: b.data.source('x', { format: 'bed', uri: 'x.bed' }),
    });
    expect(codes(b)).toEqual(['seq.lang.tooltip.content', 'seq.lang.tooltip.content']);
  });

  it('a field its template names that the table lacks is said (EC-251); `{}` reads `field`', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.tooltip('a', { target: s, rows: { start: [0], end: [1], note: ['x'] } });
    b.views.tooltip('b', { target: s, rows: { start: [0], end: [1], note: ['x'] }, field: 'note' });
    const notes = b.data.ranges('n', { target: s, rows: { start: [0], end: [1], label: ['a'] } });
    b.views.tooltip('c', { on: notes, text: '{label}: {note}' });
    const issues = b.check();
    expect(issues.map((i) => [i.code, i.address])).toEqual([
      ['seq.lang.tooltip.field', '/views/a'],
      ['seq.lang.tooltip.field', '/views/c'],
    ]);
  });

  it('a table off its schema is said (Seq*’s schemas here, MVS’s at resolve: resolve’s tooltip tables, EC-252)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' });
    b.views.tooltip('a', { target: s, rows: { start: [0], tooltip: ['x'] } as never });
    b.views.tooltip('b', {
      target: s,
      schema: 'chain',
      rows: { label_asym_id: ['A'], beg_label_seq_id: [3], tooltip: ['x'] } as never,
    });
    b.views.tooltip('c', { target: s, schema: 'residue', rows: { start: [0], end: [1], tooltip: ['x'] } });
    expect(b.check().map((i) => [i.code, i.address])).toEqual([['seq.lang.tooltip.schema', '/views/a']]);
  });

  it('an `on` naming nothing is said, that tooltip inert (EC-250)', () => {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.tooltip('a', { on: ['/data/s/(range :from 1 :to 2)', '/data/gone'], text: 'x' });
    const issues = b.check();
    expect(issues.map((i) => i.code)).toEqual(['seq.lang.tooltip.on']);
    expect(issues[0].path?.slice(-2)).toEqual(['on', 1]);
  });

  it('a template past the bound, or one not an f-string, is refused (EC-256)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.tooltip('long', { target: s, rows: { start: [0], end: [1], tooltip: ['x'] }, text: 'x'.repeat(10_000) });
    b.views.tooltip('bad', { target: s, rows: { start: [0], end: [1], tooltip: ['x'] }, text: '{tooltip' });
    expect(codes(b)).toEqual(['seq.lang.tooltip.text', 'seq.lang.tooltip.text']);
    expect(checkSeqSpec(b.build(), specBounds({ tooltipText: 20_000 })).map((i) => i.address)).toEqual(['/views/bad']);
  });

  it('a CIF category is said once: not read yet (EC-252)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' });
    const cif = b.data.source('cif', { format: 'mmcif', uri: '1tup.cif' });
    b.views.tooltip('t', { target: s, source: cif, schema: 'auth_residue', block: 0, category: 'notes' });
    expect(b.check().map((i) => [i.code, i.severity])).toEqual([['seq.lang.tooltip.category', 'info']]);
  });

  it('a template’s fields, `{}` as its field (MVS’s text_format and field_name)', () => {
    expect(templateFields('{} at {start:d}: {label} {{x}}', 'note')).toEqual(['note', 'start', 'label']);
    for (const bad of ['}', '{a', '{a{b}}']) expect(typeof templateFields(bad)).toBe('string');
  });

  it('a view exports its items’ titles as their names or their cards', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.view('v', { export: { titles: 'card' } }, (v) => v.section('a', { of: s }));
    expect(b.check()).toEqual([]);
    const bad = createSeqSpec();
    bad.views.view('v', { export: { titles: 'all' as never } });
    expect(bad.check().map((i) => i.code)).toContain('repr.schema.param.enum');
  });
});
