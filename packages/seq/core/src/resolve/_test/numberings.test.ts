/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Declared numberings (_spec/README.md, Language, Numberings; "Resolve", Numberings): a sequence with Kabat and IMGT
// numberings as `numbering` children, read by the query forms (a label in each, a range across insertion codes, IMGT
// out of index order, a missing end snapped, an unknown label hinted), printed back by `addressFor`, checked by the
// language (a label list's length, a duplicate, an unknown scheme word accepted), and a change to the sequence's
// letters reaching its numbering's check.

import { applyPatch } from '@molstar/core/repr/patch';
import { describe, expect, it } from 'vitest';
import { checkSeqChanges } from '@molstar/seqx-core/lang/changes';
import { checkSeqSpec, createSeqSpec, SeqLanguage, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { addressFor, type QueryResult, rangesLocation, resolveAddressIn } from '@molstar/seqx-core/query';
import { resolveSpec } from '../model.js';
import { antibodySpec, H } from './antibody.js';

const root = antibodySpec();
const model = resolveSpec(root);

function at(address: string) {
  const r = resolveAddressIn(root, model, address);
  const loc = r.result?.type === 'location' ? r.result.location : undefined;
  return { ranges: loc ? Array.from(loc.parts[0]?.ranges ?? []) : undefined, codes: r.issues.map((i) => i.code), r };
}

const location = (start: number, end: number): QueryResult => ({
  type: 'location',
  location: rangesLocation('/data/h', [start, end]),
});

describe('declared numberings: the model', () => {
  it('lists them beside the built-ins, and holds each with its scheme', () => {
    expect(model.systems.get('/data/h')?.numberings).toEqual(['index', 'ordinal', 'kabat', 'imgt']);
    expect(model.numberingOf('/data/h', 'kabat')).toMatchObject({ address: '/data/h/kabat', scheme: 'kabat' });
    expect(model.numberings.get('/data/h/imgt')?.numbering?.labelAt(16)).toBe('111.1');
    expect(model.system('/data/h')?.numberingIds()).toEqual(['kabat', 'imgt']);
    // a position with no label
    expect(model.system('/data/h')?.numbering('kabat')?.labelAt(0)).toBeUndefined();
    expect(model.issues).toEqual([]);
  });
});

describe('declared numberings: the query forms', () => {
  it('(at :label …) in each numbering', () => {
    expect(at('/data/h/(at :label 82B :numbering kabat)').ranges).toEqual([12, 13]);
    expect(at('/data/h/(at :label 82 :ins C :numbering kabat)').ranges).toEqual([13, 14]);
    expect(at('/data/h/(at :label "111.1" :numbering imgt)').ranges).toEqual([16, 17]);
    expect(at('/data/h/(at :label 1 :numbering kabat)').ranges).toEqual([1, 2]);
    expect(at('/data/h/(at :label 1 :numbering imgt)').ranges).toEqual([0, 1]);
    // the built-ins as they were
    expect(at('/data/h/(at :label 1)').ranges).toEqual([0, 1]);
  });

  it('a range across insertion codes', () => {
    expect(at('/data/h/(range :from 82 :to 83 :numbering kabat)').ranges).toEqual([10, 15]);
    expect(at('/data/h/(range :from 82A :to 82C :numbering kabat)').ranges).toEqual([11, 14]);
    expect(at('/data/h/(range :from 100 :to 100B :numbering kabat)').ranges).toEqual([15, 18]);
  });

  it('IMGT out of index order: ends found by label, the range in index order', () => {
    expect(at('/data/h/(range :from "111.1" :to "112.1" :numbering imgt)').ranges).toEqual([16, 18]);
    expect(at('/data/h/(range :from 112 :to 113 :numbering imgt)').ranges).toEqual([18, 20]);
    // 112.1 is after 112 by label, before it by index: still a range, in index order
    expect(at('/data/h/(range :from "112.1" :to 112 :numbering imgt)').ranges).toEqual([17, 19]);
    // 112 to 112.1: backwards in index order
    expect(at('/data/h/(range :from 112 :to "112.1" :numbering imgt)').codes).toContain('seq.query.range.backwards');
  });

  it('a missing end snaps inward with an info; an unknown single label is reported with the labels near it', () => {
    const snapped = at('/data/h/(range :from 50 :to 60 :numbering kabat)');
    expect(snapped.ranges).toEqual([6, 10]);
    expect(snapped.codes).toEqual(['seq.query.label.snapped']);
    // past the labels each side: the system's own ends
    expect(at('/data/h/(range :from 110 :numbering imgt)').ranges).toEqual([15, 20]);
    const unknown = at('/data/h/(at :label 99 :numbering kabat)');
    expect(unknown.r.status).toBe('failed');
    expect(unknown.codes).toEqual(['seq.query.label.unknown']);
    expect(unknown.r.issues[0].hint).toMatch(/83 100/);
  });

  it('an unknown numbering lists the system’s', () => {
    const r = at('/data/h/(at :label 1 :numbering chothia)');
    expect(r.codes).toEqual(['seq.query.numbering.unknown']);
    expect(r.r.issues[0].hint).toBe('numberings: index, ordinal, kabat, imgt');
  });
});

describe('declared numberings: an address printed back', () => {
  it('in the numbering asked, reading back to the same positions', () => {
    expect(addressFor(root, model, location(12, 13), { numbering: 'kabat' })).toBe(
      '/data/h/(at :label 82B :numbering kabat)',
    );
    expect(addressFor(root, model, location(10, 15), { numbering: 'kabat' })).toBe(
      '/data/h/(range :from 82 :numbering kabat :to 83)',
    );
    expect(addressFor(root, model, location(16, 18), { numbering: 'imgt' })).toBe(
      '/data/h/(range :from "111.1" :numbering imgt :to "112.1")',
    );
  });

  it('in ordinals where a label is missing there', () => {
    expect(addressFor(root, model, location(0, 3), { numbering: 'kabat' })).toBe('/data/h/(range :from 1 :to 3)');
  });
});

describe('declared numberings: the language', () => {
  it('a valid spec has no issues; a scheme is an open word', () => {
    expect(checkSeqSpec(root)).toEqual([]);
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'protein', letters: 'ME' }, (s) =>
      s.numbering('house', { method: 'labels', labels: ['A1', 'A2'], scheme: 'my-lab' }),
    );
    expect(b.check()).toEqual([]);
  });

  it('a label list of another length is an error at its labels; a duplicate a warning at its second place', () => {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' }, (s) => {
      s.numbering('short', { method: 'labels', labels: ['1', '2'] });
      s.numbering('twice', { method: 'labels', labels: ['1', '2', '2', null] });
    });
    const issues = b.check();
    expect(issues.map((i) => [i.code, i.address, i.path])).toEqual([
      ['seq.lang.numbering.length', '/data/s/short', ['children', 0, 'children', 0, 'children', 0, 'params', 'labels']],
      [
        'seq.lang.numbering.duplicate',
        '/data/s/twice',
        ['children', 0, 'children', 0, 'children', 1, 'params', 'labels', 2],
      ],
    ]);
  });

  it('a change to the sequence’s letters reaches its numbering’s check', () => {
    const r = applyPatch(
      root,
      { ops: [{ op: 'set', at: ['data', 'h'], key: 'letters', value: H.slice(0, 10) }] },
      { language: SeqLanguage, validate: SeqValidateOptions },
    );
    if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
    const issues = checkSeqChanges(root, r.root, r.changes);
    expect(issues.map((i) => [i.code, i.address])).toEqual([
      ['seq.lang.numbering.length', '/data/h/kabat'],
      ['seq.lang.numbering.length', '/data/h/imgt'],
    ]);
  });
});
