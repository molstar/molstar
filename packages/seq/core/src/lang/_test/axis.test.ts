/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { checkSeqChanges } from '../check-changes.js';
import { checkSeqSpec } from '../check.js';
import { SeqLanguage, SeqValidateOptions } from '../language.js';

// The axis's params (_spec/README.md, Language): a section's direction, pieces, collapsed and pad; a letters layer's
// complement and code; the references inside them checked by kind.

function spec(section: Record<string, unknown>, layer: Record<string, unknown> = {}) {
  const b = createSeqSpec();
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: 'MEEPQSDP' });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  b.data.map('codon', { from: p53, to: gene, method: 'codon', exons: { start: [2, 20], end: [14, 32] }, strand: '-' });
  b.data.ranges('exons', { target: gene, rows: { start: [2, 20], end: [14, 32] } });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('t', {}).layer('l', { type: 'letters' });
    const g = v.section('gene', { of: gene, follow: protein, ...(section as object) });
    g.track('t', {}).layer('l', { type: 'letters', ...(layer as object) });
  });
  return b.build();
}

const codes = (s: ReturnType<typeof spec>) => checkSeqSpec(s).map((i) => [i.code, i.path?.slice(-2)]);

describe('the axis params', () => {
  it('take each direction, pieces form, pad and the letters params', () => {
    for (const direction of ['genomic', 'reversed', 'follow', { along: '/data/p53' }, { along: '/data/codon' }])
      expect(codes(spec({ direction }))).toEqual([]);
    for (const pieces of ['whole', 'coding', '/data/exons', '/data/codon', ['/data/exons', '/data/codon']])
      expect(codes(spec({ pieces, collapsed: pieces === 'whole' ? 'coding' : pieces, pad: 3 }))).toEqual([]);
    expect(codes(spec({}, { complement: 'always', code: 'three' }))).toEqual([]);
  });

  it('refuse values outside them', () => {
    expect(checkSeqSpec(spec({ direction: 'up' })).length).toBeGreaterThan(0);
    expect(checkSeqSpec(spec({ pieces: 'introns' })).length).toBeGreaterThan(0);
    expect(checkSeqSpec(spec({ pad: -1 })).length).toBeGreaterThan(0);
    expect(checkSeqSpec(spec({}, { complement: 'sometimes' })).length).toBeGreaterThan(0);
  });

  it('check the kinds their references name', () => {
    expect(codes(spec({ direction: { along: '/views/main/protein' } }))).toEqual([
      ['seq.lang.ref.target', ['direction', 'along']],
    ]);
    expect(codes(spec({ pieces: ['/data/exons', '/data/p53'] }))).toEqual([['seq.lang.ref.target', ['pieces', 1]]]);
    expect(codes(spec({ collapsed: '/data/tp53' }))).toEqual([['seq.lang.ref.target', ['params', 'collapsed']]]);
  });

  it('are followed by a change: a node `along` or `pieces` names, removed, is reported for the section', () => {
    const before = spec({ direction: { along: '/data/codon' }, pieces: ['/data/exons'] });
    for (const at of [
      ['data', 'exons'],
      ['data', 'codon'],
    ]) {
      const r = applyPatch(before, { ops: [{ op: 'remove', at }] } as never, {
        language: SeqLanguage,
        validate: SeqValidateOptions,
        allowDangling: true,
      });
      if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
      const issues = checkSeqChanges(before, r.root, r.changes);
      expect(issues.filter((i) => i.address === '/views/main/gene').map((i) => i.code)).toContain(
        'seq.lang.ref.target',
      );
    }
  });
});
