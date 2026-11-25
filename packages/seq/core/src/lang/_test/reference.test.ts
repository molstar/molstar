/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A values layer's `reference` param (language 0.11.0; _spec/README.md "Render", Layer types, `values`): up to four
// lines `{ at, label?, color?, style? }`, each checked by the description and by the layer's own rule for a color.

import { parseJson, printJson } from '@molstar/core/repr/codecs/json';
import { validateAndNormalize } from '@molstar/core/repr/schema';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { checkSeqSpec } from '../check.js';
import { SeqLanguage, SeqLanguageVersion } from '../language.js';

function withReference(reference: unknown) {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
  const v = b.data.values('v', { target: s, values: Float32Array.from([1, 2, 3, 4, 5, 6, 7, 8]) });
  b.views.view('main', (view) =>
    view
      .section('s', { of: s })
      .track('v')
      .layer('values', { type: 'values', data: v, ...(reference === undefined ? {} : { reference }) } as never),
  );
  return b.build();
}

const codes = (reference: unknown) => checkSeqSpec(withReference(reference)).map((i) => [i.code, i.severity]);

describe('a values layer’s reference lines (language 0.11.0)', () => {
  it('is a param of the layer, in the version 0.11.0', () => {
    expect(SeqLanguageVersion).toBe('0.11.0');
    expect(Object.keys(SeqLanguage.kinds.layer.params)).toContain('reference');
  });

  it('takes up to four lines, each at a value, with a label, a color and a style', () => {
    expect(codes(undefined)).toEqual([]);
    expect(codes([])).toEqual([]);
    expect(codes([{ at: 7.3 }])).toEqual([]);
    expect(
      codes([
        { at: 7.3, label: 'genome-wide', style: 'dashed' },
        { at: 5, label: 'suggestive', style: 'dotted', color: '#d62728' },
        { at: 0, style: 'solid', color: 'token:msui-text-muted' },
        { at: 3, color: 'light-dark(#333333, #dddddd)' },
      ]),
    ).toEqual([]);
  });

  it('refuses a fifth line, a line with no `at`, an `at` that is no number, a style outside the words', () => {
    const line = { at: 1 };
    expect(codes([line, line, line, line, line]).length).toBeGreaterThan(0);
    expect(codes([{ label: 'x' }]).length).toBeGreaterThan(0);
    expect(codes([{ at: 'high' }]).length).toBeGreaterThan(0);
    expect(codes([{ at: 1, style: 'wavy' }]).length).toBeGreaterThan(0);
    expect(codes([{ at: 1, label: 4 }]).length).toBeGreaterThan(0);
  });

  it('refuses a color that is no color or token, and a theme form in its place (a line is one colour)', () => {
    expect(codes([{ at: 1, color: 'reddish' }])).toEqual([['seq.lang.theme.color', 'error']]);
    expect(codes([{ at: 1, color: 'token:Not A Token' }])).toEqual([['seq.lang.theme.color', 'error']]);
    expect(codes([{ at: 1, color: { kind: 'uniform' } }]).length).toBeGreaterThan(0);
  });

  it('round trips through JSON', () => {
    const reference = [
      { at: 7.3, label: 'genome-wide' },
      { at: 5, label: 'suggestive', style: 'dotted', color: '#d62728' },
    ];
    const spec = withReference(reference);
    const read = parseJson(printJson(spec, { language: SeqLanguage }));
    expect(read.issues).toEqual([]);
    const back = validateAndNormalize(SeqLanguage, read.node as Node);
    expect(back.issues).toEqual([]);
    expect(hashNode(back.node as Node)).toBe(hashNode(spec));
    expect(checkSeqSpec(back.node as Node)).toEqual([]);
  });
});
