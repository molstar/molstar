/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Story steps: `/snapshots/<id>`, each a stored patch in `seq.patch` with no base. The builder writes them, the checks
// read them through the slot's language, and their addresses are refs of the spec: a rename rewrites them.

import { applyPatch, patchFromNode } from '@molstar/core/repr/patch';
import { type Node, resolveAddress } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { checkSeqSpec } from '../check.js';
import { SeqLanguage, SeqPatchLanguage, SeqValidateOptions } from '../language.js';
import { P53_LETTERS } from './p53.js';

function withSteps(): Node {
  const b = createSeqSpec({ title: 'steps' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  b.views.view('main', (v) => v.section('protein', { of: p53 }).track('letters').layer('letters', { type: 'letters' }));
  b.snapshots.snapshot('zoom', {
    label: 'The start',
    patch: { ops: [{ op: 'set', at: ['views', 'main', 'protein'], key: 'range', value: [0, 10] }] },
  });
  return b.build();
}

describe('story steps', () => {
  it('are written by the builder, in the stored patches’ language, and check clean', () => {
    const spec = withSteps();
    expect(SeqPatchLanguage.id).toBe('seq.patch');
    expect(checkSeqSpec(spec)).toEqual([]);
    const step = resolveAddress(spec, ['snapshots', 'zoom']);
    expect(step?.kind).toBe('snapshot');
    expect(patchFromNode(step!.params!.patch as Node).patch?.ops).toHaveLength(1);
  });

  it('refuse a base: a step applies over whatever the spec is', () => {
    expect(() => createSeqSpec().snapshots.snapshot('s', { patch: { base: '0'.repeat(32), ops: [] } })).toThrow(
      /no base/,
    );
  });

  it('follow a rename of what they name', () => {
    const r = applyPatch(
      withSteps(),
      {
        ops: [
          {
            op: 'move',
            at: ['views', 'main', 'protein'],
            parent: ['views', 'main'],
            place: { last: true },
            rename: 'p',
          },
        ],
      },
      { language: SeqLanguage, validate: SeqValidateOptions },
    );
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    const step = patchFromNode(resolveAddress(r.root, ['snapshots', 'zoom'])!.params!.patch as Node).patch!;
    expect(step.ops[0]).toMatchObject({ op: 'set', at: ['views', 'main', 'p'] });
  });
});
