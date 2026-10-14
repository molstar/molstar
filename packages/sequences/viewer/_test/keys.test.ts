/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Keys hold every input a stage reads (VC1): a same-length residue edit reruns what reads residues and keeps the rest,
// a numbering edit reruns no letters, a mapping edit re-places the layers placed through it. Each case also ran against
// the old stack until SW1 (its parity layer, deleted).

import { Paths } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { mappedDoc, patch, refsRan, residuesDoc, viewer } from './graph-fixtures.js';

const set = (ref: string, params: Record<string, unknown>) => patch([{ op: 'set', ref, params }]);
const stateOf = (ctx: ReturnType<typeof viewer>['ctx'], ref: string) => {
  const v = ctx.graph.assembled.value?.view;
  const at = v?.resolved.layers.ref.findIndex((r) => v.set.strings.get(r) === ref) ?? -1;
  return v?.prepared[at]?.state;
};

describe('content keys', () => {
  it("a same-length residue edit reruns the letters' prepare and keeps the features layer's state", async () => {
    const { ctx, prepared } = viewer();
    await ctx.setDocument(residuesDoc());
    const [letters0, feat0] = [stateOf(ctx, 'v.letters.l'), stateOf(ctx, 'v.feat.l')];
    expect((await ctx.applyPatch(set('x', { value: 'MKVLATGHEW' }))).ok).toBe(true);
    expect(refsRan(ctx).filter((r) => r.startsWith('#prepare:'))).toEqual(['#prepare:v.letters.l']);
    expect(prepared.get('v.letters.l')).toBe(2);
    expect(stateOf(ctx, 'v.letters.l')).not.toBe(letters0);
    expect(stateOf(ctx, 'v.feat.l')).toBe(feat0); // the same object: `content: false`
    expect(refsRan(ctx)).not.toContain('v.letters.l'); // the layer stage reads no residues
  });

  it('a numbering edit rebuilds the molecule but not its residues: no letters prepare reruns', async () => {
    const { ctx, prepared } = viewer();
    await ctx.setDocument(residuesDoc());
    await ctx.applyPatch(set('x.auth', { first: 201 }));
    expect(refsRan(ctx)).toContain('x');
    expect(prepared.get('v.letters.l')).toBe(1);
  });

  it('an edit of a mapping re-places the layers placed through it', async () => {
    const { ctx, prepared } = viewer();
    await ctx.setDocument(mappedDoc());
    const forward = () => {
      const p = ctx.graph.assembled.value?.view?.resolved.paths[0];
      return p ? Array.from({ length: p.length }, (_, i) => Paths.at(p, i)) : [];
    };
    expect(forward().slice(0, 3)).toEqual([0, 1, 2]);
    await ctx.applyPatch(set('m', { segments: [[0, 2, 4]] }));
    expect(refsRan(ctx)).toContain('v.site.l');
    expect(prepared.get('v.site.l')).toBe(2);
    expect(forward().slice(0, 3)).toEqual([2, 3, 4]);
  });
});
