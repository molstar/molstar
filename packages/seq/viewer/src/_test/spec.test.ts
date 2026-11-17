/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashNode } from '@molstar/core/repr/tree';
import { HistoryModel } from '@molstar/ui/history/history-model';
import { describe, expect, it } from 'vitest';
import { SpecModel, SpecRunOptions } from '../models/spec.js';
import { proteinOnGene } from './fixture.js';

const person = { kind: 'person' };

describe('SpecModel', () => {
  it('holds the root, its revision as its hash, and the static model of what is shown', () => {
    const spec = new SpecModel(proteinOnGene());
    expect(spec.revision.value).toBe(hashNode(spec.root.value));
    expect(spec.model.value.systems.has('/data/g')).toBe(true);
    // an equal tree is no change, and keeps the revision
    expect(spec.apply(proteinOnGene())).toBe(false);
    expect(spec.apply(proteinOnGene({ title: 'other' }))).toBe(true);
    expect(spec.revision.value).toBe(hashNode(proteinOnGene({ title: 'other' })));
  });

  it('refuses a change made against another revision (EC-81)', () => {
    const spec = new SpecModel(proteinOnGene());
    const before = spec.revision.value;
    spec.apply(proteinOnGene({ title: 'a' }));
    expect(spec.apply(proteinOnGene({ title: 'b' }), { against: before })).toBe('stale-revision');
    expect(spec.apply(proteinOnGene({ title: 'b' }), { against: spec.revision.value })).toBe(true);
  });

  it('shows a preview over the root, which cancels to nothing or commits as one change', () => {
    const spec = new SpecModel(proteinOnGene());
    const root = spec.root.value;
    spec.setPreview(proteinOnGene({ strand: '-' }));
    expect(spec.shown.value).not.toBe(root);
    expect(spec.root.value).toBe(root);
    expect(spec.model.value.spec).toBe(spec.shown.value);
    spec.cancelPreview();
    expect(spec.shown.value).toBe(root);
    spec.setPreview(proteinOnGene({ strand: '-' }));
    expect(spec.commitPreview()).toBe(true);
    expect(spec.preview.value).toBeUndefined();
    expect(spec.root.value).not.toBe(root);
  });

  it('undo through the history restores the spec’s root, one step a batch', async () => {
    const spec = new SpecModel(proteinOnGene());
    const history = new HistoryModel(spec);
    const r0 = spec.root.value;
    await history.run('Retitle', person, async () => spec.apply(proteinOnGene({ title: 'one' })), SpecRunOptions);
    const r1 = spec.root.value;
    // a preview during a gesture ends in one step
    spec.setPreview(proteinOnGene({ title: 'dragging' }));
    await history.run('Resize', person, async () => spec.commitPreview(), SpecRunOptions);
    expect(spec.root.value).not.toBe(r1);
    expect((await history.undo()).ok).toBe(true);
    expect(spec.root.value).toBe(r1);
    expect((await history.undo()).ok).toBe(true);
    expect(spec.root.value).toBe(r0);
    expect(spec.revision.value).toBe(hashNode(r0));
    // nothing changed: no step
    const none = await history.run('Nothing', person, async () => spec.apply(r0), SpecRunOptions);
    expect(none.step).toBeUndefined();
    expect((await history.undo()).ok).toBe(false);
    await history.redo();
    expect(spec.root.value).toBe(r1);
    // a refused change fails and is no step
    const refused = await history.run('Stale', person, async () => spec.apply(r0, { against: 'old' }), SpecRunOptions);
    expect(refused.failed).toBe(true);
    history.dispose();
  });

  it('a change outside the history makes its earlier steps unsafe: undo never restores over it', async () => {
    const spec = new SpecModel(proteinOnGene());
    const history = new HistoryModel(spec);
    spec.attachHistory(history);
    await history.run(
      'Retitle',
      person,
      async (context) => spec.apply(proteinOnGene({ title: 'one' }), { context }),
      SpecRunOptions,
    );
    // inside a run nothing is invalidated
    expect((await history.undo()).ok).toBe(true);
    await history.redo();
    const outside = proteinOnGene({ title: 'outside' });
    spec.apply(outside);
    const undo = await history.undo();
    expect(undo.ok).toBe(false);
    expect(!undo.ok && undo.code).toBe('unsafe');
    expect(spec.root.value).toBe(outside);
    history.dispose();
  });
});
