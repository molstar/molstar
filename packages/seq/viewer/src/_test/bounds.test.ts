/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bounds as the viewer's composition holds them (_spec/README.md "Configuration"): items with the security
// table's defaults, the patch's nesting check and the store's table bounds reading them.

import { StoryPart } from '../story.js';
import { applyPatch, type Patch, patchToNode } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { DefaultViewerConfig, viewerConfig, viewerConfigOf, viewerItem } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { StoreModel } from '../models/store.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer } from '../viewer.js';
import { proteinOnGene } from './fixture.js';
import { installFrames } from './frames.js';

const hideLetters = (): Patch => ({
  ops: [{ op: 'set', at: ['views', 'main', 'protein', 'letters', 'l'], key: 'hidden', value: true }],
});

describe('the bounds items', () => {
  it('have the security table’s defaults', () => {
    expect(DefaultViewerConfig.bounds).toEqual({
      maxBytes: 64 * 1024 * 1024,
      maxDepth: 256,
      maxQueryDepth: 64,
      maxColumnLength: 2 ** 31 - 1,
      maxNodeIssues: 100,
      maxPlaces: 10,
      themeDepth: 4,
      maxBlockValues: 2 ** 24,
    });
  });

  it('are set as items of a viewer spec', () => {
    const config = viewerConfigOf([
      [viewerItem('bounds.maxDepth'), 4],
      [viewerItem('bounds.maxBlockValues'), 100],
    ]);
    expect(config.bounds).toMatchObject({ maxDepth: 4, maxBlockValues: 100, maxBytes: 64 * 1024 * 1024 });
  });
});

describe('a patch under the nesting bound', () => {
  it('is applied at the bound and refused one level past it, as the whole spec’s check does', () => {
    // spec / views / main / protein / letters / l: the layer is the sixth level, the root the zeroth
    const at = new SpecModel(
      proteinOnGene(),
      {},
      viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxDepth: 5 } }),
    );
    const ok = at.applyPatch(hideLetters());
    expect(ok.ok && ok.changed).toBe(true);
    const past = new SpecModel(
      proteinOnGene(),
      {},
      viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxDepth: 4 } }),
    );
    const before = past.head.value.revision;
    const refused = past.applyPatch(hideLetters());
    expect(refused.ok).toBe(false);
    expect(refused.issues.map((i) => i.code)).toContain('seq.lang.bound.depth');
    expect(past.head.value.revision).toBe(before);
  });

  it('is applied by default', () => {
    const spec = new SpecModel(proteinOnGene());
    expect(spec.applyPatch(hideLetters()).ok).toBe(true);
  });
});

describe('the store’s table bounds', () => {
  it('are the composition’s', () => {
    const store = new StoreModel(viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxBlockValues: 7 } }));
    expect(store.store.settings).toMatchObject({ maxBlockValues: 7, maxColumnLength: 2 ** 31 - 1 });
  });
});

describe('a whole root under the bounds', () => {
  const root = proteinOnGene();
  const deep = (n: number): Node => {
    let node: Node = { kind: 'ext:note', id: 'c' };
    for (let i = 1; i < n; i++) node = { kind: 'ext:note', id: 'c', children: [node] };
    return { ...root, children: [...(root.children ?? []), node] };
  };
  const withQuery = (n: number): Node => ({
    ...root,
    children: [
      ...(root.children ?? []),
      { kind: 'ext:note', id: 'q', params: { to: `/data/${'(a :b '.repeat(n)}x${')'.repeat(n)}` } as never },
    ],
  });
  const bigText = (n: number): Node => ({
    ...root,
    children: [...(root.children ?? []), { kind: 'ext:note', id: 'big', params: { s: 'x'.repeat(n) } as never }],
  });

  it('SpecModel.apply refuses what is past a bound, before it hashes, with the issues and the root as it was', () => {
    const spec = new SpecModel(root, {}, viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxBytes: 10_000 } }));
    const head = spec.head.value;
    expect(spec.apply(deep(20_000))).toBe('out-of-bounds');
    expect(spec.boundsIssues.value.map((i) => i.code)).toEqual(['seq.lang.bound.depth']);
    expect(spec.apply(bigText(20_000))).toBe('out-of-bounds');
    expect(spec.boundsIssues.value.map((i) => i.code)).toEqual(['seq.lang.bound.size']);
    expect(spec.apply(withQuery(65))).toBe('out-of-bounds');
    expect(spec.boundsIssues.value.map((i) => i.code)).toEqual(['seq.lang.bound.query_depth']);
    expect(spec.head.value).toBe(head);
    // one within the bounds is taken, and the issues are cleared
    expect(spec.apply(withQuery(64))).toBe(true);
    expect(spec.boundsIssues.value).toEqual([]);
  });

  it('a viewer’s load rejects with the issues', async () => {
    const frames = installFrames();
    try {
      const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
      await viewer.load(root);
      for (const bad of [deep(20_000), withQuery(65)]) {
        const e = await viewer.load(bad).then(
          () => undefined,
          (x: unknown) => x as { issues: readonly { code: string }[] },
        );
        expect(e?.issues.map((i) => i.code)).toEqual([expect.stringMatching(/^seq\.lang\.bound\./)]);
      }
      expect(viewer.spec.root.value).toBe(viewer.spec.head.value.root);
      viewer.retire();
    } finally {
      frames.restore();
    }
  });
});

describe('previews and story steps under the bounds', () => {
  const stepRoot = (): Node => {
    // two ops: a track, then a layer in it: past a depth of 4 only together
    const step: Patch = {
      ops: [
        {
          op: 'insert',
          parent: ['views', 'main', 'protein'],
          place: { last: true },
          node: { kind: 'track', id: 'note' },
        },
        {
          op: 'insert',
          parent: ['views', 'main', 'protein', 'note'],
          place: { last: true },
          node: { kind: 'layer', id: 'n', params: { type: 'letters' } },
        },
      ],
    };
    const r = applyPatch(proteinOnGene(), {
      ops: [
        {
          op: 'insert',
          parent: [],
          place: { last: true },
          node: {
            kind: 'snapshots',
            id: 'snapshots',
            children: [{ kind: 'snapshot', id: 'deep', params: { label: 'Deep', patch: patchToNode(step) } }],
          },
        },
      ],
    });
    if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
    return r.root;
  };

  it('a stored step whose ops build past maxDepth does not apply; within it, it does', () => {
    const tight = new SpecModel(
      stepRoot(),
      {},
      viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxDepth: 4 } }),
    );
    tight.attachStory(StoryPart);
    expect(tight.enterStep('deep')).toBeUndefined();
    expect(tight.step.value).toMatchObject({ id: 'deep', applies: false });
    expect(tight.issues.value.map((i) => i.code)).toEqual(['seq.viewer.step.stale']);
    const loose = new SpecModel(stepRoot());
    loose.attachStory(StoryPart);
    expect(loose.enterStep('deep')).toBeUndefined();
    expect(loose.step.value).toMatchObject({ id: 'deep', applies: true });
  });

  it('a preview that builds past maxDepth is dropped with an issue', () => {
    const spec = new SpecModel(
      proteinOnGene(),
      {},
      viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxDepth: 4 } }),
    );
    const preview: Patch = {
      ops: [
        {
          op: 'insert',
          parent: ['views', 'main', 'protein'],
          place: { last: true },
          node: { kind: 'track', id: 'note', children: [{ kind: 'layer', id: 'n', params: { type: 'letters' } }] },
        },
      ],
    };
    expect(spec.setPreview(preview).ok).toBe(false);
    expect(spec.preview.value).toBeUndefined();
  });
});
