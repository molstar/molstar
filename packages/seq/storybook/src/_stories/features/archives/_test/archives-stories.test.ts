/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { archiveReader } from '@molstar/seqx-core/data';
import { ExampleViews } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { HostileArchives, printRefusal, saveView } from '../archive-model.js';

// The Archives stories, headless: a view saved and drawn from its archive alone as the story composes it (the default
// viewer spec, the archive's reader as the assets port), and the refusals the story prints.

const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
});

const view = (dataset: string, id?: string) => ExampleViews.find((v) => v.dataset === dataset && (!id || v.id === id))!;

describe('Features/Archives, headless', () => {
  it.each([
    ['plasmid', undefined],
    ['p53', 'protein-on-gene'],
  ] as const)('%s: saved, listed and drawn from the archive alone, as from the fixtures', async (dataset, id) => {
    const v = view(dataset, id);
    const saved = await saveView(v.build(), nodeFixtureAssets().reader, { keep: (u) => u.startsWith('generated:') });
    expect(saved.listing.split('\n')[0]).toMatch(/^\d+ bytes as zip$/);
    expect(saved.listing).toContain('index.json');
    expect(saved.listing).toContain('fixtures/');
    const open = async (
      spec: ReturnType<typeof v.build>,
      assets: () => Parameters<typeof createSeqViewer>[0]['ports'],
    ) => {
      const base = DefaultSeqViewerSpec();
      const viewer = createSeqViewer({
        ...base,
        ui: undefined,
        registry: [...(base.registry ?? []), ...(v.registry?.() ?? [])],
        ports: assets(),
      });
      viewers.push(viewer);
      viewer.resize(720);
      await viewer.load(spec);
      await viewer.settle();
      viewer.draw();
      return viewer.renderSvg({ idPrefix: 'a-' });
    };
    const fromFixtures = await open(v.build(), () => ({ assets: () => nodeFixtureAssets() }));
    const fromArchive = await open(saved.archive.spec, () => ({
      assets: () => ({ reader: archiveReader(saved.archive) }),
    }));
    expect(fromArchive.issues).toEqual([]);
    expect(fromArchive.settled).toBe(true);
    expect(fromArchive.svg).toBe(fromFixtures.svg);
  });

  it.each(HostileArchives.map((h) => [h.id, h] as const))(
    'EC-91: %s is refused, as the story prints it',
    async (id, h) => {
      const { read, text } = await printRefusal(h);
      expect(read.archive).toBeUndefined();
      expect(text).toContain('refused');
      expect(text).toContain(id === 'above-root' || id === 'absolute' ? 'seq.archive.name' : 'seq.archive.bounds');
    },
  );
});
