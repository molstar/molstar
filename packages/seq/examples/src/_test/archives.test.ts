/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { AssetsPort } from '@molstar/seqx-core/context';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { archiveReader, readArchive, SeqRuntime, writeArchive } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets } from '../node-assets.js';
import { type ExampleView, ExampleViews } from '../views.js';

// Every dataset saved as an archive and opened again (EC-91; _spec/README.md "Data", Archives): each view's spec with
// the files its sources name, read through the Node adapter that refuses all but the fixtures; the opened archive's
// reader serves the same bytes, and a view drawn from the archive alone is the view drawn from the fixtures.

/** A location a view’s own runtime serves (a generated chromosome): no file to hold. */
const generated = (uri: string) => uri.startsWith('generated:');

const signal = new AbortController().signal;

const urisOf = (node: Node): string[] => [
  ...(node.kind === 'source'
    ? (['uri', 'index'] as const).flatMap((k) =>
        typeof node.params?.[k] === 'string' ? [node.params[k] as string] : [],
      )
    : []),
  ...(node.children ?? []).flatMap(urisOf),
];

async function saved(view: ExampleView) {
  const spec = view.build();
  const assets = nodeFixtureAssets().reader;
  const w = await writeArchive(spec, { assets, keep: generated });
  expect(w.issues.map((i) => i.message)).toEqual([]);
  const r = await readArchive(w.bytes!);
  expect(r.issues).toEqual([]);
  return { spec, assets, archive: r.archive!, size: w.bytes!.byteLength };
}

/** The view drawn through the fixtures' adapter and through its archive's reader alone: the same SVG, settled. */
async function expectSameDrawing(view: ExampleView) {
  const { spec, archive } = await saved(view);
  const draw = async (s: Node, port: () => AssetsPort) => {
    const ctx = createSeqContext({
      registry: [{ kind: 'format-runtime', runtime: SeqRuntime }, ...(view.registry?.({ delayMs: 0 }) ?? [])],
      ports: { assets: port },
    });
    const r = await ctx.renderSvg(s, { size: { width: 720, height: 320 } });
    ctx.retire();
    return r;
  };
  const fixtures = await draw(spec, () => nodeFixtureAssets());
  const archived = await draw(archive.spec, () => ({ reader: archiveReader(archive) }));
  expect(archived.issues, `${view.dataset}/${view.id}`).toEqual(fixtures.issues);
  expect(archived.settled).toBe(true);
  expect(archived.svg, `${view.dataset}/${view.id}`).toBe(fixtures.svg);
}

const datasets = [...new Set(ExampleViews.map((v) => v.dataset))];

describe('every dataset saved and opened as an archive (EC-91)', () => {
  it.each(datasets)('%s: each view’s files are in the archive, byte for byte', async (dataset) => {
    for (const view of ExampleViews.filter((v) => v.dataset === dataset)) {
      const { spec, assets, archive } = await saved(view);
      const before = urisOf(spec);
      const after = urisOf(archive.spec);
      expect(after).toHaveLength(before.length);
      const reader = archiveReader(archive);
      for (let i = 0; i < before.length; i++) {
        if (generated(before[i])) {
          expect(after[i]).toBe(before[i]);
          continue;
        }
        const original = await assets.readBytes(before[i], signal);
        expect(await reader.readBytes(after[i], signal), `${view.id}: ${before[i]}`).toEqual(original);
      }
    }
  });

  it('draws a view from its archive alone as it is drawn from the fixtures', async () => {
    await expectSameDrawing(ExampleViews.find((v) => v.dataset === 'plasmid')!);
  });

  it('draws every view from its archive as from the fixtures', { tags: ['medium'] }, async () => {
    // a view with a runtime of its own (`generated:`, a profile named in a `value`) reads files its sources do not name,
    // which no archive holds: the writer cannot know them (the host keeps those, `keep`)
    for (const view of ExampleViews) if (!view.registry) await expectSameDrawing(view);
  });
});
