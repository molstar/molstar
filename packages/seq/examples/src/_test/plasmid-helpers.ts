/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import type { Node } from '@molstar/core/repr/tree';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import {
  attachSources,
  type Frame,
  type FrameInput,
  bundledTextMetrics,
  LightTokens,
  Renderer,
} from '@molstar/seqx-core/render';
import { resolveSpec, sourceLengths } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';
import { PlasmidFiles } from '../plasmid/derive.js';

// What the Plasmid's tests share: the published record's letters (read from the committed file), and a view settled
// against the fixtures.

/** The letters of the committed `seq` file. */
export function seqLetters(file: string): string {
  const data = JSON.parse(readFileSync(`${PackageFolder}${file}`, 'utf8'));
  return data.params.pairs[0].answer.tables.letters.letters as string;
}

export const publishedLetters = () => seqLetters(PlasmidFiles.published);

export interface Settled {
  readonly renderer: Renderer;
  frame(input: FrameInput): Promise<Frame>;
}

/** A spec resolved and its sources opened against the fixtures; `frame` runs until the data has arrived. */
export async function settled(spec: Node): Promise<Settled> {
  const runtimes = [SeqRuntime, RecordsRuntime];
  const assets = nodeFixtureAssets().reader;
  const store = new DataStore();
  // the sources that give a system its length first, then the model resolved again with their lengths
  let model = resolveSpec(spec, { runtimes });
  const sources = new Set([...model.systems.values()].flatMap((s) => (s.source ? [s.source] : [])));
  if (sources.size) await attachSources(model, store, { assets, only: sources });
  const lengths = sourceLengths(model, (s) => store.coverage(s));
  if (lengths.size) model = resolveSpec(spec, { runtimes, lengths });
  const view = resolveViews(model)[0];
  await attachSources(model, store, { assets });
  const renderer = new Renderer({ model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() });
  return {
    renderer,
    async frame(input) {
      let f = renderer.frame(input);
      for (let i = 0; i < 4; i++) {
        store.want(f.wants);
        await store.whenSettled();
        const next = renderer.frame(input);
        if (!next.rebuilt) break;
        f = next;
      }
      return renderer.frame(input);
    },
  };
}
