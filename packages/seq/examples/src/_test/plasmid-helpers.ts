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
import { resolveSpec } from '@molstar/seqx-core/resolve';
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
  const model = resolveSpec(spec, { runtimes: [SeqRuntime, RecordsRuntime] });
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store, { assets: nodeFixtureAssets().reader });
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
