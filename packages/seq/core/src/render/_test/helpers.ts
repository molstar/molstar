/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Test helpers: a settled frame of a spec (resolved, its sources attached, its wants answered), and a pass's
// primitives counted per track (each track's data is one clip block).

import type { Node } from '@molstar/core/repr/tree';
import { DataStore, type FormatRuntime } from '@molstar/seqx-core/data';
import { resolveSpec, type StaticModel, sourceLengths } from '@molstar/seqx-core/resolve';
import { type ResolvedView, resolveViews } from '@molstar/seqx-core/view';
import {
  attachSources,
  bundledTextMetrics,
  type DrawLists,
  type Frame,
  type FrameInput,
  LightTokens,
  type PassName,
  Renderer,
  type RenderConfigOverrides,
  renderConfig,
  walkPass,
} from '../index.js';

export interface Settled {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly store: DataStore;
  readonly renderer: Renderer;
  /** A frame once its wants have settled. */
  frame(input: FrameInput): Promise<Frame>;
}

export async function settle(
  spec: Node,
  runtimes: readonly FormatRuntime[] = [],
  config: RenderConfigOverrides = {},
): Promise<Settled> {
  const store = new DataStore();
  // the sources that give a system its length first, then the model resolved again with their lengths
  let model = resolveSpec(spec, { runtimes });
  const sources = new Set([...model.systems.values()].flatMap((s) => (s.source ? [s.source] : [])));
  if (sources.size) await attachSources(model, store, { only: sources });
  const lengths = sourceLengths(model, (s) => store.coverage(s));
  if (lengths.size) model = resolveSpec(spec, { runtimes, lengths });
  const view = resolveViews(model)[0];
  await attachSources(model, store);
  const renderer = new Renderer({
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    config: renderConfig(config),
  });
  return {
    model,
    view,
    store,
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

export interface BlockCount {
  rects: number;
  letters: number;
  texts: number;
  images: number;
  paths: number;
}

/** A pass's primitives per clip block (a track's data), in order. */
export function blocks(lists: DrawLists, pass: PassName = 'data'): BlockCount[] {
  const out: BlockCount[] = [];
  let cur: BlockCount | undefined;
  const add = (k: keyof BlockCount, n = 1) => {
    if (cur) cur[k] += n;
  };
  walkPass(lists.passes[pass], {
    rect: () => add('rects'),
    path: () => add('paths'),
    glyphs: (l, i) => add('letters', l.count[i]),
    text: () => add('texts'),
    image: () => add('images'),
    clip: () => {
      cur = { rects: 0, letters: 0, texts: 0, images: 0, paths: 0 };
      out.push(cur);
    },
    unclip: () => {
      cur = undefined;
    },
  });
  return out;
}
