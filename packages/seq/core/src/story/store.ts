/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The story's encoder (_spec/README.md "Stories", Storing): states in, steps out. A full state and a patch over the
// spec are keyframes, a patch `from` the step before a delta; `smaller` picks per step the smallest form within
// `story.maxChain` deltas of a keyframe; a step kept `full` stays a full state, one kept `spec` a patch over it.

import { diff, type Patch, patchToNode } from '@molstar/core/repr/patch';
import type { Node, Value } from '@molstar/core/repr/tree';
import { DefaultStoryConfig, type Easing, type StoryConfig } from './config.js';
import { stripStory } from './plan.js';

/** A step to store: its state, and what is kept beside it. */
export interface StepInput {
  readonly id: string;
  readonly root: Node;
  readonly label?: string;
  readonly keep?: 'full' | 'spec';
  readonly view?: Value;
  readonly duration?: number;
  readonly easing?: Easing;
  readonly linger?: number;
}

export interface StoreOptions {
  readonly store?: StoryConfig['store'];
  readonly maxChain?: number;
}

const sizeOf = (x: unknown) => JSON.stringify(x, (_k, v) => (ArrayBuffer.isView(v) ? Array.from(v as any) : v)).length;

const noBase = (p: Patch): Patch => ({ ops: p.ops });

/**
 * The `snapshot` nodes storing the states over `spec`, in order: each as `story.store` says (a patch over the spec by
 * default). A step's `stateOf` gives its state back (the story's own container aside).
 */
export function storeSteps(spec: Node, steps: readonly StepInput[], options: StoreOptions = {}): Node[] {
  const mode = options.store ?? DefaultStoryConfig.store;
  const maxChain = Math.max(0, options.maxChain ?? DefaultStoryConfig.maxChain);
  const base = stripStory(spec);
  const out: Node[] = [];
  let prev: { id: string; root: Node; chain: number } | undefined;
  for (const s of steps) {
    const root = stripStory(s.root);
    const kept: Record<string, Value> = {};
    for (const k of ['label', 'keep', 'view', 'duration', 'easing', 'linger'] as const)
      if (s[k] !== undefined) kept[k] = s[k] as Value;
    const overSpec = () => ({ patch: patchToNode(noBase(diff(base, root))) as Value });
    let content: Record<string, Value>;
    let chain = 0;
    const form = s.keep === 'full' ? 'full' : s.keep === 'spec' ? 'patch' : mode;
    if (form === 'full') content = { state: root as Value };
    else if (form === 'patch') content = overSpec();
    else {
      const candidates: { content: Record<string, Value>; chain: number }[] = [
        { content: overSpec(), chain: 0 },
        { content: { state: root as Value }, chain: 0 },
      ];
      if (prev && prev.chain < maxChain)
        candidates.push({
          content: { from: `/snapshots/${prev.id}`, patch: patchToNode(noBase(diff(prev.root, root))) as Value },
          chain: prev.chain + 1,
        });
      let best = candidates[0];
      let bestSize = sizeOf(best.content);
      for (const c of candidates.slice(1)) {
        const size = sizeOf(c.content);
        if (size < bestSize) {
          best = c;
          bestSize = size;
        }
      }
      content = best.content;
      chain = best.chain;
    }
    out.push({ kind: 'snapshot', id: s.id, params: { ...kept, ...content } });
    prev = { id: s.id, root, chain };
  }
  return out;
}

/** The spec with its story replaced by `steps` (the `snapshots` container last, as the builder writes it). */
export function withSteps(spec: Node, steps: readonly Node[]): Node {
  const rest = stripStory(spec).children ?? [];
  return {
    ...spec,
    children: steps.length ? [...rest, { kind: 'snapshots', id: 'snapshots', children: steps }] : rest,
  };
}
