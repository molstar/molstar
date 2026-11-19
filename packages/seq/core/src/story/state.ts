/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A step's state (_spec/README.md "Stories", A step's state): its full state, or the spec (or the step it starts
// `from`) with its patch applied; computed on demand, cached by the spec's revision and the step. A patch that no
// longer applies whole gives the ops that still apply, flagged.

import { applyPatch, type ApplyOptions, type Patch, patchFromNode } from '@molstar/core/repr/patch';
import { hashNode, type Issue, issueOf, type Node, printAddress, resolveAddress } from '@molstar/core/repr/tree';
import { LRUCache } from '@molstar/core/util/lru-cache';
import { SeqLanguage, SeqValidateOptions } from '../lang/language.js';
import { DefaultStoryConfig, easingOf, type StepTiming, type StoryConfig } from './config.js';
import { StoryIssues } from './issues.js';

/** A section's viewport in a step, as positions (never offsets: offsets change meaning when the axis does). */
export interface StepViewport {
  readonly start: number;
  readonly end: number;
  readonly top?: number;
  readonly scale?: number;
  readonly bands?: Readonly<Record<string, { readonly bandTop?: number; readonly rowScale?: number }>>;
  readonly y?: { readonly start: number; readonly end: number };
}

/** What a step shows: the view (its address; the first when absent) and each section's viewport by address. */
export interface StepView {
  readonly view?: string;
  readonly sections: Readonly<Record<string, StepViewport>>;
}

export interface StepState {
  readonly id: string;
  /** `/snapshots/<id>`. */
  readonly address: string;
  readonly label: string;
  /** How it is stored: a patch (over the spec, or a delta `from` a step) or a full state. */
  readonly form: 'patch' | 'full';
  /** The step it starts from (a delta's). */
  readonly from?: string;
  /** The state: a spec tree. */
  readonly root: Node;
  readonly view?: StepView;
  readonly timing: StepTiming;
  /** False when some of its ops (or a start's) no longer apply: the state is the ops that do. */
  readonly applies: boolean;
  /** Its own ops that did not apply. */
  readonly skipped: number;
  /** Deltas between it and its keyframe (0: a keyframe). */
  readonly chain: number;
  readonly issues: readonly Issue[];
}

/** Where `stateOf` keeps what it computed: by the spec's revision and the step. */
export interface StoryCache {
  readonly lru: LRUCache<StepState>;
}

export function createStoryCache(capacity = 64): StoryCache {
  return { lru: LRUCache.create(capacity) };
}

export interface StateOfOptions {
  /** How a patch applies (the language, its slots, refs followed, dangling allowed: a step is never refused). */
  readonly apply?: ApplyOptions;
  readonly cache?: StoryCache;
  readonly config?: StoryConfig;
}

/** A step applies with the language and its slots, refs followed, nothing else refused. */
export const StepApplyOptions: ApplyOptions = {
  language: SeqLanguage,
  validate: SeqValidateOptions,
  allowDangling: true,
};

/** The ids of the spec's steps, in order. */
export function storySteps(spec: Node): readonly string[] {
  const c = resolveAddress(spec, ['snapshots']);
  const out: string[] = [];
  for (const n of c?.children ?? []) if (n.kind === 'snapshot' && n.id !== undefined) out.push(n.id);
  return out;
}

/** A stored patch node read once. */
const Patches = new WeakMap<Node, Patch | undefined>();

function patchOf(node: Node): Patch | undefined {
  const p = node.params?.patch;
  if (!p || typeof p !== 'object' || Array.isArray(p) || !('kind' in p)) return undefined;
  if (!Patches.has(p as Node)) Patches.set(p as Node, patchFromNode(p as Node).patch);
  return Patches.get(p as Node);
}

const num = (x: unknown) => (typeof x === 'number' && Number.isFinite(x) ? x : undefined);

/** A step's `view` param read: malformed fields are left out. */
export function readStepView(x: unknown): StepView | undefined {
  if (!x || typeof x !== 'object' || Array.isArray(x)) return undefined;
  const v = x as Record<string, unknown>;
  const sections: Record<string, StepViewport> = {};
  const ss = v.sections;
  if (ss && typeof ss === 'object' && !Array.isArray(ss))
    for (const [address, s] of Object.entries(ss as Record<string, any>)) {
      const start = num(s?.start);
      const end = num(s?.end);
      if (start === undefined || end === undefined) continue;
      const top = num(s.top);
      const scale = num(s.scale);
      const y = s.y && num(s.y.start) !== undefined && num(s.y.end) !== undefined ? s.y : undefined;
      let bands: Record<string, { bandTop?: number; rowScale?: number }> | undefined;
      if (s.bands && typeof s.bands === 'object')
        for (const [t, b] of Object.entries(s.bands as Record<string, any>)) {
          const bandTop = num(b?.bandTop);
          const rowScale = num(b?.rowScale);
          if (bandTop === undefined && rowScale === undefined) continue;
          (bands ??= {})[t] = {
            ...(bandTop !== undefined ? { bandTop } : {}),
            ...(rowScale !== undefined ? { rowScale } : {}),
          };
        }
      sections[address] = {
        start,
        end,
        ...(top !== undefined ? { top } : {}),
        ...(scale !== undefined ? { scale } : {}),
        ...(bands ? { bands } : {}),
        ...(y ? { y: { start: y.start, end: y.end } } : {}),
      };
    }
  return { ...(typeof v.view === 'string' ? { view: v.view } : {}), sections };
}

/** A step's timing: its own, else the configuration's. */
export function timingOf(node: Node, config: StoryConfig = DefaultStoryConfig): StepTiming {
  return {
    duration: num(node.params?.duration) ?? config.transition.duration,
    easing: easingOf(node.params?.easing, config.transition.easing),
    linger: num(node.params?.linger) ?? config.linger,
  };
}

/** Applies a patch; when it fails whole, op by op, keeping those that apply. */
export function applyLoose(root: Node, patch: Patch, options: ApplyOptions): { root: Node; skipped: number } {
  const r = applyPatch(root, patch, options);
  if (r.ok) return { root: r.root, skipped: 0 };
  let skipped = 0;
  for (const op of patch.ops) {
    const one = applyPatch(root, { ops: [op] }, options);
    if (one.ok) root = one.root;
    else skipped++;
  }
  return { root, skipped };
}

const idOfAddress = (a: unknown) => (typeof a === 'string' ? /^\/snapshots\/([^/]+)$/.exec(a)?.[1] : undefined);

/**
 * A step's state over `spec` (the spec shown: its root with a preview): undefined when the spec has no such step.
 * Cached by the spec's revision and the step when a cache is given.
 */
export function stateOf(spec: Node, id: string, options: StateOfOptions = {}): StepState | undefined {
  return compute(spec, id, options, new Set());
}

function compute(spec: Node, id: string, options: StateOfOptions, visiting: Set<string>): StepState | undefined {
  const cache = options.cache?.lru;
  const key = cache ? `${hashNode(spec)}\u0000${id}` : '';
  if (cache) {
    const hit = LRUCache.get(cache, key);
    if (hit) return hit;
  }
  const node = resolveAddress(spec, ['snapshots', id]);
  if (node?.kind !== 'snapshot') return undefined;
  const address = printAddress(['snapshots', id]);
  const config = options.config ?? DefaultStoryConfig;
  const label = typeof node.params?.label === 'string' ? node.params.label : id;
  const view = readStepView(node.params?.view);
  const base = { id, address, label, timing: timingOf(node, config), ...(view ? { view } : {}) };
  const issues: Issue[] = [];
  let out: StepState;
  const full = node.params?.state;
  if (full && typeof full === 'object' && !Array.isArray(full) && 'kind' in full) {
    out = { ...base, form: 'full', root: full as Node, applies: true, skipped: 0, chain: 0, issues };
  } else {
    const patch = patchOf(node);
    let start = spec;
    let chain = 0;
    let applies = true;
    const fromAddress = node.params?.from;
    const fromId = idOfAddress(fromAddress);
    if (fromAddress !== undefined) {
      visiting.add(id);
      const prior =
        fromId !== undefined && !visiting.has(fromId) ? compute(spec, fromId, options, visiting) : undefined;
      visiting.delete(id);
      // a start that cannot be reached, here or further back, is this step's issue too
      const reached = prior && !prior.issues.some((i) => i.code === 'seq.story.step.from');
      if (prior) {
        start = prior.root;
        chain = prior.chain + 1;
        applies = prior.applies;
      }
      if (!reached) {
        applies = false;
        applies = false;
        issues.push(
          issueOf(StoryIssues, 'seq.story.step.from', `The step ${label} starts from no step it can reach`, {
            address,
            hint: StoryIssues.codes['seq.story.step.from'].hint,
          }),
        );
      }
    }
    if (!patch) {
      issues.push(
        issueOf(StoryIssues, 'seq.story.step.form', `The step ${label} holds no patch nor state`, {
          address,
          hint: StoryIssues.codes['seq.story.step.form'].hint,
        }),
      );
      out = { ...base, form: 'patch', root: start, applies: false, skipped: 0, chain, issues };
    } else {
      const r = applyLoose(start, patch, options.apply ?? StepApplyOptions);
      if (r.skipped)
        issues.push(
          issueOf(StoryIssues, 'seq.story.step.stale', `${r.skipped} of the step ${label}’s ops no longer apply`, {
            address,
            hint: StoryIssues.codes['seq.story.step.stale'].hint,
          }),
        );
      out = {
        ...base,
        form: 'patch',
        ...(fromId !== undefined ? { from: fromId } : {}),
        root: r.root,
        applies: applies && r.skipped === 0,
        skipped: r.skipped,
        chain,
        issues,
      };
    }
  }
  if (cache) LRUCache.set(cache, key, out);
  return out;
}
