/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hasErrors, type Issue, issueOf, type Node } from '@molstar/core/repr/tree';
import {
  type AssetReader,
  DataStore,
  type FormatRuntime,
  type OpenedSource,
  SourceOpenError,
  type Want,
} from '@molstar/seqx-core/data';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { type ResolveConfig, resolveSpec, type StaticModel } from '@molstar/seqx-core/resolve';
import { type ResolvedView, resolveView, resolveViews } from '@molstar/seqx-core/view';
import type { Viewport } from './axis.js';
import { type RenderConfig, renderConfig } from './config.js';
import { DrawLists } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { layerSource } from './layer-data.js';
import type { Emphasis } from './layer-type.js';
import { Renderer } from './renderer.js';
import { drawListsToSvg } from './svg.js';
import { bundledTextMetrics, type TextMetrics } from './text.js';
import { builtInTokens, type SeqTokens } from './tokens.js';

// Headless rendering (_spec/README.md "Render", Export): a spec to an SVG string, asynchronously. It resolves the spec,
// opens the sources the view reads, wants what the export's size and viewport need through a want set of its own
// (never the live view's, EC-125), waits for those entries to settle or a timeout (then the export says it is
// unsettled), and draws once more with the view settled (open domains computed, lanes packed).

/** The store's sources a view reads: its placed layers' data, and its sections' systems' sources (for a length). */
export function viewSources(model: StaticModel, view: ResolvedView): Set<string> {
  const out = new Set<string>();
  for (const s of view.sections) {
    const source = model.systems.get(s.axis.system)?.source;
    if (source) out.add(source);
    for (const t of s.tracks)
      for (const l of t.layers) {
        const src = layerSource(model, l);
        if (src) out.add(src.source);
      }
  }
  return out;
}

/** An asset reader that reads nothing: a headless render reads no file or URL unless given a reader (security). */
export const NoAssets: AssetReader = {
  readText: (uri) => Promise.reject(new Error(`This render reads no assets: ${uri}`)),
};

export interface AttachOptions {
  readonly assets?: AssetReader;
  readonly signal?: AbortSignal;
  /** The store's sources to open (all the model names when absent). */
  readonly only?: ReadonlySet<string>;
}

/**
 * Opens the model's sources into the store: inline data through its `seq` runtime, a source node through the runtime
 * bound to it. One already attached is kept. What cannot be opened is an issue, and its layers stay pending.
 */
export async function attachSources(
  model: StaticModel,
  store: DataStore,
  options: AttachOptions = {},
): Promise<Issue[]> {
  const ctx = { assets: options.assets ?? NoAssets, signal: options.signal ?? new AbortController().signal };
  const opens = new Map<string, () => Promise<OpenedSource>>();
  for (const plan of [
    ...[...model.systems.values()].map((s) => s.access),
    ...[...model.annotations.values()].map((a) => a.access),
  ])
    if (plan.kind === 'inline') opens.set(plan.source.address, () => plan.runtime.open(plan.source, ctx));
  for (const b of model.sources.values()) {
    const runtime = b.runtime;
    if (runtime) opens.set(b.address, () => runtime.open(b.entry, ctx));
  }
  const issues: Issue[] = [];
  await Promise.all(
    [...opens].map(async ([address, open]) => {
      if ((options.only && !options.only.has(address)) || store.coverage(address)) return;
      try {
        store.attach(address, await open());
      } catch (e) {
        if (e instanceof SourceOpenError) issues.push(...e.issues.map((i) => ({ ...i, address })));
        else
          issues.push(
            issueOf(
              RenderIssues,
              'seq.render.source',
              `The source could not be opened: ${String((e as Error)?.message ?? e)}`,
              {
                address,
                hint: RenderIssues.codes['seq.render.source'].hint,
              },
            ),
          );
      }
    }),
  );
  return issues;
}

export interface RenderSvgOptions {
  /** CSS px. */
  readonly size: { readonly width: number; readonly height: number };
  readonly viewport?: Viewport;
  /** The format runtimes the application composed (inline data always reads through `seq`). */
  readonly runtimes?: readonly FormatRuntime[];
  /** The view's address; the spec's first view when absent. */
  readonly view?: string;
  /** The built-in token table of a theme (`light` by default), or a table. */
  readonly theme?: 'light' | 'dark';
  readonly tokens?: SeqTokens;
  /** The text metrics port: the bundled table by default. */
  readonly text?: TextMetrics;
  readonly config?: Partial<RenderConfig>;
  readonly resolve?: Partial<ResolveConfig>;
  /** What sources may read: nothing by default. */
  readonly assets?: AssetReader;
  /** How long the export waits for its data, ms: 10 s by default. */
  readonly timeoutMs?: number;
  readonly dpr?: number;
  readonly emphasis?: Emphasis;
  /** A store to read through (a viewer's): the export wants through a want set of its own. A new one by default. */
  readonly store?: DataStore;
}

export interface RenderSvgResult {
  readonly svg: string;
  /** Whether everything the export needed arrived before the timeout. */
  readonly settled: boolean;
  /** The spec's, the view's, the sources' and the render's issues. */
  readonly issues: readonly Issue[];
  readonly lists: DrawLists;
}

function wantsKey(wants: readonly Want[]): string {
  return wants.map((w) => `${w.source}|${w.system}|${w.start}|${w.end}|${w.level}|${w.priority}`).join(';');
}

/** A spec drawn to an SVG string once the data its view needs has settled (or the timeout passed). */
export async function renderSvg(spec: Node, options: RenderSvgOptions): Promise<RenderSvgResult> {
  const issues: Issue[] = [...checkSeqSpec(spec)];
  const blank = (settled: boolean): RenderSvgResult => {
    const lists = new DrawLists();
    lists.width = options.size.width;
    lists.height = options.size.height;
    return { svg: drawListsToSvg(lists), settled, issues, lists };
  };
  if (hasErrors(issues)) return blank(true);
  const model = resolveSpec(spec, {
    ...(options.runtimes ? { runtimes: options.runtimes } : {}),
    config: options.resolve ?? {},
  });
  issues.push(...model.issues);
  const view = options.view ? resolveView(model, options.view) : resolveViews(model)[0];
  if (!view) {
    issues.push(
      issueOf(RenderIssues, 'seq.render.view', 'The spec has no such view', {
        ...(options.view ? { address: options.view } : {}),
        hint: RenderIssues.codes['seq.render.view'].hint,
      }),
    );
    return blank(true);
  }
  issues.push(...view.issues);
  const own = !options.store;
  const store = options.store ?? new DataStore();
  const controller = new AbortController();
  const wanter = store.wanter();
  try {
    issues.push(
      ...(await attachSources(model, store, {
        ...(options.assets ? { assets: options.assets } : {}),
        signal: controller.signal,
        only: viewSources(model, view),
      })),
    );
    const tokens = options.tokens ?? builtInTokens(options.theme ?? 'light');
    const config = renderConfig(options.config);
    const renderer = new Renderer({
      model,
      view,
      data: store,
      tokens,
      text: options.text ?? bundledTextMetrics(),
      config,
    });
    const frameInput = {
      size: options.size,
      ...(options.viewport ? { viewport: options.viewport } : {}),
      dpr: options.dpr ?? 1,
      ...(options.emphasis ? { emphasis: options.emphasis } : {}),
    };
    const timeout = options.timeoutMs ?? 10_000;
    const until = Date.now() + timeout;
    let settled = true;
    let last = '';
    // what a layout wants can follow from what arrived (a length, a follower's range): ask until it holds
    for (let round = 0; round < 8; round++) {
      const f = renderer.frame({ ...frameInput, settled: false });
      const key = wantsKey(f.wants);
      if (key === last) break;
      last = key;
      wanter.want(f.wants);
      if (!(await wanter.whenSettled(Math.max(0, until - Date.now())))) {
        settled = false;
        break;
      }
    }
    const frame = renderer.frame({ ...frameInput, settled: true });
    issues.push(...frame.layout.issues);
    if (!settled)
      issues.push(
        issueOf(RenderIssues, 'seq.render.unsettled', 'The export did not settle before its timeout', {
          hint: RenderIssues.codes['seq.render.unsettled'].hint,
        }),
      );
    issues.push(...store.issues.value);
    const domains: string[] = [];
    for (const s of frame.layout.sections)
      for (const t of s.tracks)
        for (const l of t.layers) {
          const legend = l.type.legend?.(l.lay);
          if (legend?.domain) domains.push(`${l.key}=${legend.domain[0]},${legend.domain[1]}`);
        }
    const svg = drawListsToSvg(frame.lists, {
      fonts: tokens.fonts,
      fade: config.fade,
      metadata: { 'seq-settled': String(settled), ...(domains.length ? { 'seq-domains': domains.join(';') } : {}) },
    });
    return { svg, settled, issues, lists: frame.lists };
  } finally {
    wanter.dispose();
    controller.abort();
    if (own) store.retire();
  }
}
