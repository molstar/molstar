/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import type { Coverage, DataEntry, Want } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { Viewport } from './axis.js';
import { DefaultRenderConfig, type RenderConfig } from './config.js';
import { DrawLists, EmphasisPasses } from './draw-lists.js';
import { EmptyOverlay, type Emphasis, type LayerType, type Overlay, type StoreReader } from './layer-type.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import { BuiltInLayerTypes, drawView, LayoutState, layoutView, type ViewLayout } from './tracks.js';

// Frames keyed by what they read (_spec/README.md "Render", The key): the model's hash, the view, the size, the
// viewport, the device pixel ratio, the token table's version, the text metrics' version, the overlay's version, whether the view settled, and the
// versions of the store entries the last layout read. An equal key reuses the layout and the lists and redraws only the
// emphasis passes (a hover); anything else lays out again into the same buffers.

export interface RendererOptions {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly data: StoreReader;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  readonly config?: RenderConfig;
  readonly types?: Readonly<Record<string, LayerType<any, any>>>;
}

export interface FrameInput {
  readonly size: { readonly width: number; readonly height: number };
  readonly viewport?: Viewport;
  readonly dpr?: number;
  readonly settled?: boolean;
  readonly overlay?: Overlay;
  readonly emphasis?: Emphasis;
}

export interface Frame {
  readonly lists: DrawLists;
  readonly layout: ViewLayout;
  /** The draw lists' key. */
  readonly key: string;
  /** Whether this frame laid out again (else only the emphasis passes were drawn). */
  readonly rebuilt: boolean;
  /** What the layout needs from the store. */
  readonly wants: readonly Want[];
}

export class Renderer {
  tokens: SeqTokens;
  readonly lists = new DrawLists();
  private readonly state = new LayoutState();
  private readonly config: RenderConfig;
  private readonly configKey: string;
  private layout: ViewLayout | undefined;
  private baseKey = '';
  private readonly reads = new Map<string, number>();
  private readonly coverages = new Map<string, Coverage | undefined>();

  constructor(readonly options: RendererOptions) {
    this.tokens = options.tokens;
    this.config = options.config ?? DefaultRenderConfig;
    this.configKey = (hashString(JSON.stringify(this.config)) >>> 0).toString(16);
  }

  /** Whether the entries and coverages the last layout read are as they were. */
  private readsHold(): boolean {
    for (const [k, v] of this.reads) if ((this.options.data.get(k)?.version ?? -1) !== v) return false;
    for (const [s, c] of this.coverages) if (this.options.data.coverage(s) !== c) return false;
    return true;
  }

  frame(input: FrameInput): Frame {
    const { model, view } = this.options;
    const dpr = input.dpr ?? 1;
    const vp = input.viewport ? `${input.viewport.start}:${input.viewport.end}` : 'all';
    const overlay = input.overlay ?? EmptyOverlay;
    const base = [
      model.hash,
      view.address,
      `${input.size.width}x${input.size.height}@${dpr}`,
      vp,
      this.tokens.version,
      this.options.text.version,
      overlay.version,
      input.settled ? 's' : 'm',
      this.configKey,
    ].join('|');
    if (this.layout && base === this.baseKey && this.readsHold()) {
      drawView(this.layout, this.lists, input.emphasis, EmphasisPasses);
      return { lists: this.lists, layout: this.layout, key: this.lists.key, rebuilt: false, wants: this.layout.wants };
    }
    this.reads.clear();
    this.coverages.clear();
    const data = this.options.data;
    const reader: StoreReader = {
      coverage: (s) => {
        const c = data.coverage(s);
        this.coverages.set(s, c);
        return c;
      },
      keysOf: (w) => data.keysOf(w),
      get: (k): DataEntry | undefined => {
        const e = data.get(k);
        this.reads.set(k, e?.version ?? -1);
        return e;
      },
    };
    const layout = layoutView(
      {
        model,
        view,
        size: input.size,
        ...(input.viewport ? { viewport: input.viewport } : {}),
        dpr,
        tokens: this.tokens,
        text: this.options.text,
        config: this.config,
        data: reader,
        overlay,
        settled: input.settled ?? false,
        types: this.options.types ?? BuiltInLayerTypes,
      },
      this.state,
    );
    this.layout = layout;
    this.baseKey = base;
    drawView(layout, this.lists, input.emphasis);
    let sig = '';
    for (const [k, v] of [...this.reads].sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0))) sig += `${k}=${v};`;
    this.lists.key = `${base}|${(hashString(sig) >>> 0).toString(16)}`;
    return { lists: this.lists, layout, key: this.lists.key, rebuilt: true, wants: layout.wants };
  }
}
