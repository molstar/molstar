/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Issue, Node } from '@molstar/core/repr/tree';
import { createSeqContext, type SeqContext } from '@molstar/seqx-core/context/index';
import {
  type Emphasis,
  type Frame,
  type FrameInput,
  type Hit as RenderHit,
  Renderer,
  type RenderSvgOptions,
  type RenderSvgResult,
  type Viewport as RenderViewport,
  hitView,
} from '@molstar/seqx-core/render';
import { type ResolvedView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DomHost } from '@molstar/ui/dom';
import { bindingTable, DefaultActions, type BindingTable, BindingPresets } from '@molstar/seqx-viewer/bindings';
import type { CanvasBackend } from '@molstar/seqx-viewer/canvas/index';
import { type ViewerConfig, viewerConfigOf } from '@molstar/seqx-viewer/config';
import { HoverModel, SpecModel, StoreModel, ViewModel } from '@molstar/seqx-viewer/models';
import { ViewerInput } from './input.js';
import type { BehaviorContext, SeqViewerSpec } from './viewer-spec.js';
import { mountCanvas } from './mount.js';

// The viewer (_spec/README.md "Composition"): composed from a spec by `createSeqViewer`, as a Mol* plugin is: the
// context (registry, configuration, ports) first, then the models, the bindings and the behaviors, then the UI. Frames
// and time are `DomHost`'s. The frame: the renderer lays the view out (or, for a hover, redraws the emphasis passes
// alone) and the Canvas2D backend draws what changed.

const NeverAborted = new AbortController().signal;

type Mutable<T> = { -readonly [K in keyof T]: T[K] };

/** A command a behavior provides, with the context it runs in. */
export type ViewerCommand = CommandDefinition<any>;

export class SeqViewer extends ReactiveModel {
  readonly context: SeqContext;
  readonly config: ViewerConfig;
  readonly spec: SpecModel;
  readonly view: ViewModel;
  readonly hover: HoverModel;
  readonly store: StoreModel;
  readonly bindings: BindingTable;
  readonly input: ViewerInput;
  /** The behaviors' handles by provider name. */
  readonly behaviors = new Map<string, any>();
  /** The commands the behaviors provide. */
  readonly commands: ViewerCommand[] = [];
  /** The layout's issues of the last frame (a color it cannot read, a layer it cannot draw). */
  readonly layoutIssues: ReadonlyState<readonly Issue[]>;

  /** The canvas's size in CSS px (the data's width and the content's height by default). */
  width = 0;
  height = 0;
  /** The device pixel ratio of the last frame. */
  dpr = 1;
  /** Frames drawn and layouts made: what a test or a profile reads. */
  readonly counts = { frames: 0, layouts: 0 };

  private readonly _issues = new State<readonly Issue[]>([]);
  private readonly _emphasis: Mutable<Emphasis> = {};
  private readonly _emphasisHit = {} as Mutable<RenderHit>;
  private readonly _size = { width: 0, height: 0 };
  private readonly _frameInput = { size: this._size } as Mutable<FrameInput>;
  private readonly _viewportOut = { start: 0, end: 0 };
  private _backend: CanvasBackend | undefined;
  private _renderer: Renderer | undefined;
  private _rendererFor: unknown[] = [];
  private _wanter: { want(w: readonly any[]): void; dispose(): void } | undefined;
  private _frame: Frame | undefined;
  private _scheduled = false;
  private _quiet = true;
  private _settleTimer: unknown;
  private _syncing: Promise<void> = Promise.resolve();
  private _contentHeight = 0;
  private _baseDirty = true;
  private _hoverVersion = -1;
  private _retired = false;

  constructor(
    viewerSpec: SeqViewerSpec,
    target?: HTMLElement,
    private readonly fixedDpr?: number,
  ) {
    super();
    this.config = viewerConfig(viewerSpec);
    // the headless pieces first: the registry, the configuration and the ports
    this.context = createSeqContext({
      ...(viewerSpec.registry ? { registry: viewerSpec.registry } : {}),
      ...(viewerSpec.config ? { config: viewerSpec.config } : {}),
      ...(viewerSpec.ports ? { ports: viewerSpec.ports } : {}),
    });
    this.registerPermanentDispose(() => this.context.retire());
    this.layoutIssues = this._issues.asReadonly;

    this.spec = new SpecModel(createSeqSpec().build(), {
      runtimes: this.context.runtimes,
      config: this.context.resolveConfig,
    });
    this.view = new ViewModel(this.spec, this.config);
    this.hover = new HoverModel(this.spec, this.view, this.config);
    this.store = new StoreModel(this.config);
    for (const m of [this.store, this.hover, this.view, this.spec]) this.registerPermanentDispose(() => m.retire());

    const preset = viewerSpec.bindings?.preset ?? BindingPresets[this.config.bindings.preset];
    const extraActions = viewerSpec.bindings?.actions ?? [];
    this.bindings = bindingTable(preset, viewerSpec.bindings?.overrides ?? [], {
      actions: [...DefaultActions, ...extraActions],
    });

    // behaviors, after the models; one of each name
    const context: BehaviorContext = {
      spec: this.spec,
      view: this.view,
      hover: this.hover,
      store: this.store,
      config: this.config,
      bindings: this.bindings,
    };
    for (const b of viewerSpec.behaviors) {
      if (this.behaviors.has(b.provider.name)) throw new Error(`The behavior ${b.provider.name} is listed twice`);
      const handle = b.provider.create(context, b.params);
      this.behaviors.set(b.provider.name, handle);
      if (Array.isArray(handle?.commands)) this.commands.push(...handle.commands);
      if (typeof handle?.dispose === 'function') this.registerPermanentDispose(() => handle.dispose());
    }
    this.input = new ViewerInput(this);

    // what moves the picture asks for a frame
    this.listenPermanently(this.spec.model, () => this._sync());
    this.listenPermanently(this.hover.version, () => this.schedule());
    this.listenPermanently(this.store.store.version, () => this.schedule());
    this.listenPermanently(this.store.store.settled, () => this.schedule());
    this.listenPermanently(this.view.resolved, () => {
      this._followViewports();
      this._applyWidths();
      this.schedule();
    });
    this._followViewports();
    this.registerPermanentDispose(() => {
      this._retired = true;
      DomHost.clearTimeout(this._settleTimer);
      this._wanter?.dispose();
      this._viewportListeners.forEach((u) => u());
      this._viewportListeners.length = 0;
    });

    if (target && viewerSpec.ui?.canvas) this.registerPermanentDispose(mountCanvas(this, target));
  }

  /** Whether `retire` has run: a frame asked for before it does nothing. */
  get retired(): boolean {
    return this._retired;
  }

  // ---- the document ----

  /**
   * Shows a spec: sets the root and opens its sources. Resolves when they are open (their data arrives after: `settle`).
   * `view`: the view's address (the first by default).
   */
  async load(root: Node, options: { view?: string } = {}): Promise<void> {
    this.spec.apply(root);
    this.view.show(options.view);
    await this._sync();
  }

  private _sync(): Promise<void> {
    if (this._retired) return this._syncing;
    this._applyWidths();
    this._syncing = this._syncing.then(() =>
      this._retired ? undefined : this.store.sync(this.spec.model.value, this.context.assets),
    );
    this._syncing = this._syncing.then(() => this.schedule());
    return this._syncing;
  }

  /**
   * Resolves once the picture is complete: the sources open, the data the view wants arrived, the view settled. A
   * viewer with no size yet is given the default one. For tests, exports and stories that wait for the picture.
   */
  async settle(timeoutMs = 10_000): Promise<boolean> {
    await this._syncing;
    const until = DomHost.now() + timeoutMs;
    this._quiet = true;
    for (let round = 0; round < 8; round++) {
      const frame = this.draw();
      if (!frame) return true;
      const left = Math.max(0, until - DomHost.now());
      if (!(await this.store.store.whenSettled(left))) return false;
      if (!this.draw()?.rebuilt && !frame.rebuilt) break;
    }
    this.draw();
    return true;
  }

  // ---- size and frames ----

  /** Sets the canvas width (CSS px) and, for a fixed height, the height; the sections' widths follow. */
  resize(width: number, height?: number): void {
    if (this._retired) return;
    this.width = Math.max(0, Math.round(width));
    if (height !== undefined) this.height = Math.max(0, Math.round(height));
    this._applyWidths();
    this.schedule();
  }

  /** Attaches the Canvas2D backend (mounted UI), or none: the viewer then lays out and hit-tests without drawing. */
  useBackend(backend: CanvasBackend | undefined): void {
    this._backend = backend;
    this._baseDirty = true;
    this.schedule();
  }

  /** Drops what was measured: a web font arrived. */
  invalidateText(): void {
    this.context.text.clear();
    this.schedule();
  }

  /** The data box's left edge: the page's inset, the label column and its gap. */
  get dataX(): number {
    const c = this.context.renderConfig;
    return c.inset + c.labelWidth + c.labelGap;
  }

  private _applyWidths(): void {
    const sections = this.view.resolved.value?.sections;
    if (!sections || this.width <= 0) return;
    const dataW = Math.max(1, this.width - this.dataX - this.context.renderConfig.inset);
    for (const s of sections) this.view.setWidth(s.address, dataW);
  }

  private _viewportListeners: (() => void)[] = [];

  /** A frame is asked for by each section's viewport changing: listened to, again when the view's sections change. */
  private _followViewports(): void {
    for (const u of this._viewportListeners) u();
    this._viewportListeners.length = 0;
    for (const s of this.view.resolved.value?.sections ?? [])
      this._viewportListeners.push(
        this.view.viewport(s.address).listen(() => {
          this._moved();
          this.schedule();
        }),
      );
  }

  /** The view moved: it is not settled until it rests for `canvas.settleMs`. */
  private _moved(): void {
    if (this._retired) return;
    this._quiet = false;
    DomHost.clearTimeout(this._settleTimer);
    this._settleTimer = DomHost.setTimeout(() => {
      this._settleTimer = undefined;
      this._quiet = true;
      this.schedule();
    }, this.config.canvas.settleMs);
  }

  /** Asks for a frame: one a `requestAnimationFrame`, whatever asked. */
  schedule(): void {
    if (this._scheduled || this._retired) return;
    this._scheduled = true;
    DomHost.requestAnimationFrame(() => {
      this._scheduled = false;
      if (!this._retired) this.draw();
    });
  }

  /** The window the renderer draws: the first section that leads (the renderer takes one window for all). */
  private _viewport(resolved: ResolvedView): RenderViewport | undefined {
    for (const s of resolved.sections) {
      if (this.view.isFollowing(s.address)) continue;
      const v = this.view.get(s.address);
      if (!v || v.pending) return undefined;
      this._viewportOut.start = v.start;
      this._viewportOut.end = v.end;
      return this._viewportOut;
    }
    return undefined;
  }

  /** The viewport the last frame was drawn with, a copy (for the SVG and tests). */
  get viewport(): RenderViewport | undefined {
    const resolved = this.view.resolved.value;
    const v = resolved ? this._viewport(resolved) : undefined;
    return v ? { start: v.start, end: v.end } : undefined;
  }

  /**
   * Draws a frame now: lays the view out when anything it reads changed (else only the emphasis passes are redrawn, a
   * hover), asks the store for what the layout wants, and draws what changed. Undefined when there is nothing to draw
   * (no view, no width).
   */
  draw(): Frame | undefined {
    if (this._retired) return undefined;
    const model = this.spec.model.value;
    const resolved = this.view.resolved.value;
    if (!resolved || this.width <= 0) return undefined;
    const rr = this._rendererFor;
    if (!this._renderer || rr[0] !== model || rr[1] !== resolved) {
      this._renderer = new Renderer({
        model,
        view: resolved,
        data: this.store.store,
        tokens: this.context.tokens,
        text: this.context.text,
        config: this.context.renderConfig,
      });
      rr[0] = model;
      rr[1] = resolved;
      this._wanter?.dispose();
      this._wanter = this.store.store.wanter();
    }
    const renderer = this._renderer;
    const dpr =
      this.fixedDpr ??
      (this.config.canvas.dpr === 'device' ? (globalThis.devicePixelRatio ?? 1) : this.config.canvas.dpr);
    this.dpr = dpr;

    const { height: h } = this.config.canvas;
    const fixed = h !== 'content';
    const size = this._size;
    size.width = this.width;
    size.height = fixed ? h : Math.min(this.config.canvas.maxHeight, Math.max(this._contentHeight, 32));
    this._updateEmphasis();
    const input = this._frameInput;
    const vp = this._viewport(resolved);
    if (vp) input.viewport = vp;
    else delete input.viewport;
    input.dpr = dpr;
    input.settled = this._quiet && this.store.store.settled.value;
    input.emphasis = this._emphasis;

    let frame = renderer.frame(input);
    if (!fixed) {
      const last = frame.layout.sections.at(-1);
      const content = last ? Math.ceil(last.y + last.height + this.context.renderConfig.inset) : 32;
      if (content !== this._contentHeight) {
        this._contentHeight = content;
        size.height = Math.min(this.config.canvas.maxHeight, Math.max(content, 32));
        frame = renderer.frame(input);
      }
    }
    this.height = size.height;
    if (frame.rebuilt) {
      this.counts.layouts++;
      this._wanter?.want(frame.wants);
      this._issues.set(frame.layout.issues);
    }
    this._frame = frame;

    // a resting pointer is hit again once a frame the view moves under it (EC-143)
    if (frame.rebuilt && this.config.pointer.rehit && this.input.rehit()) {
      this._updateEmphasis();
      frame = renderer.frame(input);
      this._frame = frame;
    }

    const backend = this._backend;
    if (backend) {
      this._baseDirty = backend.resize(size.width, size.height, dpr) || this._baseDirty;
      backend.draw(frame.lists, frame.rebuilt || this._baseDirty);
      this._baseDirty = false;
    }
    this.counts.frames++;
    return frame;
  }

  /** The last frame: its layout is what hit-testing reads. */
  get frame(): Frame | undefined {
    return this._frame;
  }

  /** The hover as the emphasis passes draw it: one reused record, mutated when the hover changes. */
  private _updateEmphasis(): void {
    const v = this.hover.version.value;
    if (v === this._hoverVersion) return;
    this._hoverVersion = v;
    const r = this.hover.current;
    if (r.owner === undefined || r.layer === undefined || r.system === undefined) {
      this._emphasis.hover = undefined;
      return;
    }
    const h = this._emphasisHit;
    h.layer = r.layer;
    h.section = r.section === undefined ? '' : r.section.slice(r.section.lastIndexOf('/') + 1);
    h.system = r.system;
    (h as { position: number | undefined }).position = Number.isNaN(r.position) ? undefined : r.position;
    (h as { element: unknown }).element = r.element;
    (h as { extent: unknown }).extent = r.extent;
    this._emphasis.hover = h;
  }

  /** What is under a point of the last frame's layout. */
  hit(x: number, y: number): RenderHit | undefined {
    return this._frame ? hitView(this._frame.layout, x, y) : undefined;
  }

  // ---- commands ----

  /** Runs a command a behavior provides; the message of its result, or why it could not run. */
  runCommand(id: string, context: unknown): string {
    const c = this.commands.find((x) => x.id === id);
    if (!c) return `No command ${id}`;
    const can = c.when?.(context);
    if (can !== undefined && can !== true) return typeof can === 'string' ? can : 'Not available';
    const r = c.run({}, context as never, NeverAborted) as CommandResult;
    return r.message ?? '';
  }

  // ---- export ----

  /**
   * The picture as an SVG string: the same view, size and window through the context's ports and this viewer's store
   * (data already loaded is shared). At device pixel ratio 1, whatever the screen's: an export is device independent.
   */
  renderSvg(options: Partial<RenderSvgOptions> = {}): Promise<RenderSvgResult> {
    const resolved = this.view.resolved.value;
    const vp = resolved ? this._viewport(resolved) : undefined;
    return this.context.renderSvg(this.spec.shown.value, {
      size: { width: this.width, height: this.height },
      ...(resolved ? { view: resolved.address } : {}),
      ...(vp ? { viewport: { start: vp.start, end: vp.end } } : {}),
      store: this.store.store,
      ...options,
    });
  }
}

function viewerConfig(spec: SeqViewerSpec): ViewerConfig {
  return viewerConfigOf(spec.config);
}

/**
 * Makes a viewer from a spec, as a Mol* plugin is made: its context first (registry, configuration, ports), then the
 * models, the bindings and the behaviors, then the UI mounted into `target` when there is one and the spec has a UI.
 * Without a target it is headless: models, commands, frames, SVG. `retire()` disposes what the viewer made, its owned
 * ports included, and leaves a host's.
 */
export function createSeqViewer(spec: SeqViewerSpec, target?: HTMLElement): SeqViewer {
  return new SeqViewer(spec, target);
}
