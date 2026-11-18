/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { type Issue, issueOf, type Node } from '@molstar/core/repr/tree';
import { createSeqContext, type SeqContext } from '@molstar/seqx-core/context';
import {
  type Emphasis,
  type EmphasisHighlight,
  type Frame,
  type FrameInput,
  type Hit as RenderHit,
  Renderer,
  type RenderSvgOptions,
  type RenderSvgResult,
  type Viewport as RenderViewport,
  hitView,
  resolveColor,
} from '@molstar/seqx-core/render';
import { type Location, locationOfResult } from '@molstar/seqx-core/query';
import type { SeqTokens } from '@molstar/seqx-core/render';
import { type ResolvedView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { cancelFrameCallback, DomHost, requestFrameCallback } from '@molstar/ui/dom';
import {
  bindingTable,
  DefaultActions,
  type BindingTable,
  BindingPresets,
  SelectActions,
} from '@molstar/seqx-viewer/bindings';
import type { CanvasBackend } from '@molstar/seqx-viewer/canvas';
import { type ViewerConfig, viewerConfigOf } from '@molstar/seqx-viewer/config';
import {
  AddressModel,
  HighlightModel,
  HoverModel,
  isNothing,
  Projector,
  SelectionModel,
  SpecModel,
  StoreModel,
  ViewModel,
} from '@molstar/seqx-viewer/models';
import { AddressingProvider, type LinkRun, SelectProvider } from '@molstar/seqx-viewer/behaviors';
import { ViewerInput } from './input.js';
import { uiHistory } from './ui-history.js';
import { ViewerIssues } from './issues.js';
import type { BehaviorContext, BehaviorProvider, SeqViewerSpec } from './viewer-spec.js';
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
  /** The selection (session state): what the selection pass draws. */
  readonly selection: SelectionModel;
  /** Named highlight channels (session state). */
  readonly highlight: HighlightModel;
  /** Addresses in the spec shown: resolved, acted on, linked; the refs; the selection and the focus as addresses. */
  readonly addresses: AddressModel;
  readonly bindings: BindingTable;
  readonly input: ViewerInput;
  /** The behaviors' handles by provider name. */
  readonly behaviors = new Map<string, any>();
  /** The commands the behaviors provide. */
  readonly commands: ViewerCommand[] = [];
  /** The layout's issues of the last frame (a color it cannot read, a layer it cannot draw). */
  readonly layoutIssues: ReadonlyState<readonly Issue[]>;
  /** The viewer's own: a gesture that needs a behavior the spec does not list, once each. */
  readonly behaviorIssues: ReadonlyState<readonly Issue[]>;
  /**
   * Moves once each time the picture comes to rest: when the frame laid out settled is drawn, `canvas.settleMs` after the
   * view's last move; or `canvas.settleMs` after the last layout that changed what is drawn while the view was still
   * (data arriving, the content's height): what follows the picture waits for it.
   */
  readonly rested: ReadonlyState<number>;

  /** Moves to the token table's version each time the tokens change (`setTokens`): what follows the picture redraws. */
  readonly tokensVersion: ReadonlyState<number>;
  /** The token table the frames and the SVG are drawn with: the context's, or the root's CSS when mounted. */
  tokens: SeqTokens;

  /** The canvas's size in CSS px (the data's width and the content's height by default). */
  width = 0;
  height = 0;
  /** The device pixel ratio of the last frame. */
  dpr = 1;
  /** Frames drawn and layouts made: what a test or a profile reads. */
  readonly counts = { frames: 0, layouts: 0 };

  /** Moves when the issues do, not with each layout: what follows them (a story's SVG export) runs once a change. */
  private readonly _issues = new State<readonly Issue[]>([], sameIssues);
  private readonly _behaviorIssues = new State<readonly Issue[]>([], sameIssues);
  private readonly _rested = new State(0);
  private readonly _tokensVersion = new State(0);
  private readonly _handles = new Map<BehaviorProvider<any, any>, unknown>();
  private _frameId: number | undefined;
  private _rehitting = false;
  private readonly _emphasis: Mutable<Emphasis> = {};
  private readonly _projector: Projector;
  /** The emphasis's projected record, reused: its fields are the cached projections. */
  private readonly _projected: { hover?: Location; selection?: Location; hoverSection?: string } = {};
  private readonly _cursor = { system: '', position: 0 };
  /** The highlights as drawn, made again only when a channel, the tokens, the view or the model changed. */
  private _highlights: readonly EmphasisHighlight[] = [];
  private readonly _highlightsFor: unknown[] = [];
  private readonly _emphasisHit = {} as Mutable<RenderHit>;
  /** Section ids by address (the last step), kept: a hover makes no string. */
  private readonly _sectionIds = new Map<string, string>();
  private readonly _size = { width: 0, height: 0 };
  private readonly _frameInput = { size: this._size } as Mutable<FrameInput>;
  private readonly _viewportOut = { start: 0, end: 0 };
  /** The window the last frame drew (NaN before one): where a jump's motion starts. */
  private readonly _drawn = { start: Number.NaN, end: Number.NaN };
  /** A jump's motion (_spec/README.md "Composition", Motion): from the window drawn, since a time, at a jump count. */
  private _motion: { start: number; end: number; t0: number; jumps: number } | undefined;
  /** The platform asks for reduced motion (the mounted UI's media query): every jump is a cut. */
  reducedMotion = false;
  private _backend: CanvasBackend | undefined;
  private _renderer: Renderer | undefined;
  private _rendererFor: unknown[] = [];
  private _wanter: { want(w: readonly any[]): void; dispose(): void } | undefined;
  private _frame: Frame | undefined;
  private _quiet = true;
  private _settleTimer: unknown;
  private _restTimer: unknown;
  /** The view came to rest: `rested` moves once the frame that lays it out settled is drawn. */
  private _restDue = false;
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
    this.behaviorIssues = this._behaviorIssues.asReadonly;
    this.rested = this._rested.asReadonly;
    this.tokensVersion = this._tokensVersion.asReadonly;
    this.tokens = this.context.tokens;

    this.spec = new SpecModel(
      createSeqSpec().build(),
      { runtimes: this.context.runtimes, config: this.context.resolveConfig },
      this.config,
    );
    this.view = new ViewModel(this.spec, this.config);
    this.hover = new HoverModel(this.spec, this.view, this.config);
    this.store = new StoreModel(this.config);
    // an element selected by an address and then toggled is drawn by its extent, resolved over what is shown
    this.selection = new SelectionModel(
      (els) =>
        locationOfResult(this.spec.shown.value, this.addresses.queryModel(), { type: 'elements', elements: els })
          .location,
      { extendScope: this.config.selection.extendScope },
    );
    this._projector = new Projector(this.spec, this.view);
    this.highlight = new HighlightModel();
    this.addresses = new AddressModel(
      this.spec,
      this.view,
      this.store,
      this.hover,
      this.selection,
      this.highlight,
      this.config,
    );
    for (const m of [this.addresses, this.highlight, this.selection, this.store, this.hover, this.view, this.spec])
      this.registerPermanentDispose(() => m.retire());
    // the host's history, through Seq*'s port: the spec takes part in it (none headless unless the spec gives one)
    const history = viewerSpec.history ?? (target && viewerSpec.ui?.history ? uiHistory() : undefined);
    if (history) this.spec.attachHistory(history);

    const preset = viewerSpec.bindings?.preset ?? BindingPresets[this.config.bindings.preset];
    const extraActions = viewerSpec.bindings?.actions ?? [];
    // the selection's actions only with the behavior that runs them: else its bindings are reported, never dead
    const select = viewerSpec.behaviors.some((b) => b.provider === SelectProvider) ? SelectActions : [];
    this.bindings = bindingTable(preset, viewerSpec.bindings?.overrides ?? [], {
      actions: [...DefaultActions, ...select, ...extraActions],
    });

    // behaviors, after the models; one of each name
    const context: BehaviorContext = {
      spec: this.spec,
      view: this.view,
      hover: this.hover,
      store: this.store,
      selection: this.selection,
      highlight: this.highlight,
      addresses: this.addresses,
      config: this.config,
      bindings: this.bindings,
    };
    try {
      for (const b of viewerSpec.behaviors) {
        if (this.behaviors.has(b.provider.name)) throw new Error(`The behavior ${b.provider.name} is listed twice`);
        const handle = b.provider.create(context, b.params);
        this.behaviors.set(b.provider.name, handle);
        this._handles.set(b.provider, handle);
        if (Array.isArray(handle?.commands)) this.commands.push(...handle.commands);
        if (typeof handle?.dispose === 'function') this.registerPermanentDispose(() => handle.dispose());
      }
    } catch (e) {
      // what was made is not left behind
      this.retire();
      throw e;
    }
    this.input = new ViewerInput(this);

    // what moves the picture asks for a frame
    this.listenPermanently(this.spec.model, () => this._sync());
    this.listenPermanently(this.hover.version, () => this.schedule());
    this.listenPermanently(this.selection.value, () => {
      this._holdSelection();
      this.schedule();
    });
    this.listenPermanently(this.selection.cursor, () => this.schedule());
    this.listenPermanently(this.highlight.version, () => this.schedule());
    // a source-backed system's length arrives with its source's coverage: the view's axes are placed then
    this.view.setLengthSource((system) => {
      const s = this.spec.model.value.systems.get(system);
      if (!s?.source) return undefined;
      const info = this.store.store
        .coverage(s.source)
        ?.systems.find((x) => x.name === s.name || x.aliases?.includes(s.name));
      return info?.length;
    });
    this.listenPermanently(this.store.store.version, () => {
      this.view.refreshLengths();
      this.schedule();
    });
    this.listenPermanently(this.store.store.settled, () => this.schedule());
    // a jump animates from the window drawn to the model's (the model is there at once: hits read the frame drawn)
    this.listenPermanently(this.view.jumps, (n) => this._startMotion(n));
    this.listenPermanently(this.view.resolved, () => {
      this._followViewports();
      this._applyWidths();
      this.schedule();
    });
    this._followViewports();
    if (target && viewerSpec.ui?.canvas) this.registerPermanentDispose(mountCanvas(this, target));
    // registered last, so it runs first: everything after it (the UI's teardown, the behaviors, the models) is disposed
    // while the viewer is retired, and none of it asks for a frame
    this.registerPermanentDispose(() => {
      this._retired = true;
      if (this._frameId !== undefined) cancelFrameCallback(this._frameId);
      this._frameId = undefined;
      DomHost.clearTimeout(this._settleTimer);
      this._settleTimer = undefined;
      DomHost.clearTimeout(this._restTimer);
      this._restTimer = undefined;
      this._wanter?.dispose();
      for (const u of this._viewportListeners) u();
      this._viewportListeners.length = 0;
    });
  }

  /**
   * The handle of a behavior by its provider (what the input calls). One the spec does not list is reported once as an
   * issue (`behaviorIssues`) and the gesture does nothing.
   */
  behavior<H>(provider: BehaviorProvider<any, H>): H | undefined {
    if (this._handles.has(provider)) return this._handles.get(provider) as H;
    if (!this._behaviorIssues.value.some((i) => i.address === provider.name))
      this._behaviorIssues.set([
        ...this._behaviorIssues.value,
        issueOf(
          ViewerIssues,
          'seq.viewer.behavior.missing',
          `The viewer has no ${provider.name} behavior: its gestures do nothing`,
          { address: provider.name, hint: ViewerIssues.codes['seq.viewer.behavior.missing'].hint },
        ),
      ]);
    return undefined;
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
    // a new root syncs the sources through the model's own listener: once
    const refused = this.spec.apply(root);
    if (refused === 'out-of-bounds')
      throw Object.assign(new Error('The spec is past the configured bounds'), {
        issues: this.spec.boundsIssues.value,
      });
    this.view.show(options.view);
    // a new document: the view's history was the old one's (an edit of the same spec keeps it)
    this.view.clearHistory();
    await this._syncing;
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

  /**
   * Sets the token table (a theme switch). A table equal to the current one, by its version, is no change; else the
   * renderer keeps its layout key's token version, the backend redraws and a frame is asked for.
   */
  setTokens(tokens: SeqTokens): boolean {
    if (this._retired || tokens.version === this.tokens.version) return false;
    this.tokens = tokens;
    if (this._renderer) this._renderer.tokens = tokens;
    this._baseDirty = true;
    this._tokensVersion.set(tokens.version);
    this.schedule();
    return true;
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
          // a move by hand ends a jump's motion (a jump's own move is followed by its count moving)
          if (this._motion && this._motion.jumps === this.view.jumps.value) this._motion = undefined;
          this._moved();
          this.schedule();
        }),
      );
  }

  /** The view moved: it is not settled until it rests for `canvas.settleMs`. */
  private _moved(): void {
    if (this._retired) return;
    this._quiet = false;
    this._restDue = false;
    DomHost.clearTimeout(this._settleTimer);
    this._settleTimer = DomHost.setTimeout(() => {
      this._settleTimer = undefined;
      this._quiet = true;
      this._restDue = true;
      this.schedule();
    }, this.config.canvas.settleMs);
  }

  /** The picture changed while the view is quiet (data arrived, the height followed): it rests `canvas.settleMs` later. */
  private _changed(): void {
    if (this._retired || !this._quiet) return;
    DomHost.clearTimeout(this._restTimer);
    this._restTimer = DomHost.setTimeout(() => {
      this._restTimer = undefined;
      if (this._quiet && !this._retired) this._rested.set(this._rested.value + 1);
    }, this.config.canvas.settleMs);
  }

  /** Asks for a frame: one a frame (`@molstar/ui/dom`'s frame callbacks, on `DomHost`), whatever asked; `retire` cancels it. */
  schedule(): void {
    // a re-hit's own hover change is drawn by the frame that made it
    if (this._frameId !== undefined || this._retired || this._rehitting) return;
    this._frameId = requestFrameCallback(() => {
      this._frameId = undefined;
      if (!this._retired) this.draw();
    });
  }

  /** A jump moved the view: its motion starts from the window drawn, unless motion is cut. */
  private _startMotion(jumps: number): void {
    const { jumpMs, motion } = this.config.view;
    if (this._retired || motion === 'cut' || this.reducedMotion || !(jumpMs > 0) || Number.isNaN(this._drawn.start)) {
      this._motion = undefined;
      return;
    }
    this._motion = { start: this._drawn.start, end: this._drawn.end, t0: DomHost.now(), jumps };
    this.schedule();
  }

  /**
   * Stops a jump's motion where it is drawn (_spec/README.md "Composition", Motion): the leading section takes the window
   * drawn now (a move, no step of the history) and a frame is drawn, so the input that stopped it reads one picture: what
   * it hits is what the model shows. What every pointer and key input does first.
   */
  stopMotion(): void {
    const m = this._motion;
    const resolved = this.view.resolved.value;
    if (!m || !resolved) return;
    const vp = this._viewport(resolved, DomHost.now());
    this._motion = undefined;
    if (!vp) return;
    const lead = resolved.sections.find((s) => !this.view.isFollowing(s.address));
    if (lead) this.view.setRange(lead.address, vp.start, vp.end);
    this.draw();
  }

  /** Whether a jump's motion is drawing. */
  get moving(): boolean {
    return this._motion !== undefined;
  }

  /** The window the last frame drew (during a jump's motion, between where it started and the model's), a copy. */
  get drawnViewport(): RenderViewport | undefined {
    return Number.isNaN(this._drawn.start) ? undefined : { start: this._drawn.start, end: this._drawn.end };
  }

  /**
   * The window the renderer draws: the first section that leads (the renderer takes one window for all); during a
   * jump's motion (`at` the time), the window between where it started and the model's, eased in and out, its width
   * moving geometrically and its centre along it so the content under the eye moves smoothly. Undefined past the end.
   */
  private _viewport(resolved: ResolvedView, at?: number): RenderViewport | undefined {
    for (const s of resolved.sections) {
      if (this.view.isFollowing(s.address)) continue;
      const v = this.view.get(s.address);
      if (!v || v.pending) return undefined;
      const out = this._viewportOut;
      out.start = v.start;
      out.end = v.end;
      const m = this._motion;
      if (at !== undefined && m) {
        const t = (at - m.t0) / this.config.view.jumpMs;
        if (t >= 1 || !(t >= 0)) this._motion = undefined;
        else {
          const e = t < 0.5 ? 4 * t * t * t : 1 - (-2 * t + 2) ** 3 / 2;
          const w0 = m.end - m.start;
          const w1 = v.end - v.start;
          const w = w0 > 0 && w1 > 0 ? w0 * (w1 / w0) ** e : w0 + (w1 - w0) * e;
          const c = (m.start + m.end) / 2 + ((v.start + v.end) / 2 - (m.start + m.end) / 2) * e;
          out.start = c - w / 2;
          out.end = c + w / 2;
        }
      }
      return out;
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
        tokens: this.tokens,
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
    const vp = this._viewport(resolved, DomHost.now());
    if (vp) input.viewport = vp;
    else delete input.viewport;
    if (vp) {
      this._drawn.start = vp.start;
      this._drawn.end = vp.end;
    }
    // a motion draws a frame a frame until it ends; the view is not settled meanwhile
    if (this._motion) {
      this._moved();
      this.schedule();
    }
    input.dpr = dpr;
    input.settled = this._quiet && this.store.store.settled.value;
    input.emphasis = this._emphasis;
    // tier 2: the store's dynamic nodes (its version is part of the frame's key)
    input.overlay = this.store.overlay(model, resolved, this.context);

    let frame = renderer.frame(input);
    let rebuilt = frame.rebuilt;
    if (!fixed) {
      const last = frame.layout.sections.at(-1);
      const content = last ? Math.ceil(last.y + last.height + this.context.renderConfig.inset) : 32;
      if (content !== this._contentHeight) {
        this._contentHeight = content;
        size.height = Math.min(this.config.canvas.maxHeight, Math.max(content, 32));
        frame = renderer.frame(input);
        rebuilt ||= frame.rebuilt;
      }
    }
    this.height = size.height;
    if (rebuilt) {
      this.counts.layouts++;
      if (!this._restDue) this._changed();
      this._wanter?.want(frame.wants);
      this._issues.set(frame.layout.issues);
      this.store.holdShown(frame.layout.shown);
    }
    this._frame = frame;

    // a resting pointer is hit again once a frame the view moves under it (EC-143)
    // (what that frame laid out is still to be drawn: `rebuilt` is the layout's, not the emphasis redraw's)
    if (rebuilt && this.config.pointer.rehit) {
      this._rehitting = true;
      const changed = this.input.rehit();
      this._rehitting = false;
      if (changed) {
        this._updateEmphasis();
        frame = renderer.frame(input);
        this._frame = frame;
      }
    }

    const backend = this._backend;
    if (backend) {
      backend.fonts.setStacks(renderer.tokens.fonts);
      this._baseDirty = backend.resize(size.width, size.height, dpr) || this._baseDirty;
      backend.draw(frame.lists, rebuilt || frame.base || this._baseDirty);
      this._baseDirty = false;
    }
    this.counts.frames++;
    if (this._restDue) {
      // the settled picture is drawn: it rests now, once (a layout of this frame waits for no second rest)
      this._restDue = false;
      DomHost.clearTimeout(this._restTimer);
      this._restTimer = undefined;
      this._rested.set(this._rested.value + 1);
    }
    return frame;
  }

  /** The last frame: its layout is what hit-testing reads. */
  get frame(): Frame | undefined {
    return this._frame;
  }

  /**
   * The emphasis the passes draw: the selection's Location and its projection onto the linked sections, the cursor, the
   * highlight channels (each with its projection), and the hover with its projection. Reused records, mutated when what
   * they show changes; the projections are cached (`Projector`, `HoverModel.projectedAll`), so a redraw makes none.
   */
  private _updateEmphasis(): void {
    const e = this._emphasis;
    const p = this._projected;
    const sel = this.selection.value.value;
    e.elements = sel.elements.keys.length ? sel.elements : undefined;
    if (isNothing(sel)) e.selection = p.selection = undefined;
    else {
      e.selection = sel.shown;
      p.selection = this._projector.of(sel.shown);
    }
    // the cursor shows once the keyboard moved it (focus-visible: a click places it unseen)
    const c = this.selection.cursor.value;
    if (c?.keyboard) {
      this._cursor.system = c.system;
      this._cursor.position = c.position;
      e.cursor = this._cursor;
    } else e.cursor = undefined;
    this._updateHighlights();
    e.projected = p;
    p.hover = this.hover.projectedAll();
    p.hoverSection = this.hover.current.section;
    const v = this.hover.version.value;
    if (v === this._hoverVersion) return;
    this._hoverVersion = v;
    const r = this.hover.current;
    if (r.owner === undefined || r.layer === undefined || r.system === undefined) {
      e.hover = undefined;
      return;
    }
    const h = this._emphasisHit;
    h.layer = r.layer;
    h.section = r.section === undefined ? '' : this._sectionId(r.section);
    h.system = r.system;
    (h as { position: number | undefined }).position = Number.isNaN(r.position) ? undefined : r.position;
    (h as { element: unknown }).element = r.element;
    (h as { extent: unknown }).extent = r.extent;
    e.hover = h;
  }

  /** The highlight channels as drawn: each its Location on every section it shows on, in its channel's color. */
  private _updateHighlights(): void {
    const k = this._highlightsFor;
    const version = this.highlight.version.value;
    const view = this.view.resolved.value;
    const model = this.spec.model.value;
    if (k[0] === version && k[1] === this.tokens.version && k[2] === view && k[3] === model) {
      this._emphasis.highlights = this._highlights.length ? this._highlights : undefined;
      return;
    }
    k[0] = version;
    k[1] = this.tokens.version;
    k[2] = view;
    k[3] = model;
    const out: EmphasisHighlight[] = [];
    const { color, colors } = this.config.highlight;
    for (const name of this.highlight.names.value) {
      const h = this.highlight.get(name);
      if (!h || isNothing(h)) continue;
      const location = this._projector.of(h.shown, { own: true }) ?? h.shown;
      const c =
        resolveColor(Object.hasOwn(colors, name) ? colors[name] : color, this.tokens) ??
        resolveColor(color, this.tokens) ??
        this.tokens.colors['seq-cat-5'];
      out.push({ location, color: c });
    }
    this._highlights = out;
    this._emphasis.highlights = out.length ? out : undefined;
  }

  /** The selection's dynamic nodes held while it is selected (`selection.holdNodes`). */
  private _holdSelection(): void {
    if (!this.config.selection.holdNodes || this._retired) return;
    const sel = this.selection.value.value;
    const addresses: string[] = [];
    for (const k of sel.elements.keys) addresses.push(k.owner);
    for (const part of sel.shown.parts) addresses.push(part.system);
    this.store.holdSession('selection', this.spec.model.value, addresses);
  }

  /** A section's id, its address's last step (kept: a hover allocates nothing once warm). */
  private _sectionId(address: string): string {
    let id = this._sectionIds.get(address);
    if (id === undefined) this._sectionIds.set(address, (id = address.slice(address.lastIndexOf('/') + 1)));
    return id;
  }

  /** What is under a point of the last frame's layout. */
  hit(x: number, y: number): RenderHit | undefined {
    return this._frame ? hitView(this._frame.layout, x, y) : undefined;
  }

  // ---- commands ----

  /**
   * Runs a command a behavior provides; the message of its result, or why it could not run. A command that runs in
   * turn (an edit through the host's history) says nothing here: `execute` waits for its message.
   */
  runCommand(id: string, context: unknown): string {
    const r = this.execute(id, context);
    return typeof r === 'string' ? r : '';
  }

  /** Runs a command; its message, at once for a command that runs at once, else when it has run. */
  execute(id: string, context: unknown): string | Promise<string> {
    const c = this.commands.find((x) => x.id === id);
    if (!c) return `No command ${id}`;
    const can = c.when?.(context);
    if (can !== undefined && can !== true) return typeof can === 'string' ? can : 'Not available';
    const r = c.run({}, context as never, NeverAborted);
    return r instanceof Promise ? r.then((x: CommandResult) => x.message ?? '') : (r.message ?? '');
  }

  /**
   * Runs a link's actions (`!seq-focus=/refs/dbd`), each through the action its command runs, in order: undefined when
   * the text is no link or the viewer has no `addressing` behavior (reported in `behaviorIssues`).
   */
  runLink(text: string): LinkRun | undefined {
    if (this._retired) return undefined;
    return this.behavior(AddressingProvider)?.runLink(text);
  }

  // ---- export ----

  /** The session as an export draws it: the selection (its address, else its Location) and the highlight channels. */
  private _svgEmphasis(): Pick<RenderSvgOptions, 'selection' | 'emphasis'> {
    this._updateHighlights();
    const sel = this.selection.value.value;
    const highlights = this._emphasis.highlights;
    const location = sel.address === undefined && !isNothing(sel) ? sel.shown : undefined;
    return {
      ...(sel.address !== undefined ? { selection: sel.address } : {}),
      ...(location || highlights
        ? { emphasis: { ...(location ? { selection: location } : {}), ...(highlights ? { highlights } : {}) } }
        : {}),
    };
  }

  /**
   * The picture as an SVG string: the same view, size and window through the context's ports and this viewer's store
   * (data already loaded is shared), at the device pixel ratio of the last frame: the density steps are in device
   * pixels, so this is the picture the canvas shows (`dpr: 1` for a device independent export).
   */
  renderSvg(options: Partial<RenderSvgOptions> = {}): Promise<RenderSvgResult> {
    const resolved = this.view.resolved.value;
    const vp = resolved ? this._viewport(resolved) : undefined;
    return this.context.renderSvg(this.spec.shown.value, {
      size: { width: this.width, height: this.height },
      ...(resolved ? { view: resolved.address } : {}),
      ...(vp ? { viewport: { start: vp.start, end: vp.end } } : {}),
      store: this.store.store,
      dpr: this.dpr,
      tokens: this.tokens,
      // what the canvas shows of the session: the selection, as its address when it has one (the export resolves it
      // over what it holds, and projects it onto the linked sections), else as its Location; the highlight channels
      ...this._svgEmphasis(),
      ...options,
    });
  }
}

/** The same issues in the same order, by code, address and message. */
function sameIssues(a: readonly Issue[], b: readonly Issue[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++)
    if (a[i].code !== b[i].code || a[i].address !== b[i].address || a[i].message !== b[i].message) return false;
  return true;
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
