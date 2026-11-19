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
  BuiltInLayouts,
  type Emphasis,
  type EmphasisHighlight,
  type BandViewport,
  type Frame,
  type FrameInput,
  type Hit as RenderHit,
  type Placement,
  Renderer,
  type RenderSvgOptions,
  type RenderSvgResult,
  type Viewport as RenderViewport,
  hitView,
  railSections,
  rowsRailTracks,
  resolveColor,
  sectionLayoutName,
  textMinScale,
} from '@molstar/seqx-core/render';
import { type Location, locationOfResult } from '@molstar/seqx-core/query';
import type { SeqTokens } from '@molstar/seqx-core/render';
import {
  type Frame as StoryFrame,
  frameRoot,
  frameView,
  plan as planTransition,
  type PlanOptions,
  sample,
  type StoryState,
  type StepViewport,
  type TransitionPlan,
  type ViewState,
} from '@molstar/seqx-core/story';
import type { ResolvedNavigator, ResolvedView } from '@molstar/seqx-core/view';
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
  axisTotal,
  HighlightModel,
  HoverModel,
  isNothing,
  nearestOffsetOf,
  Projector,
  SelectionModel,
  SpecModel,
  StoreModel,
  ViewModel,
} from '@molstar/seqx-viewer/models';
import {
  AddressingProvider,
  type LinkRun,
  NavigatorProvider,
  windowIdle,
  RowsActions,
  RowsProvider,
  SelectProvider,
} from '@molstar/seqx-viewer/behaviors';
import type { AnimationRequest, SeqAnimation, SeqAnimationParticipant } from '@molstar/seqx-viewer/ports';
import { frameAnimation } from './animation.js';
import { ViewerInput } from './input.js';
import { highlightColorOf, LegendModel } from './legends.js';
import { uiHistory } from './ui-history.js';
import { ViewerIssues } from './issues.js';
import type { BehaviorContext, BehaviorProvider, SeqViewerSpec } from './viewer-spec.js';
import { mountCanvas } from './mount.js';
import { type NavigatorDef, NavigatorModel } from './navigator.js';

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
  /** The legends (derived on settle), the hovered entry and the reader's color-blind substitution. */
  readonly legends: LegendModel;
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
  /** The section layouts by name: one map the frames draw with and the view's moves navigate by. */
  readonly layouts = BuiltInLayouts;
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
  /** The window the last frame drew (NaN before one). */
  private readonly _drawnWindow = { start: Number.NaN, end: Number.NaN };
  /** Each leading section's viewport the last frame drew, by address: where a motion starts. */
  private readonly _drawnViews = new Map<string, StepViewport>();
  /**
   * The frame slot (_spec/README.md "Composition", Motion): the plan and frame shown over the state, its viewports by
   * section; a request on the animation port running it (a jump's, with its count).
   */
  private _shownFrame: { plan: TransitionPlan; frame: StoryFrame; views?: ViewState } | undefined;
  private _request: { request: AnimationRequest; plan: TransitionPlan; jumps?: number } | undefined;
  private readonly _drawnFrames = new State(0);
  /** Bumped after every frame drawn: a navigator redraws over it. */
  private readonly _drawn = new State(0);
  readonly drawn: ReadonlyState<number> = this._drawn.asReadonly;
  /**
   * The ruler rails shown (_spec/README.md "Composition", Navigators): the reader's session (the embedded layout's
   * corner button), never the spec; `navigator.form: none` draws none whatever it says.
   */
  private readonly _railsShown = new State(true);
  readonly railsShown: ReadonlyState<boolean> = this._railsShown.asReadonly;
  /** The sections that draw a rail and which can move: kept while they hold (the frame's key reads them by identity). */
  private _rails: Record<string, true> | undefined;
  private _railsFor: [unknown, boolean, number] = [undefined, false, 0];
  /** Bumped when a host's navigator comes or goes: the rails it replaces are found again. */
  private _hostedVersion = 0;
  private _railsActive: Record<string, boolean> = {};
  /** The rails' keys that are rows tracks' (a rows rail). */
  private _rowsRails: ReadonlySet<string> = new Set();
  /** The navigators' models by id: the view's declared ones as asked for, and a host's (`navigator(def)`). */
  private readonly _navigators = new Map<string, NavigatorModel>();
  /** The spec the last frame drew, the step shown, and whether `playStep` is entering one (no second transition). */
  private _drawnRoot: Node | undefined;
  private _shownStep: string | undefined;
  private _stepping = false;
  /** The application's clock (the animation port): Seq*'s own motions are plans requested of it. */
  private readonly _animation: SeqAnimation;
  /** The participant face of the animation port: what the application's clock drives. */
  readonly animation: SeqAnimationParticipant;
  /**
   * The platform asks for reduced motion (the mounted UI's media query): under `view.motion: platform` every jump and
   * story step is a cut.
   */
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
      false,
    );
    this.view = new ViewModel(this.spec, this.config);
    // a section's moves are its layout's arithmetic over its placement in the last frame (_spec/README.md "ViewModel")
    this.view.setLayoutSource({
      placement: (section) => this.placementOf(section),
      layouts: this.layouts,
      minScale: textMinScale(this.context.renderConfig),
    });
    this.hover = new HoverModel(this.spec, this.view, this.config);
    this.store = new StoreModel(this.config);
    // an element selected by an address and then toggled is drawn by its extent, resolved over what is shown
    this.selection = new SelectionModel(
      (els) =>
        locationOfResult(this.spec.shown.value, this.addresses.queryModel(), { type: 'elements', elements: els })
          .location,
      { extendScope: this.config.selection.extendScope, acrossBreaks: this.config.selection.acrossBreaks },
    );
    this._projector = new Projector(this.spec, this.view);
    this.highlight = new HighlightModel();
    this.legends = new LegendModel(this.config.theme.cvdSafe);
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
    const rows = viewerSpec.behaviors.some((b) => b.provider === RowsProvider) ? RowsActions : [];
    this.bindings = bindingTable(preset, viewerSpec.bindings?.overrides ?? [], {
      actions: [...DefaultActions, ...select, ...rows, ...extraActions],
      keys: this.config.view.keys,
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
      legends: this.legends,
      layout: () => this._frame?.layout,
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
      this._revealLinked();
      this.schedule();
    });
    this.listenPermanently(this.selection.cursor, () => this.schedule());
    this.listenPermanently(this.highlight.version, () => this.schedule());
    // the reader's substitution redraws (its frames are keyed by it); the hover moves the legend's entry
    this.listenPermanently(this.legends.cvdSafe, () => this.schedule());
    // a band's own scroll, row zoom and folds (a rows band's)
    this.listenPermanently(this.view.bands, () => this.schedule());
    this.listenPermanently(this.hover.version, () => this.legends.updateHover(this._frame?.layout, this.hover.current));
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
    // the animation port (_spec/README.md "Ports", Animation): the application's clock, the page's frames by default
    this.animation = {
      progress: (plan, t) => this._showFrame(plan, sample(plan, t)),
      show: (plan, frame) => this._showFrame(plan, frame),
      clear: () => this._clearFrame(),
      draw: () => {
        this.draw();
      },
      drawn: this._drawnFrames.asReadonly,
      durationOf: (plan) => plan.duration,
    };
    this._animation = (viewerSpec.animation ?? frameAnimation())(this.animation);
    this.registerPermanentDispose(() => {
      this._request?.request.cancel();
      this._request = undefined;
      if (this._animation.disposeWithViewer !== false) this._animation.dispose?.();
    });
    // a jump animates from the windows drawn to the model's (the model is there at once: hits read the frame drawn)
    this.listenPermanently(this.view.jumps, (n) => this._startMotion(n));
    // an edit, an undo or a load under a frame that changes the spec ends it: the frame was made from the old root
    this.listenPermanently(this.spec.head, () => {
      const f = this._shownFrame;
      if (!f || f.plan.from.root === f.plan.to.root) return;
      this._endMotion();
      this._clearFrame();
    });
    // a story step entered by a command or a host is a transition too (its view is a snapshot's: a jump)
    this.listenPermanently(this.spec.step, () => this._stepChanged());
    this.listenPermanently(this.view.resolved, () => {
      this._followViewports();
      this._applyWidths();
      this.schedule();
    });
    this._followViewports();
    if (target && viewerSpec.ui?.toolbar) this.registerPermanentDispose(viewerSpec.ui.toolbar(this, target));
    if (target && viewerSpec.ui?.canvas) this.registerPermanentDispose(mountCanvas(this, target));
    if (target && viewerSpec.ui?.navigators) this.registerPermanentDispose(viewerSpec.ui.navigators(this, target));
    if (target && viewerSpec.ui?.railToggle) this.registerPermanentDispose(viewerSpec.ui.railToggle(this, target));
    if (target && viewerSpec.ui?.legends) this.registerPermanentDispose(viewerSpec.ui.legends(this, target));
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
    // a new document: the story step shown was the old one's (and a full-state step refuses a whole root)
    this.spec.leaveStep();
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
          // a move by hand ends a motion, a story's too, and clears a frame a host shows (a jump's own move is
          // followed by its count moving; a step's view is placed before its motion starts)
          if (this._request && (this._request.jumps === undefined || this._request.jumps === this.view.jumps.value))
            this._endMotion();
          else if (!this._request && this._shownFrame) this._clearFrame();
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

  /**
   * Whether a motion, a jump's or a story's, is a cut: `view.motion` (`cut`; `platform`: the platform's reduced motion;
   * `animate`: never), one item for both (_spec/README.md, Configuration).
   */
  private get _cut(): boolean {
    const motion = this.config.view.motion;
    return motion === 'cut' || (motion === 'platform' && this.reducedMotion);
  }

  /** Each leading section's viewport in the model now (what a jump goes to), by address. */
  private _modelViews(): ViewState {
    const sections: Record<string, StepViewport> = {};
    for (const s of this.view.resolved.value?.sections ?? []) {
      if (this.view.isFollowing(s.address)) continue;
      const v = this.view.get(s.address);
      if (!v || v.pending) continue;
      // the bands' top and zoom (display rows, px a row): what a step tweens; the folds are the session's, never a step's
      const own = this.view.bandsOf(s.address);
      let bands: Record<string, { bandTop: number; rowScale?: number }> | undefined;
      for (const track in own ?? {}) {
        const b = own![track];
        (bands ??= {})[track] = { bandTop: b.bandTop, ...(b.rowScale !== undefined ? { rowScale: b.rowScale } : {}) };
      }
      sections[s.address] = {
        start: v.start,
        end: v.end,
        ...(v.top ? { top: v.top } : {}),
        ...(v.scale !== undefined ? { scale: v.scale } : {}),
        ...(bands ? { bands } : {}),
      };
    }
    return { sections };
  }

  /** What a plan of this viewer needs besides its states: the cyclic sections' totals (a window's shortest arc). */
  private _planOptions(): Pick<PlanOptions, 'cyclic' | 'config' | 'resolveColor'> {
    const cyclic: Record<string, number> = {};
    for (const s of this.view.resolved.value?.sections ?? [])
      if (this.view.wraps(s.address)) {
        const total = axisTotal(s.axis);
        if (total) cyclic[s.address] = total;
      }
    // a token blends through the color it draws in (the viewer's token table, light or dark as shown)
    const tokens = this.tokens;
    const color = (c: string) => {
      const r = resolveColor(c, tokens);
      // RGBA: an opaque color as hex; one with alpha is no color a plan blends
      return r === undefined || (r & 0xff) !== 0xff ? undefined : `#${(r >>> 8).toString(16).padStart(6, '0')}`;
    };
    return { cyclic, config: this.config.story, resolveColor: color };
  }

  /**
   * A jump moved the view (_spec/README.md "Composition", Motion): a plan from every leading section's viewport drawn
   * (window, letter size) to the model's, requested of the animation port over `view.jumpMs`, unless motion is cut.
   */
  private _startMotion(jumps: number): void {
    const { jumpMs } = this.config.view;
    this._endMotion();
    if (this._retired || this._cut || !(jumpMs > 0) || !this._drawnViews.size) return;
    const from: ViewState = { sections: Object.fromEntries(this._drawnViews) };
    const shown = this.spec.shown.value;
    const p = planTransition(
      { root: shown, view: from },
      { root: shown, view: this._modelViews() },
      { ...this._planOptions(), duration: jumpMs },
    );
    if (p.channels.length) this._run(p, jumps);
  }

  /** Runs a plan on the application's clock, shown from its start now; the frame left at its end. */
  private _run(plan: TransitionPlan, jumps?: number): void {
    this._request?.request.cancel();
    this._showFrame(plan, sample(plan, 0));
    const request = this._animation.request(plan, {
      duration: plan.duration,
      done: () => {
        if (this._request?.plan !== plan) return;
        this._request = undefined;
        this._clearFrame();
      },
    });
    this._request = { request, plan, ...(jumps !== undefined ? { jumps } : {}) };
  }

  /** Ends a running motion where the model is (no frame left over it). */
  private _endMotion(): void {
    if (!this._request) return;
    this._request.request.cancel();
    this._request = undefined;
    this._clearFrame();
  }

  /** Shows a frame over the state (the frame slot): its spec when the plan changes the spec, its viewports. */
  private _showFrame(plan: TransitionPlan, frame: StoryFrame): void {
    if (this._retired) return;
    // a plan between two equal roots (a jump's) shows the state itself
    this.spec.setFrame(plan.from.root === plan.to.root ? undefined : frameRoot(plan, frame));
    const views = frameView(plan, frame);
    this._shownFrame = { plan, frame, ...(views ? { views } : {}) };
    this._moved();
    this.schedule();
  }

  /** Leaves the frame: the state is shown and drawn again. */
  private _clearFrame(): void {
    if (!this._shownFrame && this.spec.frame === undefined) return;
    this._shownFrame = undefined;
    this.spec.setFrame(undefined);
    this.schedule();
  }

  /**
   * Stops a motion (a jump's, a story step's) or a host's frame where it is drawn (_spec/README.md "Composition",
   * Motion): each leading section takes the viewport drawn now (a move, no step of the history), the spec shown is the
   * state again, and a frame is drawn, so the input that stopped it reads one picture: what it hits is what the model
   * shows. What every pointer and key input does first.
   */
  stopMotion(): void {
    const r = this._request;
    if (!r && !this._shownFrame) return;
    const drawn = [...this._drawnViews];
    r?.request.cancel();
    this._request = undefined;
    this._clearFrame();
    for (const [address, v] of drawn) this.view.setViewport(address, v);
    this.draw();
  }

  /**
   * Enters a story step, or leaves the step shown (`undefined`), as a transition (_spec/README.md "Composition",
   * Stories): navigation, no commit, no history step. The step's state is shown, its view's sections placed (positions
   * to offsets, exactly), and a plan from what was drawn to it is made, timed by the step. `run` (true): the plan runs
   * on the animation port; false: its start is shown and the application drives it (`animation.progress`, then
   * `animation.clear()`). A cut under reduced motion. The issue when the step is unknown.
   */
  playStep(id: string | undefined, options: { readonly run?: boolean } = {}): { plan?: TransitionPlan; issue?: Issue } {
    this._endMotion();
    const from = this._fromState(this.spec.shown.value);
    this._stepping = true;
    try {
      if (id === undefined) this.spec.leaveStep();
      else {
        const issue = this.spec.enterStep(id);
        if (issue) return { issue };
      }
    } finally {
      this._stepping = false;
    }
    this._shownStep = this.spec.step.value?.id;
    return { plan: this._toStep(from, options.run ?? true) };
  }

  /** What a transition starts from: a root and the viewports drawn last (the model's before a frame). */
  private _fromState(root: Node): StoryState {
    const views = this._drawnViews.size ? { sections: Object.fromEntries(this._drawnViews) } : this._modelViews();
    return { root, view: { ...views, ...this._viewAddress() } };
  }

  /** A step entered elsewhere (the `step.enter` command, a host): the same transition, from what was drawn. */
  private _stepChanged(): void {
    const id = this.spec.step.value?.id;
    if (this._stepping || this._retired || id === this._shownStep) return;
    this._shownStep = id;
    this._endMotion();
    this._toStep(this._fromState(this._drawnRoot ?? this.spec.shown.value), true);
  }

  /** The step shown (its view placed) and the plan to it from `from`, run or shown at its start. */
  private _toStep(from: StoryState, run: boolean): TransitionPlan {
    const st = this.spec.step.value?.state;
    if (st?.view) {
      if (st.view.view !== undefined) this.view.show(st.view.view);
      const resolved = this.view.resolved.value;
      for (const [address, v] of Object.entries(st.view.sections)) {
        const axis = resolved?.sections.find((s) => s.address === address)?.axis;
        const s = axis && nearestOffsetOf(axis, v.start);
        const e = axis && nearestOffsetOf(axis, v.end);
        if (s === undefined || e === undefined) continue;
        this.view.setViewport(address, {
          start: s,
          end: e,
          ...(v.top !== undefined ? { top: v.top } : {}),
          ...(v.scale !== undefined ? { scale: v.scale } : {}),
        });
        // a step's bands: their top and zoom placed in the view model, the session's folds kept
        for (const track in v.bands ?? {}) {
          const b = v.bands![track];
          this.view.setBand(track, {
            ...(b.bandTop !== undefined ? { bandTop: b.bandTop } : {}),
            ...(b.rowScale !== undefined ? { rowScale: b.rowScale } : {}),
          });
        }
      }
    }
    const to = { root: this.spec.shown.value, view: { ...this._modelViews(), ...this._viewAddress() } };
    const timing = st?.timing;
    const p = planTransition(from, to, {
      ...this._planOptions(),
      reducedMotion: this._cut,
      ...(timing ? { duration: timing.duration, easing: timing.easing } : {}),
    });
    if (this._cut || !p.channels.length) {
      this.schedule();
      return p;
    }
    if (run) this._run(p);
    else this._showFrame(p, sample(p, 0));
    return p;
  }

  /** The view shown's address, when one is shown by address. */
  private _viewAddress(): { view?: string } {
    const v = this.view.view.value;
    return v !== undefined ? { view: v } : {};
  }

  /** Whether a motion (a jump's, a story step's) is drawing. */
  get moving(): boolean {
    return this._request !== undefined;
  }

  /** The window the last frame drew (during a motion, between where it started and the model's), a copy. */
  get drawnViewport(): RenderViewport | undefined {
    return Number.isNaN(this._drawnWindow.start)
      ? undefined
      : { start: this._drawnWindow.start, end: this._drawnWindow.end };
  }

  /** Each leading section's viewport the last frame drew, by address (during a motion, the frame's), copies. */
  get drawnViewports(): Readonly<Record<string, StepViewport>> {
    return Object.fromEntries([...this._drawnViews].map(([k, v]) => [k, { ...v }]));
  }

  /**
   * A leading section's viewport as drawn now: the frame's (an animation's), else the model's. `eye` false: a jump's
   * motion is the eye's way only (the SVG and `viewport` read the model, the destination), a story's frame is shown.
   */
  private _shownViewport(address: string, eye = true): StepViewport | undefined {
    const jump =
      this._request !== undefined && this._request.jumps !== undefined && this._request.plan === this._shownFrame?.plan;
    const f = eye || !jump ? this._shownFrame?.views?.sections[address] : undefined;
    if (f) return f;
    const v = this.view.get(address);
    return v && !v.pending ? v : undefined;
  }

  /**
   * The window the renderer draws: the first section that leads (the renderer takes one window for all), as the frame
   * slot shows it during a motion. Undefined when it is not placed.
   */
  private _viewport(resolved: ResolvedView, eye = true): RenderViewport | undefined {
    for (const s of resolved.sections) {
      if (this.view.isFollowing(s.address)) continue;
      const w = this._shownViewport(s.address, eye);
      if (!w) return undefined;
      const out = this._viewportOut;
      out.start = w.start;
      out.end = w.end;
      return out;
    }
    return undefined;
  }

  /**
   * Every leading section's viewport (the frame's `viewports`, by address), as its layout reads it (a window; a text's
   * anchor and letter size; a top, a band's scroll), as the frame slot shows it during a motion: each section's own,
   * field by field. Computed once a frame and kept while equal, so the frame's key holds; undefined when no section is
   * placed.
   */
  private _viewports(resolved: ResolvedView, eye = true): Readonly<Record<string, RenderViewport>> | undefined {
    const kept = this._viewportsKept;
    let n = 0;
    let same = kept !== undefined;
    for (const s of resolved.sections) {
      // a grid's own y window: its entry even when x follows (the frame reads its `y`)
      const y = this.view.get(s.address)?.y;
      if (this.view.isFollowing(s.address) && !y) continue;
      const w = this._shownViewport(s.address, eye);
      if (!w) continue;
      const k = kept?.[s.address];
      if (
        !k ||
        k.start !== w.start ||
        k.end !== w.end ||
        k.scale !== w.scale ||
        k.top !== (w.top || undefined) ||
        k.y !== y ||
        k.bands !== this._bandsFor(s.address, w)
      )
        same = false;
      n++;
    }
    if (!n) return (this._viewportsKept = undefined);
    if (same && Object.keys(kept!).length === n) return kept;
    const out: Record<string, RenderViewport> = {};
    for (const s of resolved.sections) {
      const y = this.view.get(s.address)?.y;
      if (this.view.isFollowing(s.address) && !y) continue;
      const w = this._shownViewport(s.address, eye);
      if (!w) continue;
      const bands = this._bandsFor(s.address, w);
      out[s.address] = {
        start: w.start,
        end: w.end,
        ...(w.scale !== undefined ? { scale: w.scale } : {}),
        // `top`: the frame's during a motion, else the model's (when set)
        ...(w.top ? { top: w.top } : {}),
        // `bands`: the ViewModel's (scroll, row zoom, folds), a shown frame's `bandTop` and `rowScale` over them
        ...(bands ? { bands } : {}),
        ...(y ? { y } : {}),
      };
    }
    this._viewportsKept = out;
    return out;
  }

  private _viewportsKept: Readonly<Record<string, RenderViewport>> | undefined;

  /**
   * A section's bands as a frame reads them (_spec/README.md "ViewModel", Bands): the ViewModel's (a rows band's scroll,
   * row zoom and folds, session state), a shown story frame's `bandTop` and `rowScale` merged over them, the session's
   * folds kept. Kept while both are the same objects, so the frame's key holds.
   */
  private _bandsFor(section: string, w: StepViewport): RenderViewport['bands'] | undefined {
    const own = this.view.bandsOf(section);
    const framed = this._shownFrame ? w.bands : undefined;
    if (!framed) return own;
    const kept = this._bandsMerged.get(section);
    if (kept && kept[0] === own && kept[1] === framed) return kept[2];
    const out: Record<string, BandViewport> = { ...own };
    for (const track of Object.keys(framed)) {
      const f = framed[track];
      const o = own?.[track];
      out[track] = {
        bandTop: f.bandTop ?? o?.bandTop ?? 0,
        ...((f.rowScale ?? o?.rowScale) !== undefined ? { rowScale: f.rowScale ?? o?.rowScale } : {}),
        ...(o?.folds !== undefined ? { folds: o.folds } : {}),
      };
    }
    this._bandsMerged.set(section, [own, framed, out]);
    return out;
  }

  private readonly _bandsMerged = new Map<
    string,
    [RenderViewport['bands'] | undefined, NonNullable<StepViewport['bands']>, RenderViewport['bands']]
  >();

  /** A section's placement in the last frame (its strips, what its layout's arithmetic reads), by address. */
  placementOf(section: string): Placement | undefined {
    const sections = this._frame?.layout.sections;
    if (!sections) return undefined;
    for (let i = 0; i < sections.length; i++) if (sections[i].section.address === section) return sections[i].placement;
    return undefined;
  }

  /** The viewport the last frame was drawn with, a copy (for the SVG and tests). */
  get viewport(): RenderViewport | undefined {
    const resolved = this.view.resolved.value;
    const v = resolved ? this._viewport(resolved, false) : undefined;
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
        layouts: this.layouts,
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
    // every leading section its own viewport, read by its layout (_spec/README.md "Composition"); the first's once
    const vp = this._viewport(resolved);
    const viewports = this._viewports(resolved);
    if (viewports) input.viewports = viewports;
    else delete input.viewports;
    if (vp) {
      this._drawnWindow.start = vp.start;
      this._drawnWindow.end = vp.end;
    }
    this._updateRails(resolved);
    if (this._rails) {
      input.rails = this._rails;
      input.railsActive = this._railsActive;
    } else {
      delete input.rails;
      delete input.railsActive;
    }
    // what each leading section showed: where the next motion starts
    this._drawnViews.clear();
    if (viewports) for (const k in viewports) this._drawnViews.set(k, viewports[k]);
    // a motion is not settled while it runs
    if (this._request) this._moved();
    input.dpr = dpr;
    input.settled = this._quiet && this.store.store.settled.value;
    input.emphasis = this._emphasis;
    input.cvdSafe = this.legends.cvdSafe.value;
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
      // what the models read of the placement (a text's lines in view) follows the frame
      this.view.placed();
      if (!this._restDue) this._changed();
      this._wanter?.want(frame.wants);
      this._issues.set(frame.layout.issues);
      // legends on settle only (a pan's frames derive none), compared structurally
      if (input.settled) this.legends.update(frame.layout);
      this.store.holdShown(frame.layout.shown);
    }
    this._frame = frame;
    // the layout placed the view (its width): a rail that went idle or live with it is drawn so now (the window pass
    // alone: which rails can move is outside the key)
    const wasActive = this._railsActive;
    this._updateRails(resolved);
    if (this._rails && this._railsActive !== wasActive) {
      input.railsActive = this._railsActive;
      frame = renderer.frame(input);
      this._frame = frame;
    }

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
    this._drawn.set(this._drawn.value + 1);
    if (!this._shownFrame) this._drawnRoot = this.spec.shown.value;
    if (this._shownFrame) this._drawnFrames.set(this._drawnFrames.value + 1);
    if (this._restDue) {
      // the settled picture is drawn: it rests now, once (a layout of this frame waits for no second rest)
      this._restDue = false;
      DomHost.clearTimeout(this._restTimer);
      this._restTimer = undefined;
      this._rested.set(this._rested.value + 1);
    }
    return frame;
  }

  /**
   * The rails this frame draws (_spec/README.md "Composition", Navigators): every linear section no shown navigator
   * binds, and every rows track no navigator binds along its rows (a rows rail), while `navigator.form` is `rail` and
   * the reader shows them; each active when the view can move it (the `navigator` behavior listed, not idle: idle is
   * the whole axis at the zoom limit, a follower under `ignore`; a rows rail's every row shown, the layout's).
   */
  private _updateRails(resolved: ResolvedView): void {
    const on = this.config.navigator.form === 'rail' && this._railsShown.value;
    const f = this._railsFor;
    if (f[0] !== resolved || f[1] !== on || f[2] !== this._hostedVersion) {
      f[0] = resolved;
      f[1] = on;
      f[2] = this._hostedVersion;
      // a host's navigator (`navigator(def)`) replaces the rail on what it binds, as a declared one does
      const hosted: ResolvedNavigator[] = [];
      for (const m of this._navigators.values()) {
        const n = m.isHosted ? m.navigator : undefined;
        if (n) hosted.push(n);
      }
      // and every rows layer's rows rail, by its track (a navigator over the section's rows replaces it)
      const rows = on ? rowsRailTracks(resolved, hosted) : [];
      this._rowsRails = new Set(rows);
      this._rails = on
        ? Object.fromEntries([...railSections(resolved, hosted), ...rows].map((a) => [a, true as const]))
        : undefined;
      this._railsActive = {};
    }
    const rails = this._rails;
    if (!rails) return;
    const behavior = this._handles.has(NavigatorProvider);
    let changed = false;
    // a rows rail can move while the behavior is there (idle while every row shows: the layout knows, `whole`)
    const can = (a: string) => behavior && (this._rowsRails.has(a) || !windowIdle(this.view, a, this.config));
    for (const a in rails) if (this._railsActive[a] !== can(a)) changed = true;
    if (!changed) return;
    const next: Record<string, boolean> = {};
    for (const a in rails) next[a] = can(a);
    this._railsActive = next;
  }

  /** Shows or hides the ruler rails (the reader's session: the embedded layout's corner button). */
  setRailsShown(shown: boolean): void {
    if (this._railsShown.value === shown) return;
    this._railsShown.set(shown);
    this.schedule();
  }

  /** The emphasis the last frame drew (the selection, the highlights, the cursor): what a navigator's marks read. */
  get emphasis(): Emphasis {
    return this._emphasis;
  }

  /**
   * A navigator's model by its id (a declared node of the view's, by id or address, or a host's): made once and kept
   * with the viewer. Undefined when there is none of that id.
   */
  navigatorModel(id: string): NavigatorModel | undefined {
    const kept = this._navigators.get(id);
    if (kept) return kept;
    const view = this.view.resolved.value;
    if (!view?.navigators.some((n) => n.id === id || n.address === id)) return undefined;
    const m = new NavigatorModel(this, id);
    this._navigators.set(id, m);
    this.registerPermanentDispose(() => m.retire());
    return m;
  }

  /**
   * A host's navigator (_spec/README.md "Composition", Navigators): a definition of the `navigator` node's shape, never
   * saved in the spec (a Mol* plugin's strip over a user's spec); mounted by its id like a declared one.
   */
  navigator(def: NavigatorDef): NavigatorModel {
    this._navigators.get(def.id)?.retire();
    const m = new NavigatorModel(this, def.id, def);
    this._navigators.set(def.id, m);
    this._hostedVersion++;
    this.registerPermanentDispose(() => m.retire());
    this.schedule();
    return m;
  }

  /** The view's renderer (its image cache and held domains: a navigator's picture shares them). */
  get renderer(): Renderer | undefined {
    return this._renderer;
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
    // a window's part under the pointer (a rail's, a ring's) lights it; a window being dragged grows its rail
    e.drag = this.input?.dragHit;
    const wh = this.input?.windowHit;
    if (wh) {
      e.hover = wh;
      this._hoverVersion = -1;
      return;
    }
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
    // a cell of two axes (a grid's, a rows cell): its second position, what the crosshair draws
    (h as { ySystem: string | undefined }).ySystem = r.ySystem;
    (h as { yPosition: number | undefined }).yPosition = Number.isNaN(r.yPosition) ? undefined : r.yPosition;
    // the clade hover of a rows layer reads it (_spec/README.md "Render", The rows layer, Hover)
    (h as { target: unknown }).target = r.target;
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
    for (const name of this.highlight.names.value) {
      const h = this.highlight.get(name);
      if (!h || isNothing(h)) continue;
      const location = this._projector.of(h.shown, { own: true }) ?? h.shown;
      out.push({ location, color: highlightColorOf(this, name), ...(h.style ? { style: h.style } : {}) });
    }
    this._highlights = out;
    this._emphasis.highlights = out.length ? out : undefined;
  }

  /** The selection's dynamic nodes held while it is selected (`selection.holdNodes`). */
  /**
   * A selection made elsewhere (another section, a link, a structure) scrolls each section that shows none of it to its
   * first part (_spec/README.md "Composition", `reveal.linked`): text sections by default, as a move by hand (no
   * step of the history). Not the section the selection is being made in (a drag's, the cursor's): it shows it already.
   */
  private _revealLinked(): void {
    const mode = this.config.reveal.linked;
    const resolved = this.view.resolved.value;
    const shown = this.selection.value.value.shown;
    if (mode === 'none' || !resolved || this._retired || !shown.parts.length) return;
    const own = this.input.selectingIn ?? this.selection.cursor.value?.section;
    let loc: Location | undefined;
    for (const s of resolved.sections) {
      if (s.address === own || (mode === 'text' && s.layout !== 'text')) continue;
      const v = this.view.get(s.address);
      if (!v || v.pending) continue;
      loc ??= this._projector.of(shown, { own: true }) ?? shown;
      const r = loc.parts.find((p) => p.system === s.axis.system)?.ranges;
      if (!r?.length) continue;
      let first = Number.POSITIVE_INFINITY;
      let shownNow = false;
      for (let i = 0; i + 1 < r.length && !shownNow; i += 2) {
        const a = nearestOffsetOf(s.axis, r[i]);
        const b = nearestOffsetOf(s.axis, r[i + 1]);
        if (a === undefined || b === undefined) continue;
        const lo = Math.min(a, b);
        const hi = Math.max(a, b);
        if (hi > v.start && lo < v.end) shownNow = true;
        first = Math.min(first, lo);
      }
      if (!shownNow && Number.isFinite(first)) this.view.reveal(s.address, first, first + 1);
    }
  }

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
  hit(x: number, y: number, touch = false): RenderHit | undefined {
    return this._frame ? hitView(this._frame.layout, x, y, touch) : undefined;
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
    const viewports = resolved ? this._viewports(resolved, false) : undefined;
    // a text section exports every line (the figure): the picture is as tall as it, not the canvas's scrolled box
    const figure = resolved?.sections.some((s) => sectionLayoutName(s, resolved) === 'text') ?? false;
    return this.context.renderSvg(this.spec.shown.value, {
      size: { width: this.width, height: figure ? 'auto' : this.height },
      ...(resolved ? { view: resolved.address } : {}),
      ...(viewports ? { viewports } : {}),
      store: this.store.store,
      dpr: this.dpr,
      tokens: this.tokens,
      // the navigators as the parts draw them (an included one's zoom strip, its direction)
      navigator: { zoomStrip: this.config.navigator.zoomStrip, direction: this.config.navigator.direction },
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
