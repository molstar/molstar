/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Listenable, ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Issue, Node } from '@molstar/core/repr/tree';
import {
  composeSeqContext,
  type SeqContext,
  type SeqExtension,
  type SeqRegistryEntry,
} from '@molstar/seqx-core/context';
import {
  type Emphasis,
  type Frame,
  type Hit as RenderHit,
  type Placement,
  type Renderer,
  type RendererOptions,
  type RenderSvgOptions,
  type RenderSvgResult,
  type Viewport as RenderViewport,
  hitView,
  textMinScale,
} from '@molstar/seqx-core/render';
import { labelColumnWidth } from '@molstar/seqx-core/render/tracks';
import { locationOfResult } from '@molstar/seqx-core/query';
import { type Lengths, sourceLengths } from '@molstar/seqx-core/resolve';
import { mapsEqual } from '@molstar/core/util/map';
import type { SeqTokens } from '@molstar/seqx-core/render';
import type { StepViewport, TransitionPlan } from '@molstar/seqx-core/story';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DomHost } from '@molstar/ui/dom/host';
import { bindingTable, DefaultActions, type BindingTable, BindingPresets } from '@molstar/seqx-viewer/bindings';
import type { CanvasBackend } from '@molstar/seqx-viewer/canvas';
import { type ViewerConfig, viewerConfigOf } from '@molstar/seqx-viewer/config';
import {
  AddressModel,
  HighlightModel,
  HoverModel,
  SelectionModel,
  SpecModel,
  sameIssues,
  StoreModel,
  ViewModel,
} from '@molstar/seqx-viewer/models';
import type { AddressingBehavior, LinkRun } from '@molstar/seqx-viewer/behaviors';
import type { SeqAnimationParticipant } from '@molstar/seqx-viewer/ports';
import { frameAnimation } from './animation.js';
import { ViewerInput } from './input.js';
import { LegendModel } from './legends.js';
import { type HiddenTrackEntry, hiddenTracksState } from './hidden-tracks.js';
import { modelIssuesState } from './model-issues.js';
import type {
  BehaviorContext,
  BehaviorHandle,
  BehaviorProvider,
  SeqViewerExtension,
  SeqViewerSpec,
} from './viewer-spec.js';
import type { HeaderGeometry, SeqViewerRegistryEntry } from './header-entries.js';
import type { ViewerExtensions } from './viewer-extend.js';
import { mountCanvas } from './mount.js';
import type { NavigatorDef, NavigatorModel } from './navigator.js';
import { behaviorOf, createBehaviors, listedActions } from './viewer-behaviors.js';
import { FrameClock } from './viewer-clock.js';
import { ViewerEmphasis } from './viewer-emphasis.js';
import { ViewerFrame } from './viewer-frame.js';
import { ViewerViewports } from './viewer-viewports.js';
import { ViewerMotion } from './viewer-motion.js';
import { ViewerRails } from './viewer-rails.js';
import { holdSelection, revealLinked } from './viewer-selection.js';

// The viewer (_spec/README.md "Composition"): composed from a spec by `createSeqViewer`, as a Mol* plugin is: the
// context (registry, configuration, ports) first, then the models, the bindings and the behaviors, then the UI. Frames
// and time are `DomHost`'s. This module composes and wires; each concern is a collaborator beside it, given what it
// reads: the behaviors (`viewer-behaviors.ts`), the frame clock (`viewer-clock.ts`), motion, the frame slot
// and the viewports a frame reads as it shows them (`viewer-motion.ts`), the rails and navigators
// (`viewer-rails.ts`), the emphasis (`viewer-emphasis.ts`), what a selection moves (`viewer-selection.ts`) and the
// frame itself (`viewer-frame.ts`): the renderer lays the view out (or, for a hover, redraws the emphasis passes
// alone) and the Canvas2D backend draws what changed.

const NeverAborted = new AbortController().signal;

/** A command a behavior provides, with the context it runs in. */
export type ViewerCommand = CommandDefinition<any>;

export class SeqViewer extends ReactiveModel {
  readonly context: SeqContext;
  /** The spec's registry as listed: the core's entries and the viewer's own kinds (header entries). */
  readonly registry: readonly SeqViewerRegistryEntry[];
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
  /**
   * Aborted when the viewer retires (_spec/README.md "Composition", Retire), before its models are: what a helper over
   * the viewer (`linkStructure`) ends with.
   */
  readonly retiredSignal: AbortSignal;
  /** The binding table: the spec's, and a layer for each extension's bindings (`extend` makes it again). */
  bindings: BindingTable;
  readonly input: ViewerInput;
  /** The behaviors' handles by provider name. */
  readonly behaviors = new Map<string, any>();
  /** The commands the behaviors provide. */
  readonly commands: ViewerCommand[] = [];
  /** The layout's issues of the last frame (a color it cannot read, a layer it cannot draw). */
  readonly layoutIssues: ReadonlyState<readonly Issue[]>;
  /** The model's issues as the viewer says them: one of an address still pending held back (`model-issues.ts`). */
  readonly modelIssues: ReadonlyState<readonly Issue[]>;
  /**
   * Each section's hidden tracks from the spec, in its order (address, label, section); its change is the event. The
   * view draws nothing for them: the application lists and restores them (`mountHiddenTracks`, `track.show`).
   */
  readonly hiddenTracks: ReadonlyState<readonly HiddenTrackEntry[]>;
  /** The viewer's own: a gesture that needs a behavior the spec does not list, once each. */
  readonly behaviorIssues: ReadonlyState<readonly Issue[]>;
  /**
   * Moves once each time the picture comes to rest: when the frame laid out settled is drawn, `canvas.settleMs` after the
   * view's last move; or `canvas.settleMs` after the last layout that changed what is drawn while the view was still
   * (data arriving, the content's height): what follows the picture waits for it.
   */
  readonly rested: ReadonlyState<number>;
  /**
   * Moves (a new record each time) after each command has run, whether or not it changed anything: what says what a
   * key did when the state it moves stays (the announcements: a cursor key at the end of its axis).
   */
  readonly commandRan: ReadonlyState<{ readonly id: string; readonly n: number }>;

  /** Moves to the token table's version each time the tokens change (`setTokens`): what follows the picture redraws. */
  readonly tokensVersion: ReadonlyState<number>;
  /** The token table the frames and the SVG are drawn with: the context's, or the root's CSS when mounted. */
  tokens: SeqTokens;

  /** The canvas's size in CSS px (the data's width and the content's height by default). */
  width = 0;
  height = 0;
  /** The device pixel ratio of the last frame. */
  dpr = 1;
  /**
   * The section layouts by name, the spec's registry's (and its extensions'): one map the frames draw with and the
   * view's moves navigate by (a section naming another is laid out as tracks, and said).
   */
  get layouts(): RendererOptions['layouts'] {
    return this.context.layouts;
  }
  /** Frames drawn and layouts made: what a test or a profile reads. */
  readonly counts = { frames: 0, layouts: 0 };
  /** Bumped after every frame drawn: a navigator redraws over it. */
  readonly drawn: ReadonlyState<number>;
  /**
   * The ruler rails shown (_spec/README.md "Composition", Navigators): the reader's session (the embedded layout's
   * corner button), never the spec; `navigator.form: none` draws none whatever it says.
   */
  readonly railsShown: ReadonlyState<boolean>;
  /** The participant face of the animation port: what the application's clock drives. */
  readonly animation: SeqAnimationParticipant;
  /**
   * The platform asks for reduced motion (the mounted UI's media query): under `view.motion: platform` every jump and
   * story step is a cut.
   */
  reducedMotion = false;

  /** Moves when the issues do, not with each layout: what follows them (a story's SVG export) runs once a change. */
  private readonly _issues = new State<readonly Issue[]>([], sameIssues);
  private readonly _behaviorIssues = new State<readonly Issue[]>([], sameIssues);
  private readonly _ran = new State<{ readonly id: string; readonly n: number }>({ id: '', n: 0 });
  private readonly _tokensVersion = new State(0);
  /**
   * The collaborators, each a module beside this one given what it reads (the frame calls them without allocating):
   * one frame at a time and the rest, the motion, the frame slot and the viewports a frame reads, the rails and
   * navigators, the emphasis, the frame.
   */
  private readonly _clock: FrameClock;
  private readonly _motion: ViewerMotion;
  private readonly _rails: ViewerRails;
  private readonly _emphasis: ViewerEmphasis;
  private readonly _frame: ViewerFrame;
  private readonly _viewports: ViewerViewports;
  private _retired = false;
  /**
   * What `extend` (its chunk) reads: the extensions it registered with their counts, the bindings made with them, the
   * behaviors' context.
   * @internal
   */
  readonly _ext: ViewerExtensions;
  private _syncing: Promise<void> = Promise.resolve();
  private _viewportListeners: (() => void)[] = [];

  constructor(viewerSpec: SeqViewerSpec, target?: HTMLElement, fixedDpr?: number) {
    super();
    // (what the composition calls often, named once)
    const listen = <T>(source: Listenable<T>, fn: (value: T) => void) => this.listenPermanently(source, fn);
    const onRetire = (dispose: () => void) => this.registerPermanentDispose(dispose);
    const ended = new AbortController();
    this.retiredSignal = ended.signal;
    const schedule = () => this.schedule();
    const config = (this.config = viewerConfigOf(viewerSpec.config));
    // the headless pieces first: the registry, the configuration and the ports
    // the registry as listed: what the canvas and the SVG draw (_spec/README.md "Composition", Injected, opt in)
    // extensions (_spec/README.md "Composition", Extensions): the core's checks and entries the context's, the viewer's
    // kinds after the spec's own
    const exts = (viewerSpec.extensions ?? []) as readonly SeqViewerExtension[];
    // each one's viewer parts added by its maker (`defineViewerExtension`; `defineExtension` refuses viewer parts)
    for (const e of exts) viewerSpec = e.composeViewer?.(viewerSpec) ?? viewerSpec;
    this.registry = viewerSpec.registry ?? [];
    const context = (this.context = composeSeqContext({
      // the core's kinds (the viewer's own, the header entries, are its behaviors')
      ...(viewerSpec.registry
        ? { registry: viewerSpec.registry.filter((e): e is SeqRegistryEntry => e.kind !== 'header-action') }
        : {}),
      extensions: exts,
      ...(viewerSpec.config ? { config: viewerSpec.config } : {}),
      ...(viewerSpec.ports ? { ports: viewerSpec.ports } : {}),
    }));
    onRetire(() => context.retire());
    const clock = (this._clock = new FrameClock(this));
    const rails = (this._rails = new ViewerRails(this, viewerSpec.navigators, this._behaviorIssues, onRetire));
    this.layoutIssues = this._issues.asReadonly;
    this.behaviorIssues = this._behaviorIssues.asReadonly;
    this.rested = clock.rested;
    this.commandRan = this._ran.asReadonly;
    this.tokensVersion = this._tokensVersion.asReadonly;
    this.railsShown = rails.shown;
    this.tokens = context.tokens;

    const spec = (this.spec = new SpecModel(
      createSeqSpec().build(),
      { runtimes: context.runtimes, config: context.resolveConfig },
      config,
      false,
      this._lengths.asReadonly,
      viewerSpec.story ? { story: viewerSpec.story } : {},
    ));
    const view = (this.view = new ViewModel(spec, config, context.kinds));
    // a section's moves are its layout's arithmetic over its placement in the last frame (_spec/README.md "ViewModel")
    view.setLayoutSource({
      placement: (section) => this.placementOf(section),
      // the context's, as `extend` leaves them
      get layouts() {
        return context.layouts;
      },
      minScale: textMinScale(context.renderConfig),
    });
    const hover = (this.hover = new HoverModel(spec, view, config));
    const store = (this.store = new StoreModel(config, { removed: (gone) => this.selection.dropGone(gone) }));
    this.modelIssues = modelIssuesState(this);
    this.hiddenTracks = hiddenTracksState(this);
    // an element selected by an address and then toggled is drawn by its extent, resolved over what is shown
    const selection = (this.selection = new SelectionModel(
      (els) => locationOfResult(spec.shown.value, this.addresses.queryModel(), { type: 'keys', keys: els }).location,
      { extendScope: config.selection.extendScope, acrossBreaks: config.selection.acrossBreaks },
    ));
    const emphasis = (this._emphasis = new ViewerEmphasis(this));
    const highlight = (this.highlight = new HighlightModel());
    const legends = (this.legends = new LegendModel(config.theme.cvdSafe));
    const addresses = (this.addresses = new AddressModel(
      spec,
      view,
      store,
      hover,
      selection,
      highlight,
      config,
      context.renderConfig.maxNumberingPositions,
    ));
    for (const m of [addresses, highlight, selection, store, hover, view, spec]) onRetire(() => m.retire());
    // the host's history, through Seq*'s port: the spec takes part in it (none headless unless the spec gives one)
    const history = viewerSpec.history ?? (target && viewerSpec.ui?.history ? viewerSpec.ui.history() : undefined);
    if (history) spec.attachHistory(history);

    const bind = (vs: SeqViewerSpec) =>
      bindingTable(vs.bindings?.preset ?? BindingPresets[config.bindings.preset], vs.bindings?.overrides ?? [], {
        actions: [...DefaultActions, ...listedActions(vs.behaviors), ...(vs.bindings?.actions ?? [])],
        keys: config.view.keys,
      });
    const bindings = (this.bindings = bind(viewerSpec));

    // behaviors, after the models; one of each name
    const behaviorContext: BehaviorContext = {
      spec,
      view,
      hover,
      store,
      selection,
      highlight,
      addresses,
      config,
      bindings,
      legends,
      layout: () => this._frame?.frame?.layout,
    };
    this._ext = {
      extended: new Map(),
      spec: viewerSpec,
      bind,
      context: behaviorContext,
      onRetire,
      create: createBehaviors,
    };
    try {
      createBehaviors(viewerSpec.behaviors, behaviorContext, this.behaviors, this.commands, onRetire);
      this.input = new ViewerInput(this);
    } catch (e) {
      // what was made is not left behind
      this.retire();
      throw e;
    }

    // what moves the picture asks for a frame
    listen(spec.model, (model) => {
      // what the model no longer opens goes at once: no frame before the sync reads an edited source's old file
      store.detachStale(model);
      this._sync();
    });
    listen(hover.version, schedule);
    listen(selection.value, () => {
      holdSelection(this);
      revealLinked(this, emphasis.projector);
      this.schedule();
    });
    listen(selection.cursor, schedule);
    listen(highlight.version, schedule);
    // the reader's substitution redraws (its frames are keyed by it); the hover moves the legend's entry
    listen(legends.cvdSafe, schedule);
    // a band's own scroll, row zoom and folds (a rows band's)
    listen(view.bands, schedule);
    listen(hover.version, () => legends.updateHover(this._frame?.frame?.layout, hover.current));
    listen(store.store.version, () => {
      this._readLengths();
      this.schedule();
    });
    listen(store.store.settled, schedule);
    // the animation port (_spec/README.md "Ports", Animation): the application's clock, the page's frames by default
    const motion = (this._motion = new ViewerMotion(
      this,
      clock,
      viewerSpec.story,
      viewerSpec.animation ?? frameAnimation(),
    ));
    const viewports = (this._viewports = new ViewerViewports(view, motion));
    this.animation = motion.participant;
    onRetire(() => motion.dispose());
    const frame = (this._frame = new ViewerFrame(
      this,
      fixedDpr,
      clock,
      motion,
      viewports,
      rails,
      emphasis,
      this._issues,
    ));
    this.drawn = frame.drawn;
    // a jump animates from the windows drawn to the model's (the model is there at once: hits read the frame drawn)
    listen(view.jumps, (n) => motion.jumped(n));
    // an edit, an undo or a load under a frame that changes the spec ends it: the frame was made from the old root
    listen(spec.head, () => motion.headMoved());
    // a story step entered by a command or a host is a transition too (its view is a snapshot's: a jump)
    listen(spec.step, () => motion.stepChanged());
    let resolvedFor = spec.model.value;
    listen(view.resolved, () => {
      this._followViewports();
      this._applyWidths();
      // another view of the same model: what it asks the model to build (a ranges layer's links' ends, core's _spec
      // "Render", Layer types: ranges) is opened by a sync, as the model listener's is for a new model
      if (spec.model.value === resolvedFor) void this._sync();
      resolvedFor = spec.model.value;
      this.schedule();
    });
    this._followViewports();
    // what reads the commands of every behavior and the frame (the header entries' checks), before the UI
    try {
      for (const h of this.behaviors.values() as Iterable<BehaviorHandle | undefined>) h?.composed?.(this);
    } catch (e) {
      this.retire();
      throw e;
    }
    const ui = viewerSpec.ui;
    if (target && ui?.toolbar) onRetire(ui.toolbar(this, target));
    if (target && ui?.canvas) onRetire(mountCanvas(this, target));
    if (target && ui?.navigators) onRetire(ui.navigators(this, target));
    if (target && ui?.railToggle) onRetire(ui.railToggle(this, target));
    if (target && ui?.legends) onRetire(ui.legends(this, target));
    // registered last, so it runs first: everything after it (the UI's teardown, the behaviors, the models) is disposed
    // while the viewer is retired, and none of it asks for a frame
    onRetire(() => {
      this._retired = true;
      // a helper over the viewer ends first, while the models it clears still stand
      ended.abort();
      clock.dispose();
      frame.dispose();
      for (const u of this._viewportListeners) u();
      this._viewportListeners.length = 0;
    });
  }

  /**
   * The handle of a behavior by its provider or its name (what the input calls). One the spec does not list, or one of a
   * built-in's name whose handle lacks a method the viewer calls (`BehaviorShapes`: a host's provider under that name
   * with another shape), is reported once as an issue (`behaviorIssues`) and the gesture does nothing.
   */
  behavior<H>(provider: BehaviorProvider<any, H> | string): H | undefined {
    // by name: what the input asks for names no provider object, so a viewer bundles only the behaviors it lists
    return behaviorOf<H>(this.behaviors, typeof provider === 'string' ? provider : provider.name, this._behaviorIssues);
  }

  /**
   * The headers' geometry (_spec/README.md "Behaviors", Headers): each shown header's boxes in the canvas's CSS px,
   * published when a layout moves them; undefined without the `headers` behavior.
   */
  get headers(): ReadonlyState<readonly HeaderGeometry[]> | undefined {
    return this.behaviors.get('headers')?.geometry;
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
    // an attach moves no entry: the lengths the opened sources give are read once they are in
    this._syncing = this._syncing.then(() => this._readLengths());
    this._syncing = this._syncing.then(() => this.schedule());
    return this._syncing;
  }

  /**
   * A source-backed system's length arrives with its source's coverage: the model is resolved again with it
   * (_spec/README.md "SpecModel"); an equal set of lengths is no change.
   */
  private _readLengths(): void {
    if (this._retired) return;
    // a coverage only from the source the model names as it is now (an edited source's old file gives none)
    const model = this.spec.model.value;
    this._lengths.set(sourceLengths(model, (s) => this.store.openedCoverage(model, s)));
  }

  /** The source-backed systems' lengths the spec's model is resolved with (`_readLengths`). */
  private readonly _lengths = new State<Lengths>(new Map(), mapsEqual);

  /**
   * Resolves once the picture is complete: the sources open, the behaviors' lazy chunks landed, the data the view
   * wants arrived, the view settled. A viewer with no size yet is given the default one. For tests, exports and
   * stories that wait for the picture.
   */
  async settle(timeoutMs = 10_000): Promise<boolean> {
    await this._syncing;
    for (const h of this.behaviors.values() as Iterable<BehaviorHandle | undefined>) await h?.loaded?.();
    const until = DomHost.now() + timeoutMs;
    this._clock.quiet = true;
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
    this._frame.useBackend(backend);
    this.schedule();
  }

  /**
   * Sets the token table (a theme switch). A table equal to the current one, by its version, is no change; else the
   * renderer keeps its layout key's token version, the backend redraws and a frame is asked for.
   */
  setTokens(tokens: SeqTokens): boolean {
    if (this._retired || tokens.version === this.tokens.version) return false;
    this.tokens = tokens;
    this._frame.setTokens(tokens);
    this._tokensVersion.set(tokens.version);
    this.schedule();
    return true;
  }

  /** Drops what was measured: a web font arrived. */
  invalidateText(): void {
    this.context.text.clear();
    this._column[0] = undefined;
    // the label column measured again: the data box's width with it
    this._applyWidths();
    this.schedule();
  }

  /** The data box's left edge: the page's inset, the label column (measured over the view's labels) and its gap. */
  get dataX(): number {
    const c = this.context.renderConfig;
    const view = this.view.resolved.value;
    if (this._column[0] !== view) {
      this._column[0] = view;
      this._column[1] = view ? labelColumnWidth(view, c, this.context.text) : c.labelWidth;
    }
    return c.inset + (this._column[1] as number) + c.labelGap;
  }
  /** The label column's width and the resolved view it was measured for. */
  private readonly _column: [unknown, number] = [undefined, 0];

  private _applyWidths(): void {
    const sections = this.view.resolved.value?.sections;
    if (!sections || this.width <= 0) return;
    const dataW = Math.max(1, this.width - this.dataX - this.context.renderConfig.inset);
    for (const s of sections) this.view.setWidth(s.address, dataW);
  }

  /** A frame is asked for by each section's viewport changing: listened to, again when the view's sections change. */
  private _followViewports(): void {
    for (const u of this._viewportListeners) u();
    this._viewportListeners.length = 0;
    for (const s of this.view.resolved.value?.sections ?? [])
      this._viewportListeners.push(
        this.view.viewport(s.address).listen(() => {
          // a move by hand ends a motion, a story's too, and clears a frame a host shows
          this._motion.handMoved();
          this._clock.moved();
          this.schedule();
        }),
      );
  }

  /** Asks for a frame: one a frame (`@molstar/ui/dom`'s frame callbacks, on `DomHost`), whatever asked; `retire` cancels it. */
  schedule(): void {
    // a re-hit's own hover change is drawn by the frame that made it
    this._clock.schedule();
  }

  /**
   * Stops a motion (a jump's, a story step's) or a host's frame where it is drawn (_spec/README.md "Composition",
   * Motion): each leading section takes the viewport drawn now (a move, no step of the history), the spec shown is the
   * state again, and a frame is drawn, so the input that stopped it reads one picture: what it hits is what the model
   * shows. What every pointer and key input does first.
   */
  stopMotion(): void {
    this._motion.stop();
  }

  /**
   * Enters a story step, or leaves the step shown (`undefined`), as a transition (_spec/README.md "Composition",
   * Stories): navigation, no commit, no history step. The step's state is shown, its view's sections placed (positions
   * to offsets, exactly), and a plan from what was drawn to it is made, timed by the step. `run` (true): the plan runs
   * on the animation port; false: its start is shown and the application drives it (`animation.progress`, then
   * `animation.clear()`). A cut under reduced motion, and with no `story` part listed (the step is not entered then: the
   * issue says so). The issue when the step is unknown.
   */
  playStep(id: string | undefined, options: { readonly run?: boolean } = {}): { plan?: TransitionPlan; issue?: Issue } {
    return this._motion.playStep(id, options);
  }

  /** Whether a motion (a jump's, a story step's) is drawing. */
  get moving(): boolean {
    return this._motion.moving;
  }

  /** The window the last frame drew (during a motion, between where it started and the model's), a copy. */
  get drawnViewport(): RenderViewport | undefined {
    const w = this._motion.drawnWindow;
    return Number.isNaN(w.start) ? undefined : { start: w.start, end: w.end };
  }

  /** Each leading section's viewport the last frame drew, by address (during a motion, the frame's), copies. */
  get drawnViewports(): Readonly<Record<string, StepViewport>> {
    return Object.fromEntries([...this._motion.drawnViews].map(([k, v]) => [k, { ...v }]));
  }

  /** A section's placement in the last frame (its strips, what its layout's arithmetic reads), by address. */
  placementOf(section: string): Placement | undefined {
    return this._frame?.placementOf(section);
  }

  /**
   * The window the renderer draws (`ViewerViewports.window`): the first leading section's, as the frame slot shows it
   * during a motion; a frame reads it once.
   */
  private _viewport(resolved: ResolvedView, eye = true): RenderViewport | undefined {
    return this._viewports.window(resolved, eye);
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
    return this._frame.draw();
  }

  /** Shows or hides the ruler rails (the reader's session: the embedded layout's corner button). */
  setRailsShown(shown: boolean): void {
    this._rails.setShown(shown);
  }

  /** The emphasis the last frame drew (the selection, the highlights, the cursor): what a navigator's marks read. */
  get emphasis(): Emphasis {
    return this._emphasis.emphasis;
  }

  /**
   * A navigator's model by its id (a declared node of the view's, by id or address, or a host's): made once and kept
   * with the viewer. Undefined when there is none of that id, or the spec lists no navigators (`navigators`).
   */
  navigatorModel(id: string): NavigatorModel | undefined {
    return this._rails.model(id);
  }

  /**
   * A host's navigator (_spec/README.md "Composition", Navigators): a definition of the `navigator` node's shape, never
   * saved in the spec (a Mol* plugin's strip over a user's spec); mounted by its id like a declared one. Throws when the
   * spec lists no navigators (`navigators`): a host that shows one composes them.
   */
  navigator(def: NavigatorDef): NavigatorModel {
    return this._rails.hosted(def);
  }

  /** The view's renderer (its image cache and held domains: a navigator's picture shares them). */
  get renderer(): Renderer | undefined {
    return this._renderer;
  }

  /** The renderer the frame made (made again when the model or the view changes). */
  private get _renderer(): Renderer | undefined {
    return this._frame.renderer;
  }

  /** The last frame: its layout is what hit-testing reads. */
  get frame(): Frame | undefined {
    return this._frame.frame;
  }

  /** What is under a point of the last frame's layout. */
  hit(x: number, y: number, touch = false): RenderHit | undefined {
    const frame = this._frame.frame;
    return frame ? hitView(frame.layout, x, y, touch, this._emphasis.hoverHit()) : undefined;
  }

  // ---- extensions ----

  /**
   * Registers an extension at run time (_spec/README.md "Composition", Extensions) and resolves to its undo: the
   * context's checks and entries (nothing changed on an error), then its behaviors, commands and bindings; the view is
   * resolved again with the new layer kinds, so a layer drawn as its fallback is drawn as its type, and as the fallback
   * again after the undo (behaviors disposed, entries removed). The same descriptor again is counted. A lazy chunk:
   * registering late is rare, and a composition that never does bundles none of it.
   */
  async extend(extension: SeqExtension | SeqViewerExtension): Promise<() => void> {
    const { extendViewer } = await import('./viewer-extend.js');
    return extendViewer(this, extension as SeqViewerExtension, this._ext);
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
    this._ran.set({ id, n: this._ran.value.n + 1 });
    return r instanceof Promise ? r.then((x: CommandResult) => x.message ?? '') : (r.message ?? '');
  }

  /**
   * Runs a link's actions (`!seq-focus=/refs/dbd`), each through the action its command runs, in order: undefined when
   * the text is no link or the viewer has no `addressing` behavior (reported in `behaviorIssues`).
   */
  runLink(text: string): LinkRun | undefined {
    if (this._retired) return undefined;
    return this.behavior<AddressingBehavior>('addressing')?.runLink(text);
  }

  // ---- export ----

  /**
   * The picture as an SVG string: the same view, size and window through the context's ports and this viewer's store
   * (data already loaded is shared), at the device pixel ratio of the last frame: the density steps are in device
   * pixels, so this is the picture the canvas shows (`dpr: 1` for a device independent export).
   */
  renderSvg(options: Partial<RenderSvgOptions> = {}): Promise<RenderSvgResult> {
    const resolved = this.view.resolved.value;
    const viewports = resolved ? this._viewports.all(resolved, false) : undefined;
    // a section that scrolls inside its box (a text's, as its layout placed it) exports every line (the figure): the
    // picture is as tall as it, not the canvas's scrolled box; a text laid out as tracks (no text layout listed) is not
    const figure = resolved?.sections.some((s) => this.placementOf(s.address)?.clip === true) ?? false;
    const navigators = this._rails.part;
    return this.context.renderSvg(this.spec.shown.value, {
      size: { width: this.width, height: figure ? 'auto' : this.height },
      ...(resolved ? { view: resolved.address } : {}),
      ...(viewports ? { viewports } : {}),
      store: this.store.store,
      dpr: this.dpr,
      tokens: this.tokens,
      // the navigators as the parts draw them (an included one's zoom strip, its direction)
      navigator: { zoomStrip: this.config.navigator.zoomStrip, direction: this.config.navigator.direction },
      // the navigators part: a viewer that lists none exports none (and says so when asked)
      ...(navigators ? { navigatorExport: navigators.svgExport } : {}),
      // what the canvas shows of the session: the selection, as its address when it has one (the export resolves it
      // over what it holds, and projects it onto the linked sections), else as its Location; the highlight channels
      ...this._emphasis.svgEmphasis(),
      // the headers' states as the canvas draws them (a figure lays out no slots)
      ...(this.behaviors.get('headers') ? { headers: this.behaviors.get('headers').partsOf } : {}),
      ...options,
    });
  }
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
