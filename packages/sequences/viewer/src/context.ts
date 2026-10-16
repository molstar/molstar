/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headless context (seqspec:viewer#the-context, seqspec:viewer#data-mode): all state of one viewer. VC1 builds
// config, the registry, events, the scheduler, the graph, the document, the history and the views; later phases add the
// managers in their order, each subscribing through the context's guarded fan-out (`Fanout`).

import { issueOf, type Loaders, type Patch, type SeqViewSpecDoc, validate } from '@molstar/seq-core/format';
import type { Issue, KeyResolver } from '@molstar/seq-core/model';
import { type ActionManager, createActions, nodePath } from './actions.js';
import { type BehaviorManager, Behaviors } from './behaviors.js';
import { type CommandManager, Commands, TargetRangeError } from './commands.js';
import { type ConfigManager, createConfigManager, SeqViewerConfig } from './config.js';
import { dataKeys, dataMode, dataViews, type SeqDataInput, type SetDataOptions } from './data.js';
import {
  type DocumentManager,
  type DocumentPort,
  Documents,
  type HistoryManager,
  type Routed,
  type ViewManager,
  Views,
} from './document.js';
import { Focus, type FocusManager } from './focus.js';
import { type Assembled, Graph, type SeqGraph, type Skin } from './graph.js';
import { cutHeaderText } from './header-text.js';
import { Highlight, type HighlightManager } from './highlight.js';
import { Hover, type HoverManager, type PointerState } from './hover.js';
import { InputDispatcher, type SeqInput } from './input.js';
import { Interaction, type SeqInteraction } from './interaction.js';
import { type LabelManager, Labels } from './labels.js';
import { DefaultRenderProps, type LayoutManager, Layouts, mergeProps, type RenderProps } from './layout.js';
import type { SeqLociAny } from './lineage.js';
import { SetChannels } from './location-sets.js';
import { snapshotHints } from './prefetch.js';
import { type Coordinates, isCoordinates } from './projection.js';
import { type RevealOptions, revealTarget } from './reveal.js';
import { type RefManager, Refs } from './refs.js';
import { type SeqViewerRegistry, type SeqViewerRegistryEntry, ViewerRegistry } from './registry.js';
import { Results, type ResultsManager } from './results.js';
import { manualScheduler, type Scheduler, type TrackedScheduler, trackedScheduler } from './scheduler.js';
import { type SelectionManager, Selections } from './selection.js';
import { type SnapshotManager, Snapshots } from './snapshots.js';
import type { SeqViewerSpec } from './spec.js';
import {
  event,
  type Live,
  live,
  type Origin,
  type SeqEvent,
  type SeqEventEmitter,
  type SeqState,
  state,
} from './state.js';
import { type StreamManager, Streams } from './streams.js';
import { ViewerSystems } from './systems.js';
import { monospaceMeasurer, TextManager, type TextManagerApi, type TextMeasurer } from './text.js';
import { type IssueManager, Issues, type TaskManager, Tasks, Ui, type UiManager } from './ui.js';
import { type ScrollGroupManager, ScrollGroups, sectionsOf, type ViewportManager, Viewports } from './viewport.js';

export interface LogEntry {
  readonly level: 'info' | 'warning' | 'error';
  readonly message: string;
  readonly data?: unknown;
}
export type SeqIssue = Issue & { readonly ref: string | null };

export interface SeqViewerPorts {
  /** The format's I/O port; default: refuses every request (`svs.source.load`). */
  loaders: Loaders;
  /** Who commits proposed patches; default: the context itself, with history. */
  document: DocumentPort;
  /** Frames and time; default: manual (the root touches no DOM global; a browser host passes `rafScheduler(…)`). */
  scheduler: Scheduler;
  /** Canvas widths (seqspec:viewer#text-measuring); default: a monospace estimate. */
  text: TextMeasurer;
  log?: (entry: LogEntry) => void;
}
export interface GraphResult {
  readonly ok: boolean;
  readonly issues: readonly Issue[];
}
export interface PatchEvent {
  readonly patch: Patch;
  readonly label: string;
  readonly routed: Routed;
  readonly origin: Origin;
}

/** What is not input (input is `interaction`'s streams, seqspec:interaction#streams). */
export interface SeqViewerEvents {
  /** A new issue from a node, the engine, a registration, a behavior or a command. */
  readonly issue: SeqEvent<SeqIssue>;
  readonly log: SeqEvent<LogEntry>;
  /** Every proposal (seqspec:viewer#the-document). */
  readonly patch: SeqEvent<PatchEvent>;
}

/** The light skin, until the engine's skin reader or a host sets one. */
export const DefaultSkin: Skin = {
  dark: false,
  surface: 'rgb(255 255 255)',
  surfaceAlt: 'rgb(246 247 249)',
  ink: 'rgb(27 31 36)',
  inkQuiet: 'rgb(91 99 112)',
  inkFaint: 'rgb(196 200 207)',
  grid: 'rgb(229 231 235)',
  accent: 'rgb(181 69 44)',
  hover: 'rgb(43 127 255 / 0.14)',
  highlight: 'rgb(242 180 0 / 0.22)',
  selection: 'rgb(43 127 255 / 0.1)',
  selectionBorder: 'rgb(43 127 255 / 0.6)',
  unobserved: 'rgb(214 217 222)',
  pending: 'rgb(226 222 214)',
  grip: 'rgb(124 131 142)',
  gripHover: 'rgb(27 31 36)',
  window: 'rgb(91 99 112 / 0.1)',
  conflict: 'rgb(192 38 45)',
  focusRing: 'rgb(43 127 255)',
  fontMono: 'ui-monospace, SFMono-Regular, Menlo, monospace',
  fontUi: 'ui-sans-serif, system-ui, sans-serif',
};
/** The dark skin where a token is not set: skin.css's dark block (the skin reader's fallback; the light one is `DefaultSkin`). */
export const DarkSkin: Skin = {
  ...DefaultSkin,
  dark: true,
  surface: 'rgb(24 26 30)',
  surfaceAlt: 'rgb(34 37 43)',
  ink: 'rgb(230 232 235)',
  inkQuiet: 'rgb(156 163 175)',
  inkFaint: 'rgb(75 80 90)',
  grid: 'rgb(48 52 60)',
  accent: 'rgb(240 128 96)',
  hover: 'rgb(110 160 255 / 0.2)',
  highlight: 'rgb(250 204 21 / 0.28)',
  selection: 'rgb(110 160 255 / 0.14)',
  selectionBorder: 'rgb(110 160 255 / 0.7)',
  unobserved: 'rgb(60 64 72)',
  pending: 'rgb(56 53 48)',
  grip: 'rgb(140 147 158)',
  gripHover: 'rgb(240 242 245)',
  window: 'rgb(156 163 175 / 0.12)',
  conflict: 'rgb(255 123 114)',
  focusRing: 'rgb(110 160 255)',
};
const refusing: Loaders = {
  load: (r) => Promise.reject(new Error(`the host passed no loaders (${r.uri ?? r.hostFile ?? 'a file'})`)),
};

export interface SeqViewerContext {
  readonly spec: SeqViewerSpec;
  readonly config: ConfigManager;
  readonly registry: SeqViewerRegistry;
  /** seqspec:viewer#registration: both halves, conflicts first, an idempotent undo. */
  register(entry: SeqViewerRegistryEntry | readonly SeqViewerRegistryEntry[]): () => void;
  readonly graph: SeqGraph;
  readonly document: DocumentManager;
  readonly history: HistoryManager;
  readonly views: ViewManager;
  /** The layout's and engine's sizes (seqspec:viewer#render-props): `spec.render`, then `setRenderProps`. */
  readonly renderProps: SeqState<RenderProps>;
  readonly text: TextManagerApi;
  readonly ui: UiManager;
  readonly tasks: TaskManager;
  /** The graph's issues by node, then the others (registration, behaviors, snapshots, the engine). */
  readonly issues: IssueManager;
  readonly layout: LayoutManager;
  /** Data held by region, loaded by view (seqspec:sources#the-viewers-side). */
  readonly streams: StreamManager;
  /**
   * Input for consumers (seqspec:interaction#streams): streams of fresh payloads after the behaviors ran, observe-only,
   * scoped by `on`; built only while observed.
   */
  readonly interaction: SeqInteraction;
  /** Requests by registered id (seqspec:viewer#commands): `dispatch` from code, `run` from plain data. */
  readonly commands: CommandManager;
  /**
   * Where the engine (or a test) hands input records; `events`, the base events behaviors claim from
   * (seqspec:interaction#input-events), for behavior authors.
   */
  readonly input: SeqInput;
  readonly behaviors: BehaviorManager;
  /** Written by the input dispatcher. */
  readonly pointer: Live<PointerState>;
  readonly hover: HoverManager;
  readonly highlight: HighlightManager;
  /**
   * Coordinate systems at the edge (seqspec:libraries#the-system-registry): core's registry over the shown set (a
   * host's systems with their maps, the set's own by key) plus its reactive `version`.
   */
  readonly systems: ViewerSystems;
  readonly labels: LabelManager;
  readonly viewports: ViewportManager;
  readonly scrollGroups: ScrollGroupManager;
  readonly focus: FocusManager;
  /**
   * Brings a target (loci, `Coordinates`, or a list: their hull) into view in every section that shows it, or
   * `opts.section`, changing no selection, focus or zoom: text scrolls its lines in, tracks pan at their cell; only
   * where it is not shown already (seqspec:interaction#focus). False: no section shows it. The context's, not a
   * manager's: it holds no state (a focus is state, a reveal a move) and needs both the projector (the viewports) and
   * the text's lines (the layout).
   */
  reveal(target: SeqLociAny | Coordinates | readonly (SeqLociAny | Coordinates)[], opts?: RevealOptions): boolean;
  readonly selection: SelectionManager;
  readonly refs: RefManager;
  /** Results an application found, lit and stepped (seqspec:interaction#results); bookmarks step through it too. */
  readonly results: ResultsManager;
  readonly snapshots: SnapshotManager;
  /** Registered actions by name (seqspec:viewer#actions). */
  readonly actions: ActionManager;
  readonly skin: SeqState<Skin>;
  /** The port, tracked: `dispose` clears every timer and frame set through it, and refuses later ones. */
  readonly scheduler: Scheduler;
  readonly events: SeqViewerEvents;
  readonly disposed: boolean;
  /** Validated, then the base (by `diff` over the previous one); resolves when the graph settles. */
  setDocument(doc: SeqViewSpecDoc, opts?: { origin?: Origin }): Promise<GraphResult>;
  /** Data without a document: a set, or build inputs (one node per entry), with views by the builder. */
  setData(data: SeqDataInput, opts?: SetDataOptions): Promise<GraphResult>;
  /** `snapshot`: a local patch's tag (default: the snapshot shown). */
  applyPatch(
    patch: Patch,
    opts?: { layer?: 'base' | 'local' | 'preview'; snapshot?: string | null; origin?: Origin },
  ): Promise<GraphResult>;
  /** A host override of the skin until `setSkin(null)`, which shows the skin the reader last read. */
  setSkin(skin: Skin | null): void;
  /** The skin reader's: remembered, and shown unless a host's `setSkin` holds (seqspec:interaction#the-skin). */
  publishSkin(skin: Skin): void;
  setRenderProps(props: Partial<RenderProps>): void;
  /** Replaces the text port: every prepared text is stale (`text.version`). */
  setTextMeasurer(m: TextMeasurer): void;
  /** False pauses a covered or hidden viewer: animations jump to their end, a press in progress cancels, the player pauses. */
  setActive(on: boolean): void;
  /**
   * An issue from outside the managers (the engine's `seqv.renderer.error`, a host's part): reported and logged, at its
   * registered severity (an error when unregistered; seqspec:format#validation).
   */
  report(code: string, message: string): void;
  dispose(): void;
}

/**
 * The context's own subscriptions, in dependency order, each guarded: a throw is reported as `seqv.state.error` and the
 * rest still run. A new shown document is two passes, so every manager keeps or makes a section's objects (`sync`, kept
 * across a removal and return within one update) before any emits (`emit`: `sections` after the layout and ui sync).
 * One graph emission feeds `assembled`, axes first.
 */
interface Fanout {
  readonly sync: ((doc: SeqViewSpecDoc | null) => void)[];
  readonly emit: (() => void)[];
  readonly assembled: ((a: Assembled | null) => void)[];
}

class Context implements SeqViewerContext {
  readonly spec: SeqViewerSpec;
  readonly ports: SeqViewerPorts;
  readonly scheduler: TrackedScheduler;
  readonly config: ConfigManager;
  readonly registry: ViewerRegistry;
  readonly issueEvent: SeqEventEmitter<SeqIssue>;
  readonly logEvent: SeqEventEmitter<LogEntry>;
  readonly patchEvent: SeqEventEmitter<PatchEvent>;
  readonly events: SeqViewerEvents;
  readonly graph: Graph;
  readonly document: Documents;
  readonly views = new Views();
  readonly renderProps = state<RenderProps>(DefaultRenderProps);
  readonly text: TextManager;
  readonly ui = new Ui();
  readonly tasks: Tasks;
  readonly issues = new Issues();
  readonly pointer = live<PointerState>({ section: -1, x: 0.5, y: 0.5, inside: false });
  // A header's hover stands for what the `headerHover` behavior says (seqspec:interaction#hover), in the view shown.
  readonly hover = new Hover(
    () => this.viewports?.projector ?? null,
    () => this.graph?.assembled.value?.view ?? null,
    () => this.systems?.version.value ?? 0,
  );
  readonly highlight = new Highlight((m) => this.report('seqv.state.error', m));
  readonly systems = new ViewerSystems();
  readonly labels: Labels;
  readonly input: InputDispatcher;
  readonly behaviors: Behaviors;
  readonly commands: Commands;
  /** The public streams (`interaction`), fed by the dispatcher. */
  readonly out: Interaction;
  readonly layout: Layouts;
  readonly streams: Streams;
  readonly viewports: Viewports;
  readonly scrollGroups: ScrollGroups;
  readonly focus: Focus;
  readonly selection = new Selections(() => this.viewports?.projector ?? null);
  readonly refs: Refs;
  readonly results: Results;
  readonly snapshots: Snapshots;
  readonly actions: ActionManager;
  readonly skin = state(DefaultSkin);
  /** A host's `setSkin` holds the skin; the reader's last, shown again by `setSkin(null)`. */
  private overridden = false;
  private read: Skin = DefaultSkin;
  readonly fanout: Fanout = { sync: [], emit: [], assembled: [] };
  /** Showing `setData`'s views, not a document. */
  private data = false;
  /** In data mode, what refs resolve against: the incoming set or a resolver over the inputs. */
  private keys: KeyResolver | undefined;
  disposed = false;

  constructor(spec: SeqViewerSpec, ports: Partial<SeqViewerPorts>) {
    this.spec = spec;
    this.document = new Documents({
      svs: () => this.registry.svs,
      set: () => this.keys,
      base: () => this.ports.loaders.base,
      settle: (issues, ok) => this.settle(issues, ok, true),
      port: ports.document ?? null,
      historyLimit: () => this.config.get(SeqViewerConfig.History.Limit),
      report: (i) => this.issueEvent.emit({ ...i, ref: null }),
      patched: (e) => this.patchEvent.emit(e),
    });
    this.ports = {
      loaders: ports.loaders ?? refusing,
      document: this.document.port,
      scheduler: ports.scheduler ?? manualScheduler(),
      text: ports.text ?? monospaceMeasurer,
      log: ports.log,
    };
    this.scheduler = trackedScheduler(this.ports.scheduler);
    this.config = createConfigManager(spec.config);
    this.registry = new ViewerRegistry();
    const onError = (e: unknown) => this.log('error', 'an event subscriber threw', e);
    this.issueEvent = event(onError);
    this.logEvent = event(onError);
    this.patchEvent = event(onError);
    this.events = {
      issue: this.issueEvent,
      log: this.logEvent,
      patch: this.patchEvent,
    };
    this.graph = new Graph({
      registry: this.registry,
      loaders: this.ports.loaders,
      skin: this.skin,
      colorSchemes: this.ui.colorSchemes,
      colorRanges: this.ui.colorRanges,
      doc: () => this.document.shown.value,
      currentView: () => this.views.current.value,
      onIssue: (i) => this.issueEvent.emit(i),
      issueLimit: () => this.config.get(SeqViewerConfig.Issues.PerNode),
      // Read when a node runs, after the stream manager below exists.
      dynamic: () => this.streams.nodes,
    });
    this.tasks = new Tasks(this.scheduler, this.config);
    this.graph.issues.subscribe((l) => this.issues.setGraph(l));
    this.issueEvent.subscribe((i) => this.issues.add(i));
    // VC3 and LY1: render props, text, viewports, the layout (it builds the axes and hands them to the viewports),
    // scroll groups, focus.
    this.renderProps.set(mergeProps(spec.render));
    this.text = new TextManager(this.ports.text, () => this.config.get(SeqViewerConfig.Text.CacheEntries));
    // A disposed or paused (`setActive(false)`) viewer is inactive.
    const active = () => !this.disposed && this.ui.active.value;
    this.viewports = new Viewports({ scheduler: this.scheduler, config: this.config, active });
    this.viewports.projector.systems = this.systems;
    this.layout = new Layouts({
      registry: this.registry,
      config: this.config,
      renderProps: this.renderProps,
      text: this.text,
      ui: this.ui,
      viewports: this.viewports,
      scheduler: this.scheduler,
      report: (code, text) => this.report(code, text),
      // A still pointer's hover follows a rebuild (the dispatcher leaves the keyboard's, and its own hit's, alone).
      rebuilt: (i) => this.input?.rehit(i),
      // And a moved view's, once a frame: the hover follows what is under the pointer now.
      moved: (i) => this.input?.moved(i),
      stream: (layer) => this.streams.get(layer),
    });
    // LY4: the stream manager, which reads the layout and the viewports.
    this.streams = new Streams({
      registry: this.registry,
      config: this.config,
      scheduler: this.scheduler,
      layout: this.layout,
      viewports: this.viewports,
      active,
      has: (ref) => nodePath(this.document.shown.value, ref) !== null,
      expanded: (track) => this.ui.expanded.value.has(track),
      report: (code, text) => this.report(code, text),
      nodesChanged: () => void this.run(),
      held: (out) => {
        const s = this.selection.current.value;
        for (const l of [s.seq, s.boundaries, s.items]) heldBy(l, out);
        for (const e of this.highlight.channels.value.values()) if (!isCoordinates(e.target)) heldBy(e.target, out);
      },
    });
    this.scrollGroups = new ScrollGroups(this.viewports, (code, text) => this.report(code, text));
    this.focus = new Focus(this.viewports, onError);
    const f = this.fanout;
    f.sync.push(
      (doc) => this.viewports.sync(doc),
      (doc) => this.layout.sync(sectionsOf(doc)),
      (doc) => this.ui.keep(sectionsOf(doc)),
    );
    f.emit.push(() => this.viewports.publish());
    const view = (a: Assembled | null) => a?.view?.resolved ?? null;
    // The view's managers update when the view does: a set that changes under the same view (a dynamic node no section
    // reads, ST1b-b) lays nothing out again.
    let shown: Assembled['view'] | undefined;
    let moved = true;
    f.assembled.push(
      (a) => {
        moved = a?.view !== shown;
        shown = a?.view;
      },
      (a) => moved && this.layout.update(a?.view ?? null),
      (a) => moved && this.streams.update(a?.view ?? null),
      (a) => moved && this.scrollGroups.update(view(a), a?.view?.set ?? null),
      (a) => moved && this.focus.update(view(a), a?.view?.set ?? null),
      (a) => moved && this.ui.update(view(a), a?.view?.set ?? null),
    );
    // VC2: refs, selection (carried to every set, lineage.ts), snapshots, actions.
    this.refs = new Refs({ set: () => this.graph.assembled.value?.set ?? null, doc: () => this.document.shown.value });
    f.assembled.push((a) => this.selection.update(a?.set ?? null));

    this.snapshots = new Snapshots({
      documents: this.document,
      views: this.views,
      scheduler: this.scheduler,
      update: () => this.run().then(() => undefined),
      report: (i) => this.issueEvent.emit({ ...i, ref: null }),
      enter: (t) => {
        this.focus.transition = t;
      },
    });
    this.actions = createActions(this);
    // LS1: the shown sections' `highlights` as channels, set again only when an entry or the set changes. LS2: results.
    const sets = new SetChannels(this.highlight);
    this.results = new Results(this.highlight, this.focus, () => this.graph.assembled.value?.set ?? null);
    // IN3: the shown set's own systems, then the held targets carried to every set; then LS1's.
    f.assembled.push(
      (a) => this.systems.update(a?.set ?? null),
      (a) => this.hover.update(a?.set ?? null),
      (a) => this.highlight.update(a?.set ?? null),
      (a) => sets.update(a),
      (a) => this.results.update(a?.set ?? null),
    );
    // IN2: labels, the input dispatcher, the behaviors (added by `createSeqViewer`, after the registry).
    this.labels = new Labels({
      registry: this.registry,
      hover: this.hover.live,
      view: () => this.graph.assembled.value?.view ?? null,
      stream: (layer) => this.streams.get(layer),
      streamVersion: () => this.streams.version.value,
      letters: (data, s, at, n, out) => this.streams.lettersAt?.(data, s, at, n, out) ?? 0,
      cutText: (h) => cutHeaderText(this, h),
    });
    this.out = new Interaction({
      sectionRef: (i) => this.layout.sectionRef(i),
      assembled: () => this.graph.assembled.value,
      now: () => this.scheduler.now(),
    });
    this.input = new InputDispatcher({
      out: this.out,
      layout: this.layout,
      pointer: this.pointer,
      hover: this.hover,
      config: this.config,
      scheduler: this.scheduler,
      active,
      setDrag: (d) => this.ui.setDrag(d),
      onError: (e) => this.report('seqv.behavior.error', `a behavior's subscriber threw: ${message(e)}`),
    });
    // CM1: commands, before the behaviors that handle them.
    this.commands = new Commands(this.registry, (id, e) =>
      e instanceof TargetRangeError
        ? this.report('seqv.target.range', `command "${id}": ${message(e)}`)
        : this.report('seqv.command.error', `command "${id}" failed: ${message(e)}`),
    );
    this.behaviors = new Behaviors(this, (name, e) =>
      this.report('seqv.behavior.error', `behavior "${name}" threw and was disposed: ${message(e)}`),
    );
    this.watch(this.selection.current, () => this.streams.rehold());
    this.watch(this.highlight.channels, () => this.streams.rehold());
    // Regions arrived or left: what is under a still pointer may be new (an item), hit again at the next frame.
    this.watch(this.streams.version, () => this.input.moved(-1));
    // ST6: the focus regions of the snapshots next to the one shown are prefetched (`Stream.Prefetch`).
    const hint = () => this.streams.hint(snapshotHints(this.document.base.value, this.snapshots.current.value));
    this.watch(this.document.base, hint);
    this.watch(this.snapshots.current, hint);
    this.watch(this.graph.assembled, (a) => {
      for (const fn of f.assembled) this.guard(fn, a);
    });
    const update = () => void this.run();
    this.watch(this.document.shown, (doc) => {
      this.views.update(doc);
      for (const fn of f.sync) this.guard(fn, doc);
      for (const fn of f.emit) this.guard(fn, undefined);
      update();
    });
    for (const s of [
      this.views.current,
      this.registry.version,
      this.skin,
      this.ui.colorSchemes,
      this.ui.colorRanges,
    ] as SeqState<unknown>[])
      this.watch(s, update);
  }

  get interaction(): SeqInteraction {
    return this.out.api;
  }

  /** A subscription of the context's own, guarded. */
  watch<T>(source: SeqState<T>, fn: (v: T) => void) {
    return source.subscribe((v) => this.guard(fn, v));
  }

  private guard<T>(fn: (v: T) => void, v: T) {
    try {
      fn(v);
    } catch (e) {
      this.report('seqv.state.error', `a viewer subscriber threw: ${message(e)}`);
    }
  }

  reveal(target: SeqLociAny | Coordinates | readonly (SeqLociAny | Coordinates)[], opts?: RevealOptions): boolean {
    this.check();
    return revealTarget(this.viewports, this.layout, target, opts);
  }
  register(entry: SeqViewerRegistryEntry | readonly SeqViewerRegistryEntry[]) {
    this.check();
    return this.registry.register(entry);
  }

  get history(): HistoryManager {
    return this.document.history;
  }

  async setDocument(doc: SeqViewSpecDoc, opts: { origin?: Origin } = {}): Promise<GraphResult> {
    this.check();
    // An archive's document: its relative URIs are valid against the base its loaders carry, as `compile` reads them
    // (seqspec:sources#archives).
    const issues = validate(doc, this.registry.svs, { base: this.ports.loaders.base });
    if (issues.some((i) => i.severity === 'error')) return { ok: false, issues };
    this.issues.clear();
    // Its warnings are listed until the next document, as a base patch's are (an inline style's bad key,
    // seqspec:format#sections).
    for (const i of issues) this.issueEvent.emit({ ...i, ref: null });
    if (this.data) {
      this.data = false;
      this.keys = undefined;
      this.graph.setData(null);
      this.document.reset();
    }
    // A new document (not an echo of this viewer's commits): every layer its own color again
    // (seqspec:interaction#color-schemes-at-runtime).
    if (this.document.replace(doc, opts.origin, issues)) this.ui.clearColorSchemes();
    return this.settle(issues, true);
  }

  /** The views are validated against the incoming data before the run; invalid ones change nothing. */
  async setData(data: SeqDataInput, opts?: SetDataOptions): Promise<GraphResult> {
    this.check();
    const keys = dataKeys(data);
    const doc = dataViews(data, opts);
    const issues = validate(doc, this.registry.svs, { set: keys });
    if (issues.some((i) => i.severity === 'error')) return { ok: false, issues };
    this.issues.clear();
    // Its warnings are listed until the next document, as a base patch's are (an inline style's bad key).
    for (const i of issues) this.issueEvent.emit({ ...i, ref: null });
    if (!this.data) this.document.reset();
    this.data = true;
    this.keys = keys;
    this.graph.setData(dataMode(data));
    if (this.document.replace(doc, undefined, issues)) this.ui.clearColorSchemes();
    return this.settle(issues, true);
  }

  async applyPatch(
    patch: Patch,
    opts: { layer?: 'base' | 'local' | 'preview'; snapshot?: string | null; origin?: Origin } = {},
  ) {
    this.check();
    const d = this.document;
    const r =
      opts.layer === 'preview'
        ? d.previewPatch(patch)
        : opts.layer === 'local'
          ? d.addLocal(patch, opts.snapshot)
          : d.applyBase(patch, opts.origin, true);
    return this.settle(r.issues, r.ok, opts.layer === 'preview');
  }

  setSkin(skin: Skin | null) {
    this.check();
    this.overridden = skin !== null;
    this.skin.set(skin ?? this.read);
  }

  publishSkin(skin: Skin) {
    this.check();
    this.read = skin;
    if (!this.overridden) this.skin.set(skin);
  }

  setRenderProps(props: Partial<RenderProps>) {
    this.check();
    this.renderProps.set(mergeProps(this.renderProps.value, props));
  }

  setTextMeasurer(m: TextMeasurer) {
    this.check();
    this.text.setMeasurer(m);
  }

  setActive(on: boolean) {
    this.check();
    const was = this.ui.active.value;
    this.ui.active.set(on);
    if (!on) {
      this.viewports.finish();
      this.input.cancel();
    } else if (!was) this.streams.resume();
    this.snapshots.setActive(on);
  }

  /** A graph update, counted by the tasks manager when there is something to run. */
  private run() {
    const work = this.graph.update();
    return this.document.shown.value || this.graph.assembled.value ? this.tasks.track(work) : work;
  }

  /**
   * `preview`: its `svs.snapshot.invalid` warnings are returned, and listed in `issues.current` until the next preview
   * replaces them; an invalid preview changes nothing, its list included.
   */
  private async settle(issues: readonly Issue[], ok: boolean, preview = false): Promise<GraphResult> {
    if (preview && ok)
      this.issues.setPreview(issues.flatMap((i) => (i.code === 'svs.snapshot.invalid' ? [{ ...i, ref: null }] : [])));
    await this.run();
    return { ok, issues: [...issues, ...this.graph.issues.value] };
  }

  report(code: string, text: string) {
    const issue: SeqIssue = { ...issueOf(this.registry.svs, code, '', text), ref: null };
    this.issueEvent.emit(issue);
    this.log(issue.severity, text, issue);
  }

  log(level: LogEntry['level'], text: string, data?: unknown) {
    const entry: LogEntry = { level, message: text, data };
    try {
      this.ports.log?.(entry);
    } catch {
      // A host's log that throws loses that entry, never the caller's work.
    }
    this.logEvent.emit(entry);
  }

  check() {
    if (this.disposed) throw new Error('SeqViewerContext: used after dispose');
  }

  /** Disposal, the parts built so far: timers, the graph's work, states and events, the registry. */
  dispose() {
    if (this.disposed) return;
    this.behaviors.dispose();
    this.input.dispose();
    this.out.dispose();
    this.scheduler.clearAll();
    this.graph.dispose();
    this.disposed = true;
    this.scrollGroups.dispose();
    this.layout.dispose();
    this.streams.dispose();
    this.commands.dispose();
    this.viewports.dispose();
    this.focus.dispose();
    this.ui.dispose();
    this.tasks.dispose();
    this.issues.dispose();
    this.renderProps.complete();
    this.snapshots.dispose();
    this.selection.dispose();
    this.labels.dispose();
    this.pointer.complete();
    this.hover.live.complete();
    this.results.dispose();
    this.highlight.dispose();
    this.systems.version.complete();
    this.document.dispose();
    this.views.dispose();
    this.skin.complete();
    this.issueEvent.complete();
    this.logEvent.complete();
    this.patchEvent.complete();
    this.registry.dispose();
  }
}

const message = (e: unknown) => (e instanceof Error ? e.message : String(e));

/** Builds the managers, registers `spec.registry`, then creates `spec.behaviors` in order. */
export function createSeqViewer(spec: SeqViewerSpec, ports: Partial<SeqViewerPorts> = {}): SeqViewerContext {
  const ctx = new Context(spec, ports);
  if (spec.registry) ctx.register(spec.registry);
  for (const b of spec.behaviors) ctx.behaviors.add(b.behavior, b.defaultParams);
  return ctx;
}

/**
 * The keys a selection's or highlight's loci name (units and their molecules, items' annotations): what it holds of the
 * stream. The hover holds nothing: what it points at is drawn, held by its placement or its region.
 */
function heldBy(l: SeqLociAny | null, out: Set<string>) {
  if (l?.kind === 'item-loci')
    for (const e of l.elements) {
      const a = l.set.annotations[e.annotation];
      if (a) out.add(l.set.strings.get(a.key));
    }
  if (l?.kind !== 'seq-loci' && l?.kind !== 'boundary-loci') return;
  const set = l.set;
  for (const e of l.elements) {
    out.add(set.strings.get(set.units.key[e.unit]));
    const m = set.molecules[set.units.molecule[e.unit]];
    if (m) out.add(set.strings.get(m.key));
  }
}
