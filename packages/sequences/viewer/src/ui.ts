/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// UI state, tasks and issues (seqspec:viewer#state): runtime state behaviors, actions and hosts set, never in the
// document. LY1 ported the part the layout reads (groups and split tracks expanded, collapsed clades, a section's
// preset at runtime), IN2 the pointer's (the cursor per section, the drag, the window part lit), IN4 the keyboard's
// cursor; UX1 adds the legend's expanded groups and `active`, with the tasks and the issues.

import type { ResolvedView } from "@molstar/seq-core/format";
import { type SeqSet, SortedArray } from "@molstar/seq-core/model";
import { sameList } from "./axis.js";
import { type ConfigManager, SeqViewerConfig } from "./config.js";
import type { SeqIssue } from "./context.js";
import { severityRank } from "./graph.js";
import type { Scheduler } from "./scheduler.js";
import { type Live, live, type SeqState, state, type WritableSeqState } from "./state.js";

export interface DragState {
  readonly section: number;
  readonly owner: string;
  readonly region: number;
}
/**
 * The window part lit: under the pointer, or held by a window's drag from its press. `kind` as windows.ts `WindowKind`
 * (0 the navigator, 1 the vertical overview, 2 the scrollbar; -1 none), `part` a `WindowPart`: an end (its grip), `pan`
 * (the body: every grip), `page` (outside the window: none).
 */
export interface WindowHover {
  section: number;
  kind: number;
  part: number;
}

/**
 * The keyboard's place in one section (seqspec:viewer#behaviors), in the block coordinate system, so a filter, a zoom
 * or a relayout keeps it on its position. `track` null: a group's header (`displayRow` is its group). `row` is the
 * alignment's data row, `displayRow` the layout's row (a rows layer's display row, a split track's value row, the
 * wrapped line; -1: none), which a data row cannot give back without a search.
 */
export interface KeyCursor {
  readonly section: string;
  readonly track: string | null;
  readonly header: boolean;
  readonly block: number;
  readonly index: number;
  readonly row: number;
  readonly displayRow: number;
}

/** The `hoverRow` behavior's params (seqspec:interaction#hover). */
export interface RowHoverParams {
  /** `entry`: a data hover lights its track's box (a split row's, a rows track's row, a renderer's `rowHover`). */
  readonly band: "entry" | "off";
  /** The band's strength (0 to 1) over the data area and over the header. */
  readonly data: number;
  readonly header: number;
}
/** The `hoverRow` behavior's params per section (null: no band there). */
export type RowHoverOf = (section: number) => RowHoverParams | null;

export interface UiManager {
  /** Expanded group refs, and split tracks' refs (their rows shown). */
  readonly expanded: SeqState<ReadonlySet<string>>;
  /** Alignment ref → its collapsed tree nodes. */
  readonly clades: SeqState<ReadonlyMap<string, SortedArray>>;
  /** Section → a preset at runtime. */
  readonly style: SeqState<ReadonlyMap<string, string>>;
  /**
   * Layer → a color scheme at runtime (seqspec:interaction#color-schemes-at-runtime); kept across patches, cleared by a
   * new document.
   */
  readonly colorSchemes: SeqState<ReadonlyMap<string, string>>;
  /** Layer → a color range ([lo, hi], wider or narrower than its data's) at runtime; kept and cleared as the schemes. */
  readonly colorRanges: SeqState<ReadonlyMap<string, readonly [number, number]>>;
  /** A CSS cursor name; per section. */
  cursor(section: string): SeqState<string>;
  readonly drag: SeqState<DragState | null>;
  /** The window part lit (one pointer: one record for all sections); reused, its version bumps on a change. */
  readonly window: Live<WindowHover>;
  /** One per viewer, in the section whose surface has keyboard focus; null: none. Set by the `keyboard` behavior. */
  readonly keyCursor: SeqState<KeyCursor | null>;
  /** The `hoverRow` behavior's, which the overlay draws a data hover's row band by; null: none. A drawing only. */
  readonly rowHover: SeqState<RowHoverOf | null>;
  /**
   * Legend groups expanded past `UI.LegendItems`, by group key `section/trackRef`; never in snapshots
   * (seqspec:interaction#the-legend).
   */
  readonly legendOpen: SeqState<ReadonlySet<string>>;
  /**
   * The legend open (true) or closed, as the person's toggle or `view.legend` left it; null: as the UI's `legend`
   * option says. Kept across patches and a new document (a host resets it with `setLegendExpanded(null)`).
   */
  readonly legendExpanded: SeqState<boolean | null>;
  /** False: a covered or hidden viewer (`ctx.setActive`); the viewports, streams and input read it. */
  readonly active: SeqState<boolean>;
  /** A group, or a split track's rows. */
  toggleGroup(ref: string): void;
  toggleClade(alignment: string, node: number): void;
  toggleLegend(group: string): void;
  setLegendExpanded(open: boolean | null): void;
  setStyle(section: string, preset: string | null): void;
  /** null: the layer's own color again. */
  setColorScheme(layer: string, scheme: string | null): void;
  /** null: the layer's own range again. */
  setColorRange(layer: string, range: readonly [number, number] | null): void;
  setCursor(section: string, cursor: string): void;
  /** The cursor behavior's, with the cursor: bumps `window` only on a change. */
  setWindow(section: number, kind: number, part: number): void;
  setKeyCursor(cursor: KeyCursor | null): void;
  setRowHover(of: RowHoverOf | null): void;
}

export class Ui implements UiManager {
  readonly expanded = state<ReadonlySet<string>>(new Set());
  readonly clades = state<ReadonlyMap<string, SortedArray>>(new Map());
  readonly style = state<ReadonlyMap<string, string>>(new Map());
  readonly colorSchemes = state<ReadonlyMap<string, string>>(new Map());
  readonly colorRanges = state<ReadonlyMap<string, readonly [number, number]>>(new Map());
  readonly drag = state<DragState | null>(null);
  readonly window = live<WindowHover>({ section: -1, kind: -1, part: 0 });
  readonly keyCursor = state<KeyCursor | null>(null);
  readonly rowHover = state<RowHoverOf | null>(null);
  readonly legendOpen = state<ReadonlySet<string>>(new Set());
  readonly legendExpanded = state<boolean | null>(null);
  readonly active = state(true);
  private readonly cursors = new Map<string, WritableSeqState<string>>();
  /** Cursors of sections that left, until the update ends; the sections of the last `keep`. */
  private readonly leaving = new Map<string, WritableSeqState<string>>();
  private sections: readonly string[] = [];
  private pending = false;
  private closed = false;
  /** A group's `collapsed` param when last seen: a change resets its expansion. */
  private readonly seen = new Map<string, number>();

  /**
   * Groups start from their `collapsed` param the first time their ref is seen, and again when a patch changes it; a
   * split track (its rows, seqspec:format#layer-types) likewise from its resolved `expanded`.
   */
  update(view: ResolvedView | null, set: SeqSet | null) {
    let next: Set<string> | null = null;
    const see = (id: number, collapsed: number) => {
      const ref = (set as SeqSet).strings.get(id);
      if (this.seen.get(ref) === collapsed) return;
      this.seen.set(ref, collapsed);
      next ??= new Set(this.expanded.value);
      if (collapsed) next.delete(ref);
      else next.add(ref);
    };
    const [g, t] = [view?.groups, view?.tracks];
    for (let i = 0; g && set && i < g.count; i++) see(g.ref[i], g.collapsed[i]);
    for (let i = 0; t && set && i < t.count; i++) if (t.split[i] >= 0) see(t.ref[i], 1 - t.expanded[i]);
    if (next) this.expanded.set(next, "graph");
  }

  toggleGroup(ref: string) {
    const next = new Set(this.expanded.value);
    if (!next.delete(ref)) next.add(ref);
    this.expanded.set(next, "api");
  }
  toggleClade(alignment: string, node: number) {
    const held = Array.from(this.clades.value.get(alignment) ?? []);
    const nodes = held.includes(node) ? held.filter((n) => n !== node) : [...held, node];
    this.clades.set(new Map(this.clades.value).set(alignment, SortedArray.ofUnsortedArray(Int32Array.from(nodes))));
  }
  setLegendExpanded(open: boolean | null) {
    if (this.legendExpanded.value !== open) this.legendExpanded.set(open);
  }
  toggleLegend(group: string) {
    const next = new Set(this.legendOpen.value);
    if (!next.delete(group)) next.add(group);
    this.legendOpen.set(next);
  }
  setStyle(section: string, preset: string | null) {
    const next = new Map(this.style.value);
    if (preset === null) next.delete(section);
    else next.set(section, preset);
    this.style.set(next);
  }

  setColorScheme(layer: string, scheme: string | null) {
    if ((this.colorSchemes.value.get(layer) ?? null) === scheme) return;
    const next = new Map(this.colorSchemes.value);
    if (scheme === null) next.delete(layer);
    else next.set(layer, scheme);
    this.colorSchemes.set(next);
  }
  setColorRange(layer: string, range: readonly [number, number] | null) {
    const held = this.colorRanges.value.get(layer);
    if (range === null ? held === undefined : held?.[0] === range[0] && held[1] === range[1]) return;
    const next = new Map(this.colorRanges.value);
    if (range === null) next.delete(layer);
    else next.set(layer, range);
    this.colorRanges.set(next);
  }
  /** A new document: every layer its own color and range again. */
  clearColorSchemes() {
    if (this.colorSchemes.value.size > 0) this.colorSchemes.set(new Map());
    if (this.colorRanges.value.size > 0) this.colorRanges.set(new Map());
  }

  cursor(section: string): SeqState<string> {
    const c = this.cursors.get(section) ?? state("default");
    this.cursors.set(section, c);
    return c;
  }
  setCursor(section: string, cursor: string) {
    (this.cursor(section) as WritableSeqState<string>).set(cursor);
  }
  /** The input dispatcher's: the drag's owner while one runs. */
  setDrag(drag: DragState | null) {
    this.drag.set(drag);
  }
  setWindow(section: number, kind: number, part: number) {
    const w = this.window.current;
    const none = kind < 0 || part === 0;
    const s = none ? -1 : section;
    const k = none ? -1 : kind;
    const p = none ? 0 : part;
    if (w.section === s && w.kind === k && w.part === p) return;
    w.section = s;
    w.kind = k;
    w.part = p;
    this.window.bump();
  }

  setKeyCursor(cursor: KeyCursor | null) {
    this.keyCursor.set(cursor);
  }
  setRowHover(of: RowHoverOf | null) {
    this.rowHover.set(of);
  }

  /**
   * Per-section states of sections that left the shown document, with their legend groups, are completed and dropped,
   * but only after the update: a section that leaves and returns within one keeps its objects.
   */
  keep(sections: readonly string[]) {
    this.sections = sections;
    for (const ref of sections) {
      const back = this.leaving.get(ref);
      if (back && this.leaving.delete(ref)) this.cursors.set(ref, back);
    }
    for (const [ref, c] of this.cursors) {
      if (sections.includes(ref)) continue;
      this.cursors.delete(ref);
      this.leaving.set(ref, c);
    }
    // A group's key is `<section>/<track>`, the track a template instance's `<track>/<node ref>`: a section's ref (an
    // author ref) holds no `/`, so it is the text before the first.
    const orphan = [...this.legendOpen.value].some((k) => !sections.includes(k.slice(0, k.indexOf("/"))));
    if (!this.pending && (this.leaving.size > 0 || orphan)) {
      this.pending = true;
      void Promise.resolve().then(this.drop);
    }
  }
  private readonly drop = () => {
    this.pending = false;
    if (this.closed) return;
    for (const c of this.leaving.values()) c.complete();
    this.leaving.clear();
    const open = [...this.legendOpen.value].filter((k) => this.sections.includes(k.slice(0, k.indexOf("/"))));
    if (open.length !== this.legendOpen.value.size) this.legendOpen.set(new Set(open));
  };

  dispose() {
    this.closed = true;
    for (const c of [...this.cursors.values(), ...this.leaving.values()]) c.complete();
    this.cursors.clear();
    this.leaving.clear();
    for (const s of [
      this.expanded,
      this.clades,
      this.style,
      this.colorSchemes,
      this.colorRanges,
      this.drag,
      this.keyCursor,
      this.rowHover,
      this.legendOpen,
      this.legendExpanded,
      this.active,
    ])
      (s as WritableSeqState<unknown>).complete();
    this.window.complete();
  }
}

export interface TaskManager {
  /** True only after `Tasks.BusyDelayMs` of continuous work. */
  readonly busy: SeqState<boolean>;
  readonly progress: SeqState<{ readonly pending: number; readonly done: number }>;
}

/** Counts graph updates requested and settled; `busy` turns on after the delay, off when all are done. */
export class Tasks implements TaskManager {
  readonly busy = state(false);
  readonly progress = state({ pending: 0, done: 0 });
  private timer = -1;
  private readonly scheduler: Scheduler;
  private readonly config: ConfigManager;

  constructor(scheduler: Scheduler, config: ConfigManager) {
    this.scheduler = scheduler;
    this.config = config;
  }

  track<T>(work: Promise<T>): Promise<T> {
    const { pending, done } = this.progress.value;
    if (pending === done) {
      const delay = this.config.get(SeqViewerConfig.Tasks.BusyDelayMs);
      this.timer = this.scheduler.setTimer(() => this.busy.set(true), delay);
    }
    this.progress.set({ pending: pending + 1, done });
    const end = () => {
      const p = this.progress.value;
      if (p.done + 1 < p.pending) return this.progress.set({ pending: p.pending, done: p.done + 1 });
      this.scheduler.clearTimer(this.timer);
      this.busy.set(false);
      this.progress.set({ pending: 0, done: 0 });
    };
    work.then(end, end);
    return work;
  }

  dispose() {
    this.busy.complete();
    this.progress.complete();
  }
}

export interface IssueManager {
  /** The graph's issues by node, then the others (registration, behaviors, snapshots, the engine), grouped by ref. */
  readonly current: SeqState<readonly SeqIssue[]>;
}

/**
 * Issues not from a node are kept until the document is replaced, each once, at most `Limit`: past it the least severe
 * goes first, the oldest among equals (a flood of warnings never evicts an error); the shown preview's snapshot
 * warnings until the next preview.
 */
export class Issues implements IssueManager {
  static readonly Limit = 200;
  readonly current = state<readonly SeqIssue[]>([], sameList);
  private graph: readonly SeqIssue[] = [];
  private other: SeqIssue[] = [];
  private preview: readonly SeqIssue[] = [];

  setGraph(issues: readonly SeqIssue[]) {
    this.graph = issues;
    this.publish();
  }
  add(issue: SeqIssue) {
    if (issue.ref !== null) return;
    const same = (i: SeqIssue) => i.code === issue.code && i.path === issue.path && i.message === issue.message;
    if (this.other.some(same)) return;
    const next = [...this.other, issue];
    if (next.length > Issues.Limit) {
      let out = 0;
      for (let k = 1; k < next.length; k++) if (severityRank(next[k]) > severityRank(next[out])) out = k;
      if (out === next.length - 1) return;
      next.splice(out, 1);
    }
    this.other = next;
    this.publish();
  }
  setPreview(issues: readonly SeqIssue[]) {
    if (issues.length === 0 && this.preview.length === 0) return;
    this.preview = issues;
    this.publish();
  }
  clear() {
    this.other = [];
    this.preview = [];
    this.publish();
  }
  private publish() {
    const all = [...this.graph, ...this.other, ...this.preview];
    this.current.set(
      all.sort((a, b) => (a.ref ?? "￿").localeCompare(b.ref ?? "￿")),
      "graph",
    );
  }
  dispose() {
    this.current.complete();
  }
}
