/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  LimitRowView,
  LimitsData,
  StepView,
  ToolGroupView,
  ToolsData,
  ToolView,
  TraceData,
  TurnView,
} from '@molstar/ai-trace';
import { keyedSelector, observableState, type ReadonlyState } from '@molstar/prism-ui';
import { combineLatest, debounceTime, distinctUntilChanged, map, merge, type Observable, skip, take } from 'rxjs';

/** The Inspector's three views (analysis/viewer-mode.md §5). */
export type InspectorView = 'trace' | 'tools' | 'limits';

/** The overview strip's scale: every step the same width, or widths that follow how long each took. */
export type OverviewScale = 'steps' | 'time';

/**
 * A row of the trace: a turn's header, its request, or one of its steps (a model call, a tool call; a nested run's
 * steps under the call that ran them). Keys: `turn:<id>`, `req:<id>`, `step:<id>`.
 */
export type TraceRow =
  | { readonly kind: 'turn'; readonly key: string; readonly turn: TurnView; readonly index: number }
  | { readonly kind: 'request'; readonly key: string; readonly turn: TurnView; readonly index: number }
  | {
      readonly kind: 'step';
      readonly key: string;
      readonly turn: TurnView;
      readonly index: number;
      readonly step: StepView;
      /** It ran under another call (a tool's helper): indented under it. */
      readonly nested: boolean;
    };

/** A row of the Tools view: a group's header (core, a skill, the single tools) or a tool. Keys: `group:<id>`, `tool:<name>`. */
export type ToolsRow =
  | { readonly kind: 'group'; readonly key: string; readonly group: ToolGroupView }
  | { readonly kind: 'tool'; readonly key: string; readonly group: ToolGroupView; readonly tool: ToolView };

/** The detail pane's tab, per kind of row: a model call, a tool call (the trace), a tool's definition (Tools). */
export interface DetailTabs {
  model: 'saw' | 'answer' | 'timing';
  tool: 'result' | 'args' | 'definition' | 'timing';
  definition: 'input' | 'schema';
}

/** What the Inspector shows, from a session (or a mock). */
export interface InspectorSources {
  trace: ReadonlyState<TraceData | undefined>;
  tools: ReadonlyState<ToolsData>;
  limits: ReadonlyState<LimitsData>;
  /**
   * Keeps a limit (in its row's unit), or drops it back to its default (undefined); the limits follow. Without it the
   * limits are shown as they are.
   */
  setLimit?(key: string, value: number | undefined): void;
  /** How long a burst of trace changes (a running turn's records) settles before the rows re-read it, in ms (80). */
  settle?: number;
}

const NoRows: readonly never[] = [];

/** The first value at once, then the latest once a burst has settled: re-read while a turn runs, not per token. */
function settled<T>(source: Observable<T>, ms: number): Observable<T> {
  return ms > 0 ? merge(source.pipe(take(1)), source.pipe(skip(1), debounceTime(ms))) : source;
}

/**
 * The agent's Inspector (Trace, Tools, Limits): which view, the selected trace row and tool, the overview's scale, the
 * detail tabs; and the jumps between them (a tool call's definition to its tool, a page focus). The data is the
 * session's, as `@molstar/ai-trace` prepares it; the mocks feed fixed data.
 */
export class InspectorModel {
  readonly trace: ReadonlyState<TraceData | undefined>;
  readonly tools: ReadonlyState<ToolsData>;
  readonly limits: ReadonlyState<LimitsData>;
  /** Limits can be changed here (the session keeps them). */
  readonly editable: boolean;

  readonly #view = observableState<InspectorView>('trace');
  readonly #scale = observableState<OverviewScale>('steps');
  readonly #tabs = observableState<DetailTabs>({ model: 'saw', tool: 'result', definition: 'input' });
  readonly #selected = observableState<string | undefined>(undefined);
  readonly #tool = observableState<string | undefined>(undefined);
  readonly #setLimit: InspectorSources['setLimit'];

  /** The trace's turns, re-read once a burst of changes settles (the overview strip). */
  readonly turns: Observable<readonly TurnView[]>;
  /** The trace's rows (re-read once a burst of changes settles): a row is kept while its turn and step are. */
  readonly rows: Observable<readonly TraceRow[]>;
  /** The selected row's key: the chosen one while it is there, else the latest row. */
  readonly selected: Observable<string | undefined>;
  readonly row: Observable<TraceRow | undefined>;
  readonly toolRows: Observable<readonly ToolsRow[]>;
  /** The selected group or tool's key: the chosen one while it is there, else the first group. */
  readonly selectedTool: Observable<string | undefined>;
  readonly toolRow: Observable<ToolsRow | undefined>;
  /** Whether a row is the selected one (only the two rows a change concerns hear of it). */
  readonly isSelected: (key: string) => Observable<boolean>;
  readonly isSelectedTool: (key: string) => Observable<boolean>;

  constructor(sources: InspectorSources) {
    this.trace = sources.trace;
    this.tools = sources.tools;
    this.limits = sources.limits;
    this.#setLimit = sources.setLimit;
    this.editable = !!sources.setLimit;

    const trace = settled(sources.trace, sources.settle ?? 80);
    this.turns = trace.pipe(map((t) => t?.turns ?? NoRows));
    this.rows = trace.pipe(map((t) => this.#rowsOf(t)));
    const pick = <R extends { key: string }>(rows: readonly R[], key: string | undefined, fallback: R | undefined) =>
      rows.find((r) => r.key === key) ?? fallback;
    const rowState = combineLatest([this.rows, this.#selected.asReadonly]).pipe(
      map(([rows, key]) => pick(rows, key, rows.at(-1))),
      distinctUntilChanged(),
    );
    this.row = rowState;
    this.selected = rowState.pipe(
      map((r) => r?.key),
      distinctUntilChanged(),
    );
    this.toolRows = sources.tools.pipe(map((t) => this.#toolRowsOf(t)));
    const toolState = combineLatest([this.toolRows, this.#tool.asReadonly]).pipe(
      map(([rows, key]) => pick(rows, key, rows[0])),
      distinctUntilChanged(),
    );
    this.toolRow = toolState;
    this.selectedTool = toolState.pipe(
      map((r) => r?.key),
      distinctUntilChanged(),
    );
    const rowSelector = keyedSelector(this.selected);
    const toolSelector = keyedSelector(this.selectedTool);
    this.isSelected = (key) => rowSelector.is(key);
    this.isSelectedTool = (key) => toolSelector.is(key);
  }

  get view(): ReadonlyState<InspectorView> {
    return this.#view.asReadonly;
  }

  get scale(): ReadonlyState<OverviewScale> {
    return this.#scale.asReadonly;
  }

  get tabs(): ReadonlyState<DetailTabs> {
    return this.#tabs.asReadonly;
  }

  setView(view: InspectorView): void {
    this.#view.set(view);
  }

  setScale(scale: OverviewScale): void {
    this.#scale.set(scale);
  }

  setTab<K extends keyof DetailTabs>(kind: K, tab: DetailTabs[K]): void {
    this.#tabs.set({ ...this.#tabs.value, [kind]: tab });
  }

  /** Selects a trace row by key (`turn:<id>`, `req:<id>`, `step:<id>`). */
  select(key: string): void {
    this.#selected.set(key);
  }

  /** Moves the trace's selection by `delta` rows (j/k, the arrows), stopping at the ends. */
  move(delta: number): void {
    this.#selected.set(step(this.#rowsOf(this.trace.value), this.#selected.value, delta, 'last'));
  }

  /** Selects a group (`group:<id>`) or a tool (`tool:<name>`) in Tools. */
  selectTool(key: string): void {
    this.#tool.set(key);
  }

  moveTool(delta: number): void {
    const rows = this.#toolRowsOf(this.tools.value);
    this.#tool.set(step(rows, this.#tool.value, delta, 'first'));
  }

  /** From a tool call's definition to the tool in Tools. */
  openTool(name: string): void {
    this.#tool.set(`tool:${name}`);
    this.#view.set('tools');
  }

  /**
   * Sets a limit, in its row's unit, or (undefined) back to its default. A value the row cannot take (not a number, not
   * above zero, over its ceiling, not one of its choices) is refused: false.
   */
  setLimit(key: string, value: number | undefined): boolean {
    const row = findLimit(this.limits.value, key);
    if (!row || !row.editable || !this.#setLimit) return false;
    if (value === undefined) {
      this.#setLimit(key, undefined);
      return true;
    }
    if (!Number.isFinite(value)) return false;
    if (row.choices) {
      if (!row.choices.some((c) => c.value === value)) return false;
    } else {
      if (row.step) value = Math.round(value / row.step) * row.step;
      if (value <= 0 || (row.ceiling !== undefined && value > row.ceiling)) return false;
    }
    this.#setLimit(key, value === row.default ? undefined : value);
    return true;
  }

  /**
   * Applies a page focus (`ShellViewContext.focus`): a view (`trace`, `tools`, `limits`), a turn (`turn:<id>`; the
   * latest when it is not in the trace), a step (`step:<id>`), a tool (`tool:<name>`) or a group (`group:<id>`).
   * Returns whether it was understood.
   */
  focus(text: string): boolean {
    if (text === 'trace' || text === 'tools' || text === 'limits') {
      this.#view.set(text);
      return true;
    }
    const colon = text.indexOf(':');
    const kind = text.slice(0, colon);
    const id = text.slice(colon + 1);
    if (colon < 0 || !id) return false;
    switch (kind) {
      case 'turn': {
        const turns = this.trace.value?.turns ?? [];
        const turn = turns.find((t) => t.id === id) ?? turns.at(-1);
        if (turn) this.#selected.set(`turn:${turn.id}`);
        this.#view.set('trace');
        return true;
      }
      case 'step':
        this.#selected.set(`step:${id}`);
        this.#view.set('trace');
        return true;
      case 'tool':
        this.openTool(id);
        return true;
      case 'group':
        this.#tool.set(`group:${id}`);
        this.#view.set('tools');
        return true;
    }
    return false;
  }

  // The rows of each trace read once; a row kept from the previous read while its turn and step are the same records
  // (ai-trace's cache keeps them across re-reads), so the list renders again only what changed.
  readonly #rowsRead = new WeakMap<TraceData, readonly TraceRow[]>();
  #rowRecords = new Map<string, TraceRow>();

  #rowsOf(data: TraceData | undefined): readonly TraceRow[] {
    if (!data) return NoRows;
    const read = this.#rowsRead.get(data);
    if (read) return read;
    const previous = this.#rowRecords;
    const records = new Map<string, TraceRow>();
    const keep = (row: TraceRow) => {
      const old = previous.get(row.key);
      const same =
        old &&
        old.kind === row.kind &&
        old.turn === row.turn &&
        old.index === row.index &&
        (old.kind !== 'step' || (row.kind === 'step' && old.step === row.step));
      const kept = same ? old : row;
      records.set(row.key, kept);
      return kept;
    };
    const rows = data.turns.flatMap((turn, index): TraceRow[] => [
      keep({ kind: 'turn', key: `turn:${turn.id}`, turn, index }),
      keep({ kind: 'request', key: `req:${turn.id}`, turn, index }),
      ...turn.steps.map((step) =>
        keep({ kind: 'step', key: `step:${step.id}`, turn, index, step, nested: step.under !== undefined }),
      ),
    ]);
    this.#rowRecords = records;
    this.#rowsRead.set(data, rows);
    return rows;
  }

  readonly #toolRowsRead = new WeakMap<ToolsData, readonly ToolsRow[]>();

  #toolRowsOf(data: ToolsData): readonly ToolsRow[] {
    let rows = this.#toolRowsRead.get(data);
    if (!rows) {
      rows = data.groups.flatMap((group): ToolsRow[] => [
        { kind: 'group', key: `group:${group.id}`, group },
        ...group.tools.map((tool): ToolsRow => ({ kind: 'tool', key: `tool:${tool.name}`, group, tool })),
      ]);
      this.#toolRowsRead.set(data, rows);
    }
    return rows;
  }
}

/** The key `delta` rows from `key` (or from the row shown selected while none is: the `from` end), stopping at the ends. */
function step(
  rows: readonly { key: string }[],
  key: string | undefined,
  delta: number,
  from: 'first' | 'last',
): string | undefined {
  if (!rows.length) return key;
  let i = rows.findIndex((r) => r.key === key);
  if (i < 0) i = from === 'first' ? 0 : rows.length - 1;
  return rows[Math.max(0, Math.min(rows.length - 1, i + delta))].key;
}

export function findLimit(limits: LimitsData, key: string): LimitRowView | undefined {
  for (const g of limits.groups) {
    const row = g.rows.find((r) => r.key === key);
    if (row) return row;
  }
  return undefined;
}

/** How many tools and skills the Tools view lists. */
export function toolCounts(tools: ToolsData): { tools: number; skills: number } {
  let n = 0;
  for (const g of tools.groups) n += g.tools.length;
  return { tools: n, skills: tools.groups.filter((g) => g.kind === 'skill').length };
}
