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
} from '@molstar/ai/trace';
import { keyedSelector, observableState, type ReadonlyState } from '@molstar/prism-ui';
import type { KeyedValue } from '@molstar/prism-ui/reactive/keyed-selector';
import { derive } from '@molstar/prism-ui/reactive/state';
import { debounceTime, merge, Observable, skip, take } from 'rxjs';

/** The Inspector's three views (packages/ai/_spec/trace.md § The Inspector). */
export type InspectorView = 'trace' | 'tools' | 'limits';

/** The overview strip's scale: every step the same width, or widths that follow how long each took. */
export type OverviewScale = 'steps' | 'time';

/** What a request row shows of its turn: the person's words and focus, not the turn's progress. */
export type TurnRequest = Pick<TurnView, 'id' | 'prompt' | 'focus'>;

/**
 * A row of the trace: a turn's header, its request, or one of its steps (a model call, a tool call; a nested run's
 * steps under the call that ran them). Keys: `turn:<id>`, `req:<id>`, `step:<id>`. A row holds what it shows and no
 * more, so it stays the same object (and its view stays) while that is the same: a step row while its step is, a
 * request row while the words are, however the turn around them grows.
 */
export type TraceRow =
  | { readonly kind: 'turn'; readonly key: string; readonly turn: TurnView; readonly index: number }
  | { readonly kind: 'request'; readonly key: string; readonly turn: TurnRequest; readonly index: number }
  | {
      readonly kind: 'step';
      readonly key: string;
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

/** Whether two lists hold the same items in the same order. */
function sameItems<T>(a: readonly T[], b: readonly T[]): boolean {
  return a === b || (a.length === b.length && a.every((x, i) => x === b[i]));
}

/**
 * `source` as it settles: the first value at once, then the latest once a burst of changes has settled (`ms`), followed
 * only while something listens (a state, to `derive` from: the rows are re-read while a turn runs, not per token).
 */
function settledState<T>(source: ReadonlyState<T>, ms: number): ReadonlyState<T> {
  if (ms <= 0) return source;
  let connected = 0;
  let latest: T;
  const settled = new Observable<T>((subscriber) => {
    connected += 1;
    const sub = merge(source.pipe(take(1)), source.pipe(skip(1), debounceTime(ms))).subscribe((v) => {
      latest = v;
      subscriber.next(v);
    });
    return () => {
      connected -= 1;
      sub.unsubscribe();
    };
  });
  // Unfollowed there is no burst to wait for: the source's own value.
  return Object.defineProperty(settled, 'value', {
    get: () => (connected > 0 ? latest : source.value),
  }) as ReadonlyState<T>;
}

/**
 * The agent's Inspector (Trace, Tools, Limits): which view, the selected trace row and tool, the overview's scale, the
 * detail tabs; and the jumps between them (a tool call's definition to its tool, a page focus). The data is the
 * session's, as `@molstar/ai/trace` prepares it; the mocks feed fixed data.
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
  readonly turns: ReadonlyState<readonly TurnView[]>;
  /** The trace's rows (re-read once a burst of changes settles): a row is kept while what it shows is the same. */
  readonly rows: ReadonlyState<readonly TraceRow[]>;
  /** The selected row's key: the chosen one while it is there, else the latest row. */
  readonly selected: ReadonlyState<string | undefined>;
  readonly row: ReadonlyState<TraceRow | undefined>;
  readonly toolRows: ReadonlyState<readonly ToolsRow[]>;
  /** The selected group or tool's key: the chosen one while it is there, else the first group. */
  readonly selectedTool: ReadonlyState<string | undefined>;
  readonly toolRow: ReadonlyState<ToolsRow | undefined>;
  /** Whether a row is the selected one (told only to the two rows a change concerns). */
  readonly isSelected: (key: string) => KeyedValue<boolean>;
  readonly isSelectedTool: (key: string) => KeyedValue<boolean>;

  constructor(sources: InspectorSources) {
    this.trace = sources.trace;
    this.tools = sources.tools;
    this.limits = sources.limits;
    this.#setLimit = sources.setLimit;
    this.editable = !!sources.setLimit;

    const trace = settledState(sources.trace, sources.settle ?? 80);
    // The same turns and rows in a new array (a re-read that changed none of them) are not a change.
    this.turns = derive(trace, (t) => t?.turns ?? NoRows, sameItems);
    this.rows = derive(trace, (t) => this.#rowsOf(t));
    const pick = <R extends { key: string }>(rows: readonly R[], key: string | undefined, fallback: R | undefined) =>
      rows.find((r) => r.key === key) ?? fallback;
    this.row = derive([this.rows, this.#selected.asReadonly], (rows, key) => pick(rows, key, rows.at(-1)));
    this.selected = derive(this.row, (r) => r?.key);
    this.toolRows = derive(sources.tools, (t) => this.#toolRowsOf(t));
    this.toolRow = derive([this.toolRows, this.#tool.asReadonly], (rows, key) => pick(rows, key, rows[0]));
    this.selectedTool = derive(this.toolRow, (r) => r?.key);
    const rowSelector = keyedSelector(this.selected);
    const toolSelector = keyedSelector(this.selectedTool);
    this.isSelected = (key) => rowSelector.key(key);
    this.isSelectedTool = (key) => toolSelector.key(key);
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

  // The rows of each trace read once, and of each turn while the turn is the same record (the trace cache of
  // `@molstar/ai/trace` keeps the views of what did not change): a re-read of a long conversation builds the rows of
  // the turn that changed, and a row of it is kept while what it shows is.
  readonly #rowsRead = new WeakMap<TraceData, readonly TraceRow[]>();
  #turnRows = new Map<string, TurnRows>();
  #parts: readonly (readonly TraceRow[])[] = [];
  #rowsAll: readonly TraceRow[] = NoRows;

  #rowsOf(data: TraceData | undefined): readonly TraceRow[] {
    if (!data) return NoRows;
    const read = this.#rowsRead.get(data);
    if (read) return read;
    const previous = this.#turnRows;
    const turnRows = new Map<string, TurnRows>();
    const parts: (readonly TraceRow[])[] = [];
    for (const [index, turn] of data.turns.entries()) {
      const old = previous.get(turn.id);
      const made = old?.turn === turn && old.index === index ? old : turnRowsOf(turn, index, old);
      turnRows.set(turn.id, made);
      parts.push(made.rows);
    }
    // The same rows in the same turns: the list they were last in.
    const before = this.#parts;
    const same = parts.length === before.length && parts.every((p, i) => p === before[i]);
    const rows = same ? this.#rowsAll : parts.flat();
    this.#turnRows = turnRows;
    this.#parts = parts;
    this.#rowsAll = rows;
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

/** The rows of one turn as last built, with the turn and its place they are of. */
interface TurnRows {
  readonly turn: TurnView;
  readonly index: number;
  readonly rows: readonly TraceRow[];
  readonly byKey: ReadonlyMap<string, TraceRow>;
}

/** A turn's rows: its header, its request and its steps; a row of the same turn built before stays while it shows the same. */
function turnRowsOf(turn: TurnView, index: number, old: TurnRows | undefined): TurnRows {
  const byKey = new Map<string, TraceRow>();
  const keep = (row: TraceRow): TraceRow => {
    const was = old?.byKey.get(row.key);
    const same =
      was !== undefined &&
      was.kind === row.kind &&
      was.index === row.index &&
      (row.kind === 'request'
        ? was.kind === 'request' && was.turn.prompt === row.turn.prompt && was.turn.focus === row.turn.focus
        : row.kind === 'step'
          ? was.kind === 'step' && was.step === row.step && was.nested === row.nested
          : was.kind === 'turn' && was.turn === row.turn);
    const kept = same ? (was as TraceRow) : row;
    byKey.set(row.key, kept);
    return kept;
  };
  const request: TurnRequest = { id: turn.id, prompt: turn.prompt, ...(turn.focus ? { focus: turn.focus } : {}) };
  const rows = [
    keep({ kind: 'turn', key: `turn:${turn.id}`, turn, index }),
    keep({ kind: 'request', key: `req:${turn.id}`, turn: request, index }),
    ...turn.steps.map((step) =>
      keep({ kind: 'step', key: `step:${step.id}`, index, step, nested: step.under !== undefined }),
    ),
  ];
  return { turn, index, rows, byKey };
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
