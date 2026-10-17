/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  canonicalJson,
  estimatedTokens,
  type Json,
  type JsonObject,
  type Message,
  LABEL_GUIDANCE,
  LABEL_GUIDANCE_UPDATE,
  LABEL_REQUEST,
  type ModelFinish,
  REQUEST_MARKER,
} from '../messages/index.js';
import {
  type AttemptRecord,
  type CallRecord,
  type Execution,
  type HarnessEvent,
  type HarnessReader,
  type Limits,
  MAIN,
  type StepEvent,
  type TurnSummary,
  type Usage,
} from '../harness/index.js';
import type {
  AttemptView,
  BudgetKind,
  BudgetView,
  ExecutionView,
  LoadedToolView,
  MilestoneView,
  ModelStepView,
  PrefixBreak,
  PreambleView,
  StepView,
  TokensView,
  ToolStepView,
  TraceData,
  TurnStatus,
  TurnView,
  WindowEntryView,
  WindowPart,
} from './types.js';

/**
 * The agent trace (harness.md §11, §13; trace.md): the event log, heads and blobs of one harness as
 * data for the Inspector. Read-only and safe while a turn runs: it shows what has been committed. Everything that
 * belongs to a host (how a request is worded, what a tool call changed, prices) is a hook of the `TraceProfile`, with
 * a default.
 */

/**
 * The labels (`EntryLabel.kind`) of the entries a session and its toolkit write, by the names `session/`
 * and `tools/` give them. An entry is read by its label, never by its text.
 */
export const LABELS = {
  request: LABEL_REQUEST,
  guidance: LABEL_GUIDANCE,
  guidanceUpdate: LABEL_GUIDANCE_UPDATE,
  definitions: 'tools:definitions',
  skill: 'tools:skill',
  entry: 'tools:entry',
} as const;

/** The turn's request as the person wrote it. */
export interface TraceRequest {
  readonly prompt: string;
  readonly focus?: string;
}

/** A tool call, as the profile's `toolCall` hook sees it. */
export interface TraceToolCall {
  readonly callId: string;
  /** The call's key (harness.md §8). */
  readonly key: string;
  readonly toolId: string;
  readonly name: string;
  readonly arguments: JsonObject;
  readonly status: ToolStepView['status'];
  /** The result's text, or the last failure's message. */
  readonly text?: string;
}

/** What the host knows of a tool call that the log does not say. */
export interface TraceToolNotes {
  /** A short label for the call. */
  readonly title?: string;
  /** It only read. */
  readonly read?: boolean;
  /** What it changed: an id to cross-highlight and the label to show. */
  readonly changed?: readonly { readonly id: string; readonly label: string }[];
}

export interface TraceProfile {
  /** A turn's request from the last user message of its begin step. Default: `defaultRequestOf`. */
  readonly requestOf?: (message: Message) => TraceRequest;
  /** Which part of the window a message is. Default: `defaultWindowPart`. */
  readonly windowPart?: (message: Message) => WindowPart;
  /** What the host knows of a tool call. Default: nothing. */
  readonly toolCall?: (call: TraceToolCall) => TraceToolNotes | undefined;
  /** What a model call cost, in the host's credits. Default: nothing is priced. */
  readonly credits?: (model: string, usage: Usage) => number | undefined;
  /**
   * The budget bar's ceilings by kind. Default: the limits the latest turn ran under (the turn's configuration, as
   * extended).
   */
  readonly ceilings?: Partial<Record<BudgetKind, number>>;
  /** Tool ids to the names shown. Default: the name in the tool's definition. */
  readonly toolNames?: ReadonlyMap<string, string>;
  /** The model's context window, for scale. */
  readonly window?: number;
  /** Kept across reads of one harness, so a re-read fetches only what is new. */
  readonly cache?: TraceCache;
}

/**
 * What the trace caches between reads of one harness: heads, entries and blobs are immutable, and the event log is
 * append-only (seq order), so a re-read fetches only the events after the last one it has. What it made of them is
 * kept too (`memo`, trace.md, Incremental reads): a call whose records, and the calls it depends on, are the same
 * records keeps its view, a finished turn keeps its own, so a re-read maps only the changed tail and the rows built
 * on the views keep what did not change.
 */
export interface TraceCache {
  readonly entries: Map<string, Message>;
  readonly heads: Map<string, readonly string[]>;
  readonly blobs: Map<string, Json>;
  /**
   * The events read so far, and the harness they are of (another harness reads its log from the start), with the
   * events grouped by turn as far as `grouped` of them.
   */
  log?: TraceLog;
  /** What was made of the records (views by call and by turn), of the harness `log` is of. */
  memo: TraceMemo;
}

export interface TraceLog {
  readonly harness: string;
  readonly events: HarnessEvent[];
  /** How many of `events` are in `byTurn`, `ended` and `origin`. */
  grouped: number;
  readonly byTurn: Map<string, HarnessEvent[]>;
  /** The summaries the log itself holds (the event that ended a turn). */
  readonly ended: Map<string, TurnSummary>;
  /** The time of the first step event: times are from it. */
  origin?: number;
}

/** What a read made of a model call, valid while the records it came from, and the calls it follows, are the same. */
interface ModelRead {
  readonly call: CallRecord;
  readonly origin: number;
  readonly limit: number | undefined;
  /** The calls it followed: the run's previous call, the run's previous request, the latest request of any run. */
  readonly after: FlowAt;
  readonly view: ModelStepView;
  /** The ids of the messages it saw. */
  readonly ids: readonly string[];
  readonly shape: RequestShape | undefined;
  /** The tool names its window defines. */
  readonly learned: readonly (readonly [string, string])[];
  /** The tools its window loads in full (a call of the turn's own run only). */
  readonly loads: readonly string[];
  /** The preamble its window gives, made only for the call that is first. */
  readonly preamble: PreambleView | undefined;
}

/** What a read of a tool call kept: the host's notes it was made with, to see them change. */
interface ToolRead {
  readonly call: CallRecord;
  readonly origin: number;
  readonly input: TraceToolCall;
  readonly notes: TraceToolNotes | undefined;
  readonly view: ToolStepView;
}

/** What a read of a turn kept, to see that nothing it is made of has changed. */
interface TurnRead {
  readonly view: TurnView;
  readonly count: number;
  readonly last: HarnessEvent | undefined;
  readonly status: TurnStatus;
  readonly summary: string;
  readonly origin: number;
  readonly begin: StepEvent | undefined;
  readonly request: TraceRequest;
  /** The calls it is made of, in order, as read. */
  readonly steps: readonly StepView[];
  /** The flow into it, and out of it. */
  readonly enter: FlowAt & { readonly preamble: boolean; readonly index: number };
  readonly exit: FlowAt;
  /** What it adds to what the reads carry: tool names, the tools first loaded in it, the preamble it makes. */
  readonly learned: ReadonlyMap<string, string>;
  readonly loads: ReadonlySet<string>;
  readonly firstLoads: readonly string[];
  readonly preamble: PreambleView | undefined;
  /** Every call in it was read in full (a missing blob is read again, not kept). */
  readonly whole: boolean;
}

export interface TraceMemo {
  readonly models: Map<string, ModelRead>;
  readonly tools: Map<string, ToolRead>;
  readonly turns: Map<string, TurnRead>;
  /** A message's entry in a window: the same in every call that sees it. */
  readonly windowEntries: Map<string, WindowEntryView>;
  /** The tool definitions in `main`'s head, as far as the head was read. */
  main?: { readonly ids: readonly string[]; readonly definitions: readonly DefinitionPart[] };
  /** The limits the latest turn ran under, for as many events as it had. */
  limits?: { readonly turnId: string; readonly count: number; readonly limits: Limits };
  /** The tools loaded and declared, as the last read gave them (the same arrays while they are the same lists). */
  loaded?: readonly LoadedToolView[];
  declared?: readonly string[];
}

type DefinitionPart = Extract<Message['parts'][number], { kind: 'tool-definition' }>;

/** What a call follows: the model calls whose windows and requests it is compared with. */
interface FlowAt {
  readonly window: ModelRead | undefined;
  readonly request: ModelRead | undefined;
  readonly latest: ModelRead | undefined;
}

const newMemo = (): TraceMemo => ({
  models: new Map(),
  tools: new Map(),
  turns: new Map(),
  windowEntries: new Map(),
});

export function createTraceCache(): TraceCache {
  return { entries: new Map(), heads: new Map(), blobs: new Map(), memo: newMemo() };
}

const PAGE = 500;
const MAX_PAGES = 200;
const SUMMARIES = 50;

const KIND_NAMES: Record<BudgetKind, string> = {
  modelCalls: 'model calls',
  toolCalls: 'tool calls',
  uncachedInput: 'uncached input',
  cacheRead: 'cache reads',
  cacheWrite: 'cache writes',
  output: 'output',
  reasoning: 'reasoning',
};
const KINDS = Object.keys(KIND_NAMES) as BudgetKind[];

/** A turn's state (the mapping is documented at `TurnStatus`). */
export function turnStatusOf(
  status: TurnSummary['status'] | undefined,
  execution: Pick<Execution, 'turn'>,
  turnId: string,
): TurnStatus {
  if (status === 'completed') return 'done';
  if (status === 'failed') return 'failed';
  if (status === 'stopped' || status === 'reset') return 'stopped';
  const turn = execution.turn?.id === turnId ? execution.turn : undefined;
  if (turn === undefined) return 'working';
  if (turn.park !== undefined) return 'waiting';
  return turn.pending.kind === 'model' ? 'thinking' : 'working';
}

/** A usage as the trace shows it. */
export function tokensOf(
  u: Pick<Usage, 'uncachedInput' | 'cacheRead' | 'cacheWrite' | 'output'> &
    Partial<Pick<Usage, 'reasoning' | 'estimated'>>,
  credits?: number,
): TokensView {
  return {
    input: u.uncachedInput + u.cacheRead + u.cacheWrite,
    output: u.output,
    ...(u.cacheRead ? { cached: u.cacheRead } : {}),
    ...(u.cacheWrite ? { written: u.cacheWrite } : {}),
    ...(u.reasoning ? { reasoning: u.reasoning } : {}),
    ...(u.estimated ? { estimated: true } : {}),
    ...(credits !== undefined ? { credits } : {}),
  };
}

/** The text parts of a message, joined. */
function textParts(m: Message): string {
  return m.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');
}

/**
 * The request is the prompt and focus in the label of the turn's request entry (`session:request`); a message without
 * that label records no request.
 */
export function defaultRequestOf(message: Message): TraceRequest {
  const data = message.label?.kind === LABELS.request ? message.label.data : undefined;
  return {
    prompt: (typeof data?.prompt === 'string' ? data.prompt.trim() : '') || '(no request recorded)',
    ...(typeof data?.focus === 'string' && data.focus ? { focus: data.focus } : {}),
  };
}

/**
 * A system message of definitions only is the tools, other system text the guidance; a request entry with something
 * before its request (the host's state, the Focus line, notes) is the state, any other message the conversation.
 */
export function defaultWindowPart(message: Message): WindowPart {
  switch (message.label?.kind) {
    case LABELS.definitions:
      return 'tools';
    case LABELS.guidance:
    case LABELS.guidanceUpdate:
      return 'guidance';
  }
  if (message.role === 'system')
    return message.parts.some((p) => p.kind === 'tool-definition') && !message.parts.some((p) => p.kind === 'text')
      ? 'tools'
      : 'guidance';
  if (message.label?.kind === LABELS.request) {
    // The entry's text ends with `REQUEST_MARKER` and the prompt the label holds: what precedes is the state.
    const prompt = message.label.data?.prompt;
    const text = textParts(message);
    const end = typeof prompt === 'string' ? text.length - prompt.length - REQUEST_MARKER.length : 0;
    if (end > 0 && text.slice(0, end).trim() !== '') return 'state';
  }
  return 'conversation';
}

const isObject = (v: Json | undefined): v is JsonObject => !!v && typeof v === 'object' && !Array.isArray(v);

/** What the harness's own `reason` says: a message, or the name of the reason. */
function reasonText(reason: Json | undefined): string | undefined {
  if (typeof reason === 'string') return reason;
  if (!isObject(reason)) return undefined;
  if (typeof reason.message === 'string') return reason.message;
  return typeof reason.reason === 'string' ? `The turn failed: ${reason.reason}.` : undefined;
}

/** The limits a turn ran under: as extended by its last budget decision, else as configured, else those it holds now. */
function limitsOf(turnEvents: readonly HarnessEvent[], execution: Execution, turnId: string): Limits {
  if (execution.turn?.id === turnId) return execution.turn.budget.limits;
  let begun: Limits | undefined;
  for (let i = turnEvents.length - 1; i >= 0; i -= 1) {
    const e = turnEvents[i] as HarnessEvent;
    if (e.kind !== 'step') continue;
    if (e.limits !== undefined) return e.limits;
    if (e.begun) begun = e.begun.configuration.limits;
  }
  return begun ?? {};
}

/** The latest turn against its budget: the kind nearest its ceiling. */
function budgetView(
  usage: Usage | undefined,
  ceilings: Partial<Record<BudgetKind, number>>,
  credits: number | undefined,
): BudgetView {
  let shown: { kind: BudgetKind; consumed: number; ceiling: number } = {
    kind: 'uncachedInput',
    consumed: 0,
    ceiling: ceilings.uncachedInput ?? 0,
  };
  let share = -1;
  for (const kind of KINDS) {
    const ceiling = ceilings[kind];
    const consumed = usage?.[kind] ?? 0;
    if (ceiling !== undefined && ceiling > 0 && consumed / ceiling > share) {
      share = consumed / ceiling;
      shown = { kind, consumed, ceiling };
    }
  }
  return {
    ...shown,
    dimension: KIND_NAMES[shown.kind],
    scope: 'latest turn',
    ...(credits !== undefined ? { credits } : {}),
  };
}

/** A request as the cache sees it: its parts (provider-given digests) and its model. */
interface RequestShape {
  readonly prefix: readonly string[];
  readonly model: string;
  readonly configuration: string;
}

interface Evidence {
  readonly head?: string | null;
  readonly configuration?: { readonly target?: { readonly model?: string } };
  readonly prefix?: unknown;
}

function requestShape(evidence: Evidence | undefined): RequestShape | undefined {
  if (!Array.isArray(evidence?.prefix)) return undefined;
  return {
    prefix: evidence.prefix as string[],
    model: evidence.configuration?.target?.model ?? '',
    configuration: canonicalJson((evidence.configuration ?? null) as Json),
  };
}

/**
 * Why a request does not extend the previous one, or undefined when it does: the parts are the tools, the system
 * text, then each message, and the provider's cache reaches as far as they agree (harness.md §11).
 */
function breakOf(previous: RequestShape | undefined, now: RequestShape): PrefixBreak | undefined {
  if (previous === undefined) return undefined;
  const of = previous.prefix.length;
  if (previous.model !== now.model) return { reason: `another model (${now.model})`, kept: 0, of };
  if (previous.configuration !== now.configuration) return { reason: "the model's settings changed", kept: 0, of };
  let kept = 0;
  while (kept < of && previous.prefix[kept] === now.prefix[kept]) kept++;
  if (kept === of) return undefined;
  const reason =
    kept === 0
      ? 'its tools changed'
      : kept === 1
        ? 'its system text changed'
        : `the conversation changed at message ${kept - 1} (an edit or a compaction)`;
  return { reason, kept, of };
}

function attemptView(
  a: AttemptRecord,
  ordinal: number,
  origin: number,
  previous: AttemptRecord | undefined,
): AttemptView {
  const waited = previous === undefined ? 0 : a.startedAt - previous.endedAt;
  return {
    ordinal,
    start: a.startedAt - origin,
    end: a.endedAt - origin,
    status: a.outcome,
    ...(a.failure ? { failure: a.failure.message } : {}),
    ...(waited > 0 ? { waitBefore: { reason: 'retry back-off', ms: waited } } : {}),
  };
}

function referenceText(p: { reference: { uri: string; name?: string } }): string {
  return `[${p.reference.name ?? p.reference.uri}]`;
}

function callText(name: string, args: Json): string {
  const text = isObject(args)
    ? Object.entries(args)
        .map(([k, v]) => `${k}=${JSON.stringify(v)}`)
        .join(', ')
    : JSON.stringify(args);
  return `${name}(${text})`;
}

/** The text of a message, as the Inspector shows it. */
export function messageText(m: Message, nameOf: (id: string) => string = (id) => id): string {
  return m.parts
    .map((p) =>
      p.kind === 'text'
        ? p.text
        : p.kind === 'reference'
          ? referenceText(p)
          : p.kind === 'tool-call'
            ? callText(nameOf(p.tool.id), p.arguments)
            : p.kind === 'tool-result'
              ? p.parts.map((q) => (q.kind === 'text' ? q.text : referenceText(q))).join('')
              : `${p.definition.name}: ${p.definition.description}`,
    )
    .join('');
}

/** The tools a context has loaded: a definition in full (a deferred one only declares it). */
function definitionsIn(messages: readonly Message[]) {
  return messages.flatMap((m) =>
    m.parts.flatMap((p) => (p.kind === 'tool-definition' && m.role === 'system' ? [p] : [])),
  );
}

function preambleOf(head: readonly Message[]): PreambleView {
  const system = head.find((m) => m.role === 'system' && m.parts.some((p) => p.kind === 'text'));
  const guidance = system ? textParts(system) : '';
  const loaded = definitionsIn(head).filter((p) => !p.deferred);
  return {
    guidance,
    guidanceTokens: estimatedTokens(guidance),
    tools: loaded.map((p) => ({
      id: p.definition.tool.id,
      name: p.definition.name,
      description: p.definition.description,
    })),
    toolsTokens: loaded.length ? estimatedTokens(JSON.stringify(loaded.map((p) => p.definition))) : 0,
  };
}

/** Tool definitions shown as one entry (they are many and do not change within a conversation). */
function collapseTools(entries: readonly WindowEntryView[]): readonly WindowEntryView[] {
  const tools = entries.filter((e) => e.role === 'tools');
  if (tools.length < 2) return entries;
  const at = entries.indexOf(tools[0] as WindowEntryView);
  const merged: WindowEntryView = {
    role: 'tools',
    part: 'tools',
    title: `${tools.length} tool definition messages`,
    text: tools.map((t) => t.text).join('\n'),
    tokens: tools.reduce((a, t) => a + t.tokens, 0),
    ...(tools.some((t) => t.fresh) ? { fresh: true } : {}),
  };
  const rest = entries.filter((e) => e.role !== 'tools');
  return [...rest.slice(0, at), merged, ...rest.slice(at)];
}

/** The list read before while it is the same list, else the one read now. */
function kept<T>(before: readonly T[] | undefined, now: readonly T[], same: (a: T, b: T) => boolean): readonly T[] {
  return before && before.length === now.length && before.every((x, i) => same(x, now[i] as T)) ? before : now;
}

/** Whether two readings of the host's notes on a tool call are the same. */
function sameNotes(a: TraceToolNotes | undefined, b: TraceToolNotes | undefined): boolean {
  if (a === b) return true;
  if (a?.title !== b?.title || a?.read !== b?.read) return false;
  const x = a?.changed ?? [];
  const y = b?.changed ?? [];
  return x.length === y.length && x.every((c, i) => c.id === y[i]?.id && c.label === y[i]?.label);
}

const sameFlow = (a: FlowAt, b: FlowAt) => a.window === b.window && a.request === b.request && a.latest === b.latest;

/**
 * Reads the trace of one harness: its turns with their model and tool calls (attempts, usage, what each model call
 * saw, where its prompt cache broke), the budget, the preamble and the tools loaded. Times are milliseconds from the
 * first step event.
 *
 * With a `cache` a read is incremental (trace.md, Incremental reads): the events are grouped as they come, a call
 * whose records are the ones it was read from, and that follows the same calls, keeps its view, and a turn that ended
 * and is unchanged keeps its own. Only the changed tail (the running turn's new calls) is mapped, so a read costs what
 * changed, not calls × messages, and the views of what did not change are the same objects.
 */
export async function readTrace(reader: HarnessReader, profile: TraceProfile = {}): Promise<TraceData> {
  const requestOf = profile.requestOf ?? defaultRequestOf;
  const partOf = profile.windowPart ?? defaultWindowPart;
  const cache = profile.cache ?? createTraceCache();

  const execution = await reader.execution();
  if (cache.log?.harness !== reader.id) {
    cache.log = { harness: reader.id, events: [], grouped: 0, byTurn: new Map(), ended: new Map() };
    cache.memo = newMemo();
  }
  const log = cache.log;
  const memo = cache.memo;
  const { events } = log;
  for (let page = 0; page < MAX_PAGES; page += 1) {
    const got = await reader.events({ after: events.at(-1)?.seq ?? 0, limit: PAGE });
    // Two reads at once over one cache: what the other appended meanwhile is not taken twice.
    for (const e of got) if (e.seq > (events.at(-1)?.seq ?? 0)) events.push(e);
    if (got.length < PAGE) break;
  }
  // The events by turn, and the summaries the log holds (the event that ended a turn): only the new ones are looked at.
  for (; log.grouped < events.length; log.grouped += 1) {
    const e = events[log.grouped] as HarnessEvent;
    if (log.origin === undefined && e.kind === 'step') log.origin = e.at;
    if (e.turnId === undefined) continue;
    if (e.summary !== undefined) log.ended.set(e.turnId, e.summary);
    const list = log.byTurn.get(e.turnId);
    if (list) list.push(e);
    else log.byTurn.set(e.turnId, [e]);
  }
  // A turn's summary: from the log, the active turn's (and the latest ended ones') from the reader.
  const fromReader = new Map<string, TurnSummary>();
  for (const t of await reader.turns({ limit: SUMMARIES })) fromReader.set(t.id, t);
  const summaryOf = (turnId: string) => fromReader.get(turnId) ?? log.ended.get(turnId);
  const origin = log.origin ?? 0;

  const messagesOf = async (ids: readonly string[]) => {
    const missing = [...new Set(ids)].filter((id) => !cache.entries.has(id));
    if (missing.length)
      for (const [i, m] of (await reader.entries(missing)).entries()) cache.entries.set(missing[i] as string, m);
    return ids.map((id) => ({ id, message: cache.entries.get(id) as Message }));
  };
  const headIds = async (headId: string | undefined): Promise<readonly string[]> => {
    if (headId === undefined) return [];
    const known = cache.heads.get(headId);
    if (known) return known;
    const [head] = await reader.heads([headId]).catch(() => []);
    if (head === undefined) return [];
    cache.heads.set(headId, head.entries);
    return head.entries;
  };
  const headOf = async (headId: string | undefined) => messagesOf(await headIds(headId));
  const blob = async (ref: string | undefined): Promise<Json | undefined> => {
    if (ref === undefined) return undefined;
    if (cache.blobs.has(ref)) return cache.blobs.get(ref);
    const [body] = await reader.blobs([ref]).catch(() => []);
    if (body !== undefined) cache.blobs.set(ref, body);
    return body;
  };

  // Tool names, learned from the definitions in the heads read so far.
  const names = new Map<string, string>();
  const nameOf = (id: string) => profile.toolNames?.get(id) ?? names.get(id) ?? id;

  // A message's window entry is the same in every call that sees it; only `fresh` differs.
  const windowEntry = (id: string, m: Message): WindowEntryView => {
    const known = memo.windowEntries.get(id);
    if (known) return known;
    const entry = describe(m, partOf(m), nameOf);
    memo.windowEntries.set(id, entry);
    return entry;
  };

  let preamble: PreambleView | undefined;
  // What the runs have seen so far (the turn's run is ''; a nested run's is its call's id): each run's previous model
  // call (its window, to mark what is fresh), each run's previous request and the latest of any run (where the
  // provider's cache broke: a nested run's first request follows its parent's latest).
  const flow = {
    window: new Map<string, ModelRead>(),
    request: new Map<string, ModelRead>(),
    latest: undefined as ModelRead | undefined,
  };
  const flowAt = (key: string, under: string | undefined): FlowAt => ({
    window: flow.window.get(key),
    request: flow.request.get(key),
    latest: under !== undefined ? flow.latest : undefined,
  });
  // The turn (1-based) whose model call first saw each tool loaded; the first call's are "from the start".
  const loadedFirst = new Map<string, number | undefined>();
  const views: TurnView[] = [];

  /** A model call read now (its window mapped): the one place a read costs the size of a window. */
  const readModel = async (call: CallRecord, entering: FlowAt, under: string | undefined) => {
    const attempts = call.attempts.map((a, i) => attemptView(a, i + 1, origin, call.attempts[i - 1]));
    const evidence = (await blob(call.request)) as Evidence | undefined;
    const ids = await headIds(evidence?.head ?? undefined);
    const head = await messagesOf(ids);
    const messages = head.map((e) => e.message);
    const defined = definitionsIn(messages);
    const learned = defined.map((p) => [p.definition.tool.id, p.definition.name] as const);
    for (const [id, name] of learned) names.set(id, name);
    const before = new Set(entering.window?.ids ?? []);
    const entries = head.map((e) => {
      const entry = windowEntry(e.id, e.message);
      return before.has(e.id) ? entry : { ...entry, fresh: true };
    });
    const reply = (await blob(call.result)) as Message | undefined;
    const usage = call.attempts.at(-1)?.usage;
    const model = evidence?.configuration?.target?.model;
    const credits = usage && model ? profile.credits?.(model, usage) : undefined;
    const shape = requestShape(evidence);
    const prefixBreak = shape
      ? breakOf(entering.request?.shape ?? (under !== undefined ? entering.latest?.shape : undefined), shape)
      : undefined;
    const view = {
      kind: 'model',
      id: call.id,
      callId: call.id,
      attempts,
      ...(model ? { model } : {}),
      ...(usage ? { tokens: tokensOf(usage, credits) } : {}),
      ...(under !== undefined ? { under } : {}),
      ...(prefixBreak ? { prefixBreak } : {}),
      ...modelOutput(reply, nameOf),
      window: {
        entries: collapseTools(entries),
        ...(profile.window !== undefined ? { limit: profile.window } : {}),
      },
    } satisfies ModelStepView;
    const read: ModelRead = {
      call,
      origin,
      limit: profile.window,
      after: entering,
      view,
      ids,
      shape,
      learned,
      loads:
        under === undefined ? [...new Set(defined.filter((p) => !p.deferred).map((p) => p.definition.tool.id))] : [],
      preamble: under === undefined && preamble === undefined ? preambleOf(messages) : undefined,
    };
    const whole =
      evidence !== undefined &&
      (evidence.head == null || ids.length > 0) &&
      (call.result === undefined || reply !== undefined);
    return { read, whole };
  };

  /** A tool call: the view kept while the call is the same record and the host's notes on it are the same. */
  const readTool = async (call: CallRecord, under: string | undefined) => {
    const known = memo.tools.get(call.id);
    if (known?.call === call && known.origin === origin) {
      const notes = profile.toolCall?.(known.input);
      if (sameNotes(notes, known.notes)) return { read: known, whole: true };
      const read: ToolRead = { ...known, notes, view: toolView(call, under, known.input, notes) };
      memo.tools.set(call.id, read);
      return { read, whole: true };
    }
    const request = (await blob(call.request)) as { toolId?: string; arguments?: Json } | undefined;
    const result = (await blob(call.result)) as
      | { status: 'succeeded' | 'failed'; parts: readonly { kind: string; text?: string }[] }
      | undefined;
    const last = call.attempts.at(-1);
    const text = result?.parts.map((p) => p.text ?? '').join('') ?? last?.failure?.message;
    const toolId = request?.toolId ?? call.key;
    const status: ToolStepView['status'] =
      result?.status ?? (last === undefined ? 'running' : last.outcome === 'aborted' ? 'aborted' : 'failed');
    const input: TraceToolCall = {
      callId: call.id,
      key: call.key,
      toolId,
      name: nameOf(toolId),
      arguments: isObject(request?.arguments) ? request.arguments : {},
      status,
      ...(text ? { text } : {}),
    };
    const notes = profile.toolCall?.(input);
    const read: ToolRead = { call, origin, input, notes, view: toolView(call, under, input, notes) };
    const whole = request !== undefined && (call.result === undefined || result !== undefined);
    if (whole) memo.tools.set(call.id, read);
    return { read, whole };
  };

  const toolView = (
    call: CallRecord,
    under: string | undefined,
    input: TraceToolCall,
    notes: TraceToolNotes | undefined,
  ): ToolStepView => ({
    kind: 'tool',
    id: call.id,
    callId: call.id,
    toolId: input.toolId,
    name: input.name,
    args: input.arguments,
    status: input.status,
    attempts: call.attempts.map((a, i) => attemptView(a, i + 1, origin, call.attempts[i - 1])),
    ...(notes?.title ? { title: notes.title } : {}),
    ...(notes?.read !== undefined ? { read: notes.read } : {}),
    ...(input.text ? { result: input.text } : {}),
    ...(notes?.changed?.length ? { changed: notes.changed } : {}),
    ...(under !== undefined ? { under } : {}),
  });

  /** The turn's request: the begin step's last user message. */
  const requestFor = async (begin: StepEvent | undefined): Promise<TraceRequest> => {
    const move = begin?.moves.find((m) => m.pointer === MAIN) ?? begin?.moves[0];
    const head = await headOf(move?.to);
    const user = head.findLast((e) => e.message.role === 'user');
    return user ? requestOf(user.message) : { prompt: '(no request recorded)' };
  };

  /** What a model call adds to what the reads carry (its tool names, the tools first loaded, the preamble). */
  const carry = (read: ModelRead, key: string, under: string | undefined, effects?: Effects) => {
    for (const [id, name] of read.learned) {
      names.set(id, name);
      effects?.learned.set(id, name);
    }
    if (under === undefined) {
      const first = preamble === undefined;
      if (first && read.preamble) {
        preamble = read.preamble;
        if (effects) {
          effects.preamble = read.preamble;
          effects.firstLoads = read.loads;
        }
      }
      for (const id of read.loads) {
        effects?.loads.add(id);
        if (!loadedFirst.has(id)) loadedFirst.set(id, first ? undefined : views.length + 1);
      }
    }
    flow.window.set(key, read);
    if (read.shape) {
      flow.request.set(key, read);
      flow.latest = read;
    }
  };

  for (const [turnId, turnEvents] of log.byTurn) {
    const summary = summaryOf(turnId);
    const status = turnStatusOf(summary?.status, execution, turnId);
    const known = memo.turns.get(turnId);
    const last = turnEvents.at(-1);
    const summaryKey = () => (summary === undefined ? '' : JSON.stringify(summary));
    const enter = {
      ...flowAt('', undefined),
      latest: flow.latest,
      preamble: preamble !== undefined,
      index: views.length,
    };

    // A turn whose records and summary are the ones it was read from (every turn but the running one, and it between
    // two records): its view as it was, and what it carries on to the turns after it.
    if (
      known?.whole &&
      known.status === status &&
      known.count === turnEvents.length &&
      known.last === last &&
      known.origin === origin &&
      known.enter.preamble === enter.preamble &&
      known.enter.index === enter.index &&
      sameFlow(known.enter, enter) &&
      known.summary === summaryKey()
    ) {
      for (const [id, name] of known.learned) names.set(id, name);
      if (!enter.preamble && known.preamble) {
        preamble = known.preamble;
        for (const id of known.firstLoads) if (!loadedFirst.has(id)) loadedFirst.set(id, undefined);
      }
      for (const id of known.loads) if (!loadedFirst.has(id)) loadedFirst.set(id, views.length + 1);
      if (known.exit.window) flow.window.set('', known.exit.window);
      if (known.exit.request) flow.request.set('', known.exit.request);
      flow.latest = known.exit.latest;
      views.push(known.view);
      continue;
    }

    const stepEvents = turnEvents.filter((e): e is StepEvent => e.kind === 'step');
    const begin = stepEvents.find((e) => e.stepKind === 'begin');
    // A call is recorded once in the step that ran it; a parked step's later event may repeat it. The calls come in the
    // order they began (a nested run's follow the call that ran them).
    const found = new Map<string, { call: CallRecord; at: number; order: number }>();
    for (const e of stepEvents)
      for (const call of e.calls) {
        const was = found.get(call.id);
        found.set(call.id, { call, at: e.at, order: was?.order ?? found.size });
      }
    const calls = [...found.values()]
      .sort(
        (a, b) =>
          (a.call.attempts[0]?.startedAt ?? a.at) - (b.call.attempts[0]?.startedAt ?? b.at) || a.order - b.order,
      )
      .map((c) => c.call);

    const effects: Effects = { learned: new Map(), loads: new Set(), firstLoads: [], preamble: undefined };
    const steps: StepView[] = [];
    let whole = true;
    for (const call of calls) {
      const under = call.parentCallId;
      const key = under ?? '';
      if (call.kind === 'model') {
        const entering = flowAt(key, under);
        const kept = memo.models.get(call.id);
        const top = under === undefined;
        // The same record, following the same calls: its view as it was.
        let read =
          kept?.call === call &&
          kept.origin === origin &&
          kept.limit === profile.window &&
          sameFlow(kept.after, entering) &&
          (!top || preamble !== undefined || kept.preamble !== undefined)
            ? kept
            : undefined;
        if (read === undefined) {
          const made = await readModel(call, entering, under);
          read = made.read;
          if (made.whole) memo.models.set(call.id, read);
          else whole = false;
        }
        carry(read, key, under, effects);
        steps.push(read.view);
      } else {
        const made = await readTool(call, under);
        if (!made.whole) whole = false;
        steps.push(made.read.view);
      }
    }

    const ended = turnEvents.findLast((e) => (e.kind === 'step' ? e.end !== undefined : e.summary !== undefined));
    const request = known && begin !== undefined && known.begin === begin ? known.request : await requestFor(begin);
    const key = summaryKey();
    const exit: FlowAt = { window: flow.window.get(''), request: flow.request.get(''), latest: flow.latest };
    // The same turn as the last read, call for call: that view again.
    const same =
      known !== undefined &&
      known.status === status &&
      known.count === turnEvents.length &&
      known.last === last &&
      known.origin === origin &&
      known.summary === key &&
      known.request === request &&
      known.steps.length === steps.length &&
      known.steps.every((s, i) => s === steps[i]);
    const view: TurnView = same
      ? known.view
      : {
          id: turnId,
          prompt: request.prompt,
          ...(request.focus ? { focus: request.focus } : {}),
          status,
          start: (begin?.at ?? origin) - origin,
          ...(summary?.endedAt !== undefined ? { end: summary.endedAt - origin } : {}),
          tokens: summary ? tokensOf(summary.usage, sumCredits(steps)) : { input: 0, output: 0 },
          steps,
          milestones: turnEvents.map((e) => milestone(e, origin)),
          ...(status === 'failed'
            ? { failure: reasonText(ended?.kind === 'step' ? ended.end?.reason : undefined) ?? 'The turn failed.' }
            : {}),
        };
    if (whole)
      memo.turns.set(turnId, {
        view,
        count: turnEvents.length,
        last,
        status,
        summary: key,
        origin,
        begin,
        request,
        steps,
        enter,
        exit,
        learned: effects.learned,
        loads: effects.loads,
        firstLoads: effects.firstLoads,
        preamble: effects.preamble,
        whole,
      });
    views.push(view);
  }

  // The latest context: what it has loaded and declares. The definitions in `main`'s head are read once: a head that
  // extends the one read before is read from where that ended.
  const mainIds = await headIds(execution.main);
  let definitions: readonly DefinitionPart[];
  if (memo.main?.ids === mainIds) definitions = memo.main.definitions;
  else {
    const before = memo.main;
    const extends_ =
      before !== undefined && before.ids.length <= mainIds.length && before.ids.every((id, i) => id === mainIds[i]);
    const tail = extends_ ? mainIds.slice(before.ids.length) : mainIds;
    const found = definitionsIn((await messagesOf(tail)).map((e) => e.message));
    definitions = extends_ ? [...before.definitions, ...found] : found;
    memo.main = { ids: mainIds, definitions };
  }
  for (const p of definitions) names.set(p.definition.tool.id, p.definition.name);
  const loadedIds = [...new Set(definitions.filter((p) => !p.deferred).map((p) => p.definition.tool.id))];
  const loaded = kept(
    memo.loaded,
    loadedIds.map((id): LoadedToolView => {
      const turn = loadedFirst.get(id);
      return { id, name: nameOf(id), ...(turn !== undefined ? { turn } : {}) };
    }),
    (a, b) => a.id === b.id && a.name === b.name && a.turn === b.turn,
  );
  const declared = kept(
    memo.declared,
    [
      ...new Set(
        definitions
          .filter((p) => p.deferred && !loadedIds.includes(p.definition.tool.id))
          .map((p) => p.definition.tool.id),
      ),
    ],
    (a, b) => a === b,
  );
  memo.loaded = loaded;
  memo.declared = declared;

  const latest = views.at(-1);
  let ceilings = profile.ceilings;
  if (ceilings === undefined && latest !== undefined) {
    // The limits of the latest turn, scanned once for as many events as it has.
    const turnEvents = log.byTurn.get(latest.id) ?? [];
    const kept = memo.limits;
    const limits =
      kept?.turnId === latest.id && kept.count === turnEvents.length && execution.turn?.id !== latest.id
        ? kept.limits
        : limitsOf(turnEvents, execution, latest.id);
    memo.limits = { turnId: latest.id, count: turnEvents.length, limits };
    ceilings = Object.fromEntries(
      KINDS.flatMap((kind) => (limits[kind] === undefined ? [] : [[kind, limits[kind]]])),
    ) as Partial<Record<BudgetKind, number>>;
  }

  return {
    harnessId: reader.id,
    execution: executionView(execution),
    budget: budgetView(summaryOf(latest?.id ?? '')?.usage, ceilings ?? {}, latest?.tokens.credits),
    preamble: preamble ?? { guidance: '', guidanceTokens: 0, tools: [], toolsTokens: 0 },
    loaded,
    declared,
    turns: views,
  };
}

/** What the model calls of a turn add to what the reads carry on to the turns after it. */
interface Effects {
  readonly learned: Map<string, string>;
  readonly loads: Set<string>;
  firstLoads: readonly string[];
  preamble: PreambleView | undefined;
}

function executionView(execution: Execution): ExecutionView {
  return {
    seq: execution.seq,
    headId: execution.main,
    ...(execution.turn !== undefined ? { pendingStep: execution.turn.pending.kind } : {}),
    ...(execution.turn?.park !== undefined ? { parked: true } : {}),
  };
}

function milestone(e: HarnessEvent, origin: number): MilestoneView {
  const compacted = e.kind === 'step' ? e.compaction : undefined;
  const detail = compacted?.length
    ? compacted.map((c) => `compacted ${c.pointer} to ${c.sizeAfter.inputTokens} tokens`).join(', ')
    : undefined;
  return {
    position: e.seq,
    at: e.at - origin,
    kind: e.kind === 'step' ? `${e.stepKind} ${e.ending}` : e.control,
    ...(detail ? { detail } : {}),
  };
}

/** A turn's credits: its model calls' (nested runs' included); undefined when none is priced. */
function sumCredits(steps: readonly StepView[]): number | undefined {
  const priced = steps.flatMap((s) =>
    s.kind === 'model' && s.tokens?.credits !== undefined ? [s.tokens.credits] : [],
  );
  return priced.length ? priced.reduce((a, b) => a + b, 0) : undefined;
}

/** What a model call answered: its text and the calls it asked for, and how that ended. */
function modelOutput(
  result: Message | undefined,
  nameOf: (id: string) => string,
): Pick<ModelStepView, 'reply' | 'toolCalls' | 'finish'> {
  if (result === undefined) return {};
  const text = textParts(result).trim();
  const calls = result.parts.flatMap((p) => (p.kind === 'tool-call' ? [callText(nameOf(p.tool.id), p.arguments)] : []));
  const finish: ModelFinish = calls.length ? 'tool-calls' : 'complete';
  return { ...(text ? { reply: text } : {}), ...(calls.length ? { toolCalls: calls } : {}), finish };
}

/** One message of the model's input, as the Inspector shows it. Tokens are estimates (chars / 4). */
function describe(m: Message, part: WindowPart, nameOf: (id: string) => string): WindowEntryView {
  const text = messageText(m, nameOf);
  if (m.role === 'system') {
    const definitions = definitionsIn([m]);
    if (definitions.length > 0 && !m.parts.some((p) => p.kind === 'text'))
      return {
        role: 'tools',
        part,
        // A deferred one is loaded later where it is disclosed.
        title: definitions.map((d) => (d.deferred ? `${d.definition.name} (declared)` : d.definition.name)).join(', '),
        text,
        // A declared tool shows by its name only (what a provider sends at most).
        tokens: estimatedTokens(
          JSON.stringify(definitions.map((d) => (d.deferred ? d.definition.name : d.definition))),
        ),
      };
    return { role: 'system', part, title: 'Guidance', text, tokens: estimatedTokens(text) };
  }
  if (m.role === 'user' && m.parts.some((p) => p.kind === 'tool-result'))
    return {
      role: 'tool',
      part,
      title: 'tool results',
      text,
      tokens: estimatedTokens(text),
      ...(m.parts.some((p) => p.kind === 'tool-result' && p.status === 'failed') ? { failed: true } : {}),
    };
  if (m.role === 'assistant')
    return {
      role: 'assistant',
      part,
      title: m.parts.some((p) => p.kind === 'tool-call') ? 'Tool calls' : 'Reply',
      text,
      tokens: estimatedTokens(text),
    };
  return {
    role: 'user',
    part,
    title: part === 'state' ? 'State + request' : 'Request',
    text,
    tokens: estimatedTokens(text),
  };
}
