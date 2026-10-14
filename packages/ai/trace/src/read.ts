/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { canonicalJson, type Json, type JsonObject, type Message, type ModelFinish } from '@molstar/ai-agent';
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
} from '@molstar/ai-harness';
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
 * The agent trace (harness.md §11, §13; analysis/viewer-mode.md §5): the event log, heads and blobs of one harness as
 * data for the Inspector. Read-only and safe while a turn runs: it shows what has been committed. Everything that
 * belongs to a host (how a request is worded, what a tool call changed, prices) is a hook of the `TraceProfile`, with
 * a default.
 */

/** Opens the person's words in a request. The same text as `@molstar/ai-session`'s, which this package does not import. */
export const REQUEST_MARKER = 'Request: ';
/** Opens the Focus line of a request. The same text as `@molstar/ai-session`'s. */
export const FOCUS_MARKER = 'Focus: ';

/**
 * The labels (`EntryLabel.kind`) of the entries a session and its toolkit write, by the names `@molstar/ai-session`
 * and `@molstar/ai-tools` give them. A labelled entry is read by its label; one written before labels, by its text.
 */
export const LABELS = {
  request: 'session:request',
  guidance: 'session:guidance',
  guidanceUpdate: 'session:guidance-update',
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

/** The model-visible calls the trace caches between reads: heads, entries and blobs are immutable. */
export interface TraceCache {
  readonly entries: Map<string, Message>;
  readonly heads: Map<string, readonly string[]>;
  readonly blobs: Map<string, Json>;
}

export function createTraceCache(): TraceCache {
  return { entries: new Map(), heads: new Map(), blobs: new Map() };
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
  /** The clock the lease is measured against (default `Date.now`). */
  readonly now?: () => number;
  /** Kept across reads of one harness, so a re-read fetches only what is new. */
  readonly cache?: TraceCache;
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

/** A text's size in tokens, as the trace estimates it (about four characters a token). */
export const estimatedTokens = (text: string): number => Math.max(1, Math.round(text.length / 4));

/** A turn's state (the mapping is documented at `TurnStatus`). */
export function turnStatusOf(
  status: TurnSummary['status'] | undefined,
  execution: Pick<Execution, 'runner' | 'turn'>,
  turnId: string,
): TurnStatus {
  if (status === 'completed') return 'done';
  if (status === 'failed') return 'failed';
  if (status === 'stopped' || status === 'reset') return 'stopped';
  const turn = execution.turn?.id === turnId ? execution.turn : undefined;
  if (turn === undefined) return 'working';
  if (turn.park !== undefined || turn.blocked !== undefined || execution.runner.phase === 'none') return 'waiting';
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

/** The index of a marker that opens the text or a paragraph, or -1. */
function paragraphAt(text: string, marker: string): number {
  if (text.startsWith(marker)) return 0;
  const at = text.indexOf(`\n\n${marker}`);
  return at < 0 ? -1 : at + 2;
}

/** The rest of the line that starts with a marker, if any line does. */
function lineAfter(text: string, marker: string): string | undefined {
  for (const line of text.split('\n')) if (line.startsWith(marker)) return line.slice(marker.length).trim();
  return undefined;
}

/**
 * The request is the text after the request marker (the first that opens a paragraph, so a prompt that contains
 * the marker is kept whole), the focus the Focus line before it.
 */
export function defaultRequestOf(message: Message): TraceRequest {
  const data = message.label?.kind === LABELS.request ? message.label.data : undefined;
  if (typeof data?.prompt === 'string')
    return {
      prompt: data.prompt.trim() || '(no request recorded)',
      ...(typeof data.focus === 'string' && data.focus ? { focus: data.focus } : {}),
    };
  const text = textParts(message);
  const at = paragraphAt(text, REQUEST_MARKER);
  const prompt = (at >= 0 ? text.slice(at + REQUEST_MARKER.length) : text).trim() || '(no request recorded)';
  const focus = at > 0 ? lineAfter(text.slice(0, at), FOCUS_MARKER) : undefined;
  return { prompt, ...(focus ? { focus } : {}) };
}

/**
 * A system message of definitions only is the tools, other system text the guidance; a user message with something
 * before the request marker (the host's state, the Focus line, notes) is the state, any other message the
 * conversation.
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
  if (message.role === 'user' && !message.parts.some((p) => p.kind === 'tool-result')) {
    const text = textParts(message);
    const at = paragraphAt(text, REQUEST_MARKER);
    if (at > 0 && text.slice(0, at).trim() !== '') return 'state';
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
  const steps = turnEvents.filter((e): e is StepEvent => e.kind === 'step');
  return (
    steps.findLast((e) => e.limits !== undefined)?.limits ??
    steps.find((e) => e.begun)?.begun?.configuration.limits ??
    {}
  );
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

/** A request as the cache sees it: its parts (adapter-given digests) and its model. */
interface RequestShape {
  readonly segments: readonly string[];
  readonly model: string;
  readonly configuration: string;
}

interface Evidence {
  readonly head?: string | null;
  readonly configuration?: { readonly target?: { readonly model?: string } };
  readonly segments?: unknown;
}

function requestShape(evidence: Evidence | undefined): RequestShape | undefined {
  if (!Array.isArray(evidence?.segments)) return undefined;
  return {
    segments: evidence.segments as string[],
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
  const of = previous.segments.length;
  if (previous.model !== now.model) return { reason: `another model (${now.model})`, kept: 0, of };
  if (previous.configuration !== now.configuration) return { reason: "the model's settings changed", kept: 0, of };
  let kept = 0;
  while (kept < of && previous.segments[kept] === now.segments[kept]) kept++;
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

/**
 * Reads the trace of one harness: its turns with their model and tool calls (attempts, usage, what each model call
 * saw, where its prompt cache broke), the budget, the preamble and the tools loaded. Times are milliseconds from the
 * first step event.
 */
export async function readTrace(reader: HarnessReader, profile: TraceProfile = {}): Promise<TraceData> {
  const requestOf = profile.requestOf ?? defaultRequestOf;
  const partOf = profile.windowPart ?? defaultWindowPart;
  const cache = profile.cache ?? createTraceCache();
  const now = profile.now ?? Date.now;

  const execution = await reader.execution();
  const events: HarnessEvent[] = [];
  for (let page = 0; page < MAX_PAGES; page += 1) {
    const got = await reader.events({ after: events.at(-1)?.seq ?? 0, limit: PAGE });
    events.push(...got);
    if (got.length < PAGE) break;
  }
  // A turn's summary: from the log (the event that ended it), the active turn's from the reader.
  const summaries = new Map<string, TurnSummary>();
  for (const e of events) if (e.turnId !== undefined && e.summary !== undefined) summaries.set(e.turnId, e.summary);
  for (const t of await reader.turns({ limit: SUMMARIES })) summaries.set(t.id, t);
  const byTurn = new Map<string, HarnessEvent[]>();
  for (const e of events) {
    if (e.turnId === undefined) continue;
    const list = byTurn.get(e.turnId);
    if (list) list.push(e);
    else byTurn.set(e.turnId, [e]);
  }
  const origin = events.find((e) => e.kind === 'step')?.at ?? 0;

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
  const learn = (messages: readonly Message[]) => {
    for (const p of definitionsIn(messages)) names.set(p.definition.tool.id, p.definition.name);
  };
  const nameOf = (id: string) => profile.toolNames?.get(id) ?? names.get(id) ?? id;

  // A message's window entry is the same in every call that sees it; only `fresh` differs.
  const entryCache = new Map<string, WindowEntryView>();
  const windowEntry = (id: string, m: Message): WindowEntryView => {
    const known = entryCache.get(id);
    if (known) return known;
    const entry = describe(m, partOf(m), nameOf);
    entryCache.set(id, entry);
    return entry;
  };

  let preamble: PreambleView | undefined;
  // What each run's previous model call saw (the turn's: ""; a nested run's: its call's id).
  const previousWindow = new Map<string, Set<string>>();
  // Each run's previous request, and the latest of any run: a nested run's first request follows its parent's latest.
  const previousRequest = new Map<string, RequestShape>();
  let latestRequest: RequestShape | undefined;
  // The turn (1-based) whose model call first saw each tool loaded; the first call's are "from the start".
  const loadedFirst = new Map<string, number | undefined>();
  const views: TurnView[] = [];

  for (const [turnId, turnEvents] of byTurn) {
    const stepEvents = turnEvents.filter((e): e is StepEvent => e.kind === 'step');
    const begin = stepEvents.find((e) => e.stepKind === 'begin');
    // A call is recorded once in the step that ran it; a parked step's later event may repeat it. The calls come in the
    // order they began (a nested run's follow the call that ran them).
    const found = new Map<string, { call: CallRecord; at: number; order: number }>();
    for (const e of stepEvents)
      for (const call of e.calls) {
        const known = found.get(call.id);
        found.set(call.id, { call, at: e.at, order: known?.order ?? found.size });
      }
    const calls = [...found.values()]
      .sort(
        (a, b) =>
          (a.call.attempts[0]?.startedAt ?? a.at) - (b.call.attempts[0]?.startedAt ?? b.at) || a.order - b.order,
      )
      .map((c) => c.call);

    const steps: StepView[] = [];
    for (const call of calls) {
      const attempts = call.attempts.map((a, i) => attemptView(a, i + 1, origin, call.attempts[i - 1]));
      const under = call.parentCallId;
      if (call.kind === 'model') {
        const evidence = (await blob(call.request)) as Evidence | undefined;
        const head = await headOf(evidence?.head ?? undefined);
        learn(head.map((e) => e.message));
        const before = previousWindow.get(under ?? '') ?? new Set<string>();
        const entries = head.map((e) => {
          const entry = windowEntry(e.id, e.message);
          return before.has(e.id) ? entry : { ...entry, fresh: true };
        });
        previousWindow.set(under ?? '', new Set(head.map((e) => e.id)));
        if (under === undefined) {
          const first = preamble === undefined;
          preamble ??= preambleOf(head.map((e) => e.message));
          for (const p of definitionsIn(head.map((e) => e.message)))
            if (!p.deferred && !loadedFirst.has(p.definition.tool.id))
              loadedFirst.set(p.definition.tool.id, first ? undefined : views.length + 1);
        }
        const reply = (await blob(call.result)) as Message | undefined;
        const usage = call.attempts.at(-1)?.usage;
        const model = evidence?.configuration?.target?.model;
        const credits = usage && model ? profile.credits?.(model, usage) : undefined;
        const shape = requestShape(evidence);
        const prefixBreak = shape
          ? breakOf(previousRequest.get(under ?? '') ?? (under !== undefined ? latestRequest : undefined), shape)
          : undefined;
        if (shape) {
          previousRequest.set(under ?? '', shape);
          latestRequest = shape;
        }
        steps.push({
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
        } satisfies ModelStepView);
      } else steps.push({ ...(await toolStep(call, attempts)), ...(under !== undefined ? { under } : {}) });
    }

    const summary = summaries.get(turnId);
    const ended = turnEvents.findLast((e) => (e.kind === 'step' ? e.end !== undefined : e.summary !== undefined));
    const status = turnStatusOf(summary?.status, execution, turnId);
    const request = await requestFor(begin);
    views.push({
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
    });
  }

  async function toolStep(call: CallRecord, attempts: AttemptView[]): Promise<ToolStepView> {
    const request = (await blob(call.request)) as { toolId?: string; arguments?: Json } | undefined;
    const result = (await blob(call.result)) as
      | { status: 'succeeded' | 'failed'; parts: readonly { kind: string; text?: string }[] }
      | undefined;
    const last = call.attempts.at(-1);
    const text = result?.parts.map((p) => p.text ?? '').join('') ?? last?.failure?.message;
    const toolId = request?.toolId ?? call.key;
    const status: ToolStepView['status'] =
      result?.status ?? (last === undefined ? 'running' : last.outcome === 'aborted' ? 'aborted' : 'failed');
    const args = isObject(request?.arguments) ? request.arguments : {};
    const name = nameOf(toolId);
    const notes = profile.toolCall?.({
      callId: call.id,
      key: call.key,
      toolId,
      name,
      arguments: args,
      status,
      ...(text ? { text } : {}),
    });
    return {
      kind: 'tool',
      id: call.id,
      callId: call.id,
      toolId,
      name,
      args,
      status,
      attempts,
      ...(notes?.title ? { title: notes.title } : {}),
      ...(notes?.read !== undefined ? { read: notes.read } : {}),
      ...(text ? { result: text } : {}),
      ...(notes?.changed?.length ? { changed: notes.changed } : {}),
    };
  }

  /** The turn's request: the begin step's last user message. */
  async function requestFor(begin: StepEvent | undefined) {
    const move = begin?.moves.find((m) => m.pointer === MAIN) ?? begin?.moves[0];
    const head = await headOf(move?.to);
    const user = head.findLast((e) => e.message.role === 'user');
    return user ? requestOf(user.message) : { prompt: '(no request recorded)' };
  }

  // The latest context: what it has loaded and declares.
  const main = await headOf(execution.main);
  learn(main.map((e) => e.message));
  const mainDefinitions = definitionsIn(main.map((e) => e.message));
  const loadedIds = [...new Set(mainDefinitions.filter((p) => !p.deferred).map((p) => p.definition.tool.id))];
  const loaded: LoadedToolView[] = loadedIds.map((id) => {
    const turn = loadedFirst.get(id);
    return { id, name: nameOf(id), ...(turn !== undefined ? { turn } : {}) };
  });
  const declared = [
    ...new Set(
      mainDefinitions
        .filter((p) => p.deferred && !loadedIds.includes(p.definition.tool.id))
        .map((p) => p.definition.tool.id),
    ),
  ];

  const latest = views.at(-1);
  const ceilings =
    profile.ceilings ??
    (latest !== undefined
      ? (Object.fromEntries(
          KINDS.flatMap((kind) => {
            const limit = limitsOf(byTurn.get(latest.id) ?? [], execution, latest.id)[kind];
            return limit === undefined ? [] : [[kind, limit]];
          }),
        ) as Partial<Record<BudgetKind, number>>)
      : {});

  return {
    harnessId: reader.id,
    execution: executionView(execution, now),
    budget: budgetView(summaries.get(latest?.id ?? '')?.usage, ceilings, latest?.tokens.credits),
    preamble: preamble ?? { guidance: '', guidanceTokens: 0, tools: [], toolsTokens: 0 },
    loaded,
    declared,
    turns: views,
  };
}

function executionView(execution: Execution, now: () => number): ExecutionView {
  const runner = execution.runner;
  return {
    seq: execution.seq,
    headId: execution.main,
    runner: runner.phase,
    ...(runner.phase === 'live' && execution.turn !== undefined
      ? {
          lease: {
            worker: runner.runnerId,
            renewsInMs: Math.max(0, runner.leaseUntil - now()),
            leaseMs: runner.leaseUntil - runner.signalAt,
          },
        }
      : {}),
    ...(execution.turn !== undefined ? { pendingStep: execution.turn.pending.kind } : {}),
    ...(execution.turn?.park !== undefined ? { parked: true } : {}),
  };
}

function milestone(e: HarnessEvent, origin: number): MilestoneView {
  const compacted = e.kind === 'step' ? e.compaction : undefined;
  const detail = compacted?.length
    ? compacted.map((c) => `compacted ${c.pointer} to ${c.sizeAfter.inputTokens} tokens`).join(', ')
    : e.kind === 'step' && e.consumed.length
      ? `${e.consumed.length} items`
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
