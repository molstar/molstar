/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  Provider,
  ProviderBindings,
  ModelCapabilities,
  ModelConfiguration,
  ReasoningParameters,
  ToolDefinition,
} from '../messages/index.js';
import {
  type Clock,
  createHarness,
  type DriveContext,
  type DriveHandlers,
  emptyHead,
  type Harness,
  HarnessError,
  type HarnessOptions,
  type HarnessReader,
  type HeadEdit,
  type LastCompaction,
  MAIN,
  type Observation,
  type ToolBindings,
  type TurnConfiguration,
  type TurnOutcome,
  type TurnSummary,
  wallClock,
} from '../harness/index.js';
import {
  type AfterWrite,
  createToolBindings,
  type DisclosureApi,
  NestedTools,
  RepeatedReads,
  type Toolkit,
  type ToolRegistration,
} from '../tools/index.js';
import {
  type CompactionDomain,
  type CompactionSettings,
  compactionFor,
  runCompaction,
  toolkitDomain,
} from './compaction/index.js';
import { type CatchUp, guidanceLabel, LABEL_GUIDANCE, LABEL_GUIDANCE_UPDATE, turnInput } from './compose.js';
import { type BudgetPolicy, createDrive, type DriveHost, type DrivePolicy, turnBudget } from './drive.js';
import { newTurnState, type TurnState } from './state.js';
import type {
  ApprovalDecision,
  ModelChoice,
  PendingApproval,
  ResumeOptions,
  StatusEvent,
  ToolEvent,
  TurnHooks,
  TurnRequest,
  TurnResult,
} from './types.js';

/** All the model said in the turn: each top-level model call's reply, in order, separated by blank lines. */
const turnTextOf = (state: TurnState): string =>
  [...state.replies.values()]
    .map((text) => text.trim())
    .filter((text) => text !== '')
    .join('\n\n');

/** The session's thinking: one for every model, or per model. */
export type SessionThinking =
  | ModelChoice['thinking']
  | ((target: { readonly provider: string; readonly model: string }) => ReasoningParameters | undefined);

/** The reasoning parameters for a turn: the choice's, else the session's; "model" sends none. */
function reasoningFor(choice: ModelChoice, thinking: SessionThinking): { reasoning?: ReasoningParameters } {
  const chosen = choice.thinking ?? (typeof thinking === 'function' ? thinking(choice) : thinking);
  return chosen === undefined || chosen === 'model' ? {} : { reasoning: chosen };
}

/** Where the session finds models: a provider by its id (a lookup), and the model it runs by default. */
export interface ModelSource<M extends ModelChoice = ModelChoice> {
  /** The provider of an id (a target's `provider`); a lookup, so providers a person adds while the conversation is open are found. */
  providerFor(provider: string): Provider | undefined;
  default: M;
}

export interface AgentSessionOptions<Api extends DisclosureApi, M extends ModelChoice = ModelChoice>
  extends DrivePolicy {
  models: ModelSource<M>;
  /** The tools and their groups: the core is loaded from the first turn, the rest declared and loaded on demand. */
  toolkit: Toolkit<Api>;
  /** The host's modules; the session adds the `disclosure` module (what the conversation has loaded) when it has none. */
  api: Api;
  /**
   * The tools a conversation starts with: the core (default; the rest is disclosed on demand) or all of them (evals
   * compare the two).
   */
  startTools?: 'core' | 'all';
  /** System guidance (the first message of a conversation). */
  guidance: string;
  clock?: Clock;
  /**
   * Compaction (session.md §4): sizes to override, or false for none. Default: from each turn's model's
   * context window (`compactionFor`: at 80% of it).
   */
  compaction?: Partial<CompactionSettings> | false;
  /** What the host adds to compaction (the memory, pins, more wish sources); the rest comes from the toolkit. */
  compactionDomain?: Partial<CompactionDomain>;
  /**
   * Tool calls a turn may make before the drive stops it; default 200. The turn starts with at most 60 and is
   * extended (automatically, by that much again) until this ceiling.
   */
  maxToolCalls?: number;
  /**
   * The model's thinking when the turn's choice gives none, e.g. `{ mode: "adaptive", effort: "medium" }`, or per model
   * (a host passes the catalog's default: `defaultReasoning` of `providers/models.ts`); "model" or
   * undefined leaves it to the model.
   */
  thinking?: SessionThinking;
  /** Output tokens per model call, thinking included; default 16384. A reply cut short at the limit is asked for
   * once more with twice the limit. */
  maxOutputTokens?: number;
  /** The host's state version, for repeated reads (`RepeatedReads`): the same read while it is unchanged answers in a line. */
  version?: () => unknown;
  /** Runs after a call that changed something, before its result is recorded (where a host saves). */
  afterWrite?: AfterWrite;
  /** What a model call cost, in the host's unit (the harness keeps it in the usage). */
  price?: HarnessOptions['price'];
  /** The host's own decisions, over the session's policies (which are over `defaultDrive`); 'pass' leaves one to them. */
  drive?: Partial<DriveHandlers>;
  /** Heads the request's notes (default: "Notes since your last turn:"). */
  notesHeading?: string;
}

/** What `activeTurn` holds while a turn is being composed, before the harness gives its id. */
const COMPOSING = 'composing';

/**
 * One conversation over one in-memory harness (session.md §1–§2): each prompt is a turn of the harness, run
 * under this session's drive. One turn at a time; stop() cancels the running (or parked) turn.
 */
export class AgentSession<Api extends DisclosureApi, M extends ModelChoice = ModelChoice> {
  private readonly options: AgentSessionOptions<Api, M>;
  private readonly clock: Clock;
  /**
   * Tools by id, as the toolkit is now (a dynamic toolkit: decision 25 of the Prism viewer's spec,
   * `packages/prism/viewer/_spec/README.md`).
   */
  private registration(toolId: string): ToolRegistration<Api> | undefined {
    return this.options.toolkit.registration(toolId);
  }

  /** What compaction keeps, as the toolkit is now; the tools it no longer has are dropped (the eviction hint). */
  private get domain(): CompactionDomain {
    const { toolkit } = this.options;
    return { ...toolkitDomain(toolkit), evicted: () => toolkit.removed(), ...this.options.compactionDomain };
  }
  private made: Harness | undefined;
  /** Everything of the running turn, made afresh for each turn and resume. */
  private state: TurnState = newTurnState({});
  /** The running turn: `COMPOSING` until the harness has begun it, then its turn id. */
  private activeTurn: string | undefined;
  /** stop() came while the turn was being composed: it does not begin. */
  private stopAsked = false;
  /** Answers a repeated read briefly (per turn). */
  private readonly repeats: RepeatedReads | undefined;
  /** The person's answers to calls that needed approval, by call id (kept across resumes). */
  private readonly decided = new Map<string, ApprovalDecision>();
  /** No turn has ever run in this harness (or it was forgotten): the first turn carries guidance and tool definitions. */
  private fresh = true;

  /** Limits the providers gave, by model ("provider:model"): the one source of a model's window and output maximum. */
  private readonly known = new Map<string, ModelCapabilities>();
  /** The conversation's size after the last model call, in tokens, as the provider counted it. */
  private contextTokens: number | undefined;
  /** The running turn's model's output maximum (a widened retry stays within it). */
  private outputCap: number | undefined;

  constructor(options: AgentSessionOptions<Api, M>) {
    this.options = options;
    this.clock = options.clock ?? wallClock;
    // The disclosure tools are setting up, not reading: exempt by the toolkit's names for them.
    const { names } = options.toolkit;
    this.repeats = options.version
      ? new RepeatedReads({ version: options.version, exempt: [names.skills, names.tools] })
      : undefined;
  }

  get busy(): boolean {
    return this.activeTurn !== undefined;
  }

  private providerFor(provider: string): Provider | undefined {
    return this.options.models.providerFor(provider);
  }

  /**
   * The default model's context window in tokens: what its provider reported (its catalog, or the host's
   * `modelLimits` for a custom endpoint); 200k until the first call has asked.
   */
  contextWindow(): number {
    const { provider, model } = this.options.models.default;
    return this.windowOf(provider, model);
  }

  private windowOf(provider: string, model: string): number {
    return this.known.get(`${provider}:${model}`)?.contextWindowTokens ?? 200_000;
  }

  /** A model's limits: the provider's `capabilities()`, the one source (a model it does not know fails the turn). */
  private async capabilitiesOf(found: Provider, configuration: ModelConfiguration): Promise<ModelCapabilities> {
    const { provider, model } = configuration.target;
    const capabilities = await found.capabilities(configuration);
    this.known.set(`${provider}:${model}`, capabilities);
    return capabilities;
  }

  /** Compaction for a model's window, or none. */
  private compactionSettings(window: number): CompactionSettings | undefined {
    const o = this.options.compaction;
    return o === false ? undefined : compactionFor(window, o ?? {});
  }

  /** The turn's budget for a model (`turnBudget`), from its window's compaction trigger. */
  private budgetFor(configuration: ModelConfiguration): BudgetPolicy {
    const { provider, model } = configuration.target;
    const window = this.windowOf(provider, model);
    const trigger = this.compactionSettings(window)?.triggerInputTokens ?? Math.floor(window * 0.8);
    return turnBudget(trigger, this.options.maxToolCalls ?? 200);
  }

  /** The budget policy for a model (the default's, or a choice's): its initial limits and ceilings. */
  budget(choice?: M): BudgetPolicy {
    return this.budgetFor(this.modelConfiguration(choice));
  }

  /** The harness's reader, for reading the conversation (a trace, the Inspector). */
  client(): HarnessReader {
    return this.harness;
  }

  /** Whether the next turn is the conversation's first (no head yet): what a host's preamble depends on. */
  async isFresh(): Promise<boolean> {
    return this.fresh;
  }

  /** The messages of the conversation's committed context. */
  private async head(client: HarnessReader) {
    const execution = await client.execution();
    const [head] = await client.heads([execution.main]);
    return head === undefined ? [] : await client.entries(head.entries);
  }

  /** The tools a conversation starts with and the ones declared for later. */
  private toolsAtStart(): { starting: ToolDefinition[]; deferring: ToolDefinition[] } {
    const { toolkit, startTools } = this.options;
    const starting = startTools === 'all' ? toolkit.allDefinitions() : toolkit.coreDefinitions();
    // The rest are declared from the start, loaded when disclosed.
    const deferring = toolkit.allDefinitions().filter((d) => !starting.some((t) => t.tool.id === d.tool.id));
    return { starting, deferring };
  }

  /**
   * What the conversation has not seen: tools added or changed since its first turn (a dynamic toolkit changes between
   * turns; a changed definition replaces the old one): the starting tools and any it loaded later, loaded; the others
   * declared (deferred); and the guidance, when the context does not hold its current text. Undefined when nothing is
   * missing.
   */
  private async catchUp(
    client: HarnessReader,
    starting: readonly ToolDefinition[],
    deferring: readonly ToolDefinition[],
    guidance: string,
  ): Promise<CatchUp | undefined> {
    const seen = new Map<string, string>();
    const loaded = new Set<string>();
    const text = (d: ToolDefinition) => JSON.stringify([d.name, d.description, d.inputSchema]);
    // The current guidance is in the context when an entry labelled with its fingerprint is.
    const { key } = guidanceLabel(guidance);
    let current = false;
    for (const message of await this.head(client)) {
      const label = message.label;
      if ((label?.kind === LABEL_GUIDANCE || label?.kind === LABEL_GUIDANCE_UPDATE) && label.key === key)
        current = true;
      for (const part of message.parts) {
        if (part.kind === 'tool-definition') {
          seen.set(part.definition.tool.id, text(part.definition));
          if (!part.deferred) loaded.add(part.definition.tool.id);
        }
      }
    }
    const all = this.options.toolkit.allDefinitions();
    const changed = (d: ToolDefinition) => seen.get(d.tool.id) !== text(d);
    const isStarting = (d: ToolDefinition) => starting.some((s) => s.tool.id === d.tool.id);
    const tools = all.filter((d) => changed(d) && (isStarting(d) || loaded.has(d.tool.id)));
    const deferred = deferring.filter((d) => changed(d) && !isStarting(d) && !loaded.has(d.tool.id));
    return tools.length || deferred.length || !current ? { tools, deferred, guidance: !current } : undefined;
  }

  /** The harness, made on first use. */
  private get harness(): Harness {
    this.made ??= this.setUp();
    return this.made;
  }

  private setUp(): Harness {
    // Which tools and skills the conversation has: read from its committed context.
    const { toolkit } = this.options;
    const disclosure = {
      loaded: async () =>
        (await this.head(this.harness)).flatMap((m) =>
          m.parts.flatMap((p) => (p.kind === 'tool-definition' && !p.deferred ? [p.definition.tool.id] : [])),
        ),
      skills: async () => toolkit.groupsIn(await this.head(this.harness)),
    };
    const given = this.options.api;
    const api = given.disclosure ? given : ({ ...given, disclosure } as Api);
    const nested = new NestedTools<Api, M>((choice) => this.modelConfiguration(choice));
    const tools = this.observedTools(
      createToolBindings(toolkit, api, {
        ...(this.repeats ? { repeats: this.repeats } : {}),
        nested,
        ...(this.options.afterWrite ? { afterWrite: this.options.afterWrite } : {}),
      }),
    );
    const models: ProviderBindings = {
      resolve: async (target) => {
        const found = this.providerFor(target.provider);
        if (!found) throw new Error(`No provider for ${target.provider}.`);
        // Delegates rather than spreads: a provider made by a class keeps its methods on the prototype.
        return {
          capabilities: (configuration) => this.capabilitiesOf(found, configuration),
          project: (input) => found.project(input),
          invoke: (input) => found.invoke(input),
        };
      },
    };
    const harness = createHarness({
      models,
      tools,
      drive: createDrive(this.driveHost(), this.options.drive),
      clock: this.clock,
      ...(this.options.price ? { price: this.options.price } : {}),
    });
    harness.observe().subscribe((o) => this.observe(o));
    return harness;
  }

  /** What the drive reads and writes of this session. */
  private driveHost(): DriveHost<Api> {
    return {
      toolkit: this.options.toolkit,
      api: this.options.api,
      policy: this.options,
      decided: this.decided,
      turn: () => this.state,
      registration: (toolId) => this.registration(toolId),
      budgetFor: (configuration) => this.budgetFor(configuration),
      outputCap: () => this.outputCap,
      compaction: (context, when, sized) => this.compaction(context, when, sized),
      status: (event) => this.status(event),
      shown: (event) => this.shown(event),
    };
  }

  private toolSet(): Pick<TurnConfiguration, 'tools'> {
    // A tool marked approval runs only with the person's yes (the permission decision); the tools lazy skills will
    // bring are in it too (a turn's set is fixed when it begins).
    return { tools: this.options.toolkit.toolSet() };
  }

  // ------------------------------------------------------------------ what the person sees

  private status(event: StatusEvent): void {
    this.state.hooks.onStatus?.(event);
  }

  /** Tells the host the turn has begun (once per turn): where it takes its checkpoint. */
  private begin(turnId: string, resumed = false): void {
    if (this.state.announced === turnId) return;
    this.state.announced = turnId;
    this.state.hooks.onTurn?.({ kind: 'begun', turnId, ...(resumed ? { resumed } : {}) });
  }

  /** An event of a call in progress (started, progress): to the host, not kept. */
  private emit(event: ToolEvent): void {
    this.state.hooks.onTool?.(event);
  }

  /**
   * An event the person sees and the result keeps; once: a resumed step asks again about its denials and repairs,
   * which the person has seen.
   */
  private shown(event: Omit<ToolEvent, 'startedAt'> & { startedAt?: number }): void {
    const seen = this.state.events.some(
      (e) => e.id !== undefined && e.id === event.id && e.status === event.status && e.text === event.text,
    );
    if (seen) return;
    const now = this.clock.now();
    const done: ToolEvent = {
      ...event,
      startedAt: event.startedAt ?? (event.id !== undefined ? this.state.starts.get(event.id) : undefined) ?? now,
      endedAt: event.endedAt ?? now,
    };
    this.state.events.push(done);
    this.state.hooks.onTool?.(done);
  }

  /** Wraps tool handlers so the host sees each call start and finish. */
  private observedTools(inner: ToolBindings): ToolBindings {
    return {
      resolve: (definition) => {
        const handler = inner.resolve(definition);
        if (handler === undefined) return undefined;
        const registration = this.registration(definition.tool.id);
        const name = registration?.definition.name ?? definition.tool.id;
        const read = registration?.read === true;
        return {
          invoke: (context) => {
            // A helper's own calls (a nested run under a tool's call) are the trace's, not the person's.
            if (context.path.length > 0) return handler.invoke(context);
            const session = this;
            return (async function* () {
              const id = context.callId;
              const args = context.arguments;
              // The host's checkpoint comes before any change a call makes, whatever reached it first.
              session.begin(context.turnId);
              const { starts } = session.state;
              if (!starts.has(id)) starts.set(id, session.clock.now());
              const startedAt = starts.get(id) ?? session.clock.now();
              session.emit({ id, name, args, status: 'started', read, startedAt });
              for await (const frame of handler.invoke(context)) {
                if (frame.kind === 'progress')
                  session.emit({ id, name, args, status: 'progress', progress: String(frame.value), read, startedAt });
                if (frame.kind === 'result') {
                  const text = frame.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');
                  const done: ToolEvent = { id, name, args, status: frame.status, text, read, startedAt };
                  session.state.retrying.delete(id);
                  // A transient failure may be retried; it shows only if the retries run out.
                  if (frame.status === 'failed' && frame.transient === true) session.state.retrying.set(id, done);
                  else session.shown(done);
                }
                yield frame;
              }
            })();
          },
        };
      },
    };
  }

  private observe(o: Observation): void {
    const state = this.state;
    if (o.kind === 'model.delta' && 'text' in o.delta && o.path.length === 0) {
      const text = (state.previews.get(o.attemptId) ?? '') + o.delta.text;
      state.previews.set(o.attemptId, text);
      state.replies.set(o.callId, text);
      // Whether the text says anything yet, from the delta alone: `trim()` of the whole text would flatten its rope
      // at every delta (an answer's length squared), and a text with something in it keeps it as it grows.
      if (!state.said.has(o.attemptId) && o.delta.text.trim()) state.said.add(o.attemptId);
      if (state.said.has(o.attemptId)) {
        state.lastText = text;
        state.hooks.onText?.(text);
      }
    } else if (o.kind === 'decision.opened' && o.path.length === 0) {
      if (o.decisionKind === 'permission' || o.decisionKind === 'budget') this.status({ kind: o.decisionKind });
    } else if (o.kind === 'decision.answered' && o.path.length === 0) {
      // A reply cut short and asked for again: the new call's reply replaces it.
      if (o.decisionKind === 'model-result' && o.answer === 'replace' && state.lastModelCall !== undefined)
        state.replies.delete(state.lastModelCall);
    } else if (o.kind === 'attempt.started' && o.path.length === 0) {
      const tools = state.planned.some((c) => c.callId === o.callId);
      if (!tools) {
        // A retried attempt replaces what the failed one said.
        state.lastModelCall = o.callId;
        state.replies.set(o.callId, '');
      }
      this.status({ kind: tools ? 'tools' : 'model' });
    } else if (o.kind === 'event' && o.event.kind === 'step') {
      const step = o.event;
      for (const call of step.calls)
        if (call.kind === 'model' && call.parentCallId === undefined) {
          const u = call.attempts.at(-1)?.usage;
          if (u) this.contextTokens = u.uncachedInput + u.cacheRead + u.cacheWrite + u.output;
        }
      if (step.stepKind === 'tools' && step.ending === 'settled') {
        // The calls that never ran (after a failed edit): skipped, apart from the failure that stopped them (the
        // model reads the skip reason with each).
        for (const c of state.planned)
          if (!state.events.some((e) => e.id === c.callId))
            this.shown({
              id: c.callId,
              name: this.registration(c.toolId)?.definition.name ?? c.toolId,
              args: c.arguments,
              status: 'skipped',
              text: state.skipShown ?? 'Not run.',
              read: this.registration(c.toolId)?.read === true,
            });
        state.planned = [];
        state.skipShown = undefined;
      }
      if (step.next?.kind === 'tools')
        state.planned = step.next.plan.map((c) => ({ callId: c.callId, toolId: c.toolId, arguments: c.arguments }));
    }
  }

  /** A compaction's edit when the context needs one (runCompaction), else undefined. */
  private async compaction(
    context: DriveContext,
    when: 'force' | 'trigger' | 'soft',
    sized?: { size: { inputTokens: number }; lastCompaction?: LastCompaction },
  ): Promise<HeadEdit | undefined> {
    const model = context.configuration?.model ?? this.modelConfiguration(undefined);
    const settings = this.compactionSettings(this.windowOf(model.target.provider, model.target.model));
    if (settings === undefined) return undefined;
    const state = this.state;
    return await runCompaction(context, {
      settings,
      domain: this.domain,
      model,
      when,
      ...(sized !== undefined ? { size: sized.size.inputTokens } : {}),
      ...(sized?.lastCompaction !== undefined || state.compactedSize !== undefined
        ? { lastSize: Math.max(sized?.lastCompaction?.sizeAfter.inputTokens ?? 0, state.compactedSize ?? 0) }
        : {}),
      onPlan: (plan) => {
        state.compactions.push(plan.report);
        state.compactedSize = plan.report.estimatedTokens;
        this.status({ kind: 'compacting' });
      },
    });
  }

  // ------------------------------------------------------------------ turns

  /** Runs one turn to completion (or stop/failure/park). Throws only on programming errors. */
  async turn(request: TurnRequest<M>, hooks: TurnHooks = {}): Promise<TurnResult> {
    if (this.activeTurn) throw new Error('A turn is already running.');
    this.activeTurn = COMPOSING;
    this.stopAsked = false;
    this.reset(hooks);
    try {
      const harness = this.harness;
      let configuration: TurnConfiguration;
      try {
        configuration = await this.configurationFor(request.model);
      } catch (e) {
        return {
          status: 'failed',
          text: '',
          turnText: '',
          error: e instanceof Error ? e.message : String(e),
          tools: [],
        };
      }
      const { starting, deferring } = this.toolsAtStart();
      const { guidance } = this.options;
      const catchUp = this.fresh ? undefined : await this.catchUp(harness, starting, deferring, guidance);
      const messages = turnInput({
        prompt: request.prompt,
        ...(request.state !== undefined ? { state: request.state } : {}),
        ...(request.focus !== undefined ? { focus: request.focus } : {}),
        ...(request.notes?.length ? { notes: request.notes } : {}),
        ...(request.preamble ? { preamble: request.preamble } : {}),
        ...(request.attachments ? { attachments: request.attachments } : {}),
        ...(this.options.notesHeading ? { notesHeading: this.options.notesHeading } : {}),
        first: this.fresh,
        tools: starting,
        deferred: deferring,
        guidance,
        ...(catchUp ? { catchUp } : {}),
      });
      if (this.stopAsked) return { status: 'stopped', text: '', turnText: '', tools: [] };
      const handle = harness.turn({ messages, configuration });
      this.fresh = false;
      this.activeTurn = handle.turnId;
      // Before the turn's first model call, so before any change it makes.
      this.begin(handle.turnId);
      return await this.outcome(await handle.done);
    } catch (e) {
      return {
        status: 'failed',
        text: this.state.lastText,
        turnText: turnTextOf(this.state),
        error: e instanceof Error ? e.message : String(e),
        tools: this.state.events,
      };
    } finally {
      this.activeTurn = undefined;
    }
  }

  private reset(hooks: TurnHooks): void {
    this.state = newTurnState(hooks);
    this.repeats?.reset();
  }

  /** The session's work is over: a turn still running or parked ends as stopped. */
  dispose(): void {
    this.made?.stop();
  }

  /**
   * Back to an empty conversation: `main` returns to the empty head, so the next turn is a first turn again
   * (guidance, tools, then its request). A parked turn ends as reset. Throws while a turn runs.
   */
  async forget(): Promise<void> {
    if (this.activeTurn) throw new Error('A turn is running: stop it first.');
    this.harness.reset({ [MAIN]: emptyHead().id });
    this.fresh = true;
    this.decided.clear();
    this.contextTokens = undefined;
    this.repeats?.reset();
  }

  /**
   * Resumes a turn that was parked for approval, with the person's decisions. Returns undefined when a turn runs or
   * none is parked.
   */
  async resume({ decisions = {}, hooks = {} }: ResumeOptions = {}): Promise<TurnResult | undefined> {
    if (this.activeTurn) return undefined;
    // The parked step's calls: its re-run reports their attempts as tools, and the ones that never run as skipped.
    const pending = (await this.harness.execution()).turn?.pending;
    if (this.activeTurn) return undefined;
    let handle: ReturnType<Harness['resume']>;
    try {
      handle = this.harness.resume();
    } catch (e) {
      if (e instanceof HarnessError) return undefined;
      throw e;
    }
    // The loop starts after this tick: the decisions and the hooks are in place for its first decision.
    for (const [callId, d] of Object.entries(decisions)) this.decided.set(callId, d);
    // What the person saw before the park stays the turn's (a call that ended then is not run again).
    const { events, starts } = this.state;
    this.reset(hooks);
    this.state.events = [...events];
    for (const [id, at] of starts) this.state.starts.set(id, at);
    if (pending?.kind === 'tools')
      this.state.planned = pending.plan.map((c) => ({ callId: c.callId, toolId: c.toolId, arguments: c.arguments }));
    this.activeTurn = handle.turnId;
    try {
      this.begin(handle.turnId, true);
      return await this.outcome(await handle.done);
    } catch (e) {
      return {
        status: 'failed',
        text: this.state.lastText,
        turnText: turnTextOf(this.state),
        error: e instanceof Error ? e.message : String(e),
        tools: this.state.events,
      };
    } finally {
      this.activeTurn = undefined;
    }
  }

  /** The result of a turn that ended or parked. */
  private async outcome(outcome: TurnOutcome): Promise<TurnResult> {
    const state = this.state;
    const client = this.harness;
    if ('parked' in outcome)
      return {
        status: 'awaiting-approval',
        approval: await this.awaiting(),
        text: state.lastText,
        turnText: turnTextOf(state),
        tools: state.events,
      };
    const summary = (await client.turns({ limit: 16 })).find((t) => t.id === outcome.turnId);
    const status = outcome.end.status;
    const error =
      status === 'failed' ? ((summary && (await this.endReason(client, summary))) ?? 'The turn failed.') : undefined;
    if (!state.lastText && status === 'completed') {
      const last = (await this.head(client)).findLast((m) => m.role === 'assistant');
      state.lastText = last?.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('') ?? '';
    }
    const turnText = turnTextOf(state) || state.lastText;
    const u = summary?.usage;
    return {
      status,
      text: state.lastText,
      turnText,
      ...(error !== undefined ? { error } : {}),
      tools: state.events,
      ...(u !== undefined
        ? {
            usage: {
              input: u.uncachedInput + u.cacheRead + u.cacheWrite,
              output: u.output,
              cached: u.cacheRead,
              written: u.cacheWrite,
              modelCalls: u.modelCalls,
            },
          }
        : {}),
      ...(this.contextTokens !== undefined ? { context: this.contextTokens } : {}),
      ...(state.compactions.length ? { compactions: [...state.compactions] } : {}),
    };
  }

  /** Why a turn ended, in words, from the event that ended it. */
  private async endReason(client: HarnessReader, summary: TurnSummary): Promise<string | undefined> {
    if (summary.lastStep === undefined) return undefined;
    const ending = (await client.steps([summary.lastStep])).findLast((e) => e.end !== undefined);
    const reason = ending?.end?.reason as { reason?: string; message?: string; kind?: string } | undefined;
    if (reason === undefined) return undefined;
    return reason.message ?? (reason.kind !== undefined ? `${reason.reason}: ${reason.kind}` : reason.reason);
  }

  /** The calls a parked turn waits on: its pending calls that need approval and are not decided. */
  private async awaiting(): Promise<PendingApproval[]> {
    const pending = (await this.harness.execution()).turn?.pending;
    if (pending?.kind !== 'tools') return [];
    const set = this.options.toolkit.toolSet();
    const asked = new Set(set.enabled.filter((id) => !set.preauthorized.includes(id)));
    return pending.plan
      .filter((c) => asked.has(c.toolId) && !this.decided.has(c.callId))
      .map((c) => ({ callId: c.callId, tool: c.toolId, args: c.arguments }));
  }

  /**
   * A model's configuration: the output the session asks for, clamped to the model's maximum once its provider has
   * given it (a nested run on the turn's model runs within it too).
   */
  private modelConfiguration(choice: M | undefined): ModelConfiguration {
    const chosen = choice ?? this.options.models.default;
    const { provider, model } = chosen;
    const wanted = this.options.maxOutputTokens ?? 16_384;
    const cap = this.known.get(`${provider}:${model}`)?.maxOutputTokens;
    return {
      target: { provider, model },
      parameters: {
        maxOutputTokens: cap !== undefined ? Math.min(wanted, cap) : wanted,
        ...reasoningFor(chosen, this.options.thinking),
      },
    };
  }

  /**
   * The configuration for this turn's model (the session's own, or the one chosen): its limits resolved, the output it
   * asks for clamped to the model's maximum, its budget. The limits are asked without an output: a provider checks the
   * output it is given against the model's maximum, and the clamp is what keeps it within.
   */
  private async configurationFor(choice: M | undefined): Promise<TurnConfiguration> {
    const chosen = choice ?? this.options.models.default;
    const found = this.providerFor(chosen.provider);
    if (!found) throw new Error(`No provider for ${chosen.provider}.`);
    const { maxOutputTokens: _, ...unsized } = this.modelConfiguration(choice).parameters;
    const capabilities = await this.capabilitiesOf(found, {
      target: { provider: chosen.provider, model: chosen.model },
      parameters: unsized,
    });
    this.outputCap = capabilities.maxOutputTokens;
    const configured = this.modelConfiguration(choice);
    return { ...this.toolSet(), model: configured, limits: this.budgetFor(configured).initial };
  }

  /** Cancels the running (or parked) turn. Changes already applied stay (the host can undo the turn). */
  async stop(): Promise<void> {
    if (this.activeTurn === COMPOSING) this.stopAsked = true;
    else this.made?.stop();
  }
}
