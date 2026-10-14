/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  Json,
  ModelAdapter,
  ModelBindings,
  ModelCapabilities,
  ModelConfiguration,
  ReasoningParameters,
  ToolDefinition,
} from '@molstar/ai-agent';
import { createInProcessSupervisor, createMemoryStorage } from '@molstar/ai-components';
import {
  type Clock,
  createHarnessHost,
  createRunner,
  type DriveContext,
  type DriveHandlers,
  type Execution,
  emptyHead,
  type HarnessClient,
  type HarnessStorage,
  type HeadEdit,
  type Item,
  type LastCompaction,
  MAIN,
  type Observation,
  type RunnerOptions,
  type ToolBindings,
  type TurnConfiguration,
  type TurnSummary,
} from '@molstar/ai-harness';
import {
  type AfterWrite,
  createToolBindings,
  type DisclosureApi,
  NestedTools,
  RepeatedReads,
  type Toolkit,
  type ToolRegistration,
} from '@molstar/ai-tools';
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
  StatusEvent,
  ToolEvent,
  TurnHooks,
  TurnRequest,
  TurnResult,
} from './types.js';

/** The session's thinking: one for every model, or per model. */
export type SessionThinking =
  | ModelChoice['thinking']
  | ((target: { readonly provider: string; readonly model: string }) => ReasoningParameters | undefined);

/** The reasoning parameters for a turn: the choice's, else the session's; "model" sends none. */
function reasoningFor(choice: ModelChoice, thinking: SessionThinking): { reasoning?: ReasoningParameters } {
  const chosen = choice.thinking ?? (typeof thinking === 'function' ? thinking(choice) : thinking);
  return chosen === undefined || chosen === 'model' ? {} : { reasoning: chosen };
}

/** Where the session finds models: an adapter per provider (a lookup), and the model it runs by default. */
export interface ModelSource<M extends ModelChoice = ModelChoice> {
  /** The adapter for a provider; a lookup, so providers a person adds while the conversation is open are found. */
  adapterFor(provider: string): ModelAdapter | undefined;
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
  /** Where the harness keeps its records; default: in memory only. Pass one over IndexedDB to resume across
   * reloads and restarts, or bring your own `HarnessStorage`. Several conversations may share one. */
  storage?: HarnessStorage;
  /** Stable harness id, so durable storage finds this conversation's history. */
  harnessId?: string;
  clock?: Clock;
  /** Ownership lease; a reopened tab waits at most this long to resume a turn. */
  leaseMs?: number;
  /** The default model's context window in tokens; default: what its adapter reports. */
  contextWindow?: number;
  /**
   * Compaction (tools-and-session.md §6): sizes to override, or false for none. Default: from each turn's model's
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
   * (a host passes the catalog's default: `defaultReasoning` of `@molstar/ai-model-adapters/model-catalog`); "model" or
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
  price?: RunnerOptions['price'];
  /** The host's own decisions, over the session's policies (which are over `defaultDrive`); 'pass' leaves one to them. */
  drive?: Partial<DriveHandlers>;
  /** Heads the request's notes (default: "Notes since your last turn:"). */
  notesHeading?: string;
}

interface Harness {
  readonly id: string;
  readonly client: HarnessClient;
  readonly supervisor: ReturnType<typeof createInProcessSupervisor>;
}

const wallClock: Clock = {
  now: () => Date.now(),
  sleep: (ms, signal) =>
    new Promise<void>((resolve, reject) => {
      if (signal?.aborted) return reject(signal.reason);
      const timer = setTimeout(resolve, ms);
      signal?.addEventListener(
        'abort',
        () => {
          clearTimeout(timer);
          reject(signal.reason);
        },
        { once: true },
      );
    }),
};

const uuid = () => globalThis.crypto.randomUUID();

/** What `activeTurn` holds while a resume or a forget reads the harness, before a turn id is known (or for none). */
const RESUMING = 'resuming';
const FORGETTING = 'forgetting';
const claimed = (active: string | undefined) => active === RESUMING || active === FORGETTING;

/**
 * One conversation over one harness (tools-and-session.md §4). Prompts go to its harness as items; a runner (started
 * by an in-process supervisor) works each turn under this session's drive. One turn at a time; stop() cancels the
 * running turn.
 */
export class AgentSession<Api extends DisclosureApi, M extends ModelChoice = ModelChoice> {
  private readonly options: AgentSessionOptions<Api, M>;
  private readonly clock: Clock;
  private readonly registrations: ReadonlyMap<string, ToolRegistration<Api>>;
  private readonly domain: CompactionDomain;
  private harness: Harness | undefined;
  /** Everything of the running turn, made afresh for each turn and resume. */
  private state: TurnState = newTurnState({});
  /** The running turn: its prompt's key until it has begun, then its turn id. */
  private activeTurn: string | undefined;
  /** Answers a repeated read briefly (per turn). */
  private readonly repeats: RepeatedReads | undefined;
  /** The person's answers to calls that needed approval, by call id (kept across resumes). */
  private readonly decided = new Map<string, ApprovalDecision>();
  /** No turn has ever run in this harness (or it was forgotten): the first turn carries guidance and tool definitions. */
  private fresh = true;

  /** Context windows by model ("provider:model"), as the turns' choices gave them. */
  private readonly windows = new Map<string, number>();
  /** Limits the adapters gave, by model. */
  private readonly known = new Map<string, ModelCapabilities>();
  /** The conversation's size after the last model call, in tokens, as the provider counted it. */
  private contextTokens: number | undefined;
  /** The running turn's model's output maximum (a widened retry stays within it). */
  private outputCap: number | undefined;
  /** Output limits by model ("provider:model"), for models the adapters do not know. */
  private readonly outputs = new Map<string, number>();

  constructor(options: AgentSessionOptions<Api, M>) {
    this.options = options;
    this.clock = options.clock ?? wallClock;
    this.registrations = new Map(options.toolkit.tools.map((r) => [r.definition.tool.id, r]));
    this.domain = { ...toolkitDomain(options.toolkit), ...options.compactionDomain };
    this.repeats = options.version ? new RepeatedReads({ version: options.version }) : undefined;
    this.remember(options.models.default);
  }

  get busy(): boolean {
    return this.activeTurn !== undefined;
  }

  /** The window and output limit a choice gives, kept for the model it names. */
  private remember(choice: ModelChoice): void {
    const key = `${choice.provider}:${choice.model}`;
    if (choice.contextWindow) this.windows.set(key, choice.contextWindow);
    if (choice.maxOutput) this.outputs.set(key, choice.maxOutput);
  }

  private adapterFor(provider: string): ModelAdapter | undefined {
    return this.options.models.adapterFor(provider);
  }

  /** A model's context window as the application gives it (its catalog, a custom model's), if it does. */
  private givenWindow(provider: string, model: string): number | undefined {
    const preferred = this.options.models.default;
    if (provider === preferred.provider && model === preferred.model && this.options.contextWindow)
      return this.options.contextWindow;
    return this.windows.get(`${provider}:${model}`);
  }

  /** The default model's context window in tokens (the application's, else the adapter's last word). */
  contextWindow(): number {
    const { provider, model } = this.options.models.default;
    return this.windowOf(provider, model);
  }

  private windowOf(provider: string, model: string): number {
    return this.givenWindow(provider, model) ?? this.known.get(`${provider}:${model}`)?.contextWindowTokens ?? 200_000;
  }

  /**
   * A model's limits: the adapter's for a model it knows, with the application's window when it gives one; for a
   * model the adapter does not know (a custom one), the application's window and output limit.
   */
  private async capabilitiesOf(adapter: ModelAdapter, configuration: ModelConfiguration): Promise<ModelCapabilities> {
    const { provider, model } = configuration.target;
    const window = this.givenWindow(provider, model);
    let known: ModelCapabilities | undefined;
    try {
      known = await adapter.capabilitiesFor(configuration);
      this.known.set(`${provider}:${model}`, known);
    } catch (e) {
      if (window === undefined)
        throw new Error(
          `The adapter does not know ${provider} model "${model}": give its context window in the model options. (${e instanceof Error ? e.message : String(e)})`,
        );
    }
    if (known && window === undefined) return known;
    if (known && window !== undefined)
      return {
        ...known,
        contextWindowTokens: window,
        ...(known.maxInputTokens !== undefined ? { maxInputTokens: Math.min(known.maxInputTokens, window) } : {}),
      };
    return {
      contextWindowTokens: window ?? 200_000,
      maxOutputTokens: this.outputs.get(`${provider}:${model}`) ?? 32_768,
      outputLimitIncludesReasoning: true,
      appendsTools: false,
      cacheScope: false,
    };
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

  /** The harness's client, for reading the conversation (a trace). */
  async client(): Promise<HarnessClient> {
    return (await this.ensure()).client;
  }

  /** Whether the next turn is the conversation's first (no head yet): what a host's preamble depends on. */
  async isFresh(): Promise<boolean> {
    await this.ensure();
    return this.fresh;
  }

  /** The messages of the conversation's committed context. */
  private async head(client: HarnessClient) {
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
   * What a conversation that began earlier has not seen: tools added or changed since (a conversation keeps its
   * harness across releases; a changed definition replaces the old one): the starting tools and any it loaded later,
   * loaded; the others declared (deferred); and the guidance, when its text changed. Undefined when nothing is
   * missing.
   */
  private async catchUp(
    client: HarnessClient,
    starting: readonly ToolDefinition[],
    deferring: readonly ToolDefinition[],
    guidance: string,
  ): Promise<CatchUp | undefined> {
    const seen = new Map<string, string>();
    const loaded = new Set<string>();
    const text = (d: ToolDefinition) => JSON.stringify([d.name, d.description, d.inputSchema]);
    // The current guidance is in the context when an entry labelled with its fingerprint is (one written before labels:
    // when a text ends with it).
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
        } else if (label === undefined && part.kind === 'text' && part.text.endsWith(guidance)) current = true;
      }
    }
    const all = this.options.toolkit.allDefinitions();
    const changed = (d: ToolDefinition) => seen.get(d.tool.id) !== text(d);
    const isStarting = (d: ToolDefinition) => starting.some((s) => s.tool.id === d.tool.id);
    const tools = all.filter((d) => changed(d) && (isStarting(d) || loaded.has(d.tool.id)));
    const deferred = deferring.filter((d) => changed(d) && !isStarting(d) && !loaded.has(d.tool.id));
    return tools.length || deferred.length || !current ? { tools, deferred, guidance: !current } : undefined;
  }

  private async ensure(): Promise<Harness> {
    if (this.harness) return this.harness;
    const { clock } = this;
    const storage = this.options.storage ?? createMemoryStorage();
    const leaseMs = this.options.leaseMs ?? 10_000;
    const host = createHarnessHost(storage, { clock, leaseMs, launchMs: leaseMs });
    const id = this.options.harnessId ?? `session-${uuid()}`;
    await host.create(id, { key: `create-${id}` });
    const client = host.client(id);
    this.fresh = (await this.head(client)).length === 0;
    // Which tools and skills the conversation has: read from its committed context.
    const { toolkit } = this.options;
    const disclosure = {
      loaded: async () =>
        (await this.head(client)).flatMap((m) =>
          m.parts.flatMap((p) => (p.kind === 'tool-definition' && !p.deferred ? [p.definition.tool.id] : [])),
        ),
      skills: async () => toolkit.groupsIn(await this.head(client)),
    };
    const given = this.options.api;
    const api = given.disclosure ? given : ({ ...given, disclosure } as Api);
    const nested = new NestedTools<Api, M>((choice) => this.modelConfiguration(choice));
    const tools = this.observedTools(
      createToolBindings(toolkit.tools, api, {
        ...(this.repeats ? { repeats: this.repeats } : {}),
        nested,
        ...(this.options.afterWrite ? { afterWrite: this.options.afterWrite } : {}),
      }),
    );
    const models: ModelBindings = {
      resolve: async (target) => {
        const found = this.adapterFor(target.provider);
        if (!found) throw new Error(`No adapter for ${target.provider}.`);
        return { ...found, capabilitiesFor: (configuration) => this.capabilitiesOf(found, configuration) };
      },
    };
    const port = host.runnerPort(id);
    const drive = createDrive(this.driveHost(), this.options.drive);
    const supervisor = createInProcessSupervisor({
      client,
      clock,
      startRunner: () => {
        const runner = createRunner({
          port,
          runnerId: `agent-session-${uuid()}`,
          models,
          tools,
          drive,
          configure: (items) => this.configure(items),
          clock,
          ...(this.options.price ? { price: this.options.price } : {}),
        });
        runner.observe().subscribe((o) => this.observe(o));
        return runner;
      },
    });
    this.harness = { id, client, supervisor };
    return this.harness;
  }

  /** What the drive reads and writes of this session. */
  private driveHost(): DriveHost<Api> {
    return {
      toolkit: this.options.toolkit,
      api: this.options.api,
      policy: this.options,
      decided: this.decided,
      turn: () => this.state,
      registration: (toolId) => this.registrations.get(toolId),
      budgetFor: (configuration) => this.budgetFor(configuration),
      outputCap: () => this.outputCap,
      compaction: (context, when, sized) => this.compaction(context, when, sized),
      status: (event) => this.status(event),
      shown: (event) => this.shown(event),
    };
  }

  /** A turn's configuration, as its item carries it (else the session's own). */
  private configure(items: readonly Item[]): TurnConfiguration {
    const given = (items[0]?.options as { configuration?: TurnConfiguration } | undefined)?.configuration;
    if (given !== undefined) return given;
    const model = this.modelConfiguration(undefined);
    return { ...this.toolSet(), model, limits: this.budgetFor(model).initial };
  }

  private toolSet(): Pick<TurnConfiguration, 'tools'> {
    const all = [...this.registrations.values()];
    // A tool marked approval runs only with the person's yes (the permission decision).
    return {
      tools: {
        enabled: all.map((r) => r.definition.tool.id),
        preauthorized: all.filter((r) => !r.approval).map((r) => r.definition.tool.id),
      },
    };
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

  /** An event the person sees and the result keeps. */
  private shown(event: Omit<ToolEvent, 'startedAt'> & { startedAt?: number }): void {
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
        const registration = this.registrations.get(definition.tool.id);
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
      if (text.trim()) {
        state.lastText = text;
        state.hooks.onText?.(text);
      }
    } else if (o.kind === 'decision.opened' && o.path.length === 0) {
      if (o.decisionKind === 'permission' || o.decisionKind === 'budget') this.status({ kind: o.decisionKind });
    } else if (o.kind === 'attempt.started' && o.path.length === 0) {
      this.status({ kind: state.planned.some((c) => c.callId === o.callId) ? 'tools' : 'model' });
    } else if (o.kind === 'event' && o.event.kind === 'step') {
      const step = o.event;
      // The step that begins this prompt's turn: its id is known, before the first model call.
      if (state.key !== undefined && step.begun?.items.includes(state.key)) this.begin(step.turnId);
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
              name: this.registrations.get(c.toolId)?.definition.name ?? c.toolId,
              args: c.arguments,
              status: 'skipped',
              text: state.skipShown ?? 'Not run.',
              read: this.registrations.get(c.toolId)?.read === true,
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
    // At a turn's start, the coming turn's model (its prompt carries the configuration).
    const coming = (context.inbox[0]?.options as { configuration?: TurnConfiguration } | undefined)?.configuration;
    const model = context.configuration?.model ?? coming?.model ?? this.modelConfiguration(undefined);
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
    const key = `prompt-${Date.now().toString(36)}-${uuid()}`;
    this.activeTurn = key;
    this.reset(hooks, key);
    try {
      const harness = await this.ensure();
      let configuration: TurnConfiguration;
      try {
        configuration = await this.configurationFor(request.model);
      } catch (e) {
        return { status: 'failed', text: '', error: e instanceof Error ? e.message : String(e), tools: [] };
      }
      const { starting, deferring } = this.toolsAtStart();
      const { guidance } = this.options;
      const catchUp = this.fresh ? undefined : await this.catchUp(harness.client, starting, deferring, guidance);
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
      const submitted = await harness.client.submit({
        key,
        messages,
        options: { configuration } as unknown as Json,
      });
      this.fresh = false;
      if (submitted.launch !== undefined) harness.supervisor.ensure(harness.id, submitted.launch);
      return await this.outcome(harness, { key });
    } catch (e) {
      return {
        status: 'failed',
        text: this.state.lastText,
        error: e instanceof Error ? e.message : String(e),
        tools: this.state.events,
      };
    } finally {
      this.activeTurn = undefined;
    }
  }

  private reset(hooks: TurnHooks, key?: string): void {
    this.state = newTurnState(hooks, key);
    this.repeats?.reset();
  }

  /**
   * Stops looking after this conversation's turns (its supervisor): the session it served is gone, as a closed tab's
   * is, and another (a reload, the server's next session) takes them up. A runner at work finishes its step; one
   * whose lease runs out is relaunched by whoever supervises next.
   */
  dispose(): void {
    this.harness?.supervisor.close();
  }

  /** Whether the harness holds a turn that was interrupted or parked (resume() would continue it). */
  async interrupted(): Promise<boolean> {
    if (this.activeTurn) return false;
    const { client } = await this.ensure();
    return (await client.execution()).turn !== undefined;
  }

  /**
   * Back to an empty conversation: `main` returns to the empty head, so the next turn is a first turn again
   * (guidance, tools, then its request). A turn that was parked or interrupted ends as reset. Throws while a turn runs.
   */
  async forget(): Promise<void> {
    if (this.activeTurn) throw new Error('A turn is running: stop it first.');
    this.activeTurn = FORGETTING;
    try {
      const { client } = await this.ensure();
      const empty = await emptyHead();
      await client.reset({ key: `forget-${uuid()}`, pointers: { [MAIN]: empty.id } });
      this.fresh = true;
      this.decided.clear();
      this.contextTokens = undefined;
      this.repeats?.reset();
    } finally {
      this.activeTurn = undefined;
    }
  }

  /**
   * Resumes a turn that was parked for approval (with the person's decisions) or interrupted (the tab closed
   * mid-turn: a runner takes it over once the old owner's lease runs out). Returns undefined when there is no such
   * turn.
   */
  async resume(
    hooks: TurnHooks = {},
    decisions: Readonly<Record<string, ApprovalDecision>> = {},
  ): Promise<TurnResult | undefined> {
    if (this.activeTurn) return undefined;
    // Claimed before the first await, as `turn` does, so a turn, a resume or a forget cannot start meanwhile.
    this.activeTurn = RESUMING;
    try {
      for (const [callId, d] of Object.entries(decisions)) this.decided.set(callId, d);
      const harness = await this.ensure();
      const turn = (await harness.client.execution()).turn;
      if (turn === undefined) return undefined;
      this.reset(hooks);
      this.activeTurn = turn.id;
      return await this.resumed(harness, turn);
    } finally {
      this.activeTurn = undefined;
    }
  }

  private async resumed(harness: Harness, turn: NonNullable<Execution['turn']>): Promise<TurnResult> {
    try {
      this.begin(turn.id, true);
      const key = `resume-${turn.id}-${Date.now().toString(36)}`;
      const { launch } =
        turn.park !== undefined ? await harness.client.resume({ key, turnId: turn.id }) : await harness.client.check();
      if (launch !== undefined) harness.supervisor.ensure(harness.id, launch);
      return await this.outcome(harness, { turnId: turn.id });
    } catch (e) {
      return {
        status: 'failed',
        text: this.state.lastText,
        error: e instanceof Error ? e.message : String(e),
        tools: this.state.events,
      };
    }
  }

  /** Waits until the turn (begun by the prompt `key`, or `turnId`) ends or parks. */
  private outcome(harness: Harness, want: { key?: string; turnId?: string }): Promise<TurnResult> {
    let turnId = want.turnId;
    return new Promise<TurnResult>((resolve, reject) => {
      let busy = false;
      let again = false;
      let done = false;
      let unsubscribe: () => void = () => undefined;
      const finish = (f: () => void) => {
        done = true;
        unsubscribe();
        f();
      };
      const check = async () => {
        if (busy) {
          again = true;
          return;
        }
        busy = true;
        try {
          do {
            again = false;
            const execution = await harness.client.execution();
            if (want.key !== undefined && execution.inbox.some((i) => i.key === want.key)) continue;
            if (turnId === undefined) {
              const begun = (await harness.client.turns({ limit: 8 })).find((t) =>
                t.startedBy.includes(want.key ?? ''),
              );
              turnId = begun?.id;
              if (turnId === undefined) continue;
              if (this.activeTurn === want.key) this.activeTurn = turnId;
              this.begin(turnId);
            }
            const result = await this.settled(harness.client, execution, turnId);
            if (result !== undefined) return finish(() => resolve(result));
          } while (again);
        } catch (e) {
          finish(() => reject(e));
        } finally {
          busy = false;
        }
      };
      unsubscribe = harness.client.watch().subscribe(
        () => {
          if (!done) void check();
        },
        (e) => {
          if (!done) finish(() => reject(e));
        },
      );
    });
  }

  /** The turn's result once it ended or parked; undefined while it runs. */
  private async settled(client: HarnessClient, execution: Execution, turnId: string): Promise<TurnResult | undefined> {
    const state = this.state;
    if (execution.turn?.id === turnId) {
      if (execution.turn.park === undefined) return undefined;
      return {
        status: 'awaiting-approval',
        approval: this.awaiting(execution),
        text: state.lastText,
        tools: state.events,
      };
    }
    const summary = (await client.turns({ limit: 16 })).find((t) => t.id === turnId);
    if (summary === undefined || summary.status === 'active') return undefined;
    const status = summary.status === 'completed' ? 'completed' : summary.status === 'failed' ? 'failed' : 'stopped';
    const error = status === 'failed' ? ((await this.endReason(client, summary)) ?? 'The turn failed.') : undefined;
    if (!state.lastText && status === 'completed') {
      const last = (await this.head(client)).findLast((m) => m.role === 'assistant');
      state.lastText = last?.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('') ?? '';
    }
    const u = summary.usage;
    return {
      status,
      text: state.lastText,
      ...(error !== undefined ? { error } : {}),
      tools: state.events,
      usage: {
        input: u.uncachedInput + u.cacheRead + u.cacheWrite,
        output: u.output,
        cached: u.cacheRead,
        written: u.cacheWrite,
        modelCalls: u.modelCalls,
      },
      ...(this.contextTokens !== undefined ? { context: this.contextTokens } : {}),
      ...(state.compactions.length ? { compactions: [...state.compactions] } : {}),
    };
  }

  /** Why a turn ended, in words, from the event that ended it. */
  private async endReason(client: HarnessClient, summary: TurnSummary): Promise<string | undefined> {
    if (summary.lastStep === undefined) return undefined;
    const ending = (await client.steps([summary.lastStep])).findLast((e) => e.end !== undefined);
    const reason = ending?.end?.reason as { reason?: string; message?: string; kind?: string } | undefined;
    if (reason === undefined) return undefined;
    return reason.message ?? (reason.kind !== undefined ? `${reason.reason}: ${reason.kind}` : reason.reason);
  }

  /** The calls a parked turn waits on: its pending calls that need approval and are not decided. */
  private awaiting(execution: Execution): PendingApproval[] {
    const pending = execution.turn?.pending;
    if (pending?.kind !== 'tools') return [];
    const asked = new Set([...this.registrations.values()].filter((r) => r.approval).map((r) => r.definition.tool.id));
    return pending.plan
      .filter((c) => asked.has(c.toolId) && !this.decided.has(c.callId))
      .map((c) => ({ callId: c.callId, tool: c.toolId, args: c.arguments }));
  }

  private modelConfiguration(choice: M | undefined): ModelConfiguration {
    const chosen = choice ?? this.options.models.default;
    const { provider, model } = chosen;
    return {
      target: { provider, model },
      parameters: {
        maxOutputTokens: this.options.maxOutputTokens ?? 16_384,
        ...reasoningFor(chosen, this.options.thinking),
      },
    };
  }

  /**
   * The configuration for this turn's model (the session's own, or the one chosen): its limits resolved, the output it
   * asks for clamped to the model's maximum, its budget.
   */
  private async configurationFor(choice: M | undefined): Promise<TurnConfiguration> {
    const chosen = choice ?? this.options.models.default;
    const adapter = this.adapterFor(chosen.provider);
    if (!adapter) throw new Error(`No adapter for ${chosen.provider}.`);
    this.remember(chosen);
    const model = this.modelConfiguration(choice);
    const capabilities = await this.capabilitiesOf(adapter, model);
    this.outputCap = capabilities.maxOutputTokens;
    const wanted = this.options.maxOutputTokens ?? 16_384;
    const maxOutputTokens =
      capabilities.maxOutputTokens !== undefined ? Math.min(wanted, capabilities.maxOutputTokens) : wanted;
    const configured = { ...model, parameters: { ...model.parameters, maxOutputTokens } };
    return { ...this.toolSet(), model: configured, limits: this.budgetFor(configured).initial };
  }

  /** Cancels the running (or parked) turn. Changes already applied stay (the host can undo the turn). */
  async stop(): Promise<void> {
    if (!this.harness) return;
    const { client } = this.harness;
    const key = `stop-${Date.now().toString(36)}-${uuid()}`;
    try {
      const active = this.activeTurn;
      if (active?.startsWith('prompt-')) await client.stop({ key, item: active });
      else {
        const turnId = (claimed(active) ? undefined : active) ?? (await client.execution()).turn?.id;
        if (turnId !== undefined) await client.stop({ key, turnId });
      }
    } catch {
      // A turn that ended meanwhile needs no stop.
    }
  }
}
