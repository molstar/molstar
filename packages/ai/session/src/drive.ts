/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ModelConfiguration, ReasoningParameters, ToolDefinition } from '@molstar/ai-agent';
import {
  type Decisions,
  type Drive,
  type DriveContext,
  type DriveHandlers,
  defaultDrive,
  drive,
  type HeadEdit,
  type LastCompaction,
  type Limits,
} from '@molstar/ai-harness';
import { type Toolkit, type ToolRegistration, toolSchedule } from '@molstar/ai-tools';
import type { TurnState } from './state.js';
import type { ApprovalDecision, PendingApproval, StatusEvent, ToolEvent } from './types.js';

/** The session's decisions over `defaultDrive` (tools-and-session.md §5). */

/** What the model reads for a call the person did not allow. */
export const DENIED = 'The person did not allow this call.';

/** A turn's budget policy: the limits it starts with, and how far extensions may take them. */
export interface BudgetPolicy {
  readonly initial: Limits;
  readonly ceiling: Limits;
}

/** A turn's budget ceilings by kind of token, as the trace shows them (also read with no session: headless). */
export function ceilingsOf(budget: BudgetPolicy): Partial<Record<'uncachedInput' | 'cacheRead' | 'output', number>> {
  const { uncachedInput, cacheRead, output } = budget.ceiling;
  return {
    ...(uncachedInput !== undefined ? { uncachedInput } : {}),
    ...(cacheRead !== undefined ? { cacheRead } : {}),
    ...(output !== undefined ? { output } : {}),
  };
}

/**
 * A turn's budget, in tokens by kind (the same for every model; credits are shown, never a limit): backstops against
 * a runaway turn, scaled with the compaction trigger (`trigger`, from the model's window), as a turn resends its
 * context on every call: tool calls, uncached input (what a cache miss or a new model costs), cache reads (the
 * context read again, cheap but not free), output (thinking included). The drive extends the short kind toward its
 * ceiling.
 */
export function turnBudget(trigger: number, maxToolCalls = 200): BudgetPolicy {
  return {
    initial: {
      toolCalls: Math.min(maxToolCalls, 60),
      uncachedInput: Math.max(500_000, 4 * trigger),
      cacheRead: Math.max(5_000_000, 40 * trigger),
      output: 200_000,
    },
    ceiling: {
      toolCalls: maxToolCalls,
      uncachedInput: Math.max(2_000_000, 20 * trigger),
      cacheRead: Math.max(25_000_000, 200 * trigger),
      output: 1_000_000,
    },
  };
}

/**
 * An automatic budget extension: the short dimension grows by its initial amount (or to what the attempt needs, if
 * more), capped at the ceiling. Undefined when the ceiling does not allow it.
 */
export function extendBudget(payload: Decisions['budget'][0], policy: BudgetPolicy): Limits | undefined {
  const d = payload.dimension;
  if (d === 'endBy') return undefined;
  const limit = payload.budget.limits[d];
  if (limit === undefined) return undefined;
  const used = (payload.budget.consumed as unknown as Record<string, number | undefined>)[d] ?? 0;
  const need = payload.need ?? 0;
  const ceiling = policy.ceiling[d] ?? Number.POSITIVE_INFINITY;
  if (used + need > ceiling) return undefined;
  return { [d]: Math.min(Math.max(limit + (policy.initial[d] ?? 0), used + need), ceiling) };
}

const LOW_THINKING: ReasoningParameters = { mode: 'adaptive', effort: 'low' };

/** Gemini models that take reasoning parameters (the adapter refuses them for earlier ones). */
const GEMINI_THINKING = /^gemini-(?:[3-9]|[1-9]\d)/;

/**
 * The reasoning parameters for a model: `'model'` leaves it to the model; otherwise the given ones; by default
 * Google's low (low passed every eval, 2–5× faster: more thinking made Gemini read and re-read) and the provider's
 * own elsewhere.
 */
export function thinkingFor(
  provider: string,
  model: string,
  thinking: ReasoningParameters | 'model' | undefined,
): { reasoning?: ReasoningParameters } {
  if (thinking === 'model') return {};
  if (provider === 'google' && !GEMINI_THINKING.test(model)) return {};
  if (thinking) return { reasoning: thinking };
  return provider === 'google' ? { reasoning: LOW_THINKING } : {};
}

/** The providers' names in a failure's words (the host's own table replaces it). */
export const DEFAULT_PROVIDER_NAMES: Readonly<Record<string, string>> = {
  google: 'Google',
  anthropic: 'Anthropic',
  openai: 'OpenAI',
};

/**
 * A model failure in plain words: a context that cannot fit, a request the model cannot take, a missing model or a
 * refused key say what to do; a limit that asks for a wait no turn can make says so.
 */
export function plainFailure(
  failure: { code: string; message: string },
  model?: ModelConfiguration,
  providerNames: Readonly<Record<string, string>> = DEFAULT_PROVIDER_NAMES,
): string {
  if (failure.code === 'capacity')
    return "The conversation no longer fits this model's context window, even condensed. Choose a model with a larger window, or start a new conversation.";
  if (failure.code === 'unrepresentable')
    return `This model cannot take the request as it is (${failure.message}). Try another model.`;
  const who = model === undefined ? 'The provider' : (providerNames[model.target.provider] ?? model.target.provider);
  const wait = waitOf(failure.message);
  if (wait !== undefined && wait > MAX_WAIT_SECONDS)
    return `${who} asks to wait ${waitText(wait)} before the next call (${failure.message})`;
  const status = Number(/\(HTTP (\d{3})\)$/.exec(failure.message)?.[1]);
  if (status === 404) return `No model "${model?.target.model ?? ''}" at ${who}: check the model id in the options.`;
  if (status === 401 || status === 403) return `${who} refused the key (${failure.message}).`;
  return failure.message;
}

/** The longest wait a provider may ask for that the turn makes (seconds). */
const MAX_WAIT_SECONDS = 60;

/** The wait a provider's message asks for, in seconds ("try again in 1h2m3.5s", "in 20s", "in 337ms"). */
function waitOf(message: string): number | undefined {
  const m = /try again in ((?:\d+(?:\.\d+)?(?:ms|[hms]))+)/i.exec(message);
  if (!m) return undefined;
  const per: Record<string, number> = { h: 3600, m: 60, s: 1, ms: 0.001 };
  let seconds = 0;
  for (const [, n, unit] of (m[1] ?? '').matchAll(/(\d+(?:\.\d+)?)(ms|[hms])/gi))
    seconds += Number(n) * (per[(unit ?? 's').toLowerCase()] ?? 1);
  return seconds;
}

/** A wait in words: "63 hours", "5 minutes". */
function waitText(seconds: number): string {
  if (seconds >= 5400) return `${Math.round(seconds / 3600)} hours`;
  if (seconds >= 90) return `${Math.round(seconds / 60)} minutes`;
  return `${Math.round(seconds)} seconds`;
}

// ------------------------------------------------------------------ the drive

/** The policies the drive reads (the session's options). */
export interface DrivePolicy {
  /** Rules for calls that need approval, decided without asking; undefined asks the person (the turn parks). */
  approve?(call: PendingApproval): ApprovalDecision | undefined;
  /** Max model retries on transient errors (rate limits, network): 3. */
  retries?: number;
  /** The back-off step: 1500 ms. */
  retryDelayMs?: number;
  /** Max retries of a tool call that failed transiently: 2. */
  toolRetries?: number;
  /** The back-off step: 1000 ms. */
  toolRetryDelayMs?: number;
  providerNames?: Readonly<Record<string, string>>;
}

/** What the drive needs of the session that owns it. */
export interface DriveHost<Api> {
  readonly toolkit: Toolkit<Api>;
  readonly api: Api;
  readonly policy: DrivePolicy;
  /** The person's answers to calls that needed approval, by call id (resume). */
  readonly decided: ReadonlyMap<string, ApprovalDecision>;
  /** The running turn. */
  turn(): TurnState;
  registration(toolId: string): ToolRegistration<Api> | undefined;
  budgetFor(configuration: ModelConfiguration): BudgetPolicy;
  /** The running turn's model's output maximum (a widened retry stays within it). */
  outputCap(): number | undefined;
  compaction(
    context: DriveContext,
    when: 'force' | 'trigger' | 'soft',
    sized?: { size: { inputTokens: number }; lastCompaction?: LastCompaction },
  ): Promise<HeadEdit | undefined>;
  status(event: StatusEvent): void;
  /** Shows the person a call that ended without running its handler, or that gave up. */
  shown(event: Omit<ToolEvent, 'startedAt'> & { startedAt?: number }): void;
}

/**
 * The session's drive: its policies over `defaultDrive`, with the host's own handlers over them (a handler of the
 * host that answers 'pass' leaves the decision to the session's).
 */
export function createDrive<Api>(host: DriveHost<Api>, overrides: DriveHandlers = {}): Drive {
  const rule = (p: PendingApproval) =>
    host.registration(p.tool)?.approvedBy?.(p.args, host.api) ? ('allow' as const) : undefined;
  const reads = (toolId: string) => host.registration(toolId)?.read === true;

  const permission: DriveHandlers['permission'] = ({ payload, path }) => {
    const p = { callId: payload.call.callId, tool: payload.call.toolId, args: payload.call.arguments };
    const decided = host.decided.get(p.callId) ?? host.policy.approve?.(p) ?? rule(p);
    // Ask: the turn parks (durably) until the person answers and the turn resumes.
    if (decided === undefined) return { park: { note: p as unknown as Json } };
    if (decided === 'deny') {
      // Told to the model and the turn goes on; the person sees a failed call.
      if (path.length === 0)
        host.shown({
          id: p.callId,
          name: nameOf(host, p.tool),
          args: p.args,
          status: 'failed',
          text: DENIED,
          read: reads(p.tool),
        });
      return { deny: DENIED };
    }
    return 'allow';
  };

  const failure: DriveHandlers['failure'] = async ({ payload, path }, context) => {
    const f = payload.failure;
    const attempt = payload.attempt ?? 0;
    const top = path.length === 0;
    const turn = host.turn();
    if (payload.call !== undefined && payload.call.toolId !== '') {
      const { callId, toolId } = payload.call;
      const registration = host.registration(toolId);
      // Most tool errors are deterministic (validation) and go back to the model at once. A transient one (network,
      // 5xx) is retried a few times first, but only when running it again cannot repeat a change: the attempt made
      // none, or the tool is safe to run again (tools-and-session.md §5).
      if (
        f.transient &&
        attempt <= (host.policy.toolRetries ?? 2) &&
        (f.effect === 'none' || registration?.idempotent === true)
      )
        return { retry: { afterMs: (host.policy.toolRetryDelayMs ?? 1000) * attempt } };
      const gaveUp = turn.retrying.get(callId);
      if (gaveUp && top) {
        turn.retrying.delete(callId);
        host.shown(gaveUp);
      }
      // A failed call that is not a read (an edit, a view command, a disclosure) ends the step: the calls after it
      // were written before its result.
      if (reads(toolId)) return { report: {} };
      const why = f.message.split('\n')[0]?.trim() ?? '';
      const skipReason = `Tool execution skipped: ${toolId} failed${why ? ` (${why.length > 160 ? `${why.slice(0, 159)}…` : why})` : ''}, so the calls after it in this response did not run. Fix it, then make the ones still needed again.`;
      if (top) turn.skipShown = `Not run: ${toolId} failed before it.`;
      // Every call after it (the policy for now; the harness can skip only those that depend on it).
      return { report: { skip: { text: skipReason } } };
    }
    if (payload.attempt === undefined) {
      // Projection: a context that no longer fits is compacted when it can be.
      if (f.code === 'capacity') {
        const edit = context.depth === 0 ? await host.compaction(context, 'force') : undefined;
        if (edit !== undefined) return { compact: edit };
      }
      return {
        fail: { reason: f.code, message: plainFailure(f, context.configuration?.model, host.policy.providerNames) },
      };
    }
    const wait = waitOf(f.message);
    const retryable = f.transient && (wait === undefined || wait <= MAX_WAIT_SECONDS);
    if (retryable && attempt <= (host.policy.retries ?? 3)) {
      host.status({ kind: 'retrying', detail: f.message });
      return { retry: { afterMs: (host.policy.retryDelayMs ?? 1500) * attempt } };
    }
    return {
      fail: { reason: f.code, message: plainFailure(f, context.configuration?.model, host.policy.providerNames) },
    };
  };

  const modelResult: DriveHandlers['model-result'] = ({ payload }, context) => {
    // Cut short at the output limit (e.g. long thinking): once, ask again with twice the room.
    if (context.depth > 0) return 'accept';
    const model = context.configuration?.model;
    const limit = model?.parameters.maxOutputTokens;
    const cap = host.outputCap();
    const room = limit ? Math.min(limit * 2, cap ?? limit * 2) : undefined;
    const turn = host.turn();
    if (payload.finish === 'output-limit' && !turn.widened && model && limit && room && room > limit) {
      turn.widened = true;
      host.status({ kind: 'widened' });
      return { replace: { ...model, parameters: { ...model.parameters, maxOutputTokens: room } } };
    }
    return 'accept';
  };

  const budget: DriveHandlers['budget'] = ({ payload }, context) => {
    // A nested run (a tool's helper) keeps the limits its tool gave it.
    if (context.depth > 0) return { fail: { reason: 'budget', message: 'The helper reached its budget.' } };
    const policy = context.configuration === undefined ? undefined : host.budgetFor(context.configuration.model);
    const extension = policy === undefined ? undefined : extendBudget(payload, policy);
    if (extension === undefined) return { fail: { reason: 'budget', message: 'The turn reached its budget.' } };
    const said = Object.entries(extension)
      .map(([d, v]) => `${d} → ${v}`)
      .join(', ');
    host.status({ kind: 'budget-extended', detail: said });
    return { extend: extension };
  };

  const repair: DriveHandlers['repair'] = ({ payload, path }) => {
    const { call, reason } = payload;
    const top = path.length === 0;
    const { names } = host.toolkit;
    const registration = host.registration(call.toolId);
    const name = registration?.definition.name ?? call.toolId;
    const unknown = `There is no tool "${call.toolId}". Load tools with ${names.skills} or ${names.tools}; their descriptions list what there is.`;
    let text: string;
    let load: ToolDefinition | undefined;
    if (reason === 'invalid-arguments') text = `${call.toolId} takes its arguments as one JSON object.`;
    else if (reason === 'unknown-tool' || (reason === 'not-loaded' && registration === undefined)) text = unknown;
    else if (reason === 'not-loaded' && registration !== undefined) {
      // A call to a tool the conversation has not loaded: answered failed and its definition loaded, so the model
      // calls it again.
      const place = host.toolkit.placementOf(call.toolId);
      load = registration.definition;
      const where =
        place.kind === 'group' && !place.extra
          ? ` It belongs to the ${place.group.id} skill: ${names.skills}(["${place.group.id}"]) gives its guidance.`
          : place.kind === 'entry'
            ? ` It comes with ${place.entry.id}: ${names.skills}(["${place.entry.id}"]) gives its text.`
            : '';
      text = `${name} was not loaded; it is now. Call it again.${where}`;
    } else {
      // Not enabled, or no handler: nothing to load, and the rest of the response is not affected.
      text = `${name} cannot run now.`;
    }
    const skipText =
      load !== undefined
        ? `Not run: ${load.name} was not loaded yet, so the rest of this response did not run; call again what is still needed.`
        : undefined;
    if (top) {
      host.shown({ id: call.callId, name, args: call.arguments, status: 'failed', text, read: reads(call.toolId) });
      if (skipText !== undefined) host.turn().skipShown = skipText;
    }
    return {
      report: {
        parts: [{ kind: 'text', text }],
        ...(load !== undefined
          ? { contributions: host.toolkit.definitionEntries([load]), skip: { text: skipText ?? 'Not run.' } }
          : {}),
      },
    };
  };

  return drive(
    overrides,
    drive(
      {
        permission,
        schedule: ({ payload }) => ({
          // Edits one at a time in the model's order; reads between them together.
          stages: toolSchedule(
            payload.calls.map((c) => ({ callId: c.callId, tool: { id: c.toolId } })),
            reads,
          ),
        }),
        failure,
        'model-result': modelResult,
        budget,
        repair,
        // Mid-turn compaction once the context reaches the trigger. Queued prompts are not steered in: each prompt is
        // its own turn.
        boundary: async ({ payload }, context) => {
          const edit = await host.compaction(context, 'trigger', payload);
          return edit === undefined ? 'proceed' : { compact: edit };
        },
        // One prompt per turn; from the soft threshold, the turn begins with a compaction.
        'next-turn': async ({ payload }, context) => {
          const item = context.inbox[0];
          if (item === undefined) return 'release';
          const compact = await host.compaction(context, 'soft', payload);
          return { begin: { keys: [item.key], ...(compact !== undefined ? { compact } : {}) } };
        },
      },
      defaultDrive(),
    ),
  );
}

const nameOf = <Api>(host: DriveHost<Api>, toolId: string) => host.registration(toolId)?.definition.name ?? toolId;
