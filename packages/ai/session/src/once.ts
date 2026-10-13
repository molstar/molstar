/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  DisclosureApi,
  InputSpec,
  RunInput,
  ArgsFrom,
  Toolkit,
  ToolRegistration,
  ToolResult,
  ToolSpec,
} from '@molstar/ai-tools';
import { defineTool } from '@molstar/ai-tools';
import { AgentSession, type ModelSource, type SessionThinking } from './session.js';
import type { ModelChoice } from './types.js';

/**
 * One-off runs (tools-and-session.md §7): a throw-away session on a fresh in-memory harness, nothing reaching any
 * conversation, ended by one terminal call: a tool whose call settles the run with its value (the command menu's
 * assistant proposing a filled command). The run is raced against a timeout and a signal.
 */

/** What a run settles with: a terminal tool's call sets it, and the run stops there. */
export class Settlement<V> {
  private result: { value: V } | undefined;
  private onSettle: (() => void) | undefined;

  get settled(): boolean {
    return this.result !== undefined;
  }

  get value(): V | undefined {
    return this.result?.value;
  }

  /** The first call settles; later ones (a response with several terminal calls) change nothing. */
  settle(value: V): void {
    if (this.result !== undefined) return;
    this.result = { value };
    this.onSettle?.();
  }

  /** @internal Starts a run over this settlement: cleared, and `stop` called once it settles. */
  begin(stop: () => void): () => void {
    this.result = undefined;
    this.onSettle = stop;
    return () => {
      this.onSettle = undefined;
    };
  }
}

/** A terminal tool's outcome: the value the run settles with, or a failure that goes back to the model. */
export type TerminalResult<V> = { ok: true; value: V; did?: string } | Extract<ToolResult, { ok: false }>;

export interface TerminalSpec<Api, K extends keyof Api, I extends InputSpec, V>
  extends Omit<ToolSpec<Api, K, I>, 'run' | 'read' | 'idempotent' | 'approval' | 'approvedBy'> {
  /** Checks what the model gave and returns the value; a failure is sent back to the model within the same run. */
  run(input: RunInput<Api, K, ArgsFrom<I>>): TerminalResult<V> | Promise<TerminalResult<V>>;
}

/**
 * A tool for one host's API whose call settles the run: `const terminal = terminalToolFor<Api>()`, then
 * `terminal(settlement, { id, input, modules, run })`. A `run` that returns `{ ok: true, value }` settles it (the
 * session stops, the value is the run's); one that fails answers the model, which may call again.
 */
export function terminalToolFor<Api>() {
  return <V, const K extends keyof Api, const I extends InputSpec>(
    settlement: Settlement<V>,
    spec: TerminalSpec<Api, K, I, V>,
  ): ToolRegistration<Api> =>
    defineTool<Api, K, I>({
      ...spec,
      run: async (input) => {
        const r = await spec.run(input);
        if (!r.ok) return r;
        settlement.settle(r.value);
        return { ok: true, did: r.did ?? 'Done.' };
      },
    });
}

export interface RunOnceOptions<Api extends DisclosureApi, V, M extends ModelChoice = ModelChoice> {
  models: ModelSource<M>;
  /** The model for this run; default: `models.default`. */
  model?: M;
  guidance: string;
  /** The tools, all loaded from the start; the terminal ones among them settle `settlement`. */
  toolkit: Toolkit<Api>;
  api: Api;
  prompt: string;
  /** What the terminal tools settle; the run's `value` is what it holds. */
  settlement?: Settlement<V>;
  signal?: AbortSignal;
  /** Gives up after this long (default 45 s). */
  timeoutMs?: number;
  /** Tool calls the run may make (default 8). */
  maxToolCalls?: number;
  /** Model retries on transient errors (default 1). */
  retries?: number;
  thinking?: SessionThinking;
}

export interface OnceResult<V> {
  /**
   * settled: a terminal call ended it (`value`); completed: the model answered in words (`text`); failed (`error`);
   * timeout and aborted: it was cut short before either.
   */
  status: 'settled' | 'completed' | 'failed' | 'timeout' | 'aborted';
  value?: V;
  /** What the model said last (an answer in words is its reason for not calling a terminal tool). */
  text: string;
  error?: string;
}

/** Runs one prompt on a throw-away session; never throws for a model or tool failure, which are the status. */
export async function runOnce<Api extends DisclosureApi, V = never, M extends ModelChoice = ModelChoice>(
  options: RunOnceOptions<Api, V, M>,
): Promise<OnceResult<V>> {
  const session = new AgentSession<Api, M>({
    models: options.models,
    toolkit: options.toolkit,
    api: options.api,
    guidance: options.guidance,
    startTools: 'all',
    maxToolCalls: options.maxToolCalls ?? 8,
    retries: options.retries ?? 1,
    ...(options.thinking ? { thinking: options.thinking } : {}),
  });
  const settlement = options.settlement;
  const stop = () => void session.stop();
  const end = settlement?.begin(stop);
  let timedOut = false;
  const onAbort = stop;
  options.signal?.addEventListener('abort', onAbort, { once: true });
  const timer = setTimeout(() => {
    timedOut = true;
    stop();
  }, options.timeoutMs ?? 45_000);
  try {
    if (options.signal?.aborted) return { status: 'aborted', text: '' };
    const result = await session.turn({ prompt: options.prompt, ...(options.model ? { model: options.model } : {}) });
    if (settlement?.settled) return { status: 'settled', value: settlement.value as V, text: result.text };
    if (options.signal?.aborted) return { status: 'aborted', text: result.text };
    if (timedOut) return { status: 'timeout', text: result.text };
    if (result.status === 'failed')
      return { status: 'failed', text: result.text, error: result.error ?? 'The run failed.' };
    return { status: 'completed', text: result.text };
  } finally {
    clearTimeout(timer);
    options.signal?.removeEventListener('abort', onAbort);
    end?.();
    session.dispose();
  }
}
