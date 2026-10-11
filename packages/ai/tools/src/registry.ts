/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ToolDefinition } from '@molstar/ai-agent';
import type { ToolBindings, ToolFrame, ToolHandler } from '@molstar/ai-harness';
import type { Nest, NestedBinding, NestedRuns } from './nest.js';
import type { RepeatedReads } from './repeats.js';
import { type Effect, type ToolOutcome, toFrame } from './result.js';

/** Which call this is, for a tool that must answer a re-run of it (tools-and-session.md §1). */
export interface CallIdentity {
  readonly callId: string;
  /** The call's key: the same for a re-run of the same call. */
  readonly key: string;
  /** `possible`: an earlier attempt may have run (the harness is re-running it). */
  readonly prior: 'none' | 'possible';
}

/** The modules of the host's API a tool declared, and only those. */
export type Projection<Api, K extends keyof Api> = Readonly<Pick<Api, K>>;

export interface ToolRegistration<Api> {
  readonly definition: ToolDefinition;
  readonly modules: readonly (keyof Api)[];
  /** Of its modules, those it uses when the host has them and does without otherwise. */
  readonly optional?: readonly (keyof Api)[];
  /** Only reads: its result depends on its arguments and the host's state alone (see RepeatedReads). */
  readonly read?: boolean;
  /** Safe to run again after an attempt that may have changed something: the drive may retry it (harness.md §7). */
  readonly idempotent?: boolean;
  /**
   * Runs only with the person's approval of each call (or a rule they set): not preauthorized, so the harness stops at
   * its permission step.
   */
  readonly approval?: boolean;
  /** A rule that approves a call without asking. */
  readonly approvedBy?: (args: Json, api: Api) => boolean;
  /** Resolves the handler against the full API, projected to the declared modules. */
  handler(api: Api, options?: HandlerOptions<Api>): ToolHandler;
}

/** What binds a handler beyond the tool itself. */
export interface HandlerOptions<Api> {
  readonly repeats?: RepeatedReads;
  readonly nested?: NestedRuns<Api>;
  readonly afterWrite?: AfterWrite;
}

/** Runs after a call that changed something, before its result is recorded (tools-and-session.md §1). */
export type AfterWrite = (call: { toolId: string; callId: string }) => Promise<void>;

export interface RegisterOptions<Api> {
  read?: boolean;
  idempotent?: boolean;
  optional?: readonly (keyof Api)[];
  approval?: boolean;
  approvedBy?: (args: Json, api: Api) => boolean;
}

/** What a registered tool's `execute` gets. */
export interface ExecuteInput<Api, K extends keyof Api> {
  args: Record<string, Json>;
  api: Projection<Api, K>;
  signal: AbortSignal;
  /** What a long call is doing ("RCSB FATCAT-flexible: running, 6 s"), for the person watching. */
  progress: (text: string) => void;
  /** Runs a helper under this call, when a turn made it (absent when called directly). */
  nest?: Nest<Api>;
  call: CallIdentity;
}

const isRecord = (v: Json): v is Record<string, Json> => typeof v === 'object' && v !== null && !Array.isArray(v);

/**
 * Registers a tool: its definition, the API modules it may use, and how it executes. The execute function only sees the
 * modules it declared, frozen. The handler pumps the call: progress frames while it runs, then one result frame whose
 * `effect` is the outcome's own or derived (tools-and-session.md §1):
 *
 * | the call                              | effect                       |
 * | ------------------------------------- | ---------------------------- |
 * | success of a read                     | none                         |
 * | success of a write                    | applied                      |
 * | the tool's own `{ ok: false }`        | none                         |
 * | a throw                               | none for a read, else unknown|
 * | the `afterWrite` hook failing         | applied, transient           |
 */
export function registered<Api, const K extends keyof Api>(
  definition: ToolDefinition,
  modules: readonly K[],
  execute: (input: ExecuteInput<Api, K>) => ToolOutcome | Promise<ToolOutcome>,
  options: RegisterOptions<Api> = {},
): ToolRegistration<Api> {
  const read = options.read === true;
  return Object.freeze({
    definition,
    modules,
    ...(options.optional?.length ? { optional: options.optional } : {}),
    ...(read ? { read: true } : {}),
    ...(options.idempotent ? { idempotent: true } : {}),
    ...(options.approval ? { approval: true } : {}),
    ...(options.approvedBy ? { approvedBy: options.approvedBy } : {}),
    handler(api: Api, handlerOptions: HandlerOptions<Api> = {}): ToolHandler {
      const { repeats, nested, afterWrite } = handlerOptions;
      const projected = Object.freeze(Object.fromEntries(modules.map((k) => [k, api[k]]))) as Projection<Api, K>;
      return {
        async *invoke(input) {
          const args = isRecord(input.arguments) ? input.arguments : {};
          const nest: Nest<Api> | undefined = nested ? (request) => nested.run(input, request) : undefined;
          const call: CallIdentity = { callId: input.callId, key: input.key, prior: input.prior };
          const repeat = read ? repeats?.check(definition, args) : undefined;
          if (!read) repeats?.acted(definition);
          // Progress frames go out while the call runs (the harness passes them to observers as preview deltas, never
          // into the journal); the result follows them.
          const pending: ToolFrame[] = [];
          let wake: (() => void) | undefined;
          let done = false;
          const progress = (text: string) => {
            if (done) return;
            pending.push({ kind: 'progress', value: text });
            wake?.();
          };
          let outcome: ToolOutcome | undefined;
          const running = (async () => {
            try {
              outcome =
                repeat?.outcome ??
                (await execute({
                  args,
                  api: projected,
                  signal: input.signal,
                  progress,
                  call,
                  ...(nest ? { nest } : {}),
                }));
            } finally {
              done = true;
              wake?.();
            }
          })();
          let failure: unknown;
          running.catch((e: unknown) => {
            failure = e ?? new Error('The tool failed.');
          });
          while (!done || pending.length) {
            const next = pending.shift();
            if (next) {
              yield next;
              continue;
            }
            await new Promise<void>((resolve) => {
              wake = resolve;
              if (done || pending.length) resolve();
            });
          }
          await running.catch(() => undefined);
          let result: ToolOutcome;
          let effect: Effect;
          if (failure === undefined && outcome) {
            result = outcome;
            effect = outcome.effect ?? (outcome.ok && !read ? 'applied' : 'none');
            // A call that changed something: kept before the harness records it as done, or failed (and retried when
            // the effect allows) when it cannot be kept.
            if (outcome.ok && !read && !repeat?.outcome && effect !== 'none' && afterWrite) {
              try {
                await afterWrite({ toolId: definition.tool.id, callId: input.callId });
              } catch (e) {
                result = {
                  ok: false,
                  code: 'tool.transient',
                  message: `The change could not be saved: ${e instanceof Error ? e.message : String(e)}`,
                  transient: true,
                };
                effect = 'applied';
              }
            }
            if (result.ok && !repeat?.outcome) repeat?.remember();
            const nudge = result.ok && repeat ? repeats?.nudge() : undefined;
            if (result.ok && nudge) result = { ...result, text: `${result.text}\n\n${nudge}` };
          } else {
            const e = failure;
            const transient = (e as { transient?: unknown } | null)?.transient === true;
            result = {
              ok: false,
              code: transient ? 'tool.transient' : 'tool.exception',
              message: e instanceof Error ? e.message : String(e),
              ...(transient ? { transient } : {}),
            };
            effect = read ? 'none' : 'unknown';
          }
          yield toFrame(result, effect);
        },
      };
    },
  });
}

export interface BindingOptions<Api> {
  readonly repeats?: RepeatedReads;
  readonly nested?: NestedBinding<Api>;
  readonly afterWrite?: AfterWrite;
}

export function createToolBindings<Api>(
  registrations: readonly ToolRegistration<Api>[],
  api: Api,
  options: BindingOptions<Api> = {},
): ToolBindings {
  const byId = new Map(registrations.map((r) => [r.definition.tool.id, r]));
  return {
    resolve: (definition) => byId.get(definition.tool.id)?.handler(api, options) ?? options.nested?.resolve(definition),
  };
}
