/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '../messages/index.js';
import { type ArgsFrom, type InputSpec, inputOf } from './input.js';
import type { Nest } from './nest.js';
import { type CallIdentity, type Projection, registered, type ToolRegistration } from './registry.js';
import { failed, render, type ToolResult } from './result.js';
import { normalize } from './shape.js';

/** What a tool's `run` gets (tools-and-session.md §1). */
export interface RunInput<Api, K extends keyof Api, A, M = unknown> {
  /** Its arguments, normalized and checked, typed as the input says. */
  args: A;
  /** The modules it declared, frozen. */
  api: Projection<Api, K>;
  signal: AbortSignal;
  /** What a long call is doing, for the person watching (never the model). */
  progress: (text: string) => void;
  /** Runs a helper under this call, when a turn made it (NestRequest). */
  nest?: Nest<Api, M>;
  /** Which call this is: its key and whether an earlier attempt may have run, for a tool that must answer a re-run. */
  call: CallIdentity;
}

export interface ToolSpec<Api, K extends keyof Api, I extends InputSpec, M = unknown> {
  id: string;
  description: string;
  /** A field shape (tools-and-session.md §2). */
  input: I;
  modules: readonly K[];
  /** Of its modules, those it does without where the host lacks them. */
  optional?: readonly K[];
  /** Only reads: its result depends on its arguments and the host's state alone. */
  read?: boolean;
  /** Safe to run again after an attempt that may have changed something (the drive retries only these). */
  idempotent?: boolean;
  /** Runs only with the person's approval. */
  approval?: boolean;
  /** Approves a call without asking. */
  approvedBy?: (args: Json, api: Api) => boolean;
  run(input: RunInput<Api, K, ArgsFrom<I>, M>): ToolResult | Promise<ToolResult>;
}

/**
 * A tool: its input, the API modules it uses (it gets only those), and what it does. Its arguments arrive normalized
 * and checked; arguments it does not take are left out, with a warning naming the ones it takes. `Api` is the host's
 * object of modules; to give it once for many tools, use `toolFor<Api>()`.
 */
export function defineTool<Api, const K extends keyof Api, const I extends InputSpec, M = unknown>(
  spec: ToolSpec<Api, K, I, M>,
): ToolRegistration<Api> {
  const input = inputOf(spec.input);
  const fields = input.shape.fields;
  const takes = Object.keys(fields).join(', ') || 'no arguments';
  const split = (args: Readonly<Record<string, unknown>>) => ({
    known: Object.fromEntries(Object.entries(args).filter(([k]) => Object.hasOwn(fields, k))),
    unknown: Object.keys(args).filter((k) => !Object.hasOwn(fields, k)),
  });
  const registration = registered<Api, K>(
    { tool: { id: spec.id }, name: spec.id, description: spec.description, inputSchema: input.jsonSchema },
    spec.modules,
    async ({ args, api, signal, progress, nest, call }) => {
      // Arguments it does not take are left out, and the result says so (the call still does what it can: a small
      // model's stray field should not cost it the call).
      const { known, unknown } = split(args);
      const parsed = input.parse(normalize(input.shape, known));
      if (!parsed.ok) return { ok: false, code: 'tool.arguments', message: parsed.message };
      const r = await spec.run({
        args: parsed.value as ArgsFrom<I>,
        api,
        signal,
        progress,
        call,
        ...(nest ? { nest: nest as Nest<Api, M> } : {}),
      });
      if (!unknown.length || !r.ok) return render(r);
      const ignored = `Ignored ${unknown.join(', ')}: ${spec.id} takes ${takes}.`;
      return render({ ...r, warnings: [...(r.warnings ?? []), ignored] });
    },
    {
      ...(spec.read ? { read: true } : {}),
      ...(spec.idempotent ? { idempotent: true } : {}),
      ...(spec.optional ? { optional: spec.optional } : {}),
      ...(spec.approval ? { approval: true } : {}),
      ...(spec.approvedBy ? { approvedBy: spec.approvedBy } : {}),
    },
  );
  DIRECT.set(registration, async (args, api, signal, call = DIRECT_CALL) => {
    const { known, unknown } = split(args);
    if (unknown.length) return failed('tool.arguments', `${spec.id} takes ${takes}, not ${unknown.join(', ')}.`);
    const parsed = input.parse(normalize(input.shape, known));
    if (!parsed.ok) return failed('tool.arguments', parsed.message);
    const projected = Object.freeze(Object.fromEntries(spec.modules.map((k) => [k, api[k]]))) as unknown as Projection<
      Api,
      K
    >;
    return spec.run({ args: parsed.value as ArgsFrom<I>, api: projected, signal, progress: () => {}, call });
  });
  return registration;
}

/** `defineTool` for one host: `const tool = toolFor<ViewerApi>()`, then `tool({ id, input, modules, run })`. */
export function toolFor<Api, M = unknown>() {
  return <const K extends keyof Api, const I extends InputSpec>(spec: ToolSpec<Api, K, I, M>): ToolRegistration<Api> =>
    defineTool<Api, K, I, M>(spec);
}

const DIRECT_CALL: CallIdentity = { callId: 'direct', key: 'direct', prior: 'none' };

export type DirectRun<Api> = (
  args: Readonly<Record<string, unknown>>,
  api: Api,
  signal: AbortSignal,
  call?: CallIdentity,
) => Promise<ToolResult>;
const DIRECT = new WeakMap<object, DirectRun<never>>();

/**
 * A tool called from code (a skill that writes code against the tools): its arguments checked as a call's are, its
 * result as data, not rendered. Only tools made by `defineTool` can be called so.
 */
export function directRun<Api>(registration: ToolRegistration<Api>): DirectRun<Api> | undefined {
  return DIRECT.get(registration) as DirectRun<Api> | undefined;
}
