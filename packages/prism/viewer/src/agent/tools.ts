/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import {
  type DisclosureApi,
  failed,
  type FieldShape,
  type ToolRegistration,
  type ToolResult,
  toolFor,
} from '@molstar/ai-tools';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ApiResult, Changed, OperationDefinition, OperationDescription, ViewerApi } from '../api/area.js';
import { paramsShape } from './params.js';
import { type AgentPolicy, mayAsk } from './policies.js';

// The agent's tools (agent.md §2): thin adapters over operations. A tool takes flat fields in the model's words, maps
// them to an operation's values and runs it with the agent's origin (one run path: the history, the checkpoint and the
// run log are a person's too), and renders the result in one order: what it did, what changed, warnings, what it read,
// the next step.

/** What the agent's tools get of the viewer: the API bound to the turn that runs (`agent`, with its turn). */
export interface AgentViewer {
  /** The API as the current turn's agent calls it. */
  readonly api: ViewerApi;
  /** The agent's previous finished turn (`undo_turn`'s). */
  readonly lastTurn: string | undefined;
  /** The operation's description: its params as they are now (registries change), and what this viewer supports. */
  describe(id: string): OperationDescription | undefined;
}

/** The host's modules for the tools (`ai-tools`: each tool gets only what it declares). */
export interface ViewerAgentApi extends DisclosureApi {
  readonly viewer: AgentViewer;
  /** The scene in words, for the reads (`account.ts`). */
  readonly account: AgentAccountModule;
}

export interface AgentAccountModule {
  /** The scene's outline: structures, components, representations, volumes, measurements; at most `max` lines. */
  outline(max?: number, detail?: 'summary' | 'components' | 'full'): string;
  /** What is selected and focused, in words and as addresses. */
  focus(): string;
}

export const tool = toolFor<ViewerAgentApi>();

/** A tool over an operation (agent.md §2). */
export interface OperationToolSpec {
  /** The tool's id, in the model's words (`load_entry`, `color`). */
  readonly id: string;
  /** The operation it runs. */
  readonly operation: string;
  /** What the model reads: when to use it, what it does; the operation's title when not given. */
  readonly description?: string;
  /** The operation's params offered, by name (default: all that JSON can carry). */
  readonly fields?: readonly string[];
  /** Values set for the model (`mode: 'add'`). */
  readonly fixed?: Readonly<Record<string, unknown>>;
  /** Field formats said in their descriptions (`target`: a target string, targets.md). */
  readonly formats?: Readonly<Record<string, string>>;
}

/**
 * A tool over an operation. Its input is the operation's params as field shapes (`paramsShape`), described now: a
 * registry's choices are what they are when it is made; the toolkit makes it again between turns when they change
 * (`ViewerToolkit.refresh`, agent.md §2).
 */
export function operationTool(
  viewer: AgentViewer,
  spec: OperationToolSpec,
  policy?: AgentPolicy,
): ToolRegistration<ViewerAgentApi> {
  const described = viewer.describe(spec.operation);
  if (!described) throw new Error(`No operation '${spec.operation}' for the tool '${spec.id}'.`);
  const { definition } = described;
  const offered = Object.fromEntries(
    Object.entries(described.params).filter(
      ([name, p]) => !p.isHidden && !(name in (spec.fixed ?? {})) && (!spec.fields || spec.fields.includes(name)),
    ),
  );
  // Targets (a `PD.Value` with a string default: `TargetParam`) are strings to the model (agent.md §2).
  const asText = Object.fromEntries(
    Object.entries(offered).map(([name, p]) =>
      p.type === 'value' && typeof p.defaultValue === 'string'
        ? [name, PD.Text(p.defaultValue, { label: p.label, description: p.description })]
        : [name, p],
    ),
  );
  const shape = paramsShape(asText, { formats: spec.formats ?? targetFormats(offered) });
  // Every call that may ask goes through the policy, decided at call time on the values the operation will run with
  // (agent.md §4: destructive, upload and what no undo restores ask; network by address; a change while the turn
  // cannot be undone asks once). A read or a view change that reaches no address never asks.
  const asks = mayAsk(definition);
  const valuesOf = (args: Json) => ({ ...shape.toValues((args ?? {}) as Record<string, unknown>), ...spec.fixed });
  return tool({
    id: spec.id,
    // What this viewer cannot do at all is said; what cannot be done now (nothing loaded) is the call's answer.
    description: describe(spec, definition, described.supported),
    input: shape.shape,
    modules: ['viewer'],
    read: !!definition.read,
    // A read or an already-so answer is safe to repeat; a change is not.
    idempotent: !!definition.read,
    approval: asks,
    ...(asks
      ? {
          approvedBy: (args: Json) => {
            // Values the model's arguments do not make (a bad field) are the person's to see: asked.
            try {
              return policy ? policy.approves(definition, valuesOf(args)) : false;
            } catch {
              return false;
            }
          },
        }
      : {}),
    run: async ({ args, api, signal }) => {
      const wrong = shape.check((args ?? {}) as Record<string, unknown>);
      if (wrong) return failed('invalid', wrong);
      const values = valuesOf(args);
      const result = await api.viewer.api.run(spec.operation, values, { signal });
      return toolResult(result);
    },
  });
}

/** The operation's words for the model: the tool's description, and why it cannot run now when it cannot. */
function describe(spec: OperationToolSpec, def: OperationDefinition<any, any>, supported: true | string): string {
  const text = spec.description ?? `${def.title}.`;
  const view = def.history === false && !def.read ? ' It changes the view only: the structures stay as they are.' : '';
  return supported === true ? `${text}${view}` : `${text}${view} (Not in this viewer: ${supported})`;
}

/** Text fields that take targets (`TargetParam`: a PD.Value with a string default) say so. */
function targetFormats(params: PD.Params): Record<string, string> {
  const out: Record<string, string> = {};
  for (const [name, p] of Object.entries(params)) {
    if (p.type === 'value' && typeof p.defaultValue === 'string') out[name] = 'target';
  }
  return out;
}

/** An operation's result as the model reads it (ai-tools `render`: did, changed, warnings, body, next). */
export function toolResult(r: ApiResult<unknown>): ToolResult {
  if (!r.ok) return failed(r.code, r.message, r.hint);
  const changed = changedWords(r.changed);
  return {
    ok: true,
    did: r.unchanged ? `${r.did.replace(/\.$/, '')}: already so.` : r.did,
    ...(changed.length ? { changed } : {}),
    ...(r.warnings?.length ? { warnings: r.warnings } : {}),
    ...(r.next ? { next: r.next } : {}),
    // Effects: nothing changed answers as a read does (the drive may repeat it).
    ...(r.unchanged ? { effect: 'none' as const } : {}),
  };
}

function changedWords(changed: Changed | undefined): string[] {
  if (!changed) return [];
  const out: string[] = [];
  if (changed.nodes?.length) out.push(`${changed.nodes.length} scene item${changed.nodes.length === 1 ? '' : 's'}`);
  if (changed.selection) out.push('the selection');
  if (changed.view) out.push('the view');
  if (changed.sequence) out.push('the sequence');
  return out;
}

/** A field shape's type, for a hand-written tool (re-exported for the viewer's own tools). */
export type { FieldShape };
