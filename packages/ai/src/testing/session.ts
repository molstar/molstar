/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message, ProviderBindings } from '../messages/index.js';
import {
  createHarness,
  type Drive,
  defaultDrive,
  type Harness,
  type HarnessEvent,
  type HarnessOptions,
  type Observation,
  type ToolBindings,
  type TurnConfiguration,
  type TurnHandle,
  type TurnOutcome,
  type TurnSummary,
} from '../harness/index.js';
import { SimulatedClock } from './clock.js';
import { simulatedIds } from './ids.js';
import { declareTools, user } from './messages.js';
import { type ModelStep, type ScriptedModel, type ScriptedModelOptions, scriptedModel } from './model.js';
import { type ScriptedTools, scriptedTools, type ToolRun, type ToolSpec } from './tools.js';

export interface SessionSetup {
  readonly model?: readonly ModelStep[];
  readonly modelOptions?: Omit<ScriptedModelOptions, 'clock'>;
  readonly tools?: Readonly<Record<string, ToolSpec | ToolRun>>;
  /** Tools the first prompt declares and loads (default: every scripted tool). */
  readonly declare?: readonly string[];
  /** Tools the first prompt declares without loading. */
  readonly deferred?: readonly string[];
  /** Merged over the default: the scripted model, every tool enabled and preauthorized, no limits. */
  readonly configuration?: Partial<TurnConfiguration>;
  readonly drive?: Drive;
  readonly limits?: HarnessOptions['limits'];
  /** The model bindings, from the scripted model's (default: those). */
  models?(model: ScriptedModel): ProviderBindings;
  /** The tool bindings (default: the scripted tools'; their definitions are declared either way). */
  readonly toolBindings?: ToolBindings;
}

export interface SendOptions {
  /** This turn's configuration (default: the session's). */
  readonly configuration?: TurnConfiguration;
  /** Run the turn until it ends or parks, in virtual time (default true); else only begin it. */
  readonly wait?: boolean;
}

export interface Session {
  readonly harnessId: string;
  readonly clock: SimulatedClock;
  readonly model: ScriptedModel;
  readonly tools: ScriptedTools;
  readonly harness: Harness;
  /** The turns' configuration (`setup.configuration` over the default). */
  readonly configuration: TurnConfiguration;
  /** Every observation, in order, since the session was made. */
  readonly observations: Observation[];
  /** Begins a turn with a prompt (the first carries the tool declarations); runs it unless `wait` is false. */
  send(
    prompt: string | readonly Message[],
    options?: SendOptions,
  ): Promise<{ handle: TurnHandle; outcome?: TurnOutcome }>;
  /** Waits for a turn's outcome in virtual time. */
  run(handle: TurnHandle): Promise<TurnOutcome>;
  /** Resumes the parked turn and runs it to its outcome. */
  resume(): Promise<TurnOutcome>;
  events(): Promise<HarnessEvent[]>;
  /** The messages of `main`'s current head. */
  main(): Promise<Message[]>;
  turns(): Promise<TurnSummary[]>;
}

let sessions = 0;

/** A harness over a scripted model and scripted tools, on a simulated clock (testing.md). */
export function createSession(setup: SessionSetup = {}): Session {
  const clock = new SimulatedClock();
  const ids = simulatedIds();
  const harnessId = `harness-${++sessions}`;
  const model = scriptedModel(setup.model ?? [], { ...setup.modelOptions, clock });
  const tools = scriptedTools(setup.tools ?? {});
  const toolIds = tools.definitions.map((definition) => definition.tool.id);
  const configuration: TurnConfiguration = {
    model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
    tools: { enabled: toolIds, preauthorized: toolIds },
    limits: {},
    ...setup.configuration,
  };
  const declared = (setup.declare ?? toolIds).map((id) => tools.definition(id));
  const deferred = (setup.deferred ?? []).map((id) => tools.definition(id));
  const harness = createHarness({
    id: harnessId,
    models: setup.models?.(model) ?? model.bindings,
    tools: setup.toolBindings ?? tools.bindings,
    drive: setup.drive ?? defaultDrive(),
    configure: () => configuration,
    clock,
    ids,
    ...(setup.limits !== undefined ? { limits: setup.limits } : {}),
  });
  const observations: Observation[] = [];
  harness.observe().subscribe((observation) => observations.push(observation));
  let first = true;

  const run = (handle: TurnHandle) => clock.runUntil(handle.done);

  return {
    harnessId,
    clock,
    model,
    tools,
    harness,
    configuration,
    observations,
    run,
    async send(prompt, options = {}) {
      const messages: Message[] = typeof prompt === 'string' ? [user(prompt)] : [...prompt];
      if (first) {
        first = false;
        if (deferred.length > 0) messages.unshift(declareTools(deferred, true));
        if (declared.length > 0) messages.unshift(declareTools(declared));
      }
      const handle = harness.turn({
        messages,
        ...(options.configuration !== undefined ? { configuration: options.configuration } : {}),
      });
      if (options.wait === false) return { handle };
      return { handle, outcome: await run(handle) };
    },
    resume: () => run(harness.resume()),
    async events() {
      const out: HarnessEvent[] = [];
      for (;;) {
        const page = await harness.events({ after: out.at(-1)?.seq ?? 0, limit: 500 });
        out.push(...page);
        if (page.length < 500) return out;
      }
    },
    async main() {
      const execution = await harness.execution();
      const [head] = await harness.heads([execution.main]);
      return head === undefined ? [] : await harness.entries(head.entries);
    },
    turns: () => harness.turns({ limit: 100 }),
  };
}
