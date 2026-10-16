/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message, ModelBindings } from '@molstar/ai-agent';
import {
  createHarnessHost,
  createRunner,
  type Drive,
  defaultDrive,
  type HarnessClient,
  type HarnessEvent,
  type HarnessHost,
  type HarnessStorage,
  type HostOptions,
  type LaunchKey,
  type Observation,
  type Runner,
  type RunnerExit,
  type RunnerOptions,
  type RunnerPort,
  type ToolBindings,
  type TurnConfiguration,
  type TurnSummary,
} from '@molstar/ai-harness';
import { SimulatedClock } from './testing/clock.js';
import { simulatedIds } from './testing/ids.js';
import { declareTools, user } from './testing/messages.js';
import { type ModelStep, type ScriptedModel, type ScriptedModelOptions, scriptedModel } from './testing/model.js';
import { type ScriptedTools, scriptedTools, type ToolRun, type ToolSpec } from './testing/tools.js';

/** What a backend gives the suite. */
export interface ConformanceFixture {
  readonly name: string;
  /** A fresh, empty storage for one case. */
  storage(): HarnessStorage | Promise<HarnessStorage>;
}

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
  readonly runnerLimits?: RunnerOptions['limits'];
  readonly host?: Omit<HostOptions, 'clock'>;
}

export interface SendOptions {
  readonly key?: string;
  readonly options?: Json;
  readonly interrupt?: { readonly turnId: string };
  /** Run a runner if the submit launched one (default true). */
  readonly run?: boolean;
  readonly runner?: Runner;
}

export interface Session {
  readonly harnessId: string;
  readonly clock: SimulatedClock;
  readonly model: ScriptedModel;
  readonly tools: ScriptedTools;
  readonly host: HarnessHost;
  readonly client: HarnessClient;
  readonly port: RunnerPort;
  /** A runner over this session's port; its observations collect in `observations`. */
  runner(overrides?: Partial<RunnerOptions>): Runner & { readonly observations: Observation[] };
  /** Submits a prompt; if that launched a runner, runs one until it exits. */
  send(prompt: string | readonly Message[], options?: SendOptions): Promise<{ launch?: LaunchKey; exit?: RunnerExit }>;
  /** Runs a runner for a launch key until it exits, in virtual time. */
  run(launch: LaunchKey, runner?: Runner): Promise<RunnerExit>;
  events(): Promise<HarnessEvent[]>;
  /** The messages of `main`'s current head. */
  main(): Promise<Message[]>;
  turns(): Promise<TurnSummary[]>;
}

let sessions = 0;

export async function createSession(fixture: ConformanceFixture, setup: SessionSetup = {}): Promise<Session> {
  const clock = new SimulatedClock();
  const ids = simulatedIds();
  const harnessId = `harness-${++sessions}`;
  // Storage work in flight is work in progress for the simulated clock (slow backends).
  const inner = await fixture.storage();
  const storage: HarnessStorage = {
    transact: (harnessId, run) => clock.track(() => inner.transact(harnessId, run)),
    putBlobs: (harnessId, blobs) => clock.track(() => inner.putBlobs(harnessId, blobs)),
    getBlobs: (harnessId, refs) => clock.track(() => inner.getBlobs(harnessId, refs)),
    watch: (harnessId) => inner.watch(harnessId),
    drop: (harnessId) => clock.track(() => inner.drop(harnessId)),
  };
  const host = createHarnessHost(storage, { ...setup.host, clock });
  await host.create(harnessId, { key: `create-${harnessId}` });
  const client = host.client(harnessId);
  const port = host.runnerPort(harnessId);
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
  let first = true;
  let runners = 0;

  const runner: Session['runner'] = (overrides = {}) => {
    const observations: Observation[] = [];
    const made = createRunner({
      port,
      runnerId: `runner-${++runners}`,
      models: model.bindings satisfies ModelBindings,
      tools: tools.bindings satisfies ToolBindings,
      drive: setup.drive ?? defaultDrive(),
      configure: () => configuration,
      clock,
      ids,
      ...(setup.runnerLimits !== undefined ? { limits: setup.runnerLimits } : {}),
      ...overrides,
    });
    made.observe().subscribe((observation) => observations.push(observation));
    return Object.assign(made, { observations });
  };

  const run: Session['run'] = (launch, given) => clock.runUntil((given ?? runner()).run(launch));

  return {
    harnessId,
    clock,
    model,
    tools,
    host,
    client,
    port,
    runner,
    run,
    async send(prompt, options = {}) {
      const messages: Message[] = typeof prompt === 'string' ? [user(prompt)] : [...prompt];
      if (first) {
        first = false;
        if (deferred.length > 0) messages.unshift(declareTools(deferred, true));
        if (declared.length > 0) messages.unshift(declareTools(declared));
      }
      const submitted = await client.submit({
        key: options.key ?? ids.next('key'),
        messages,
        ...(options.options !== undefined ? { options: options.options } : {}),
        ...(options.interrupt !== undefined ? { interrupt: options.interrupt } : {}),
      });
      if (submitted.launch === undefined || options.run === false) return submitted;
      return { launch: submitted.launch, exit: await run(submitted.launch, options.runner) };
    },
    async events() {
      const out: HarnessEvent[] = [];
      for (;;) {
        const page = await client.events({ after: out.at(-1)?.seq ?? 0, limit: 500 });
        out.push(...page);
        if (page.length < 500) return out;
      }
    },
    async main() {
      const execution = await client.execution();
      const [head] = await client.heads([execution.main]);
      return head === undefined ? [] : await client.entries(head.entries);
    },
    turns: () => client.turns({ limit: 100 }),
  };
}
