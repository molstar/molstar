/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type AgentModel,
  type AgentTurn,
  type AIProvider,
  type AISettingsModel,
  type AISetup,
  type CommandRun,
  type HistoryItem,
  type KeyScope,
  type RunSource,
  type RunStatus,
} from '../../model/agent.js';
import { observableState } from '../../reactive/observable-state.js';
import type { CommandMenuAssistant } from '../../views/command-menu/command-menu.js';

/** One scripted step of a mock turn. */
interface Step {
  command: string;
  title: string;
  /** How long it runs (ms); with progress, it reports it. */
  ms: number;
  read?: boolean;
  progress?: boolean;
  approval?: string;
  detail?: string;
  fail?: string;
}

interface Script {
  thinking: string;
  steps: Step[];
  answer: string;
  retries?: number;
}

const Agent: RunSource = { kind: 'agent', name: 'Mol* assistant' };
const Person: RunSource = { kind: 'person' };

/** Scripts picked by words in the request; anything else gets the active-site one. */
function scriptFor(prompt: string): Script {
  const p = prompt.toLowerCase();
  if (/superpose|align|load|download|fetch/.test(p)) {
    return {
      thinking: 'I will fetch 1W0E and superpose it on 1TQN by their A chains.',
      retries: 1,
      steps: [
        { command: 'scene.get', title: 'Read the scene', ms: 300, read: true, detail: '1 structure: 1TQN' },
        {
          command: 'data.download-structure',
          title: 'Download 1W0E from PDBe',
          ms: 1600,
          progress: true,
          approval: 'Fetches data from www.ebi.ac.uk',
          detail: '1W0E: 1 model, 3,898 atoms',
        },
        {
          command: 'superpose.chains',
          title: 'Superpose 1W0E on 1TQN (chain A)',
          ms: 900,
          detail: 'RMSD 0.38 Å over 412 Cα',
        },
        { command: 'camera.reset', title: 'Reset the camera', ms: 300 },
      ],
      answer: 'Loaded 1W0E and superposed it on 1TQN by chain A (RMSD 0.38 Å over 412 Cα atoms).',
    };
  }
  const color = /\b(red|blue|green|yellow|orange|purple|white|black|gr[ae]y)\b/.exec(p)?.[1];
  if (/\b(colou?r|paint)\b/.test(p) && color) {
    return {
      thinking: `Coloring every representation of every structure ${color}, uniformly.`,
      steps: [
        { command: 'scene.get', title: 'Read the scene', ms: 300, read: true, detail: '1TQN: 4 representations' },
        {
          command: 'structure.set-color',
          title: `Color all representations ${color} (uniform)`,
          ms: 1200,
          progress: true,
          detail: '4 representations in 1TQN',
        },
      ],
      answer: `Colored all 4 representations of 1TQN ${color}. Undo turn brings the colors back.`,
    };
  }
  if (/fail|broken|error/.test(p)) {
    return {
      thinking: 'Trying to color by the custom property.',
      steps: [
        { command: 'scene.get', title: 'Read the scene', ms: 300, read: true },
        {
          command: 'structure.set-color',
          title: 'Color by B-factor (custom range)',
          ms: 600,
          fail: 'No B-factor values in the selection',
        },
      ],
      answer: 'The selection has no B-factor values, so I could not color by them. Try a different theme?',
    };
  }
  return {
    thinking: 'Looking at the heme pocket: I will select the residues around HEM 508 and give them a component.',
    steps: [
      { command: 'scene.get', title: 'Read the scene', ms: 300, read: true, detail: '1TQN: polymer, ligand, water' },
      {
        command: 'select.residues',
        title: 'Select residues within 5 Å of HEM 508',
        ms: 500,
        detail: '23 residues in chain A',
      },
      { command: 'structure.component.add', title: 'Add component “Active site” from the selection', ms: 700 },
      {
        command: 'structure.add-representation',
        title: 'Show “Active site” as ball & stick',
        ms: 900,
        progress: true,
      },
      { command: 'camera.focus', title: 'Focus the view on “Active site”', ms: 400 },
    ],
    answer:
      'Made an “Active site” component (23 residues around HEM 508), shown as ball & stick, and focused the view on it.',
  };
}

/** The providers as Stories untold has them (any-model-adapters), none with a key yet. */
const Providers: readonly AIProvider[] = [
  {
    id: 'anthropic',
    label: 'Anthropic',
    family: 'anthropic',
    models: ['claude-sonnet-5-5', 'claude-opus-5-5', 'claude-haiku-5-5'],
  },
  { id: 'google', label: 'Gemini', family: 'google', models: ['gemini-3.8-flash', 'gemini-3.8-pro'] },
  { id: 'openai', label: 'OpenAI', family: 'openai', models: ['gpt-6-luna', 'gpt-6-luna-mini'] },
];

/** Keys and the model choice, in memory (the real one keeps keys encrypted in IndexedDB). */
export class MockAISettings implements AISettingsModel {
  readonly #providers = observableState<readonly AIProvider[]>(Providers);
  readonly #setup = observableState<AISetup | undefined>(undefined);
  readonly providers = this.#providers.asReadonly;
  readonly setup = this.#setup.asReadonly;

  constructor(
    private readonly post: (text: string, kind?: 'info' | 'success' | 'warning' | 'error') => void,
    configured?: boolean,
  ) {
    if (configured) {
      this.#patch('anthropic', { key: { from: 'here', scope: 'browser' } });
      this.#setup.set({ provider: 'anthropic', model: 'claude-sonnet-5-5' });
    }
  }

  setKey(provider: string, key: string, scope: KeyScope) {
    if (!key.trim()) return;
    this.#patch(provider, { key: { from: 'here', scope } });
    const p = this.#providers.value.find((x) => x.id === provider)!;
    if (!this.#setup.value) this.#setup.set({ provider, model: p.models[0] });
    const kept = { reload: 'until the page reloads', session: 'for this session', browser: 'in this browser' }[scope];
    this.post(`Key saved for ${p.label}, kept ${kept}`, 'success');
  }

  removeKey(provider: string) {
    this.#patch(provider, { key: undefined });
    if (this.#setup.value?.provider === provider) {
      const other = this.#providers.value.find((p) => p.key);
      this.#setup.set(other ? { provider: other.id, model: other.models[0] } : undefined);
    }
  }

  choose(setup: AISetup) {
    this.#setup.set(setup);
  }

  async refreshModels(provider?: string) {
    await new Promise((r) => setTimeout(r, 700));
    // What the providers' /models would add (a new release since the list was made).
    const more: Record<string, string[]> = {
      anthropic: ['claude-fable-5-1'],
      google: ['gemini-3.8-flash-lite'],
      openai: ['gpt-6-luna-nano'],
    };
    this.#providers.set(
      this.#providers.value.map((p) =>
        p.key && (!provider || p.id === provider) && more[p.id]
          ? { ...p, models: [...new Set([...p.models, ...more[p.id]])] }
          : p,
      ),
    );
    this.post('Model lists updated from the providers', 'success');
  }

  async listModels(endpoint: { baseURL: string; key: string }) {
    await new Promise((r) => setTimeout(r, 700));
    if (!/^https?:\/\//.test(endpoint.baseURL)) throw new Error('Not a URL');
    return ['llama-4-scout', 'qwen3-coder-next', 'mistral-large-3'];
  }

  addProvider(own: {
    label: string;
    baseURL: string;
    api: 'chat' | 'responses';
    key: string;
    scope: KeyScope;
    models: readonly string[];
  }) {
    const id = `custom:${own.label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`;
    this.#providers.set([
      ...this.#providers.value,
      {
        id,
        label: own.label,
        family: 'openai',
        baseURL: own.baseURL,
        key: own.key ? { from: 'here', scope: own.scope } : undefined,
        models: own.models,
      },
    ]);
    this.#setup.set({ provider: id, model: own.models[0] });
    this.post(`Added ${own.label}`, 'success');
  }

  #patch(id: string, patch: Partial<AIProvider>) {
    this.#providers.set(this.#providers.value.map((p) => (p.id === id ? { ...p, ...patch } : p)));
  }
}

const wait = (ms: number, signal: AbortSignal) =>
  new Promise<void>((resolve, reject) => {
    const t = setTimeout(resolve, ms);
    signal.addEventListener('abort', () => {
      clearTimeout(t);
      reject(new Error('stopped'));
    });
  });

/**
 * A scripted agent for the stories: streams its words, runs steps one by one (reads quietly, progress where the command
 * reports it), waits for approval on downloads, fails when asked to, can be stopped, and undoes a turn to its checkpoint.
 */
export class MockAgentModel implements AgentModel {
  readonly #turn = observableState<AgentTurn | undefined>(undefined);
  readonly #view = observableState<'card' | 'line' | 'none'>('none');
  readonly #history = observableState<readonly HistoryItem[]>([]);
  #serial = 0;
  #abort: AbortController | undefined;
  #approval: ((allow: boolean) => void) | undefined;

  readonly settings: MockAISettings;
  readonly turn = this.#turn.asReadonly;
  readonly view = this.#view.asReadonly;
  readonly history = this.#history.asReadonly;

  constructor(
    private readonly post: (text: string, kind?: 'info' | 'success' | 'warning' | 'error') => void,
    options: { configured?: boolean; history?: boolean } = {},
  ) {
    this.settings = new MockAISettings(post, options.configured);
    if (options.history) this.#seed();
  }

  /** A command the person ran (the command menu, a button): one history row. */
  record(title: string, command: string, detail?: string) {
    this.#release();
    const run: CommandRun = {
      id: `run-${++this.#serial}`,
      command,
      title,
      source: Person,
      status: 'done',
      detail,
      time: new Date(),
    };
    this.#history.set([...this.#history.value, { kind: 'run', run }]);
  }

  ask(prompt: string) {
    if (this.#abort) this.stop();
    this.#release();
    const script = scriptFor(prompt);
    const turn: AgentTurn = {
      id: `turn-${++this.#serial}`,
      prompt,
      source: Agent,
      status: 'thinking',
      text: '',
      runs: [],
      checkpoint: false,
      time: new Date(),
    };
    this.#show(turn);
    void this.#play(turn.id, script);
  }

  stop() {
    this.#abort?.abort();
  }

  approve(runId: string, allow: boolean) {
    const t = this.#turn.value;
    if (!t?.runs.some((r) => r.id === runId && r.status === 'waiting')) return;
    this.#approval?.(allow);
  }

  undo(turnId: string) {
    const item = this.#history.value.find((h) => h.kind === 'turn' && h.turn.id === turnId);
    if (!item || item.kind !== 'turn' || !item.turn.checkpoint || item.turn.undone) return;
    this.#update(turnId, (t) => ({ ...t, undone: true, checkpoint: false }));
    this.post(`Undid “${item.turn.prompt}”`, 'success');
  }

  /** The next action releases the checkpoint before the last turn (one snapshot at most is kept). */
  #release() {
    if (!this.#history.value.some((h) => h.kind === 'turn' && h.turn.checkpoint)) return;
    this.#history.set(
      this.#history.value.map((h) =>
        h.kind === 'turn' && h.turn.checkpoint ? { kind: 'turn', turn: { ...h.turn, checkpoint: false } } : h,
      ),
    );
    const latest = this.#turn.value;
    if (latest?.checkpoint) this.#turn.set({ ...latest, checkpoint: false });
  }

  dismiss() {
    if (this.#view.value === 'card') this.#view.set('line');
  }

  hide() {
    this.#view.set('none');
  }

  reopen() {
    if (this.#turn.value) this.#view.set('card');
  }

  // The turn on show and its row in the history are one record, replaced on every change.
  #show(turn: AgentTurn) {
    this.#turn.set(turn);
    this.#view.set('card');
    this.#history.set([...this.#history.value, { kind: 'turn', turn }]);
  }

  #update(id: string, change: (t: AgentTurn) => AgentTurn) {
    let next: AgentTurn | undefined;
    this.#history.set(
      this.#history.value.map((h) => {
        if (h.kind !== 'turn' || h.turn.id !== id) return h;
        next = change(h.turn);
        return { kind: 'turn', turn: next };
      }),
    );
    if (next && this.#turn.value?.id === id) this.#turn.set(next);
  }

  #setRun(turnId: string, runId: string, patch: Partial<CommandRun>) {
    this.#update(turnId, (t) => ({ ...t, runs: t.runs.map((r) => (r.id === runId ? { ...r, ...patch } : r)) }));
  }

  async #play(id: string, script: Script) {
    const abort = new AbortController();
    this.#abort = abort;
    const started = Date.now();
    let usage = { input: 2400, output: 0 };
    try {
      // Thinking: the words stream in.
      for (const word of script.thinking.split(' ')) {
        await wait(45, abort.signal);
        usage = { input: usage.input, output: usage.output + 2 };
        this.#update(id, (t) => ({ ...t, text: t.text ? `${t.text} ${word}` : word, usage }));
      }
      if (script.retries) {
        await wait(400, abort.signal);
        this.#update(id, (t) => ({ ...t, retries: script.retries }));
      }
      this.#update(id, (t) => ({ ...t, status: 'working' }));
      for (const step of script.steps) {
        const run: CommandRun = {
          id: `run-${++this.#serial}`,
          command: step.command,
          title: step.title,
          source: Agent,
          status: 'running',
          read: step.read,
          approval: !!step.approval,
          time: new Date(),
        };
        // The first change takes the checkpoint.
        this.#update(id, (t) => ({ ...t, checkpoint: t.checkpoint || !step.read, runs: [...t.runs, run] }));
        if (step.approval) {
          this.#setRun(id, run.id, { status: 'waiting', detail: step.approval });
          this.#update(id, (t) => ({ ...t, status: 'waiting' }));
          const allow = await new Promise<boolean>((resolve, reject) => {
            this.#approval = resolve;
            abort.signal.addEventListener('abort', () => reject(new Error('stopped')));
          });
          this.#approval = undefined;
          if (!allow) {
            this.#setRun(id, run.id, { status: 'denied', detail: 'Not allowed' });
            this.#update(id, (t) => ({
              ...t,
              status: 'done',
              text: 'You did not allow the download, so I stopped there.',
              durationMs: Date.now() - started,
            }));
            this.#abort = undefined;
            return;
          }
          this.#setRun(id, run.id, { status: 'running', detail: undefined });
          this.#update(id, (t) => ({ ...t, status: 'working' }));
        }
        const t0 = Date.now();
        if (step.progress) {
          for (let i = 1; i <= 10; i++) {
            await wait(step.ms / 10, abort.signal);
            this.#setRun(id, run.id, { progress: i / 10 });
          }
        } else {
          await wait(step.ms, abort.signal);
        }
        const status: RunStatus = step.fail ? 'failed' : 'done';
        this.#setRun(id, run.id, {
          status,
          detail: step.fail ?? step.detail,
          progress: undefined,
          durationMs: Date.now() - t0,
        });
        usage = { input: usage.input + 900, output: usage.output + 60 };
        this.#update(id, (t) => ({ ...t, usage }));
        if (step.fail) break;
      }
      await wait(200, abort.signal);
      this.#update(id, (t) => ({
        ...t,
        status: script.steps.some((s) => s.fail) ? 'failed' : 'done',
        text: script.answer,
        durationMs: Date.now() - started,
      }));
    } catch {
      this.#update(id, (t) => ({
        ...t,
        status: 'stopped',
        runs: t.runs.map((r) => (r.status === 'running' || r.status === 'waiting' ? { ...r, status: 'stopped' } : r)),
        durationMs: Date.now() - started,
      }));
    } finally {
      if (this.#abort === abort) this.#abort = undefined;
    }
  }

  /** A history to look at: the person's commands and a finished turn between them. */
  #seed() {
    const at = (min: number) => new Date(Date.now() - min * 60_000);
    const run = (title: string, command: string, min: number, extra: Partial<CommandRun> = {}): CommandRun => ({
      id: `run-${++this.#serial}`,
      command,
      title,
      source: Person,
      status: 'done',
      time: at(min),
      ...extra,
    });
    const turn: AgentTurn = {
      id: `turn-${++this.#serial}`,
      prompt: 'show the active site around the heme',
      source: Agent,
      status: 'done',
      text: 'Made an “Active site” component (23 residues around HEM 508), shown as ball & stick, and focused the view on it.',
      runs: [
        { ...run('Read the scene', 'scene.get', 6), source: Agent, read: true, durationMs: 120 },
        { ...run('Select residues within 5 Å of HEM 508', 'select.residues', 6), source: Agent, detail: '23 residues' },
        { ...run('Add component “Active site” from the selection', 'structure.component.add', 6), source: Agent },
        { ...run('Show “Active site” as ball & stick', 'structure.add-representation', 6), source: Agent },
        { ...run('Focus the view on “Active site”', 'camera.focus', 6), source: Agent },
      ],
      // Released: commands came after it.
      checkpoint: false,
      usage: { input: 6100, output: 240 },
      time: at(6),
      durationMs: 4200,
    };
    this.#history.set([
      { kind: 'run', run: run('Download 1TQN from PDBe', 'data.download-structure', 9, { detail: '1 model' }) },
      { kind: 'run', run: run('Apply preset: Polymer & ligand', 'preset.polymer-and-ligand', 8) },
      { kind: 'turn', turn },
      { kind: 'run', run: run('Color the selection red', 'selection.color', 3) },
      { kind: 'run', run: run('Take snapshot “Overview”', 'snapshot.add', 1) },
    ]);
  }
}

/**
 * The command menu's assistant over the mock agent: a few kinds of words read into one command (its sentence opens
 * filled), anything else handed to the agent as a turn.
 */
export function mockAssistant(agent: MockAgentModel, openSetup: () => void): CommandMenuAssistant {
  const colors = /\b(red|blue|green|yellow|orange|purple|white|black|gr[ae]y)\b/;
  return {
    ready: () => !!agent.settings.setup.value,
    propose(text) {
      const t = text.toLowerCase();
      const color = colors.exec(t);
      // The selection, in a color: one command. Everything, all, or more than one thing: a turn for the agent.
      if (/\b(colou?r|paint|make)\b/.test(t) && color && !/\b(and|everything|all|each|every)\b/.test(t))
        return { command: 'selection.color', tail: color[1] };
      const id = /\b(?:fetch|download|load|open)\s+(?:pdb\s+)?([0-9][a-z0-9]{3})\b/.exec(t);
      if (id && !/\b(and|superpose|align)\b/.test(t)) return { command: 'data.download-structure', tail: id[1] };
      const focus = /^(?:focus|zoom)(?: in)?(?: on)?\s+(.+)$/.exec(t);
      if (focus) return { command: 'focus.target', tail: focus[1] };
      return undefined;
    },
    ask: (text) => agent.ask(text),
    setup: openSetup,
  };
}
