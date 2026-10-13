/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui';

/** Who ran a command (design/agent.md §4): shown in the history, so an agent's work is never mistaken for the person's. */
export type RunSource =
  | { kind: 'person' }
  | { kind: 'assistant' }
  | { kind: 'agent'; name: string }
  | { kind: 'script'; name?: string };

export type RunStatus = 'running' | 'waiting' | 'done' | 'failed' | 'stopped' | 'denied';

/** One command run (a tool call, for an agent): what the history and a turn's card list. */
export interface CommandRun {
  id: string;
  /** The command's id (the tool's name). */
  command: string;
  /** What it does, in words, with its arguments ("Color the selection red", "Download 1W0E from PDBe"). */
  title: string;
  source: RunSource;
  status: RunStatus;
  /** A read (no effect): shown quieter, needs no approval or checkpoint. */
  read?: boolean;
  /** 0..1 while running, when known. */
  progress?: number;
  /** Its outcome ("12 residues selected") or error; while waiting, why it needs approval. */
  detail?: string;
  /** Leaves the page or cannot be undone: it waited for the person's Allow. */
  approval?: boolean;
  time: Date;
  durationMs?: number;
}

export type TurnStatus = 'thinking' | 'working' | 'waiting' | 'done' | 'stopped' | 'failed';

/**
 * An agent turn: one request and everything done for it (design/agent.md §5). A checkpoint is taken before its first
 * change; undoing the turn restores it (data, camera, selection), however many commands ran.
 */
export interface AgentTurn {
  id: string;
  prompt: string;
  source: RunSource;
  status: TurnStatus;
  /** The model's words so far (streamed), and its last word at the end. */
  text: string;
  runs: readonly CommandRun[];
  /**
   * Its checkpoint is held (the turn changed something, and nothing came after it): Undo turn restores it. The next
   * action (the person's command, another turn) releases it, so one snapshot at most is ever kept.
   */
  checkpoint: boolean;
  /** Undone: back at its checkpoint. */
  undone?: boolean;
  /** Model calls retried (a failed request, a rate limit), shown so a slow turn is explained. */
  retries?: number;
  /** Tokens in and out, when the provider reports them. */
  usage?: { input: number; output: number };
  time: Date;
  durationMs?: number;
}

/** The history: the person's commands one by one, an agent's grouped by turn. */
export type HistoryItem = { kind: 'run'; run: CommandRun } | { kind: 'turn'; turn: AgentTurn };

/**
 * How long a key is kept here: until the page reloads (memory only), for this session (the tab), or in this browser
 * (until removed). Kept ones are stored scrambled (encrypted with a key this browser keeps), so the storage read as it
 * is gives no working key; that does not guard against code running in the page.
 */
export type KeyScope = 'reload' | 'session' | 'browser';

/**
 * A model provider (any-model-adapters: Anthropic, Google, OpenAI, or an OpenAI-compatible endpoint of one's own, e.g.
 * OpenRouter, Ollama, vLLM). Keys go only to their provider.
 */
export interface AIProvider {
  id: string;
  label: string;
  /** The adapter family: how it is called. */
  family: 'anthropic' | 'google' | 'openai';
  /** One's own: an OpenAI-compatible endpoint (https, or http on localhost). */
  baseURL?: string;
  /** Its key: entered here (kept for a scope), or from the server's environment; none yet when absent. */
  key?: { from: 'here'; scope: KeyScope } | { from: 'server' };
  /** Models offered for it (the default first). */
  models: readonly string[];
}

/** The model the assistant uses (from the next request on). */
export interface AISetup {
  provider: string;
  model: string;
}

/** Setting up the assistant's model: providers, their keys, the choice, one's own endpoints. */
export interface AISettingsModel {
  providers: ReadonlyState<readonly AIProvider[]>;
  /** The choice; undefined while no provider has a key. */
  setup: ReadonlyState<AISetup | undefined>;
  setKey(provider: string, key: string, scope: KeyScope): void;
  removeKey(provider: string): void;
  choose(setup: AISetup): void;
  /** Asks the providers with a key (or one) which models they offer now (their `/models`), and offers those. */
  refreshModels(provider?: string): Promise<void>;
  /** Asks an endpoint not added yet which models it offers (OpenAI-compatible `/models`). */
  listModels(endpoint: { baseURL: string; key: string }): Promise<readonly string[]>;
  /** One's own OpenAI-compatible endpoint (name, URL, which API, key, its models). */
  addProvider(own: {
    label: string;
    baseURL: string;
    api: 'chat' | 'responses';
    key: string;
    scope: KeyScope;
    models: readonly string[];
  }): void;
}

/**
 * The assistant and agents at work (over any-harness and the command registry): asking, the running turn, the
 * history of every run, undoing a turn back to its checkpoint, and the model setup.
 */
export interface AgentModel {
  /** The model setup; without a model, asking opens it. */
  settings: AISettingsModel;
  /** The latest turn: running, waiting for the person, or the last one finished. */
  turn: ReadonlyState<AgentTurn | undefined>;
  /**
   * How it shows, under the status line's tasks: its card (asking opens it), a line (Esc or ✕ closes the card to it; a
   * click reopens it), or nothing (the line's ✕; it stays in History).
   */
  view: ReadonlyState<'card' | 'line' | 'none'>;
  /** Newest last. */
  history: ReadonlyState<readonly HistoryItem[]>;
  /** Hands a request to the agent (several steps, a turn of its own). */
  ask(prompt: string): void;
  /** Stops the running turn (what ran stays, and can be undone with the turn). */
  stop(): void;
  /** Answers a run waiting for approval. */
  approve(runId: string, allow: boolean): void;
  /** Restores the checkpoint before the turn (only the latest turn has one, until the next action). */
  undo(turnId: string): void;
  /** Closes the turn's card to its line (a running turn goes on). */
  dismiss(): void;
  /** Takes the line away too (the turn stays in History). */
  hide(): void;
  /** Opens the latest turn's card again. */
  reopen(): void;
}
