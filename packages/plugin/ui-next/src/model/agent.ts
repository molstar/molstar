/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '../reactive/readonly-state.js';

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
  /** A checkpoint was taken (the turn changed something) and can be restored. */
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
 * A model provider (any-model-adapters: Anthropic, Google, OpenAI, or an OpenAI-compatible endpoint of one's own, e.g.
 * OpenRouter, Ollama, vLLM). Keys stay in this browser (encrypted at rest) and go only to that provider.
 */
export interface AIProvider {
  id: string;
  label: string;
  /** The adapter family: how it is called. */
  family: 'anthropic' | 'google' | 'openai';
  /** One's own: an OpenAI-compatible endpoint (https, or http on localhost). */
  baseURL?: string;
  /** Where the key comes from: this browser, or the server's environment (none: no key yet). */
  key?: 'browser' | 'server';
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
  setKey(provider: string, key: string): void;
  removeKey(provider: string): void;
  choose(setup: AISetup): void;
  /** One's own OpenAI-compatible endpoint (name, URL, which API, key). */
  addProvider(own: { label: string; baseURL: string; api: 'chat' | 'responses'; key: string; model: string }): void;
}

/**
 * The assistant and agents at work (over any-harness and the command registry): asking, the running turn, the
 * history of every run, undoing a turn back to its checkpoint, and the model setup.
 */
export interface AgentModel {
  /** The model setup; without a model, asking opens it. */
  settings: AISettingsModel;
  /** The turn on show: running, waiting for the person, or just finished (until dismissed). */
  turn: ReadonlyState<AgentTurn | undefined>;
  /** Newest last. */
  history: ReadonlyState<readonly HistoryItem[]>;
  /** Hands a request to the agent (several steps, a turn of its own). */
  ask(prompt: string): void;
  /** Stops the running turn (what ran stays, and can be undone with the turn). */
  stop(): void;
  /** Answers a run waiting for approval. */
  approve(runId: string, allow: boolean): void;
  /** Restores the checkpoint before the turn (says what else it discards: the person's runs after it). */
  undo(turnId: string): void;
  /** Hides the finished turn's card (it stays in the history). */
  dismiss(): void;
}
