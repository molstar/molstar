/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ReasoningParameters, Reference } from '@molstar/ai-agent';
import type { CompactionReport } from './compaction/index.js';

/** What the hooks and results of a session carry (tools-and-session.md §4). */

export interface ToolEvent {
  /** The call it is about: a started call, its progress and its end share it. */
  id?: string;
  name: string;
  args: Json;
  /** skipped: never ran, as a call before it in the same response failed (the drive stopped there). */
  status: 'started' | 'progress' | 'succeeded' | 'failed' | 'skipped';
  text?: string;
  /** What a long call is doing ("RCSB FATCAT-flexible: running, 6 s"), while it runs. */
  progress?: string;
  /** The tool only reads (its registration says so). */
  read?: boolean;
  /** When the call started (the first attempt), by the session's clock. */
  startedAt: number;
  /** When it ended (succeeded, failed or skipped). */
  endedAt?: number;
}

/** What the session is doing, for a status line. */
export interface StatusEvent {
  kind: 'model' | 'tools' | 'permission' | 'budget' | 'retrying' | 'widened' | 'budget-extended' | 'compacting';
  /** retrying: why; budget-extended: what grew. */
  detail?: string;
}

/** Where a turn is, for a host that must act at a turn's edge. */
export interface TurnEvent {
  /**
   * The turn has begun and its id is known: before its first model call, so before any change it makes (where a
   * host takes its checkpoint). `resumed`: the turn began earlier and this session took it up.
   */
  kind: 'begun';
  turnId: string;
  resumed?: boolean;
}

export interface TurnHooks {
  /** The assistant's text so far (streaming). */
  onText?(text: string): void;
  onTool?(event: ToolEvent): void;
  onStatus?(event: StatusEvent): void;
  onTurn?(event: TurnEvent): void;
}

/** A call that waits for the person's approval (a tool marked `approval`). */
export interface PendingApproval {
  callId: string;
  tool: string;
  args: Json;
}
export type ApprovalDecision = 'allow' | 'deny';

/** How a parked or interrupted turn goes on (`AgentSession.resume`). */
export interface ResumeOptions {
  /** The person's answers to the calls that wait, by call id; none: an interrupted turn is taken up as it is. */
  readonly decisions?: Readonly<Record<string, ApprovalDecision>>;
  readonly hooks?: TurnHooks;
}

/** What a prompt carries beyond its words (tools-and-session.md §4); each part is the host's text. */
export interface TurnRequest<M extends ModelChoice = ModelChoice> {
  prompt: string;
  /** The state the model sees with every prompt: compact, ids first. */
  state?: string;
  /** What the person points at; the request has a Focus line only when this is given. */
  focus?: string;
  /** What changed since the agent's last turn that the person did (a list of lines). */
  notes?: readonly string[];
  /** Ahead of the state (a memory, a digest of earlier turns): the host gives it when the conversation needs it. */
  preamble?: string;
  /** Files sent with the request: references that carry their bytes. */
  attachments?: readonly Reference[];
  /** The model for this turn; default: the session's own. */
  model?: M;
}

export interface TurnResult {
  /**
   * awaiting-approval: the turn is parked at its permission decision (resume() with decisions goes on). reset: a reset
   * ended it (`forget()`, or another client's).
   */
  status: 'completed' | 'stopped' | 'failed' | 'reset' | 'awaiting-approval';
  /** awaiting-approval: the calls the person is asked about. */
  approval?: PendingApproval[];
  /** The last thing the model said (its last model call with text). */
  text: string;
  /**
   * All the model said in the turn, one reply per model call (a retried attempt counts once, a reply asked for again
   * with more room replaces the one cut short), in order, separated by blank lines. Since this session took the turn
   * up, for a resumed one.
   */
  turnText: string;
  error?: string;
  tools: ToolEvent[];
  /** The turn's tokens: input (with the parts read from and written to the provider's prompt cache,
   * where reported) and output; model calls counted. */
  usage?: { input: number; output: number; cached: number; written: number; modelCalls: number };
  /** The conversation's size after the turn, in tokens (the provider's count of its last call). */
  context?: number;
  /** Compactions during the turn (and the one before it, at the soft threshold): what each kept. */
  compactions?: CompactionReport[];
}

/** Which model runs a turn, and how much it thinks. */
export interface ModelChoice {
  provider: string;
  model: string;
  /** The model's thinking, over `AgentSessionOptions.thinking`; "model" leaves it to the model. */
  thinking?: ReasoningParameters | 'model';
  /** Its context window in tokens (the application's catalog, or what the person gave for a custom model). */
  contextWindow?: number;
  /** Its output limit in tokens, for a model the adapter does not know (default 32k). */
  maxOutput?: number;
}
