/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject, ModelFinish } from '../messages/index.js';

/**
 * What the Inspector shows (trace.md), prepared from the harness records (harness.md §3, §13): plain
 * data, no functions, so a view can be built on mock data and a host can keep it in a store. Times are milliseconds
 * from the session's first step event; tokens are the model's own counts (estimates where marked).
 */

/**
 * A turn's state. A turn that ended maps from the harness's status: `completed` is `done`, `failed` is `failed`,
 * `stopped` and `reset` are `stopped`. A turn still active is read from the Execution (harness.md §4–§6): `waiting`
 * when it is parked (a permission answered `park`); else by the step it will run next: `thinking` for a model step,
 * `working` for a tools step.
 */
export type TurnStatus = 'thinking' | 'working' | 'waiting' | 'done' | 'stopped' | 'failed';

/** The kinds of count a turn's budget limits (harness.md §11), `cost` and `endBy` aside. */
export type BudgetKind =
  | 'modelCalls'
  | 'toolCalls'
  | 'uncachedInput'
  | 'cacheRead'
  | 'cacheWrite'
  | 'output'
  | 'reasoning';

export interface TraceData {
  readonly harnessId: string;
  readonly execution: ExecutionView;
  /** The latest turn against its budget: the kind nearest its ceiling. */
  readonly budget: BudgetView;
  /** The guidance and tool definitions every model call starts with. */
  readonly preamble: PreambleView;
  readonly turns: readonly TurnView[];
  /**
   * The tools the conversation's latest context has loaded (a definition in full), each with the turn that first
   * showed it to a model (1-based; absent: from the start).
   */
  readonly loaded: readonly LoadedToolView[];
  /** The tools it declares but has not loaded (a deferred definition: a name only), by id. */
  readonly declared: readonly string[];
}

export interface ExecutionView {
  /** The event log's last `seq`. */
  readonly seq: number;
  /** The head `main` points at. */
  readonly headId: string;
  /** The step the harness will run or resume next. */
  readonly pendingStep?: 'model' | 'tools';
  /** The active turn is parked, waiting for an answer. */
  readonly parked?: boolean;
}

/** A turn's budget: the kind of count nearest its ceiling, against that ceiling. */
export interface BudgetView {
  readonly consumed: number;
  /** 0 when no ceiling is known. */
  readonly ceiling: number;
  readonly kind?: BudgetKind;
  /** The kind in words: "uncached input". */
  readonly dimension?: string;
  /** What is counted: "latest turn". */
  readonly scope?: string;
  readonly credits?: number;
}

export interface PreambleView {
  readonly guidance: string;
  readonly guidanceTokens: number;
  readonly tools: readonly ToolDefinitionView[];
  readonly toolsTokens: number;
}

export interface ToolDefinitionView {
  readonly id: string;
  readonly name: string;
  readonly description: string;
}

export interface LoadedToolView {
  readonly id: string;
  readonly name: string;
  readonly turn?: number;
}

/**
 * A call's or turn's tokens. `cached` and `written` are parts of `input` (read from and written to the provider's
 * prompt cache), where the provider reports them.
 */
export interface TokensView {
  readonly input: number;
  readonly output: number;
  readonly cached?: number;
  readonly written?: number;
  /** Part of `output`, where the provider reports it. */
  readonly reasoning?: number;
  /** The harness counted them itself (the provider reported nothing). */
  readonly estimated?: boolean;
  /** What it cost, in the host's credits (a display unit, never a limit). */
  readonly credits?: number;
}

export interface TurnView {
  readonly id: string;
  /** The request as the person wrote it. */
  readonly prompt: string;
  /** What the request said was in focus. */
  readonly focus?: string;
  readonly status: TurnStatus;
  readonly start: number;
  readonly end?: number;
  readonly tokens: TokensView;
  /** Every call of the turn, nested runs' included, in the order they began. */
  readonly steps: readonly StepView[];
  /** Sparse harness milestones (events). */
  readonly milestones: readonly MilestoneView[];
  readonly failure?: string;
}

export interface MilestoneView {
  readonly position: number;
  readonly at: number;
  readonly kind: string;
  readonly detail?: string;
}

export interface AttemptView {
  /** 1-based. */
  readonly ordinal: number;
  readonly start: number;
  readonly end?: number;
  readonly status: 'running' | 'succeeded' | 'failed' | 'aborted';
  readonly failure?: string;
  /** Streaming model calls: when the first output arrived, where a host records it. */
  readonly firstOutput?: number;
  /** A wait before this attempt (the back-off of a retry). */
  readonly waitBefore?: { readonly reason: string; readonly ms: number };
}

interface StepBase {
  readonly id: string;
  readonly callId: string;
  readonly attempts: readonly AttemptView[];
  /** The call it ran under (a tool's helper, compaction's summary: a nested run), when it did. */
  readonly under?: string;
}

export interface ModelStepView extends StepBase {
  readonly kind: 'model';
  /** The model that was called. */
  readonly model?: string;
  /** The last attempt's usage. */
  readonly tokens?: TokensView;
  /** The model's text, if it wrote any. */
  readonly reply?: string;
  /** The tools it asked for, as "name(args)". */
  readonly toolCalls?: readonly string[];
  /**
   * When its request does not extend the previous one (so the provider's cache stops there): why, and how many of
   * the previous request's parts it kept.
   */
  readonly prefixBreak?: PrefixBreak;
  /** How its answer ended: `tool-calls` or `complete` (the log does not keep the other endings). */
  readonly finish?: ModelFinish;
  /** Its input: what it saw. */
  readonly window?: WindowView;
}

export interface PrefixBreak {
  readonly reason: string;
  readonly kept: number;
  readonly of: number;
}

export interface ToolStepView extends StepBase {
  readonly kind: 'tool';
  readonly toolId: string;
  readonly name: string;
  /** A short label for the call, from the host's profile ("Color chain A"). */
  readonly title?: string;
  /** It only read: from the host's profile. */
  readonly read?: boolean;
  readonly args: JsonObject;
  readonly status: 'running' | 'succeeded' | 'failed' | 'aborted';
  readonly result?: string;
  /** What this call changed, from the host's profile, for cross-highlighting. */
  readonly changed?: readonly { readonly id: string; readonly label: string }[];
}

export type StepView = ModelStepView | ToolStepView;

/** The parts of what a model sees: the guidance, the tool definitions, the host's state, the conversation. */
export type WindowPart = 'guidance' | 'tools' | 'state' | 'conversation';

export interface WindowEntryView {
  readonly role: 'system' | 'tools' | 'user' | 'assistant' | 'tool';
  readonly part: WindowPart;
  readonly title: string;
  readonly text: string;
  /** An estimate (characters / 4). */
  readonly tokens: number;
  /** Appended since the run's previous model call. */
  readonly fresh?: boolean;
  readonly failed?: boolean;
}

export interface WindowView {
  readonly entries: readonly WindowEntryView[];
  /** The model's context limit, for scale, where the host gave it. */
  readonly limit?: number;
}

/* The Tools view. */

/** Where a tool loads: with the first turn, with its skill, or by name. */
export type ToolGroupKind = 'core' | 'skill' | 'single';

/**
 * A tool in the conversation: `loaded` in full in its latest context, `deferred` (declared by name, loaded when
 * called or disclosed), or `offered` (not in the context yet).
 */
export type ToolState = 'offered' | 'deferred' | 'loaded';

export interface ToolView {
  readonly id: string;
  readonly name: string;
  readonly description: string;
  /** Its input, as the JSON Schema the model sees. */
  readonly input: JsonObject;
  /** The API modules it declares (it gets only those). */
  readonly modules?: readonly string[];
  /** It only reads: its result depends on its arguments and the host's state alone. */
  readonly read: boolean;
  /** It needs the person's approval before it runs. */
  readonly approval: boolean;
  /** The turn's tool set includes it. */
  readonly enabled: boolean;
  readonly state: ToolState;
  /** The turn that first showed it loaded to a model (1-based; absent: from the start, or not loaded). */
  readonly loadedAtTurn?: number;
  /** Its definition's size in tokens (an estimate). */
  readonly tokens: number;
  /** One of a skill's extras: its guidance names it, and it loads by name or when called. */
  readonly extra?: boolean;
  /** Why it cannot be used here. */
  readonly unavailable?: string;
}

export interface ToolGroupView {
  /** "core", "single", or the skill's id. */
  readonly id: string;
  readonly kind: ToolGroupKind;
  readonly title: string;
  /** A skill's one line: what it is for. */
  readonly summary?: string;
  /** A skill's guidance, exactly as the model receives it when it loads the skill. */
  readonly guidance?: string;
  readonly guidanceTokens?: number;
  /** Why a skill cannot be used here (a module this host does not provide). */
  readonly unavailable?: string;
  readonly tools: readonly ToolView[];
}

export interface ToolsData {
  readonly groups: readonly ToolGroupView[];
}

/* The Limits view. */

/** One limit, in the unit people read (tokens, MB, seconds, %). */
export interface LimitRowView {
  readonly key: string;
  readonly label: string;
  /** What it bounds, in a phrase. */
  readonly note?: string;
  readonly unit: string;
  readonly value: number;
  readonly default: number;
  /** What it counts now (the latest turn's consumption), drawn against `value`. */
  readonly used?: number;
  /** The most the host lets `value` take. */
  readonly ceiling?: number;
  /** The smallest step the field takes (default 1). */
  readonly step?: number;
  /** A limit picked from a few levels (the model's thinking level): each value with its name. */
  readonly choices?: readonly { readonly value: number; readonly label: string }[];
  readonly editable: boolean;
}

/** A group of limits ("Turn budget", "Context", "Web"). */
export interface LimitsGroupView {
  readonly id: string;
  readonly title: string;
  /** Where they are kept, and when a change applies. */
  readonly note?: string;
  readonly rows: readonly LimitRowView[];
}

export interface LimitsData {
  readonly groups: readonly LimitsGroupView[];
}
