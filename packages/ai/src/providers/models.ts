/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReasoningParameters } from '../messages/index.js';

/**
 * The model catalog (plan.md §2.5): what the providers' own models are, as data, with no dependencies, so a host's
 * model list reads it without loading a provider. The providers take their limits and quirks from it; a host takes the
 * labels, the reasoning default and the prices. Endpoints of the host's own (a gateway, a plan) are the host's.
 */

export interface ModelLimits {
  readonly contextWindowTokens: number;
  /** An input cap for providers whose window does not include the output. */
  readonly maxInputTokens?: number;
  readonly maxOutputTokens: number;
  readonly defaultMaxOutputTokens?: number;
}

export interface ModelReasoning {
  /** The modes the model is known to take (`off` aside). */
  readonly modes: readonly ('adaptive' | 'budget')[];
  /** The adaptive efforts known to work, least first. */
  readonly efforts?: readonly string[];
  /** What a host sends when the person chose nothing; absent, nothing (the provider's own). */
  readonly default?: ReasoningParameters;
}

/** A list price in dollars per million tokens, as of a date. */
export interface ModelPrice {
  /** The date the price was read (YYYY-MM-DD). */
  readonly asOf: string;
  readonly source: string;
  readonly input: number;
  readonly cachedInput: number;
  /** Cache writes, where the provider bills them; else as input. */
  readonly cacheWrite?: number;
  readonly output: number;
  /** A prompt over `tokens` costs `factor` times as much, all of it. */
  readonly longPrompt?: { readonly tokens: number; readonly factor: number };
}

export interface ModelQuirks {
  /** Tools loaded later are appended where they load (Anthropic `tool_addition`, OpenAI `additional_tools`). */
  readonly appendsTools?: boolean;
}

/**
 * An entry: a family (`family: true`, matched by prefix at a boundary: `gpt-6` matches `gpt-6.1-sol`, not `gpt-60`)
 * or a model (matched as itself and its dated snapshots: `claude-haiku-4-5-20251001`). A model's fields are over its
 * families', the longer prefix over the shorter, so a new version of a family works unchanged (limits, reasoning,
 * quirks) and only a listed model has a label and a price.
 */
export interface CatalogEntry {
  readonly provider: string;
  readonly id: string;
  readonly family?: true;
  readonly label?: string;
  /** A sentence: what it is like to use. */
  readonly note?: string;
  readonly limits?: ModelLimits;
  readonly reasoning?: ModelReasoning;
  readonly quirks?: ModelQuirks;
  readonly price?: ModelPrice;
}

/** What the catalog knows of one model id, its entries merged. */
export interface ModelInfo {
  readonly provider: string;
  readonly id: string;
  /** Whether a model entry matched (else only families: a label and a price are not known). */
  readonly listed: boolean;
  readonly label?: string;
  readonly note?: string;
  readonly limits?: ModelLimits;
  readonly reasoning?: ModelReasoning;
  readonly quirks: ModelQuirks;
  readonly price?: ModelPrice;
}

// The efforts Stories' evals ran (2026-10); the APIs may take more (the providers check what the wire takes).
const EFFORTS = ['low', 'medium', 'high'] as const;
const LOW: ReasoningParameters = { mode: 'adaptive', effort: 'low' };

const GEMINI_WINDOW: ModelLimits = {
  contextWindowTokens: 1_114_112,
  maxInputTokens: 1_048_576,
  maxOutputTokens: 65_536,
};
const CLAUDE_4_5: ModelLimits = { contextWindowTokens: 200_000, maxOutputTokens: 64_000 };
const CLAUDE_5: ModelLimits = { contextWindowTokens: 1_000_000, maxOutputTokens: 128_000 };
// Adaptive thinking with an effort (Anthropic's models API, 2026-10-08).
const CLAUDE_5_REASONING: ModelReasoning = { modes: ['adaptive'], efforts: EFFORTS };

export const MODEL_CATALOG: readonly CatalogEntry[] = [
  // Google. Reasoning parameters from Gemini 3 (`google/provider.ts`). Low by default: in Stories' evals it did as
  // well as more thinking, 2–5× faster (more made Gemini read and re-read).
  { provider: 'google', id: 'gemini-2.5', family: true, limits: GEMINI_WINDOW },
  {
    provider: 'google',
    id: 'gemini-3',
    family: true,
    limits: GEMINI_WINDOW,
    reasoning: { modes: ['adaptive', 'budget'], efforts: ['minimal', ...EFFORTS], default: LOW },
  },
  {
    provider: 'google',
    id: 'gemini-3.8-flash',
    label: 'Gemini 3.8 Flash',
    note: 'Fast.',
    // Through 2026; doubles in 2027.
    price: { asOf: '2026-09-30', source: 'ofox.ai', input: 0.75, cachedInput: 0.075, output: 3.75 },
  },

  // Anthropic.
  { provider: 'anthropic', id: 'claude-haiku-4-5', family: true, limits: CLAUDE_4_5, reasoning: { modes: ['budget'] } },
  {
    provider: 'anthropic',
    id: 'claude-sonnet-4-5',
    family: true,
    limits: CLAUDE_4_5,
    reasoning: { modes: ['budget'] },
  },
  { provider: 'anthropic', id: 'claude-haiku-5', family: true, limits: CLAUDE_5, reasoning: CLAUDE_5_REASONING },
  { provider: 'anthropic', id: 'claude-sonnet-5', family: true, limits: CLAUDE_5, reasoning: CLAUDE_5_REASONING },
  { provider: 'anthropic', id: 'claude-opus-5', family: true, limits: CLAUDE_5, reasoning: CLAUDE_5_REASONING },
  {
    provider: 'anthropic',
    id: 'claude-haiku-5-5',
    label: 'Claude Haiku 5.5',
    note: 'Fast.',
    price: {
      asOf: '2026-10-08',
      source: 'anthropic.com',
      input: 0.1,
      cachedInput: 0.01,
      cacheWrite: 0.125,
      output: 0.5,
      longPrompt: { tokens: 100_000, factor: 5 },
    },
  },
  {
    provider: 'anthropic',
    id: 'claude-sonnet-5-5',
    label: 'Claude Sonnet 5.5',
    note: 'Careful; slower.',
    quirks: { appendsTools: true },
    // Sonnet 5's, taken as Sonnet 5.5's.
    price: { asOf: '2026-09-30', source: 'sentra.app', input: 2, cachedInput: 0.2, cacheWrite: 2.5, output: 10 },
  },
  {
    provider: 'anthropic',
    id: 'claude-opus-5-5',
    label: 'Claude Opus 5.5',
    note: 'Most capable; slowest.',
    price: { asOf: '2026-09-30', source: 'sentra.app', input: 4, cachedInput: 0.2, cacheWrite: 5, output: 20 },
  },

  // OpenAI (the Responses API; another endpoint knows none of these).
  {
    provider: 'openai',
    id: 'gpt-6',
    family: true,
    limits: { contextWindowTokens: 1_050_000, maxOutputTokens: 128_000 },
    reasoning: { modes: ['adaptive'], efforts: EFFORTS },
  },
  {
    provider: 'openai',
    id: 'gpt-6-luna',
    label: 'GPT-6 Luna',
    note: 'Fast.',
    quirks: { appendsTools: true },
    // Cached input assumed a tenth.
    price: { asOf: '2026-09-30', source: 'openrouter.ai', input: 0.1, cachedInput: 0.01, output: 0.5 },
  },
  {
    provider: 'openai',
    id: 'gpt-6.1-sol',
    label: 'GPT-6.1 Sol',
    note: 'Most capable GPT; slower.',
    price: { asOf: '2026-09-30', source: 'openai.com', input: 2, cachedInput: 0.1, cacheWrite: 2.5, output: 10 },
  },
];

const BOUNDARY = /^[-.:/_@]/;
// A dated snapshot of a model: -20251001, -2025-10-01.
const SNAPSHOT = /^-(?:\d{8}|\d{4}-\d{2}-\d{2})$/;

function matches(entry: CatalogEntry, id: string): boolean {
  if (id === entry.id) return true;
  if (!id.startsWith(entry.id)) return false;
  const rest = id.slice(entry.id.length);
  return entry.family ? BOUNDARY.test(rest) : SNAPSHOT.test(rest);
}

/** What the catalog knows of a model; undefined when no entry matches. `catalog` defaults to the built-in one. */
export function modelInfo(
  provider: string,
  id: string,
  catalog: readonly CatalogEntry[] = MODEL_CATALOG,
): ModelInfo | undefined {
  const matched = catalog
    .filter((e) => e.provider === provider && matches(e, id))
    // Families first, the shorter prefix first; the model last.
    .sort((a, b) => Number(!a.family) - Number(!b.family) || a.id.length - b.id.length);
  if (matched.length === 0) return undefined;
  let info: ModelInfo = { provider, id, listed: false, quirks: {} };
  for (const { provider: _, id: __, family, quirks, ...fields } of matched) {
    info = { ...info, ...fields, quirks: { ...info.quirks, ...quirks }, listed: info.listed || !family };
  }
  return info;
}

/** The listed models, for a host's model list (`provider` narrows it). */
export function listedModels(provider?: string, catalog: readonly CatalogEntry[] = MODEL_CATALOG): ModelInfo[] {
  return catalog
    .filter((e) => !e.family && (provider === undefined || e.provider === provider))
    .map((e) => modelInfo(e.provider, e.id, catalog))
    .filter((m): m is ModelInfo => m !== undefined);
}

/** The reasoning a host sends for a model when the person chose nothing (`AgentSessionOptions.thinking`). */
export function defaultReasoning(provider: string, id: string): ReasoningParameters | undefined {
  return modelInfo(provider, id)?.reasoning?.default;
}

/** One model call's tokens, as a price reads them: input split into uncached, cache reads and cache writes. */
export interface PricedUsage {
  readonly uncachedInput: number;
  readonly cacheRead: number;
  readonly cacheWrite: number;
  readonly output: number;
}

/** What one model call cost in dollars at a list price. */
export function dollarsOf(price: ModelPrice, usage: PricedUsage): number {
  const prompt = usage.uncachedInput + usage.cacheRead + usage.cacheWrite;
  const factor = price.longPrompt && prompt > price.longPrompt.tokens ? price.longPrompt.factor : 1;
  return (
    (factor *
      (usage.uncachedInput * price.input +
        usage.cacheRead * price.cachedInput +
        usage.cacheWrite * (price.cacheWrite ?? price.input) +
        usage.output * price.output)) /
    1e6
  );
}
