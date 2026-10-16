/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type CatalogEntry, defaultReasoning, dollarsOf, listedModels, MODEL_CATALOG, modelInfo } from '../models.js';

describe('the model catalog', () => {
  it('knows a listed model from its entry over its family', () => {
    const sonnet = modelInfo('anthropic', 'claude-sonnet-5-5');
    expect(sonnet).toMatchObject({
      listed: true,
      label: 'Claude Sonnet 5.5',
      limits: { contextWindowTokens: 1_000_000, maxOutputTokens: 128_000 },
      quirks: { appendsTools: true },
    });
    expect(sonnet?.price?.asOf).toMatch(/^\d{4}-\d{2}-\d{2}$/);
  });

  it('knows a new version of a family by the longest prefix, without a label or a price', () => {
    const next = modelInfo('google', 'gemini-3.9-flash');
    expect(next).toMatchObject({ listed: false, limits: { maxOutputTokens: 65_536 } });
    expect(next?.label).toBeUndefined();
    expect(next?.price).toBeUndefined();
    // A newer Sonnet 5 is a Sonnet 5, without 5.5's tool additions.
    expect(modelInfo('anthropic', 'claude-sonnet-5-7')?.quirks.appendsTools).toBeUndefined();
  });

  it('matches at a boundary, and a model as itself and its dated snapshots', () => {
    expect(modelInfo('openai', 'gpt-6.1-sol')?.listed).toBe(true);
    expect(modelInfo('openai', 'gpt-60')).toBeUndefined();
    expect(modelInfo('anthropic', 'claude-haiku-4-5-20251001')?.limits?.contextWindowTokens).toBe(200_000);
    // A variant is not the model: Luna's price and tool additions are not Luna Mini's.
    const mini = modelInfo('openai', 'gpt-6-luna-mini');
    expect(mini).toMatchObject({ listed: false, quirks: {} });
    expect(mini?.price).toBeUndefined();
    expect(modelInfo('google', 'claude-sonnet-5-5')).toBeUndefined();
  });

  it('lists every model entry, each with a label, a note and limits', () => {
    const listed = listedModels();
    expect(listed.length).toBe(MODEL_CATALOG.filter((e) => !e.family).length);
    for (const m of listed) {
      expect(m.label).toBeTruthy();
      expect(m.note).toBeTruthy();
      expect(m.limits).toBeDefined();
    }
    expect(listedModels('openai').map((m) => m.id)).toEqual(['gpt-6-luna', 'gpt-6.1-sol']);
  });

  it('gives the reasoning default: low on Gemini 3, none elsewhere', () => {
    expect(defaultReasoning('google', 'gemini-3.8-flash')).toEqual({ mode: 'adaptive', effort: 'low' });
    expect(defaultReasoning('google', 'gemini-2.5-flash')).toBeUndefined();
    expect(defaultReasoning('anthropic', 'claude-sonnet-5-5')).toBeUndefined();
  });

  it('takes a catalog of the host’s own', () => {
    const own: CatalogEntry[] = [
      {
        provider: 'einfra',
        id: 'glm-5',
        family: true,
        limits: { contextWindowTokens: 1_048_576, maxOutputTokens: 131_072 },
      },
      { provider: 'einfra', id: 'glm-5.3', label: 'GLM 5.3', note: 'Open.' },
    ];
    expect(modelInfo('einfra', 'glm-5.3', own)).toMatchObject({ listed: true, label: 'GLM 5.3' });
    expect(modelInfo('einfra', 'glm-5.3')).toBeUndefined();
  });

  it('prices a call, cache reads and writes apart, and a long prompt at its factor', () => {
    const usage = { uncachedInput: 1_000_000, cacheRead: 1_000_000, cacheWrite: 1_000_000, output: 1_000_000 };
    expect(dollarsOf({ asOf: '2026-01-01', source: 'x', input: 1, cachedInput: 0.1, output: 5 }, usage)).toBeCloseTo(
      7.1,
    );
    const tiered = {
      asOf: '2026-01-01',
      source: 'x',
      input: 1,
      cachedInput: 0.1,
      cacheWrite: 1.25,
      output: 5,
      longPrompt: { tokens: 100_000, factor: 5 },
    };
    expect(dollarsOf(tiered, { uncachedInput: 50_000, cacheRead: 0, cacheWrite: 0, output: 0 })).toBeCloseTo(0.05);
    expect(dollarsOf(tiered, { uncachedInput: 50_000, cacheRead: 60_000, cacheWrite: 0, output: 0 })).toBeCloseTo(0.28);
  });
});
