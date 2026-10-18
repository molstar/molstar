/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it, vi } from 'vitest';
import { assistant } from '../spec.js';
import { preferredModel } from '../viewer/ui.js';

// `assistant(options)`, the spec's `agent.assistant` (viewer-spec.md §2.4): a lazy loader, the model a host prefers.

describe('the assistant in a spec', () => {
  it('imports only types: the assistant loads when the viewer calls the loader', () => {
    const source = readFileSync(new URL('../spec.ts', import.meta.url), 'utf8');
    const imports = source.match(/^import\b.*$/gm) ?? [];
    expect(imports.length).toBeGreaterThan(0);
    for (const line of imports) expect(line).toMatch(/^import type /);
  });

  it('gives the viewer the assistant and the agent, models made when it loads', async () => {
    const provided = vi.fn(async () => []);
    const load = assistant({ provided, provider: 'google' });
    expect(provided).not.toHaveBeenCalled();
    const m = await load();
    expect(provided).toHaveBeenCalledTimes(1);
    expect(typeof m.createAssistant).toBe('function');
    expect(typeof m.createViewerAgent).toBe('function');
  });

  it("fills in the host's model: the provider's first listed, the catalog's provider, the fallback", () => {
    expect(preferredModel({})).toEqual({ provider: 'anthropic', model: 'claude-sonnet-5-5' });
    expect(preferredModel({ provider: 'anthropic' })).toEqual({ provider: 'anthropic', model: 'claude-sonnet-5-5' });
    expect(preferredModel({ model: 'gemini-3.8-flash' })).toEqual({ provider: 'google', model: 'gemini-3.8-flash' });
    expect(preferredModel({ provider: 'google' })).toEqual({ provider: 'google', model: 'gemini-3.8-flash' });
    expect(
      preferredModel({
        provider: 'scripted',
        provided: [{ id: 'scripted', label: 'S', models: ['s1'], provider: {} as never }],
      }),
    ).toEqual({ provider: 'scripted', model: 's1' });
  });
});
