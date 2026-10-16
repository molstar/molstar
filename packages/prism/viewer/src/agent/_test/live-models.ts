/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import type { Provider } from '@molstar/ai';
import { createAnthropicProvider, createGoogleProvider } from '@molstar/ai/providers';

// The models live runs use (`*.live.ts`, `pnpm test:live`): Haiku and Gemini Flash. The keys come from the repository's
// `.env.local` and are never printed; a run skips a model without its key.

function keys(): Record<string, string> {
  const file = path.resolve(import.meta.dirname, '../../../../../../.env.local');
  if (!fs.existsSync(file)) return {};
  const out: Record<string, string> = {};
  for (const line of fs.readFileSync(file, 'utf8').split('\n')) {
    const m = /^\s*([A-Z0-9_]+)\s*=\s*"?([^"\n]*)"?\s*$/.exec(line);
    if (m) out[m[1]] = m[2];
  }
  return out;
}

export const env = keys();

export interface LiveModel {
  readonly provider: string;
  readonly model: string;
  readonly key: string;
  readonly create: () => Provider;
}

export const LiveModels: readonly LiveModel[] = [
  {
    provider: 'anthropic',
    model: 'claude-haiku-5-5',
    key: 'ANTHROPIC_API_KEY',
    create: () => createAnthropicProvider({ resolveApiKey: () => env.ANTHROPIC_API_KEY }),
  },
  {
    provider: 'google',
    model: 'gemini-3.8-flash',
    key: 'GEMINI_API_KEY',
    create: () => createGoogleProvider({ resolveApiKey: () => env.GEMINI_API_KEY }),
  },
];
