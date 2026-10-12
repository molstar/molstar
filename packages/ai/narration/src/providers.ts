/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { geminiNarration } from './gemini.js';
import { openaiNarration } from './openai.js';
import type { Fetch, NarrationProvider } from './provider.js';

/** The providers there are keys for (as the app names providers), in the order offered. */
export function narrationProviders(keys: { google?: string; openai?: string }, fetch?: Fetch): NarrationProvider[] {
  return [
    ...(keys.google ? [geminiNarration({ apiKey: keys.google, ...(fetch ? { fetch } : {}) })] : []),
    ...(keys.openai ? [openaiNarration({ apiKey: keys.openai, ...(fetch ? { fetch } : {}) })] : []),
  ];
}
