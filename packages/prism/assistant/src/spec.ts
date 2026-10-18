/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PrismViewerAgentSpec } from '@molstar/prism-viewer/spec';
import type { AssistantOptions, ProvidedModels } from './viewer.js';

// `@molstar/prism-assistant/spec`: what a page's spec writes for the assistant (viewer-spec.md §2.4, agent.md §1). Types
// only and one dynamic import: the assistant, the agent and their dependencies (the AI SDK, zod) stay out of the page's
// startup bundle and load with the assistant's chunk.

/** The assistant's options as a spec gives them. */
export interface AssistantSpecOptions extends Omit<AssistantOptions, 'provided'> {
  /**
   * Models the host brings with their providers; as a function, made when the assistant loads, so the provider's code is
   * not in the page's startup bundle either.
   */
  readonly provided?: readonly ProvidedModels[] | (() => Promise<readonly ProvidedModels[]>);
}

/**
 * The spec's `agent.assistant`: `agent: { assistant: assistant({ provider: 'google', model: 'gemini-3.8-flash' }) }`.
 * Returns the lazy loader; nothing of the assistant is imported until the viewer calls it (the Assistant tab opened, the
 * Ask row used).
 */
export function assistant(options: AssistantSpecOptions = {}): NonNullable<PrismViewerAgentSpec['assistant']> {
  return async () => {
    const m = await import('./viewer.js');
    const { provided, ...rest } = options;
    const models = typeof provided === 'function' ? await provided() : provided;
    return {
      createViewerAgent: m.createViewerAgent,
      createAssistant: m.createAssistantWith({ ...rest, ...(models ? { provided: models } : {}) }),
    };
  };
}
