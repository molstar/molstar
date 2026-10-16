/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/prism-assistant/viewer`: what the viewer's spec loads lazily (`agent.assistant`, viewer-spec.md §2.4).
export { type AssistantOptions, createAssistant, createAssistantWith } from './viewer/ui.js';
export type { ProvidedModels, SettingsPrefs } from './viewer/settings-model.js';
// The agent comes with the assistant (the viewer's `agent.assistant` module, viewer-spec.md §2.4): the viewer does not
// import its code itself.
export { createViewerAgent } from '@molstar/prism-viewer/agent/session';
export { createViewerAgentModel, historyItemsOf } from './viewer/agent-model.js';
export { createAISettings, localPrefs } from './viewer/settings-model.js';
