/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The assistant, opt-in: an app gets it only by depending on this package (nothing in @molstar/prism-ui imports it).
// It plugs into the UI through its ports: tabs (Assistant, History in place of Log), the command menu's `assistant`,
// the canvas's agent slot (the turn card).

export * from './model.js';
export { AgentTurnCard, RunList, RunRow, UndoTurn } from './views/agent-turn.js';
export { AISetup } from './views/ai-setup.js';
export { AssistantTab } from './views/assistant-tab.js';
export { HistoryTab } from './views/history-tab.js';
