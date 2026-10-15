/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The agent's Inspector (Trace, Tools, Limits), a page over the whole app: its own entry
// (`@molstar/prism-assistant/inspector`), kept out of the package's index so the Assistant's tab and card do not pull it
// in. The ways in (the page id, the commands) are in the index (`inspector-access`).

export { InspectorPageId, inspectorCommands, type InspectorFocus } from './inspector-access.js';
export { LimitsView, parseLimit } from './inspector/limits.js';
export { inspectLinks } from './inspector/links.js';
export * from './inspector/model.js';
export { InspectorBody, InspectorHeader, inspectorPage } from './inspector/page.js';
export { BudgetMeter, Facts, Json, Meter, TextTabs, Usage } from './inspector/parts.js';
export { fieldsOf, type SchemaField, ToolsView } from './inspector/tools.js';
export { TraceView } from './inspector/trace.js';
