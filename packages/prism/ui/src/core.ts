/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/prism-ui/core`: the DOM-free part of the library (spec library.md §1), what a headless app imports: reactive
// state, commands, forms, the managers, the shell's model, the canvas ink. The package's index adds the DOM runtime,
// components and views.

export * from './reactive/combine-behaviors.js';
export * from './reactive/disposable-group.js';
export * from './reactive/keyed-selector.js';
export * from './reactive/keyed-state.js';
export * from './reactive/work-queue.js';
export * from './reactive/map-state.js';
export * from './reactive/observable-state.js';
export * from './reactive/readonly-state.js';
export * from './reactive/transaction.js';

export * from './model/canvas-theme.js';

export * from './managers/history.js';
export * from './managers/navigation.js';
export * from './managers/notifications.js';
export {
  DefaultDockSizes,
  type DockLayout,
  defaultViewLayout,
  MinDockSizes,
  type StoredViewLayout,
  type ViewLayout,
  type ViewLayoutSeed,
  ViewStateManager,
} from './managers/view-state.js';

export * from './shell/shell.js';
export * from './shell/theme.js';

export * from './commands/items.js';
export * from './commands/keymap.js';
export * from './commands/keys.js';
export * from './commands/manager.js';
export * from './commands/types.js';

// Types only (erased): what the DOM-free contracts name.
export type { IconName } from './components/icon/icon.js';
export type { ViewRenderer } from './dom/view.js';
