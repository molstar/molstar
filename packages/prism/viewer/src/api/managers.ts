/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  CommandManager,
  HistoryManager,
  NavigationManager,
  NotificationsManager,
  UiTransactionCoordinator,
  ViewStateManager,
} from '@molstar/prism-ui/core';
import type { ViewerCheckpoint } from '../sources/checkpoint.js';

// The viewer's managers (api.md §11): generic and DOM-free, from prism-ui, composed by the viewer as
// `PluginContext` has its managers.

/** What commands see: whether anything is loaded, the viewer's API bound to the person. */
export interface ViewerCommandContext {
  readonly empty: boolean;
}

export interface ViewerManagers {
  readonly commands: CommandManager<ViewerCommandContext>;
  /** Optimistic UI transactions: an edit shows at once and is absorbed when the source reports it, or rolled back. */
  readonly transactions: UiTransactionCoordinator;
  readonly history: HistoryManager<ViewerCheckpoint>;
  readonly navigation: NavigationManager;
  readonly viewState: ViewStateManager;
  readonly notifications: NotificationsManager;
}
