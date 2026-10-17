/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  CommandDefinition,
  NavigationManager,
  ShellModel,
  StoredNavigation,
  StoredViewLayout,
  ViewStateManager,
} from '@molstar/prism-ui';
import { merge, skip } from 'rxjs';
import type { ViewerCommandContext } from '../api/managers.js';
import type { PrismStore, PrismViewerSpec } from '../spec.js';

// The layout kept in the preferences (journal.md §6, shell.md §3): the shell's layout (`viewState.save()`: docks, sizes,
// float, pin, hidden tabs) and what is open in the docks (`navigation.save()`, without the card and the page: they are
// summoned, not kept), restored at start and written a moment after each change. The host's parts of the layout (the
// mode, the compact rail, the buttons) stay the spec's: a kept layout never overrides what the host asked for.

/** The preferences' key the layout is kept under. */
export const LayoutKey = 'layout';

/** What is kept: both managers' own versioned forms; either of another version is dropped (journal.md §7). */
export interface KeptLayout {
  readonly view: StoredViewLayout;
  readonly navigation: StoredNavigation;
}

/** Ms after the last change before the layout is written (a drag of the splitter writes once). */
export const KeepLayoutDelayMs = 500;

/**
 * Where the layout is kept, if anywhere: the spec's `persist.preferences` unless `persist.layout` is off. By default an
 * embedded layout (`ui.layout.mode: 'embedded'`, or `compactRail`) is not kept: it is the host page's, not the person's.
 */
export function layoutStore(spec: PrismViewerSpec): PrismStore | undefined {
  const persist = spec.persist;
  if (!persist?.preferences) return undefined;
  const layout = spec.ui?.layout;
  const embedded = layout?.mode === 'embedded' || !!layout?.compactRail;
  return (persist.layout ?? !embedded) ? persist.preferences : undefined;
}

/**
 * Restores the kept layout into the managers, before the tabs register (a tab lands where it was left). A value of
 * another shape or version is ignored (false); tabs that are not registered stay in the layout.
 */
export function restoreLayout(
  managers: { readonly viewState: ViewStateManager; readonly navigation: NavigationManager },
  store: PrismStore,
): boolean {
  let kept: unknown;
  try {
    kept = store.get(LayoutKey);
  } catch {
    return false;
  }
  if (!kept || typeof kept !== 'object') return false;
  const { view, navigation } = kept as Partial<KeptLayout>;
  if (!view || typeof view !== 'object' || !view.layout || typeof view.layout !== 'object') return false;
  const host = managers.viewState.layout.value;
  const restored = managers.viewState.restore({
    ...view,
    layout: { ...view.layout, mode: host.mode, compactRail: host.compactRail, buttons: host.buttons },
  });
  if (!restored) return false;
  if (navigation && typeof navigation === 'object' && navigation.navigation) {
    managers.navigation.restore({
      ...navigation,
      navigation: { docks: navigation.navigation.docks, card: undefined, page: undefined },
    });
  }
  return true;
}

/** What the managers keep now. */
export function keptLayout(managers: {
  readonly viewState: ViewStateManager;
  readonly navigation: NavigationManager;
}): KeptLayout {
  const { version, navigation } = managers.navigation.save();
  return {
    view: managers.viewState.save(),
    navigation: { version, navigation: { docks: navigation.docks, card: undefined, page: undefined } },
  };
}

/**
 * Writes the layout to the store a moment after each change (`delay`), only when it differs from what was written;
 * once more when the page is hidden (a reload) and when stopped. Returns the stop.
 */
export function keepLayout(
  managers: { readonly viewState: ViewStateManager; readonly navigation: NavigationManager },
  store: PrismStore,
  delay = KeepLayoutDelayMs,
): () => void {
  let written = JSON.stringify(keptLayout(managers));
  let timer: ReturnType<typeof setTimeout> | undefined;
  const write = () => {
    if (timer !== undefined) clearTimeout(timer);
    timer = undefined;
    const value = keptLayout(managers);
    const text = JSON.stringify(value);
    if (text === written) return;
    written = text;
    try {
      store.set(LayoutKey, value);
    } catch {
      // Not kept: a store that fails breaks nothing (journal.md §6).
    }
  };
  const schedule = () => {
    if (timer !== undefined) clearTimeout(timer);
    timer = setTimeout(write, delay);
  };
  // Their current values are `written`'s: only what follows is a change.
  const changes = merge(managers.viewState.layout.pipe(skip(1)), managers.navigation.state.pipe(skip(1))).subscribe(
    schedule,
  );
  const onHide = () => {
    if (timer !== undefined) write();
  };
  const win = typeof window !== 'undefined' ? window : undefined;
  win?.addEventListener('pagehide', onHide);
  return () => {
    win?.removeEventListener('pagehide', onHide);
    changes.unsubscribe();
    if (timer !== undefined) write();
  };
}

/** The layout's commands: reset it, show the hidden tabs (shell.md §3). */
export function layoutCommands(shell: ShellModel): CommandDefinition<ViewerCommandContext>[] {
  return [
    {
      id: 'view.reset-layout',
      title: 'Reset layout',
      alias: ['reset layout', 'reset panels', 'default layout'],
      category: 'View',
      icon: 'rotate-ccw',
      run: () => {
        shell.resetLayout();
        return { message: 'The layout is back to its defaults.' };
      },
    },
    {
      id: 'view.show-tabs',
      title: 'Show hidden tabs',
      alias: ['show tabs', 'unhide tabs', 'hidden tabs'],
      category: 'View',
      icon: 'eye',
      when: () => (shell.hiddenTabs.value.length ? true : 'No tab is hidden'),
      run: () => {
        const titles = shell.hiddenTabs.value.map((t) => t.title);
        shell.showAllTabs();
        return { message: `Shown: ${titles.join(', ')}.` };
      },
    },
  ];
}
