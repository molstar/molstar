/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { KeyBindings, KeymapStore } from '@molstar/prism-ui/core';
import type { PrismViewerSpec, PrismStore } from '../spec.js';

// The keymap's layers in the viewer (commands.md §6, viewer-spec.md §2.5): the app's keys are the spec's
// `ui.keybindings`; the person's remaps are kept with the preferences (`persist.preferences`, key `keybindings`), so
// they come back after a reload as the theme does. Without preferences they last for the session.

/** The preferences' key of the person's remapped keys (the Storage view names it). */
export const KeybindingsKey = 'keybindings';

/** The spec's `ui.keybindings` as the keymap's app layer (by command id; `null` removes a key). */
export function appKeybindings(spec: PrismViewerSpec): KeyBindings {
  return Object.fromEntries((spec.ui?.keybindings ?? []).map((k) => [k.command, k.key]));
}

/** The person's remaps in the preferences store; none without one. */
export function keymapStore(preferences: PrismStore | undefined): KeymapStore | undefined {
  if (!preferences) return undefined;
  return {
    load: () => preferences.get(KeybindingsKey) as KeyBindings | undefined,
    save: (bindings) => {
      if (Object.keys(bindings).length) preferences.set(KeybindingsKey, bindings);
      else preferences.remove(KeybindingsKey);
    },
  };
}
