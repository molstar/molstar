/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StorageEntry, StorageSource } from '@molstar/prism-viewer/storage/model';
import type { KeyScope, KeyStore } from '../keys.js';
import type { AISettingsModel } from '../model.js';

// The assistant's storage source (journal.md §7): the stored keys, by provider and scope, never the key. Forget removes
// one through the settings model, so its setup forgets it too. Domain rules are not listed: the assistant keeps none
// yet (the agent's "Allow this domain" lasts for the page, agent/policies.ts).

const Scopes: Readonly<Record<KeyScope, string>> = {
  browser: 'Kept in this browser (localStorage, scrambled)',
  session: 'Kept for this tab (sessionStorage, scrambled)',
  reload: 'In memory until the page reloads',
};

export function keysStorageSource(
  keys: Pick<KeyStore, 'list' | 'remove'>,
  settings?: Pick<AISettingsModel, 'providers' | 'removeKey'>,
): StorageSource {
  const labelOf = (provider: string) => settings?.providers.value.find((p) => p.id === provider)?.label ?? provider;
  return {
    id: 'assistant-keys',
    title: 'Assistant keys',
    where: 'localStorage and sessionStorage',
    async list() {
      return (await keys.list()).map(
        ({ provider, scope }): StorageEntry => ({ id: provider, label: labelOf(provider), detail: Scopes[scope] }),
      );
    },
    async forget(provider) {
      if (settings) settings.removeKey(provider);
      else await keys.remove(provider);
    },
  };
}
