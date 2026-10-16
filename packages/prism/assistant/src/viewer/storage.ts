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
import type { DomainRules } from '@molstar/prism-viewer/agent/policies';

// The assistant's storage sources (journal.md §7): the stored keys, by provider and scope, never the key (Forget removes
// one through the settings model, so its setup forgets it too); and the domains the person allowed the agent to reach
// (agent.md §4), for the session or kept in this browser, each with Forget.

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

/** The domains kept in this browser, from the prefs' text (a JSON list); anything else: none. */
export function keptDomains(text: string | undefined): string[] {
  try {
    const v = text ? JSON.parse(text) : [];
    return Array.isArray(v) ? v.filter((d): d is string => typeof d === 'string' && /^[a-z0-9.-]+$/i.test(d)) : [];
  } catch {
    return [];
  }
}

/** The allowed domains in the Storage view: each with where it lasts, and Forget. */
export function domainsStorageSource(rules: Pick<DomainRules, 'list' | 'remove'>): StorageSource {
  return {
    id: 'assistant-domains',
    title: 'Allowed domains',
    where: 'localStorage',
    async list() {
      return rules.list().map(
        ({ domain, scope }): StorageEntry => ({
          id: domain,
          label: domain,
          detail: scope === 'browser' ? 'The assistant reaches it without asking, in this browser' : 'For this session',
        }),
      );
    },
    async forget(domain) {
      rules.remove(domain);
    },
  };
}
