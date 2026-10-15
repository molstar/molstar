/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MessageKind } from '@molstar/prism-ui/core';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ApiResult, ViewerApi } from '../../api/area.js';
import type { RemoteState, RemoteStatesModel } from '../../model/home.js';
import type { RemoteStateEntry } from './operations.js';

// The view model of the Remote states card (model/home.ts `RemoteStatesModel`), over the `session.remote.*` operations
// as the person: it holds the list and the server's address (view state of the session: no operation sets the
// plugin's server) and acts through the API only. DOM-free.

export interface RemoteStatesDeps {
  readonly api: ViewerApi;
  notify(text: string, kind?: MessageKind): void;
}

/** "Sun, 11 Oct 2026 09:05" (UTC), as the card shows dates. */
export const formatDate = (timestamp: number) =>
  timestamp ? new Date(timestamp).toUTCString().slice(0, 22) : undefined;

const stateOf = (e: RemoteStateEntry): RemoteState => ({
  id: e.id,
  label: e.name || formatDate(e.timestamp) || e.id,
  date: e.name ? formatDate(e.timestamp) : undefined,
});

export function createRemoteStatesModel(deps: RemoteStatesDeps): RemoteStatesModel {
  const serverUrl = observableState('');
  const states = observableState<readonly RemoteState[]>([]);
  const loading = observableState(false);
  let entries: readonly RemoteStateEntry[] = [];
  // The latest refresh wins: a slow answer from a server the person has since left is dropped.
  let token = 0;

  /** One operation at a time shows as loading; a failure is posted, not thrown. */
  const run = async (work: () => Promise<unknown>) => {
    loading.set(true);
    try {
      await work();
    } finally {
      loading.set(false);
    }
  };
  const server = () => serverUrl.value.trim();
  const say = (r: ApiResult) =>
    deps.notify(r.ok ? r.did : r.hint ? `${r.message} ${r.hint}` : r.message, r.ok ? 'success' : 'error');

  const refresh = async () => {
    const mine = ++token;
    loading.set(true);
    try {
      const r = await deps.api.session.remoteList?.({ server: server() });
      if (mine !== token || !r) return;
      if (!r.ok) {
        entries = [];
        states.set([]);
        say(r);
        return;
      }
      entries = r.value?.entries ?? [];
      states.set(entries.map(stateOf));
      // The server that answered: what the viewer's option said when the field was empty.
      if (!server() && r.value) serverUrl.set(r.value.server);
    } finally {
      if (mine === token) loading.set(false);
    }
  };

  const idOf = (id: string) => entries.find((e) => e.id === id)?.id ?? id;

  return {
    serverUrl: serverUrl.asReadonly,
    setServerUrl(url) {
      serverUrl.set(url.trim());
      void refresh();
    },
    states: states.asReadonly,
    loading: loading.asReadonly,
    refresh: () => void refresh(),
    open(id) {
      void run(async () => {
        const r = await deps.api.session.remoteOpen!({ entry: idOf(id), server: server() });
        say(r);
        return r;
      });
    },
    upload(name) {
      void run(async () => {
        const r = await deps.api.session.remoteUpload!({ name, server: server() });
        say(r);
        if (r.ok) await refresh();
        return r;
      });
    },
    remove(id) {
      void run(async () => {
        const r = await deps.api.session.remoteRemove!({ entry: idOf(id), server: server() });
        say(r);
        if (r.ok) await refresh();
        return r;
      });
    },
  };
}
