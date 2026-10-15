/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { urlCombine } from '@molstar/core/util/url';
import { PluginConfig } from '@molstar/plugin/config';
import { type ApiContext, ApiResult, memberName, type OperationDefinition } from '../../api/area.js';
import { messageOf } from '../../api/run.js';

// The `remote-states` feature's operations (operations.md §2, the `session` row: `remote.list`, `remote.open`,
// `remote.upload`, `remote.remove`; parity.md §1): Mol*'s state server (the classic State tab's Remote States,
// `RemoteStateSnapshots`) over the plugin's snapshot manager, as operations of the core's `session` area. The server is
// a parameter (empty: the viewer's `state-server` option, which is what the plugin's `CurrentServer` holds), so the
// views can point at another without a setting of the plugin's. Failures name the server in the hint.

const network = { modules: ['network'] as const, approval: 'network' as const };

/** An entry of the server's list, with the addresses that get and remove it. */
export interface RemoteStateEntry {
  readonly id: string;
  readonly name: string;
  readonly description: string;
  /** Milliseconds since the epoch, as the server saved it. */
  readonly timestamp: number;
  readonly isSticky: boolean;
  readonly url: string;
  readonly removeUrl: string;
}

/** What `session.remote.list` gives. */
export interface RemoteStateList {
  /** The server that answered. */
  readonly server: string;
  readonly entries: readonly RemoteStateEntry[];
}

const ServerParam = PD.Text('', {
  description: "The state server's address; empty: the viewer's state server",
});

export const ListParams = { server: ServerParam };

export const OpenParams = {
  entry: PD.Text('', { description: 'The state to open: its id, or its name' }),
  server: ServerParam,
};

export const UploadParams = {
  name: PD.Text('', { description: 'The name to save the current state under' }),
  description: PD.Text(''),
  playOnLoad: PD.Boolean(false, { description: 'Play the snapshots when the state is opened' }),
  server: ServerParam,
};

export const RemoveParams = {
  entry: PD.Text('', { description: 'The state to remove: its id, or its name' }),
  server: ServerParam,
};

function serverOf(ctx: ApiContext, given: string): string {
  const server = given.trim() || ctx.plugin.config.get(PluginConfig.State.CurrentServer) || '';
  return server.replace(/\/+$/, '');
}

const failure = (what: string, server: string, e: unknown) =>
  ApiResult.fail(
    'failed',
    `Could not ${what}: ${messageOf(e)}`,
    `Check the state server (${server}) and the network connection.`,
  );

/** The server's entries, the sticky ones first and then by time (as the classic list). */
async function fetchEntries(ctx: ApiContext, server: string): Promise<RemoteStateEntry[]> {
  const plugin = ctx.plugin;
  await plugin.initialized;
  const json = (await plugin.runTask(plugin.fetch({ url: urlCombine(server, 'list'), type: 'json' }))) ?? [];
  if (!Array.isArray(json)) throw new Error('the server did not answer with a list of states');
  return (json as Partial<RemoteStateEntry>[])
    .filter((e) => typeof e?.id === 'string' || typeof e?.id === 'number')
    .map((e) => ({
      id: String(e.id),
      name: e.name ?? '',
      description: e.description ?? '',
      timestamp: Number(e.timestamp) || 0,
      isSticky: !!e.isSticky,
      url: urlCombine(server, `get/${e.id}`),
      removeUrl: urlCombine(server, `remove/${e.id}`),
    }))
    .sort((a, b) => (a.isSticky === b.isSticky ? a.timestamp - b.timestamp : a.isSticky ? -1 : 1));
}

/** The id first, then the name; a name that several entries share is ambiguous. */
function pick(entries: readonly RemoteStateEntry[], key: string, server: string) {
  const wanted = key.trim();
  if (!wanted) return ApiResult.fail('invalid', 'Which state?', 'Give its id or its name.');
  const byId = entries.find((e) => e.id === wanted);
  if (byId) return byId;
  const named = entries.filter((e) => e.name.toLowerCase() === wanted.toLowerCase());
  if (named.length === 1) return named[0];
  if (named.length > 1) {
    return ApiResult.fail(
      'ambiguous',
      `Several states are named '${wanted}'.`,
      `Use an id: ${named.map((e) => e.id).join(', ')}.`,
    );
  }
  return ApiResult.fail(
    'not-found',
    `There is no state '${wanted}' on the server.`,
    `Server: ${server}. List the states with session.remote.list.`,
  );
}

const list: OperationDefinition<typeof ListParams, RemoteStateList> = {
  id: 'session.remote.list',
  title: 'List the remote states',
  params: ListParams,
  read: true,
  ...network,
  async run({ server: given }, ctx) {
    const server = serverOf(ctx, given);
    try {
      const entries = await fetchEntries(ctx, server);
      return ApiResult.ok(
        `${entries.length} state${entries.length === 1 ? '' : 's'} on the server.`,
        {},
        { value: { server, entries } },
      );
    } catch (e) {
      return failure('list the remote states', server, e);
    }
  },
};

const open: OperationDefinition<typeof OpenParams, void> = {
  id: 'session.remote.open',
  title: 'Open a remote state',
  params: OpenParams,
  ...network,
  // Destructive wins over network (api/area.ts).
  approval: 'destructive',
  async run({ entry: key, server: given }, ctx) {
    const server = serverOf(ctx, given);
    const plugin = ctx.plugin;
    let name: string;
    try {
      const found = pick(await fetchEntries(ctx, server), key, server);
      if (ApiResult.isFailure(found)) return found;
      name = found.name || found.id;
      // The server's entry holds the snapshot as `data` (`Snapshots.Fetch` of the plugin's state commands).
      const json = await plugin.runTask(plugin.fetch({ url: found.url, type: 'json' }));
      const snapshot = json?.data;
      if (!snapshot || !Array.isArray(snapshot.entries)) throw new Error(`'${name}' is not a saved state`);
      // A lazy extension whose transformers the state uses loads first, as `session.open` does.
      const transformers = [...JSON.stringify(snapshot).matchAll(/"transformer"\s*:\s*"([^"]+)"/g)].map((m) => m[1]);
      await ctx.extensions.ensureFor({ transformers });
      let failed: string | undefined;
      const sub = plugin.events.log.subscribe((e) => {
        if (e.type === 'error') failed ??= String(e.message);
      });
      try {
        await plugin.managers.snapshot.setStateSnapshot(snapshot);
      } finally {
        sub.unsubscribe();
      }
      if (failed) throw new Error(failed);
    } catch (e) {
      return failure(`open '${key}'`, server, e);
    }
    return ApiResult.ok(`Opened '${name}'.`, { view: true, sequence: true });
  },
};

const upload: OperationDefinition<typeof UploadParams, void> = {
  id: 'session.remote.upload',
  title: 'Upload the state',
  params: UploadParams,
  history: false,
  undoable: false,
  ...network,
  // It sends the state out (agent.md §4).
  approval: 'upload',
  when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
  async run({ name, description, playOnLoad, server: given }, ctx) {
    const server = serverOf(ctx, given);
    try {
      const body = JSON.stringify(
        // The person's Snapshots list is theirs: the scene goes into the upload only (keepList).
        await ctx.plugin.managers.snapshot.getStateSnapshot({ name, description, playOnLoad, keepList: true }),
      );
      const response = await fetch(
        urlCombine(server, `set?name=${encodeURIComponent(name)}&description=${encodeURIComponent(description)}`),
        {
          method: 'POST',
          mode: 'cors',
          referrer: 'no-referrer',
          headers: { 'Content-Type': 'application/json; charset=utf-8' },
          body,
        },
      );
      if (!response.ok) throw new Error(`the server answered ${response.status} ${response.statusText}`.trim());
    } catch (e) {
      return failure('upload the state', server, e);
    }
    return ApiResult.ok(`Uploaded the state${name ? ` '${name}'` : ''}.`);
  },
};

const remove: OperationDefinition<typeof RemoveParams, void> = {
  id: 'session.remote.remove',
  title: 'Remove a remote state',
  params: RemoveParams,
  undoable: false,
  ...network,
  // Destructive wins over network (api/area.ts).
  approval: 'destructive',
  async run({ entry: key, server: given }, ctx) {
    const server = serverOf(ctx, given);
    try {
      const found = pick(await fetchEntries(ctx, server), key, server);
      if (ApiResult.isFailure(found)) return found;
      if (found.isSticky) {
        return ApiResult.fail(
          'blocked',
          `'${found.name || found.id}' is a pinned state of the server.`,
          'The server keeps pinned states.',
        );
      }
      const response = await fetch(found.removeUrl);
      if (!response.ok) throw new Error(`the server answered ${response.status} ${response.statusText}`.trim());
      return ApiResult.ok(`Removed '${found.name || found.id}' from the server.`);
    } catch (e) {
      return failure(`remove '${key}'`, server, e);
    }
  },
};

export const remoteStatesOperations = (): { area: string; name: string; operation: OperationDefinition<any, any> }[] =>
  [list, open, upload, remove].map((operation) => ({ area: 'session', name: memberName(operation.id), operation }));
