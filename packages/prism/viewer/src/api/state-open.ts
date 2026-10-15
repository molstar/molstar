/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { urlCombine } from '@molstar/core/util/url';
import { PluginConfig } from '@molstar/plugin/config';
import type { ApiContext } from './area.js';

// Opening a state of Mol*'s state server by address, shared by `load.snapshot` (by id, the loaders feature) and
// `session.remote.open` (the remote-states feature); features do not import each other, so the one implementation
// lives here. It throws what went wrong: the callers say it in their own words and hint.

/** The state server the viewer points at (its `state-server` option), without a trailing slash. */
export function currentStateServer(ctx: ApiContext): string {
  return (ctx.plugin.config.get(PluginConfig.State.CurrentServer) ?? '').replace(/\/+$/, '');
}

/** The transformers a state's JSON names: a lazy extension that owns one loads before the state is applied. */
export function transformersOf(text: string): string[] {
  return [...text.matchAll(/"transformer"\s*:\s*"([^"]+)"/g)].map((m) => m[1]);
}

/** The address of a state on a server: `<server>/get/<id>`. */
export const stateAddress = (server: string, id: string) => urlCombine(server, `get/${encodeURIComponent(id)}`);

/**
 * Fetches a state from its address (the server's entry holds the snapshot as `data`, `Snapshots.Fetch` of the plugin's
 * state commands) and applies it: lazy extensions owning its transformers load first, as `session.open` does, and an
 * error the plugin logs while applying (it reports and carries on) is thrown. `label` names it in the message when it
 * is not a state.
 */
export async function openStateFrom(ctx: ApiContext, url: string, label: string): Promise<void> {
  const plugin = ctx.plugin;
  await plugin.initialized;
  const json = await plugin.runTask(plugin.fetch({ url, type: 'json' }));
  const snapshot = json?.data;
  if (!snapshot || !Array.isArray(snapshot.entries)) throw new Error(`'${label}' is not a saved state`);
  await ctx.extensions.ensureFor({ transformers: transformersOf(JSON.stringify(snapshot)) });
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
}
