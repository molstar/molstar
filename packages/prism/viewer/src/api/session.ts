/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ApiContext } from './area.js';
import { ApiArea, ApiResult } from './area.js';
import { asText, contentOf, extensionOf, InputParams, inputModuleProblem, inputOf, isZip } from './input.js';
import { OptionalParam } from './params.js';
import { messageOf } from './run.js';
import { labelOf, loadedStructures } from './scope.js';
import { currentStateServer, openStateFrom, stateAddress, transformersOf } from './state-open.js';

// `session` (operations.md §2): `open` (`.molx`, `.molj`, from a url, data, a file, or a state server's id) and `clear`
// now; `save` with `export-session`, the remote operations with `remote-states` (plan phase 3). Both are history steps:
// Undo brings the scene back. The loaders' `load.snapshot` is `open` under another name (`openSession`), for the app's
// `snapshot-url` and `snapshot-id` URL options.

/** A session file's record (journal.md §5): `record.jsonl` in a session's zip, a `record` field in a state. */
function recordLinesOf(extra: { files?: Record<string, Uint8Array>; fields?: Record<string, unknown> } | undefined) {
  const file = extra?.files?.['record.jsonl'];
  if (file) return new TextDecoder().decode(file);
  const field = extra?.fields?.record;
  return Array.isArray(field) ? field.map((e) => JSON.stringify(e)).join('\n') : '';
}

/** The nodes an opened session could not make: how many, the first one's error, and whether every top one failed. */
function failedNodes(ctx: ApiContext): { count: number; first: string; all: boolean } {
  const state = ctx.plugin.state.data;
  let count = 0;
  let first = '';
  state.cells.forEach((cell) => {
    if (cell.status !== 'error') return;
    count++;
    first ||= cell.errorText ?? 'an error';
  });
  const top = state.tree.children.get(StateTransform.RootRef);
  const all = count > 0 && top.size > 0 && top.toArray().every((ref) => state.cells.get(ref)?.status === 'error');
  return { count, first, all };
}

export type SessionType = 'auto' | 'molx' | 'molj';

/** `session.open`'s params (and `load.snapshot`'s): one of url, data, file (operations.md §1, "Inputs") or id. */
export const SessionOpenParams = {
  ...InputParams({
    url: 'The address of a session (.molx) or a state (.molj)',
    data: "A state's JSON (.molj), or a session's bytes (.molx; 'base64,…' as text)",
    file: 'A session (.molx) or a state (.molj) file',
    accept: '.molx,.molj,.json',
  }),
  id: OptionalParam(
    PD.Text('', { description: "A state on the state server (the viewer's state-server option), by id" }),
  ),
  type: PD.Select<SessionType>(
    'auto',
    [
      ['auto', 'From the name, else the content'],
      ['molx', 'Session (.molx)'],
      ['molj', 'State (.molj)'],
    ],
    { description: "The kind: by default the name's extension (.molx, .molj), else the content (a zip is a session)" },
  ),
};

/** What an opened session holds: the names of its structures. */
export interface Opened {
  readonly loaded: readonly string[];
}

/** The kind a name's or an address's extension says (`.molx`; `.molj` or `.json`); a query or fragment is ignored. */
export function sessionTypeOf(name: string): 'molj' | 'molx' | undefined {
  const ext = extensionOf(name);
  return ext === 'molx' ? 'molx' : ext === 'molj' || ext === 'json' ? 'molj' : undefined;
}

/** The address a call reaches (the agent's policy, agent.md §4); an id is the configured state server's. */
export const sessionAddresses = ({ url }: { url?: string }): string[] => (url?.trim() ? [url.trim()] : []);

/**
 * Opens a session or a state (`session.open`, and `load.snapshot`, which is the same): from a url, data or a file
 * (`.molx` or `.molj`, by `type`, else the name, else the content), or a state server's state by id. It replaces the
 * scene. Only the open itself failing is fatal; a node of the session that cannot be made (a server out of reach) is a
 * warning beside what did open, and an error logged by other work changes nothing.
 */
export async function openSession(
  values: PD.ValuesFor<typeof SessionOpenParams>,
  ctx: ApiContext,
): Promise<ApiResult<Opened>> {
  const id = values.id?.trim() ?? '';
  const input = inputOf(values, { verb: 'Open', or: { name: 'id', given: !!id } });
  if (ApiResult.isFailure(input)) return input;
  const unavailable = inputModuleProblem(input ?? { kind: 'id' }, ctx, 'Open a session');
  if (unavailable) return unavailable;
  let name: string;
  let extra: Awaited<ReturnType<typeof ctx.plugin.managers.snapshot.open>> | undefined;
  if (!input) {
    const server = currentStateServer(ctx);
    if (!server) {
      return ApiResult.fail(
        'unavailable',
        `This viewer has no state server to open state ${id} from.`,
        "Set its state-server option, or give the state's url, data or file.",
      );
    }
    name = `state ${id}`;
    try {
      await openStateFrom(ctx, stateAddress(server, id), id);
    } catch (e) {
      return ApiResult.fail(
        'failed',
        `Could not open ${name}: ${messageOf(e)}`,
        'Check the id and the network connection.',
      );
    }
  } else {
    let content: string | Uint8Array;
    try {
      content = await contentOf(ctx, input);
    } catch (e) {
      return ApiResult.fail(
        'failed',
        `Could not load ${input.kind === 'url' ? input.url : input.name}: ${messageOf(e)}`,
        'Check the address and the network connection.',
      );
    }
    const kind =
      values.type !== 'auto' ? values.type : (sessionTypeOf(input.name) ?? (isZip(content) ? 'molx' : 'molj'));
    name = input.name || (kind === 'molx' ? 'the session' : 'the state');
    // A lazy extension whose transformers the state uses loads first (a .molj names them; a .molx is zipped).
    if (kind === 'molj') await ctx.extensions.ensureFor({ transformers: transformersOf(asText(content)) });
    // Mol*'s snapshot manager reads the kind from the file's name.
    const base = input.name.replace(/\.[^.]*$/, '') || 'session';
    const file = new File([content as BlobPart], `${base}.${kind}`);
    let opened = false;
    let logged: string | undefined;
    const subs = [
      ctx.plugin.managers.snapshot.events.opened.subscribe(() => {
        opened = true;
      }),
      ctx.plugin.events.log.subscribe((e) => {
        if (e.type === 'error') logged = String(e.message);
      }),
    ];
    try {
      await ctx.plugin.initialized;
      extra = await ctx.plugin.managers.snapshot.open(file);
    } finally {
      for (const s of subs) s.unsubscribe();
    }
    if (!opened) return ApiResult.fail('failed', `Could not open ${name}: ${logged ?? 'it is not a session'}`);
  }
  const failed = failedNodes(ctx);
  if (failed.all) {
    return ApiResult.fail(
      'failed',
      `Could not open ${name}: nothing in it could be made (${failed.first}).`,
      'A state names its data; data opened from local files is only in a session (.molx).',
    );
  }
  // A record saved with it (journal.md §5): shown in History, read-only, beside the opened scene.
  const lines = recordLinesOf(extra);
  const shown = !!lines && ctx.journal.show(name, lines);
  const names = loadedStructures(ctx).map((s) => labelOf(ctx, s));
  const warnings = failed.count
    ? [`${failed.count} part${failed.count === 1 ? '' : 's'} of it could not be made: ${failed.first}`]
    : undefined;
  return ApiResult.ok(
    `Opened ${name}.${shown ? ' Its record is in History, read-only.' : ''}`,
    { view: true, sequence: true },
    { value: { loaded: names }, ...(warnings ? { warnings } : {}) },
  );
}

export class SessionApi extends ApiArea {
  readonly id = 'session';

  readonly open = this.operation({
    id: 'session.open',
    title: 'Open a session',
    params: SessionOpenParams,
    // The modules are the input's, checked per call: a url and an id need the network, a file access to files.
    // It replaces the scene: destructive, which wins over network (api/area.ts).
    approval: 'destructive',
    addresses: sessionAddresses,
    run: openSession,
  });

  readonly clear = this.operation({
    id: 'session.clear',
    title: 'Clear the session',
    params: {},
    approval: 'destructive',
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run(_, ctx) {
      const tree = ctx.plugin.state.data.tree;
      const removed = tree.children.get(StateTransform.RootRef).size;
      await ctx.plugin.state.data.build().delete(StateTransform.RootRef).commit();
      // A cleared scene is small: per-step undo is back on (sequence.md §7.4).
      ctx.managers.history.setCheckpointing(true);
      return ApiResult.ok(`Cleared the session (${removed} item${removed === 1 ? '' : 's'}).`, { view: true });
    },
  });
}
