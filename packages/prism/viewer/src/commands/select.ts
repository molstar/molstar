/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition, CommandGroup } from '@molstar/prism-ui/core';
import type { Observable } from 'rxjs';
import type { ApiFailure, ViewerApi } from '../api/area.js';
import type { ViewerCommandContext } from '../api/managers.js';
import type { ResolvedTarget, TargetResolver } from '../api/targets.js';
import type { LociPreview, QueryEntry, SelectQuerySource } from '../sources/select-queries.js';
import { SelectCommandGroups } from './select-groups.js';

// The `select` group (commands.md §5, api.md §9 generated commands): one command per selection query, each running
// `select.set` with that query as its target and the group's shared mode. By query: the plugin's selection-query
// registry by category (Type, Structure property, …; an extension's queries too); in the loaded structures: their
// entities, chains, ligands and non-standard residues, and elements, made again when the structures change
// (`sources/select-queries.ts`). A highlighted row lights its atoms in 3D (the viewer's preview sink,
// `lociPreviewSink`) and says how much it would select.

type Ctx = ViewerCommandContext;

export { SelectCommandGroups, SelectGroupId } from './select-groups.js';

/** "412 residues · 3,301 atoms" (one residue, or atoms alone when the counts are the same). */
export function countLine(residues: number, atoms: number): string {
  if (!atoms) return 'nothing';
  const n = (x: number) => x.toLocaleString('en-US');
  const r = `${n(residues)} residue${residues === 1 ? '' : 's'}`;
  return residues === atoms ? `${n(atoms)} atom${atoms === 1 ? '' : 's'}` : `${r} · ${n(atoms)} atoms`;
}

export interface SelectCommandsOptions {
  readonly queries: SelectQuerySource;
  /** The API as the person. */
  readonly api: ViewerApi;
  readonly targets: TargetResolver;
  /** Emits when the loaded structures change (`sources.structures`). */
  readonly structures: Observable<unknown>;
  /** Registers commands (the viewer's `registerCommands`: overrides applied); returns the removal. */
  readonly register: (defs: readonly CommandDefinition<Ctx>[]) => () => void;
  readonly registerGroups: (groups: readonly CommandGroup[]) => () => void;
}

/**
 * The `select` group's commands, kept current: made again when the loaded structures change (the registry's are read
 * then too, so an extension's queries come in with the next change). Returns the removal of all of it.
 */
export function watchSelectCommands(o: SelectCommandsOptions): () => void {
  const { api, queries } = o;
  const dropGroups = o.registerGroups(SelectCommandGroups);
  let drop: (() => void) | undefined;
  let dropDynamicGroups: (() => void) | undefined;
  let made = '';
  /** Resolved targets of this structure set, by command id: the preview and the count resolve once. */
  let cache = new Map<string, Promise<ResolvedTarget | ApiFailure>>();

  const ensure = (entry: QueryEntry) => queries.ensure(entry);
  const resolve = (entry: QueryEntry) => {
    let r = cache.get(entry.id);
    if (!r) {
      r = ensure(entry).then(() => o.targets.resolve({ molql: entry.expression }));
      cache.set(entry.id, r);
    }
    return r;
  };
  const command = (entry: QueryEntry): CommandDefinition<Ctx> => ({
    id: entry.id,
    title: entry.title,
    category: 'Select',
    group: entry.group,
    icon: entry.icon,
    ...(entry.alias?.length ? { alias: entry.alias } : {}),
    when: (ctx) =>
      ctx.empty ? 'Nothing is loaded.' : (api.describe('select.set')?.available ?? 'Selecting is not available.'),
    preview: () =>
      ({
        kind: 'loci',
        resolve: async () => {
          const r = await resolve(entry);
          return r.ok ? r.parts.map((p) => p.loci) : [];
        },
      }) satisfies LociPreview,
    count: async (_args, _ctx, signal) => {
      const r = await resolve(entry);
      if (signal.aborted || !r.ok) return undefined;
      let residues = 0;
      let atoms = 0;
      for (const p of r.parts) {
        residues += p.selected.residues;
        atoms += p.selected.atoms;
      }
      return countLine(residues, atoms);
    },
    async run(args) {
      await ensure(entry);
      const mode = typeof args.mode === 'string' ? args.mode : 'set';
      const r = await api.run('select.set', { target: { molql: entry.expression }, mode });
      if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
      // The query's MolQL is long: the operation is said with the query's name.
      return { message: r.did, ops: `select.set({ target: ${JSON.stringify(entry.title)} (query), mode: "${mode}" })` };
    },
  });

  const update = () => {
    cache = new Map();
    const registry = queries.list();
    const entries = registry.entries;
    // Made again only when what is listed changed (a structure's own change keeps the same rows).
    const key = `${registry.groups.map((g) => g.id).join()}|${entries.map((e) => `${e.id}:${e.title}`).join()}`;
    if (key === made) return;
    made = key;
    drop?.();
    dropDynamicGroups?.();
    dropDynamicGroups = o.registerGroups(registry.groups);
    drop = o.register(entries.map(command));
  };
  const sub = o.structures.subscribe(() => update());
  return () => {
    sub.unsubscribe();
    drop?.();
    dropDynamicGroups?.();
    dropGroups();
    cache.clear();
  };
}
