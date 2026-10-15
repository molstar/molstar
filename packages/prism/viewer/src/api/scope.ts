/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { findStructure, structureNames } from '../sources/naming.js';
import { type ApiContext, type ApiFailure, ApiResult } from './area.js';

// Which structures an operation acts on (targets.md §1, "the structure"): a name (a state ref, a label, an entry id),
// else the one the person works with: the focused one, else the only one loaded; operations that act on many take
// every loaded one. Shared by the scene areas.

export function loadedStructures(ctx: ApiContext): readonly StructureRef[] {
  return ctx.plugin.managers.structure.hierarchy.current.structures;
}

/** The structure's name for people (`sources/naming.ts`): its entry, told apart when two share one. */
export function labelOf(ctx: ApiContext, s: StructureRef): string {
  return structureNames(loadedStructures(ctx)).get(s.cell.transform.ref) ?? s.cell.transform.ref;
}

function byName(ctx: ApiContext, name: string): StructureRef | undefined {
  return findStructure(loadedStructures(ctx), name);
}

function notFound(ctx: ApiContext, name: string): ApiFailure {
  const list = loadedStructures(ctx);
  return ApiResult.fail(
    'not-found',
    `There is no structure '${name}'.`,
    list.length ? `Loaded: ${list.map((x) => labelOf(ctx, x)).join(', ')}.` : 'Nothing is loaded.',
  );
}

/** One structure: the named one, else the focused one, else the only one loaded. */
export function oneStructure(ctx: ApiContext, name: string | undefined): StructureRef | ApiFailure {
  const list = loadedStructures(ctx);
  if (!list.length) return ApiResult.fail('not-found', 'Nothing is loaded.', 'Load a structure first.');
  if (name) return byName(ctx, name) ?? notFound(ctx, name);
  if (list.length === 1) return list[0];
  const focus = ctx.plugin.managers.structure.focus.current?.loci.structure.root;
  const focused = focus && list.find((s) => s.cell.obj?.data.root === focus);
  if (focused) return focused;
  return ApiResult.fail(
    'ambiguous',
    `Which structure? ${list.length} are loaded.`,
    `Name one: ${list.map((x) => labelOf(ctx, x)).join(', ')}.`,
  );
}

/** The named structure, or every loaded one. */
export function someStructures(ctx: ApiContext, name: string | undefined): readonly StructureRef[] | ApiFailure {
  const list = loadedStructures(ctx);
  if (!list.length) return ApiResult.fail('not-found', 'Nothing is loaded.', 'Load a structure first.');
  if (!name) return list;
  const s = byName(ctx, name);
  return s ? [s] : notFound(ctx, name);
}

export function isFailure(x: unknown): x is ApiFailure {
  return !!x && typeof x === 'object' && (x as { ok?: unknown }).ok === false;
}
