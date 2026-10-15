/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';

// What a structure is called for people, the model and the API (targets.md §1, "a label"): its entry ("1CBS"), the
// plugin's label of it added only when another loaded structure has the same entry ("1CBS · Assembly 1"). The plugin
// labels structures by kind ("Assembly 1", "Model"), which does not tell them apart.

export function entryOf(s: StructureRef): string {
  const model = s.model?.cell.obj?.data;
  return model?.entryId || s.model?.cell.obj?.label || s.cell.obj?.label || s.cell.transform.ref;
}

/** Each structure's name, by its state ref. */
export function structureNames(list: readonly StructureRef[]): ReadonlyMap<string, string> {
  const counts = new Map<string, number>();
  for (const s of list) counts.set(entryOf(s), (counts.get(entryOf(s)) ?? 0) + 1);
  const names = new Map<string, string>();
  for (const s of list) {
    const entry = entryOf(s);
    const own = s.cell.obj?.label;
    names.set(s.cell.transform.ref, (counts.get(entry) ?? 0) > 1 && own ? `${entry} · ${own}` : entry);
  }
  return names;
}

/** Finds a structure by state ref, name, entry, or the plugin's label (case-insensitive). */
export function findStructure(list: readonly StructureRef[], name: string): StructureRef | undefined {
  const lower = name.trim().toLowerCase();
  const names = structureNames(list);
  return (
    list.find((s) => s.cell.transform.ref === name) ??
    list.find((s) => names.get(s.cell.transform.ref)?.toLowerCase() === lower) ??
    list.find((s) => entryOf(s).toLowerCase() === lower) ??
    list.find((s) => (s.cell.obj?.label ?? '').toLowerCase() === lower)
  );
}
