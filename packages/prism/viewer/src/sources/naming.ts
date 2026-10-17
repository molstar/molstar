/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
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

/**
 * The structure the person is working with (targets.md §1, "the structure"): the focused one, else the one holding the
 * selection when only one does. Undefined when neither says; then an operation takes the only one loaded, or (a target)
 * every one, or fails `ambiguous`.
 */
export function workingStructure(
  plugin: PluginContext,
  list: readonly StructureRef[] = plugin.managers.structure.hierarchy.current.structures,
): { readonly structure: StructureRef; readonly why: 'focused' | 'selected' } | undefined {
  const focus = plugin.managers.structure.focus.current?.loci;
  if (focus) {
    const root = focus.structure.root;
    const s = list.find((x) => x.cell.obj?.data === focus.structure || x.cell.obj?.data.root === root);
    if (s) return { structure: s, why: 'focused' };
  }
  const selected = structuresWithSelection(plugin, list);
  return selected.length === 1 ? { structure: selected[0], why: 'selected' } : undefined;
}

/** The structures that hold some of the selection. */
export function structuresWithSelection(plugin: PluginContext, list: readonly StructureRef[]): StructureRef[] {
  return list.filter((s) => {
    const data = s.cell.obj?.data;
    if (!data) return false;
    const loci = plugin.managers.structure.selection.getLoci(data);
    return StructureElement.Loci.is(loci) && !StructureElement.Loci.isEmpty(loci);
  });
}

/**
 * The loaded structures' names by ref, made once while the structures and the state stay the same (the hierarchy's list
 * and the state tree are new objects at every change; a name changes only with one). Made per call, a result naming
 * every structure (a preset over all) was structures × structures.
 */
const namesCache = new WeakMap<object, { list: unknown; tree: unknown; names: ReadonlyMap<string, string> }>();
export function loadedNames(plugin: PluginContext): ReadonlyMap<string, string> {
  const list = plugin.managers.structure.hierarchy.current.structures;
  const tree = plugin.state.data.tree;
  const hit = namesCache.get(plugin);
  if (hit && hit.list === list && hit.tree === tree) return hit.names;
  const names = structureNames(list);
  namesCache.set(plugin, { list, tree, names });
  return names;
}
