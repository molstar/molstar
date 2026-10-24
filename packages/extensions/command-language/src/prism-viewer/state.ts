/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type { ApiAreaContext } from '@molstar/prism-viewer/api/area';
import { labelOf, loadedStructures } from '@molstar/prism-viewer/api/scope';
import { DefaultVmdRep, type ObjectInfo, type ScriptEnv } from '../language/env.js';
import type { RepType, Sel, VmdRep } from '../language/ir.js';

// The session's state in the viewer (_spec §5.1): owner-keyed custom state, so Undo's checkpoints and saved sessions
// keep it with the tree. Object names and VMD molecule ids by structure ref, the named selections, VMD's top molecule
// and rep defaults, and which components the console made for which PyMOL representation (so `show` grows the same
// component). Tcl variables and atomselect handles are the console's, not this.

export const StateKey = 'command-language';

export interface CommandLanguageState {
  readonly version: 1;
  /** Object name and molecule id by structure ref. */
  readonly objects: Readonly<Record<string, { readonly name: string; readonly molecule: number }>>;
  readonly nextMolecule: number;
  /** VMD's top molecule, by structure ref. */
  readonly top?: string;
  readonly names: Readonly<Record<string, Sel>>;
  readonly repDefaults: VmdRep;
  /** The components the console made for PyMOL's representations: component ref to its type. */
  readonly display: Readonly<Record<string, RepType>>;
}

export const EmptyState: CommandLanguageState = {
  version: 1,
  objects: {},
  nextMolecule: 0,
  names: {},
  repDefaults: DefaultVmdRep,
  display: {},
};

/** The state's holder: one per viewer, registered with the plugin's custom state. */
export class StateHolder {
  value: CommandLanguageState = EmptyState;
  readonly #unregister: () => void;

  constructor(plugin: PluginContext) {
    this.#unregister = plugin.state.registerCustomState(StateKey, {
      get: () => (this.value === EmptyState ? undefined : this.value),
      set: (v) => {
        const s = v as CommandLanguageState | undefined;
        this.value = s?.version === 1 ? s : EmptyState;
      },
    });
  }

  update(patch: Partial<CommandLanguageState>): void {
    this.value = { ...this.value, ...patch };
  }

  dispose(): void {
    this.#unregister();
  }
}

const sanitize = (s: string) => s.replace(/[^A-Za-z0-9_]/g, '_').replace(/^_+|_+$/g, '') || 'obj';

/**
 * Each loaded structure as an object: the name and molecule id the state gives it, else one made now (its entry, in
 * load order; ids after the state's). The same for a preview and a run while the scene stays the same.
 */
export function objectsOf(
  ctx: Pick<ApiAreaContext, 'plugin'>,
  state: CommandLanguageState,
): (ObjectInfo & { readonly ref: string })[] {
  const out: (ObjectInfo & { ref: string })[] = [];
  const taken = new Set<string>();
  let next = state.nextMolecule;
  const list = loadedStructures(ctx);
  // Known ones first, so a made name never takes one the state has.
  for (const s of list) {
    const known = state.objects[s.cell.transform.ref];
    if (known) taken.add(known.name.toLowerCase());
  }
  for (const s of list) {
    const ref = s.cell.transform.ref;
    const known = state.objects[ref];
    if (known) {
      out.push({ ref, name: known.name, molecule: known.molecule, ...contentsOf(s) });
      continue;
    }
    let name = sanitize(labelOf(ctx, s).toLowerCase());
    for (let i = 2; taken.has(name.toLowerCase()); i++) name = `${sanitize(labelOf(ctx, s).toLowerCase())}_${i}`;
    taken.add(name.toLowerCase());
    out.push({ ref, name, molecule: next++, ...contentsOf(s) });
  }
  return out;
}

/** For completion: the structure's author chains and its ligands' residue names (a bounded walk). */
function contentsOf(s: StructureRef): { chains: string[]; ligands: string[] } {
  const data = s.cell.obj?.data;
  if (!data) return { chains: [], ligands: [] };
  const chains = new Set<string>();
  const ligands = new Set<string>();
  for (const unit of data.units) {
    if (!Unit.isAtomic(unit)) continue;
    const { chains: ch, atoms, chainAtomSegments } = unit.model.atomicHierarchy;
    const first = unit.elements[0];
    chains.add(ch.auth_asym_id.value(chainAtomSegments.index[first]));
    if (chains.size > 64) break;
    const entities = unit.model.entities;
    const eIdx = entities.getEntityIndex(ch.label_entity_id.value(chainAtomSegments.index[first]));
    if (entities.data.type.value(eIdx) === 'non-polymer' && ligands.size < 64) {
      ligands.add(atoms.label_comp_id.value(first));
    }
  }
  return { chains: [...chains].sort(), ligands: [...ligands].sort() };
}

/** What names mean now (`ScriptHost.env`). */
export function envOf(
  ctx: Pick<ApiAreaContext, 'plugin'>,
  state: CommandLanguageState,
  files: readonly string[] = [],
): ScriptEnv & { readonly objects: readonly (ObjectInfo & { readonly ref: string })[] } {
  const objects = objectsOf(ctx, state);
  const top = state.top ? objects.find((o) => o.ref === state.top)?.name : undefined;
  return {
    objects,
    names: new Map(Object.entries(state.names)),
    ...(top ? { top } : {}),
    repDefaults: state.repDefaults,
    nextMolecule: Math.max(state.nextMolecule, ...objects.map((o) => o.molecule + 1)),
    files,
  };
}
