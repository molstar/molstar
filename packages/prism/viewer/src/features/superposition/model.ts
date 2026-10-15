/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { MessageKind } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { combineLatest, type Subscription } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { SuperpositionEntry, SuperpositionModel } from '../../model/superposition.js';
import type { PickedStructure } from '../../sources/selection.js';
import type { ViewerSources } from '../../sources/viewer-sources.js';
import type { SuperposeMethod } from './operations.js';

// The view model of the Superpose section (state.md: "the algorithms move out of the view into the API"): the form, what
// the method would use (the selection of each structure, or every structure for UniProt) and Superpose, which calls
// `superpose.*` as the person. Over the people's API and the sources only (no plugin). DOM-free.

/** The method and its options, as the classic panel offers them. */
export const FormParams = {
  method: PD.Select<SuperposeMethod>(
    'chains',
    [
      ['chains', 'Chains'],
      ['tm-align', 'TM-align'],
      ['db', 'UniProt (SIFTS)'],
      ['ligands', 'Ligands'],
      ['atoms', 'Atoms'],
    ],
    { isEssential: true, description: 'What the structures are superposed by.' },
  ),
  alignSequences: PD.Boolean(true, {
    label: 'Align sequences',
    isEssential: true,
    description: 'Guide the superposition by a sequence alignment of the chains.',
    hideIf: (p: { method: SuperposeMethod }) => p.method !== 'chains',
  }),
  traceOnly: PD.Boolean(true, {
    label: 'Trace only',
    isEssential: true,
    description: 'Use only CA (and equivalent) atoms.',
    hideIf: (p: { method: SuperposeMethod }) => p.method !== 'chains' && p.method !== 'tm-align' && p.method !== 'db',
  }),
};

const Hints: Record<SuperposeMethod, string> = {
  chains: 'Select one chain (or residues within it) in each of two or more structures.',
  'tm-align': 'Select one chain in each of two or more structures; best when sequence identity is low.',
  ligands: 'Select one ligand in each of two or more structures.',
  atoms: 'Select the same atoms, in the same order, in each structure.',
  db: 'Superposes every structure by the residues they share in UniProt.',
};

export interface SuperpositionDeps {
  readonly api: ViewerApi;
  readonly sources: ViewerSources;
  notify(text: string, kind?: MessageKind): void;
}

export interface Superposition extends SuperpositionModel {
  /**
   * The section shows (true) or is folded away (false): the entries follow the structures and the selection only
   * while it shows (prism AGENTS.md, hidden UI). Not active, the model holds no subscription.
   */
  setActive(active: boolean): void;
  dispose(): void;
}

/** The fewest atoms that fix an orientation (as the operations). */
const MinAtoms = 3;

const count = (n: number, what: string) => `${n} ${what}${n === 1 ? '' : 's'}`;

/** Why a method cannot use this selection, or undefined. Mirrors what `superpose.*` checks of each target. */
export function problemOf(
  method: SuperposeMethod,
  p: PickedStructure,
  reference?: PickedStructure,
): string | undefined {
  switch (method) {
    case 'chains':
    case 'tm-align':
      if (p.polymerAtoms === 0) return 'no polymer in it: select a chain';
      if (p.polymerChains > 1) return `its polymer is in ${count(p.polymerChains, 'chain')}: select one`;
      return undefined;
    case 'ligands':
      return p.ligand ? undefined : 'not one ligand residue';
    case 'atoms':
      if (p.atoms < MinAtoms) return `${count(p.atoms, 'atom')}: ${MinAtoms} are needed`;
      if (reference && p !== reference && p.atoms !== reference.atoms) {
        return `${count(p.atoms, 'atom')}, the reference has ${reference.atoms}`;
      }
      return undefined;
    default:
      return undefined;
  }
}

export function createSuperpositionModel(deps: SuperpositionDeps): Superposition {
  const form = new ParamsFormModel(FormParams);
  const available = observableState(false);
  const entries = observableState<readonly SuperpositionEntry[]>([]);
  const ready = observableState(false);
  const hint = observableState('');
  const blocked = observableState<string | undefined>(undefined);
  let sub: Subscription | undefined;
  let busy = false;

  const method = () => form.values.value.method;

  // The entries follow the method, the structures and the selection: what each structure's selection holds comes from
  // the interaction source (no runs of the API for a read), flagged where the method would refuse it.
  const refresh = () => {
    const structures = deps.sources.structures.value;
    const m = method();
    available.set(structures.length > 1);
    let list: SuperpositionEntry[];
    if (m === 'db') {
      list = structures.map((s) => ({ id: s.ref, label: `${s.label} · all chains` }));
    } else {
      const picked = deps.sources.interaction.picked.value;
      list = picked.map((p) => {
        const problem = problemOf(m, p, picked[0]);
        return { id: p.ref, label: `${p.label} · ${p.line}`, ...(problem ? { problem } : {}) };
      });
    }
    const flagged = list.find((e) => e.problem);
    entries.set(list);
    blocked.set(flagged ? `${flagged.label.split(' · ')[0]}: ${flagged.problem}.` : undefined);
    ready.set(list.length >= 2 && !flagged);
    hint.set(list.length === 1 ? `${Hints[m]} One so far.` : Hints[m]);
  };
  refresh();

  const setActive = (active: boolean) => {
    if (active === !!sub) return;
    if (!active) {
      sub?.unsubscribe();
      sub = undefined;
      return;
    }
    // The first emission of each (the current value) is the refresh itself.
    sub = combineLatest([form.values, deps.sources.structures, deps.sources.interaction.picked]).subscribe(refresh);
  };

  const superpose = async () => {
    if (busy || !ready.value) return;
    busy = true;
    try {
      const v = form.values.value;
      const ids = entries.value.map((e) => e.id);
      const api = deps.api.superpose;
      if (!api) return;
      const scoped = ids.map((structure) => ({ target: 'selection', structure }));
      const r =
        v.method === 'chains'
          ? await api.chains({ targets: scoped, alignSequences: v.alignSequences, traceOnly: v.traceOnly })
          : v.method === 'tm-align'
            ? await api.tmAlign({ targets: scoped, traceOnly: v.traceOnly })
            : v.method === 'ligands'
              ? await api.ligands({ targets: scoped })
              : v.method === 'atoms'
                ? await api.atoms({ targets: scoped })
                : await api.db({ structures: ids, traceOnly: v.traceOnly });
      // A failure says itself where the person looks (the notification), with what to do.
      deps.notify(r.ok ? r.did : r.hint ? `${r.message} ${r.hint}` : r.message, r.ok ? 'success' : 'error');
    } finally {
      busy = false;
    }
  };

  return {
    available: available.asReadonly,
    form,
    entries: entries.asReadonly,
    ready: ready.asReadonly,
    hint: hint.asReadonly,
    blocked: blocked.asReadonly,
    focus: (id) => void deps.api.focus.set({ target: 'selection', structure: id }),
    superpose: () => void superpose(),
    setActive,
    dispose() {
      setActive(false);
    },
  };
}
