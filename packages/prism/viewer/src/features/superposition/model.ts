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
import { combineLatest, debounceTime } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { SuperpositionEntry, SuperpositionModel } from '../../model/superposition.js';
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

const SelectionDelayMs = 80;

export interface SuperpositionDeps {
  readonly api: ViewerApi;
  readonly sources: ViewerSources;
  notify(text: string, kind?: MessageKind): void;
}

export interface Superposition extends SuperpositionModel {
  dispose(): void;
}

export function createSuperpositionModel(deps: SuperpositionDeps): Superposition {
  const form = new ParamsFormModel(FormParams);
  const available = observableState(false);
  const entries = observableState<readonly SuperpositionEntry[]>([]);
  const ready = observableState(false);
  const hint = observableState('');
  let disposed = false;
  let token = 0;
  let busy = false;

  const method = () => form.values.value.method;

  // The entries follow the method, the structures and the selection (short debounce: a pick fires several changes).
  const refresh = async () => {
    const mine = ++token;
    const structures = deps.sources.structures.value;
    const m = method();
    available.set(structures.length > 1);
    let list: SuperpositionEntry[];
    if (m === 'db') {
      list = structures.map((s) => ({ id: s.ref, label: `${s.label} · all chains` }));
    } else {
      const r = await deps.api.select.get();
      if (disposed || mine !== token) return;
      list =
        r.ok && r.value
          ? r.value.map((v) => ({ id: v.structure.ref, label: `${v.structure.label} · ${v.selected.line}` }))
          : [];
    }
    entries.set(list);
    ready.set(list.length >= 2);
    hint.set(list.length === 1 ? `${Hints[m]} One so far.` : Hints[m]);
  };

  const sub = combineLatest([form.values, deps.sources.structures, deps.sources.interaction.selection])
    .pipe(debounceTime(SelectionDelayMs))
    .subscribe(() => void refresh());
  void refresh();

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
    focus: (id) => void deps.api.focus.set({ target: 'selection', structure: id }),
    superpose: () => void superpose(),
    dispose() {
      disposed = true;
      token++;
      sub.unsubscribe();
    },
  };
}
