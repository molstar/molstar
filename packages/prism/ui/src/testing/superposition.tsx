/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { combineLatest, map, type Observable } from 'rxjs';
import { domView } from '../dom/view.js';
import { ParamsFormModel } from '../forms/params-form.js';
import type { SuperpositionEntry, SuperpositionModel } from '../model/superposition.js';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import type { ShellModel } from '../shell/shell.js';
import { Superpose } from '../views/structure/superpose.js';
import type { MockScene } from './scene.js';

const chainBased = (p: any) => p.method !== 'chains' && p.method !== 'uniprot';

// Mol*'s superposition methods and options (StructureSuperpositionParams).
const Params = {
  method: PD.Select(
    'chains',
    [
      ['chains', 'Chains'],
      ['tm-align', 'TM-align'],
      ['uniprot', 'UniProt (SIFTS)'],
      ['ligands', 'Ligands'],
      ['atoms', 'Atoms'],
    ],
    { isEssential: true, description: 'What the structures are superposed by.' },
  ),
  alignSequences: PD.Boolean(true, {
    label: 'Align sequences',
    isEssential: true,
    description: 'Guide the superposition by a sequence alignment of the chains.',
    hideIf: (p: any) => p.method !== 'chains',
  }),
  traceOnly: PD.Boolean(true, {
    label: 'Trace only',
    isEssential: true,
    description: 'Use only CA (and equivalent) atoms.',
    hideIf: chainBased,
  }),
};

const Hints: Record<string, string> = {
  chains: 'Select one chain (or residues within it) in each of two or more structures.',
  'tm-align': 'Select one chain in each of two or more structures; best when sequence identity is low.',
  ligands: 'Select one ligand in each of two or more structures.',
  atoms: 'Select the same atoms, in the same order, in each structure.',
};

function follow<T>(source: Observable<T>, initial: T): ReadonlyState<T> {
  const state = observableState(initial);
  source.subscribe((v) => state.set(v));
  return state.asReadonly;
}

/** Superposition for the stories: the CYP3A4 structures, with chains selected in two of them. */
export function mockSuperposition(scene: MockScene, post: (text: string) => void): SuperpositionModel {
  const form = new ParamsFormModel(Params);
  const structures = scene.roots.pipe(map((r) => r.filter((n) => n.kind === 'structure')));
  const method = form.values.pipe(map((v: any) => v.method as string));
  const selections: Record<string, SuperpositionEntry[]> = {
    chains: [
      { id: 'c1', label: '1TQN · Chain A' },
      { id: 'c2', label: '1W0E · Chain A' },
    ],
    'tm-align': [
      { id: 't1', label: '1TQN · Chain A' },
      { id: 't2', label: '1W0E · Chain A' },
    ],
    ligands: [{ id: 'l1', label: '1TQN · HEM 508' }],
    atoms: [],
  };
  const entries = combineLatest([method, structures]).pipe(
    map(([m, s]): readonly SuperpositionEntry[] =>
      m === 'uniprot' ? s.map((n) => ({ id: n.id, label: `${n.label.value} · all chains` })) : (selections[m] ?? []),
    ),
  );
  const ready = combineLatest([method, entries]).pipe(map(([m, e]) => (m === 'atoms' ? e.length >= 1 : e.length >= 2)));
  const hint = combineLatest([method, entries]).pipe(
    map(([m, e]) => (e.length === 1 ? `${Hints[m] ?? ''} One so far.` : (Hints[m] ?? ''))),
  );
  const count = structures.pipe(map((s) => s.length));
  return {
    available: follow(count.pipe(map((c) => c > 1)), false),
    form,
    entries: follow(entries, []),
    ready: follow(ready, false),
    hint: follow(hint, ''),
    focus: (id) =>
      post(
        `Focused on ${
          Object.values(selections)
            .flat()
            .find((e) => e.id === id)?.label ?? id
        }`,
      ),
    superpose: () => post(`Superposed by ${(form.values.value as any).method} · RMSD 0.38 Å`),
  };
}

/** Superposing is done with the selection: a Superpose section of the Selection tab, with two or more structures. */
export function registerMockSuperposition(shell: ShellModel, model: SuperpositionModel) {
  shell.registerSections([
    {
      id: 'superpose',
      slot: 'selection',
      title: 'Superpose',
      icon: 'superpose',
      keywords: ['superposition', 'align', 'tm-align', 'rmsd', 'uniprot', 'sifts'],
      order: 3,
      when: () => model.available.value,
      render: domView(() => <Superpose model={model} />),
    },
  ]);
}
