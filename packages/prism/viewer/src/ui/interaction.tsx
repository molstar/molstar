/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mapState } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type { SelectionModel } from '../model/selection.js';
import type { InteractionSources } from '../sources/interaction.js';
import type { CanvasMode, Granularity, ModeSwitchProps } from '../views/canvas/mode-switch.js';

// The canvas's interaction chrome over the viewer (plan 2.3): the mode switch and the picking level (Mol*'s selection
// mode and granularity, `select.mode` and `select.granularity`), and the selection bar (what is selected, and what to do
// with it: `focus.set`, `component.add`, `select.clear`). Every action is an operation, as the person.

/** Mol*'s granularities, labelled and drawn as picking levels. */
export const Granularities: readonly Granularity[] = [
  { value: 'element', label: 'Atom', icon: 'pick-atom' },
  { value: 'residue', label: 'Residue', icon: 'pick-residue' },
  { value: 'chain', label: 'Chain', icon: 'pick-chain' },
  { value: 'entity', label: 'Entity', icon: 'pick-entity' },
  { value: 'model', label: 'Model', icon: 'pick-model' },
  { value: 'operator', label: 'Operator', icon: 'copy' },
  { value: 'structure', label: 'Structure', icon: 'structure' },
  { value: 'elementInstances', label: 'Atom', icon: 'pick-atom', group: 'All instances' },
  { value: 'residueInstances', label: 'Residue', icon: 'pick-residue', group: 'All instances' },
  { value: 'chainInstances', label: 'Chain', icon: 'pick-chain', group: 'All instances' },
];

export function modeSwitchProps(interaction: InteractionSources, api: ViewerApi): ModeSwitchProps {
  return {
    mode: mapState(interaction.selectionMode, (on): CanvasMode => (on ? 'select' : 'interact'), Object.is),
    setMode: (m) => void api.select.mode({ on: m === 'select' }),
    granularity: interaction.granularity,
    granularities: Granularities,
    setGranularity: (g) => void api.run('select.granularity', { granularity: g }),
  };
}

export function selectionModel(interaction: InteractionSources, api: ViewerApi): SelectionModel {
  return {
    summary: interaction.selection,
    actions: [
      { id: 'focus', label: 'Focus', icon: 'focus', run: () => void api.focus.set({ target: 'selection' }) },
      {
        id: 'component',
        label: 'Component',
        icon: 'layers',
        run: () => void api.component.add({ target: 'selection', label: 'Selection' }, { copy: true }),
      },
    ],
    clear: () => void api.select.clear(),
  };
}

export { ModeSwitch, PickLevel } from '../views/canvas/mode-switch.js';
export { SelectionBar } from '../views/canvas/selection-bar.js';
