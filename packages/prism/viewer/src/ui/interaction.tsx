/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mapState, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerApi } from '../api/area.js';
import type { SelectionAction, SelectionModel, SelectionPart, SelectionSummaryModel } from '../model/selection.js';
import type { InteractionSources } from '../sources/interaction.js';
import type { PickedStructure } from '../sources/selection.js';
import { sameItems } from '../sources/state.js';
import type { CanvasMode, Granularity, ModeSwitchProps } from '../views/canvas/mode-switch.js';

// The canvas's interaction chrome over the viewer: the mode switch and the picking level (Mol*'s selection
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

/**
 * The selection bar's menu: Focus, Show in structure tree (`showInTree`: the Structure panel with its tree filtered by
 * the selection), Make a component, the extensions' rows (`extra`, by order), and the Selection tab for the rest
 * (`openPanel`).
 */
export function selectionModel(
  interaction: InteractionSources,
  api: ViewerApi,
  o: {
    openPanel?: () => void;
    showInTree?: () => void;
    extra?: ReadonlyState<readonly SelectionAction[]>;
  } = {},
): SelectionModel {
  const own: SelectionAction[] = [
    {
      id: 'focus',
      label: 'Focus',
      hint: 'Bring the selection into view',
      icon: 'focus',
      order: 0,
      run: () => void api.focus.set({ target: 'selection' }),
    },
    ...(o.showInTree
      ? [
          {
            id: 'show-in-tree',
            label: 'Show in structure tree',
            hint: 'Only the rows that hold it',
            icon: 'list-tree' as const,
            order: 1,
            run: o.showInTree,
          },
        ]
      : []),
    // Without the structure feature there is no component area to add to (api.md §6).
    ...(api.area('component')
      ? [
          {
            id: 'component',
            label: 'Make a component',
            hint: 'A copy of it, as its own component',
            icon: 'layers' as const,
            order: 2,
            run: () => void api.component.add({ target: 'selection', label: 'Selection' }, { copy: true }),
          },
        ]
      : []),
    ...(o.openPanel
      ? [
          {
            id: 'more',
            label: 'More in the Selection panel…',
            hint: 'Color, measure, superpose',
            icon: 'mouse-pointer-click' as const,
            order: 1000,
            run: o.openPanel,
          },
        ]
      : []),
  ];
  const byOrder = (a: SelectionAction, b: SelectionAction) => (a.order ?? 500) - (b.order ?? 500);
  const actions = o.extra
    ? derive(o.extra, (extra) => [...own, ...extra].sort(byOrder))
    : observableState<readonly SelectionAction[]>(own).asReadonly;
  return {
    summary: interaction.selection,
    actions,
    clear: () => void api.select.clear(),
    focused: mapState(interaction.focus, (f) => f !== undefined, Object.is),
    unfocus: () => void api.focus.clear(),
  };
}

/** Where in a structure the selection is: "chain A · 60 atoms"; with several structures, "chains A, B · 8 residues · 120 atoms". */
function partDetail(p: PickedStructure, several: boolean): string {
  const chains =
    p.chains.length === 0
      ? ''
      : p.chains.length === 1
        ? `chain ${p.chains[0]}`
        : p.chains.length <= 4
          ? `chains ${p.chains.join(', ')}`
          : `${p.chains.length} chains`;
  // One structure: the summary above counts the residues already.
  const residues = several && p.residues > 1 ? `${p.residues.toLocaleString('en-US')} residues` : '';
  const atoms = `${p.atoms.toLocaleString('en-US')} atom${p.atoms === 1 ? '' : 's'}`;
  return [chains, residues, atoms].filter(Boolean).join(' · ');
}

const samePart = (a: SelectionPart, b: SelectionPart) => a.id === b.id && a.label === b.label && a.detail === b.detail;

/**
 * The selection at the top of the Selection tab: the summary, a line per structure it is in (the interaction source's
 * `picked`, followed only while the tab shows), Focus (`focus.set`) and Clear (`select.clear`).
 */
export function selectionSummaryModel(interaction: InteractionSources, api: ViewerApi): SelectionSummaryModel {
  return {
    summary: interaction.selection,
    parts: mapState(
      interaction.picked,
      (list) => list.map((p) => ({ id: p.ref, label: p.label, detail: partDetail(p, list.length > 1) })),
      (a, b) => sameItems(a, b, samePart),
    ),
    selectMode: interaction.selectionMode,
    pickLevel: mapState(interaction.granularity, (g) => pickLevelLabel(g), Object.is),
    startSelecting: () => void api.select.mode({ on: true }),
    focus: () => void api.focus.set({ target: 'selection' }),
    clear: () => void api.select.clear(),
  };
}

/** A granularity as the picking level menu says it ("Residue", "Chain (all instances)"). */
function pickLevelLabel(value: string): string {
  const g = Granularities.find((x) => x.value === value);
  if (!g) return value;
  return g.group ? `${g.label} (${g.group.toLowerCase()})` : g.label;
}

export { ModeSwitch, modeMenuRows, PickLevel } from '../views/canvas/mode-switch.js';
export { SelectionBar } from '../views/canvas/selection-bar.js';
