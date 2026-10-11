/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { domView } from '../../dom/view.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { SelectionModel } from '../../model/selection.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel } from '../../shell/shell.js';

// Mol*'s selection theme actions (the selection toolbar's Theme panel): what to apply to the selected part of which
// representations.
const SelectionThemeParams = {
  action: PD.MappedStatic(
    'color',
    {
      color: PD.Group({ color: PD.Color(0x0090ff as any, { isEssential: true }) }),
      'reset-color': PD.EmptyGroup(),
      transparency: PD.Group({ value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { isEssential: true }) }),
      emissive: PD.Group({ value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { isEssential: true }) }),
      material: PD.Group({
        metalness: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }, { isEssential: true }),
        roughness: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { isEssential: true }),
      }),
      'reset-material': PD.EmptyGroup(),
      clipping: PD.Group({
        groups: PD.Select('1', [
          ['1', 'Group 1'],
          ['2', 'Group 2'],
          ['3', 'Group 3'],
        ]),
      }),
      wiggle: PD.Group({ amplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, { isEssential: true }) }),
    },
    {
      label: 'Action',
      isEssential: true,
      options: [
        ['color', 'Color'],
        ['reset-color', 'Reset color'],
        ['transparency', 'Transparency'],
        ['emissive', 'Emissive'],
        ['material', 'Material'],
        ['reset-material', 'Reset material'],
        ['clipping', 'Clipping'],
        ['wiggle', 'Wiggle'],
      ] as [any, string][],
    },
  ),
  representation: PD.Select(
    'all',
    [
      ['all', 'All'],
      ['cartoon', 'Cartoon'],
      ['ball-and-stick', 'Ball & Stick'],
      ['molecular-surface', 'Molecular Surface'],
    ],
    { label: 'Representations', isEssential: true },
  ),
};

/**
 * The selection (12 residues of 1TQN) with Mol*'s selection actions; Color and Theme open the Selection card (in the
 * panel's place) with the action preset.
 */
export function mockSelection(shell: ShellModel, post: (text: string) => void): SelectionModel {
  const summary = observableState<string | undefined>('12 residues');
  const form = new ParamsFormModel(SelectionThemeParams);
  shell.registerOverlays([
    {
      id: 'selection',
      title: 'Selection',
      render: domView(() => (
        <div style="display: grid; gap: 8px; padding: 4px 12px 12px">
          <ParamsForm form={form} />
          <div style="display: flex">
            <Button
              variant="primary"
              block
              icon={<Icon name="paintbrush" />}
              onClick={() => {
                const v = form.values.value as any;
                post(`Applied ${v.action.name} to ${summary.value} (${v.representation})`);
              }}
            >
              Apply
            </Button>
          </div>
        </div>
      )),
    },
  ]);
  const openWith = (action?: string) => {
    if (action)
      form.set('action', {
        name: action,
        params: PD.getDefaultValues((SelectionThemeParams.action.map(action) as PD.Group<any>).params),
      } as any);
    shell.openOverlay('selection');
  };
  return {
    summary: summary.asReadonly,
    actions: [
      { id: 'color', label: 'Color', icon: 'palette', run: () => openWith('color') },
      { id: 'theme', label: 'Transparency, material, …', icon: 'paintbrush', run: () => openWith() },
      {
        id: 'component',
        label: 'Create component',
        icon: 'layers',
        run: () => post(`Component from ${summary.value}`),
      },
      { id: 'focus', label: 'Focus', icon: 'focus', run: () => post(`Focused on ${summary.value}`) },
    ],
    clear: () => summary.set(undefined),
  };
}
