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
import { SectionList } from '../../views/shell/section-list.js';

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
  // The Selection card: what to do with it, as sections (like Export): apply a theme action (color, transparency,
  // material, …), or make a component of it.
  const component = new ParamsFormModel({
    label: PD.Text('Selection', { label: 'Label', isEssential: true }),
    representation: PD.Select(
      'ball-and-stick',
      [
        ['cartoon', 'Cartoon'],
        ['ball-and-stick', 'Ball & Stick'],
        ['spacefill', 'Spacefill'],
        ['molecular-surface', 'Molecular Surface'],
      ],
      { label: 'Representation', isEssential: true },
    ),
  });
  shell.registerOverlays([
    {
      id: 'selection',
      title: 'Selection',
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="selection" focus={ctx.focus} />),
    },
  ]);
  shell.registerSections([
    {
      id: 'theme',
      slot: 'selection',
      title: 'Color, transparency, material…',
      icon: 'palette',
      order: 0,
      render: domView(() => (
        <div style="display: grid; gap: 8px">
          <ParamsForm form={form} />
          <div style="display: flex; justify-content: flex-end">
            <Button
              variant="primary"
              size="sm"
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
    {
      id: 'component',
      slot: 'selection',
      title: 'Create component',
      icon: 'layers',
      order: 1,
      render: domView(() => (
        <div style="display: grid; gap: 8px">
          <ParamsForm form={component} />
          <div style="display: flex; justify-content: flex-end">
            <Button
              variant="primary"
              size="sm"
              icon={<Icon name="plus" />}
              onClick={() => post(`Component "${(component.values.value as any).label}" from ${summary.value}`)}
            >
              Create
            </Button>
          </div>
        </div>
      )),
    },
  ]);
  return {
    summary: summary.asReadonly,
    // One way in: the card offers what to do (color, theme, component); focus stays one click.
    actions: [
      { id: 'edit', label: 'Color, theme, component…', icon: 'paintbrush', run: () => shell.openOverlay('selection') },
      { id: 'focus', label: 'Focus', icon: 'focus', run: () => post(`Focused on ${summary.value}`) },
    ],
    clear: () => summary.set(undefined),
  };
}
