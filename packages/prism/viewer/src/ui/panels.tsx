/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  ActionFooter,
  type CanvasToolbarItem,
  domView,
  onCleanup,
  ParamsForm,
  SectionList,
  type ShellModel,
  type ShellOverlay,
  type ShellSection,
  type ShellTab,
} from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { map, type Observable } from 'rxjs';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type { SelectionSummaryModel } from '../model/selection.js';
import type { StyleLayerInfo } from '../sources/style-layers.js';
import { SelectionSummary } from '../views/selection/selection-summary.js';
import { StyleLayers } from '../views/selection/style-layers.js';
import { createStyleLayersModel } from './style-layers.js';

// The core's shared homes for features' sections (viewer-spec.md §3.1, the `ViewerSlot`s): the Selection tab
// (`selection`: what to do with the selection; Measure and Superpose come with their features) and the Export card
// (`export`: image and session with their features, geometry, models, video with Mol*'s extensions). The Export
// button shows once something exports.

/** The Selection tab: what is selected, then its slot's sections, folded until one is opened (the selection bar opens one). */
export function selectionTab(selection: SelectionSummaryModel): ShellTab {
  return {
    id: 'selection',
    title: 'Selection',
    icon: 'mouse-pointer-click',
    order: 1.2,
    render: domView(({ ctx }) => (
      <div style="display: grid; gap: var(--msp-space-2)">
        {/* What the sections below act on. */}
        <SelectionSummary model={selection} />
        <SectionList shell={ctx.shell} slot="selection" focus={ctx.focus} collapsed />
      </div>
    )),
  };
}

/**
 * The core's sections in the Selection tab: style what is selected (and the style layers applied, `sources.layers`),
 * make a component of it.
 */
export function selectionSections(
  api: ViewerApi,
  selection: ReadonlyState<string | undefined>,
  layers: ReadonlyState<readonly StyleLayerInfo[]>,
  /** Shows a component just made where it is: the Structure panel, its row picked. */
  showComponent?: (ref: string) => void,
): ShellSection[] {
  // Nothing to act on without a selection: the buttons say so.
  const none = selection.pipe(map((s) => !s));
  return [
    {
      id: 'selection-style',
      slot: 'selection',
      title: 'Color, transparency, glow',
      icon: 'palette',
      keywords: ['overpaint', 'color', 'transparency', 'emissive', 'glow', 'layers', 'clear'],
      order: 0,
      render: domView(() => {
        const host = createStyleLayersModel(api, layers);
        onCleanup(() => host.dispose());
        return (
          <div style="display: grid; grid-template-columns: minmax(0, 1fr)">
            <StyleSelection api={api} none={none} />
            <StyleLayers model={host.model} />
          </div>
        );
      }),
    },
    {
      id: 'selection-component',
      slot: 'selection',
      title: 'Make a component',
      icon: 'layers',
      order: 1,
      render: domView(() => (
        <ComponentOfSelection api={api} none={none} {...(showComponent ? { show: showComponent } : {})} />
      )),
    },
  ];
}

function StyleSelection(props: { api: ViewerApi; none: Observable<boolean> }) {
  const { api } = props;
  const themes = (api.describe('representation.color')?.params.theme as PD.Select<string> | undefined)?.options ?? [];
  const form = new ParamsFormModel({
    action: PD.MappedStatic(
      'color',
      {
        color: PD.Group({ color: PD.Color(Color(0xe0523a)) }, { isFlat: true }),
        theme: PD.Group(
          {
            theme: PD.Select(
              themes.find((t) => t[0])?.[0] ?? '',
              themes.filter((t) => t[0]),
              { label: 'Color by' },
            ),
          },
          { isFlat: true },
        ),
        transparency: PD.Group(
          {
            alpha: PD.Numeric(
              0.5,
              { min: 0, max: 1, step: 0.05 },
              { label: 'Transparency', description: '0: opaque, 1: hidden' },
            ),
          },
          {
            isFlat: true,
          },
        ),
        glow: PD.Group({ strength: PD.Numeric(0.5, { min: 0, max: 1, step: 0.05 }) }, { isFlat: true }),
      },
      { label: 'Change' },
    ),
  });
  const busy = observableState(false);
  const apply = async () => {
    const a = (form.values.value as { action: { name: string; params: Record<string, any> } }).action;
    busy.set(true);
    try {
      if (a.name === 'color') {
        await api.representation.color({ target: 'selection', color: Color.toHexStyle(a.params.color) });
      } else if (a.name === 'theme') {
        await api.representation.color({ target: 'selection', theme: a.params.theme });
      } else if (a.name === 'transparency') {
        await api.representation.transparency({ target: 'selection', alpha: a.params.alpha });
      } else {
        await api.representation.emissive({ target: 'selection', strength: a.params.strength });
      }
    } finally {
      busy.set(false);
    }
  };
  return (
    <div style="display: grid; gap: var(--msp-space-3)">
      <ParamsForm form={form} />
      <ActionFooter
        icon="paintbrush"
        label="Apply"
        busy={busy.asReadonly}
        disabled={props.none}
        onReset={() => form.reset()}
        onApply={apply}
      />
    </div>
  );
}

function ComponentOfSelection(props: { api: ViewerApi; none: Observable<boolean>; show?: (ref: string) => void }) {
  const { api } = props;
  // The representation types, as `representation.set` lists them (`component.add` takes one by id); or none, to add one
  // later (v5's "< Create Later >").
  const types = (api.describe('representation.set')?.params.type as PD.Select<string> | undefined)?.options ?? [];
  const choices = [['none', 'None (add one later)'] as const, ...types] as PD.Select<string>['options'];
  const form = new ParamsFormModel({
    label: PD.Text('Selection', { label: 'Name' }),
    representation: PD.Select(types.find((t) => t[0] === 'ball-and-stick')?.[0] ?? types[0]?.[0] ?? '', choices),
  });
  const busy = observableState(false);
  const apply = async () => {
    const v = form.values.value as { label: string; representation: string };
    busy.set(true);
    try {
      const r = await api.component.add({ target: 'selection', label: v.label, representation: v.representation });
      // Where it is: the Structure panel, its row picked (as a person looks for what they just made).
      const made = r.ok ? (r.value as readonly { ref: string }[] | undefined)?.[0] : undefined;
      if (made) props.show?.(made.ref);
    } finally {
      busy.set(false);
    }
  };
  return (
    <div style="display: grid; gap: var(--msp-space-3)">
      <ParamsForm form={form} />
      <ActionFooter
        icon="plus"
        label="Make"
        busy={busy.asReadonly}
        disabled={props.none}
        onReset={() => form.reset()}
        onApply={apply}
      />
    </div>
  );
}

/** The Export card and its button: there once a section is in the `export` slot. */
export function exportCard(shell: ShellModel): { overlay: ShellOverlay; toolbar: CanvasToolbarItem } {
  return {
    overlay: {
      id: 'export',
      title: 'Export',
      // Folded until opened, but the image (`open`): the card's first use.
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="export" focus={ctx.focus} collapsed />),
    },
    toolbar: {
      id: 'export',
      // Last in the scene's group, under Canvas, Animation and the Story (David, 2026-10-15).
      group: 'scene',
      order: 100,
      icon: 'download',
      title: 'Export…',
      overlay: 'export',
      visible: shell.sections('export').pipe(map((list) => list.length > 0)),
    },
  };
}
