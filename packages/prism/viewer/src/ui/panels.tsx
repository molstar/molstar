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

// The core's shared homes for features' sections (viewer-spec.md §3.1, the `ViewerSlot`s): the Selection tab
// (`selection`: what to do with the selection; Measure and Superpose come with their features) and the Export card
// (`export`: image and session with their features, geometry, models, video with Mol*'s extensions). The Export
// button shows once something exports.

/** The Selection tab: its slot's sections, folded until one is opened (the selection bar opens one). */
export function selectionTab(selection: ReadonlyState<string | undefined>): ShellTab {
  return {
    id: 'selection',
    title: 'Selection',
    icon: 'mouse-pointer-click',
    order: 1.2,
    render: domView(({ ctx }) => (
      <div style="display: grid; gap: var(--msp-space-2)">
        {/* What the sections below act on. */}
        <p style="margin: 0; padding: var(--msp-space-2) var(--msp-space-3); color: var(--msp-color-text-muted)">
          {selection.pipe(
            map((s) => (s ? `Selected: ${s}` : 'Nothing is selected: pick on the canvas or in the sequence.')),
          )}
        </p>
        <SectionList shell={ctx.shell} slot="selection" focus={ctx.focus} collapsed />
      </div>
    )),
  };
}

/** The core's sections in the Selection tab: style what is selected, make a component of it. */
export function selectionSections(api: ViewerApi, selection: ReadonlyState<string | undefined>): ShellSection[] {
  // Nothing to act on without a selection: the buttons say so.
  const none = selection.pipe(map((s) => !s));
  return [
    {
      id: 'selection-style',
      slot: 'selection',
      title: 'Color, transparency, glow',
      icon: 'palette',
      order: 0,
      render: domView(() => <StyleSelection api={api} none={none} />),
    },
    {
      id: 'selection-component',
      slot: 'selection',
      title: 'Make a component',
      icon: 'layers',
      order: 1,
      render: domView(() => <ComponentOfSelection api={api} none={none} />),
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
          { alpha: PD.Numeric(0.5, { min: 0, max: 1, step: 0.05 }, { label: 'Opacity' }) },
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
        onReset={() => form.reset()}
        onApply={apply}
      />
    </div>
  );
}

function ComponentOfSelection(props: { api: ViewerApi; none: Observable<boolean> }) {
  const { api } = props;
  // The representation types, as `representation.set` lists them (`component.add` takes one by id).
  const types = (api.describe('representation.set')?.params.type as PD.Select<string> | undefined)?.options ?? [];
  const form = new ParamsFormModel({
    label: PD.Text('Selection', { label: 'Name' }),
    representation: PD.Select(types.find((t) => t[0] === 'ball-and-stick')?.[0] ?? types[0]?.[0] ?? '', types),
  });
  const busy = observableState(false);
  const apply = async () => {
    const v = form.values.value as { label: string; representation: string };
    busy.set(true);
    try {
      await api.component.add({ target: 'selection', label: v.label, representation: v.representation });
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
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="export" focus={ctx.focus} />),
    },
    toolbar: {
      id: 'export',
      group: 'export',
      icon: 'download',
      title: 'Export…',
      overlay: 'export',
      visible: shell.sections('export').pipe(map((list) => list.length > 0)),
    },
  };
}
