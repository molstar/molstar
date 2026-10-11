/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsForm } from '../../components/params/params-form.js';
import { domView } from '../../dom/view.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { SceneNode } from '../../model/scene.js';
import type { ShellModel } from '../../shell/shell.js';
import type { MockScene } from './scene.js';

/**
 * Extensions' inspector sections for structures, as their UI contributions would register them: volume streaming
 * (density for X-ray and EM entries) and the RCSB assembly symmetry (only for a symmetric assembly).
 */
export function registerMockInspectorSections(shell: ShellModel, scene: MockScene) {
  const symmetryForms = new Map<string, ParamsFormModel<any>>();
  const symmetryForm = (node: SceneNode) => {
    let form = symmetryForms.get(node.id);
    if (!form) {
      const symmetry = scene.structureInfo(node)?.symmetry ?? 'C1';
      form = new ParamsFormModel({
        symmetry: PD.Select('global', [
          ['global', `Global · ${symmetry}`],
          ['pseudo', `Pseudo · ${symmetry}`],
        ]),
        axes: PD.Boolean(true, { label: 'Show axes' }),
        color: PD.Boolean(false, { label: 'Color by symmetry' }),
      });
      symmetryForms.set(node.id, form);
    }
    return form;
  };

  shell.registerSections([
    {
      id: 'volume-streaming',
      slot: 'inspector:structure',
      title: 'Density',
      icon: 'volume',
      keywords: ['volume streaming', 'electron density', 'map', '2fo-fc', 'em'],
      order: 0,
      // Enabling it streams the entry's map around the focus; the channels appear under the structure as volumes.
      when: (ctx) => !!ctx.subject && !!scene.structureInfo(ctx.subject as SceneNode),
      render: domView(({ ctx }) => <ParamsForm form={scene.streaming(ctx.subject as SceneNode)} />),
    },
    {
      id: 'assembly-symmetry',
      slot: 'inspector:structure',
      title: 'Assembly symmetry',
      icon: 'orbit',
      keywords: ['symmetry', 'axes', 'rcsb'],
      order: 1,
      // Annotated by RCSB PDB, nothing to compute: which symmetry, its axes, coloring by it. Hidden for C1.
      when: (ctx) => !!ctx.subject && !!scene.structureInfo(ctx.subject as SceneNode)?.symmetry,
      render: domView(({ ctx }) => <ParamsForm form={symmetryForm(ctx.subject as SceneNode)} />),
    },
  ]);
}
