/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { Show } from '../../dom/control-flow.js';
import type { SceneInspection, SceneModel } from '../../model/scene.js';
import type { ShellModel } from '../../shell/shell.js';
import { SectionList } from '../shell/section-list.js';
import { classes } from './_generated/structure-tab.module.css.js';
import { SceneTree } from './scene-tree.js';

export interface StructureTabProps {
  shell: ShellModel;
  scene: SceneModel;
  /** The header's actions (presets, coloring, adding components): commands, shown as buttons. */
  onAction?: (action: 'preset' | 'color' | 'add') => void;
}

/**
 * The main tab: the scene tree, then the inspector for the picked node (its parameters, then sections that extensions
 * contribute for that kind of node: assembly symmetry, volume streaming, …).
 */
export function StructureTab(props: StructureTabProps) {
  const { shell, scene } = props;
  return (
    <div class={classes.tab}>
      <div class={classes.actions}>
        <Button size="sm" variant="ghost" icon={<Icon name="sparkles" />} onClick={() => props.onAction?.('preset')}>
          Preset
        </Button>
        <Button size="sm" variant="ghost" icon={<Icon name="palette" />} onClick={() => props.onAction?.('color')}>
          Color
        </Button>
        <Button size="sm" variant="ghost" icon={<Icon name="plus" />} onClick={() => props.onAction?.('add')}>
          Component
        </Button>
      </div>
      <SceneTree scene={scene} />
      <Show when={scene.picked}>{(node) => <Inspector shell={shell} inspection={scene.inspect(node)} />}</Show>
    </div>
  );
}

function Inspector(props: { shell: ShellModel; inspection: SceneInspection }) {
  const { node, form, slot } = props.inspection;
  return (
    <section class={classes.inspector} aria-label="Inspector">
      <header class={classes.inspectorHeader}>
        <Icon name={node.icon} size={16} />
        <span class={classes.inspectorTitle}>{node.label}</span>
        <span class={classes.inspectorKind}>{node.kind.replace('-', ' ')}</span>
      </header>
      {form ? <ParamsForm form={form} /> : null}
      <SectionList shell={props.shell} slot={slot} />
    </section>
  );
}
