/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Choices } from '../../components/choices/choices.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { Show } from '../../dom/control-flow.js';
import type { QuickStylesModel, SceneInspection, SceneModel } from '../../model/scene.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel } from '../../shell/shell.js';
import { SectionList } from '../shell/section-list.js';
import { classes } from './_generated/structure-tab.module.css.js';
import { SceneTree } from './scene-tree.js';

export interface StructureTabProps {
  shell: ShellModel;
  scene: SceneModel;
  /** The header's actions (presets, coloring, adding components): commands, shown as buttons. */
  onAction?: (action: 'preset' | 'color' | 'add') => void;
  /** Quick styles: the Preset button opens them under the header. */
  quickStyles?: QuickStylesModel;
}

/**
 * The main tab: the scene tree, then the inspector for the picked node (its parameters, then sections that extensions
 * contribute for that kind of node: assembly symmetry, volume streaming, …). The tab fills the panel: the tree scrolls
 * on its own, and the inspector is always at the bottom at a fixed height (scrolling inside; a hint while nothing is
 * picked), so picking never moves the tree.
 */
export function StructureTab(props: StructureTabProps) {
  const { shell, scene } = props;
  const stylesOpen = observableState(false);
  return (
    <div class={classes.tab}>
      <div class={classes.actions}>
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="sparkles" />}
          pressed={stylesOpen}
          onClick={() => (props.quickStyles ? stylesOpen.set(!stylesOpen.value) : props.onAction?.('preset'))}
        >
          Preset
        </Button>
        <Button size="sm" variant="ghost" icon={<Icon name="palette" />} onClick={() => props.onAction?.('color')}>
          Color
        </Button>
        <Button size="sm" variant="ghost" icon={<Icon name="plus" />} onClick={() => props.onAction?.('add')}>
          Component
        </Button>
      </div>
      {props.quickStyles ? <Show when={stylesOpen}>{() => <QuickStyles model={props.quickStyles!} />}</Show> : null}
      <div class={classes.tree}>
        <SceneTree scene={scene} />
      </div>
      <section class={classes.inspector} aria-label="Inspector">
        <Show
          when={scene.picked}
          fallback={() => <p class={classes.inspectorEmpty}>Select an item in the tree to edit its properties.</p>}
        >
          {(node) => <Inspector shell={shell} inspection={scene.inspect(node)} />}
        </Show>
      </section>
    </div>
  );
}

function Inspector(props: { shell: ShellModel; inspection: SceneInspection }) {
  const { node, form, slot } = props.inspection;
  return (
    <div class={classes.inspection}>
      <header class={classes.inspectorHeader}>
        <Icon name={node.icon} size={15} />
        <span class={classes.inspectorTitle}>{node.label}</span>
        <span class={classes.inspectorKind}>{node.kind.replace('-', ' ')}</span>
      </header>
      {props.inspection.facts?.length ? (
        <dl class={classes.facts}>
          {props.inspection.facts.map(([label, value]) => (
            <div class={classes.fact}>
              <dt>{label}</dt>
              <dd>{value}</dd>
            </div>
          ))}
        </dl>
      ) : null}
      {form ? <ParamsForm form={form} /> : null}
      <SectionList shell={props.shell} slot={slot} subject={node} flat />
    </div>
  );
}

/** One click restyles everything: a representation preset, or a visual style. */
function QuickStyles(props: { model: QuickStylesModel }) {
  const { model } = props;
  return (
    <section class={classes.quickStyles} aria-label="Quick styles">
      <Choices label="Representation" choices={model.representations} onPick={(id) => model.applyRepresentation(id)} />
      <Choices label="Style" choices={model.styles} onPick={(id) => model.applyStyle(id)} />
    </section>
  );
}
