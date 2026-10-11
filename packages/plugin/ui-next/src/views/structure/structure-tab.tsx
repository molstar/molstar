/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Choices } from '../../components/choices/choices.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { Show } from '../../dom/control-flow.js';
import type { PresetsModel, QuickStylesModel, SceneInspection, SceneModel } from '../../model/scene.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel } from '../../shell/shell.js';
import { SectionList } from '../shell/section-list.js';
import { classes } from './_generated/structure-tab.module.css.js';
import { SceneTree } from './scene-tree.js';

export interface StructureTabProps {
  shell: ShellModel;
  scene: SceneModel;
  /** The header's actions without a model here (presets, styles, adding components): commands. */
  onAction?: (action: 'preset' | 'styles' | 'add') => void;
  /** The representation presets: the Preset button opens them under the header. */
  presets?: PresetsModel;
  /** Quick styles: the Styles button opens them under the header. */
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
  // One inline block at a time under the header: the preset menu or the quick styles.
  const open = observableState<'preset' | 'styles' | undefined>(undefined);
  const toggle = (block: 'preset' | 'styles') => open.set(open.value === block ? undefined : block);
  return (
    <div class={classes.tab}>
      <div class={classes.actions}>
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="sparkles" />}
          pressed={open.asReadonly.pipe(map((o) => o === 'preset'))}
          onClick={() => (props.presets ? toggle('preset') : props.onAction?.('preset'))}
        >
          Preset
        </Button>
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="paintbrush" />}
          pressed={open.asReadonly.pipe(map((o) => o === 'styles'))}
          onClick={() => (props.quickStyles ? toggle('styles') : props.onAction?.('styles'))}
        >
          Styles
        </Button>
        <Button size="sm" variant="ghost" icon={<Icon name="plus" />} onClick={() => props.onAction?.('add')}>
          Component
        </Button>
      </div>
      <Show when={open}>
        {(block) =>
          block === 'preset' && props.presets ? (
            <PresetMenu model={props.presets} onApplied={() => open.set(undefined)} />
          ) : block === 'styles' && props.quickStyles ? (
            <QuickStyles model={props.quickStyles} />
          ) : null
        }
      </Show>
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

/**
 * The preset menu: every representation preset, grouped as the registry groups them, the one in effect marked; a click
 * applies it. The presets' options (their params) fold out below, applied with the next preset.
 */
function PresetMenu(props: { model: PresetsModel; onApplied(): void }) {
  const { model } = props;
  const optionsOpen = observableState(false);
  return (
    <section class={classes.presets} aria-label="Presets">
      {model.groups.map((group) => (
        <div class={classes.presetGroup} role="group" aria-label={group.label ?? 'Presets'}>
          {group.label ? <span class={classes.presetCaption}>{group.label}</span> : null}
          {group.presets.map((p) => (
            <button
              type="button"
              class={classes.preset}
              aria-pressed={model.current.pipe(map((c) => c === p.id))}
              title={p.description}
              onClick={() => {
                model.apply(p.id);
                props.onApplied();
              }}
            >
              <span>{p.label}</span>
              <Icon name="check" size={14} class={classes.presetCheck} />
            </button>
          ))}
        </div>
      ))}
      {model.options ? (
        <div class={classes.presetOptions}>
          <button
            type="button"
            class={classes.presetOptionsToggle}
            aria-expanded={optionsOpen}
            onClick={() => optionsOpen.set(!optionsOpen.value)}
          >
            <Icon name="sliders-horizontal" size={14} />
            Options
            <Icon name="chevron-down" size={13} class={classes.chevron} />
          </button>
          <Show when={optionsOpen}>{() => <ParamsForm form={model.options!} />}</Show>
        </div>
      ) : null}
    </section>
  );
}
