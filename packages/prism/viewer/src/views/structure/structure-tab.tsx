/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { debounceTime, map } from 'rxjs';
import { Button } from '@molstar/prism-ui';
import { Choices } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ParamsForm } from '@molstar/prism-ui';
import { ScrollArea } from '@molstar/prism-ui';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import type { PresetsModel, QuickStylesModel, SceneInspection, SceneModel } from '../../model/scene.js';
import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ShellModel } from '@molstar/prism-ui/shell/shell';
import { SectionList } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { classes } from './_generated/structure-tab.module.css.js';
import { NodeActions, SceneTree } from './scene-tree.js';

export interface StructureTabProps {
  shell: ShellModel;
  scene: SceneModel;
  /** The header's actions without a model here (presets, styles, adding components): commands. */
  onAction?: (action: 'preset' | 'styles' | 'add') => void;
  /** The representation presets: the Preset button opens them under the header. */
  presets?: PresetsModel;
  /** Quick styles: the Styles button opens them under the header. */
  quickStyles?: QuickStylesModel;
  /**
   * Opens the command menu inside a group: given, Preset and Styles open theirs there (teaching the command menu as they
   * go) instead of their inline lists.
   */
  openCommands?: (group: 'preset' | 'styles') => void;
  /**
   * The component options (Mol*'s component manager options: hydrogens, quality, lighting, material, clipping,
   * interactions, interior, animation), applied to every component; the sliders button opens them.
   */
  options?: ParamsFormModel<any>;
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
  const open = observableState<'preset' | 'styles' | 'options' | undefined>(undefined);
  const toggle = (block: 'preset' | 'styles' | 'options') => open.set(open.value === block ? undefined : block);
  // The tree's search: a field over the tree while open; closed, no filter.
  const searchOpen = observableState(false);
  const search = observableState('');
  const closeSearch = () => {
    search.set('');
    searchOpen.set(false);
  };
  return (
    <div class={classes.tab}>
      <div class={classes.actions}>
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="shapes" />}
          pressed={open.asReadonly.pipe(map((o) => o === 'preset'))}
          onClick={() =>
            props.openCommands
              ? props.openCommands('preset')
              : props.presets
                ? toggle('preset')
                : props.onAction?.('preset')
          }
        >
          <span class={classes.actionLabel}>Preset</span>
        </Button>
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="paintbrush" />}
          pressed={open.asReadonly.pipe(map((o) => o === 'styles'))}
          onClick={() =>
            props.openCommands
              ? props.openCommands('styles')
              : props.quickStyles
                ? toggle('styles')
                : props.onAction?.('styles')
          }
        >
          <span class={classes.actionLabel}>Styles</span>
        </Button>
        <Button size="sm" variant="ghost" icon={<Icon name="plus" />} onClick={() => props.onAction?.('add')}>
          <span class={classes.actionLabel}>Component</span>
        </Button>
        <Button
          size="sm"
          variant="ghost"
          class={classes.searchButton}
          icon={<Icon name="search" />}
          title="Search the scene"
          pressed={searchOpen.asReadonly}
          onClick={() => (searchOpen.value ? closeSearch() : searchOpen.set(true))}
        />
        {props.options ? (
          <Button
            size="sm"
            variant="ghost"
            class={classes.optionsButton}
            icon={<Icon name="settings" />}
            title="Component options"
            pressed={open.asReadonly.pipe(map((o) => o === 'options'))}
            onClick={() => toggle('options')}
          />
        ) : null}
      </div>
      <Show when={open} reveal>
        {(block) =>
          block === 'preset' && props.presets ? (
            <PresetMenu model={props.presets} onApplied={() => open.set(undefined)} />
          ) : block === 'styles' && props.quickStyles ? (
            <QuickStyles model={props.quickStyles} />
          ) : block === 'options' && props.options ? (
            <ScrollArea
              class={classes.options}
              contentClass={classes.optionsContent}
              role="region"
              aria-label="Component options"
            >
              <ParamsForm form={props.options} />
            </ScrollArea>
          ) : null
        }
      </Show>
      <Show when={searchOpen.asReadonly}>
        {() => (
          <div class={classes.search} ref={(el) => queueMicrotask(() => el.querySelector('input')?.focus())}>
            <TextInput
              size="sm"
              search
              placeholder="Search the scene: an entity, a chain, a component…"
              aria-label="Search the scene"
              value={search.asReadonly}
              onInput={(v) => search.set(v)}
            />
          </div>
        )}
      </Show>
      <ScrollArea edge class={classes.tree}>
        <SceneTree scene={scene} shell={shell} filter={search.asReadonly.pipe(debounceTime(120))} />
      </ScrollArea>
      <ScrollArea
        class={classes.inspector}
        contentClass={classes.inspectorContent}
        role="region"
        aria-label="Inspector"
      >
        <Show
          when={scene.picked}
          fallback={() => <p class={classes.inspectorEmpty}>Select an item in the tree to edit its properties.</p>}
        >
          {(node) => <Inspector shell={shell} inspection={scene.inspect(node)} />}
        </Show>
      </ScrollArea>
    </div>
  );
}

function Inspector(props: { shell: ShellModel; inspection: SceneInspection }) {
  const { node, form, slot } = props.inspection;
  return (
    <div class={classes.inspection}>
      <header class={classes.inspectorHeader}>
        <span class={classes.inspectorTitle}>{node.label}</span>
        {/* The kind, unless the label already says it (Measurements). */}
        <Show when={node.label.pipe(map((l) => l.toLowerCase() !== node.kind.replace('-', ' ')))}>
          {() => <span class={classes.inspectorKind}>{node.kind.replace('-', ' ')}</span>}
        </Show>
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
      {props.inspection.actions?.length ? (
        <div class={classes.inspectorActions}>
          {props.inspection.actions.map((a) => (
            <button type="button" class={classes.inspectorAction} onClick={() => a.run?.()}>
              <span>{a.label}</span>
            </button>
          ))}
        </div>
      ) : null}
      {form ? <ParamsForm form={form} /> : null}
      <SectionList shell={props.shell} slot={slot} subject={node} flat />
      {/* Nothing to set (no form or actions of its own): what can be done with it is shown here, its row menu's
          entries; with none, say so rather than leave an empty space. */}
      {!form && !props.inspection.actions?.length ? (
        (node.menu?.length ?? 0) > 0 || node.removable !== false ? (
          <NodeActions node={node} />
        ) : !props.inspection.facts?.length ? (
          <Show when={sectionsFor(props.shell, slot, node).pipe(map((has) => !has))}>
            {() => <p class={classes.nothing}>Nothing to set here.</p>}
          </Show>
        ) : null
      ) : null}
    </div>
  );
}

/** Whether any contributed section shows for a node (`inspector:<kind>` sections with their `when`). */
function sectionsFor(shell: ShellModel, slot: string, subject: unknown) {
  return shell.sections(slot).pipe(map((all) => all.some((s) => !s.when || s.when({ shell, subject }))));
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
 * applies it.
 */
function PresetMenu(props: { model: PresetsModel; onApplied(): void }) {
  const { model } = props;
  return (
    <ScrollArea class={classes.presets} contentClass={classes.presetsContent} role="region" aria-label="Presets">
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
    </ScrollArea>
  );
}
