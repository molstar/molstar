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
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { PresetsModel, QuickStylesModel, SceneInspection, SceneModel, SceneNode } from '../../model/scene.js';
import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import type { ShellModel } from '@molstar/prism-ui/shell/shell';
import { SectionList } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { VolumeLevels } from '../volume/volume-levels.js';
import { RepresentationQuick } from './representation-quick.js';
import { VolumeTransfer } from '../volume/volume-transfer.js';
import { classes } from './_generated/structure-tab.module.css.js';
import { NodeActions, SceneTree } from './scene-tree.js';
import { type SceneRowPick, SceneRowPickEvent } from './scene-row.js';

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
  // An extension's row picked (its controls in the inspector); a node picked takes its place.
  const pickedRow = observableState<SceneRowPick | undefined>(undefined);
  const picks = scene.picked.subscribe((n) => {
    if (!n || !pickedRow.value) return;
    pickedRow.value.mark(false);
    pickedRow.set(undefined);
  });
  onCleanup(() => {
    picks.unsubscribe();
    pickedRow.value?.mark(false);
  });
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
      <ScrollArea
        edge
        class={classes.tree}
        ref={(el: HTMLElement) => {
          // An extension's row picked (Membrane, Assembly symmetry): its controls in the inspector, as a node's.
          const onPick = (e: Event) => {
            e.preventDefault();
            const pick = (e as CustomEvent<SceneRowPick>).detail;
            pickedRow.value?.mark(false);
            scene.pick(undefined);
            pick.mark(true);
            pickedRow.set(pick);
          };
          el.addEventListener(SceneRowPickEvent, onPick);
          onCleanup(() => el.removeEventListener(SceneRowPickEvent, onPick));
        }}
      >
        <SceneTree scene={scene} shell={shell} filter={search.asReadonly.pipe(debounceTime(120))} />
      </ScrollArea>
      <ScrollArea
        class={classes.inspector}
        contentClass={classes.inspectorContent}
        role="region"
        aria-label="Inspector"
      >
        <Show
          when={pickedRow.asReadonly}
          fallback={() => (
            <Show
              when={scene.picked}
              fallback={() => <p class={classes.inspectorEmpty}>Select an item in the tree to edit its properties.</p>}
            >
              {(node) => <Inspector shell={shell} inspection={scene.inspect(node)} onPick={(n) => scene.pick(n)} />}
            </Show>
          )}
        >
          {(row) => <RowInspector row={row} />}
        </Show>
      </ScrollArea>
    </div>
  );
}

/** An extension's row in the inspector: its label, its line, Turn on while off, and its controls. */
function RowInspector(props: { row: SceneRowPick }) {
  const { row } = props;
  const off = row.shown ? row.shown.pipe(map((s) => !s)) : undefined;
  return (
    <div class={classes.inspection}>
      <header class={classes.inspectorHeader}>
        <span class={classes.inspectorTitle}>{row.label}</span>
        <span class={classes.inspectorKind}>{row.description ?? ''}</span>
      </header>
      {off && row.onToggle ? (
        <Show when={off}>
          {() => (
            <div class={classes.inspectorActions}>
              <button type="button" class={classes.inspectorAction} onClick={() => row.onToggle?.()}>
                <span>Turn on</span>
              </button>
            </div>
          )}
        </Show>
      ) : null}
      {row.body ? row.body() : <p class={classes.nothing}>Nothing to set here.</p>}
    </div>
  );
}

/** Rows that stand for a feature (on or off), not for something drawn (shown or hidden). */
const TurnsOn: ReadonlySet<string> = new Set(['extension', 'unit-cell']);

function Inspector(props: { shell: ShellModel; inspection: SceneInspection; onPick(node: SceneNode): void }) {
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
      {/* Off or hidden (the streaming not on, the unit cell not shown): turning it on is here too, not only its row's
          eye. */}
      <Show when={node.isHidden}>
        {() => (
          <div class={classes.inspectorActions}>
            <button type="button" class={classes.inspectorAction} onClick={() => node.toggleVisibility()}>
              <span>{TurnsOn.has(node.kind) ? 'Turn on' : 'Show'}</span>
            </button>
          </div>
        )}
      </Show>
      {props.inspection.actions?.length ? (
        <div class={classes.inspectorActions}>
          {props.inspection.actions.map((a) => (
            <button type="button" class={classes.inspectorAction} onClick={() => a.run?.()}>
              <span>{a.label}</span>
            </button>
          ))}
        </div>
      ) : null}
      {/* A component's representations: each a row to open (its inspector), its eye beside it; not only behind the
          tree's chevron. */}
      {node.kind === 'component' ? <Representations node={node} onPick={props.onPick} /> : null}
      {/* A representation's quick settings (type, color, opacity, size), above the whole form. */}
      {props.inspection.quick ? <RepresentationQuick model={props.inspection.quick} /> : null}
      {/* A map's levels (the iso values and colors people change most), above the whole form. */}
      {props.inspection.levels?.map((levels) => (
        <VolumeLevels model={levels} />
      ))}
      {/* A direct volume's transfer function (opacity and colors over the map's values). */}
      {props.inspection.transfer?.map((transfer) => (
        <VolumeTransfer model={transfer} />
      ))}
      {form ? <ParamsForm form={form} /> : null}
      <SectionList shell={props.shell} slot={slot} subject={node} flat />
      {/* Nothing to set (no form or actions of its own): what can be done with it is shown here, its row menu's
          entries; with none, say so rather than leave an empty space. */}
      {!form && !props.inspection.actions?.length ? (
        (node.menu?.length ?? 0) > 0 || node.removable !== false ? (
          <NodeActions node={node} />
        ) : !props.inspection.facts?.length ? (
          <Show when={sectionsFor(props.shell, slot, node).pipe(map((has) => !has))}>
            {() => (
              <Show when={node.isHidden.pipe(map((h) => !h))}>
                {() => <p class={classes.nothing}>Nothing to set here.</p>}
              </Show>
            )}
          </Show>
        ) : null
      ) : null}
    </div>
  );
}

/** A component's representations in its inspector: icon, name, what it is colored by, its eye; a click opens it. */
function Representations(props: { node: SceneNode; onPick(node: SceneNode): void }) {
  const reprs = props.node.children.pipe(map((c) => c.filter((n) => n.kind === 'representation')));
  return (
    <Show when={reprs.pipe(map((r) => r.length > 0))}>
      {() => (
        <section class={classes.children} aria-label="Representations">
          <span class={classes.childrenTitle}>Representations</span>
          <For each={reprs} key={(n) => n.id}>
            {(n) => (
              <div class={classes.child} classList={{ [classes.childHidden]: n.isHidden }}>
                <button type="button" class={classes.childOpen} onClick={() => props.onPick(n)}>
                  <Icon name={n.icon} size={15} class={classes.childIcon} />
                  <span class={classes.childLabel}>{n.label}</span>
                  <span class={classes.childDescription}>{n.description}</span>
                  <Icon name="chevron-right" size={14} class={classes.childChevron} />
                </button>
                <button
                  type="button"
                  class={classes.childEye}
                  aria-label={n.isHidden.pipe(map((h) => (h ? `Show ${n.label.value}` : `Hide ${n.label.value}`)))}
                  data-tip={n.isHidden.pipe(map((h) => (h ? 'Show' : 'Hide')))}
                  onClick={() => n.toggleVisibility()}
                >
                  <Show when={n.isHidden} fallback={() => <Icon name="eye" size={14} />}>
                    {() => <Icon name="eye-off" size={14} />}
                  </Show>
                </button>
              </div>
            )}
          </For>
        </section>
      )}
    </Show>
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
