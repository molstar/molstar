/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ParamsForm } from '@molstar/prism-ui';
import { ScrollArea } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { SceneLookModel } from '../../model/look.js';
import type { SceneInspection, SceneModel, SceneNode } from '../../model/scene.js';
import { keyedState } from '@molstar/prism-ui/reactive/keyed-state';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import { effect, onActive, onCleanup, requireScope, runInScope } from '@molstar/prism-ui/dom/scope';
import { on, timeout } from '@molstar/prism-ui/dom/scoped';
import type { ShellModel } from '@molstar/prism-ui/shell/shell';
import { SectionList } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { VolumeLevels } from '../volume/volume-levels.js';
import { ContentsSummaryView } from './contents-summary.js';
import { SequenceInspector } from './sequence-inspector.js';
import { RepresentationQuick } from './representation-quick.js';
import { VolumeTransfer } from '../volume/volume-transfer.js';
import { classes } from './_generated/structure-tab.module.css.js';
import { LookButton, NodeLookPopover, SceneLookPopover, sceneOverridesText } from './look-popover.js';
import { menuOf, NodeActions, RowMenu, SceneTree } from './scene-tree.js';
import { type SceneRowPick, SceneRowPickEvent } from './scene-row.js';

export interface StructureTabProps {
  shell: ShellModel;
  scene: SceneModel;
  /**
   * Opens the command menu inside a group: Preset and Styles open theirs there (teaching the command menu as they go).
   * Given, the buttons show; the tab has no inline preset or style lists.
   */
  openCommands?: (group: 'preset' | 'styles') => void;
  /**
   * The scene defaults of the look (api.md §2, "Look": hydrogens, quality, lighting, material, …), which every
   * representation follows unless it has its own; the sliders button opens them, its dot telling of overrides.
   */
  look?: SceneLookModel;
}

/**
 * The main tab: the scene tree, then the inspector for the picked node (its parameters, then sections that extensions
 * contribute for that kind of node: assembly symmetry, volume streaming, …). The tab fills the panel: the tree scrolls
 * on its own, and the inspector is always at the bottom at a fixed height (scrolling inside; a hint while nothing is
 * picked), so picking never moves the tree.
 */
export function StructureTab(props: StructureTabProps) {
  const { shell, scene } = props;
  // The tree's search: a field over the tree while open; closed, no filter.
  // An extension's row picked (its controls in the inspector); a node picked takes its place.
  const pickedRow = observableState<SceneRowPick | undefined>(undefined);
  effect(scene.picked, (n) => {
    if (!n || !pickedRow.value) return;
    pickedRow.value.mark(false);
    pickedRow.set(undefined);
  });
  // The picked row gone from the tree: its inspector goes with it (what it controls may be gone).
  let gone: (() => void) | undefined;
  onCleanup(() => {
    gone?.();
    pickedRow.value?.mark(false);
  });
  const searchOpen = observableState(false);
  const search = observableState('');
  // The tree searches what was typed once the typing pauses (120 ms): one timer, the last key's, a scoped one (it ends
  // with the tab). Paused (the tab hidden), no timer runs: what was typed is searched when the tab shows again.
  const filter = observableState('');
  const owner = requireScope('StructureTab');
  let cancelTyping: (() => void) | undefined;
  effect(search.asReadonly, (text) => {
    cancelTyping?.();
    cancelTyping = undefined;
    // Back to what the tree searches (and at the mount): no timer.
    if (text === filter.value) return;
    cancelTyping = runInScope(owner, () =>
      timeout(() => {
        cancelTyping = undefined;
        filter.set(text);
      }, 120),
    );
  });
  onActive((active) => {
    if (active) {
      if (search.value !== filter.value && !cancelTyping) filter.set(search.value);
    } else {
      cancelTyping?.();
      cancelTyping = undefined;
    }
  });
  // The search's non-text part: the tree filtered by the selection (the selection bar's Show in structure tree), a
  // chip before the text; the search field shows while either is on.
  const selectionFilter = scene.selectionFilter;
  const filterOn = selectionFilter ? selectionFilter.on : observableState(false).asReadonly;
  const searchShown = derive([searchOpen.asReadonly, filterOn], (open, on) => open || on);
  const closeSearch = () => {
    search.set('');
    searchOpen.set(false);
    selectionFilter?.set(false);
  };
  return (
    <div class={classes.tab}>
      <div class={classes.actions}>
        {props.openCommands ? (
          <Button
            size="sm"
            variant="ghost"
            icon={<Icon name="shapes" />}
            onClick={() => props.openCommands?.('preset')}
          >
            <span class={classes.actionLabel}>Preset</span>
          </Button>
        ) : null}
        {props.openCommands ? (
          <Button
            size="sm"
            variant="ghost"
            icon={<Icon name="paintbrush" />}
            onClick={() => props.openCommands?.('styles')}
          >
            <span class={classes.actionLabel}>Styles</span>
          </Button>
        ) : null}
        <Button
          size="sm"
          variant="ghost"
          class={classes.searchButton}
          icon={<Icon name="search" />}
          title="Search the scene"
          pressed={searchShown}
          onClick={() => (searchShown.value ? closeSearch() : searchOpen.set(true))}
        />
        {props.look ? (
          <LookButton
            title="Scene defaults"
            class={classes.optionsButton}
            marked={derive(props.look.overrides, sceneOverridesText)}
            popover={(anchor, close) => <SceneLookPopover model={props.look!} anchor={anchor} onClose={close} />}
          />
        ) : null}
      </div>
      <Show when={searchShown}>
        {() => (
          <div
            class={classes.search}
            ref={(el) => queueMicrotask(() => el.querySelector('input')?.focus())}
            onKeyDown={(e: KeyboardEvent) => {
              // Backspace in an empty field takes the chip, as a token field does
              if (e.key !== 'Backspace' || !filterOn.value) return;
              if ((e.target as HTMLInputElement).value === '') selectionFilter?.set(false);
            }}
          >
            <Show when={filterOn}>
              {() => (
                <span class={classes.searchChip}>
                  <Icon name="mouse-pointer-click" size={12} />
                  In the selection
                  <button
                    type="button"
                    class={classes.searchChipRemove}
                    aria-label="Show every row, not only those in the selection"
                    data-tip="Show every row"
                    onClick={() => selectionFilter?.set(false)}
                  >
                    <Icon name="x" size={12} />
                  </button>
                </span>
              )}
            </Show>
            <TextInput
              size="sm"
              search
              placeholder={derive(filterOn, (on) =>
                on
                  ? 'Narrow it: an entity, a chain, a component…'
                  : 'Search the scene: an entity, a chain, a component…',
              )}
              aria-label="Search the scene"
              value={search.asReadonly}
              onInput={(v) => search.set(v)}
            />
          </div>
        )}
      </Show>
      {/* The tree scrolls itself (virtualized: scene-tree.tsx); this holds its place and takes its rows' picks. */}
      <div
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
            // Its row or its structure gone: its inspector goes too. One watch at a time, the last cleaned up.
            const stop = pick.gone.listen(() => {
              if (pickedRow.value === pick) pickedRow.set(undefined);
            });
            gone?.();
            gone = stop;
          };
          on(el, SceneRowPickEvent, onPick);
        }}
      >
        <SceneTree scene={scene} shell={shell} filter={filter.asReadonly} />
      </div>
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
              {(node) => (
                <Inspector
                  shell={shell}
                  inspection={scene.inspect(node)}
                  updates={scene.inspection}
                  onPick={(n) => scene.pick(n)}
                />
              )}
            </Show>
          )}
        >
          {(row) => <RowInspector row={row} />}
        </Show>
      </ScrollArea>
    </div>
  );
}

/** An extension's row in the inspector: its label, its line, its eye (turn it on or off), and its controls. */
function RowInspector(props: { row: SceneRowPick }) {
  const { row } = props;
  return (
    <div class={classes.inspection}>
      <header class={classes.inspectorHeader}>
        <span class={classes.inspectorTitle}>{row.label}</span>
        <span class={classes.inspectorKind}>{row.description ?? ''}</span>
        {row.shown && row.onToggle ? (
          <VisibilityToggle hidden={derive(row.shown, (s) => !s)} feature onToggle={() => row.onToggle?.()} />
        ) : null}
      </header>
      {row.body ? row.body() : <p class={classes.nothing}>Nothing to set here.</p>}
    </div>
  );
}

/**
 * A title that can be renamed (a component): its name and a pencil; editing, a field: Enter renames, Escape leaves it
 * as it was (and does not close the panel).
 */
function RenameTitle(props: { label: ReadonlyState<string>; onRename(label: string): void }) {
  const editing = observableState(false);
  return (
    <Show
      when={editing.asReadonly}
      fallback={() => (
        <span class={classes.renameTitle}>
          <span class={classes.inspectorTitle}>{props.label}</span>
          <Button
            variant="ghost"
            size="xs"
            icon={<Icon name="pencil" size={14} />}
            title="Rename"
            onClick={() => editing.set(true)}
          />
        </span>
      )}
    >
      {() => (
        <span
          class={classes.renameField}
          ref={(el: HTMLElement) =>
            queueMicrotask(() => {
              const input = el.querySelector('input');
              input?.focus();
              input?.select();
            })
          }
          onKeyDown={(e) => {
            if (e.key !== 'Escape') return;
            e.preventDefault();
            e.stopPropagation();
            editing.set(false);
          }}
        >
          <TextInput
            size="sm"
            aria-label="Name"
            value={props.label.value}
            onSubmit={(v) => {
              const name = v.trim();
              if (name && name !== props.label.value) props.onRename(name);
              editing.set(false);
            }}
          />
        </span>
      )}
    </Show>
  );
}

/**
 * The inspector's eye at its header's end, as the tree row's: Show or Hide what is drawn, Turn on or off a feature
 * (the unit cell, the streaming, an extension's row); marked while hidden or off, so the way back is where it was.
 */
function VisibilityToggle(props: { hidden: ReadonlyState<boolean>; feature: boolean; onToggle(): void }) {
  const [on, off] = props.feature ? ['Turn off', 'Turn on'] : ['Hide', 'Show'];
  return (
    <Button
      variant="ghost"
      size="xs"
      class={classes.inspectorEye}
      icon={
        <Show when={props.hidden} fallback={() => <Icon name="eye" size={14} />}>
          {() => <Icon name="eye-off" size={14} />}
        </Show>
      }
      title={derive(props.hidden, (h) => (h ? off : on))}
      pressed={props.hidden}
      onClick={() => props.onToggle()}
    />
  );
}

/** A node's kind in words: a structure's contents node is its Entities. */
const kindName = (kind: SceneNode['kind']) => (kind === 'contents' ? 'entities' : kind.replace('-', ' '));

/** What of an inspection is built anew and followed (`Inspector`). */
const Parts = [
  'facts',
  'contents',
  'quick',
  'levels',
  'transfer',
] as const satisfies readonly (keyof SceneInspection)[];

/** Rows that stand for a feature (on or off), not for something drawn (shown or hidden). */
const TurnsOn: ReadonlySet<string> = new Set(['extension', 'unit-cell']);

function Inspector(props: {
  shell: ShellModel;
  inspection: SceneInspection;
  /** The node's inspection built anew (`SceneModel.inspection`): its facts, contents, levels and transfer follow. */
  updates?: ReadonlyState<SceneInspection | undefined>;
  onPick(node: SceneNode): void;
}) {
  const { node, form, slot } = props.inspection;
  // What is built anew follows (the form, the header and the sections stay, a gesture on them going on): its parts by
  // name, each told only when it changes (view-models.md §3); one built for another node (on its way out) is not taken.
  const parts = keyedState<keyof SceneInspection, unknown>();
  const take = (i: SceneInspection) => {
    for (const k of Parts) parts.set(k, i[k]);
  };
  take(props.inspection);
  if (props.updates)
    effect(props.updates, (i) => {
      if (i?.node === node) take(i);
    });
  const part = <K extends (typeof Parts)[number]>(key: K) =>
    parts.key(key) as ReadonlyState<SceneInspection[K] | undefined>;
  return (
    <div class={classes.inspection}>
      <header class={classes.inspectorHeader}>
        {props.inspection.rename ? (
          <RenameTitle label={node.label} onRename={props.inspection.rename} />
        ) : (
          <span class={classes.inspectorTitle}>{node.label}</span>
        )}
        {/* The kind, unless the label already says it (Measurements, Components, Entities). */}
        <Show when={derive(node.label, (l) => l.toLowerCase() !== kindName(node.kind))}>
          {() => <span class={classes.inspectorKind}>{kindName(node.kind)}</span>}
        </Show>
        {/* Its look (the representations under it): their own values over the scene defaults, a dot while it has some. */}
        {props.inspection.look ? (
          <LookButton
            title="Look"
            marked={derive(props.inspection.look.own, (o) => (o.length ? 'not the scene defaults' : undefined))}
            popover={(anchor, close) => (
              <NodeLookPopover
                model={props.inspection.look!}
                title={derive(node.label, (l) => `Look of ${l}`)}
                anchor={anchor}
                onClose={close}
              />
            )}
          />
        ) : null}
        {/* No eye where nothing shows or hides as one (Entities, Other, an emptied component). */}
        {node.hasEye === false ? null : (
          <VisibilityToggle
            hidden={node.isHidden}
            feature={TurnsOn.has(node.kind)}
            onToggle={() => node.toggleVisibility()}
          />
        )}
      </header>
      <Show when={part('facts')}>
        {(facts) =>
          facts.length ? (
            <dl class={classes.facts}>
              {facts.map(([label, value]) => (
                <div class={classes.fact}>
                  <dt>{label}</dt>
                  <dd>{value}</dd>
                </div>
              ))}
            </dl>
          ) : null
        }
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
      {/* A component's or a shape's representations: each a row to open (its inspector), its eye beside it; not only
          behind the tree's chevron. */}
      {node.kind === 'component' || node.kind === 'shape' ? (
        <Representations node={node} onPick={props.onPick} />
      ) : null}
      {/* A structure's sequence: the band's style and letters, Add track, its chains (sequence.md §8.1). */}
      {props.inspection.sequence ? <SequenceInspector model={props.inspection.sequence} /> : null}
      {/* A structure's contents: its counts, its entities to open in the tree. */}
      <Show when={part('contents')}>{(contents) => <ContentsSummaryView summary={contents} />}</Show>
      {/* A representation's quick settings (type, color, opacity, size), above the whole form. */}
      <Show when={part('quick')}>{(quick) => <RepresentationQuick model={quick} />}</Show>
      {/* A map's levels (the iso values and colors people change most), above the whole form. */}
      <Show when={part('levels')}>{(levels) => levels.map((l) => <VolumeLevels model={l} />)}</Show>
      {/* A direct volume's transfer function (opacity and colors over the map's values). */}
      <Show when={part('transfer')}>{(transfer) => transfer.map((t) => <VolumeTransfer model={t} />)}</Show>
      {form ? <ParamsForm form={form} /> : null}
      <SectionList shell={props.shell} slot={slot} subject={node} flat />
      {/* Nothing to set (no form or actions of its own): what can be done with it is shown here, its row menu's
          entries; with none, say so rather than leave an empty space. */}
      {!form && !props.inspection.actions?.length && !props.inspection.contents && !props.inspection.sequence ? (
        (node.menu?.length ?? 0) > 0 || node.removable !== false ? (
          <NodeActions node={node} />
        ) : !props.inspection.facts?.length ? (
          <Show when={derive(sectionsFor(props.shell, slot, node), (has) => !has)}>
            {() => (
              <Show when={derive(node.isHidden, (h) => !h)}>
                {() => <p class={classes.nothing}>Nothing to set here.</p>}
              </Show>
            )}
          </Show>
        ) : null
      ) : null}
    </div>
  );
}

/**
 * A component's representations in its inspector: icon, name, what it is colored by, its eye and its ⋯ (the tree row's
 * menu: Remove, its type, …); a click opens it.
 */
function Representations(props: { node: SceneNode; onPick(node: SceneNode): void }) {
  const reprs = derive(props.node.children, (c) =>
    c.filter((n) => n.kind === 'representation' || n.kind === 'shape-representation'),
  );
  return (
    <Show when={derive(reprs, (r) => r.length > 0)}>
      {() => (
        <section class={classes.children} aria-label="Representations">
          <span class={classes.childrenTitle}>Representations</span>
          <For each={reprs} by={(n) => n.id}>
            {(n) => {
              const menuOpen = observableState(false);
              let menuButton: HTMLElement | undefined;
              // In the tree too: the component's row opens, the representation's is picked (and brought into view), so
              // where it is shows.
              const open = () => {
                if (!props.node.isExpanded.value) props.node.toggleExpanded();
                props.onPick(n);
              };
              return (
                <div class={classes.child} classList={{ [classes.childHidden]: n.isHidden }}>
                  <button type="button" class={classes.childOpen} onClick={open}>
                    <Icon name={n.icon} size={15} class={classes.childIcon} />
                    <span class={classes.childLabel}>{n.label}</span>
                    <span class={classes.childDescription}>{n.description}</span>
                  </button>
                  <button
                    type="button"
                    class={classes.childEye}
                    aria-label={derive(n.isHidden, (h) => (h ? `Show ${n.label.value}` : `Hide ${n.label.value}`))}
                    data-tip={derive(n.isHidden, (h) => (h ? 'Show' : 'Hide'))}
                    onClick={() => n.toggleVisibility()}
                  >
                    <Show when={n.isHidden} fallback={() => <Icon name="eye" size={14} />}>
                      {() => <Icon name="eye-off" size={14} />}
                    </Show>
                  </button>
                  <button
                    type="button"
                    class={classes.childEye}
                    aria-label={derive(n.label, (l) => `${l}: more`)}
                    aria-haspopup="menu"
                    aria-expanded={menuOpen}
                    data-tip="More"
                    ref={(e: HTMLElement) => {
                      menuButton = e;
                    }}
                    onClick={() => menuOpen.set(!menuOpen.value)}
                  >
                    <Icon name="ellipsis" size={14} />
                  </button>
                  {/* At the row's edge, always shown: the row opens (as its label does; one tab stop, the label's). */}
                  <button type="button" class={classes.childChevron} tabIndex={-1} aria-hidden="true" onClick={open}>
                    <Icon name="chevron-right" size={14} />
                  </button>
                  <Show when={menuOpen}>
                    {() => (
                      <RowMenu
                        items={menuOf(n)}
                        anchor={menuButton!}
                        onDone={(refocus) => {
                          menuOpen.set(false);
                          if (refocus) menuButton?.focus();
                        }}
                      />
                    )}
                  </Show>
                </div>
              );
            }}
          </For>
        </section>
      )}
    </Show>
  );
}

/** Whether any contributed section shows for a node (`inspector:<kind>` sections with their `when`). */
function sectionsFor(shell: ShellModel, slot: string, subject: unknown) {
  return derive(shell.sections(slot), (all) => all.some((s) => !s.when || s.when({ shell, subject })));
}
