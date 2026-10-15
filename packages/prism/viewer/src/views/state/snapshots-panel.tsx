/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { Button } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ParamsForm } from '@molstar/prism-ui';
import { dragToReorder } from '@molstar/prism-ui/components/reorder/drag-to-reorder';
import { ScrollArea } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { SnapshotEntry, SnapshotsModel } from '../../model/snapshots.js';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { classes } from './_generated/snapshots-panel.module.css.js';

/**
 * Save state: the Structure tab's shape for snapshots. A header: Add, and what a new snapshot includes behind the
 * sliders (saving and opening a file is Export › Session & state's, and Home's Open); the snapshots as rows (a click
 * applies one; its ⋯ reorders, replaces, removes), and under them, at a fixed height, the current snapshot's own settings.
 */
export function SnapshotsPanel(props: { model: SnapshotsModel }) {
  const { model } = props;
  const open = observableState(false);
  const current = model.current.pipe(distinctUntilChanged());
  let list: HTMLElement | undefined;
  const drag = (id: string) => (e: PointerEvent) => dragToReorder(e, list, id, (index) => model.moveTo(id, index));
  return (
    <div class={classes.panel}>
      <div class={classes.actions}>
        {/* As the Structure tab's header: quiet buttons, an icon and a word each. */}
        <Button size="sm" variant="ghost" icon={<Icon name="plus" />} onClick={() => model.add()}>
          <span class={classes.actionLabel}>Add</span>
        </Button>
        <Button
          size="sm"
          variant="ghost"
          class={classes.optionsButton}
          icon={<Icon name="settings" />}
          title="Snapshot options"
          pressed={open.asReadonly}
          onClick={() => open.set(!open.value)}
        />
      </div>
      <Show when={open.asReadonly} reveal>
        {() => (
          <ScrollArea
            class={classes.block}
            contentClass={classes.blockContent}
            role="region"
            aria-label="Snapshot options"
          >
            <p class={classes.blockCaption}>What a new snapshot includes</p>
            <ParamsForm form={model.saveOptions} />
          </ScrollArea>
        )}
      </Show>
      <ScrollArea
        class={classes.list}
        contentClass={classes.listContent}
        role="listbox"
        aria-label="Snapshots"
        ref={(el: HTMLElement) => (list = el)}
      >
        <For
          each={model.entries}
          key={(e) => e.id}
          fallback={() => (
            <p class={classes.empty}>
              No snapshots yet. <strong>Add</strong> saves the current view; saved ones play back in order.
            </p>
          )}
        >
          {(entry, index) => (
            <SnapshotRow
              model={model}
              entry={entry}
              index={index}
              onGrab={drag(entry.id)}
              current={current.pipe(
                map((c) => c === entry.id),
                distinctUntilChanged(),
              )}
            />
          )}
        </For>
      </ScrollArea>
      <ScrollArea class={classes.inspector} contentClass={classes.inspectorContent} role="region" aria-label="Snapshot">
        <Show
          when={model.entries.pipe(
            map((entries) => entries.find((e) => e.id === model.current.value)),
            distinctUntilChanged(),
          )}
          fallback={() => <p class={classes.inspectorEmpty}>Apply a snapshot to edit its name, key and description.</p>}
        >
          {(entry) => <SnapshotEditor model={model} entry={entry} />}
        </Show>
      </ScrollArea>
    </div>
  );
}

/** A snapshot: its position, thumbnail, name, key and the description's first line; a click applies it. */
function SnapshotRow(props: {
  model: SnapshotsModel;
  entry: SnapshotEntry;
  index: ReadonlyState<number>;
  current: Observable<boolean>;
  /** Pointer down on the handle: starts dragging the row. */
  onGrab: (e: PointerEvent) => void;
}) {
  const { model, entry } = props;
  const menu = observableState(false);
  const firstLine = entry.description
    ?.split('\n')
    .map((l) => l.replace(/[#*`>-]/g, '').trim())
    .find((l) => l);
  return (
    <div class={classes.item} data-id={entry.id}>
      <div
        class={classes.row}
        role="option"
        tabIndex={0}
        aria-selected={props.current}
        onClick={() => model.apply(entry.id)}
        onKeyDown={(e) => {
          if (e.key === 'Enter' || e.key === ' ') {
            e.preventDefault();
            model.apply(entry.id);
          }
        }}
      >
        {/* The position; on hover, the handle that drags the row (Move up/down stay in its ⋯). */}
        <span class={classes.slot}>
          <span class={classes.index}>{props.index.pipe(map((i) => `${i + 1}`))}</span>
          <span
            class={classes.grip}
            aria-hidden="true"
            data-tip="Drag to reorder"
            onPointerDown={props.onGrab}
            onClick={(e) => e.stopPropagation()}
          >
            <Icon name="grip-vertical" size={14} />
          </span>
        </span>
        {entry.image ? <span class={classes.thumb} style={{ background: entry.image }} /> : null}
        <span class={classes.name}>{entry.name}</span>
        {entry.key ? <span class={classes.key}>#{entry.key}</span> : null}
        <span class={classes.description}>{firstLine ?? ''}</span>
        <span class={classes.rowActions}>
          <button
            type="button"
            class={classes.action}
            data-tip="Retake from the current view"
            aria-label="Retake from the current view"
            onClick={(e) => {
              e.stopPropagation();
              model.replace(entry.id);
            }}
          >
            <Icon name="camera" size={14} />
          </button>
          <button
            type="button"
            class={classes.action}
            aria-label="More"
            aria-expanded={menu}
            onClick={(e) => {
              e.stopPropagation();
              menu.set(!menu.value);
            }}
          >
            <Icon name="ellipsis" size={15} />
          </button>
        </span>
      </div>
      <Show when={menu}>
        {() => (
          <div class={classes.menu} role="menu">
            <MenuItem icon="arrow-up" label="Move up" run={() => model.move(entry.id, -1)} />
            <MenuItem icon="arrow-down" label="Move down" run={() => model.move(entry.id, 1)} />
            <MenuItem
              icon="trash-2"
              label="Remove"
              danger
              run={() => {
                menu.set(false);
                model.remove(entry.id);
              }}
            />
          </div>
        )}
      </Show>
    </div>
  );
}

function MenuItem(props: {
  icon: 'arrow-up' | 'arrow-down' | 'trash-2';
  label: string;
  danger?: boolean;
  run: () => void;
}) {
  return (
    <button
      type="button"
      role="menuitem"
      class={classes.menuItem}
      data-danger={props.danger ? '' : undefined}
      onClick={props.run}
    >
      <span>{props.label}</span>
    </button>
  );
}

/** The current snapshot's own settings: name, key, description (markdown), and how it plays. */
function SnapshotEditor(props: { model: SnapshotsModel; entry: SnapshotEntry }) {
  const { model, entry } = props;
  const options = model.entryOptions(entry.id);
  return (
    <div class={classes.inspection}>
      <header class={classes.inspectorHeader}>
        <span class={classes.inspectorTitle}>{entry.name}</span>
        <span class={classes.inspectorKind}>{entry.timestamp ?? 'snapshot'}</span>
      </header>
      <div class={classes.fields}>
        <label class={classes.field}>
          <span class={classes.label}>Name</span>
          <TextInput size="sm" value={entry.name} onSubmit={(v) => model.update(entry.id, { name: v })} />
        </label>
        <label class={classes.field}>
          <span class={classes.label}>Key</span>
          <TextInput
            size="sm"
            value={entry.key ?? ''}
            placeholder="Optional, for #links"
            onSubmit={(v) => model.update(entry.id, { key: v || undefined })}
          />
        </label>
        <label class={classes.descriptionField}>
          <span class={classes.label}>Description</span>
          <textarea
            class={classes.textarea}
            rows={4}
            placeholder="Markdown · [title](#key) links to a snapshot"
            value={entry.description ?? ''}
            onChange={(e) => model.update(entry.id, { description: e.currentTarget.value || undefined })}
          />
        </label>
      </div>
      {options ? <ParamsForm form={options} /> : null}
    </div>
  );
}
