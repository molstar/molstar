/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, MenuPopover } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ParamsForm } from '@molstar/prism-ui';
import { ScrollArea } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { keyedSelector } from '@molstar/prism-ui/reactive/keyed-selector';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { SnapshotEntry, SnapshotsModel } from '../../model/snapshots.js';
import { classes } from './_generated/snapshots-panel.module.css.js';
import { reorderRows } from '@molstar/prism-ui/components/reorder/drag-to-reorder';

/**
 * Save state: the Structure tab's shape for snapshots. A header: Add, and what a new snapshot includes behind the
 * sliders (saving and opening a file is Export › Session & state's, and Home's Open); the snapshots as rows (a click
 * applies one; its ⋯ reorders, replaces, removes), and under them, at a fixed height, the current snapshot's own settings.
 */
export function SnapshotsPanel(props: {
  model: SnapshotsModel;
  /** At the header's end, before the options (Play): what the tab adds to it. */
  actions?: () => Child;
}) {
  const { model } = props;
  const open = observableState(false);
  // Each row listens to whether it is the current one, told only when that changes (two rows a change, not all).
  const current = keyedSelector(model.current);
  let list: HTMLElement | undefined;
  // A list needs items: with none, it is the plain text saying so (accessibility.md).
  const some = derive(model.entries, (entries) => entries.length > 0);
  const panel = (
    <div class={classes.panel}>
      <div class={classes.actions}>
        {/* As the Structure tab's header: quiet buttons, an icon and a word each. */}
        <Button size="sm" variant="ghost" icon={<Icon name="plus" />} onClick={() => model.add()}>
          <span class={classes.actionLabel}>Add</span>
        </Button>
        <span class={classes.trailing}>{props.actions?.()}</span>
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
        role={derive(some, (some) => (some ? 'list' : undefined))}
        aria-label={derive(some, (some) => (some ? 'Snapshots' : undefined))}
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
          {(entry, index) => <SnapshotRow model={model} entry={entry} index={index} current={current.key(entry.id)} />}
        </For>
      </ScrollArea>
      <ScrollArea class={classes.inspector} contentClass={classes.inspectorContent} role="region" aria-label="Snapshot">
        <Show
          when={derive([model.entries, model.current], (entries, id) => entries.find((e) => e.id === id))}
          fallback={() => <p class={classes.inspectorEmpty}>Apply a snapshot to edit its name, key and description.</p>}
        >
          {(entry) => <SnapshotEditor model={model} entry={entry} />}
        </Show>
      </ScrollArea>
    </div>
  );
  if (list) reorderRows(list, `.${classes.grip}`, (id, index) => model.moveTo(id, index));
  return panel;
}

/** A snapshot: its position, thumbnail, name, key and the description's first line; a click applies it. */
function SnapshotRow(props: {
  model: SnapshotsModel;
  entry: SnapshotEntry;
  index: ReadonlyState<number>;
  current: ReadonlyState<boolean>;
}) {
  const { model, entry } = props;
  const menu = observableState(false);
  let more: HTMLElement | undefined;
  const firstLine = entry.description
    ?.split('\n')
    .map((l) => l.replace(/[#*`>-]/g, '').trim())
    .find((l) => l);
  return (
    // A list item: the row (applies the snapshot) and, beside it, its buttons (an interactive role holds no buttons).
    <div
      class={classes.item}
      role="listitem"
      data-id={entry.id}
      data-current={derive(props.current, (c) => (c ? '' : undefined))}
    >
      <div
        class={classes.row}
        role="button"
        tabIndex={0}
        aria-current={derive(props.current, (c) => (c ? 'true' : undefined))}
        aria-label={`Apply ${entry.name}`}
        aria-keyshortcuts="Alt+ArrowUp Alt+ArrowDown"
        onClick={() => model.apply(entry.id)}
        onKeyDown={(e) => {
          if (e.key === 'Enter' || e.key === ' ') {
            e.preventDefault();
            model.apply(entry.id);
          } else if (e.altKey && (e.key === 'ArrowUp' || e.key === 'ArrowDown') && e.target === e.currentTarget) {
            // The order by keyboard (the handle is for the pointer): the row moves and keeps the focus.
            e.preventDefault();
            const row = e.currentTarget as HTMLElement;
            model.move(entry.id, e.key === 'ArrowUp' ? -1 : 1);
            queueMicrotask(() => {
              if (row.isConnected && document.activeElement !== row) row.focus();
            });
          }
        }}
      >
        {/* The position; on hover, the handle that drags the row (by hand; Alt+Up and Down on the row move it by key). */}
        <span class={classes.slot}>
          <span class={classes.index}>{derive(props.index, (i) => `${i + 1}`)}</span>
          <span
            class={classes.grip}
            aria-hidden="true"
            data-tip="Drag to reorder (Alt+Up or Down on the row)"
            onClick={(e) => e.stopPropagation()}
          >
            <Icon name="grip-vertical" size={14} />
          </span>
        </span>
        {entry.image ? <span class={classes.thumb} style={{ background: entry.image }} /> : null}
        <span class={classes.name}>{entry.name}</span>
        {entry.key ? <span class={classes.key}>#{entry.key}</span> : null}
        <span class={classes.description}>{firstLine ?? ''}</span>
      </div>
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
          aria-haspopup="menu"
          aria-expanded={menu}
          ref={(el: HTMLElement) => (more = el)}
          onClick={(e) => {
            e.stopPropagation();
            menu.set(!menu.value);
          }}
        >
          <Icon name="ellipsis" size={15} />
        </button>
      </span>
      {/* A popover from ⋯ (as the tree's rows): Remove; the order is the drag handle's and Alt+Up and Down's (no Move up and down). */}
      <Show when={menu}>
        {() => (
          <MenuPopover
            anchor={more!}
            aria-label={`${entry.name}: actions`}
            onClose={(refocus) => {
              menu.set(false);
              if (refocus) more?.focus();
            }}
          >
            <MenuItem
              icon="trash-2"
              label="Remove"
              danger
              run={() => {
                menu.set(false);
                model.remove(entry.id);
              }}
            />
          </MenuPopover>
        )}
      </Show>
    </div>
  );
}

function MenuItem(props: { icon: 'trash-2'; label: string; danger?: boolean; run: () => void }) {
  return (
    <button
      type="button"
      role="menuitem"
      class={classes.menuItem}
      data-danger={props.danger ? '' : undefined}
      onClick={props.run}
    >
      {/* The column the tree's menus keep for a chevron: the label in the second. */}
      <span />
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
