/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { dragToReorder } from '../../components/reorder/drag-to-reorder.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import type { SnapshotEntry, SnapshotsModel } from '../../model/snapshots.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { classes } from './_generated/snapshots-panel.module.css.js';

/**
 * Save state: the Structure tab's shape for snapshots. A header of actions (a new snapshot, save to a file, open one; the
 * save options behind the sliders), the snapshots as rows (a click applies one; its ⋯ reorders, replaces, removes), and
 * under them, at a fixed height, the current snapshot's own settings.
 */
export function SnapshotsPanel(props: { model: SnapshotsModel }) {
  const { model } = props;
  const open = observableState<'save' | 'options' | undefined>(undefined);
  const toggle = (block: 'save' | 'options') => open.set(open.value === block ? undefined : block);
  const file = (
    <input
      type="file"
      accept=".molx,.molj,.zip,.json"
      hidden
      onChange={(e) => {
        const files = e.currentTarget.files;
        if (files?.length) model.open([...files]);
        e.currentTarget.value = '';
      }}
    />
  ) as HTMLInputElement;
  const current = model.current.pipe(distinctUntilChanged());
  let list: HTMLElement | undefined;
  const drag = (id: string) => (e: PointerEvent) => dragToReorder(e, list, id, (index) => model.moveTo(id, index));
  return (
    <div class={classes.panel}>
      <div class={classes.actions}>
        {/* Buttons that look like buttons (the control surface); Add keeps its +, the others are words. */}
        <Button size="sm" class={classes.add} icon={<Icon name="plus" />} onClick={() => model.add()}>
          <span class={classes.actionLabel}>Add</span>
        </Button>
        <Button
          size="sm"
          icon={<Icon name="download" />}
          pressed={open.asReadonly.pipe(map((o) => o === 'save'))}
          onClick={() => toggle('save')}
        >
          <span class={classes.actionLabel}>Save</span>
        </Button>
        <Button size="sm" icon={<Icon name="folder-open" />} onClick={() => file.click()}>
          <span class={classes.actionLabel}>Open</span>
        </Button>
        {file}
        <Button
          size="sm"
          variant="ghost"
          class={classes.optionsButton}
          icon={<Icon name="settings" />}
          title="Save options"
          pressed={open.asReadonly.pipe(map((o) => o === 'options'))}
          onClick={() => toggle('options')}
        />
      </div>
      <Show when={open} reveal>
        {(block) =>
          block === 'save' ? (
            <ScrollArea
              class={classes.block}
              contentClass={classes.blockContent}
              role="region"
              aria-label="Save to a file"
            >
              <SaveChoice
                icon="download"
                title="Session"
                hint=".molx · the state with its data"
                onClick={() => {
                  open.set(undefined);
                  model.download('session');
                }}
              />
              <SaveChoice
                icon="download"
                title="State"
                hint=".molj · the description; data loads from its sources"
                onClick={() => {
                  open.set(undefined);
                  model.download('state');
                }}
              />
            </ScrollArea>
          ) : (
            <ScrollArea
              class={classes.block}
              contentClass={classes.blockContent}
              role="region"
              aria-label="Save options"
            >
              <p class={classes.blockCaption}>What a new snapshot includes</p>
              <ParamsForm form={model.saveOptions} />
            </ScrollArea>
          )
        }
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

function SaveChoice(props: { icon: 'download'; title: string; hint: string; onClick: () => void }) {
  return (
    <button type="button" class={classes.choice} onClick={props.onClick}>
      <span class={classes.choiceTitle}>{props.title}</span>
      <span class={classes.choiceHint}>{props.hint}</span>
    </button>
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
