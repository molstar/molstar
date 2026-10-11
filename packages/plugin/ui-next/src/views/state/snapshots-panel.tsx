/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
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
  const drag = (id: string) => (e: PointerEvent) => startDrag(e, list, id, (index) => model.moveTo(id, index));
  return (
    <div class={classes.panel}>
      <div class={classes.actions}>
        <Button size="sm" variant="ghost" icon={<Icon name="plus" />} onClick={() => model.add()}>
          <span class={classes.actionLabel}>Snapshot</span>
        </Button>
        <Button
          size="sm"
          variant="ghost"
          icon={<Icon name="download" />}
          pressed={open.asReadonly.pipe(map((o) => o === 'save'))}
          onClick={() => toggle('save')}
        >
          <span class={classes.actionLabel}>Save</span>
        </Button>
        <Button size="sm" variant="ghost" icon={<Icon name="folder-open" />} onClick={() => file.click()}>
          <span class={classes.actionLabel}>Open</span>
        </Button>
        {file}
        <Button
          size="sm"
          variant="ghost"
          class={classes.optionsButton}
          icon={<Icon name="sliders-horizontal" />}
          title="Save options"
          pressed={open.asReadonly.pipe(map((o) => o === 'options'))}
          onClick={() => toggle('options')}
        />
      </div>
      <Show when={open}>
        {(block) =>
          block === 'save' ? (
            <section class={classes.block} aria-label="Save to a file">
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
            </section>
          ) : (
            <section class={classes.block} aria-label="Save options">
              <p class={classes.blockCaption}>What a new snapshot includes</p>
              <ParamsForm form={model.saveOptions} />
            </section>
          )
        }
      </Show>
      <div class={classes.list} role="listbox" aria-label="Snapshots" ref={(el: HTMLElement) => (list = el)}>
        <For
          each={model.entries}
          key={(e) => e.id}
          fallback={() => (
            <p class={classes.empty}>
              No snapshots yet. <strong>Snapshot</strong> saves the current view; saved ones play back in order.
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
      </div>
      <section class={classes.inspector} aria-label="Snapshot">
        <Show
          when={model.entries.pipe(
            map((entries) => entries.find((e) => e.id === model.current.value)),
            distinctUntilChanged(),
          )}
          fallback={() => <p class={classes.inspectorEmpty}>Apply a snapshot to edit its name, key and description.</p>}
        >
          {(entry) => <SnapshotEditor model={model} entry={entry} />}
        </Show>
      </section>
    </div>
  );
}

function SaveChoice(props: { icon: 'download'; title: string; hint: string; onClick: () => void }) {
  return (
    <button type="button" class={classes.choice} onClick={props.onClick}>
      <Icon name={props.icon} size={15} />
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
      <Icon name={props.icon} size={14} />
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
        <Icon name="images" size={15} />
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

/**
 * Drags a row by its handle: the row under the pointer shows where it would land (a line above or below it), release
 * moves it there. Pointer captured on the handle, so leaving the list keeps the drag.
 */
function startDrag(e: PointerEvent, list: HTMLElement | undefined, id: string, moveTo: (index: number) => void) {
  if (!list || e.button !== 0) return;
  e.preventDefault();
  e.stopPropagation();
  const handle = e.currentTarget as HTMLElement;
  handle.setPointerCapture(e.pointerId);
  const items = [...list.querySelectorAll<HTMLElement>(':scope > [data-id]')];
  const self = items.find((el) => el.dataset.id === id);
  self?.setAttribute('data-dragging', '');
  let target = -1;
  const mark = (index: number) => {
    target = index;
    for (let i = 0; i < items.length; i++) {
      const el = items[i];
      if (i === index) el.setAttribute('data-drop', 'before');
      else if (index === items.length && i === items.length - 1) el.setAttribute('data-drop', 'after');
      else el.removeAttribute('data-drop');
    }
  };
  const onMove = (ev: PointerEvent) => {
    let index = items.length;
    for (let i = 0; i < items.length; i++) {
      const r = items[i].getBoundingClientRect();
      if (ev.clientY < r.top + r.height / 2) {
        index = i;
        break;
      }
    }
    mark(index);
  };
  const onUp = () => {
    handle.removeEventListener('pointermove', onMove);
    self?.removeAttribute('data-dragging');
    const from = items.indexOf(self as HTMLElement);
    const to = target;
    mark(-1);
    // Dropping above or below itself leaves it where it is.
    if (to >= 0 && to !== from && to !== from + 1) moveTo(to > from ? to - 1 : to);
  };
  handle.addEventListener('pointermove', onMove);
  handle.addEventListener('pointerup', onUp, { once: true });
  handle.addEventListener('pointercancel', onUp, { once: true });
}
