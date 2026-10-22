/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, onActive, onCleanup } from '@molstar/prism-ui';
import { observableState } from '@molstar/prism-ui/core';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { type Observable } from 'rxjs';
import { formatBytes, type StorageEntry, type StorageModel, type StorageSourceState } from '../storage/model.js';
import { classes } from './_generated/storage.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// The Storage view (journal.md §7, README decision 32): what the viewer keeps in this browser and where (kept records,
// preferences, the light restore, the assistant's stored keys), each with its size and Forget, and Forget everything
// (asked first). A folded section of the Settings tab; it lists only while open and shown (`watch` under `onActive`):
// hidden, it reads nothing.

export function StorageSection(props: { storage: StorageModel }) {
  const open = observableState(false);
  return (
    <section class={classes.section} aria-label="Storage">
      <div class={classes.header}>
        <button
          type="button"
          class={classes.toggle}
          aria-expanded={open.asReadonly}
          onClick={() => open.set(!open.value)}
        >
          <span class={classes.title}>Storage</span>
          <Icon name="chevron-down" size={14} class={classes.chevron} />
        </button>
      </div>
      <Show when={open.asReadonly} reveal>
        {() => <StorageBody storage={props.storage} />}
      </Show>
    </section>
  );
}

function StorageBody(props: { storage: StorageModel }) {
  const { storage } = props;
  // Lists while shown; a paused (hidden) body stops watching and reads nothing.
  let stop: (() => void) | undefined;
  onActive((active) => {
    if (active) stop ??= storage.watch();
    else {
      stop?.();
      stop = undefined;
    }
  });
  onCleanup(() => stop?.());
  const busy = observableState(false);
  const confirming = observableState(false);
  const run = (work: () => Promise<void>) => {
    if (busy.value) return;
    busy.set(true);
    void work().finally(() => busy.set(false));
  };
  const empty = derive(storage.sources, (sources) => sources.every((s) => s.listed && !s.entries.length && !s.error));
  const where = storage.origin ? `Kept in this browser for ${storage.origin}.` : 'Kept in this browser.';
  return (
    <div class={classes.body}>
      <p class={classes.where}>{where} Nothing here leaves it.</p>
      <For each={storage.sources} by={(s) => s.id}>
        {(source) => (
          <Source source={source} busy={busy.asReadonly} onForget={(id) => run(() => storage.forget(source.id, id))} />
        )}
      </For>
      <div class={classes.footer}>
        <Show
          when={confirming.asReadonly}
          fallback={() => (
            <Button
              size="xs"
              variant="ghost"
              icon={<Icon name="trash-2" size={14} />}
              title="Forget everything listed here"
              disabled={busy.asReadonly}
              onClick={() => confirming.set(true)}
            >
              Forget everything
            </Button>
          )}
        >
          {() => (
            <div class={classes.confirm} role="group" aria-label="Forget everything?">
              <span class={classes.confirmText}>Forget everything listed here? It cannot be undone.</span>
              <Button
                size="xs"
                variant="danger"
                disabled={busy.asReadonly}
                onClick={() => {
                  confirming.set(false);
                  run(() => storage.forgetEverything());
                }}
              >
                Forget everything
              </Button>
              <Button size="xs" variant="ghost" onClick={() => confirming.set(false)}>
                Cancel
              </Button>
            </div>
          )}
        </Show>
        <Show when={empty}>{() => <span class={classes.caption}>Nothing is kept.</span>}</Show>
      </div>
    </div>
  );
}

function Source(props: { source: StorageSourceState; busy: Observable<boolean>; onForget(entryId: string): void }) {
  const { source } = props;
  const meta = [source.where, source.bytes !== undefined ? formatBytes(source.bytes) : undefined]
    .filter(Boolean)
    .join(' · ');
  return (
    <div class={classes.source} data-source={source.id}>
      <div class={classes.sourceHead}>
        <span class={classes.sourceTitle}>{source.title}</span>
        <span class={classes.meta}>{meta}</span>
      </div>
      {source.error ? <p class={classes.error}>{source.error}</p> : null}
      {source.listed && !source.entries.length && !source.error ? <p class={classes.none}>Nothing kept.</p> : null}
      {source.entries.map((entry) => (
        <Entry entry={entry} busy={props.busy} onForget={() => props.onForget(entry.id)} />
      ))}
    </div>
  );
}

function Entry(props: { entry: StorageEntry; busy: Observable<boolean>; onForget(): void }) {
  const { entry } = props;
  return (
    <div class={classes.entry} data-stale={entry.stale ? '' : undefined}>
      <span class={classes.label}>{entry.label}</span>
      <span class={classes.size}>{entry.bytes !== undefined ? formatBytes(entry.bytes) : ''}</span>
      <Button
        size="xs"
        variant="ghost"
        icon={<Icon name="trash-2" size={14} />}
        title={`Forget ${entry.label}`}
        disabled={props.busy}
        onClick={() => props.onForget()}
      />
      {entry.detail ? <span class={classes.detail}>{entry.detail}</span> : null}
    </div>
  );
}
