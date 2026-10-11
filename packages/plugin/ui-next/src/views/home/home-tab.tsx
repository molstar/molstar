/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, of } from 'rxjs';
import { formatChord } from '../../commands/keys.js';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { For, Show } from '../../dom/control-flow.js';
import type { HomeAction, HomeModel } from '../../model/home.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/home-tab.module.css.js';

export interface HomeTabProps {
  home: HomeModel;
  /** Opens the command palette (the hint at the bottom). */
  onOpenPalette?: () => void;
}

/**
 * The Home tab: open files (a drop area), the state actions applicable to the root as forms (download structure,
 * download density, load trajectory, …; the first one open), recently opened entries, and a pointer to the palette.
 */
export function HomeTab(props: HomeTabProps) {
  const { home } = props;
  return (
    <div class={classes.home}>
      <OpenFiles home={home} />
      <section class={classes.group} aria-label="Actions">
        <For each={home.actions} key={(a) => a.id}>
          {(action, index) => <ActionBlock action={action} open={index.value === 0} />}
        </For>
      </section>
      <Show when={home.recent.pipe(map((r) => r.length > 0))}>
        {() => (
          <section class={classes.group} aria-label="Recent">
            <h3 class={classes.caption}>Recent</h3>
            <For each={home.recent} key={(r) => r.id}>
              {(entry) => (
                <button type="button" class={classes.recent} onClick={() => entry.open()}>
                  <Icon name={entry.icon} size={16} />
                  <span class={classes.recentLabel}>{entry.label}</span>
                  <span class={classes.recentDescription}>{entry.description}</span>
                </button>
              )}
            </For>
          </section>
        )}
      </Show>
      <button type="button" class={classes.hint} onClick={() => props.onOpenPalette?.()}>
        <kbd>{formatChord({ key: 'k', mod: true })}</kbd>
        <span>Search every command</span>
      </button>
    </div>
  );
}

/** Files by picking or dropping: the whole area is the target. */
function OpenFiles(props: { home: HomeModel }) {
  const { home } = props;
  const dragging = observableState(false);
  const input = (
    <input
      type="file"
      multiple
      hidden
      onChange={(e) => {
        const files = e.currentTarget.files;
        if (files?.length) home.openFiles([...files]);
        e.currentTarget.value = '';
      }}
    />
  ) as HTMLInputElement;
  return (
    <button
      type="button"
      class={classes.drop}
      data-dragging={dragging.asReadonly.pipe(map((d) => (d ? '' : undefined)))}
      onClick={() => input.click()}
      onDragOver={(e) => {
        e.preventDefault();
        dragging.set(true);
      }}
      onDragLeave={() => dragging.set(false)}
      onDrop={(e) => {
        e.preventDefault();
        dragging.set(false);
        const files = e.dataTransfer?.files;
        if (files?.length) home.openFiles([...files]);
      }}
    >
      <Icon name="folder-open" size={20} />
      <span class={classes.dropText}>
        <span class={classes.dropTitle}>Open files</span>
        <span class={classes.dropHint}>or drop them here · {home.formats.join(', ')}</span>
      </span>
      {input}
    </button>
  );
}

/** One action: a header that folds, its form, and the apply button. */
function ActionBlock(props: { action: HomeAction; open: boolean }) {
  const { action } = props;
  const open = observableState(props.open);
  const busy = action.busy ?? of(false);
  return (
    <div class={classes.action} data-open={open.asReadonly.pipe(map((o) => (o ? '' : undefined)))}>
      <button type="button" class={classes.actionHeader} aria-expanded={open} onClick={() => open.set(!open.value)}>
        <Icon name={action.icon} size={16} />
        <span class={classes.actionTitle}>{action.title}</span>
        <span class={classes.actionDescription}>{action.description}</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      <Show when={open}>
        {() => (
          <div class={classes.actionBody}>
            {action.form ? <ParamsForm form={action.form} /> : null}
            <div class={classes.actionFooter}>
              <Button
                variant="primary"
                block
                disabled={busy}
                icon={
                  <Show when={busy} fallback={() => <Icon name={action.icon} />}>
                    {() => <Icon name="loader-circle" class={classes.spin} />}
                  </Show>
                }
                onClick={() => action.apply()}
              >
                {action.applyLabel ?? 'Apply'}
              </Button>
            </div>
          </div>
        )}
      </Show>
    </div>
  );
}
