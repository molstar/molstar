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
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import type { HomeAction, HomeModel, RemoteStatesModel } from '../../model/home.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/home-tab.module.css.js';

export interface HomeTabProps {
  home: HomeModel;
  /** Opens the command menu (the hint at the bottom). */
  onOpenCommandMenu?: () => void;
}

/**
 * The Home tab: the state actions applicable to the root as forms (download structure, open files, download density,
 * load trajectory, …; the first one open), recently opened entries, and a pointer to the command menu. Files dropped
 * anywhere on the plugin open too (the shell's drop overlay).
 */
export function HomeTab(props: HomeTabProps) {
  const { home } = props;
  return (
    <div class={classes.home}>
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
      {home.remoteStates ? <RemoteStates model={home.remoteStates} /> : null}
      <button type="button" class={classes.hint} onClick={() => props.onOpenCommandMenu?.()}>
        <kbd>{formatChord({ key: 'k', mod: true })}</kbd>
        <span>Search every command</span>
      </button>
    </div>
  );
}

/** One action: a header that folds, its form, and the apply button (Home's actions, the state tree's). */
export function ActionBlock(props: { action: HomeAction; open: boolean }) {
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
                size="sm"
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

/** Sessions on a state server: the list (newest first), the server, refresh. */
function RemoteStates(props: { model: RemoteStatesModel }) {
  const { model } = props;
  return (
    <section class={classes.group} aria-label="Remote states">
      <div class={classes.captionRow}>
        <h3 class={classes.caption}>Remote states</h3>
        <Button
          size="xs"
          variant="ghost"
          icon={
            <Show when={model.loading} fallback={() => <Icon name="rotate-ccw" />}>
              {() => <Icon name="loader-circle" class={classes.spin} />}
            </Show>
          }
          title="Refresh"
          onClick={() => model.refresh()}
        />
      </div>
      <For each={model.states} key={(r) => r.id} fallback={() => <p class={classes.none}>No states on this server.</p>}>
        {(state) => (
          <button type="button" class={classes.recent} onClick={() => model.open(state.id)}>
            <Icon name="cloud" size={16} />
            <span class={classes.recentLabel}>{state.label}</span>
            <span class={classes.recentDescription}>{state.date}</span>
          </button>
        )}
      </For>
      <label class={classes.server}>
        <span>Server</span>
        <TextInput
          size="sm"
          aria-label="State server URL"
          value={model.serverUrl}
          onSubmit={(v) => model.setServerUrl(v)}
        />
      </label>
    </section>
  );
}
