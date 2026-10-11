/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { formatChord } from '../../commands/keys.js';
import { ActionFooter } from '../../components/action-footer/action-footer.js';
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
 * load trajectory, …; the first one open), remote states when on, and a pointer to the command menu. Files dropped
 * anywhere on the plugin open too (the shell's drop overlay).
 */
export function HomeTab(props: HomeTabProps) {
  const { home } = props;
  return (
    <div class={classes.home}>
      {/* What scrolls; the hint below stays in view. */}
      <div class={classes.content}>
        <section class={classes.group} aria-label="Actions">
          <For each={home.actions} key={(a) => a.id}>
            {(action, index) => <ActionBlock action={action} open={index.value === 0} />}
          </For>
        </section>
        {home.remoteStates ? <RemoteStates model={home.remoteStates} /> : null}
      </div>
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
  return (
    <div class={classes.action} data-open={open.asReadonly.pipe(map((o) => (o ? '' : undefined)))}>
      <button type="button" class={classes.actionHeader} aria-expanded={open} onClick={() => open.set(!open.value)}>
        <span class={classes.actionTitle}>{action.title}</span>
        <span class={classes.actionDescription}>{action.description}</span>
        <Icon name="chevron-down" size={14} class={classes.chevron} />
      </button>
      <Show when={open}>
        {() => (
          <div class={classes.actionBody}>
            {action.form ? <ParamsForm form={action.form} /> : null}
            <ActionFooter
              icon={action.icon}
              label={action.applyLabel ?? 'Apply'}
              busy={action.busy}
              onApply={() => action.apply()}
              onReset={action.form ? () => action.form?.reset() : undefined}
            />
          </div>
        )}
      </Show>
    </div>
  );
}

/** Sessions on a state server: the list (newest first), refresh, and the server behind a settings button. */
function RemoteStates(props: { model: RemoteStatesModel }) {
  const { model } = props;
  const settings = observableState(false);
  return (
    <section class={classes.group} aria-label="Remote states">
      <div class={classes.captionRow}>
        <h3 class={classes.caption}>Remote states</h3>
        <span class={classes.captionActions}>
          <Button
            size="xs"
            variant="ghost"
            icon={<Icon name="settings" />}
            title="Server"
            pressed={settings.asReadonly}
            onClick={() => settings.set(!settings.value)}
          />
          <Button
            size="xs"
            variant="ghost"
            icon={
              <Show when={model.loading} fallback={() => <Icon name="refresh-cw" />}>
                {() => <Icon name="loader-circle" class={classes.spin} />}
              </Show>
            }
            title="Refresh"
            onClick={() => model.refresh()}
          />
        </span>
      </div>
      {/* The server is a setting: behind the button, in a well under the caption (always shown in the Remote states card). */}
      <Show when={settings}>
        {() => (
          <label class={classes.server}>
            <span>Server</span>
            <TextInput
              size="sm"
              aria-label="State server URL"
              value={model.serverUrl}
              onSubmit={(v) => {
                model.setServerUrl(v);
                settings.set(false);
              }}
            />
          </label>
        )}
      </Show>
      {/* A long list scrolls on its own. */}
      <div class={classes.states}>
        <For
          each={model.states}
          key={(r) => r.id}
          fallback={() => <p class={classes.none}>No states on this server.</p>}
        >
          {(state) => (
            <button type="button" class={classes.state} onClick={() => model.open(state.id)}>
              <span class={classes.stateLabel} data-tip-overflow="">
                {state.label}
              </span>
              <span class={classes.stateDate}>{state.date}</span>
            </button>
          )}
        </For>
      </div>
    </section>
  );
}
