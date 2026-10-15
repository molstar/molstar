/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, type Observable } from 'rxjs';
import { formatChord } from '@molstar/prism-ui/commands/keys';
import { ActionFooter, type ModeChoice, ModeChooser } from '@molstar/prism-ui';
import { Button } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ParamsForm } from '@molstar/prism-ui';
import { ScrollArea } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import type { HomeAction, HomeModel, RemoteStatesModel, SessionModel } from '../../model/home.js';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { classes } from './_generated/home-tab.module.css.js';

export interface HomeTabProps {
  home: HomeModel;
  /** More sections after the actions (the `home.import` slot). */
  more?: () => Child;
  /** Opens the command menu (the hint at the bottom). */
  onOpenCommandMenu?: () => void;
  /** The app's modes, when it has more than one: offered while nothing is loaded (design/modes.md). */
  modes?: { list: readonly ModeChoice[]; active: Observable<string>; choose(id: string): void };
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
      <ScrollArea class={classes.content} contentClass={classes.contentInner}>
        {home.session ? <Session session={home.session} /> : null}
        {props.modes && props.modes.list.length > 1 ? (
          <Show when={home.session?.isEmpty ?? true}>
            {() => (
              <section class={classes.group} aria-label="Mode">
                <h3 class={classes.caption}>Mode</h3>
                <ModeChooser
                  variant="rows"
                  modes={props.modes!.list}
                  active={props.modes!.active}
                  onChoose={(id) => props.modes!.choose(id)}
                />
              </section>
            )}
          </Show>
        ) : null}
        <section class={classes.group} aria-label="Actions">
          <OpenFiles home={home} />
          <For each={home.actions} key={(a) => a.id}>
            {(action, index) => <ActionBlock action={action} open={index.value === 0} />}
          </For>
        </section>
        {props.more?.()}
        {home.remoteStates ? <RemoteStates model={home.remoteStates} /> : null}
      </ScrollArea>
      <button type="button" class={classes.hint} onClick={() => props.onOpenCommandMenu?.()}>
        <kbd>{formatChord({ key: 'k', mod: true })}</kbd>
        <span>Search every command</span>
      </button>
    </div>
  );
}

/** Open files: the file picker (parity.md §1); the files open as dropped ones do. */
function OpenFiles(props: { home: HomeModel }) {
  const input = (
    <input
      type="file"
      multiple
      hidden
      onChange={(e) => {
        const el = e.currentTarget as HTMLInputElement;
        const files = [...(el.files ?? [])];
        el.value = '';
        if (files.length) props.home.openFiles(files);
      }}
    />
  ) as HTMLInputElement;
  return (
    <div>
      <Button size="sm" variant="ghost" icon={<Icon name="folder-open" />} onClick={() => input.click()}>
        Open files…
      </Button>
      {input}
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
      <Show when={open} reveal>
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

/**
 * What is loaded, while something is: a line saying what, and Clear, which asks once (it removes everything loaded;
 * saved snapshots and files stay).
 */
function Session(props: { session: SessionModel }) {
  const { session } = props;
  const confirming = observableState(false);
  return (
    <Show when={session.isEmpty.pipe(map((e) => !e))}>
      {() => (
        <section class={classes.group} aria-label="Current session">
          <div class={classes.captionRow}>
            <h3 class={classes.caption}>Current session</h3>
          </div>
          <div class={classes.session}>
            <span class={classes.sessionSummary} data-tip-overflow="">
              {session.summary}
            </span>
            <Show
              when={confirming}
              fallback={() => (
                <Button size="sm" variant="ghost" icon={<Icon name="trash-2" />} onClick={() => confirming.set(true)}>
                  Clear
                </Button>
              )}
            >
              {() => (
                <span class={classes.confirm} role="group" aria-label="Clear everything loaded?">
                  <span class={classes.confirmText}>Clear everything?</span>
                  <Button size="sm" variant="ghost" onClick={() => confirming.set(false)}>
                    Cancel
                  </Button>
                  <Button
                    size="sm"
                    variant="danger"
                    onClick={() => {
                      confirming.set(false);
                      session.clear();
                    }}
                  >
                    Clear
                  </Button>
                </span>
              )}
            </Show>
          </div>
        </section>
      )}
    </Show>
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
      <ScrollArea class={classes.states} contentClass={classes.statesContent}>
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
      </ScrollArea>
    </section>
  );
}
