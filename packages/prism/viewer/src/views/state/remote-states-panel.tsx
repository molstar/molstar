/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Button } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ScrollArea } from '@molstar/prism-ui';
import { TextInput } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { HomeAction, RemoteState, RemoteStatesModel } from '../../model/home.js';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ActionBlock } from '../home/home-tab.js';
import { classes } from './_generated/remote-states-panel.module.css.js';

/**
 * Remote states (Mol*'s state server; present only when the spec turns it on): upload the current state under a name,
 * the states on the server (a click opens one; its ⋯ downloads or removes it), and the server's address.
 */
/**
 * Uploading the current state to the state server: an action like Home's (a name, then Upload across the footer). The
 * Export card's section (the Remote states card lists what is there); nothing when the server takes no uploads.
 */
export function RemoteStatesUpload(props: { model: RemoteStatesModel }) {
  const { model } = props;
  if (!model.upload) return null;
  const form = new ParamsFormModel({ name: PD.Text('', { label: 'Name', placeholder: 'Untitled' }) });
  const upload: HomeAction = {
    id: 'upload',
    title: 'Upload the current state',
    icon: 'cloud-upload',
    form,
    applyLabel: 'Upload',
    apply: () => {
      model.upload?.(form.values.value.name.trim() || 'Untitled');
      form.reset();
    },
  };
  return (
    <section class={classes.upload} aria-label="Upload">
      <ActionBlock action={upload} open />
    </section>
  );
}

/** The states on the state server (a click opens one; ⋯ downloads or removes) and the server's address. */
export function RemoteStatesPanel(props: { model: RemoteStatesModel }) {
  const { model } = props;
  return (
    <div class={classes.panel}>
      <section class={classes.states} aria-label="On the server">
        <div class={classes.captionRow}>
          <p class={classes.caption}>On the server</p>
          <Button
            size="xs"
            variant="ghost"
            title="Refresh"
            icon={
              <Show when={model.loading} fallback={() => <Icon name="refresh-cw" />}>
                {() => <Icon name="loader-circle" class={classes.spin} />}
              </Show>
            }
            onClick={() => model.refresh()}
          />
        </div>
        <ScrollArea class={classes.list} contentClass={classes.listContent}>
          <For each={model.states} by={(s) => s.id} fallback={() => <p class={classes.empty}>No states yet.</p>}>
            {(state) => <StateRow model={model} state={state} />}
          </For>
        </ScrollArea>
      </section>
      <label class={classes.server}>
        <span>Server</span>
        <TextInput
          size="sm"
          aria-label="State server URL"
          value={model.serverUrl}
          onSubmit={(v) => model.setServerUrl(v)}
        />
      </label>
    </div>
  );
}

function StateRow(props: { model: RemoteStatesModel; state: RemoteState }) {
  const { model, state } = props;
  const menu = observableState(false);
  const hasMenu = !!(model.download || model.remove);
  return (
    <div class={classes.item}>
      <div class={classes.row}>
        <button type="button" class={classes.open} onClick={() => model.open(state.id)}>
          <span class={classes.name}>{state.label}</span>
          <span class={classes.date}>{state.date}</span>
        </button>
        {hasMenu ? (
          <button
            type="button"
            class={classes.action}
            aria-label="More"
            aria-expanded={menu}
            onClick={() => menu.set(!menu.value)}
          >
            <Icon name="ellipsis" size={15} />
          </button>
        ) : null}
      </div>
      <Show when={menu}>
        {() => (
          <div class={classes.menu} role="menu">
            {model.download ? (
              <button
                type="button"
                role="menuitem"
                class={classes.menuItem}
                onClick={() => {
                  menu.set(false);
                  model.download?.(state.id);
                }}
              >
                <span>Download</span>
              </button>
            ) : null}
            {model.remove ? (
              <button
                type="button"
                role="menuitem"
                class={classes.menuItem}
                data-danger=""
                onClick={() => {
                  menu.set(false);
                  model.remove?.(state.id);
                }}
              >
                <span>Remove</span>
              </button>
            ) : null}
          </div>
        )}
      </Show>
    </div>
  );
}
