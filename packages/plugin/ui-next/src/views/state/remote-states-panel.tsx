/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import type { RemoteState, RemoteStatesModel } from '../../model/home.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/remote-states-panel.module.css.js';

/**
 * Remote states (Mol*'s state server; present only when the spec turns it on): upload the current state under a name,
 * the states on the server (a click opens one; its ⋯ downloads or removes it), and the server's address.
 */
export function RemoteStatesPanel(props: { model: RemoteStatesModel }) {
  const { model } = props;
  const name = observableState('');
  const upload = () => {
    model.upload?.(name.value.trim() || 'Untitled');
    name.set('');
  };
  return (
    <div class={classes.panel}>
      {model.upload ? (
        <section class={classes.upload} aria-label="Upload">
          <p class={classes.caption}>Upload the current state</p>
          <TextInput
            size="sm"
            aria-label="Name"
            placeholder="Name"
            value={name.asReadonly}
            onInput={(v) => name.set(v)}
            onSubmit={upload}
          />
          <div class={classes.footer}>
            <Button size="sm" variant="primary" icon={<Icon name="cloud-upload" />} onClick={upload}>
              Upload
            </Button>
          </div>
        </section>
      ) : null}
      <section class={classes.states} aria-label="On the server">
        <div class={classes.captionRow}>
          <p class={classes.caption}>On the server</p>
          <Button
            size="sm"
            variant="ghost"
            title="Refresh"
            icon={
              <Show when={model.loading} fallback={() => <Icon name="rotate-ccw" />}>
                {() => <Icon name="loader-circle" class={classes.spin} />}
              </Show>
            }
            onClick={() => model.refresh()}
          />
        </div>
        <div class={classes.list}>
          <For each={model.states} key={(s) => s.id} fallback={() => <p class={classes.empty}>No states yet.</p>}>
            {(state) => <StateRow model={model} state={state} />}
          </For>
        </div>
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
          <Icon name="cloud" size={15} />
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
                <Icon name="download" size={14} />
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
                <Icon name="trash-2" size={14} />
                <span>Remove</span>
              </button>
            ) : null}
          </div>
        )}
      </Show>
    </div>
  );
}
