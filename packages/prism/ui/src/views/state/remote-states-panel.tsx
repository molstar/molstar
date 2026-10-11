/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { HomeAction, RemoteState, RemoteStatesModel } from '../../model/home.js';
import { observableState } from '../../reactive/observable-state.js';
import { ActionBlock } from '../home/home-tab.js';
import { classes } from './_generated/remote-states-panel.module.css.js';

/**
 * Remote states (Mol*'s state server; present only when the spec turns it on): upload the current state under a name,
 * the states on the server (a click opens one; its ⋯ downloads or removes it), and the server's address.
 */
export function RemoteStatesPanel(props: { model: RemoteStatesModel }) {
  const { model } = props;
  // Uploading is an action like Home's: a name, set defaults, and Upload across the footer.
  const form = new ParamsFormModel({ name: PD.Text('', { label: 'Name', placeholder: 'Untitled' }) });
  const upload: HomeAction | undefined = model.upload
    ? {
        id: 'upload',
        title: 'Upload the current state',
        icon: 'cloud-upload',
        form,
        applyLabel: 'Upload',
        apply: () => {
          model.upload?.(form.values.value.name.trim() || 'Untitled');
          form.reset();
        },
      }
    : undefined;
  return (
    <div class={classes.panel}>
      {upload ? (
        <section class={classes.upload} aria-label="Upload">
          <ActionBlock action={upload} open />
        </section>
      ) : null}
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
          <For each={model.states} key={(s) => s.id} fallback={() => <p class={classes.empty}>No states yet.</p>}>
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
