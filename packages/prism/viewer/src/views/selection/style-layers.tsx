/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamsForm, ScrollArea } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { StyleLayerItem, StyleLayersModel } from '../../model/style-layers.js';
import { classes } from './_generated/selection-lists.module.css.js';

/**
 * The style layers applied (under the Selection tab's style form; Mol*'s overpaint, transparency, emissive, material and
 * clipping layers, which v5 only lets go of in the state tree): a row each (the swatch or the kind, the atoms, the
 * value), its value edited inline, removed; Clear all. The list keeps a height and scrolls by itself.
 */
export function StyleLayers(props: { model: StyleLayersModel }) {
  const { model } = props;
  const count = derive(model.items, (l) => l.length);
  return (
    <Show when={derive(count, (n) => n > 0)}>
      {() => (
        <div class={classes.block}>
          <div class={classes.head}>
            <span class={classes.caption}>
              Applied <span class={classes.count}>{count}</span>
            </span>
            <button type="button" class={classes.clear} onClick={() => model.clear()}>
              Clear all
            </button>
          </div>
          <ScrollArea class={classes.list} contentClass={classes.listContent}>
            <For each={model.items} by={(i) => i.key}>
              {(item) => <LayerRow item={item} />}
            </For>
          </ScrollArea>
        </div>
      )}
    </Show>
  );
}

function LayerRow(props: { item: StyleLayerItem }) {
  const { item } = props;
  const open = observableState(false);
  const view = item.view;
  const editable = item.editable;
  return (
    <div class={classes.item}>
      <div class={classes.row} title={derive(view, (v) => v.detail)}>
        <span class={classes.lead}>
          <Show when={derive(view, (v) => v.swatch)} fallback={() => <Icon name={view.value.icon} size={13} />}>
            {(swatch) => <span class={classes.swatch} style={`background: ${swatch}`} />}
          </Show>
        </span>
        <span class={classes.label}>{derive(view, (v) => v.where)}</span>
        <span class={classes.value}>{derive(view, (v) => v.value)}</span>
        <span class={classes.actions}>
          {editable ? (
            <button
              type="button"
              class={classes.action}
              title="Edit"
              aria-label="Edit"
              aria-expanded={open}
              onClick={() => open.set(!open.value)}
            >
              <Icon name="pencil" size={13} />
            </button>
          ) : null}
          <button type="button" class={classes.action} title="Remove" aria-label="Remove" onClick={() => item.remove()}>
            <Icon name="x" size={13} />
          </button>
        </span>
      </div>
      {editable ? (
        <Show when={open} reveal>
          {() => {
            const form = item.editor();
            return form ? (
              <div class={classes.editor}>
                <ParamsForm form={form} />
              </div>
            ) : null;
          }}
        </Show>
      ) : null}
    </div>
  );
}
