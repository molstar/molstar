/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Choices } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { map } from 'rxjs';
import type { RepresentationQuickModel } from '../../model/representation-quick.js';
import { QuickColorRow, QuickRange } from '../common/quick-controls.js';
import { classes } from './_generated/representation-quick.module.css.js';

// A structure representation's quick settings (model/representation-quick.ts), above its whole form: its type and color
// theme as choices (a click applies one), one color from a row of swatches (the whole picker behind the chevron), the
// opacity and the size as sliders (a drag is one undo step).

export function RepresentationQuick(props: { model: RepresentationQuickModel }) {
  const { model } = props;
  return (
    <section class={classes.quick} aria-label="Quick settings">
      {/* The choices change with the structure (what applies) only: rebuilt when they do. */}
      <For each={model.types.pipe(map((t) => [t]))} key={(t) => t.map((c) => c.id).join()}>
        {(types) => (
          <Choices
            label="Type"
            choices={types}
            active={model.type}
            onPick={(id) => {
              model.change({ type: id });
              model.settle();
            }}
          />
        )}
      </For>
      <For each={model.themes.pipe(map((t) => [t]))} key={(t) => t.map((c) => c.id).join()}>
        {(themes) => (
          <Choices
            label="Color"
            choices={themes}
            active={model.theme}
            onPick={(id) => {
              model.change({ theme: id, color: undefined });
              model.settle();
            }}
          />
        )}
      </For>
      <Show when={model.theme.pipe(map((t) => t === 'uniform'))}>
        {() => (
          <QuickColorRow
            label="Color"
            value={model.color}
            onPick={(c) => {
              model.change({ theme: 'uniform', color: c });
              model.settle();
            }}
            onInput={(c) => model.change({ theme: 'uniform', color: c })}
          />
        )}
      </Show>
      <Show when={model.opacity.pipe(map((o) => o !== undefined))}>
        {() => (
          <QuickRange
            label="Opacity"
            min={0}
            max={1}
            step={0.01}
            value={model.opacity.pipe(map((o) => o ?? 1))}
            format={(v) => `${Math.round(v * 100)}%`}
            onInput={(v) => model.change({ opacity: v })}
            onCommit={() => model.settle()}
          />
        )}
      </Show>
      <Show when={model.size.pipe(map((s) => s !== undefined))}>
        {() => (
          <QuickRange
            label="Size"
            min={0.25}
            max={3}
            step={0.05}
            value={model.size.pipe(map((s) => s ?? 1))}
            format={(v) => `${v.toFixed(2)}×`}
            onInput={(v) => model.change({ size: v })}
            onCommit={() => model.settle()}
          />
        )}
      </Show>
    </section>
  );
}
