/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { combineLatest, map, type Observable } from 'rxjs';
import type { QuickChoice, RepresentationQuickModel } from '../../model/representation-quick.js';
import { QuickColorRow, QuickRange } from '../common/quick-controls.js';
import { classes } from './_generated/representation-quick.module.css.js';

// A structure representation's quick settings (model/representation-quick.ts), above its whole form. Each setting is a
// row: its name and the current value, then the choices drawn as what they are: the types as their icons in one segmented
// strip, the color themes as tiles painted with their colors on this structure (one color: its own swatch). One color
// opens the row of swatches (the whole picker behind its chevron); opacity and size are sliders. A drag is one undo step.

export function RepresentationQuick(props: { model: RepresentationQuickModel }) {
  const { model } = props;
  const labelOf = (list: Observable<readonly QuickChoice[]>, id: Observable<string>) =>
    combineLatest([list, id]).pipe(map(([l, i]) => l.find((c) => c.id === i)?.label ?? i));
  return (
    <section class={classes.quick} aria-label="Quick settings">
      <div class={classes.setting}>
        <div class={classes.head}>
          <span class={classes.name}>Type</span>
          <span class={classes.value}>{labelOf(model.types, model.type)}</span>
        </div>
        <div class={classes.segments} role="radiogroup" aria-label="Type">
          <For each={model.types} key={(c) => c.id}>
            {(c) => (
              <button
                type="button"
                role="radio"
                class={classes.segment}
                aria-checked={model.type.pipe(map((t) => t === c.id))}
                aria-label={c.label}
                data-tip={c.label}
                onClick={() => {
                  model.change({ type: c.id });
                  model.settle();
                }}
              >
                <Icon name={c.icon ?? 'shapes'} size={18} />
              </button>
            )}
          </For>
        </div>
      </div>
      <div class={classes.setting}>
        <div class={classes.head}>
          <span class={classes.name}>Color</span>
          <span class={classes.value}>{labelOf(model.themes, model.theme)}</span>
        </div>
        <div class={classes.tiles} role="radiogroup" aria-label="Color">
          <For each={model.themes} key={(c) => c.id}>
            {(c) => (
              <button
                type="button"
                role="radio"
                class={classes.tile}
                classList={{ [classes.unpainted]: !c.swatch && c.id !== 'uniform' }}
                style={{
                  background:
                    c.id === 'uniform'
                      ? model.color.pipe(map((v) => (v === undefined ? '' : `#${v.toString(16).padStart(6, '0')}`)))
                      : (c.swatch ?? ''),
                }}
                aria-checked={model.theme.pipe(map((t) => t === c.id))}
                aria-label={c.label}
                data-tip={c.label}
                onClick={() => {
                  model.change({ theme: c.id, color: undefined });
                  model.settle();
                }}
              >
                {c.id === 'uniform' ? <span class={classes.one} aria-hidden="true" /> : null}
              </button>
            )}
          </For>
        </div>
        <Show when={model.theme.pipe(map((t) => t === 'uniform'))} reveal>
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
      </div>
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
