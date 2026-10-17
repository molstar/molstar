/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { ActionFooter } from '@molstar/prism-ui';
import { ParamsForm } from '@molstar/prism-ui';
import { For } from '@molstar/prism-ui/dom/control-flow';
import type { SuperpositionModel } from '../../model/superposition.js';
import { Superposed } from '../selection/superposed.js';
import { classes } from './_generated/superpose.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
/**
 * Superposing structures (the Selection tab's Superpose section, after Mol*'s): the method and its options, what it
 * would superpose (one selection per structure, a click focuses one; or what to select), Superpose, and the
 * superpositions made (each removed by itself, or all).
 */
export function Superpose(props: { model: SuperpositionModel }) {
  const { model } = props;
  return (
    <div class={classes.superpose}>
      <ParamsForm form={model.form} />
      <div class={classes.entries}>
        <For each={model.entries} key={(e) => e.id} fallback={() => <p class={classes.hint}>{model.hint}</p>}>
          {(entry, index) => (
            <button
              type="button"
              class={classes.entry}
              title={entry.problem ? `${entry.label}: ${entry.problem}` : 'Focus'}
              onClick={() => model.focus(entry.id)}
            >
              <span class={classes.index}>{derive(index, (i) => `${i + 1}`)}</span>
              <span class={classes.label}>{entry.label}</span>
              {entry.problem ? <span class={classes.problem}>{entry.problem}</span> : null}
              {!entry.problem && index.value === 0 ? <span class={classes.role}>reference</span> : null}
            </button>
          )}
        </For>
      </div>
      <ActionFooter
        icon="superpose"
        label="Superpose"
        disabled={derive(model.ready, (r) => !r)}
        info={model.blocked?.pipe(map((b) => b ?? ''))}
        onReset={() => model.form.reset()}
        onApply={() => model.superpose()}
      />
      {/* The superpositions made: removed here (they go back where they were loaded), as measurements are. */}
      {model.made && model.remove ? (
        <Superposed made={model.made} remove={(id) => model.remove?.(id)} clear={model.clear} />
      ) : null}
    </div>
  );
}
