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
import { map } from 'rxjs';
import type { SelectionSummaryModel } from '../../model/selection.js';
import { classes } from './_generated/selection-lists.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
/** At most this many structures are named; the rest are counted. */
const MaxParts = 3;

/**
 * What the Selection tab's sections act on (top of the tab): the selection in Mol*'s words, as granular as it is ("8
 * residues", "REA 200 | A"), then where it is, a line per structure (its chains and atoms), with Focus and Clear. With
 * nothing selected: how to select.
 */
export function SelectionSummary(props: { model: SelectionSummaryModel }) {
  const { model } = props;
  const shown = derive(model.parts, (p) => (p.length > MaxParts ? p.slice(0, MaxParts) : p));
  const more = derive(model.parts, (p) => Math.max(0, p.length - MaxParts));
  return (
    <div class={classes.summary} role="status" aria-label="Selection">
      <Show
        when={model.summary}
        fallback={() => (
          <>
            <div class={classes.summaryHead}>
              <Icon name="mouse-pointer-click" size={15} class={classes.summaryIcon} />
              <span class={classes.summaryTitle}>Nothing selected</span>
            </div>
            {/* Guided: out of Select mode, the way in; in it, how to pick and at what level. */}
            <Show
              when={model.selectMode}
              fallback={() => (
                <div class={classes.summaryOffer}>
                  <p class={classes.summaryHint}>
                    Clicks on the canvas focus now. Switch to Select to pick atoms and residues, or pick residues in the
                    sequence.
                  </p>
                  <button type="button" class={classes.summaryStart} onClick={() => model.startSelecting()}>
                    <Icon name="mouse-pointer-click" size={14} />
                    Select on the canvas
                  </button>
                </div>
              )}
            >
              {() => (
                <p class={classes.summaryHint}>
                  Click to pick a{' '}
                  <span class={classes.summaryLevel}>{derive(model.pickLevel, (l) => l.toLowerCase())}</span>; Shift
                  extends along the chain. Or pick residues in the sequence.
                </p>
              )}
            </Show>
          </>
        )}
      >
        {(summary) => (
          <>
            <div class={classes.summaryHead}>
              <Icon name="mouse-pointer-click" size={15} class={classes.summaryIcon} />
              <span class={classes.summaryTitle}>{summary}</span>
              <button
                type="button"
                class={classes.summaryAction}
                aria-label="Focus the selection"
                data-tip="Focus"
                onClick={() => model.focus()}
              >
                <Icon name="focus" size={15} />
              </button>
              <button
                type="button"
                class={classes.summaryAction}
                aria-label="Clear the selection"
                data-tip="Clear"
                onClick={() => model.clear()}
              >
                <Icon name="x" size={15} />
              </button>
            </div>
            <For each={shown} key={(p) => p.id}>
              {(part) => (
                <div class={classes.part}>
                  <span class={classes.partName}>{part.label}</span>
                  <span class={classes.partLine}>{part.detail}</span>
                </div>
              )}
            </For>
            <Show when={more.pipe(map((n) => n > 0))}>
              {() => (
                <div class={classes.part}>
                  <span class={classes.partLine}>{more.pipe(map((n) => `+${n} more structures`))}</span>
                </div>
              )}
            </Show>
          </>
        )}
      </Show>
    </div>
  );
}
