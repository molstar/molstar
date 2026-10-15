/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, type Observable, of } from 'rxjs';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { classes } from './_generated/scene-tree.module.css.js';

export interface SceneRowProps {
  readonly label: string;
  /** The short second line ("C2 axes", "30 Å thick"). */
  readonly description?: Observable<string | undefined>;
  /** Whether what the row stands for is shown (its eye); absent: the row has no eye. */
  readonly shown?: Observable<boolean>;
  /** The eye's click: show or hide it, through the API. */
  onToggle?(): void;
  /** The row's controls, under it while it is open (rendered on first open). */
  readonly body?: () => Child;
}

/**
 * An extension's row in the scene tree (viewer-spec.md §3, `rows`), drawn as the tree draws its own: a label and a
 * description, an eye, and the controls under it when opened. Its depth comes from the tree (`--msp-tree-depth` on the
 * row's host); it acts only through what it is given, so the extension's view calls its operations.
 */
export function SceneRow(props: SceneRowProps) {
  const open = observableState(false);
  const shown = props.shown ?? of(true);
  const hidden = shown.pipe(map((s) => !s));
  const labelId = uniqueId('msp-tree-label');
  const toggleOpen = () => props.body && open.set(!open.value);
  return (
    <div
      role="treeitem"
      aria-expanded={props.body ? open.asReadonly : undefined}
      aria-labelledby={`${labelId} ${labelId}-d`}
      tabIndex={-1}
      class={classes.item}
    >
      <div
        data-tip-row=""
        data-tree-row=""
        class={classes.row}
        classList={{ [classes.hidden]: hidden }}
        onClick={toggleOpen}
      >
        <button
          type="button"
          class={classes.expander}
          classList={{ [classes.leaf]: !props.body }}
          data-tree-toggle=""
          tabIndex={-1}
          aria-label="Expand"
          onClick={(e) => {
            e.stopPropagation();
            toggleOpen();
          }}
        >
          <Icon name="chevron-right" size={14} />
        </button>
        <span class={classes.label} id={labelId} data-tip-part="">
          {props.label}
        </span>
        <span class={classes.description} id={`${labelId}-d`} data-tip-part="">
          {props.description ?? ''}
        </span>
        <span class={classes.actions}>
          {props.onToggle ? (
            <button
              type="button"
              class={classes.action}
              classList={{ [classes.keep]: hidden }}
              title={hidden.pipe(map((h) => (h ? 'Show' : 'Hide')))}
              onClick={(e) => {
                e.stopPropagation();
                props.onToggle?.();
              }}
            >
              <Show when={hidden} fallback={() => <Icon name="eye" size={14} />}>
                {() => <Icon name="eye-off" size={14} />}
              </Show>
            </button>
          ) : (
            <span class={classes.actionSlot} aria-hidden="true" />
          )}
          {/* The menu's column, so the eye lines up with the tree's. */}
          <span class={classes.actionSlot} aria-hidden="true" />
        </span>
      </div>
      {props.body ? (
        <Show when={open.asReadonly} reveal>
          {() => (
            <div class={classes.rowBody} role="group">
              {props.body!()}
            </div>
          )}
        </Show>
      ) : null}
    </div>
  );
}
