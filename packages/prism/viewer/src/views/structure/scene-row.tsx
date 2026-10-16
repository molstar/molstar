/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, type Observable, of, Subject } from 'rxjs';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { classes } from './_generated/scene-tree.module.css.js';

/** The event a row sends when clicked: the tree's host (the Structure tab) shows it in its inspector. */
export const SceneRowPickEvent = 'msp-scene-row-pick';

/** A picked row, as the inspector shows it: its label, line, eye and controls; `mark` shows it picked in the tree. */
export interface SceneRowPick {
  readonly label: string;
  readonly description?: Observable<string | undefined>;
  readonly shown?: Observable<boolean>;
  onToggle?(): void;
  readonly body?: () => Child;
  mark(on: boolean): void;
  /**
   * Emits once the row is gone (its structure removed, its `when` false, a search leaving it out): what it controls may
   * be gone too, so the inspector lets it go.
   */
  readonly gone: Observable<void>;
}

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
  /**
   * Opens what it stands for in a panel of its own (a card with room for a plot) instead of the inspector: the row's
   * click, with a panel glyph at its end.
   */
  onOpen?(): void;
}

/**
 * An extension's row in the scene tree (viewer-spec.md §3, `rows`), drawn as the tree draws its own: a label and a
 * description, an eye; picked, its controls in the inspector (`SceneRowPickEvent`; under it with no host to take them). Its depth comes from the tree (`--msp-tree-depth` on the
 * row's host); it acts only through what it is given, so the extension's view calls its operations.
 */
export function SceneRow(props: SceneRowProps) {
  const open = observableState(false);
  const shown = props.shown ?? of(true);
  const hidden = shown.pipe(map((s) => !s));
  let isHidden = false;
  const sub = hidden.subscribe((h) => {
    isHidden = h;
  });
  onCleanup(() => sub.unsubscribe());
  const labelId = uniqueId('msp-tree-label');
  let rowEl: HTMLElement | undefined;
  const gone = new Subject<void>();
  onCleanup(() => {
    gone.next();
    gone.complete();
  });
  // A click picks it, as the tree's other rows: its controls go to the inspector (the tree's host takes the event);
  // with no host to take it (a row on its own), they open under it.
  const pick = () => {
    if (props.onOpen) {
      props.onOpen();
      return;
    }
    const detail: SceneRowPick = {
      label: props.label,
      ...(props.description ? { description: props.description } : {}),
      ...(props.shown ? { shown: props.shown } : {}),
      ...(props.onToggle ? { onToggle: () => props.onToggle?.() } : {}),
      ...(props.body ? { body: props.body } : {}),
      mark: (on) => rowEl?.classList.toggle(classes.picked, on),
      gone,
    };
    const event = new CustomEvent(SceneRowPickEvent, { bubbles: true, composed: true, cancelable: true, detail });
    if (rowEl?.dispatchEvent(event) !== false) toggleOpen();
  };
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
        ref={(el) => {
          rowEl = el;
        }}
        onClick={pick}
      >
        <button
          type="button"
          class={classes.expander}
          // Its controls are in the inspector: nothing opens under it.
          classList={{ [classes.leaf]: true }}
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
                // Turned on, a row with controls is picked: they are in the inspector at once.
                if (isHidden && props.body) pick();
              }}
            >
              <Show when={hidden} fallback={() => <Icon name="eye" size={14} />}>
                {() => <Icon name="eye-off" size={14} />}
              </Show>
            </button>
          ) : (
            <span class={classes.actionSlot} aria-hidden="true" />
          )}
          {/* The menu's column, so the eye lines up with the tree's; a row that opens a panel shows it there. */}
          {props.onOpen ? (
            <span class={classes.actionSlot} title="Opens in its own panel">
              <Icon name="panel-right" size={14} />
            </span>
          ) : (
            <span class={classes.actionSlot} aria-hidden="true" />
          )}
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
