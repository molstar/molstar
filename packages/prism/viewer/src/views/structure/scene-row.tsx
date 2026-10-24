/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon, type IconName } from '@molstar/prism-ui/components/icon/icon';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import type { Listenable } from '@molstar/prism-ui/reactive/listen';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive, State } from '@molstar/prism-ui/reactive/state';
import { classes } from './_generated/scene-tree.module.css.js';

/** A row without an eye: shown, always. */
const Shown: ReadonlyState<boolean> = new State(true).asReadonly;

/** The event a row sends when clicked: the tree's host (the Structure tab) shows it in its inspector. */
export const SceneRowPickEvent = 'msp-scene-row-pick';

/** A picked row, as the inspector shows it: its label, line, eye and controls; `mark` shows it picked in the tree. */
export interface SceneRowPick {
  readonly label: string;
  readonly description?: ReadonlyState<string | undefined>;
  readonly shown?: ReadonlyState<boolean>;
  onToggle?(): void;
  readonly body?: () => Child;
  mark(on: boolean): void;
  /**
   * Tells once the row is gone (its structure removed, its `when` false, a search leaving it out; not when it only
   * scrolls out of the tree's window): what it controls may be gone too, so the inspector lets it go.
   */
  readonly gone: Listenable;
}

/**
 * The element an extension's row is mounted in, in the scene tree (`data-scene-row`): the tree tells when the row is gone
 * (its key left the rows), because the row's own scope also ends when it only scrolls out of the window.
 */
export interface SceneRowHost extends HTMLElement {
  sceneRowGone?(fn: () => void): () => void;
  /** Told by a row when its item is made (the host writes the tree's ARIA on it). */
  sceneRowSettle?(): void;
}

export interface SceneRowProps {
  readonly label: string;
  /** The short second line ("C2 axes", "30 Å thick"). */
  readonly description?: ReadonlyState<string | undefined>;
  /** Whether what the row stands for is shown (its eye); absent: the row has no eye. */
  readonly shown?: ReadonlyState<boolean>;
  /** The eye's click: show or hide it, through the API. */
  onToggle?(): void;
  /** The row's controls, under it while it is open (rendered on first open). */
  readonly body?: () => Child;
  /**
   * A few controls on the row itself, after its description (a trajectory's scrubber): they act without picking the
   * row, and stay out of the tree's keyboard (its arrows move between rows; the inspector has them all).
   */
  readonly inline?: () => Child;
  /**
   * Opens what it stands for in a panel of its own (a card with room for a plot) instead of the inspector: the row's
   * click, with a panel glyph at its end.
   */
  onOpen?(): void;
  /** Its kind's icon, before the label (a tree of its own: a kinemage, its masters). */
  readonly icon?: IconName;
  /**
   * Rows under it, a level deeper, while it is expanded (its chevron and a click toggle them): a tree of its own outside
   * the scene tree (the Kinemage card's groups and subgroups). Not in the scene tree, whose rows are flat.
   */
  readonly children?: () => Child;
  /** Expanded at first (`children`). */
  readonly expanded?: boolean;
  /** Its level outside the scene tree (the tree's host sets it there): its indent. */
  readonly depth?: number;
  /** What its click does in place of picking it (a kinemage's view: applied): an action, not an item to inspect. */
  onActivate?(): void;
}

/**
 * An extension's row in the scene tree (viewer-spec.md §3, `rows`), drawn as the tree draws its own: a label and a
 * description, an eye; picked, its controls in the inspector (`SceneRowPickEvent`; under it with no host to take them,
 * outside the tree, whose rows are of one height). Its depth comes from the tree (`--msp-tree-depth` on the
 * row's host); it acts only through what it is given, so the extension's view calls its operations.
 */
export function SceneRow(props: SceneRowProps) {
  const open = observableState(!!props.children && !!props.expanded);
  const hidden = derive(props.shown ?? Shown, (s) => !s);
  const labelId = uniqueId('msp-tree-label');
  let rowEl: HTMLElement | undefined;
  // Told once, when the row goes; a listener added after that is not kept.
  let goneListeners: Set<() => void> | undefined = new Set();
  const gone: Listenable = {
    listen(fn) {
      goneListeners?.add(fn);
      return () => void goneListeners?.delete(fn);
    },
  };
  onCleanup(() => {
    const told = goneListeners;
    goneListeners = undefined;
    if (told) for (const fn of told) fn();
  });
  // In the scene tree (its rows are of one height: virtualized), nothing opens under it.
  const inTree = () => !!rowEl?.parentElement?.closest('[data-scene-row]');
  // A click picks it, as the tree's other rows: its controls go to the inspector (the tree's host takes the event);
  // with no host to take it (a row on its own, outside the tree), they open under it (a click on the row toggles them;
  // the eye only opens them: `reveal`).
  const pick = (reveal = false) => {
    if (props.onActivate) {
      props.onActivate();
      return;
    }
    // A tree of its own: a click opens and folds its rows.
    if (props.children) {
      open.set(!open.value);
      return;
    }
    if (props.onOpen) {
      props.onOpen();
      return;
    }
    const host = rowEl?.closest<SceneRowHost>('[data-scene-row]');
    const detail: SceneRowPick = {
      label: props.label,
      ...(props.description ? { description: props.description } : {}),
      ...(props.shown ? { shown: props.shown } : {}),
      ...(props.onToggle ? { onToggle: () => props.onToggle?.() } : {}),
      ...(props.body ? { body: props.body } : {}),
      mark: (on) => rowEl?.classList.toggle(classes.picked, on),
      gone: host?.sceneRowGone ? { listen: (fn) => host.sceneRowGone!(fn) } : gone,
    };
    const event = new CustomEvent(SceneRowPickEvent, { bubbles: true, composed: true, cancelable: true, detail });
    if (rowEl?.dispatchEvent(event) !== false && !inTree()) {
      if (reveal) open.set(true);
      else toggleOpen();
    }
  };
  const toggleOpen = () => (props.children || (props.body && !inTree())) && open.set(!open.value);
  return (
    <div
      role="treeitem"
      aria-expanded={props.body || props.children ? open.asReadonly : undefined}
      aria-labelledby={`${labelId} ${labelId}-d`}
      tabIndex={-1}
      class={classes.item}
      style={props.depth !== undefined ? `--msp-tree-depth: ${props.depth}` : undefined}
      ref={(el: HTMLElement) =>
        // In its host by then (a view may render a row after the host's mount): the host writes its tree ARIA on it.
        queueMicrotask(() => el.closest<SceneRowHost>('[data-scene-row]')?.sceneRowSettle?.())
      }
    >
      <div
        data-tip-row=""
        data-tree-row=""
        class={classes.row}
        classList={{ [classes.hidden]: hidden }}
        ref={(el) => {
          rowEl = el;
        }}
        onClick={() => pick()}
      >
        <button
          type="button"
          class={classes.expander}
          // Its controls are in the inspector: nothing opens under it; rows of its own do.
          classList={{ [classes.leaf]: !props.children }}
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
        {props.icon ? <Icon name={props.icon} size={15} class={classes.kind} /> : null}
        <span class={classes.label} id={labelId} data-tip-part="">
          {props.label}
        </span>
        <span
          class={classes.description}
          classList={{ [classes.short]: !!props.inline }}
          id={`${labelId}-d`}
          data-tip-part=""
        >
          {props.description ?? ''}
        </span>
        {props.inline ? (
          // Its clicks act, they do not pick the row.
          <span class={classes.inline} onClick={(e) => e.stopPropagation()}>
            {props.inline()}
          </span>
        ) : null}
        <span class={classes.actions}>
          {props.onToggle ? (
            <button
              type="button"
              class={classes.action}
              classList={{ [classes.keep]: hidden }}
              data-tip={derive(hidden, (h) => (h ? 'Show' : 'Hide'))}
              aria-label={derive(hidden, (h) => (h ? 'Show' : 'Hide'))}
              onClick={(e) => {
                e.stopPropagation();
                props.onToggle?.();
                // Turned on, a row with controls is picked: they are in the inspector at once.
                if (hidden.value && props.body) pick(true);
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
            <span class={classes.actionSlot} data-tip="Opens in its own panel">
              <Icon name="panel-right" size={14} />
            </span>
          ) : (
            <span class={classes.actionSlot} aria-hidden="true" />
          )}
        </span>
      </div>
      {props.children ? (
        <Show when={open.asReadonly} reveal>
          {() => <div role="group">{props.children!()}</div>}
        </Show>
      ) : null}
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
