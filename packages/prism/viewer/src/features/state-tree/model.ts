/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ObservableState, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { HomeAction } from '../../model/home.js';
import type { StateTreeInspection, StateTreeModel, StateTreeNode } from '../../model/state-tree.js';
import { Debouncer } from '../../ui/debounce.js';
import type { StateApi } from './operations.js';
import type { CellDetails } from './source.js';

// The State tree tab's model (`model/state-tree.ts`) over the feature's source and its `state` area: the cells as the
// source reads them, edits and actions as the person through `state.update` / `state.action`, never the plugin.

interface Entry {
  /** What the forms were made for: a different shape (the transformer's params changed) makes new ones. */
  readonly shape: string;
  readonly updates: readonly { readonly ref: string; readonly form: ParamsFormModel<any> }[];
  readonly actions: readonly HomeAction[];
}

const shapeOf = (d: CellDetails) =>
  `${d.updates.map((u) => `${u.ref}:${Object.keys(u.params).join(',')}`).join(';')}|${d.actions.map((a) => a.name).join(',')}`;

export class LiveStateTreeModel implements StateTreeModel {
  readonly roots: ReadonlyState<readonly StateTreeNode[]>;
  readonly #selected: ObservableState<string | undefined> = observableState<string | undefined>(undefined);
  readonly selected: ReadonlyState<string | undefined> = this.#selected.asReadonly;
  readonly #cache = new Map<string, Entry>();
  /** Edits wait for the next one before they are applied (one form edit is one step). */
  readonly #edits = new Debouncer<string>((ref) => void this.#apply(ref));
  readonly #pending = new Map<string, Record<string, unknown>>();

  constructor(private readonly state: StateApi) {
    this.roots = state.tree.roots;
  }

  select(ref: string | undefined): void {
    this.#selected.set(ref);
  }

  inspect(ref: string): StateTreeInspection | undefined {
    const details = this.state.tree.details(ref);
    if (!details) return undefined;
    let entry = this.#cache.get(ref);
    if (!entry || entry.shape !== shapeOf(details)) {
      entry = this.#make(details);
      this.#cache.set(ref, entry);
    }
    return { node: details.node, updates: this.#updates(details, entry), actions: entry.actions };
  }

  toggleVisible(ref: string): void {
    const node = this.state.tree.node(ref);
    void this.state.visibility({ node: ref, visible: node?.visible === false });
  }

  remove(ref: string): void {
    void this.state.remove({ node: ref });
  }

  /** The forms take the cells' current values (the cells changed elsewhere: Undo, a command, a session). */
  refresh(refs?: ReadonlySet<string>): void {
    for (const [ref, entry] of this.#cache) {
      if (refs && !refs.has(ref) && !entry.updates.some((u) => refs.has(u.ref))) continue;
      const details = this.state.tree.details(ref);
      if (!details || entry.shape !== shapeOf(details)) {
        this.#cache.delete(ref);
        continue;
      }
      details.updates.forEach((u, i) => {
        // Not under the person's hands: an edit waiting to be applied wins.
        if (!this.#edits.has(u.ref)) entry.updates[i].form.load(u.values as never);
      });
    }
  }

  /** The selection goes when its cell does. */
  prune(): void {
    const ref = this.#selected.value;
    if (ref !== undefined && !this.state.tree.details(ref)) this.#selected.set(undefined);
  }

  /** Drops edits not applied yet (the view went). */
  cancel(): void {
    this.#edits.dispose();
    this.#pending.clear();
  }

  #updates(details: CellDetails, entry: Entry): HomeAction[] {
    return details.updates.map((u, i) => ({
      id: i === 0 ? 'update' : `update:${u.name}`,
      title: `Update ${u.name}`,
      icon: 'pencil',
      form: entry.updates[i].form,
      applyLabel: 'Update',
      apply: () => this.#apply(u.ref),
    }));
  }

  #make(details: CellDetails): Entry {
    const updates = details.updates.map((u) => ({
      ref: u.ref,
      form: new ParamsFormModel(u.params, u.values as never, {
        onChange: (values) => {
          this.#pending.set(u.ref, values);
          this.#edits.schedule(u.ref);
        },
      }),
    }));
    const target = details.node.ref;
    const actions = details.actions.map((a): HomeAction => {
      const form = Object.keys(a.params).length ? new ParamsFormModel(a.params) : undefined;
      const busy = observableState(false);
      return {
        id: a.name,
        title: a.title,
        icon: 'plus',
        ...(a.description ? { description: a.description } : {}),
        ...(form ? { form } : {}),
        applyLabel: a.trigger ?? 'Apply',
        busy: busy.asReadonly,
        apply: async () => {
          busy.set(true);
          try {
            await this.state.action({
              node: target,
              action: a.name,
              params: (form?.values.value ?? {}) as Record<string, unknown>,
            });
          } finally {
            busy.set(false);
          }
        },
      };
    });
    return { shape: shapeOf(details), updates, actions };
  }

  async #apply(ref: string): Promise<void> {
    this.#edits.cancel(ref);
    const values = this.#pending.get(ref);
    this.#pending.delete(ref);
    if (!values) return;
    const result = await this.state.update({ node: ref, params: values });
    // Not applied (invalid values, a failed update): the forms show what the cell has.
    if (!result.ok) this.refresh(new Set([ref]));
  }
}
