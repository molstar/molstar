/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { Observable } from 'rxjs';
import { BuiltInAnnotations } from '../api/annotations.js';
import type { ViewerApi } from '../api/area.js';
import type { KeptRecord } from '../journal/record.js';

/**
 * A state action the Home tab offers: one of the actions applicable to the state tree's root (download structure,
 * download density, load trajectory, open session, and whatever extensions register), with its parameters as a form.
 * Built from `plugin.state.data.actions` for the root, so a new action shows up without UI code.
 */
export interface HomeAction {
  /** The action's id (the transformer id). */
  id: string;
  title: string;
  icon: IconName;
  /** A short line under the title ("PDB, AlphaFold DB, URL…"). */
  description?: string;
  form?: ParamsFormModel<any>;
  /** The apply button's label ("Download", "Load"); "Apply" by default. */
  applyLabel?: string;
  /** Running: the button shows progress and is disabled. */
  busy?: ReadonlyState<boolean>;
  apply(): void | Promise<void>;
}

/** A session saved on a state server, to open in one click. */
export interface RemoteState {
  id: string;
  label: string;
  /** When it was saved ("Sun, 11 Oct 2026 09:05"). */
  date?: string;
}

/**
 * Sessions on a state server (Mol*'s Remote States). Optional: present only when the spec enables it (a feature flag,
 * off by default), so apps that do not run a state server never show it.
 */
export interface RemoteStatesModel {
  serverUrl: ReadonlyState<string>;
  setServerUrl(url: string): void;
  states: ReadonlyState<readonly RemoteState[]>;
  loading: ReadonlyState<boolean>;
  refresh(): void;
  open(id: string): void;
  /** Uploads the current state under a name (the Remote states panel). */
  upload?(name: string): void;
  /** Downloads a state from the server as a file. */
  download?(id: string): void;
  remove?(id: string): void;
}

export interface HomeModel {
  /** The root's applicable actions, the most used first. */
  actions: ReadonlyState<readonly HomeAction[]>;
  /** File formats the open-files action reads, for the drop area's hint ("mmCIF, PDB, CCP4, MOLX, …"). */
  formats: readonly string[];
  /** Files picked or dropped: opened with the plugin's file-format detection. */
  openFiles(files: readonly File[]): void;
  /** Sessions on a state server; absent unless the spec turns remote states on. */
  remoteStates?: RemoteStatesModel;
  /** What is loaded now, and clearing it. */
  session?: SessionModel;
  /** After a load: what the new structure carries beyond the built-ins, offered in one click (parity.md §6a). */
  annotationOffer?: AnnotationOfferModel;
  /** The records kept in earlier sessions (journal.md §3): listed while shown, one opened read-only in History. */
  keptRecords?: KeptRecordsHomeModel;
}

export interface KeptRecordsHomeModel {
  /** Newest first, this session's own aside. */
  list: ReadonlyState<readonly KeptRecord[]>;
  /** Lists now and after each keep or forget, until the returned function is called. */
  watch(): () => void;
  /** Shows it in History, read-only, beside whatever is loaded. */
  open(id: string): void;
  forget(id: string): Promise<void>;
}

/** "AF-P04637 has pLDDT": an extension's annotation a structure that just loaded carries. */
export interface AnnotationOffer {
  /** The structure's state ref (what `annotations.apply` gets) and its label (what the line says). */
  structure: string;
  label: string;
  /** The annotation's id (`plddt`, `qmean`, `validation`) and title. */
  annotation: string;
  title: string;
  description?: string;
}

export interface AnnotationOfferModel {
  /** The offer, while there is one; it follows the structures only while something watches it. */
  offer: Observable<AnnotationOffer | undefined>;
  /** Applies it (`annotations.apply`, one step) and puts the offer away. */
  accept(): void;
  dismiss(): void;
}

const BuiltInIds = new Set(BuiltInAnnotations.map((a) => a.id));

/**
 * Home's result line after a load (parity.md §6a, README decision 27; v5's "Automatic (w/ Annotation)" preset as an
 * offer): when a structure appears that carries an extension's annotation (pLDDT, QMEAN, validation; never the
 * built-ins, B-factors and secondary structure, which every file has), offer it. Each structure is asked once, by
 * `annotations.list` as the person (a read: no step); a declined or applied offer does not come back for it.
 */
export function createAnnotationOffer(
  api: ViewerApi,
  structures: ReadonlyState<readonly { readonly ref: string; readonly label: string }[]>,
): AnnotationOfferModel {
  const current = observableState<AnnotationOffer | undefined>(undefined);
  const asked = new Set<string>();
  /** Asked, and not read yet: a check superseded by a newer one leaves its structures to it. */
  let waiting: { readonly ref: string; readonly label: string }[] = [];
  let token = 0;

  const check = async (list: readonly { readonly ref: string; readonly label: string }[]) => {
    // Nothing to offer in a viewer without the `annotations` feature (viewer-spec.md §3.1).
    if (!api.area('annotations')) return;
    const shown = current.value;
    if (shown && !list.some((s) => s.ref === shown.structure)) current.set(undefined);
    const fresh = list.filter((s) => !asked.has(s.ref));
    if (!fresh.length) return;
    for (const s of fresh) asked.add(s.ref);
    waiting = [...waiting.filter((w) => list.some((s) => s.ref === w.ref)), ...fresh];
    const mine = ++token;
    // The newest first: the offer is about what just loaded.
    for (const s of [...waiting].reverse()) {
      const r = await api.annotations.list({ structure: s.ref });
      if (mine !== token) return;
      waiting = waiting.filter((w) => w !== s);
      const found = r.ok ? r.value?.[0]?.annotations.find((a) => !BuiltInIds.has(a.id)) : undefined;
      if (!found) continue;
      const provider = api.annotations.providers.value.find((p) => p.id === found.id);
      waiting = [];
      current.set({
        structure: s.ref,
        label: s.label,
        annotation: found.id,
        title: found.title,
        description: provider?.description,
      });
      return;
    }
  };

  const offer = new Observable<AnnotationOffer | undefined>((subscriber) => {
    const watch = structures.subscribe((list) => void check(list));
    const out = current.asReadonly.subscribe(subscriber);
    return () => {
      watch.unsubscribe();
      out.unsubscribe();
    };
  });

  return {
    offer,
    accept() {
      const o = current.value;
      if (!o) return;
      current.set(undefined);
      void api.annotations.apply({ annotation: o.annotation, structure: o.structure });
    },
    dismiss: () => current.set(undefined),
  };
}

/** The session as Home shows it: what is loaded (a line), and clearing it (the plugin's clear, after a confirmation). */
export interface SessionModel {
  /** Nothing loaded. */
  isEmpty: ReadonlyState<boolean>;
  /** What is loaded, in a line ("1TQN · 3,624 atoms"). */
  summary: ReadonlyState<string>;
  /** Removes everything: data, representations, snapshots' current state (not saved snapshots or files). */
  clear(): void;
}
