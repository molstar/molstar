/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
import type { ParamsFormModel } from '../forms/params-form.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

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

/** Something opened before, to open again in one click. */
export interface RecentEntry {
  id: string;
  label: string;
  /** Where it came from and when ("PDB · 2 min ago"). */
  description?: string;
  icon: IconName;
  open(): void;
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
}

export interface HomeModel {
  /** The root's applicable actions, the most used first. */
  actions: ReadonlyState<readonly HomeAction[]>;
  recent: ReadonlyState<readonly RecentEntry[]>;
  /** File formats the open-files action reads, for the drop area's hint ("mmCIF, PDB, CCP4, MOLX, …"). */
  formats: readonly string[];
  /** Files picked or dropped: opened with the plugin's file-format detection. */
  openFiles(files: readonly File[]): void;
  /** Sessions on a state server; absent unless the spec turns remote states on. */
  remoteStates?: RemoteStatesModel;
}
