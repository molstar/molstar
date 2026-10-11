/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
import type { ParamsFormModel } from '../forms/params-form.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import type { SceneNode } from './scene.js';

/** A point to measure from: one of the latest selections (Mol*'s selection additions history). */
export interface MeasurementPoint {
  id: string;
  /** "CYS 442 SG", "HEM 508 FE". */
  label: string;
}

/** What can be measured, and how many points it takes (0: the current selection as a whole). */
export interface MeasurementKind {
  id: 'distance' | 'angle' | 'dihedral' | 'label' | 'orientation' | 'plane';
  title: string;
  icon: IconName;
  points: number;
}

/**
 * Adding measurements (implemented over the measurement manager and the selection history): the latest selections, in
 * order, are the points (the first four are used); a measurement takes as many as it needs.
 */
export interface MeasurementsModel {
  /** Newest first. */
  points: ReadonlyState<readonly MeasurementPoint[]>;
  movePoint(id: string, by: -1 | 1): void;
  /** Moves it to a position (dragged by its handle). */
  movePointTo(id: string, index: number): void;
  removePoint(id: string): void;
  focusPoint(id: string): void;
  highlightPoint?(id: string | undefined): void;
  kinds: readonly MeasurementKind[];
  add(kind: MeasurementKind['id']): void;
  /** Distance unit, label text size and color (the measurement manager's options). */
  options: ParamsFormModel<any>;
  /**
   * The measurements made, to adjust after adding: the same rows as under Measurements in the Structure tab (show or
   * hide, remove), each with its style (the same form its row's inspector shows).
   */
  added: ReadonlyState<readonly SceneNode[]>;
  styleOf(node: SceneNode): ParamsFormModel<any> | undefined;
}
