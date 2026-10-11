/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
import type { ParamsFormModel } from '../forms/params-form.js';

/**
 * A group of settings by intent (Appearance, Camera, Interaction, Performance, Interface, Behaviors, Advanced): a form
 * over the plugin's params (canvas, canvas context, layout, behaviors) mapped into it, so no param is lost and new ones
 * land in Advanced. Changes apply as they are made; reset brings the group back to its defaults.
 */
export interface SettingsGroup {
  id: string;
  title: string;
  icon: IconName;
  form: ParamsFormModel<any>;
  /** Labels are long here (behaviors: "Camera Focus Loci on Canvas"): give them the row, the controls are small. */
  wideLabels?: boolean;
}

/** The Settings tab (implemented over Canvas3D, Canvas3DContext, the layout, and `plugin.state.behaviors`). */
export interface SettingsModel {
  groups: readonly SettingsGroup[];
  resetAll(): void;
}
