/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A params form's own words (styles §8): strings, or functions of what they name. The params' labels and descriptions
// are the app's (PD); these are the form's around them.

export interface ParamsFormLabels {
  /** The fold of the params that are not essential. */
  readonly advanced: string;
  /** The presets' button and row. */
  readonly presets: string;
  /** The presets' accessible name, by what they are for. */
  readonly presetsOf: (label: string) => string;
  /** Values that match no preset (a group's, a color list's). */
  readonly custom: string;
  /** The button that opens a choice's or a switch's settings. */
  readonly settingsOf: (label: string) => string;
  /** A field's help toggle. */
  readonly about: (label: string) => string;
  /** A file field: empty (one file, several), its clear button, what it holds. */
  readonly selectFile: string;
  readonly selectFiles: string;
  readonly clearFile: string;
  readonly files: (names: readonly string[]) => string;
  /** An object list: its count, its add button, an item's edit and remove buttons. */
  readonly items: (count: number) => string;
  readonly addTo: (label: string) => string;
  readonly edit: string;
  readonly remove: string;
  /** A multi-select with nothing chosen (unless the param says its own). */
  readonly none: string;
  /** A 4 × 4 matrix: the identity, any other, a cell's accessible name (row and column from 1). */
  readonly identity: string;
  readonly matrix: string;
  readonly matrixCell: (label: string, row: number, column: number) => string;
  /** A vector's number by its axis ('x', 'y', 'z'); `label` undefined when the row has none. */
  readonly axis: (label: string | undefined, axis: string) => string;
  /** An interval's end. */
  readonly intervalEnd: (label: string | undefined, end: 'from' | 'to') => string;
  /** A line graph's points, and its editor's hint. */
  readonly points: (count: number) => string;
  readonly graphHint: string;
}

export const DefaultParamsFormLabels: ParamsFormLabels = {
  advanced: 'Advanced',
  presets: 'Presets',
  presetsOf: (label) => `${label} presets`,
  custom: 'Custom',
  settingsOf: (label) => `${label} settings`,
  about: (label) => `About ${label}`,
  selectFile: 'Select a file…',
  selectFiles: 'Select files…',
  clearFile: 'Clear',
  files: (names) => (names.length === 1 ? names[0] : `${names.length} files · ${names.join(', ')}`),
  items: (n) => (n === 1 ? '1 item' : `${n} items`),
  addTo: (label) => `Add to ${label}`,
  edit: 'Edit',
  remove: 'Remove',
  none: 'None',
  identity: 'Identity',
  matrix: '4 × 4 matrix',
  matrixCell: (label, row, column) => `${label}, row ${row}, column ${column}`,
  axis: (label, axis) => (label ? `${label} ${axis}` : axis),
  intervalEnd: (label, end) => `${label ? `${label} ` : ''}${end}`,
  points: (n) => `${n} points`,
  graphHint: 'Drag a point · click to add · double-click to remove',
};
