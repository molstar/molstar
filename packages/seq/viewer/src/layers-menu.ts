/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Layers menu's two closing entries (_spec/README.md "Legends", A layer's eye; "Behaviors", Headers): a menu of
// several layers, the track header's and the legend's alike, ends with a separator, "Show all" and "Hide all". DOM-free,
// so both menus read one rule: each is one command over the layers it lists (one patch, one undo step), disabled with
// its reason when it would change nothing.

/** What a menu of layers lists of each. */
export interface LayerMenuLayer {
  readonly address: string;
  readonly hidden: boolean;
}

/** One of the two closing entries. */
export interface LayerMenuAll {
  readonly label: 'Show all' | 'Hide all';
  readonly command: 'layer.show' | 'layer.hide';
  /** The layers it changes: those in the other state, never one it would leave as it is. */
  readonly args: { readonly addresses: readonly string[] };
  /** Why it would change nothing. */
  readonly disabled?: string;
}

/** "Show all", then "Hide all", for the layers a menu lists. */
export function layerMenuAll(layers: readonly LayerMenuLayer[]): [LayerMenuAll, LayerMenuAll] {
  const hidden = layers.filter((l) => l.hidden).map((l) => l.address);
  const shown = layers.filter((l) => !l.hidden).map((l) => l.address);
  return [
    {
      label: 'Show all',
      command: 'layer.show',
      args: { addresses: hidden },
      ...(hidden.length ? {} : { disabled: 'Every layer is shown' }),
    },
    {
      label: 'Hide all',
      command: 'layer.hide',
      args: { addresses: shown },
      ...(shown.length ? {} : { disabled: 'Every layer is hidden' }),
    },
  ];
}
