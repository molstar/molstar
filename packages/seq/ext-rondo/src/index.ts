/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineExtension, type SeqExtension } from '@molstar/seqx-core/extension';
import { type RondoConfig, RondoConfigItems, rondoConfig } from './config.js';
import { CellsType, RondoSquare } from './cells.js';
import { RondoHues } from './hues.js';
import { ribbonsType } from './ribbons.js';

// The genome ring extension (_spec/README.md): `rondo`, headless (the CLI and `createSeqContext` compose it and pay
// nothing for the viewer's parts): the `rondo:ribbons` layer type, the `rondo:hues` scheme, their params and settings,
// and the builder helper `rondoRing`. The viewer's commands and binding are `./viewer`'s (`rondoViewer`).

export { DefaultRondoConfig, type RondoConfig, RondoConfigItems, rondoConfig } from './config.js';
export { CellsLay, type CellsPrep, CellsType, RondoSquare } from './cells.js';
export { hueList, RondoHues } from './hues.js';
export { ColorMode, Kind, RondoLay, type RondoPrep, stemAt } from './lay.js';
export { emphasisHeld, RibbonsType, ribbonsType } from './ribbons.js';
export { ribbonsColor, type RondoKind, type RondoRing, rondoRing, type RondoRingOptions } from './ring.js';
export { announceRondo, basesText, describeRondo, summaryRondo } from './words.js';

/** The extension's version: what a spec's `extensions: { rondo: '0.1' }` is read against. */
export const RondoVersion = '0.1.0';

/** The params `rondo:ribbons` reads (plain, on its layer), checked where the extension is composed. */
export const RondoParams = {
  'rondo:ribbons': {
    width: {
      type: 'string',
      doc: 'What a ribbon is wide by: `count` (the links, the default) or a numeric field summed (`q`).',
    },
    scale: {
      type: 'string',
      enum: ['layer', 'section'],
      doc: 'The widths’ scale: the layer’s own (`layer`, the default) or one over the section’s layers of `section` (widths compare across kinds).',
    },
    label: { type: 'string', doc: 'The layer’s name in its card and legend (a kind of link: `Translocations`).' },
    stems: { type: 'string', enum: ['level', 'none'], doc: 'Stems as trumpets (`level`, the default) or none.' },
    bundle: { type: 'number', min: 0, max: 1, doc: 'How far a ribbon is pulled toward the centre (0.7).' },
    stubs: { type: 'string', enum: ['partner', 'none'], doc: 'Links to blocks not shown: a stub a partner, or none.' },
    opens: { type: 'string', doc: 'The section a double-click opens a ribbon’s, stem’s or stub’s chromosomes in.' },
  },
  'rondo:square': {
    of_layer: {
      type: 'string',
      doc: 'The bound section’s layer whose data the cells count (its id); the first `rondo:ribbons` layer when absent.',
    },
  },
} as const;

/**
 * The extension with its settings (`rondo.*`, _spec/README.md, Configuration): a layer type reads none from the
 * context, and a descriptor with settings is composed, never registered late (`viewer.extend`), so a host's numbers
 * are given here (`rondo` is it with the defaults; `RondoConfigItems` name them).
 */
export function rondoWith(over: Partial<RondoConfig> = {}): SeqExtension {
  return defineExtension({
    id: 'rondo',
    version: RondoVersion,
    requires: '0.10',
    doc: 'A genome ring like Rondo’s: ribbons wide by count between stems, stubs to blocks not shown, 24 hues.',
    registry: [
      { kind: 'layer-type', type: ribbonsType(rondoConfig(over)) },
      { kind: 'layer-type', type: CellsType },
      { kind: 'navigator-picture', picture: RondoSquare },
      { kind: 'scheme', scheme: RondoHues },
    ],
    params: RondoParams as never,
    // its settings, readable by its layer types as the composition sets them (`rondo.*`)
    config: Object.values(RondoConfigItems),
  });
}

/** The genome ring extension, headless. */
export const rondo: SeqExtension = rondoWith();
