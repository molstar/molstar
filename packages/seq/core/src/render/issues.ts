/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The render codes (_spec/README.md "Render", Issues).

import { defineIssues } from '@molstar/core/repr/tree';

export const RenderIssues = defineIssues('seq.render', {
  'seq.render.view': {
    severity: 'error',
    description: 'The spec has no view at the address asked for (or no view at all): nothing is drawn',
    hint: 'Name a view of the spec, or add one under `views`',
  },
  'seq.render.rows.reduced': {
    severity: 'info',
    description:
      'An alignment’s rows are below a pixel and its source serves no coarse row bins: the rows in view are drawn pending',
    hint: 'Zoom in on the rows (`]`), or serve the alignment with row levels',
  },
  'seq.render.logo.profile': {
    severity: 'info',
    description: 'A logo over a lazy alignment: its profile is a source’s to serve, and none is read yet',
    hint: 'Hold the alignment, or draw the logo over a held one',
  },
  'seq.render.rows.order': {
    severity: 'info',
    description: 'A lazy alignment ordered by a field or a tree its source does not serve: drawn in data order',
    hint: 'Order a lazy alignment by `data`, or hold it',
  },
  'seq.render.color': {
    severity: 'warning',
    description: 'A layer’s color names a token the token table lacks',
    hint: 'Name a token of the table (`token:range`), or set it as `--seq-<name>`; the layer’s default color is used meanwhile',
  },
  'seq.render.theme.channel': {
    severity: 'warning',
    description:
      'A theme on a channel the layer type does not read, a value not of the channel’s kind, or a form or palette it cannot read: drawn without it',
    hint: 'Theme a channel the layer reads (its type’s description lists them), with a value of its kind',
  },
  'seq.render.theme.scheme': {
    severity: 'warning',
    description:
      'A scheme that does not fit the layer (another alphabet, a measure not available) or no scheme of the name: the layer’s default is drawn',
    hint: 'Pick a scheme the layer’s menu offers; column schemes need the rows layer’s measures',
  },
  'seq.render.theme.measure': {
    severity: 'warning',
    description: 'A rule reads a measure or a location this layer cannot: its clause does not hold',
    hint: 'Read a field of the layer’s data, `value` or `letter`; a location must resolve on the data’s system',
  },
  'seq.render.theme.color': {
    severity: 'warning',
    description:
      'A color in a theme (a `direct` field’s value, a palette’s color) does not read: it is drawn as missing',
    hint: 'Write colors as `#rrggbb`, a token or `light-dark(a, b)`',
  },
  'seq.render.scale.log': {
    severity: 'info',
    description: 'A value at or below zero on a log scale: drawn as missing (or clamped, with `non_positive: clamp`)',
    hint: 'Set the scale’s `non_positive` to `clamp` to draw them at the low end',
  },
  'seq.render.legend.unlisted': {
    severity: 'info',
    description:
      'A legend cannot list hashed values it has not painted: the values in view are listed, the others take their hashed colors',
    hint: 'List the values in the palette’s `values` (the builder does for inline data) for a complete key',
  },
  'seq.render.line-style': {
    severity: 'info',
    description: 'A ranges layer’s `line` is not a style this version draws: its items are drawn solid',
    hint: 'Use `solid`, `dashed`, `dotted`, `wavy`, `double`, `box` or `highlight`',
  },
  'seq.render.filter': {
    severity: 'warning',
    description: 'A layer’s data query is not a filter its rows can be read by: every row is drawn',
    hint: 'Filter an annotation’s items by `(item :<field> v)`, `(not :of …)`, `(and :of […])` or `(or :of […])`',
  },
  'seq.render.numbering': {
    severity: 'warning',
    description: 'A ruler’s `numbering` is not one its system has: it labels by `ordinal`',
    hint: 'Name `index`, `ordinal` or a numbering the system declares',
  },
  'seq.render.layer.type': {
    severity: 'warning',
    description: 'No type in a layer’s fallback chain is one this target draws: the layer is not drawn',
    hint: 'Give the layer a fallback this target draws',
  },
  'seq.render.layer.letters': {
    severity: 'warning',
    description: 'A letters layer shows another system’s letters through a path: letters are drawn on their own system',
    hint: 'Show the letters in a section of their own system',
  },
  'seq.render.axis.pending': {
    severity: 'info',
    description: 'A section’s system has no known length yet (its source is not open): its tracks wait',
    hint: 'The section draws once its source gives the length',
  },
  'seq.render.text.lines': {
    severity: 'warning',
    description: 'A text section holds more lines than `text.maxLines`: the lines past it are not laid out',
    hint: 'Give the section a `range`: a text is a figure, not a genome',
  },
  'seq.render.follow.leader': {
    severity: 'warning',
    description: 'A section follows one that is not in its view: it shows its own range',
    hint: 'Follow a section of the same view',
  },
  'seq.render.follow.blocks': {
    severity: 'warning',
    description:
      'A section with blocks on more than one system follows or is followed: which block’s map a follow would carry is undecided, so it shows its own range',
    hint: 'Follow one block’s section of its own, or put the blocks on one system',
  },
  'seq.render.ring.image': {
    severity: 'warning',
    description: 'A layer draws an image (a dense matrix) on a ring: images are not drawn on rings',
    hint: 'Show the layer in a tracks section',
  },
  'seq.render.follow.cycle': {
    severity: 'warning',
    description: 'Sections follow each other in a cycle: the first one met shows its own range',
    hint: 'End the chain of follows at a section that follows none',
  },
  'seq.render.source': {
    severity: 'error',
    description: 'A source could not be opened for the render: its layers stay pending',
    hint: 'Check the source’s location and the assets the render may read',
  },
  'seq.render.template.limit': {
    severity: 'warning',
    description: 'A template selects more nodes than it may make tracks for; the first ones in its order are shown',
    hint: 'Filter the template, or raise maxTemplateTracks in the configuration',
  },
  'seq.render.unsettled': {
    severity: 'warning',
    description: 'The export did not settle before its timeout: what was pending is drawn as pending',
    hint: 'Give the export a longer timeout, or check the sources that did not answer',
  },
});
