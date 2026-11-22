/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved view's codes (_spec/README.md "Resolved view", Issues).

import { defineIssues } from '@molstar/core/repr/tree';

const Compose = 'Compose the extension';

export const ViewIssues = defineIssues('seq.view', {
  'seq.view.section.system': {
    severity: 'error',
    description: 'A section’s `of` names no system of the spec; its layers are not placed',
    hint: 'Name a sequence of the spec',
  },
  'seq.view.section.complex': {
    severity: 'warning',
    description: 'A section is `of` a complex its opened source does not hold (an update dropped it); it shows nothing',
    hint: 'Name a complex the source holds',
  },
  'seq.view.block.unit': {
    severity: 'warning',
    description:
      'A `block` over a source’s complex overrides a unit id the complex never held (a mistyped id); it is left out',
    hint: 'Name a unit of the complex by its id',
  },
  'seq.view.section.blocks': {
    severity: 'warning',
    description: 'A section of assemblies has more blocks than `maxAssemblyBlocks`',
    hint: 'Keep fewer by `include: { min_length }`, or the largest first by `order: length`',
  },
  'seq.view.layer.type': {
    severity: 'warning',
    description: 'A layer’s type is not one this version draws, and it has no fallback that is',
    hint: 'Give it a `fallback`',
  },
  'seq.view.layer.extension': {
    severity: 'warning',
    description: 'What an extension this reader lacks (or has too old) draws is drawn as its fallback',
    hint: Compose,
  },
  'seq.view.navigator.picture': {
    severity: 'warning',
    description: 'A navigator’s picture is not registered: drawn as the axis',
    hint: Compose,
  },
  'seq.view.extension.version': {
    severity: 'info',
    description: 'The spec names a newer minor of an extension',
    hint: Compose,
  },
  'seq.view.layer.shape': {
    severity: 'error',
    description: 'A layer’s type does not draw the shape of its data (ranges drawn as values, values as letters)',
    hint: 'Give data of the shape the type draws, or a type that draws it',
  },
  'seq.view.layer.data': {
    severity: 'error',
    description: 'A layer’s data names nothing the static model holds: no sequence, no annotation',
    hint: 'Name a sequence or an annotation of the spec, or a query that starts at one',
  },
  'seq.view.layer.unplaced': {
    severity: 'warning',
    description: 'A layer’s data is on a system no static path joins to its section’s; it is not drawn',
    hint: 'Add a map between the systems, or show the data in a section of its own system',
  },
  'seq.view.via.conflict': {
    severity: 'warning',
    description:
      'Two sections name rival maps of one pair in their `via`: what is carried between them takes the target’s map',
    hint: 'Name the same map in both sections, or drop the `via` from one',
  },
  'seq.view.bridge.neighbours': {
    severity: 'error',
    description: 'A bridge’s sides are not the sections just above and below it in the view; it is not drawn',
    hint: 'Place the bridge between its two sections (three genomes: A, a bridge, B, a bridge, C)',
  },
  'seq.view.bridge.layout': {
    severity: 'error',
    description: 'A bridge’s side is not a row of tracks (a ring, a text, a grid, blocks in a stack); it is not drawn',
    hint: 'Join two sections in the tracks layout, their blocks in a row',
  },
  'seq.view.bridge.layer': {
    severity: 'warning',
    description:
      'A bridge draws pairs: a layer of another type is drawn as its fallbacks’ first `pairs`, else not drawn',
    hint: 'Give the bridge `pairs` layers',
  },
  'seq.view.bridge.lock': {
    severity: 'info',
    description: 'A bridge locks two sections a follow already moves together: the lock moves nothing more',
    hint: 'Drop the lock',
  },
  'seq.view.follow.unplaced': {
    severity: 'warning',
    description: 'A section follows one whose system no static path joins to its own; it does not follow',
    hint: 'Add a map between the sections’ systems',
  },
  'seq.view.section.pieces': {
    severity: 'warning',
    description:
      'A section’s `pieces` reach no piece of its system (`coding` with no codon map, a ranges node on another system); it shows `whole`',
    hint: 'Name a codon map or a ranges node on the section’s system, or follow or read along a protein for `coding`',
  },
  'seq.view.direction.unresolved': {
    severity: 'warning',
    description:
      'A section’s direction cannot be resolved (`follow` with no leader placed, `along` a system no path reaches); it is `genomic`',
    hint: 'Follow a section whose system a path joins, or read along a system that reaches this one',
  },
  'seq.view.direction.mixed': {
    severity: 'info',
    description:
      'The path a section’s direction is read along runs both ways (an inversion); the direction of most of its length is taken',
    hint: 'Spell `genomic` or `reversed` to choose',
  },
});
