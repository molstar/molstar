/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolstarSvsEntry } from './molstar/registry.js';
import { DefaultRegistry, DefaultSeqViewerSpec, SeqViewerOptionalBehaviors, SeqViewerSpec } from '@molstar/seq-viewer';
import { PrismRegistry } from './palette.js';

// The Seq* composition the Prism viewer's sequence is made from (sequence.md §1): the one place to trim it. Seq*'s
// default for now, so a SeqViewSpec document Prism opens (`load.svs`) draws in full, as UniProt's feature viewer would
// (ranges, sites, line plots, bars, coverage), streaming included until a trim pays (packages/sequences/todo.md "No slim
// build yet"); plus what Prism adds: the Mol* structure as a source, Prism's palette, and the optional reveal of what
// another view lights.

/** The letters' size in the text version, and their largest. */
const DefaultFontPx = 12;

export function PrismSequenceSpec(o: { readonly textWidth?: number } = {}): SeqViewerSpec {
  const base = DefaultSeqViewerSpec();
  return {
    ...base,
    // What another view lights comes into view (text scrolls its lines in, tracks pan at the zoom): Seq*'s optional
    // `revealHighlight`, following the link channels (`links:*`) and the 3D view's hover, which the link writes as an
    // external hover owned by `3d` (seqspec:interaction#highlight).
    behaviors: [
      ...base.behaviors,
      SeqViewerSpec.Behavior(SeqViewerOptionalBehaviors.RevealHighlight, { channels: ['links:*'], hovers: ['3d'] }),
    ],
    registry: [DefaultRegistry, { svs: MolstarSvsEntry }, { svs: PrismRegistry }],
    // The text version (Seq*'s `compact` preset leaves these to the host): whole blocks per line over the full width,
    // centered, the letters in the layer's colors and no larger than the panel text (set from the host's tokens:
    // `setFonts`); `textWidth` limits the line.
    render: {
      wrappedAlign: 'center',
      wrappedColors: true,
      unobservedGray: true,
      wrappedMaxWidth: o.textWidth ?? 0,
      wrappedFontPx: DefaultFontPx,
      maxLetterPx: DefaultFontPx,
    },
  };
}
