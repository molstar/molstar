/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { cdrRows, Data, fab, type Numbered, theDesign, wildType } from './data.js';
import type { Chain } from './numbering.js';

// The numbering schemes side by side (EC-2 to EC-8): a section a chain, a Kabat ruler and an IMGT ruler over its
// letters (each `ruler` layer's `numbering`; insertion labels such as 52A and 100C are the rulers' own), then the CDRs
// by each definition. Over the wild type (its numberings read from the files) for Kabat's insertion codes; over the
// design (inline) for IMGT's 111.1 and 112.1, out of index order (EC-7). The editing view switches which scheme a
// ruler shows: one patch of its layer's `numbering`.

export const title = 'Antibody: Kabat and IMGT side by side';

export const Chains: readonly Chain[] = ['light', 'heavy'];

export interface SchemesOptions {
  /** The instance: the wild type (read from the files) or the design. */
  readonly instance?: 'wild-type' | 'design';
  readonly title?: string;
}

/** The ruler track's layer of a section, in the scheme a patch switches (`/views/main/<chain>/ruler/ruler`). */
export const RulerLayer = (chain: Chain) => ['views', 'main', chain, 'ruler', 'ruler'];

export function build(options: SchemesOptions = {}, d: Derived = Data): Node {
  const b = createSeqSpec({ title: options.title ?? title });
  const numbered: Readonly<Record<Chain, Numbered>> = options.instance === 'design' ? theDesign(d) : wildType(d);
  const seqs = fab(b, options.instance === 'design' ? numbered : undefined);
  const cdrs = Object.fromEntries(
    Chains.flatMap((c) =>
      (['kabat', 'imgt'] as const).map((s) => [
        `${c}-${s}`,
        b.data.ranges(`${c}-cdrs-${s}`, { target: seqs[c], rows: cdrRows(c, s, numbered[c][s]) }),
      ]),
    ),
  );
  b.views.view('main', (v) => {
    for (const c of Chains) {
      const s = v.section(c, { of: seqs[c] });
      s.track('ruler', { label: 'Kabat' }).layer('ruler', { type: 'ruler', numbering: 'kabat' });
      s.track('imgt', { label: 'IMGT' }).layer('ruler', { type: 'ruler', numbering: 'imgt' });
      s.track('letters', { label: c === 'light' ? 'Light chain (A)' : 'Heavy chain (B)' }).layer('letters', {
        type: 'letters',
      });
      s.track('cdrs-kabat', { label: 'CDRs (Kabat)' }).layer('cdrs', { type: 'ranges', data: cdrs[`${c}-kabat`] });
      s.track('cdrs-imgt', { label: 'CDRs (IMGT)' }).layer('cdrs', { type: 'ranges', data: cdrs[`${c}-imgt`] });
    }
  });
  return b.build();
}

/** The design's numberings (EC-7: IMGT's 111.1 and 112.1 in CDR-H3). */
export const design = () => build({ instance: 'design', title: 'Antibody: the design’s Kabat and IMGT' });

/** The CDRs to bookmark: the wild type with both numberings, a CDR selected by its address and bookmarked. */
export const cdrs = () => build({ title: 'Antibody: CDR bookmarks' });
