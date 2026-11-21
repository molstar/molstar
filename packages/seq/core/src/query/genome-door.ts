/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome forms' door (_spec/README.md "Queries", Forms): `sequence`, `locus` and `block` are loaded the first time
// a query names one (or a spec holds an assembly: resolve's genome chunk loads them with it). Until they land the form
// waits (incomplete, so its address is `waiting`), never unknown.

import type { Resolve } from './form-kit.js';
import type { GenomeForms as Loaded } from './genome-forms.js';
import { EmptyLocation } from './location.js';

let forms: typeof Loaded | undefined;
let loading: Promise<typeof Loaded> | undefined;

/** Loads the genome forms, once. */
export function loadGenomeForms(): Promise<unknown> {
  loading ??= import('./genome-forms.js').then((m) => (forms = m.GenomeForms));
  return loading;
}

const door =
  (kind: keyof typeof Loaded): Resolve =>
  (ctx, form, prefix) => {
    if (forms) return forms[kind](ctx, form, prefix);
    void loadGenomeForms();
    ctx.complete = false;
    return kind === 'sequence' ? { type: 'nodes', nodes: [] } : { type: 'location', location: EmptyLocation };
  };

export const GenomeForms: Readonly<Record<keyof typeof Loaded, Resolve>> = {
  sequence: door('sequence'),
  locus: door('locus'),
  block: door('block'),
};
