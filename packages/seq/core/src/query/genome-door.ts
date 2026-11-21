/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome forms' door (_spec/README.md "Queries", Forms): `sequence`, `locus` and `block` are loaded the first time
// a query names one (or a spec holds an assembly: resolve's genome chunk loads them with it). Until they land the form
// waits (incomplete, so its address is `waiting`), never unknown.

import { issueOf } from '@molstar/core/repr/tree';
import type { Resolve } from './form-kit.js';
import type { GenomeForms as Loaded } from './genome-forms.js';
import { QueryIssues } from './issues.js';
import { EmptyLocation } from './location.js';

let forms: typeof Loaded | undefined;
let loading: Promise<unknown> | undefined;
let failed = false;

/**
 * Loads the genome forms, once while in flight; a failed import rejects (resolve's loader says it) and the next call
 * tries again, so the door never keeps a rejected load.
 */
export function loadGenomeForms(): Promise<unknown> {
  if (forms) return Promise.resolve();
  loading ??= import('./genome-forms.js').then(
    (m) => {
      forms = m.GenomeForms;
      failed = false;
      loading = undefined;
    },
    (e: unknown) => {
      failed = true;
      loading = undefined;
      throw e;
    },
  );
  return loading;
}

const door =
  (kind: keyof typeof Loaded): Resolve =>
  (ctx, form, prefix) => {
    if (forms) return forms[kind](ctx, form, prefix);
    // a failed load is the form's error (asked again by the next form); else it waits for the forms
    if (failed) {
      failed = false;
      ctx.issues.push(issueOf(QueryIssues, 'seq.query.form.unsupported', 'The genome forms did not load'));
      return kind === 'sequence' ? { type: 'nodes', nodes: [] } : { type: 'location', location: EmptyLocation };
    }
    loadGenomeForms().catch(() => undefined);
    ctx.complete = false;
    return kind === 'sequence' ? { type: 'nodes', nodes: [] } : { type: 'location', location: EmptyLocation };
  };

export const GenomeForms: Readonly<Record<keyof typeof Loaded, Resolve>> = {
  sequence: door('sequence'),
  locus: door('locus'),
  block: door('block'),
};
