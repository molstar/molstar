/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from '../../data/issues.js';
import type { Issue } from '@molstar/core/repr/tree';
import {
  abortError,
  type Coverage,
  type FailureKind,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
} from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';

// What a genome runtime answers (_spec/README.md "Data", Genome runtimes' internals): a region's failure by its kind, a
// source refused at open, and a source over what a whole read holds. Shared by the BED, GFF and VCF runtimes, not published
// (`./runtimes/internal/*`).

/** A region's failure of `kind`, its issue the source's code of that kind. */
export function failure(kind: FailureKind, message: string): RegionAnswer & { ok: false } {
  const code =
    kind === 'not-served'
      ? 'seq.source.not_served'
      : kind === 'malformed'
        ? 'seq.source.malformed'
        : kind === 'too-large'
          ? 'seq.source.too_large'
          : 'seq.source.failed';
  return { ok: false, kind, issue: issueOf(SourceIssues, code, message) };
}

/** Refuses a source at open (`SourceOpenError`) with one issue of the source's codes. */
export function refuse(
  address: string,
  code: 'seq.source.unreadable' | 'seq.source.malformed' | 'seq.source.too_large' | 'seq.source.no_data',
  message: string,
  hint?: string,
): never {
  throw new SourceOpenError([issueOf(SourceIssues, code, message, { address, ...(hint && { hint }) })]);
}

/** A source over what a whole read holds: a request answered from it, by system and level. */
export function wholeSource(
  coverage: Coverage,
  issues: Issue[],
  answer: (r: RegionRequest) => RegionAnswer,
): OpenedSource {
  return {
    coverage,
    issues,
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return Promise.resolve(
        !coverage.systems.some((s) => s.name === request.system)
          ? failure('not-served', `No sequence ${request.system}`)
          : !coverage.levels.includes(request.level)
            ? failure('not-served', `No level ${request.level}`)
            : answer(request),
      );
    },
    dispose() {},
  };
}
