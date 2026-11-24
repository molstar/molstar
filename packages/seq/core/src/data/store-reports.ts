/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { StoreIssues } from './issues.js';

// The data store's reports (_spec/README.md "Data", Store): each cause once per source, at most `max`, then one
// `seq.store.limit_exceeded`; cleared by source.

interface Reported {
  readonly source: string;
  readonly id: string;
  readonly issue: Issue;
}

/** The store-wide reports' source. */
export const StoreWide = '';

export class StoreReports {
  private readonly _reported: Reported[] = [];
  private _capped = false;

  constructor(
    private readonly max: number,
    private readonly publish: (issues: readonly Issue[]) => void,
  ) {}

  report(source: string, issue: Issue): void {
    const id = `${source}\u0000${issue.code}`;
    if (this._reported.some((r) => r.id === id)) return;
    if (this._reported.length >= this.max) {
      if (!this._capped) {
        this._capped = true;
        this._publish();
      }
      return;
    }
    this._reported.push({ source, id, issue });
    this._publish();
  }

  clear(source: string): void {
    const before = this._reported.length;
    for (let i = this._reported.length - 1; i >= 0; i--)
      if (this._reported[i].source === source) this._reported.splice(i, 1);
    if (this._reported.length === before) return;
    this._capped = false;
    this._publish();
  }

  private _publish(): void {
    const list = this._reported.map((r) => r.issue);
    if (this._capped) list.push(issueOf(StoreIssues, 'seq.store.limit_exceeded', `More than ${this.max} issues`));
    this.publish(list);
  }
}
