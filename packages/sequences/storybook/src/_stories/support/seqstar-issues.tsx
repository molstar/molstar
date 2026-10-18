/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/seq-core/model';
import s from './seqstar-format.module.css';

/** Issues as code, path and message; "none" when empty. */
export function IssueList({ issues }: { issues: readonly Issue[] }) {
  if (issues.length === 0) return <div className={s.none}>none</div>;
  return (
    <ul className={s.issues}>
      {issues.map((i, k) => (
        <li key={k}>
          <span className={s.code}>
            {i.code} {i.path}
          </span>{' '}
          {i.message}
        </li>
      ))}
    </ul>
  );
}
