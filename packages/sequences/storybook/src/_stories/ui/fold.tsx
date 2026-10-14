/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReactNode } from "react";
import styles from "./fold.module.css";

/**
 * What a head folds: always there, its height eased between nothing and all of it (a grid row from 0fr to 1fr, so
 * no measuring), out of reach (inert) while folded.
 */
export function Fold({ open, className, children }: { open: boolean; className?: string; children: ReactNode }) {
  return (
    <div className={styles.fold} data-open={open || undefined} inert={!open}>
      <div className={styles.inner}>
        <div className={className}>{children}</div>
      </div>
    </div>
  );
}
