/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { clsx } from "clsx";
import styles from "./switch.module.css";

export interface SwitchProps {
  /** What it turns on, beside it (and its accessible name). */
  label: string;
  checked: boolean;
  onCheckedChange: (checked: boolean) => void;
  disabled?: boolean;
  className?: string;
}

/** A small on/off switch with its label. */
export function Switch({ label, checked, onCheckedChange, disabled, className }: SwitchProps) {
  return (
    <label className={clsx(styles.root, disabled && styles.disabled, className)}>
      <span className={styles.label}>{label}</span>
      <input
        type="checkbox"
        role="switch"
        aria-checked={checked}
        className={styles.track}
        checked={checked}
        disabled={disabled}
        onChange={(e) => onCheckedChange(e.currentTarget.checked)}
      />
    </label>
  );
}
