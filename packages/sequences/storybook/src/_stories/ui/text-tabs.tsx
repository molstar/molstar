/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type KeyboardEvent, type ReactNode, useRef } from "react";
import styles from "./text-tabs.module.css";
import { Tooltip } from "./tooltip.js";

export interface TextTab {
  value: string;
  label: string;
  count?: number | string;
  /** What the tab holds, when its name alone does not say (its tooltip). */
  hint?: string;
}

export interface TextTabsProps {
  tabs: readonly TextTab[];
  value: string;
  onValueChange: (value: string) => void;
  /** Beside the tabs, at the row's end. */
  end?: ReactNode;
}

/** Quiet text tabs: a row of words, the selected one in ink with the accent underline. Arrow keys move along it. */
export function TextTabs({ tabs, value, onValueChange, end }: TextTabsProps) {
  const list = useRef<HTMLDivElement>(null);
  const onKeyDown = (e: KeyboardEvent) => {
    const step = e.key === "ArrowRight" ? 1 : e.key === "ArrowLeft" ? -1 : 0;
    if (!step) return;
    const at = tabs.findIndex((t) => t.value === value);
    const next = tabs[(at + step + tabs.length) % tabs.length];
    if (!next) return;
    e.preventDefault();
    onValueChange(next.value);
    list.current?.querySelector<HTMLElement>(`[data-value="${CSS.escape(next.value)}"]`)?.focus();
  };
  return (
    <div ref={list} role="tablist" className={styles.list} onKeyDown={onKeyDown}>
      {tabs.map((t) => {
        const count = t.count !== undefined && t.count !== "" ? t.count : undefined;
        const tab = (
          <button
            key={t.value}
            type="button"
            role="tab"
            data-value={t.value}
            aria-selected={t.value === value}
            tabIndex={t.value === value ? 0 : -1}
            className={styles.tab}
            onClick={() => onValueChange(t.value)}
          >
            <span>{t.label}</span>
            {count !== undefined && <span className={styles.count}>{count}</span>}
          </button>
        );
        return t.hint ? (
          <Tooltip key={t.value} content={t.hint}>
            {tab}
          </Tooltip>
        ) : (
          tab
        );
      })}
      {end && <span className={styles.end}>{end}</span>}
    </div>
  );
}
