/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { clsx } from 'clsx';
import type { ReactNode, Ref, UIEventHandler } from 'react';
import styles from './scroll-area.module.css';

export interface ScrollAreaProps {
  children: ReactNode;
  /** On the root: its size and place in the layout. */
  className?: string;
  /** On the viewport, the element that scrolls: padding, and layout of the content. */
  viewportClassName?: string;
  /** The element that scrolls (for scrollTop, scrollIntoView and the like). */
  viewportRef?: Ref<HTMLDivElement>;
  /** Which way it scrolls; vertical by default. */
  orientation?: 'vertical' | 'horizontal' | 'both';
  onScroll?: UIEventHandler<HTMLDivElement>;
}

/** A scrolling region with thin, quiet bars: the root takes its size from the layout, the viewport fills it. */
export function ScrollArea({
  children,
  className,
  viewportClassName,
  viewportRef,
  orientation = 'vertical',
  onScroll,
}: ScrollAreaProps) {
  return (
    <div className={clsx(styles.root, className)}>
      <div
        ref={viewportRef}
        className={clsx(styles.viewport, styles[orientation], viewportClassName)}
        {...(onScroll ? { onScroll } : {})}
      >
        {children}
      </div>
    </div>
  );
}
