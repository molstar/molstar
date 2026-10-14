/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  cloneElement,
  type FocusEvent,
  type PointerEvent,
  type ReactElement,
  type ReactNode,
  useEffect,
  useRef,
  useState,
} from "react";
import { createPortal } from "react-dom";
import styles from "./tooltip.module.css";

/** How long the pointer rests on a trigger before its hint shows, in ms. */
export const TOOLTIP_DELAY_MS = 500;

type Side = "top" | "bottom";

export interface TooltipProps {
  /** What it says: a few words, or a sentence or two. Text only. */
  content: ReactNode;
  /** The element it explains; it keeps its own handlers. */
  children: ReactElement<Record<string, unknown>>;
  /** The side it prefers; it takes the other when this one does not fit in the viewport. */
  side?: Side | "left" | "right";
  /** Where along the trigger it sits: centered, or from its start. */
  align?: "start" | "center" | "end";
}

/** A hint at an anchor's edge, in the page (fixed, so no scroll box clips it); it never takes the pointer. */
function Bubble({
  anchor,
  content,
  side,
  align,
}: { anchor: Element; content: ReactNode } & Pick<TooltipProps, "side" | "align">) {
  const ref = useRef<HTMLDivElement>(null);
  const [at, setAt] = useState<{ left: number; top: number } | null>(null);
  useEffect(() => {
    const bubble = ref.current;
    if (!bubble) return;
    const a = anchor.getBoundingClientRect();
    const b = bubble.getBoundingClientRect();
    const x = align === "start" ? a.left : align === "end" ? a.right - b.width : a.left + (a.width - b.width) / 2;
    const above = a.top - b.height - 6;
    const top =
      side === "bottom"
        ? a.bottom + 6 + b.height < innerHeight
          ? a.bottom + 6
          : above
        : above >= 0
          ? above
          : a.bottom + 6;
    setAt({ left: Math.max(4, Math.min(x, innerWidth - b.width - 4)), top });
  }, [anchor, side, align]);
  return createPortal(
    <div ref={ref} role="tooltip" className={styles.popup} style={at ?? { left: 0, top: 0, visibility: "hidden" }}>
      {content}
    </div>,
    // In the page that carries .msp-root (a story page), so it takes that page's theme and look; else the body.
    anchor.closest(".msp-root") ?? document.body,
  );
}

/**
 * A few words about the thing under the pointer: a hint, never an obstacle. It shows after the pointer rests
 * `TOOLTIP_DELAY_MS` on the trigger, and at once on keyboard focus; it hides on leaving, on a press and a scroll.
 */
export function Tooltip({ content, children, side = "top", align = "center" }: TooltipProps) {
  const [anchor, setAnchor] = useState<Element | null>(null);
  const timer = useRef<ReturnType<typeof setTimeout>>(undefined);
  const hide = () => {
    clearTimeout(timer.current);
    setAnchor(null);
  };
  useEffect(() => () => clearTimeout(timer.current), []);
  useEffect(() => {
    if (!anchor) return;
    const dismiss = () => {
      clearTimeout(timer.current);
      setAnchor(null);
    };
    const o = { capture: true, passive: true };
    window.addEventListener("pointerdown", dismiss, o);
    window.addEventListener("scroll", dismiss, o);
    window.addEventListener("wheel", dismiss, o);
    return () => {
      window.removeEventListener("pointerdown", dismiss, o);
      window.removeEventListener("scroll", dismiss, o);
      window.removeEventListener("wheel", dismiss, o);
    };
  }, [anchor]);
  const own = children.props as {
    onPointerEnter?: (e: PointerEvent) => void;
    onPointerLeave?: (e: PointerEvent) => void;
    onFocus?: (e: FocusEvent) => void;
    onBlur?: (e: FocusEvent) => void;
  };
  return (
    <>
      {cloneElement(children, {
        onPointerEnter: (e: PointerEvent<HTMLElement>) => {
          own.onPointerEnter?.(e);
          const el = e.currentTarget;
          clearTimeout(timer.current);
          timer.current = setTimeout(() => setAnchor(el), TOOLTIP_DELAY_MS);
        },
        onPointerLeave: (e: PointerEvent) => {
          own.onPointerLeave?.(e);
          hide();
        },
        onFocus: (e: FocusEvent<HTMLElement>) => {
          own.onFocus?.(e);
          if (e.currentTarget.matches(":focus-visible")) setAnchor(e.currentTarget);
        },
        onBlur: (e: FocusEvent) => {
          own.onBlur?.(e);
          hide();
        },
      })}
      {anchor && content ? (
        <Bubble anchor={anchor} content={content} side={side === "bottom" ? "bottom" : "top"} align={align} />
      ) : null}
    </>
  );
}

/** The same hint over an element it does not own: shown while there is an anchor. */
export function AnchoredTooltip({ anchor, content }: { anchor: Element | null; content: ReactNode }) {
  return anchor && content ? <Bubble anchor={anchor} content={content} side="top" align="center" /> : null;
}
