/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import * as React from 'react';
import { useEffect, useRef, useState, useCallback } from 'react';
import type { Point2D, ViewMask } from '../types';

const CLOSE_DIST_PX = 12;
const PREVIEW_COLOR = '#FF6B00';

interface Props {
    /** When set, renders this mask as a preview (pointer-events off, no drawing). */
    previewMask?: ViewMask;
    /** Called when the user successfully closes a new polygon. Only used when previewMask is unset. */
    onPolygonComplete?: (polygon: Point2D[], canvasW: number, canvasH: number) => void;
}

export function DrawingCanvas({ previewMask, onPolygonComplete }: Props) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const [inProgress, setInProgress] = useState<Point2D[]>([]);
    const [mousePos, setMousePos] = useState<Point2D | null>(null);
    const drawing = inProgress.length > 0;

    // Always-current paint function stored in a ref so the resize observer can
    // call it directly without stale closure issues.
    const paintRef = useRef<() => void>(() => {});
    paintRef.current = () => {
        const canvas = canvasRef.current;
        if (!canvas || canvas.width === 0) return;
        const ctx = canvas.getContext('2d');
        if (!ctx) return;
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Preview mode: show selected mask polygon scaled to current canvas
        if (previewMask && inProgress.length === 0) {
            const sx = canvas.width / previewMask.canvasWidth;
            const sy = canvas.height / previewMask.canvasHeight;
            const scaled = previewMask.polygon.map(([x, y]) => [x * sx, y * sy] as Point2D);
            if (previewMask.inverted) {
                drawFilledOutside(ctx, scaled, canvas.width, canvas.height, PREVIEW_COLOR, 0.2);
            } else {
                drawFilledPolygon(ctx, scaled, PREVIEW_COLOR, 0.2);
            }
            drawPolygonOutline(ctx, scaled, PREVIEW_COLOR, 2);
            scaled.forEach(p => drawVertex(ctx, p, PREVIEW_COLOR, 4));
            return;
        }

        // Drawing mode
        if (inProgress.length > 0) {
            drawPolygonOutline(ctx, inProgress, '#FF3300', 2, true);
            inProgress.forEach(p => drawVertex(ctx, p, '#FF3300'));
            if (mousePos) {
                ctx.beginPath();
                ctx.moveTo(inProgress[inProgress.length - 1][0], inProgress[inProgress.length - 1][1]);
                ctx.lineTo(mousePos[0], mousePos[1]);
                ctx.setLineDash([4, 4]);
                ctx.strokeStyle = 'rgba(255,51,0,0.7)';
                ctx.lineWidth = 1.5;
                ctx.stroke();
                ctx.setLineDash([]);
            }
            if (mousePos && inProgress.length >= 3 && distPx(mousePos, inProgress[0]) < CLOSE_DIST_PX) {
                drawVertex(ctx, inProgress[0], '#FF6B00', 8);
            }
        }
    };

    // Size the canvas and hook up the resize observer — runs once on mount.
    // Calls paintRef.current() after each resize so the canvas is always fresh.
    useEffect(() => {
        const canvas = canvasRef.current;
        if (!canvas) return;
        const resize = () => {
            canvas.width = canvas.offsetWidth;
            canvas.height = canvas.offsetHeight;
            paintRef.current();
        };
        const ro = new ResizeObserver(resize);
        ro.observe(canvas);
        resize(); // size + paint immediately
        return () => ro.disconnect();
    }, []);

    // Repaint whenever drawing state or preview changes.
    useEffect(() => { paintRef.current(); }, [inProgress, mousePos, previewMask]);

    const getPos = (e: React.MouseEvent<HTMLCanvasElement>): Point2D => {
        const r = canvasRef.current!.getBoundingClientRect();
        return [e.clientX - r.left, e.clientY - r.top];
    };

    const completePoly = useCallback(() => {
        const canvas = canvasRef.current!;
        onPolygonComplete?.(inProgress, canvas.width, canvas.height);
        setInProgress([]);
        setMousePos(null);
    }, [inProgress, onPolygonComplete]);

    const handleMouseDown = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
        e.stopPropagation();
        e.preventDefault();
        const pos = getPos(e);
        if (inProgress.length === 0) { setInProgress([pos]); return; }
        if (inProgress.length >= 3 && distPx(pos, inProgress[0]) < CLOSE_DIST_PX) { completePoly(); return; }
        setInProgress(prev => [...prev, pos]);
    }, [inProgress, completePoly]);

    const handleMouseMove = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
        if (!drawing) return;
        e.stopPropagation();
        e.preventDefault();
        setMousePos(getPos(e));
    }, [drawing]);

    const handleDblClick = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
        e.stopPropagation();
        e.preventDefault();
        if (inProgress.length >= 3) completePoly();
    }, [inProgress, completePoly]);

    const handleKeyDown = useCallback((e: KeyboardEvent) => {
        if (e.key === 'Escape') { setInProgress([]); setMousePos(null); }
        if ((e.key === 'Enter' || e.key === 'Return') && inProgress.length >= 3) completePoly();
    }, [inProgress, completePoly]);

    useEffect(() => {
        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [handleKeyDown]);

    const isPreviewOnly = !!previewMask && !drawing;

    return (
        <canvas
            ref={canvasRef}
            style={{
                position: 'absolute',
                top: 0, left: 0,
                width: '100%', height: '100%',
                pointerEvents: isPreviewOnly ? 'none' : 'auto',
                cursor: drawing ? 'crosshair' : 'cell',
                zIndex: 10,
            }}
            onMouseDown={isPreviewOnly ? undefined : handleMouseDown}
            onMouseMove={isPreviewOnly ? undefined : handleMouseMove}
            onDoubleClick={isPreviewOnly ? undefined : handleDblClick}
        />
    );
}

// ---- helpers ----------------------------------------------------------------

function distPx(a: Point2D, b: Point2D) {
    const dx = a[0] - b[0], dy = a[1] - b[1];
    return Math.sqrt(dx * dx + dy * dy);
}

function drawFilledPolygon(ctx: CanvasRenderingContext2D, pts: Point2D[], color: string, alpha: number) {
    if (pts.length < 2) return;
    ctx.beginPath();
    ctx.moveTo(pts[0][0], pts[0][1]);
    for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]);
    ctx.closePath();
    const r = parseInt(color.slice(1, 3), 16);
    const g = parseInt(color.slice(3, 5), 16);
    const b = parseInt(color.slice(5, 7), 16);
    ctx.fillStyle = `rgba(${r},${g},${b},${alpha})`;
    ctx.fill();
}

function drawPolygonOutline(ctx: CanvasRenderingContext2D, pts: Point2D[], color: string, width: number, open = false) {
    if (pts.length < 2) return;
    ctx.beginPath();
    ctx.moveTo(pts[0][0], pts[0][1]);
    for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]);
    if (!open) ctx.closePath();
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.stroke();
}

function drawVertex(ctx: CanvasRenderingContext2D, p: Point2D, color: string, r = 4) {
    ctx.beginPath();
    ctx.arc(p[0], p[1], r, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
}

function drawFilledOutside(ctx: CanvasRenderingContext2D, pts: Point2D[], w: number, h: number, color: string, alpha: number) {
    if (pts.length < 2) return;
    const r = parseInt(color.slice(1, 3), 16);
    const g = parseInt(color.slice(3, 5), 16);
    const b = parseInt(color.slice(5, 7), 16);
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, w, h);
    ctx.moveTo(pts[0][0], pts[0][1]);
    for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]);
    ctx.closePath();
    ctx.fillStyle = `rgba(${r},${g},${b},${alpha})`;
    ctx.fill('evenodd');
    ctx.restore();
}
