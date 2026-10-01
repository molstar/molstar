/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Styles shared by the mask creator and segmentor panels.
 */

import * as React from 'react';

export function SectionDivider({ label }: { label: string }) {
    return (
        <div style={{ display: 'flex', alignItems: 'center', gap: '8px', margin: '14px 0 6px' }}>
            <span style={{ color: '#999', fontSize: '12px', fontWeight: 600, textTransform: 'uppercase', letterSpacing: '0.08em', whiteSpace: 'nowrap' }}>{label}</span>
            <div style={{ flex: 1, height: '1px', background: '#333' }} />
        </div>
    );
}

export const labelStyle: React.CSSProperties = { display: 'block', color: '#888', marginTop: '4px', marginBottom: '2px', fontSize: '11px' };
export const hintStyle: React.CSSProperties = { color: '#666', fontSize: '11px', margin: '4px 0 0', fontStyle: 'italic' };
export const numInputStyle: React.CSSProperties = { width: '60px', background: '#252525', border: '1px solid #444', color: '#ddd', borderRadius: '2px', padding: '1px 4px', fontSize: '12px' };
export const btnStyle: React.CSSProperties = { border: 'none', color: '#ddd', padding: '6px 8px', cursor: 'pointer', borderRadius: '3px', fontSize: '12px' };
export const smallBtnStyle: React.CSSProperties = { ...btnStyle, padding: '2px 7px', fontSize: '11px' };
export const rowStyle: React.CSSProperties = { display: 'flex', justifyContent: 'space-between', alignItems: 'center', background: '#1e1e1e', padding: '3px 6px', marginBottom: '2px', borderRadius: '2px' };
