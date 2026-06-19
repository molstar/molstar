/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../../mol-util/color';
import { MtlFile, MtlMaterial } from './schema';

// MTL file format specification: http://www.martinreddy.net/gfx/3d/MTL.spec

/**
 * Parse a MTL file string and return a map from material name to properties.
 *
 * Only `newmtl` and `Kd` (diffuse color) directives are processed.
 * All other directives (Ka, Ks, Ns, Ni, illum, d, map_*, etc.) are ignored.
 */
export function parseMtl(data: string): MtlFile {
    const materials = new Map<string, MtlMaterial>();
    let currentName: string | null = null;
    let currentKd: Color = Color(0x808080); // default grey

    const add = () => {
        if (currentName !== null) {
            materials.set(currentName, { Kd: currentKd });
        }
    };

    const lines = data.split(/\r?\n/);
    for (const raw of lines) {
        const line = raw.trim();
        if (line.length === 0 || line.charCodeAt(0) === 35 /* '#' */) continue;

        if (line.startsWith('newmtl')) {
            add(); // Add previous material
            currentName = line.slice(6).trim();
            currentKd = Color(0x808080); // reset to default
        } else if (line.startsWith('Kd')) {
            const parts = line.slice(2).trim().split(/\s+/);
            if (parts.length >= 3) {
                const r = parseFloat(parts[0]);
                const g = parseFloat(parts[1]);
                const b = parseFloat(parts[2]);
                currentKd = Color.fromNormalizedRgb(r, g, b);
            }
        }
        // All other directives are intentionally ignored
    }

    add(); // Add last material
    return materials;
}
