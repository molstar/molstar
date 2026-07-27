/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SpacegroupData } from '../tables';
import { SyminfoEntriesByCcp4Number } from './utils';

describe('SpacegroupData aliases', () => {
    it('includes every syminfo.lib "symbol old" alias name', () => {
        const failures: string[] = [];
        for (const entry of SpacegroupData) {
            if (entry.ccp4Number === 0) continue;
            const syminfoEntry = SyminfoEntriesByCcp4Number.get(entry.ccp4Number);
            if (!syminfoEntry) {
                failures.push(`ccp4 ${entry.ccp4Number}: not present in syminfo.lib`);
                continue;
            }
            const names = new Set(entry.names);
            for (const alias of syminfoEntry.old) {
                if (!names.has(alias)) {
                    failures.push(`ccp4 ${entry.ccp4Number}: missing alias '${alias}' (has [${entry.names.join(', ')}])`);
                }
            }
        }
        expect(failures).toEqual([]);
    });

    // Regression guard: `SpacegroupNameToNumberMap` (`../tables`) only maps
    // names from the ~268 CCP4-numbered settings (see `SpacegroupEntryByName`'s
    // two-pass construction, which deliberately lets a numbered setting's name
    // win over a same-named non-numbered alternate, e.g. shared origin-choice
    // names) - so only check for collisions among numbered entries here.
    it('never assigns the same accepted name to two different CCP4-numbered spacegroups', () => {
        const failures: string[] = [];
        const nameToNumber = new Map<string, number>();
        for (const entry of SpacegroupData) {
            if (entry.ccp4Number === 0) continue;
            for (const name of entry.names) {
                const existing = nameToNumber.get(name);
                if (existing !== undefined && existing !== entry.ccp4Number) {
                    failures.push(`'${name}' claimed by both ${existing} and ${entry.ccp4Number}`);
                } else {
                    nameToNumber.set(name, entry.ccp4Number);
                }
            }
        }
        expect(failures).toEqual([]);
    });
});
