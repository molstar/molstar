/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci, EmptyLoci } from '../mol-model/loci';
import { StructureElement, Structure } from '../mol-model/structure';
import { Script } from '../mol-script/script';

export { Transparency };

interface Transparency {
    readonly loci: Loci
    readonly value: number
    readonly variant: Transparency.Variant
}

namespace Transparency {
    export type Variant = 'single' | 'multi'
    export const Empty: Transparency = { loci: EmptyLoci, value: 0, variant: 'single' };

    export function areEqual(tA: Transparency, tB: Transparency) {
        if (tA.value !== tB.value) return false;
        if (tA.variant !== tB.variant) return false;
        if (!Loci.areEqual(tA.loci, tB.loci)) return false;
        return true;
    }

    export function ofScript(script: Script, value: number, variant: Variant, structure: Structure): Transparency {
        return { loci: Script.toLoci(script, structure), value, variant };
    }

    export function ofBundle(bundle: StructureElement.Bundle, value: number, variant: Variant, structure: Structure): Transparency {
        return { loci: StructureElement.Bundle.toLoci(bundle, structure), value, variant };
    }
}