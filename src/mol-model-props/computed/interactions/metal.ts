/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * based in part on NGL (https://github.com/arose/ngl)
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { FeaturesBuilder, Features } from './features';
import { typeSymbol, compId, atomId } from '../chemistry/util';
import { Elements, isTransitionMetal, isHalogen } from '../../../mol-model/structure/model/properties/atomic/types';
import { FeatureType, FeatureGroup, InteractionType } from './common';
import { ContactProvider } from './contacts';
import { ElementSymbol, AminoAcidNames, BaseNames, ProteinBackboneAtoms, NucleicBackboneAtoms } from '../../../mol-model/structure/model/types';

export const MetalCoordinationParams = {
    distanceMax: PD.Numeric(3.0, { min: 1, max: 5, step: 0.1 }),
};
export type MetalCoordinationParams = typeof MetalCoordinationParams
export type MetalCoordinationProps = PD.Values<MetalCoordinationParams>

const IonicTypeMetals = [
    Elements.LI, Elements.NA, Elements.K, Elements.RB, Elements.CS,
    Elements.MG, Elements.CA, Elements.SR, Elements.BA, Elements.AL,
    Elements.GA, Elements.IN, Elements.TL, Elements.SC, Elements.SN,
    Elements.PB, Elements.BI, Elements.SB, Elements.HG
] as ElementSymbol[];

function addMetal(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        const element = typeSymbol(unit, i);

        let type = FeatureType.None;
        if (IonicTypeMetals.includes(element)) {
            type = FeatureType.IonicTypeMetal;
        } else if (isTransitionMetal(element) || element === Elements.ZN || element === Elements.CD) {
            type = FeatureType.TransitionMetal;
        }

        if (type) {
            builder.add(type, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
        }
    }
}

function isProteinSidechain(atomname: string) {
    return !ProteinBackboneAtoms.has(atomname);
}

function isProteinBackbone(atomname: string) {
    return ProteinBackboneAtoms.has(atomname);
}

function isNucleicBackbone(atomname: string) {
    return NucleicBackboneAtoms.has(atomname);
}

/**
 * Metal binding partners (dative bond or ionic-type interaction)
 */
function addMetalBinding(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        const element = typeSymbol(unit, i);
        const resname = compId(unit, i);
        const atomname = atomId(unit, i);
        let dative = false;
        let ionic = false;

        const isStandardAminoacid = AminoAcidNames.has(resname);
        const isStandardBase = BaseNames.has(resname);

        if (!isStandardAminoacid && !isStandardBase) {
            if (isHalogen(element) || element === Elements.O || element === Elements.S) {
                dative = true;
                ionic = true;
            } else if (element === Elements.N) {
                dative = true;
            }
        } else if (isStandardAminoacid) {
            // main chain oxygen atom or oxygen, nitrogen and sulfur from specific amino acids
            if (element === Elements.O) {
                if (['ASP', 'GLU', 'SER', 'THR', 'TYR', 'ASN', 'GLN'].includes(resname) && isProteinSidechain(atomname)) {
                    dative = true;
                    ionic = true;
                } else if (isProteinBackbone(atomname)) {
                    dative = true;
                    ionic = true;
                }
            } else if (element === Elements.S && (resname === 'CYS' || resname === 'MET')) {
                dative = true;
                ionic = true;
            } else if (element === Elements.N) {
                if (resname === 'HIS' && isProteinSidechain(atomname)) {
                    dative = true;
                }
            }
        } else if (isStandardBase) {
            // http://pubs.acs.org/doi/pdf/10.1021/acs.accounts.6b00253
            // http://onlinelibrary.wiley.com/doi/10.1002/anie.200900399/full
            if (element === Elements.O && isNucleicBackbone(atomname)) {
                dative = true;
                ionic = true;
            } else if (['N3', 'N4', 'N7'].includes(atomname)) {
                dative = true;
            } else if (['O2', 'O4', 'O6'].includes(atomname)) {
                dative = true;
                ionic = true;
            }
        }

        if (dative) {
            builder.add(FeatureType.DativeBondPartner, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
        }
        if (ionic) {
            builder.add(FeatureType.IonicTypePartner, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
        }
    }
}

function isMetalCoordination(ti: FeatureType, tj: FeatureType) {
    if (ti === FeatureType.TransitionMetal) {
        return (
            tj === FeatureType.DativeBondPartner ||
            tj === FeatureType.TransitionMetal
        );
    } else if (ti === FeatureType.IonicTypeMetal) {
        return (
            tj === FeatureType.IonicTypePartner
        );
    }
}

function testMetalCoordination(structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number): InteractionType | undefined {
    const typeA = infoA.types[infoA.feature];
    const typeB = infoB.types[infoB.feature];

    if (!isMetalCoordination(typeA, typeB) && !isMetalCoordination(typeB, typeA)) return;

    return InteractionType.MetalCoordination;
}

//

export const MetalProvider = Features.Provider([FeatureType.IonicTypeMetal, FeatureType.TransitionMetal], addMetal);
export const MetalBindingProvider = Features.Provider([FeatureType.IonicTypePartner, FeatureType.DativeBondPartner], addMetalBinding);

export const MetalCoordinationProvider: ContactProvider<MetalCoordinationParams> = {
    name: 'metal-coordination',
    params: MetalCoordinationParams,
    createTester: (props: MetalCoordinationProps) => {
        return {
            maxDistance: props.distanceMax,
            requiredFeatures: new Set([FeatureType.IonicTypeMetal, FeatureType.TransitionMetal, FeatureType.IonicTypePartner, FeatureType.DativeBondPartner]),
            getType: (structure, infoA, infoB, distanceSq) => testMetalCoordination(structure, infoA, infoB, distanceSq)
        };
    }
};