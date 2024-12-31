/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import type { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { type PluginContext } from '../../mol-plugin/context';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { ColorThemeCategory } from './categories';
import { StateSelection } from '../../mol-state/state/selection';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { QueryContext, Structure, StructureElement, StructureSelection } from '../../mol-model/structure';
import { ChainIdColorTheme, ChainIdColorThemeParams } from './chain-id';
import { EntityIdColorTheme, EntityIdColorThemeParams } from './entity-id';
import { EntitySourceColorTheme, EntitySourceColorThemeParams } from './entity-source';
import { MoleculeTypeColorTheme, MoleculeTypeColorThemeParams } from './molecule-type';
import { ModelIndexColorTheme, ModelIndexColorThemeParams } from './model-index';
import { StructureIndexColorTheme, StructureIndexColorThemeParams } from './structure-index';
import { assertUnreachable } from '../../mol-util/type-helpers';
import { ScaleLegend, TableLegend } from '../../mol-util/legend';
import { StructureLookup3DResultContext } from '../../mol-model/structure/structure/util/lookup3d';
import { StructureSelectionQueries } from '../../mol-plugin-state/helpers/structure-selection-query';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';

const Description = `Assigns a color based on structure property at a given vertex.`;

export const ExternalStructureColorThemeParams = {
    structure: PD.ValueRef<Structure>(
        (ctx: PluginContext) => {
            const structures = ctx.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure)).filter(c => c.obj?.data);
            return structures.map(v => [v.transform.ref, v.obj?.label ?? '<unknown>'] as [string, string]);
        },
        (ref, getData) => getData(ref),
    ),
    style: PD.MappedStatic('chain-id', {
        'chain-id': PD.Group(ChainIdColorThemeParams),
        'entity-id': PD.Group(EntityIdColorThemeParams),
        'entity-source': PD.Group(EntitySourceColorThemeParams),
        'molecule-type': PD.Group(MoleculeTypeColorThemeParams),
        'model-index': PD.Group(ModelIndexColorThemeParams),
        'structure-index': PD.Group(StructureIndexColorThemeParams),
    }),
    defaultColor: PD.Color(Color(0xcccccc)),
    maxDistance: PD.Numeric(8, { min: 0.1, max: 24, step: 0.1 }, { description: 'Maximum distance to search for the nearest structure element. This is done only if the approximate search fails.' }),
    approxMaxDistance: PD.Numeric(4, { min: 0, max: 12, step: 0.1 }, { description: 'Maximum distance to search for an approximately nearest structure element. This is done before the extact search.' }),
    normalOffset: PD.Numeric(0, { min: -10, max: 20, step: 0.1 }, { description: 'Offset vertex position along its normal by given amount.' }),
    backboneOnly: PD.Boolean(false),
};
export type ExternalStructureColorThemeParams = typeof ExternalStructureColorThemeParams

type ExternalStructureColorThemeProps = PD.Values<ExternalStructureColorThemeParams>

function getStyleTheme(ctx: ThemeDataContext, props: ExternalStructureColorThemeProps['style']) {
    switch (props.name) {
        case 'chain-id': return ChainIdColorTheme(ctx, props.params);
        case 'entity-id': return EntityIdColorTheme(ctx, props.params);
        case 'entity-source': return EntitySourceColorTheme(ctx, props.params);
        case 'molecule-type': return MoleculeTypeColorTheme(ctx, props.params);
        case 'model-index': return ModelIndexColorTheme(ctx, props.params);
        case 'structure-index': return StructureIndexColorTheme(ctx, props.params);
        default: assertUnreachable(props);
    }
}

export function ExternalStructureColorTheme(ctx: ThemeDataContext, props: PD.Values<ExternalStructureColorThemeParams>): ColorTheme<ExternalStructureColorThemeParams> {
    let structure: Structure | undefined;
    try {
        structure = props.structure.getValue();
    } catch {
        // .getValue() is resolved during state reconciliation => would throw from UI
    }

    // NOTE: this will currently be slow for with GPU/texture meshes due to slow iteration
    // TODO: create texture to be able to do the sampling on the GPU

    let color: LocationColor;
    let contextHash: number | undefined = undefined;
    let legend: Readonly<ScaleLegend | TableLegend> | undefined = undefined;

    const { maxDistance, approxMaxDistance: approxDistance, normalOffset, defaultColor, backboneOnly } = props;

    if (structure) {
        const styleTheme = getStyleTheme({ ...ctx, structure }, props.style);

        const lookupFirstCtx = StructureLookup3DResultContext();
        const lookupNearestCtx = StructureLookup3DResultContext();

        let s = structure;
        if (backboneOnly) {
            s = StructureSelection.unionStructure(StructureSelectionQueries.backbone.query(new QueryContext(structure)));
        }

        const position = Vec3();
        const l = StructureElement.Location.create(s);

        color = (location: Location, isSecondary: boolean): Color => {
            if (!isPositionLocation(location)) {
                return defaultColor;
            }

            // Offset the vertex position along its normal
            if (normalOffset !== 0) {
                Vec3.scaleAndAdd(position, location.position, location.normal, normalOffset);
            } else {
                Vec3.copy(position, location.position);
            }

            const [x, y, z] = position;

            if (approxDistance > 0) {
                const rf = s.lookup3d.approxNearest(x, y, z, approxDistance, lookupFirstCtx);
                if (rf.count > 0) {
                    l.unit = rf.units[0];
                    l.element = l.unit.elements[rf.indices[0]];
                    return styleTheme.color(l, isSecondary);
                }
            }

            const rn = s.lookup3d.find(x, y, z, maxDistance, lookupNearestCtx);
            if (rn.count > 0) {
                let idx = 0;
                let minD = rn.squaredDistances[0];
                for (let i = 1; i < rn.count; ++i) {
                    if (rn.squaredDistances[i] < minD) {
                        minD = rn.squaredDistances[i];
                        idx = i;
                    }
                }
                l.unit = rn.units[idx];
                l.element = l.unit.elements[rn.indices[idx]];
                return styleTheme.color(l, isSecondary);
            }

            return defaultColor;
        };

        contextHash = styleTheme.contextHash;
        legend = styleTheme.legend;
    } else {
        color = () => defaultColor;
    }

    return {
        factory: ExternalStructureColorTheme,
        granularity: 'vertex',
        preferSmoothing: true,
        color,
        props,
        contextHash,
        description: Description,
        legend,
    };
}

export const ExternalStructureColorThemeProvider: ColorTheme.Provider<ExternalStructureColorThemeParams, 'external-structure'> = {
    name: 'external-structure',
    label: 'External Structure',
    category: ColorThemeCategory.Misc,
    factory: ExternalStructureColorTheme,
    getParams: () => ExternalStructureColorThemeParams,
    defaultValues: PD.getDefaultValues(ExternalStructureColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true,
};
