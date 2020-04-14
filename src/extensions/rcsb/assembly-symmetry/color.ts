/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { AssemblySymmetryProvider, AssemblySymmetry } from './prop';
import { Color } from '../../../mol-util/color';
import { Unit, StructureElement, StructureProperties } from '../../../mol-model/structure';
import { Location } from '../../../mol-model/location';
import { ScaleLegend, TableLegend } from '../../../mol-util/legend';
import { getPalette, getPaletteParams } from '../../../mol-util/color/palette';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';

const DefaultColor = Color(0xCCCCCC);

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id;
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id;
    }
}

function clusterMemberKey(asymId: string, operList: string[]) {
    return `${asymId}-${operList.join('|')}`;
}

export const AssemblySymmetryClusterColorThemeParams = {
    ...getPaletteParams({ colorList: 'red-yellow-blue' }),
};
export type AssemblySymmetryClusterColorThemeParams = typeof AssemblySymmetryClusterColorThemeParams
export function getAssemblySymmetryClusterColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(AssemblySymmetryClusterColorThemeParams);
    return params;
}

export function AssemblySymmetryClusterColorTheme(ctx: ThemeDataContext, props: PD.Values<AssemblySymmetryClusterColorThemeParams>): ColorTheme<AssemblySymmetryClusterColorThemeParams> {
    let color: LocationColor = () => DefaultColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const assemblySymmetry = ctx.structure && AssemblySymmetryProvider.get(ctx.structure);
    const contextHash = assemblySymmetry?.version;

    const clusters = assemblySymmetry?.value?.clusters;

    if (clusters?.length && ctx.structure) {
        const clusterByMember = new Map<string, number>();
        for (let i = 0, il = clusters.length; i < il; ++i) {
            const { members } = clusters[i]!;
            for (let j = 0, jl = members.length; j < jl; ++j) {
                const asymId = members[j]!.asym_id;
                const operList = [...members[j]!.pdbx_struct_oper_list_ids || []] as string[];
                clusterByMember.set(clusterMemberKey(asymId, operList), i);
                if (operList.length === 0) {
                    operList.push('1'); // TODO hack assuming '1' is the id of the identity operator
                    clusterByMember.set(clusterMemberKey(asymId, operList), i);
                }
            }
        }
        const palette = getPalette(clusters.length, props);
        legend = palette.legend;

        const _emptyList: any[] = [];
        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                const { assembly } = location.unit.conformation.operator;
                const asymId = getAsymId(location.unit)(location);
                const cluster = clusterByMember.get(clusterMemberKey(asymId, assembly?.operList || _emptyList));
                return cluster !== undefined ? palette.color(cluster) : DefaultColor;
            }
            return DefaultColor;
        };
    }

    return {
        factory: AssemblySymmetryClusterColorTheme,
        granularity: 'instance',
        color,
        props,
        contextHash,
        description: 'Assigns chain colors according to assembly symmetry cluster membership calculated with BioJava and obtained via RCSB PDB.',
        legend
    };
}

export const AssemblySymmetryClusterColorThemeProvider: ColorTheme.Provider<AssemblySymmetryClusterColorThemeParams, AssemblySymmetry.Tag.Cluster> = {
    name: AssemblySymmetry.Tag.Cluster,
    label: 'Assembly Symmetry Cluster',
    category: ColorTheme.Category.Symmetry,
    factory: AssemblySymmetryClusterColorTheme,
    getParams: getAssemblySymmetryClusterColorThemeParams,
    defaultValues: PD.getDefaultValues(AssemblySymmetryClusterColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => AssemblySymmetry.isApplicable(ctx.structure),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? AssemblySymmetryProvider.attach(ctx, data.structure, void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && data.structure.customPropertyDescriptors.reference(AssemblySymmetryProvider.descriptor, false)
    }
};