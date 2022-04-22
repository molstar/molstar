/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, Structure, StructureSymmetry } from '../../mol-model/structure';
import { stringToWords } from '../../mol-util/string';
import { SpacegroupCell, Spacegroup } from '../../mol-math/geometry';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Vec3 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { PluginContext } from '../../mol-plugin/context';
import { Assembly, Symmetry } from '../../mol-model/structure/model/properties/symmetry';
import { PluginStateObject as SO } from '../objects';
import { ModelSymmetry } from '../../mol-model-formats/structure/property/symmetry';
import { assertUnreachable } from '../../mol-util/type-helpers';

const CommonStructureParams = {
    dynamicBonds: PD.Optional(PD.Boolean(false, { description: 'Ensure bonds are recalculated upon model changes. Also enables calculation of inter-unit bonds in water molecules and ions.' })),
};
type CommonStructureProps = PD.ValuesFor<typeof CommonStructureParams>

export namespace RootStructureDefinition {
    export function getParams(model?: Model, defaultValue?: 'auto' | 'model' | 'assembly' | 'symmetry' | 'symmetry-mates' | 'symmetry-assembly') {
        const symmetry = model && ModelSymmetry.Provider.get(model);

        const assemblyIds = symmetry ? symmetry.assemblies.map(a => [a.id, `${a.id}: ${stringToWords(a.details)}`] as [string, string]) : [];
        const showSymm = !symmetry ? true : !SpacegroupCell.isZero(symmetry.spacegroup.cell);

        const operatorOptions: [number, string][] = [];
        if (symmetry) {
            const { operators } = symmetry.spacegroup;
            for (let i = 0, il = operators.length; i < il; i++) {
                operatorOptions.push([i, `${i + 1}: ${Spacegroup.getOperatorXyz(operators[i])}`]);
            }
        }

        const asymIdsOptions: [string, string][] = [];
        if (model) {
            model.properties.structAsymMap.forEach(v => {
                const label = v.id === v.auth_id ? v.id : `${v.id} [auth ${v.auth_id}]`;
                asymIdsOptions.push([v.id, label]);
            });
        }

        const modes = {
            auto: PD.Group(CommonStructureParams),
            model: PD.Group(CommonStructureParams),
            assembly: PD.Group({
                id: PD.Optional(model
                    ? PD.Select(assemblyIds.length ? assemblyIds[0][0] : '', assemblyIds, { label: 'Asm Id', description: 'Assembly Id' })
                    : PD.Text('', { label: 'Asm Id', description: 'Assembly Id (use empty for the 1st assembly)' })),
                ...CommonStructureParams
            }, { isFlat: true }),
            'symmetry-mates': PD.Group({
                radius: PD.Numeric(5, { min: 0, max: 50, step: 1 }),
                ...CommonStructureParams
            }, { isFlat: true }),
            'symmetry': PD.Group({
                ijkMin: PD.Vec3(Vec3.create(-1, -1, -1), { step: 1 }, { label: 'Min IJK', fieldLabels: { x: 'I', y: 'J', z: 'K' } }),
                ijkMax: PD.Vec3(Vec3.create(1, 1, 1), { step: 1 }, { label: 'Max IJK', fieldLabels: { x: 'I', y: 'J', z: 'K' } }),
                ...CommonStructureParams
            }, { isFlat: true }),
            'symmetry-assembly': PD.Group({
                generators: PD.ObjectList({
                    operators: PD.ObjectList({
                        index: PD.Select(0, operatorOptions),
                        shift: PD.Vec3(Vec3(), { step: 1 }, { label: 'IJK', fieldLabels: { x: 'I', y: 'J', z: 'K' } })
                    }, e => `${e.index + 1}_${e.shift.map(a => a + 5).join('')}`, {
                        defaultValue: [] as { index: number, shift: Vec3 }[]
                    }),
                    asymIds: PD.MultiSelect([] as string[], asymIdsOptions)
                }, e => `${e.asymIds.length} asym ids, ${e.operators.length} operators`, {
                    defaultValue: [] as { operators: { index: number, shift: Vec3 }[], asymIds: string[] }[]
                }),
                ...CommonStructureParams
            }, { isFlat: true })
        };

        const options: [keyof typeof modes, string][] = [];

        if (defaultValue === 'auto') {
            options.push(['auto', 'Auto']);
        }

        options.push(['model', 'Model']);

        if (assemblyIds.length > 0) {
            options.push(['assembly', 'Assembly']);
        }

        if (showSymm) {
            options.push(['symmetry-mates', 'Symmetry Mates']);
            options.push(['symmetry', 'Symmetry (indices)']);
            options.push(['symmetry-assembly', 'Symmetry (assembly)']);
        }

        return {
            type: PD.MappedStatic(defaultValue || 'model', modes, { options })
        };
    }

    export type Params = PD.Values<ReturnType<typeof getParams>>['type']

    export function canAutoUpdate(oldParams: Params, newParams: Params) {
        if (newParams.name === 'symmetry-assembly' || (newParams.name === 'symmetry' && oldParams.name === 'symmetry')) return false;
        return true;
    }

    async function buildAssembly(plugin: PluginContext, ctx: RuntimeContext, model: Model, id?: string, props?: CommonStructureProps) {
        let asm: Assembly | undefined = void 0;

        const symmetry = ModelSymmetry.Provider.get(model);

        // if no id is specified, use the 1st assembly.
        if (!id && symmetry && symmetry.assemblies.length !== 0) {
            id = symmetry.assemblies[0].id;
        }

        if (!symmetry || symmetry.assemblies.length === 0) {
            plugin.log.warn(`Model '${model.entryId}' has no assembly, returning model structure.`);
        } else {
            asm = Symmetry.findAssembly(model, id || '');
            if (!asm) {
                plugin.log.warn(`Model '${model.entryId}' has no assembly called '${id}', returning model structure.`);
            }
        }

        const base = Structure.ofModel(model, props);
        if (!asm) {
            const label = { label: 'Model', description: Structure.elementDescription(base) };
            return new SO.Molecule.Structure(base, label);
        }

        id = asm.id;
        const s = await StructureSymmetry.buildAssembly(base, id!).runInContext(ctx);
        const objProps = { label: `Assembly ${id}`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, objProps);
    }

    async function buildSymmetry(ctx: RuntimeContext, model: Model, ijkMin: Vec3, ijkMax: Vec3, props?: CommonStructureProps) {
        const base = Structure.ofModel(model, props);
        const s = await StructureSymmetry.buildSymmetryRange(base, ijkMin, ijkMax).runInContext(ctx);
        const objProps = { label: `Symmetry [${ijkMin}] to [${ijkMax}]`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, objProps);
    }

    async function buildSymmetryMates(ctx: RuntimeContext, model: Model, radius: number, props?: CommonStructureProps) {
        const base = Structure.ofModel(model, props);
        const s = await StructureSymmetry.builderSymmetryMates(base, radius).runInContext(ctx);
        const objProps = { label: `Symmetry Mates`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, objProps);
    }

    async function buildSymmetryAssembly(ctx: RuntimeContext, model: Model, generators: StructureSymmetry.Generators, symmetry: Symmetry, props?: CommonStructureProps) {
        const base = Structure.ofModel(model, props);
        const s = await StructureSymmetry.buildSymmetryAssembly(base, generators, symmetry).runInContext(ctx);
        const objProps = { label: `Symmetry Assembly`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, objProps);
    }

    export async function create(plugin: PluginContext, ctx: RuntimeContext, model: Model, params?: Params): Promise<SO.Molecule.Structure> {
        const props = params?.params;
        const symmetry = ModelSymmetry.Provider.get(model);
        if (!symmetry || !params || params.name === 'model') {
            const s = Structure.ofModel(model, props);
            return new SO.Molecule.Structure(s, { label: 'Model', description: Structure.elementDescription(s) });
        }
        if (params.name === 'auto') {
            if (symmetry.assemblies.length === 0) {
                const s = Structure.ofModel(model, props);
                return new SO.Molecule.Structure(s, { label: 'Model', description: Structure.elementDescription(s) });
            } else {
                return buildAssembly(plugin, ctx, model, undefined, props);
            }
        }
        if (params.name === 'assembly') {
            return buildAssembly(plugin, ctx, model, params.params.id, props);
        }
        if (params.name === 'symmetry') {
            return buildSymmetry(ctx, model, params.params.ijkMin, params.params.ijkMax, props);
        }
        if (params.name === 'symmetry-mates') {
            return buildSymmetryMates(ctx, model, params.params.radius, props);
        }
        if (params.name === 'symmetry-assembly') {
            return buildSymmetryAssembly(ctx, model, params.params.generators, symmetry, props);
        }

        assertUnreachable(params);
    }
}