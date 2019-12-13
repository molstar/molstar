/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model, Structure, StructureSymmetry } from '../../../mol-model/structure';
import { stringToWords } from '../../../mol-util/string';
import { SpacegroupCell } from '../../../mol-math/geometry';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { RuntimeContext } from '../../../mol-task';
import { PluginContext } from '../../context';
import { Assembly, ModelSymmetry } from '../../../mol-model/structure/model/properties/symmetry';
import { PluginStateObject as SO } from '../objects';

export namespace ModelStructureRepresentation {
    export function getParams(model?: Model, defaultValue?: 'deposited' | 'assembly' | 'symmetry' | 'symmetry-mates') {
        const assemblyIds = model ? model.symmetry.assemblies.map(a => [a.id, `${a.id}: ${stringToWords(a.details)}`] as [string, string]) : [];
        const showSymm = !model ? true : !SpacegroupCell.isZero(model.symmetry.spacegroup.cell);

        const modes = {
            deposited: PD.EmptyGroup(),
            assembly: PD.Group({
                id: PD.Optional(model
                    ? PD.Select(assemblyIds.length ? assemblyIds[0][0] : '', assemblyIds, { label: 'Asm Id', description: 'Assembly Id' })
                    : PD.Text('', { label: 'Asm Id', description: 'Assembly Id (use empty for the 1st assembly)' }))
            }, { isFlat: true }),
            'symmetry-mates': PD.Group({
                radius: PD.Numeric(5)
            }, { isFlat: true }),
            'symmetry': PD.Group({
                ijkMin: PD.Vec3(Vec3.create(-1, -1, -1), { label: 'Min IJK', fieldLabels: { x: 'I', y: 'J', z: 'K' } }),
                ijkMax: PD.Vec3(Vec3.create(1, 1, 1), { label: 'Max IJK', fieldLabels: { x: 'I', y: 'J', z: 'K' } })
            }, { isFlat: true })
        };

        const options: [keyof typeof modes, string][] = [
            ['deposited', 'Deposited']
        ];

        if (assemblyIds) {
            options.push(['assembly', 'Assembly']);
        }

        if (showSymm) {
            options.push(['symmetry-mates', 'Symmetry Mates']);
            options.push(['symmetry', 'Symmetry (indices)']);
        }

        return {
            type: PD.MappedStatic(defaultValue || 'deposited', modes, { options })
        };
    }

    export type Params = PD.Values<ReturnType<typeof getParams>>['type']

    async function buildAssembly(plugin: PluginContext, ctx: RuntimeContext, model: Model, id?: string) {
        let asm: Assembly | undefined = void 0;

        // if no id is specified, use the 1st assembly.
        if (!id && model.symmetry.assemblies.length !== 0) {
            id = model.symmetry.assemblies[0].id;
        }

        if (model.symmetry.assemblies.length === 0) {
            if (id !== 'deposited') {
                plugin.log.warn(`Model '${model.entryId}' has no assembly, returning deposited structure.`);
            }
        } else {
            asm = ModelSymmetry.findAssembly(model, id || '');
            if (!asm) {
                plugin.log.warn(`Model '${model.entryId}' has no assembly called '${id}', returning deposited structure.`);
            }
        }

        const base = Structure.ofModel(model);
        if (!asm) {
            const label = { label: 'Deposited', description: Structure.elementDescription(base) };
            return new SO.Molecule.Structure(base, label);
        }

        id = asm.id;
        const s = await StructureSymmetry.buildAssembly(base, id!).runInContext(ctx);
        const props = { label: `Assembly ${id}`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, props);
    }

    async function buildSymmetry(ctx: RuntimeContext, model: Model, ijkMin: Vec3, ijkMax: Vec3) {
        const base = Structure.ofModel(model);
        const s = await StructureSymmetry.buildSymmetryRange(base, ijkMin, ijkMax).runInContext(ctx);
        const props = { label: `Symmetry [${ijkMin}] to [${ijkMax}]`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, props);
    }

    async function buildSymmetryMates(ctx: RuntimeContext, model: Model, radius: number) {
        const base = Structure.ofModel(model);
        const s = await StructureSymmetry.builderSymmetryMates(base, radius).runInContext(ctx);
        const props = { label: `Symmetry Mates`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, props);
    }

    export async function create(plugin: PluginContext, ctx: RuntimeContext, model: Model, params?: Params): Promise<SO.Molecule.Structure> {
        if (!params || params.name === 'deposited') {
            const s = Structure.ofModel(model);
            return new SO.Molecule.Structure(s, { label: 'Deposited', description: Structure.elementDescription(s) });
        }
        if (params.name === 'assembly') {
            return buildAssembly(plugin, ctx, model, params.params.id)
        }
        if (params.name === 'symmetry') {
            return buildSymmetry(ctx, model, params.params.ijkMin, params.params.ijkMax)
        }
        if (params.name === 'symmetry-mates') {
            return buildSymmetryMates(ctx, model, params.params.radius)
        }

        throw new Error(`Unknown represetation type: ${(params as any).name}`);
    }
}