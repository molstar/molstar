/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { PluginStateTransform } from '../base';
import { PluginStateObjects as SO } from '../objects';
// import { CartoonRepresentation, DefaultCartoonProps } from 'mol-repr/structure/representation/cartoon';
// import { BallAndStickRepresentation } from 'mol-repr/structure/representation/ball-and-stick';
import { PluginContext } from 'mol-plugin/context';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { RepresentationRegistry } from 'mol-repr/representation';

const colorThemeRegistry = new ColorTheme.Registry()
const sizeThemeRegistry = new SizeTheme.Registry()
const representationRegistry = new RepresentationRegistry()

export { CreateStructureRepresentation }
namespace CreateStructureRepresentation { export interface Params { } }
const CreateStructureRepresentation = PluginStateTransform.Create<SO.Structure, SO.StructureRepresentation3D, CreateStructureRepresentation.Params>({
    name: 'create-structure-representation',
    display: { name: 'Create 3D Representation' },
    from: [SO.Structure],
    to: [SO.StructureRepresentation3D],
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            const repr = representationRegistry.create('cartoon', { colorThemeRegistry, sizeThemeRegistry }, a.data)
            // const repr = BallAndStickRepresentation(); // CartoonRepresentation();
            await repr.createOrUpdate({ webgl: plugin.canvas3d.webgl, colorThemeRegistry, sizeThemeRegistry }, {}, {}, a.data).runInContext(ctx);
            return new SO.StructureRepresentation3D({ label: 'Visual Repr.' }, repr);
        });
    },
    update({ a, b }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            await b.data.createOrUpdate({ webgl: plugin.canvas3d.webgl, colorThemeRegistry, sizeThemeRegistry }, b.data.props, {}, a.data).runInContext(ctx);
            return Transformer.UpdateResult.Updated;
        });
    }
});