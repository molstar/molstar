/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF from 'mol-io/reader/cif'
import { FileEntity, DataEntity, UrlEntity, CifEntity, MmcifEntity, ModelEntity, StructureEntity, SpacefillEntity, AnyEntity, NullEntity, BallAndStickEntity, DistanceRestraintEntity, CartoonEntity, BackboneEntity } from './entity';
import { Model, Structure, Format } from 'mol-model/structure';

import { StateContext } from './context';
import StructureSymmetry from 'mol-model/structure/structure/symmetry';
import { SpacefillProps, SpacefillRepresentation } from 'mol-geo/representation/structure/spacefill';
import { BallAndStickProps, BallAndStickRepresentation } from 'mol-geo/representation/structure/ball-and-stick';
import { DistanceRestraintRepresentation, DistanceRestraintProps } from 'mol-geo/representation/structure/distance-restraint';
import { CartoonRepresentation, CartoonProps } from 'mol-geo/representation/structure/cartoon';
import { BackboneProps, BackboneRepresentation } from 'mol-geo/representation/structure/backbone';

type transformer<I extends AnyEntity, O extends AnyEntity, P extends {}> = (ctx: StateContext, inputEntity: I, props?: P) => Promise<O>

export interface StateTransform<I extends AnyEntity, O extends AnyEntity, P extends {}> {
    inputKind: I['kind']
    outputKind: O['kind']
    kind: string
    apply: transformer<I, O, P>
}

export namespace StateTransform {
    export function create<I extends AnyEntity, O extends AnyEntity, P extends {}>(inputKind: I['kind'], outputKind: O['kind'], kind: string, transformer: transformer<I, O, P>) {
        return { inputKind, outputKind, kind, apply: transformer }
    }
}

export type AnyTransform = StateTransform<AnyEntity, AnyEntity, {}>

export type UrlToData = StateTransform<UrlEntity, DataEntity, {}>
export const UrlToData: UrlToData = StateTransform.create('url', 'data', 'url-to-data',
    async function (ctx: StateContext, urlEntity: UrlEntity) {
        return DataEntity.ofData(ctx, await urlEntity.value.getData(), urlEntity.value.type)
    })

export type FileToData = StateTransform<FileEntity, DataEntity, {}>
export const FileToData: FileToData = StateTransform.create('file', 'data', 'file-to-data',
    async function (ctx: StateContext, fileEntity: FileEntity) {
        return DataEntity.ofData(ctx, await fileEntity.value.getData(), fileEntity.value.type)
    })

export type DataToCif = StateTransform<DataEntity, CifEntity, {}>
export const DataToCif: DataToCif = StateTransform.create('data', 'cif', 'data-to-cif',
    async function (ctx: StateContext, dataEntity: DataEntity) {
        const comp = CIF.parse(dataEntity.value.data)
        const parsed = await comp.run(ctx.log)
        if (parsed.isError) throw parsed
        return CifEntity.ofCifFile(ctx, parsed.result)
    })

export type CifToMmcif = StateTransform<CifEntity, MmcifEntity, {}>
export const CifToMmcif: CifToMmcif = StateTransform.create('cif', 'mmcif', 'cif-to-mmcif',
    async function (ctx: StateContext, cifEntity: CifEntity) {
        const frame = cifEntity.value.blocks[0];
        return MmcifEntity.ofMmcifDb(ctx, { frame, db: CIF.schema.mmCIF(frame) })
    })

export type MmcifToModel = StateTransform<MmcifEntity, ModelEntity, {}>
export const MmcifToModel: MmcifToModel = StateTransform.create('mmcif', 'model', 'mmcif-to-model',
    async function (ctx: StateContext, mmcifEntity: MmcifEntity) {
        const models = await Model.create(Format.mmCIF(mmcifEntity.value.frame, mmcifEntity.value.db)).run(ctx.log)
        return ModelEntity.ofModels(ctx, models)
    })

export interface StructureProps {
    assembly?: string
}

export type ModelToStructure = StateTransform<ModelEntity, StructureEntity, StructureProps>
export const ModelToStructure: ModelToStructure = StateTransform.create('model', 'structure', 'model-to-structure',
    async function (ctx: StateContext, modelEntity: ModelEntity, props: StructureProps = {}) {
        const model = modelEntity.value[0]
        const assembly = props.assembly
        let structure: Structure
        const assemblies = model.symmetry.assemblies
        if (assemblies.length) {
            structure = await StructureSymmetry.buildAssembly(Structure.ofModel(model), assembly || '1').run(ctx.log)
        } else {
            structure = Structure.ofModel(model)
        }
        return StructureEntity.ofStructure(ctx, structure)
    })

export type StructureCenter = StateTransform<StructureEntity, NullEntity, {}>
export const StructureCenter: StructureCenter = StateTransform.create('structure', 'null', 'structure-center',
    async function (ctx: StateContext, structureEntity: StructureEntity) {
        ctx.viewer.center(structureEntity.value.boundary.sphere.center)
        return NullEntity
    })

export type StructureToSpacefill = StateTransform<StructureEntity, SpacefillEntity, SpacefillProps>
export const StructureToSpacefill: StructureToSpacefill = StateTransform.create('structure', 'spacefill', 'structure-to-spacefill',
    async function (ctx: StateContext, structureEntity: StructureEntity, props: SpacefillProps = {}) {
        const spacefillRepr = SpacefillRepresentation()
        await spacefillRepr.create(structureEntity.value, props).run(ctx.log)
        ctx.viewer.add(spacefillRepr)
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        return SpacefillEntity.ofRepr(ctx, spacefillRepr)
    })

export type StructureToBallAndStick = StateTransform<StructureEntity, BallAndStickEntity, BallAndStickProps>
export const StructureToBallAndStick: StructureToBallAndStick = StateTransform.create('structure', 'ballandstick', 'structure-to-ballandstick',
    async function (ctx: StateContext, structureEntity: StructureEntity, props: BallAndStickProps = {}) {
        const ballAndStickRepr = BallAndStickRepresentation()
        await ballAndStickRepr.create(structureEntity.value, props).run(ctx.log)
        ctx.viewer.add(ballAndStickRepr)
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        return BallAndStickEntity.ofRepr(ctx, ballAndStickRepr)
    })

export type StructureToDistanceRestraint = StateTransform<StructureEntity, DistanceRestraintEntity, DistanceRestraintProps>
export const StructureToDistanceRestraint: StructureToDistanceRestraint = StateTransform.create('structure', 'distancerestraint', 'structure-to-distancerestraint',
    async function (ctx: StateContext, structureEntity: StructureEntity, props: DistanceRestraintProps = {}) {
        const distanceRestraintRepr = DistanceRestraintRepresentation()
        await distanceRestraintRepr.create(structureEntity.value, props).run(ctx.log)
        ctx.viewer.add(distanceRestraintRepr)
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        return DistanceRestraintEntity.ofRepr(ctx, distanceRestraintRepr)
    })
export type StructureToBackbone = StateTransform<StructureEntity, BackboneEntity, BackboneProps>
export const StructureToBackbone: StructureToBackbone = StateTransform.create('structure', 'backbone', 'structure-to-backbone',
            async function (ctx: StateContext, structureEntity: StructureEntity, props: BackboneProps = {}) {
                const backboneRepr = BackboneRepresentation()
                await backboneRepr.create(structureEntity.value, props).run(ctx.log)
                ctx.viewer.add(backboneRepr)
                ctx.viewer.requestDraw()
                console.log('stats', ctx.viewer.stats)
                return BackboneEntity.ofRepr(ctx, backboneRepr)
            })

export type StructureToCartoon = StateTransform<StructureEntity, CartoonEntity, CartoonProps>
export const StructureToCartoon: StructureToCartoon = StateTransform.create('structure', 'cartoon', 'structure-to-cartoon',
        async function (ctx: StateContext, structureEntity: StructureEntity, props: CartoonProps = {}) {
            const cartoonRepr = CartoonRepresentation()
            await cartoonRepr.create(structureEntity.value, props).run(ctx.log)
            ctx.viewer.add(cartoonRepr)
            ctx.viewer.requestDraw()
            console.log('stats', ctx.viewer.stats)
            return CartoonEntity.ofRepr(ctx, cartoonRepr)
        })

export type SpacefillUpdate = StateTransform<SpacefillEntity, NullEntity, SpacefillProps>
export const SpacefillUpdate: SpacefillUpdate = StateTransform.create('spacefill', 'null', 'spacefill-update',
    async function (ctx: StateContext, spacefillEntity: SpacefillEntity, props: SpacefillProps = {}) {
        const spacefillRepr = spacefillEntity.value
        await spacefillRepr.update(props).run(ctx.log)
        ctx.viewer.add(spacefillRepr)
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        return NullEntity
    })

export type BallAndStickUpdate = StateTransform<BallAndStickEntity, NullEntity, BallAndStickProps>
export const BallAndStickUpdate: BallAndStickUpdate = StateTransform.create('ballandstick', 'null', 'ballandstick-update',
    async function (ctx: StateContext, ballAndStickEntity: BallAndStickEntity, props: BallAndStickProps = {}) {
        const ballAndStickRepr = ballAndStickEntity.value
        await ballAndStickRepr.update(props).run(ctx.log)
        ctx.viewer.add(ballAndStickRepr)
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        return NullEntity
    })

export type DistanceRestraintUpdate = StateTransform<DistanceRestraintEntity, NullEntity, DistanceRestraintProps>
export const DistanceRestraintUpdate: DistanceRestraintUpdate = StateTransform.create('distancerestraint', 'null', 'distancerestraint-update',
    async function (ctx: StateContext, distanceRestraintEntity: DistanceRestraintEntity, props: DistanceRestraintProps = {}) {
        const distanceRestraintRepr = distanceRestraintEntity.value
        await distanceRestraintRepr.update(props).run(ctx.log)
        ctx.viewer.add(distanceRestraintRepr)
        ctx.viewer.requestDraw()
        console.log('stats', ctx.viewer.stats)
        return NullEntity
    })

export type BackboneUpdate = StateTransform<BackboneEntity, NullEntity, BackboneProps>
export const BackboneUpdate: BackboneUpdate = StateTransform.create('backbone', 'null', 'backbone-update',
            async function (ctx: StateContext, backboneEntity: BackboneEntity, props: BackboneProps = {}) {
                const backboneRepr = backboneEntity.value
                await backboneRepr.update(props).run(ctx.log)
                ctx.viewer.add(backboneRepr)
                ctx.viewer.requestDraw()
                console.log('stats', ctx.viewer.stats)
                return NullEntity
            })

export type CartoonUpdate = StateTransform<CartoonEntity, NullEntity, CartoonProps>
export const CartoonUpdate: CartoonUpdate = StateTransform.create('cartoon', 'null', 'cartoon-update',
        async function (ctx: StateContext, cartoonEntity: CartoonEntity, props: CartoonProps = {}) {
            const cartoonRepr = cartoonEntity.value
            await cartoonRepr.update(props).run(ctx.log)
            ctx.viewer.add(cartoonRepr)
            ctx.viewer.requestDraw()
            console.log('stats', ctx.viewer.stats)
            return NullEntity
        })

// composed

export type MmcifUrlToModel = StateTransform<UrlEntity, ModelEntity, {}>
export const MmcifUrlToModel: MmcifUrlToModel = StateTransform.create('url', 'model', 'url-to-model',
    async function (ctx: StateContext, urlEntity: UrlEntity) {
        const dataEntity = await UrlToData.apply(ctx, urlEntity)
        return DataToModel.apply(ctx, dataEntity)
    })

export type MmcifFileToModel = StateTransform<FileEntity, ModelEntity, {}>
export const MmcifFileToModel: MmcifFileToModel = StateTransform.create('file', 'model', 'file-to-model',
    async function (ctx: StateContext, fileEntity: FileEntity) {
        const dataEntity = await FileToData.apply(ctx, fileEntity)
        return DataToModel.apply(ctx, dataEntity)
    })

export type DataToModel = StateTransform<DataEntity, ModelEntity, {}>
export const DataToModel: DataToModel = StateTransform.create('data', 'model', 'data-to-model',
    async function getModelFromData(ctx: StateContext, dataEntity: DataEntity) {
        const cifEntity = await DataToCif.apply(ctx, dataEntity)
        const mmcifEntity = await CifToMmcif.apply(ctx, cifEntity)
        return MmcifToModel.apply(ctx, mmcifEntity)
    })

export type ModelToSpacefill = StateTransform<ModelEntity, SpacefillEntity, SpacefillProps>
export const ModelToSpacefill: ModelToSpacefill = StateTransform.create('model', 'spacefill', 'model-to-spacefill',
    async function (ctx: StateContext, modelEntity: ModelEntity, props: SpacefillProps = {}) {
        const structureEntity = await ModelToStructure.apply(ctx, modelEntity)
        // StructureToBond.apply(ctx, structureEntity, props)
        return StructureToSpacefill.apply(ctx, structureEntity, props)
    })

export type MmcifUrlToSpacefill = StateTransform<UrlEntity, SpacefillEntity, SpacefillProps>
export const MmcifUrlToSpacefill: MmcifUrlToSpacefill = StateTransform.create('url', 'spacefill', 'url-to-spacefill',
    async function (ctx: StateContext, urlEntity: UrlEntity, props: SpacefillProps = {}) {
        const modelEntity = await MmcifUrlToModel.apply(ctx, urlEntity)
        return ModelToSpacefill.apply(ctx, modelEntity, props)
    })

export type MmcifFileToSpacefill = StateTransform<FileEntity, SpacefillEntity, SpacefillProps>
export const MmcifFileToSpacefill: MmcifFileToSpacefill = StateTransform.create('file', 'spacefill', 'file-to-spacefill',
    async function (ctx: StateContext, fileEntity: FileEntity, props: SpacefillProps = {}) {
        const modelEntity = await MmcifFileToModel.apply(ctx, fileEntity)
        return ModelToSpacefill.apply(ctx, modelEntity, props)
    })