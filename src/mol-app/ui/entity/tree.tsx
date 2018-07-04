/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'

import { View } from '../view';
import { EntityTreeController } from '../../controller/entity/tree';
import { Controller } from '../../controller/controller';
import { AnyEntity, RootEntity } from 'mol-view/state/entity';
import { AnyTransform, SpacefillUpdate, UrlToData, DataToCif, FileToData, CifToMmcif, MmcifToModel, ModelToStructure, StructureToSpacefill, MmcifFileToSpacefill, StructureCenter, StructureToBallAndStick, DistanceRestraintUpdate, CartoonUpdate, BallAndStickUpdate, BackboneUpdate, MmcifUrlToSpacefill } from 'mol-view/state/transform';

function getTransforms(entity: AnyEntity): AnyTransform[] {
    const transforms: AnyTransform[] = []
    switch (entity.kind) {
        case 'root':
            transforms.push(MmcifFileToSpacefill, MmcifUrlToSpacefill)
            break;
        case 'url':
            transforms.push(UrlToData)
            break;
        case 'file':
            transforms.push(FileToData)
            break;
        case 'data':
            transforms.push(DataToCif)
            break;
        case 'cif':
            transforms.push(CifToMmcif)
            break;
        case 'mmcif':
            transforms.push(MmcifToModel)
            break;
        case 'model':
            transforms.push(ModelToStructure)
            break;
        case 'structure':
            transforms.push(StructureToSpacefill, StructureToBallAndStick, StructureCenter)
            break;
        case 'spacefill':
            transforms.push(SpacefillUpdate)
            break;
        case 'ballandstick':
            transforms.push(BallAndStickUpdate)
            break;
        case 'distancerestraint':
            transforms.push(DistanceRestraintUpdate)
            break;
        case 'backbone':
            transforms.push(BackboneUpdate)
            break;
        case 'cartoon':
            transforms.push(CartoonUpdate)
            break;
    }
    return transforms
}

export class Entity extends View<Controller<any>, {}, { entity: AnyEntity}> {
    render() {
        const entity = this.props.entity

        return <div className='molstar-entity-tree-entry'>
            <div className='molstar-entity-tree-entry-body'>
                <div className='molstar-entity-tree-entry-label-wrap'>
                    <button
                        className='molstar-entity-tree-entry-label'
                        onClick={() => {
                            console.log(entity)
                            this.controller.context.currentEntity.next(entity)
                            this.controller.context.currentTransforms.next(getTransforms(entity))
                        }}
                    >
                        <span>{entity.id} - {entity.kind}</span>
                    </button>
                </div>
            </div>
        </div>;
    }
}

export class EntityTree extends View<EntityTreeController, {}, {}> {
    render() {
        const entities: JSX.Element[] = []
        const state = this.controller.state.getValue()
        if (state) {
            state.entities.forEach(e => {
                entities.push(
                    <div key={e.id}>
                        <Entity controller={this.controller} entity={e}></Entity>
                    </div>
                )
            })
        }

        return <div className='molstar-entity-tree'>
            <div className='molstar-entity-tree-root'>
                <Entity controller={this.controller} entity={RootEntity}></Entity>
            </div>
            <div className='molstar-entity-tree-children'>
                <div>{entities}</div>
            </div>
        </div>;
    }
}