/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from 'mol-plugin/context';
import { PluginBehavior } from '../behavior';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { degToRad } from 'mol-math/misc';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { PluginStateObject as SO, PluginStateObject } from '../../state/objects';
import { StateSelection } from 'mol-state/state/selection';

// TODO this is just for testing purposes
export const Animation = PluginBehavior.create<{ play: boolean }>({
    name: 'animation',
    display: { name: 'Animation', group: 'Animation' },
    ctor: class extends PluginBehavior.Handler<{ play: boolean }> {
        private tmpMat = Mat4.identity()
        private rotMat = Mat4.identity()
        private transMat = Mat4.identity()
        private animMat = Mat4.identity()
        private transVec = Vec3.zero()
        private rotVec = Vec3.create(0, 1, 0)
        private animHandle = -1

        constructor(protected ctx: PluginContext, protected params: { play: boolean }) {
            super(ctx, params)
            this.update(params)
        }

        animate(play: boolean) {
            if (play) {
                const state = this.ctx.state.dataState
                const reprs = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Representation3D));
                const anim = (t: number) => {
                    const rad = degToRad((t / 10) % 360)
                    Mat4.rotate(this.rotMat, this.tmpMat, rad, this.rotVec)
                    for (const r of reprs) {
                        if (!SO.isRepresentation3D(r.obj)) return
                        const parent = StateSelection.findAncestorOfType(state.tree, state.cells, r.transform.ref, [PluginStateObject.Molecule.Structure])
                        if (!parent || !parent.obj) continue
                        const structure = parent.obj as PluginStateObject.Molecule.Structure

                        Vec3.negate(this.transVec, Vec3.copy(this.transVec, structure.data.boundary.sphere.center))
                        Mat4.fromTranslation(this.transMat, this.transVec)
                        Mat4.mul(this.animMat, this.rotMat, this.transMat)

                        Vec3.copy(this.transVec, structure.data.boundary.sphere.center)
                        Mat4.fromTranslation(this.transMat, this.transVec)
                        Mat4.mul(this.animMat, this.transMat, this.animMat)

                        r.obj.data.setState({ transform: this.animMat })
                        this.ctx.canvas3d.add(r.obj.data)
                        this.ctx.canvas3d.requestDraw(true)
                    }
                    this.animHandle = requestAnimationFrame(anim)
                }
                this.animHandle = requestAnimationFrame(anim)
            } else {
                cancelAnimationFrame(this.animHandle)
            }
        }

        register(): void { }

        update(p: { play: boolean }) {
            let updated = this.params.play !== p.play
            this.params.play = p.play
            if (updated) {
                this.animate(this.params.play)
            }
            return updated;
        }

        unregister() {
            cancelAnimationFrame(this.animHandle)
        }
    },
    params: () => ({
        play: PD.Boolean(false)
    })
});