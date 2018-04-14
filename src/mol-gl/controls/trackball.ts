/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mrdoob/three.js/,
 * copyright (c) 2010-2018 three.js authors. MIT License
 */

import { map, filter, scan } from 'rxjs/operators';

import { Quat, Vec2, Vec3, EPSILON } from 'mol-math/linear-algebra';
import { clamp } from 'mol-math/interpolate';
import InputObserver from 'mol-util/input/input-observer';

export const DefaultTrackballControlsProps = {
    position: Vec3.zero(),
    up: Vec3.create(0, 1, 0),
    target: Vec3.zero(),

    distance: undefined as (number|undefined),
    damping: 0.25,
    rotateSpeed: 0.28,
    zoomSpeed: 0.0075,
    pinchSpeed: 0.0075,
    panSpeed: 1.0,

    staticMoving: false,
    dynamicDampingFactor: 0.2,

    minDistance: 0,
    maxDistance: Infinity
}
export type TrackballControlsProps = Partial<typeof DefaultTrackballControlsProps>

const enum STATE {
    NONE = - 1,
    ROTATE = 0,
    ZOOM = 1,
    PAN = 2,
    TOUCH_ROTATE = 3,
    TOUCH_ZOOM_PAN = 4
}

interface Object {
    position: Vec3,
    up: Vec3,
}

interface Screen {
    left: number
    top: number
    width: number
    height: number
}

interface TrackballControls {
    update: () => void
    copyInto: (positionOut: Vec3, directionOut: Vec3, upOut: Vec3) => void

    position: Vec3
    direction: Vec3
    up: Vec3
    target: Vec3

    distance: number
    damping: number
    rotateSpeed: number
    zoomSpeed: number
    pinchSpeed: number
    panSpeed: number
}

namespace TrackballControls {
    export function create (element: Element, object: Object, props: TrackballControlsProps = {}): TrackballControls {
        const p = { ...DefaultTrackballControlsProps, ...props }

        const screen: Screen = { left: 0, top: 0, width: 0, height: 0 }

        let { rotateSpeed, zoomSpeed, panSpeed } = p
        let { staticMoving, dynamicDampingFactor } = p
        let { minDistance, maxDistance } = p

        // internals
        const target = Vec3.zero()
        const lastPosition = Vec3.zero()

        let _state = STATE.NONE
        let _prevState = STATE.NONE

        const _eye = Vec3.zero()

        const _movePrev = Vec2.zero()
        const _moveCurr = Vec2.zero()

        const _lastAxis = Vec3.zero()
        let _lastAngle = 0

        const _zoomStart = Vec2.zero()
        const _zoomEnd = Vec2.zero()

        let _touchZoomDistanceStart = 0
        let _touchZoomDistanceEnd = 0

        const _panStart = Vec2.zero()
        const _panEnd = Vec2.zero()

        // for reset
        const target0 = Vec3.clone(target)
        const position0 = Vec3.clone(object.position)
        const up0 = Vec3.clone(object.up)

        // methods
        function handleResize () {
            if ( element instanceof Document ) {
                screen.left = 0;
                screen.top = 0;
                screen.width = window.innerWidth;
                screen.height = window.innerHeight;
            } else {
                const box = element.getBoundingClientRect();
                // adjustments come from similar code in the jquery offset() function
                const d = element.ownerDocument.documentElement;
                screen.left = box.left + window.pageXOffset - d.clientLeft;
                screen.top = box.top + window.pageYOffset - d.clientTop;
                screen.width = box.width;
                screen.height = box.height;
            }
        }

        const mouseOnScreenVec2 = Vec2.zero()
        function getMouseOnScreen(pageX: number, pageY: number) {
            Vec2.set(
                mouseOnScreenVec2,
                (pageX - screen.left) / screen.width,
                (pageY - screen.top) / screen.height
            );
            return mouseOnScreenVec2;
        }

        const mouseOnCircleVec2 = Vec2.zero()
        function getMouseOnCircle(pageX: number, pageY: number) {
            Vec2.set(
                mouseOnCircleVec2,
                ((pageX - screen.width * 0.5 - screen.left) / (screen.width * 0.5)),
                ((screen.height + 2 * (screen.top - pageY)) / screen.width) // screen.width intentional
            );
            return mouseOnCircleVec2;
        }

        const rotAxis = Vec3.zero()
        const rotQuat = Quat.zero()
        const rotEyeDir = Vec3.zero()
        const rotObjUpDir = Vec3.zero()
        const rotObjSideDir = Vec3.zero()
        const rotMoveDir = Vec3.zero()

        function rotateCamera() {
            Vec3.set(rotMoveDir, _moveCurr[0] - _movePrev[0], _moveCurr[1] - _movePrev[1], 0);
            let angle = Vec3.magnitude(rotMoveDir);

            if (angle) {
                Vec3.copy(_eye, object.position)
                Vec3.sub(_eye, _eye, target)

                Vec3.normalize(rotEyeDir, Vec3.copy(rotEyeDir, _eye))
                Vec3.normalize(rotObjUpDir, Vec3.copy(rotObjUpDir, object.up))
                Vec3.normalize(rotObjSideDir, Vec3.cross(rotObjSideDir, rotObjSideDir, rotEyeDir))

                Vec3.setMagnitude(rotObjUpDir, rotObjUpDir, _moveCurr[1] - _movePrev[1])
                Vec3.setMagnitude(rotObjSideDir, rotObjSideDir, _moveCurr[0] - _movePrev[0])

                Vec3.add(rotMoveDir, Vec3.copy(rotMoveDir, rotObjUpDir), rotObjSideDir)

                Vec3.normalize(rotAxis, Vec3.cross(rotAxis, rotMoveDir, _eye))

                angle *= rotateSpeed;
                Quat.setAxisAngle(rotQuat, rotAxis, angle )

                Vec3.transformQuat(_eye, _eye, rotQuat)
                Vec3.transformQuat(object.up, object.up, rotQuat)

                Vec3.copy(_lastAxis, rotAxis)
                _lastAngle = angle;
            } else if (!staticMoving && _lastAngle) {

                _lastAngle *= Math.sqrt(1.0 - dynamicDampingFactor);
                Vec3.sub(_eye, Vec3.copy(_eye, object.position), target)
                Quat.setAxisAngle(rotQuat, _lastAxis, _lastAngle)

                Vec3.transformQuat(_eye, _eye, rotQuat)
                Vec3.transformQuat(object.up, object.up, rotQuat)
            }

            Vec2.copy(_movePrev, _moveCurr)
        }


        function zoomCamera () {
            if (_state === STATE.TOUCH_ZOOM_PAN) {
                const factor = _touchZoomDistanceStart / _touchZoomDistanceEnd
                _touchZoomDistanceStart = _touchZoomDistanceEnd;
                Vec3.scale(_eye, _eye, factor)
            } else {
                const factor = 1.0 + ( _zoomEnd[1] - _zoomStart[1] ) * zoomSpeed
                if (factor !== 1.0 && factor > 0.0) {
                    Vec3.scale(_eye, _eye, factor)
                }

                if (staticMoving) {
                    Vec2.copy(_zoomStart, _zoomEnd)
                } else {
                    _zoomStart[1] += ( _zoomEnd[1] - _zoomStart[1] ) * dynamicDampingFactor
                }
            }
        }

        const panMouseChange = Vec2.zero()
        const panObjUp = Vec3.zero()
        const panOffset = Vec3.zero()

        function panCamera() {
            Vec2.sub(panMouseChange, Vec2.copy(panMouseChange, _panEnd), _panStart)

            if ( Vec2.squaredMagnitude(panMouseChange)) {

                mouseChange.multiplyScalar( _eye.length() * _this.panSpeed );

                pan.copy( _eye ).cross( _this.object.up ).setLength( mouseChange.x );
                pan.add( objectUp.copy( _this.object.up ).setLength( mouseChange.y ) );

                _this.object.position.add( pan );
                _this.target.add( pan );

                if ( _this.staticMoving ) {

                    _panStart.copy( _panEnd );

                } else {

                    _panStart.add( mouseChange.subVectors( _panEnd, _panStart ).multiplyScalar( _this.dynamicDampingFactor ) );

                }
            }
        }

        this.checkDistances = function () {

            if ( ! _this.noZoom || ! _this.noPan ) {

                if ( _eye.lengthSq() > _this.maxDistance * _this.maxDistance ) {

                    _this.object.position.addVectors( _this.target, _eye.setLength( _this.maxDistance ) );
                    _zoomStart.copy( _zoomEnd );

                }

                if ( _eye.lengthSq() < _this.minDistance * _this.minDistance ) {

                    _this.object.position.addVectors( _this.target, _eye.setLength( _this.minDistance ) );
                    _zoomStart.copy( _zoomEnd );

                }

            }

        };

        this.update = function () {

            _eye.subVectors( _this.object.position, _this.target );

            if ( ! _this.noRotate ) {

                _this.rotateCamera();

            }

            if ( ! _this.noZoom ) {

                _this.zoomCamera();

            }

            if ( ! _this.noPan ) {

                _this.panCamera();

            }

            _this.object.position.addVectors( _this.target, _eye );

            _this.checkDistances();

            _this.object.lookAt( _this.target );

            if ( lastPosition.distanceToSquared( _this.object.position ) > EPS ) {

                _this.dispatchEvent( changeEvent );

                lastPosition.copy( _this.object.position );

            }

        };

        this.reset = function () {

            _state = STATE.NONE;
            _prevState = STATE.NONE;

            _this.target.copy( _this.target0 );
            _this.object.position.copy( _this.position0 );
            _this.object.up.copy( _this.up0 );

            _eye.subVectors( _this.object.position, _this.target );

            _this.object.lookAt( _this.target );

            _this.dispatchEvent( changeEvent );

            lastPosition.copy( _this.object.position );

        };

        // listeners

        function keydown( event ) {

            if ( _this.enabled === false ) return;

            window.removeEventListener( 'keydown', keydown );

            _prevState = _state;

            if ( _state !== STATE.NONE ) {

                return;

            } else if ( event.keyCode === _this.keys[ STATE.ROTATE ] && ! _this.noRotate ) {

                _state = STATE.ROTATE;

            } else if ( event.keyCode === _this.keys[ STATE.ZOOM ] && ! _this.noZoom ) {

                _state = STATE.ZOOM;

            } else if ( event.keyCode === _this.keys[ STATE.PAN ] && ! _this.noPan ) {

                _state = STATE.PAN;

            }

        }

        function keyup( event ) {

            if ( _this.enabled === false ) return;

            _state = _prevState;

            window.addEventListener( 'keydown', keydown, false );

        }

        function mousedown( event ) {

            if ( _this.enabled === false ) return;

            event.preventDefault();
            event.stopPropagation();

            if ( _state === STATE.NONE ) {

                _state = event.button;

            }

            if ( _state === STATE.ROTATE && ! _this.noRotate ) {

                _moveCurr.copy( getMouseOnCircle( event.pageX, event.pageY ) );
                _movePrev.copy( _moveCurr );

            } else if ( _state === STATE.ZOOM && ! _this.noZoom ) {

                _zoomStart.copy( getMouseOnScreen( event.pageX, event.pageY ) );
                _zoomEnd.copy( _zoomStart );

            } else if ( _state === STATE.PAN && ! _this.noPan ) {

                _panStart.copy( getMouseOnScreen( event.pageX, event.pageY ) );
                _panEnd.copy( _panStart );

            }

            document.addEventListener( 'mousemove', mousemove, false );
            document.addEventListener( 'mouseup', mouseup, false );

            _this.dispatchEvent( startEvent );

        }

        function mousemove( event ) {

            if ( _this.enabled === false ) return;

            event.preventDefault();
            event.stopPropagation();

            if ( _state === STATE.ROTATE && ! _this.noRotate ) {

                _movePrev.copy( _moveCurr );
                _moveCurr.copy( getMouseOnCircle( event.pageX, event.pageY ) );

            } else if ( _state === STATE.ZOOM && ! _this.noZoom ) {

                _zoomEnd.copy( getMouseOnScreen( event.pageX, event.pageY ) );

            } else if ( _state === STATE.PAN && ! _this.noPan ) {

                _panEnd.copy( getMouseOnScreen( event.pageX, event.pageY ) );

            }

        }

        function mouseup( event ) {

            if ( _this.enabled === false ) return;

            event.preventDefault();
            event.stopPropagation();

            _state = STATE.NONE;

            document.removeEventListener( 'mousemove', mousemove );
            document.removeEventListener( 'mouseup', mouseup );
            _this.dispatchEvent( endEvent );

        }

        function mousewheel( event ) {

            if ( _this.enabled === false ) return;

            if ( _this.noZoom === true ) return;

            event.preventDefault();
            event.stopPropagation();

            switch ( event.deltaMode ) {

                case 2:
                    // Zoom in pages
                    _zoomStart.y -= event.deltaY * 0.025;
                    break;

                case 1:
                    // Zoom in lines
                    _zoomStart.y -= event.deltaY * 0.01;
                    break;

                default:
                    // undefined, 0, assume pixels
                    _zoomStart.y -= event.deltaY * 0.00025;
                    break;

            }

            _this.dispatchEvent( startEvent );
            _this.dispatchEvent( endEvent );

        }

        function touchstart( event ) {

            if ( _this.enabled === false ) return;

            switch ( event.touches.length ) {

                case 1:
                    _state = STATE.TOUCH_ROTATE;
                    _moveCurr.copy( getMouseOnCircle( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY ) );
                    _movePrev.copy( _moveCurr );
                    break;

                default: // 2 or more
                    _state = STATE.TOUCH_ZOOM_PAN;
                    var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;
                    var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;
                    _touchZoomDistanceEnd = _touchZoomDistanceStart = Math.sqrt( dx * dx + dy * dy );

                    var x = ( event.touches[ 0 ].pageX + event.touches[ 1 ].pageX ) / 2;
                    var y = ( event.touches[ 0 ].pageY + event.touches[ 1 ].pageY ) / 2;
                    _panStart.copy( getMouseOnScreen( x, y ) );
                    _panEnd.copy( _panStart );
                    break;

            }

            _this.dispatchEvent( startEvent );

        }

        function touchmove( event ) {

            if ( _this.enabled === false ) return;

            event.preventDefault();
            event.stopPropagation();

            switch ( event.touches.length ) {

                case 1:
                    _movePrev.copy( _moveCurr );
                    _moveCurr.copy( getMouseOnCircle( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY ) );
                    break;

                default: // 2 or more
                    var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;
                    var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;
                    _touchZoomDistanceEnd = Math.sqrt( dx * dx + dy * dy );

                    var x = ( event.touches[ 0 ].pageX + event.touches[ 1 ].pageX ) / 2;
                    var y = ( event.touches[ 0 ].pageY + event.touches[ 1 ].pageY ) / 2;
                    _panEnd.copy( getMouseOnScreen( x, y ) );
                    break;

            }

        }

        function touchend( event ) {

            if ( _this.enabled === false ) return;

            switch ( event.touches.length ) {

                case 0:
                    _state = STATE.NONE;
                    break;

                case 1:
                    _state = STATE.TOUCH_ROTATE;
                    _moveCurr.copy( getMouseOnCircle( event.touches[ 0 ].pageX, event.touches[ 0 ].pageY ) );
                    _movePrev.copy( _moveCurr );
                    break;

            }

            _this.dispatchEvent( endEvent );

        }

        function contextmenu( event ) {

            if ( _this.enabled === false ) return;

            event.preventDefault();

        }

        this.dispose = function() {

            this.domElement.removeEventListener( 'contextmenu', contextmenu, false );
            this.domElement.removeEventListener( 'mousedown', mousedown, false );
            this.domElement.removeEventListener( 'wheel', mousewheel, false );

            this.domElement.removeEventListener( 'touchstart', touchstart, false );
            this.domElement.removeEventListener( 'touchend', touchend, false );
            this.domElement.removeEventListener( 'touchmove', touchmove, false );

            document.removeEventListener( 'mousemove', mousemove, false );
            document.removeEventListener( 'mouseup', mouseup, false );

            window.removeEventListener( 'keydown', keydown, false );
            window.removeEventListener( 'keyup', keyup, false );

        };

        this.domElement.addEventListener( 'contextmenu', contextmenu, false );
        this.domElement.addEventListener( 'mousedown', mousedown, false );
        this.domElement.addEventListener( 'wheel', mousewheel, false );

        this.domElement.addEventListener( 'touchstart', touchstart, false );
        this.domElement.addEventListener( 'touchend', touchend, false );
        this.domElement.addEventListener( 'touchmove', touchmove, false );

        window.addEventListener( 'keydown', keydown, false );
        window.addEventListener( 'keyup', keyup, false );

        this.handleResize();

        // force an update at start
        this.update();

    }
}

export default TrackballControls