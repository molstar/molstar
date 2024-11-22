import { UUID } from '../../../mol-util';
import { ColorNames } from '../../../mol-util/color/names';
import { ControlPoint } from './line-graph-component';

export interface LineGraphParams {
    height: number,
    width: number,
    padding: number,
    paddingYUnnormalized: number,
    paddingXUnnormalized: number,
    // offsetY: number,
    baseline: number
    baselineUnnormalized: number
    // real svg values (400, 600 etc.)
    roof: number,
    // [0; 1]
    roofUnnormalized: number
}


const Baseline = 35 / 2;
const Roof = 35;
const Padding = 70;
const Height = 400;
const Width = 600;

export const LineGraphParams: LineGraphParams = {
    height: Height,
    width: Width,
    padding: Padding,
    // TODO: function that caclulates this based on height etc.
    // or rather contant numbers as const vars above that
    paddingYUnnormalized: Padding / Height,
    paddingXUnnormalized: Padding / Width,
    // TODO: add offset unnormalized or
    // offsetY: 35,
    baseline: Baseline,
    baselineUnnormalized: Baseline / Height,
    roof: Roof,
    // [0; 1]
    roofUnnormalized: Roof / Height
};

export const startEndPoints: ControlPoint[] = [
    // modify this
    {
        data: {
            x: 0,
            alpha: 0,
            // adjustedAlpha: LineGraphParams.baselineUnnormalized,
        },
        id: UUID.create22(),
        color: ColorNames.black,
        index: 0,
        isTerminal: true
    },
    {
        data: {
            x: 1,
            alpha: 0,
            // adjustedAlpha: LineGraphParams.baselineUnnormalized,
        },
        id: UUID.create22(),
        color: ColorNames.black,
        index: 99999999,
        isTerminal: true
    }
];
