
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
