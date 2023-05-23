/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PresetStructureRepresentations } from "../../mol-plugin-state/builder/structure/representation-preset";
import { Color } from "../../mol-util/color";
import { CollapsableControls, PurePluginUIComponent } from "../base";
import { Button } from "../controls/common";
import { MagicWandSvg } from "../controls/icons";
import { ParamDefinition as PD } from "../../mol-util/param-definition";
import { PostprocessingParams } from "../../mol-canvas3d/passes/postprocessing";
import { PluginConfig } from "../../mol-plugin/config";
import { StructureComponentManager } from "../../mol-plugin-state/manager/structure/component";

export class StructureQuickStylesControls extends CollapsableControls {
  defaultState() {
    return {
      isCollapsed: false,
      header: "Quick Styles",
      brand: { accent: "gray" as const, svg: MagicWandSvg },
    };
  }

  renderControls() {
    return (
      <>
        <QuickStyles />
      </>
    );
  }
}

export class QuickStyles extends PurePluginUIComponent {
  state = {
    previousStyle: "",
    clicked: false,
    previousOptions: this.plugin.managers.structure.component.state.options,
  };

  leaveTimeout: number | null = null;

  // Define a method to revert to previous style
  // async revertStyle() {
  //   // Code to apply style saved in this.state.previousStyle

  //   // Check if there's a previous style saved
  //   if (!this.state.previousStyle) {
  //     this.default();
  //   } else {
  //     // Extract the properties from the previous style
  //     const { structure, postprocessing } = this.state.previousStyle;

  //     // Apply the previous structure style
  //     const { structures } = this.plugin.managers.structure.hierarchy.selection;
  //     await this.plugin.managers.structure.component.applyPreset(
  //       structures,
  //       structure
  //     );

  //     // Apply the previous postprocessing style
  //     if (this.plugin.canvas3d && postprocessing) {
  //       this.plugin.canvas3d.setProps({ postprocessing });
  //     }

  //     // Reset previousStyle state
  //     this.setState({ previousStyle: null });
  //   }
  // }

  // Add getCurrentStyle() method to return current style
  getCurrentStyle() {
    // Return the current style.
    // const structureState = this.plugin.managers.structure.component.state;
    // Similarly, you could get the current postprocessing state like this:
    // const postprocessingState = this.plugin.canvas3d ? this.plugin.canvas3d.props.postprocessing : null;
    // Now return an object that captures the current state
    // return {
    // structure: structureState,
    // postprocessing: postprocessingState,
    // };
  }
  applyStyle = (innerText: string) => {
    console.log(innerText);
    if (innerText.includes("Default")) {
      console.log("Applying Default");
      this.default();
    } else if (innerText.includes("Illustrative")) {
      console.log("Applying Illustrative");
      this.illustrative();
    } else {
      console.log("Applying stylized");
      this.stylized();
    }
  };

  handleMouseClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    console.log("Mouse Clicked");

    if (this.leaveTimeout) {
      // If timer exists, clear it to prevent defaultAfterPreview from executing
      clearTimeout(this.leaveTimeout);
      this.leaveTimeout = null;
    }
    const buttonInnerText = event.currentTarget.innerText;
    this.applyStyle(buttonInnerText);
    this.setState({ clicked: true, previousStyle: buttonInnerText }, () => {
      console.log(this.state.clicked + this.state.previousStyle); // Outputs the updated value
    });
  };

  handleMouseEnter = (event: React.MouseEvent<HTMLButtonElement>) => {
    this.setState({ clicked: false }, () => {
      console.log("Mouse Entered " + this.state.clicked); // Outputs the updated value
    });
    this.setState({ previousOptions: this.plugin.managers.structure.component.state.options }, () => {
        console.log("Mouse Entered " + this.state.previousOptions); // Outputs the updated value
    });
    this.applyStyle(event.currentTarget.innerText);
    // console.log(this.state.clicked);
  };

  handleMouseLeave = () => {
    console.log("Mouse Leave");
    console.log(this.state.previousStyle)
    if (!this.state.clicked && this.state.previousStyle) {
      // Clear the previous timeout if it exists
      console.log("Reapply previous style");
      if (this.leaveTimeout !== null) {
        clearTimeout(this.leaveTimeout);
        this.leaveTimeout = null;
      }

      // Set a timeout to delay the onMouseLeave action
      this.leaveTimeout = window.setTimeout(
        () => this.defaultAfterPreview(),
        500
      );
    }
  };

  // Add a method to apply the default style
  async defaultAfterPreview() {
    // Clear the previous timeout if it exists
    if (this.leaveTimeout) {
      clearTimeout(this.leaveTimeout);
    }

    this.applyStyle(this.state.previousStyle)



  }

  async default() {
    const { structures } = this.plugin.managers.structure.hierarchy.selection;
    const preset =
      this.plugin.config.get(
        PluginConfig.Structure.DefaultRepresentationPreset
      ) || PresetStructureRepresentations.auto.id;
    const provider =
      this.plugin.builders.structure.representation.resolveProvider(preset);
    await this.plugin.managers.structure.component.applyPreset(
      structures,
      provider
    );

    this.plugin.managers.structure.component.setOptions(
      PD.getDefaultValues(StructureComponentManager.OptionsParams)
    );

    if (this.plugin.canvas3d) {
      const p = PD.getDefaultValues(PostprocessingParams);
      this.plugin.canvas3d.setProps({
        postprocessing: { outline: p.outline, occlusion: p.occlusion },
      });
    }
  }

  async illustrative() {
    const { structures } = this.plugin.managers.structure.hierarchy.selection;
    await this.plugin.managers.structure.component.applyPreset(
      structures,
      PresetStructureRepresentations.illustrative
    );

    if (this.plugin.canvas3d) {
      this.plugin.canvas3d.setProps({
        postprocessing: {
          outline: {
            name: "on",
            params: {
              scale: 1,
              color: Color(0x000000),
              threshold: 0.25,
              includeTransparent: true,
            },
          },
          occlusion: {
            name: "on",
            params: {
              multiScale: { name: "off", params: {} },
              radius: 5,
              bias: 0.8,
              blurKernelSize: 15,
              samples: 32,
              resolutionScale: 1,
              color: Color(0x000000),
            },
          },
          shadow: { name: "off", params: {} },
        },
      });
    }
  }

  async stylized() {
    this.plugin.managers.structure.component.setOptions({
      ...this.plugin.managers.structure.component.state.options,
      ignoreLight: true,
    });

    if (this.plugin.canvas3d) {
      const pp = this.plugin.canvas3d.props.postprocessing;
      this.plugin.canvas3d.setProps({
        postprocessing: {
          outline: {
            name: "on",
            params:
              pp.outline.name === "on"
                ? pp.outline.params
                : {
                    scale: 1,
                    color: Color(0x000000),
                    threshold: 0.33,
                    includeTransparent: true,
                  },
          },
          occlusion: {
            name: "on",
            params:
              pp.occlusion.name === "on"
                ? pp.occlusion.params
                : {
                    multiScale: { name: "off", params: {} },
                    radius: 5,
                    bias: 0.8,
                    blurKernelSize: 15,
                    samples: 32,
                    resolutionScale: 1,
                    color: Color(0x000000),
                  },
          },
          shadow: { name: "off", params: {} },
        },
      });
    }
  }

  render() {
    return (
      <div className="msp-flex-row">
        <Button
          noOverflow
          title="Applies default representation preset. Set outline and occlusion effects to defaults."
          onClick={this.handleMouseClick}
          onMouseEnter={this.handleMouseEnter}
          onMouseLeave={this.handleMouseLeave}
          style={{ width: "auto" }}
        >
          Default
        </Button>
        <Button
          noOverflow
          title="Applies no representation preset. Enables outline and occlusion effects. Enables ignore-light representation parameter."
          onClick={this.handleMouseClick}
          onMouseEnter={this.handleMouseEnter}
          onMouseLeave={this.handleMouseLeave}
          style={{ width: "auto" }}
        >
          Stylized
        </Button>
        <Button
          noOverflow
          title="Applies illustrative representation preset. Enables outline and occlusion effects. Enables ignore-light parameter."
          onClick={this.handleMouseClick}
          onMouseEnter={this.handleMouseEnter}
          onMouseLeave={this.handleMouseLeave}
          style={{ width: "auto" }}
        >
          Illustrative
        </Button>
      </div>
    );
  }
}
