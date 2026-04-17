//! star-chimeric: chimeric alignment detection (Old/Mult/PEmerged).
//!
//! Ports (module → C++ source):
//! - [`segment`]  -> `ChimericSegment.{h,cpp}`
//! - [`align`]    -> `ChimericAlign.{h,cpp}` + `ChimericAlign_chimericBAMoutput.cpp`
//!                    + `ChimericAlign_chimericJunctionOutput.cpp`
//!                    + `ChimericAlign_chimericStitching.cpp`
//! - [`detection`] -> `ChimericDetection.{h,cpp}`
//!                    + `ChimericDetection_chimericDetectionMult.cpp`
//!                    + `ChimericDetection_chimericDetectionOld.cpp`
//!                    + `ChimericDetection_chimericDetectionOldOutput.cpp`
//!                    + `ReadAlign_chimericDetectionOld*.cpp`
//!                    + `ReadAlign_chimericDetectionPEmerged.cpp`
//!
//! M6 populates these modules.

pub mod align;
pub mod detection;
pub mod segment;
