use feos::pcsaft::{PcSaft, PcSaftBinaryRecord, PcSaftParameters, PcSaftRecord};
use feos::{ideal_gas::IdealGasModel, ResidualModel};
use feos_core::parameter::{BinaryRecord, Identifier, IdentifierOption, Parameter, PureRecord};
use feos_core::si::{BAR, KELVIN, KILOGRAM, METER, MOL};
use feos_core::{Contributions, DensityInitialization, EquationOfState, SolverOptions, State};
use libc::{c_char, size_t};
use ndarray::Array1;
use serde_json::Value;
use std::ffi::CStr;
use std::slice;
use std::sync::Arc;

pub mod pcsaft;

#[allow(non_camel_case_types)]
pub struct feos_equation_of_state_t(Arc<EquationOfState<IdealGasModel, ResidualModel>>);

#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_eos_from_json(json: *const c_char) -> *mut feos_equation_of_state_t {
    let c_str = unsafe {
        assert!(!json.is_null());
        CStr::from_ptr(json)
    };
    let r_str = c_str.to_str().unwrap();

    let v: Value = serde_json::from_str(r_str).unwrap();
    let residual_model_name = match &v["residual_model"] {
        Value::Null => return std::ptr::null_mut(),
        Value::String(s) => s,
        _ => return std::ptr::null_mut(),
    };

    match residual_model_name.to_lowercase().as_str() {
        "pc-saft" | "pcsaft" => {
            let pure_records: Vec<PureRecord<PcSaftRecord>> =
                serde_json::from_value(v["residual_substance_parameters"].clone()).unwrap();    
            let binary_records: Vec<BinaryRecord<Identifier, PcSaftBinaryRecord>> =
                serde_json::from_value(v["residual_binary_parameters"].clone()).unwrap_or(Vec::new());
            let binary_matrix = <PcSaftParameters as Parameter>::binary_matrix_from_records(
                &pure_records,
                &binary_records,
                IdentifierOption::Name,
            );
            let parameters =
                Arc::new(PcSaftParameters::from_records(pure_records, binary_matrix).unwrap());
            let n = parameters.sigma.len();
            let eos = Arc::new(EquationOfState {
                ideal_gas: Arc::new(IdealGasModel::NoModel(n)),
                residual: Arc::new(ResidualModel::PcSaft(PcSaft::new(parameters))),
            });
            Box::into_raw(Box::new(feos_equation_of_state_t(eos)))
        }
        _ => std::ptr::null_mut(),
    }
}

#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_eos_free(ptr: *mut feos_equation_of_state_t) {
    if !ptr.is_null() {
        unsafe {
            drop(Box::from_raw(ptr));
        }
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone)]
pub struct feos_state_t(State<EquationOfState<IdealGasModel, ResidualModel>>);

#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_state_free(ptr: *mut feos_state_t) {
    if !ptr.is_null() {
        unsafe {
            drop(Box::from_raw(ptr));
        }
    }
}

#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_state_new_npt(
    eos_ptr: *const feos_equation_of_state_t,
    t_k: f64,
    p_bar: f64,
    moles: *const f64,
    len: size_t,
    phase: *const c_char,
) -> *mut feos_state_t {
    let eos = unsafe {
        assert!(!eos_ptr.is_null());
        &*eos_ptr
    };
    let phase_str = unsafe {
        assert!(!phase.is_null());
        CStr::from_ptr(phase)
    };
    let n = unsafe {
        assert!(!moles.is_null());
        Array1::from_iter(slice::from_raw_parts(moles, len).iter().cloned()) * MOL
    };

    let density_initialization = match phase_str.to_str().unwrap() {
        "liquid" => DensityInitialization::Liquid,
        "vapor" => DensityInitialization::Vapor,
        _ => DensityInitialization::None,
    };
    let state = State::new_npt(
        &eos.0,
        t_k * KELVIN,
        p_bar * BAR,
        &n,
        density_initialization,
    )
    .unwrap();
    Box::into_raw(Box::new(feos_state_t(state)))
}

#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_state_pressure(
    state_ptr: *const feos_state_t,
    contributions: size_t,
) -> f64 {
    let state = unsafe {
        assert!(!state_ptr.is_null());
        &*state_ptr
    };
    let c = match contributions {
        0 => Contributions::IdealGas,
        1 => Contributions::Residual,
        _ => Contributions::Total,
    };
    state.0.pressure(c).convert_into(BAR)
}

#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_state_density(state_ptr: *const feos_state_t) -> f64 {
    let state = unsafe {
        assert!(!state_ptr.is_null());
        &*state_ptr
    };
    state.0.density.convert_into(MOL / (METER * METER * METER))
}

/// Mass density of the state.
/// 
/// @params state_ptr Pointer to the state.
/// @returns mass density in units of kg/m3
#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_state_mass_density(state_ptr: *const feos_state_t) -> f64 {
    let state = unsafe {
        assert!(!state_ptr.is_null());
        &*state_ptr
    };
    state
        .0
        .mass_density()
        .convert_into(KILOGRAM / (METER * METER * METER))
}

/// Test if the current state is stable.
///
/// @param state_ptr `*feos_state_t` to test stability for
///
/// @returns True if the state is stable, false otherwise.
#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_state_is_stable(state_ptr: *const feos_state_t) -> bool {
    let state = unsafe {
        assert!(!state_ptr.is_null());
        &*state_ptr
    };
    state.0.is_stable(SolverOptions::default()).unwrap()
}
