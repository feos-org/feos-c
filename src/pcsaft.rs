use feos_core::{parameter::{IdentifierOption, Parameter}, EquationOfState};
use libc::{c_char, size_t};
use std::{ffi::CStr, slice, sync::Arc};

use crate::feos_equation_of_state_t;
use feos::{pcsaft::{PcSaftParameters, PcSaft}, ideal_gas::IdealGasModel, ResidualModel};

#[allow(non_camel_case_types)]
pub struct feos_pcsaft_parameters_t(Arc<PcSaftParameters>);

#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_pcsaft_parameters_free(ptr: *mut feos_pcsaft_parameters_t) {
    if !ptr.is_null() {
        unsafe {
            drop(Box::from_raw(ptr));
        }
    }
}

/// Create PC-SAFT parameters from json files.
/// 
/// @param substances List of substances.
/// @param n number of substances. Has to match number of items in `substances`
/// @param file_pure path to json file containing PC-SAFT parameter for single substances.
/// @param file_binary path to json file containing binary PC-SAFT parameters or null pointer.
/// @param identifier_option which identifier to match the substances. One of: name, cas, iupacname, formula, inchi, or smiles.
/// 
/// @returns PC-SAFT equation of state without ideal gas model.
#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_pcsaft_parameters_from_json(
    substances: *const *const c_char,
    n: size_t,
    file_pure: *const c_char,
    file_binary: *const c_char,
    identifier_option: *const c_char,
) -> *mut feos_pcsaft_parameters_t {
    let substances_str = unsafe {
        assert!(!substances.is_null());
        slice::from_raw_parts(substances, n)
            .iter()
            .map(|&s| {
                assert!(!s.is_null());
                CStr::from_ptr(s).to_str().unwrap()
            })
            .collect::<Vec<&str>>()
    };
    // let substances_str = substances_cstr.iter().map(|s| s.to_str().unwrap()).collect::<Vec<&str>>();

    let file_pure_str = unsafe {
        assert!(!file_pure.is_null());
        CStr::from_ptr(file_pure).to_str().unwrap()
    };

    let file_binary_str = if file_binary.is_null() {
        None
    } else {
        unsafe { Some(CStr::from_ptr(file_binary).to_str().unwrap()) }
    };

    let identifier_str = unsafe {
        assert!(!identifier_option.is_null());
        CStr::from_ptr(identifier_option).to_str().unwrap()
    };
    let identifier = match identifier_str.to_lowercase().as_str() {
        "name" => IdentifierOption::Name,
        "cas" => IdentifierOption::Cas,
        "iupacname" => IdentifierOption::IupacName,
        "formula" => IdentifierOption::Formula,
        "inchi" => IdentifierOption::Inchi,
        "smiles" => IdentifierOption::Smiles,
        _ => panic!("Wrong identifier."),
    };
    Box::into_raw(Box::new(feos_pcsaft_parameters_t(Arc::new(
        PcSaftParameters::from_json(substances_str, file_pure_str, file_binary_str, identifier)
            .unwrap(),
    ))))
}


#[no_mangle]
#[allow(non_camel_case_types)]
pub extern fn feos_eos_pcsaft(parameters_ptr: *mut feos_pcsaft_parameters_t) -> *mut feos_equation_of_state_t {
    let parameters = unsafe {
        assert!(!parameters_ptr.is_null());
        &*parameters_ptr
    };
    let n = parameters.0.sigma.len();
    let eos = Arc::new(EquationOfState {
        ideal_gas: Arc::new(IdealGasModel::NoModel(n)),
        residual: Arc::new(ResidualModel::PcSaft(PcSaft::new(parameters.0.clone()))),
    });
    Box::into_raw(Box::new(feos_equation_of_state_t(eos)))
}