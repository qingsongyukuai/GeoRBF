//! Dependency-free checks for the external Surfe oracle envelope.
//!
//! This helper deliberately does not launch or discover the oracle. Normal Rust
//! builds and tests must remain independent of the frozen C++ reference.

pub const PROTOCOL: &str = "georbf-surfe-oracle";
pub const PROTOCOL_VERSION: u32 = 1;
pub const SOURCE_REPOSITORY: &str = "https://github.com/MichaelHillier/surfe.git";
pub const SOURCE_COMMIT: &str = "290dbe0ab344f4258a4935f05cad0f153f0f69a4";

const OPERATIONS: [&str; 5] = [
    "identity",
    "kernel.evaluate",
    "model.run",
    "solver.run",
    "error.probe",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProtocolViolation {
    NotOneJsonLine,
    InvalidObjectShape,
    WrongProtocol,
    WrongVersion,
    WrongSource,
    MissingRequestId,
    UnsupportedOperation,
    InvalidSuccessPayload,
    InvalidErrorPayload,
}

fn is_one_object_line(line: &str) -> bool {
    if line.is_empty() || line.contains('\n') || line.contains('\r') {
        return false;
    }

    let bytes = line.as_bytes();
    if bytes.first() != Some(&b'{') || bytes.last() != Some(&b'}') {
        return false;
    }

    let mut in_string = false;
    let mut escaped = false;
    let mut object_depth = 0_i32;
    let mut array_depth = 0_i32;
    for byte in bytes {
        if in_string {
            if escaped {
                escaped = false;
            } else if *byte == b'\\' {
                escaped = true;
            } else if *byte == b'"' {
                in_string = false;
            }
            continue;
        }

        match *byte {
            b'"' => in_string = true,
            b'{' => object_depth += 1,
            b'}' => {
                object_depth -= 1;
                if object_depth < 0 {
                    return false;
                }
            }
            b'[' => array_depth += 1,
            b']' => {
                array_depth -= 1;
                if array_depth < 0 {
                    return false;
                }
            }
            _ => {}
        }
    }

    !in_string && !escaped && object_depth == 0 && array_depth == 0
}

fn has_operation(line: &str) -> bool {
    OPERATIONS
        .iter()
        .any(|operation| line.contains(&format!("\"operation\":\"{operation}\"")))
}

pub fn validate_request_line(line: &str) -> Result<(), ProtocolViolation> {
    if !is_one_object_line(line) {
        return Err(ProtocolViolation::NotOneJsonLine);
    }
    if !line.starts_with(&format!("{{\"protocol\":\"{PROTOCOL}\",")) {
        return Err(ProtocolViolation::WrongProtocol);
    }
    if !line.contains(&format!("\"protocol_version\":{PROTOCOL_VERSION},")) {
        return Err(ProtocolViolation::WrongVersion);
    }
    if !line.contains(&format!("\"source_commit\":\"{SOURCE_COMMIT}\"")) {
        return Err(ProtocolViolation::WrongSource);
    }
    if line.contains("\"request_id\":\"\"") || !line.contains("\"request_id\":\"") {
        return Err(ProtocolViolation::MissingRequestId);
    }
    if !has_operation(line) {
        return Err(ProtocolViolation::UnsupportedOperation);
    }
    if !line.contains("\"input\":{") || !line.contains("\"evidence\":[") {
        return Err(ProtocolViolation::InvalidObjectShape);
    }
    Ok(())
}

pub fn validate_response_line(line: &str) -> Result<(), ProtocolViolation> {
    if !is_one_object_line(line) {
        return Err(ProtocolViolation::NotOneJsonLine);
    }
    if !line.starts_with(&format!("{{\"protocol\":\"{PROTOCOL}\",")) {
        return Err(ProtocolViolation::WrongProtocol);
    }
    if !line.contains(&format!("\"protocol_version\":{PROTOCOL_VERSION},")) {
        return Err(ProtocolViolation::WrongVersion);
    }
    let source = format!(
        "\"source\":{{\"repository\":\"{SOURCE_REPOSITORY}\",\"commit\":\"{SOURCE_COMMIT}\"}}"
    );
    if !line.contains(&source) {
        return Err(ProtocolViolation::WrongSource);
    }
    if line.contains("\"request_id\":\"\"") || !line.contains("\"request_id\":\"") {
        return Err(ProtocolViolation::MissingRequestId);
    }
    if !has_operation(line) {
        return Err(ProtocolViolation::UnsupportedOperation);
    }

    if line.contains("\"status\":\"ok\"") {
        if line.contains("\"result\":null") || !line.ends_with(",\"error\":null}") {
            return Err(ProtocolViolation::InvalidSuccessPayload);
        }
    } else if line.contains("\"status\":\"error\"") {
        if !line.contains("\"result\":null,\"error\":{") || line.ends_with(",\"error\":null}") {
            return Err(ProtocolViolation::InvalidErrorPayload);
        }
    } else {
        return Err(ProtocolViolation::InvalidObjectShape);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const REQUEST: &str = concat!(
        "{\"protocol\":\"georbf-surfe-oracle\",\"protocol_version\":1,",
        "\"request_id\":\"unit-1\",",
        "\"source_commit\":\"290dbe0ab344f4258a4935f05cad0f153f0f69a4\",",
        "\"operation\":\"identity\",\"input\":{},\"evidence\":[]}"
    );
    const OK_RESPONSE: &str = concat!(
        "{\"protocol\":\"georbf-surfe-oracle\",\"protocol_version\":1,",
        "\"request_id\":\"unit-1\",\"source\":{",
        "\"repository\":\"https://github.com/MichaelHillier/surfe.git\",",
        "\"commit\":\"290dbe0ab344f4258a4935f05cad0f153f0f69a4\"},",
        "\"operation\":\"identity\",\"status\":\"ok\",",
        "\"result\":{},\"error\":null}"
    );
    const ERROR_RESPONSE: &str = concat!(
        "{\"protocol\":\"georbf-surfe-oracle\",\"protocol_version\":1,",
        "\"request_id\":\"unit-2\",\"source\":{",
        "\"repository\":\"https://github.com/MichaelHillier/surfe.git\",",
        "\"commit\":\"290dbe0ab344f4258a4935f05cad0f153f0f69a4\"},",
        "\"operation\":\"error.probe\",\"status\":\"error\",",
        "\"result\":null,\"error\":{\"stage\":\"request\",",
        "\"category\":\"invalid_request\"}}"
    );

    #[test]
    fn accepts_canonical_request_and_both_response_variants() {
        assert_eq!(validate_request_line(REQUEST), Ok(()));
        assert_eq!(validate_response_line(OK_RESPONSE), Ok(()));
        assert_eq!(validate_response_line(ERROR_RESPONSE), Ok(()));
    }

    #[test]
    fn rejects_source_drift_and_multiline_output() {
        let drifted = OK_RESPONSE.replace(SOURCE_COMMIT, "main");
        assert_eq!(
            validate_response_line(&drifted),
            Err(ProtocolViolation::WrongSource)
        );
        assert_eq!(
            validate_response_line(&format!("{OK_RESPONSE}\n")),
            Err(ProtocolViolation::NotOneJsonLine)
        );
    }

    #[test]
    fn rejects_inconsistent_payload_status() {
        let invalid = OK_RESPONSE.replace("\"result\":{}", "\"result\":null");
        assert_eq!(
            validate_response_line(&invalid),
            Err(ProtocolViolation::InvalidSuccessPayload)
        );
    }

    #[test]
    fn validates_external_smoke_when_explicitly_supplied() {
        if let Ok(request) = std::env::var("SURFE_ORACLE_SMOKE_REQUEST") {
            validate_request_line(&request).expect("external oracle request must match protocol");
        }
        if let Ok(response) = std::env::var("SURFE_ORACLE_SMOKE_RESPONSE") {
            validate_response_line(&response)
                .expect("external oracle response must match protocol");
        }
    }
}
