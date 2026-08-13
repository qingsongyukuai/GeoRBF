//! Dependency-free foundations for Surfe golden-fixture parity tests.
//!
//! GeoRBF does not have a Cargo project yet.  This module therefore uses only
//! the Rust standard library and can be exercised directly with
//! `rustc --test tests/common/parity/mod.rs`.  It parses canonical fixture JSON,
//! validates the v1 envelope, and supplies the frozen numeric comparison
//! profiles.  It never discovers or launches the external C++ oracle.

pub const FIXTURE_SCHEMA: &str = "georbf-surfe-golden";
pub const FIXTURE_SCHEMA_VERSION: u32 = 1;
pub const ORACLE_PROTOCOL: &str = "georbf-surfe-oracle";
pub const ORACLE_PROTOCOL_VERSION: u32 = 1;
pub const SOURCE_REPOSITORY: &str = "https://github.com/MichaelHillier/surfe.git";
pub const SOURCE_COMMIT: &str = "290dbe0ab344f4258a4935f05cad0f153f0f69a4";

#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    Null,
    Bool(bool),
    Number(String),
    String(String),
    Array(Vec<JsonValue>),
    Object(Vec<(String, JsonValue)>),
}

impl JsonValue {
    pub fn canonical_json(&self) -> String {
        let mut output = String::new();
        self.write_canonical(&mut output);
        output
    }

    fn write_canonical(&self, output: &mut String) {
        match self {
            Self::Null => output.push_str("null"),
            Self::Bool(value) => output.push_str(if *value { "true" } else { "false" }),
            Self::Number(value) => output.push_str(value),
            Self::String(value) => write_json_string(value, output),
            Self::Array(values) => {
                output.push('[');
                for (index, value) in values.iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    value.write_canonical(output);
                }
                output.push(']');
            }
            Self::Object(fields) => {
                output.push('{');
                for (index, (name, value)) in fields.iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    write_json_string(name, output);
                    output.push(':');
                    value.write_canonical(output);
                }
                output.push('}');
            }
        }
    }
}

fn write_json_string(value: &str, output: &mut String) {
    output.push('"');
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\u{08}' => output.push_str("\\b"),
            '\u{0c}' => output.push_str("\\f"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            control if control <= '\u{1f}' => {
                output.push_str(&format!("\\u{:04x}", u32::from(control)));
            }
            other => output.push(other),
        }
    }
    output.push('"');
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JsonError {
    pub offset: usize,
    pub reason: &'static str,
}

pub fn parse_json(input: &str) -> Result<JsonValue, JsonError> {
    let mut parser = JsonParser { input, offset: 0 };
    parser.skip_whitespace();
    let value = parser.parse_value()?;
    parser.skip_whitespace();
    if parser.offset != input.len() {
        return Err(parser.error("trailing content"));
    }
    Ok(value)
}

struct JsonParser<'a> {
    input: &'a str,
    offset: usize,
}

impl JsonParser<'_> {
    fn error(&self, reason: &'static str) -> JsonError {
        JsonError {
            offset: self.offset,
            reason,
        }
    }

    fn remaining(&self) -> &str {
        &self.input[self.offset..]
    }

    fn skip_whitespace(&mut self) {
        while matches!(self.peek_byte(), Some(b' ' | b'\n' | b'\r' | b'\t')) {
            self.offset += 1;
        }
    }

    fn peek_byte(&self) -> Option<u8> {
        self.input.as_bytes().get(self.offset).copied()
    }

    fn consume_byte(&mut self, expected: u8) -> Result<(), JsonError> {
        if self.peek_byte() != Some(expected) {
            return Err(self.error("unexpected token"));
        }
        self.offset += 1;
        Ok(())
    }

    fn parse_value(&mut self) -> Result<JsonValue, JsonError> {
        match self.peek_byte() {
            Some(b'n') => self.parse_literal("null", JsonValue::Null),
            Some(b't') => self.parse_literal("true", JsonValue::Bool(true)),
            Some(b'f') => self.parse_literal("false", JsonValue::Bool(false)),
            Some(b'"') => self.parse_string().map(JsonValue::String),
            Some(b'[') => self.parse_array(),
            Some(b'{') => self.parse_object(),
            Some(b'-' | b'0'..=b'9') => self.parse_number().map(JsonValue::Number),
            Some(_) => Err(self.error("invalid JSON value")),
            None => Err(self.error("unexpected end of input")),
        }
    }

    fn parse_literal(&mut self, literal: &str, value: JsonValue) -> Result<JsonValue, JsonError> {
        if !self.remaining().starts_with(literal) {
            return Err(self.error("invalid literal"));
        }
        self.offset += literal.len();
        Ok(value)
    }

    fn parse_string(&mut self) -> Result<String, JsonError> {
        self.consume_byte(b'"')?;
        let mut output = String::new();
        loop {
            let byte = self
                .peek_byte()
                .ok_or_else(|| self.error("unterminated string"))?;
            match byte {
                b'"' => {
                    self.offset += 1;
                    return Ok(output);
                }
                b'\\' => {
                    self.offset += 1;
                    self.parse_escape(&mut output)?;
                }
                0x00..=0x1f => return Err(self.error("unescaped control character")),
                _ => {
                    let character = self
                        .remaining()
                        .chars()
                        .next()
                        .ok_or_else(|| self.error("invalid UTF-8"))?;
                    self.offset += character.len_utf8();
                    output.push(character);
                }
            }
        }
    }

    fn parse_escape(&mut self, output: &mut String) -> Result<(), JsonError> {
        let escape = self
            .peek_byte()
            .ok_or_else(|| self.error("unfinished escape"))?;
        self.offset += 1;
        match escape {
            b'"' => output.push('"'),
            b'\\' => output.push('\\'),
            b'/' => output.push('/'),
            b'b' => output.push('\u{08}'),
            b'f' => output.push('\u{0c}'),
            b'n' => output.push('\n'),
            b'r' => output.push('\r'),
            b't' => output.push('\t'),
            b'u' => {
                let first = self.parse_hex_quad()?;
                let codepoint = if (0xd800..=0xdbff).contains(&first) {
                    if !self.remaining().starts_with("\\u") {
                        return Err(self.error("missing low surrogate"));
                    }
                    self.offset += 2;
                    let second = self.parse_hex_quad()?;
                    if !(0xdc00..=0xdfff).contains(&second) {
                        return Err(self.error("invalid low surrogate"));
                    }
                    0x10000 + ((u32::from(first) - 0xd800) << 10) + (u32::from(second) - 0xdc00)
                } else if (0xdc00..=0xdfff).contains(&first) {
                    return Err(self.error("unpaired low surrogate"));
                } else {
                    u32::from(first)
                };
                output.push(char::from_u32(codepoint).ok_or_else(|| self.error("invalid escape"))?);
            }
            _ => return Err(self.error("unknown escape")),
        }
        Ok(())
    }

    fn parse_hex_quad(&mut self) -> Result<u16, JsonError> {
        let end = self.offset + 4;
        let digits = self
            .input
            .get(self.offset..end)
            .ok_or_else(|| self.error("short unicode escape"))?;
        if !digits.bytes().all(|byte| byte.is_ascii_hexdigit()) {
            return Err(self.error("invalid unicode escape"));
        }
        self.offset = end;
        u16::from_str_radix(digits, 16).map_err(|_| self.error("invalid unicode escape"))
    }

    fn parse_array(&mut self) -> Result<JsonValue, JsonError> {
        self.consume_byte(b'[')?;
        self.skip_whitespace();
        let mut values = Vec::new();
        if self.peek_byte() == Some(b']') {
            self.offset += 1;
            return Ok(JsonValue::Array(values));
        }
        loop {
            values.push(self.parse_value()?);
            self.skip_whitespace();
            match self.peek_byte() {
                Some(b',') => {
                    self.offset += 1;
                    self.skip_whitespace();
                }
                Some(b']') => {
                    self.offset += 1;
                    return Ok(JsonValue::Array(values));
                }
                _ => return Err(self.error("invalid array separator")),
            }
        }
    }

    fn parse_object(&mut self) -> Result<JsonValue, JsonError> {
        self.consume_byte(b'{')?;
        self.skip_whitespace();
        let mut fields = Vec::new();
        if self.peek_byte() == Some(b'}') {
            self.offset += 1;
            return Ok(JsonValue::Object(fields));
        }
        loop {
            if self.peek_byte() != Some(b'"') {
                return Err(self.error("object key must be a string"));
            }
            let name = self.parse_string()?;
            if fields.iter().any(|(existing, _)| existing == &name) {
                return Err(self.error("duplicate object key"));
            }
            self.skip_whitespace();
            self.consume_byte(b':')?;
            self.skip_whitespace();
            let value = self.parse_value()?;
            fields.push((name, value));
            self.skip_whitespace();
            match self.peek_byte() {
                Some(b',') => {
                    self.offset += 1;
                    self.skip_whitespace();
                }
                Some(b'}') => {
                    self.offset += 1;
                    return Ok(JsonValue::Object(fields));
                }
                _ => return Err(self.error("invalid object separator")),
            }
        }
    }

    fn parse_number(&mut self) -> Result<String, JsonError> {
        let start = self.offset;
        if self.peek_byte() == Some(b'-') {
            self.offset += 1;
        }
        match self.peek_byte() {
            Some(b'0') => self.offset += 1,
            Some(b'1'..=b'9') => {
                self.offset += 1;
                while matches!(self.peek_byte(), Some(b'0'..=b'9')) {
                    self.offset += 1;
                }
            }
            _ => return Err(self.error("invalid number integer part")),
        }
        if self.peek_byte() == Some(b'.') {
            self.offset += 1;
            let fraction_start = self.offset;
            while matches!(self.peek_byte(), Some(b'0'..=b'9')) {
                self.offset += 1;
            }
            if self.offset == fraction_start {
                return Err(self.error("empty number fraction"));
            }
        }
        if matches!(self.peek_byte(), Some(b'e' | b'E')) {
            self.offset += 1;
            if matches!(self.peek_byte(), Some(b'+' | b'-')) {
                self.offset += 1;
            }
            let exponent_start = self.offset;
            while matches!(self.peek_byte(), Some(b'0'..=b'9')) {
                self.offset += 1;
            }
            if self.offset == exponent_start {
                return Err(self.error("empty number exponent"));
            }
        }
        let token = self.input[start..self.offset].to_owned();
        let finite = token
            .parse::<f64>()
            .map_err(|_| self.error("number is not binary64"))?;
        if !finite.is_finite() {
            return Err(self.error("non-finite number must use tagged encoding"));
        }
        Ok(token)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixtureViolation {
    TopLevelFields,
    WrongSchema,
    WrongVersion,
    InvalidFixtureId,
    WrongSource,
    InvalidGeneration,
    InvalidDataset,
    InvalidComparison,
    InvalidRequest,
    InvalidExpectedResponse,
}

pub fn validate_fixture(value: &JsonValue) -> Result<(), FixtureViolation> {
    let top = exact_object(
        value,
        &[
            "schema",
            "schema_version",
            "fixture_id",
            "source",
            "generation",
            "dataset",
            "comparison",
            "request",
            "expected",
        ],
    )
    .ok_or(FixtureViolation::TopLevelFields)?;

    if string_field(top, "schema") != Some(FIXTURE_SCHEMA) {
        return Err(FixtureViolation::WrongSchema);
    }
    if number_field(top, "schema_version") != Some("1") {
        return Err(FixtureViolation::WrongVersion);
    }
    let fixture_id = string_field(top, "fixture_id").ok_or(FixtureViolation::InvalidFixtureId)?;
    if fixture_id.is_empty()
        || !fixture_id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'/' | b'_' | b'-' | b'.'))
    {
        return Err(FixtureViolation::InvalidFixtureId);
    }

    validate_source(field(top, "source").ok_or(FixtureViolation::WrongSource)?)?;
    validate_generation(field(top, "generation").ok_or(FixtureViolation::InvalidGeneration)?)?;
    validate_dataset(field(top, "dataset").ok_or(FixtureViolation::InvalidDataset)?)?;
    validate_comparison(field(top, "comparison").ok_or(FixtureViolation::InvalidComparison)?)?;
    let request = field(top, "request").ok_or(FixtureViolation::InvalidRequest)?;
    let expected = field(top, "expected").ok_or(FixtureViolation::InvalidExpectedResponse)?;
    validate_request(request)?;
    validate_expected(request, expected)?;
    Ok(())
}

fn validate_source(value: &JsonValue) -> Result<(), FixtureViolation> {
    let fields = exact_object(
        value,
        &[
            "repository",
            "commit",
            "oracle_protocol",
            "oracle_protocol_version",
        ],
    )
    .ok_or(FixtureViolation::WrongSource)?;
    if string_field(fields, "repository") != Some(SOURCE_REPOSITORY)
        || string_field(fields, "commit") != Some(SOURCE_COMMIT)
        || string_field(fields, "oracle_protocol") != Some(ORACLE_PROTOCOL)
        || number_field(fields, "oracle_protocol_version") != Some("1")
    {
        return Err(FixtureViolation::WrongSource);
    }
    Ok(())
}

fn validate_generation(value: &JsonValue) -> Result<(), FixtureViolation> {
    let fields = exact_object(value, &["command", "environment", "platform", "precision"])
        .ok_or(FixtureViolation::InvalidGeneration)?;
    if !matches!(field(fields, "command"), Some(JsonValue::Array(values)) if !values.is_empty() && values.iter().all(|value| matches!(value, JsonValue::String(_))))
    {
        return Err(FixtureViolation::InvalidGeneration);
    }
    let environment = exact_object(
        field(fields, "environment").ok_or(FixtureViolation::InvalidGeneration)?,
        &["OMP_NUM_THREADS", "LC_ALL", "TZ"],
    )
    .ok_or(FixtureViolation::InvalidGeneration)?;
    if string_field(environment, "OMP_NUM_THREADS") != Some("1")
        || string_field(environment, "LC_ALL") != Some("C")
        || string_field(environment, "TZ") != Some("UTC")
    {
        return Err(FixtureViolation::InvalidGeneration);
    }
    let platform = exact_object(
        field(fields, "platform").ok_or(FixtureViolation::InvalidGeneration)?,
        &["os", "arch", "compiler", "libc", "endianness"],
    )
    .ok_or(FixtureViolation::InvalidGeneration)?;
    if ["os", "arch", "compiler", "libc"]
        .iter()
        .any(|name| string_field(platform, name).is_none_or(str::is_empty))
        || !matches!(string_field(platform, "endianness"), Some("little" | "big"))
    {
        return Err(FixtureViolation::InvalidGeneration);
    }
    let precision = exact_object(
        field(fields, "precision").ok_or(FixtureViolation::InvalidGeneration)?,
        &[
            "default",
            "anisotropy_intermediate",
            "serialization",
            "matrix_order",
        ],
    )
    .ok_or(FixtureViolation::InvalidGeneration)?;
    if string_field(precision, "default") != Some("binary64")
        || string_field(precision, "anisotropy_intermediate") != Some("binary32")
        || string_field(precision, "serialization") != Some("max_digits10")
        || string_field(precision, "matrix_order") != Some("row_major")
    {
        return Err(FixtureViolation::InvalidGeneration);
    }
    Ok(())
}

fn validate_dataset(value: &JsonValue) -> Result<(), FixtureViolation> {
    let fields = exact_object(
        value,
        &[
            "id",
            "description",
            "coverage",
            "request_line_sha256",
            "response_line_sha256",
        ],
    )
    .ok_or(FixtureViolation::InvalidDataset)?;
    if string_field(fields, "id").is_none_or(str::is_empty)
        || string_field(fields, "description").is_none_or(str::is_empty)
        || !matches!(field(fields, "coverage"), Some(JsonValue::Array(values)) if !values.is_empty() && values.iter().all(|value| matches!(value, JsonValue::String(name) if !name.is_empty())))
        || !is_sha256(string_field(fields, "request_line_sha256"))
        || !is_sha256(string_field(fields, "response_line_sha256"))
    {
        return Err(FixtureViolation::InvalidDataset);
    }
    Ok(())
}

fn validate_comparison(value: &JsonValue) -> Result<(), FixtureViolation> {
    let fields = exact_object(
        value,
        &[
            "default_rule",
            "numeric_rules",
            "diagnostic_paths",
            "weight_policy",
            "error_message_policy",
            "acceptance",
        ],
    )
    .ok_or(FixtureViolation::InvalidComparison)?;
    if string_field(fields, "default_rule") != Some("exact")
        || !matches!(
            string_field(fields, "weight_policy"),
            Some("diagnostic_only" | "residual_feasibility_predictions")
        )
        || string_field(fields, "error_message_policy") != Some("diagnostic_only")
    {
        return Err(FixtureViolation::InvalidComparison);
    }
    let numeric_rules = match field(fields, "numeric_rules") {
        Some(JsonValue::Array(values)) => values,
        _ => return Err(FixtureViolation::InvalidComparison),
    };
    let mut patterns = Vec::new();
    for rule in numeric_rules {
        let rule_fields =
            exact_object(rule, &["path", "class"]).ok_or(FixtureViolation::InvalidComparison)?;
        let path = string_field(rule_fields, "path").ok_or(FixtureViolation::InvalidComparison)?;
        let class = string_field(rule_fields, "class")
            .and_then(ToleranceClass::from_name)
            .ok_or(FixtureViolation::InvalidComparison)?;
        if !path.starts_with('/')
            || patterns.contains(&path)
            || class == ToleranceClass::IterationIndex
        {
            return Err(FixtureViolation::InvalidComparison);
        }
        patterns.push(path);
    }
    let diagnostics = match field(fields, "diagnostic_paths") {
        Some(JsonValue::Array(values)) => values,
        _ => return Err(FixtureViolation::InvalidComparison),
    };
    if diagnostics.iter().any(|value| {
        !matches!(value, JsonValue::String(path) if path == "/expected/error/message" || path == "/expected/result/solve/weights/*")
    }) {
        return Err(FixtureViolation::InvalidComparison);
    }
    let acceptance = exact_object(
        field(fields, "acceptance").ok_or(FixtureViolation::InvalidComparison)?,
        &[
            "required_finite",
            "residual_l2_max",
            "relative_residual_max",
            "equality_residual_linf_max",
            "inequality_violation_linf_max",
            "prediction_witnesses_required",
        ],
    )
    .ok_or(FixtureViolation::InvalidComparison)?;
    if !matches!(
        field(acceptance, "required_finite"),
        Some(JsonValue::Bool(_))
    ) || !matches!(
        field(acceptance, "prediction_witnesses_required"),
        Some(JsonValue::Bool(_))
    ) || [
        "residual_l2_max",
        "relative_residual_max",
        "equality_residual_linf_max",
        "inequality_violation_linf_max",
    ]
    .iter()
    .any(|name| !is_nonnegative_number_or_null(field(acceptance, name)))
    {
        return Err(FixtureViolation::InvalidComparison);
    }
    if string_field(fields, "weight_policy") == Some("residual_feasibility_predictions")
        && (field(acceptance, "required_finite") != Some(&JsonValue::Bool(true))
            || field(acceptance, "prediction_witnesses_required") != Some(&JsonValue::Bool(true))
            || [
                "residual_l2_max",
                "relative_residual_max",
                "equality_residual_linf_max",
                "inequality_violation_linf_max",
            ]
            .iter()
            .any(|name| !matches!(field(acceptance, name), Some(JsonValue::Number(_)))))
    {
        return Err(FixtureViolation::InvalidComparison);
    }
    Ok(())
}

fn is_nonnegative_number_or_null(value: Option<&JsonValue>) -> bool {
    match value {
        Some(JsonValue::Null) => true,
        Some(JsonValue::Number(number)) => number
            .parse::<f64>()
            .is_ok_and(|number| number.is_finite() && number >= 0.0),
        _ => false,
    }
}

fn validate_request(value: &JsonValue) -> Result<(), FixtureViolation> {
    let fields = exact_object(
        value,
        &[
            "protocol",
            "protocol_version",
            "request_id",
            "source_commit",
            "operation",
            "input",
            "evidence",
        ],
    )
    .ok_or(FixtureViolation::InvalidRequest)?;
    if string_field(fields, "protocol") != Some(ORACLE_PROTOCOL)
        || number_field(fields, "protocol_version") != Some("1")
        || string_field(fields, "request_id").is_none_or(str::is_empty)
        || string_field(fields, "source_commit") != Some(SOURCE_COMMIT)
        || !matches!(
            string_field(fields, "operation"),
            Some("identity" | "kernel.evaluate" | "model.run" | "solver.run" | "error.probe")
        )
        || !matches!(field(fields, "input"), Some(JsonValue::Object(_)))
        || !matches!(field(fields, "evidence"), Some(JsonValue::Array(_)))
    {
        return Err(FixtureViolation::InvalidRequest);
    }
    Ok(())
}

fn validate_expected(request: &JsonValue, value: &JsonValue) -> Result<(), FixtureViolation> {
    let fields = exact_object(
        value,
        &[
            "protocol",
            "protocol_version",
            "request_id",
            "source",
            "operation",
            "status",
            "result",
            "error",
        ],
    )
    .ok_or(FixtureViolation::InvalidExpectedResponse)?;
    let request_fields = match request {
        JsonValue::Object(fields) => fields,
        _ => return Err(FixtureViolation::InvalidExpectedResponse),
    };
    let source = exact_object(
        field(fields, "source").ok_or(FixtureViolation::InvalidExpectedResponse)?,
        &["repository", "commit"],
    )
    .ok_or(FixtureViolation::InvalidExpectedResponse)?;
    if string_field(fields, "protocol") != Some(ORACLE_PROTOCOL)
        || number_field(fields, "protocol_version") != Some("1")
        || string_field(fields, "request_id") != string_field(request_fields, "request_id")
        || string_field(fields, "operation") != string_field(request_fields, "operation")
        || string_field(source, "repository") != Some(SOURCE_REPOSITORY)
        || string_field(source, "commit") != Some(SOURCE_COMMIT)
    {
        return Err(FixtureViolation::InvalidExpectedResponse);
    }
    match string_field(fields, "status") {
        Some("ok")
            if matches!(field(fields, "result"), Some(JsonValue::Object(_)))
                && field(fields, "error") == Some(&JsonValue::Null) =>
        {
            Ok(())
        }
        Some("error")
            if field(fields, "result") == Some(&JsonValue::Null)
                && matches!(field(fields, "error"), Some(JsonValue::Object(_))) =>
        {
            Ok(())
        }
        _ => Err(FixtureViolation::InvalidExpectedResponse),
    }
}

fn exact_object<'a>(
    value: &'a JsonValue,
    expected_names: &[&str],
) -> Option<&'a [(String, JsonValue)]> {
    let JsonValue::Object(fields) = value else {
        return None;
    };
    if fields.len() != expected_names.len()
        || fields
            .iter()
            .zip(expected_names)
            .any(|((actual, _), expected)| actual != expected)
    {
        return None;
    }
    Some(fields)
}

fn field<'a>(fields: &'a [(String, JsonValue)], name: &str) -> Option<&'a JsonValue> {
    fields
        .iter()
        .find_map(|(candidate, value)| (candidate == name).then_some(value))
}

fn string_field<'a>(fields: &'a [(String, JsonValue)], name: &str) -> Option<&'a str> {
    match field(fields, name) {
        Some(JsonValue::String(value)) => Some(value),
        _ => None,
    }
}

fn number_field<'a>(fields: &'a [(String, JsonValue)], name: &str) -> Option<&'a str> {
    match field(fields, name) {
        Some(JsonValue::Number(value)) => Some(value),
        _ => None,
    }
}

fn is_sha256(value: Option<&str>) -> bool {
    value.is_some_and(|value| {
        value.len() == 64
            && value
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    })
}

pub fn validate_matrix(value: &JsonValue) -> bool {
    let Some(fields) = exact_object(value, &["rows", "cols", "order", "data"]) else {
        return false;
    };
    let (Some(rows), Some(cols), Some("row_major"), Some(JsonValue::Array(data))) = (
        number_field(fields, "rows").and_then(parse_nonnegative_integer),
        number_field(fields, "cols").and_then(parse_nonnegative_integer),
        string_field(fields, "order"),
        field(fields, "data"),
    ) else {
        return false;
    };
    rows.checked_mul(cols) == Some(data.len())
        && data
            .iter()
            .all(|element| EncodedNumber::from_json(element).is_some())
}

pub fn validate_vector(value: &JsonValue) -> bool {
    let Some(fields) = exact_object(value, &["length", "data"]) else {
        return false;
    };
    let (Some(length), Some(JsonValue::Array(data))) = (
        number_field(fields, "length").and_then(parse_nonnegative_integer),
        field(fields, "data"),
    ) else {
        return false;
    };
    length == data.len()
        && data
            .iter()
            .all(|element| EncodedNumber::from_json(element).is_some())
}

fn parse_nonnegative_integer(value: &str) -> Option<usize> {
    if value.starts_with('-') || value.contains(['.', 'e', 'E']) {
        return None;
    }
    value.parse().ok()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerance {
    pub absolute: f64,
    pub relative: f64,
    pub exact_signed_zero: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceClass {
    KernelValue,
    FirstDerivative,
    MixedHessian,
    AnisotropyF32,
    ModifiedKernel,
    MatrixRhs,
    SolverResidual,
    SolverFeasibility,
    PredictionScalar,
    PredictionGradient,
    SolverObjective,
    IterationNumeric,
    IterationIndex,
}

impl ToleranceClass {
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "kernel_value" => Some(Self::KernelValue),
            "first_derivative" => Some(Self::FirstDerivative),
            "mixed_hessian" => Some(Self::MixedHessian),
            "anisotropy_f32" => Some(Self::AnisotropyF32),
            "modified_kernel" => Some(Self::ModifiedKernel),
            "matrix_rhs" => Some(Self::MatrixRhs),
            "solver_residual" => Some(Self::SolverResidual),
            "solver_feasibility" => Some(Self::SolverFeasibility),
            "prediction_scalar" => Some(Self::PredictionScalar),
            "prediction_gradient" => Some(Self::PredictionGradient),
            "solver_objective" => Some(Self::SolverObjective),
            "iteration_numeric" => Some(Self::IterationNumeric),
            "iteration_index" => Some(Self::IterationIndex),
            _ => None,
        }
    }

    pub fn tolerance(self) -> Option<Tolerance> {
        let (absolute, relative) = match self {
            Self::KernelValue => (1.0e-12, 1.0e-11),
            Self::FirstDerivative => (1.0e-11, 1.0e-10),
            Self::MixedHessian => (1.0e-10, 1.0e-9),
            Self::AnisotropyF32 => (2.0e-6, 2.0e-5),
            Self::ModifiedKernel => (2.0e-10, 2.0e-9),
            Self::MatrixRhs => (1.0e-11, 1.0e-10),
            Self::SolverResidual => (1.0e-10, 1.0e-8),
            Self::SolverFeasibility => (1.0e-10, 1.0e-8),
            Self::PredictionScalar => (1.0e-9, 1.0e-8),
            Self::PredictionGradient => (1.0e-8, 1.0e-7),
            Self::SolverObjective => (1.0e-8, 1.0e-7),
            Self::IterationNumeric => (1.0e-9, 1.0e-7),
            Self::IterationIndex => return None,
        };
        Some(Tolerance {
            absolute,
            relative,
            exact_signed_zero: true,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EncodedNumber {
    Finite(f64),
    NaN,
    PositiveInfinity,
    NegativeInfinity,
}

impl EncodedNumber {
    pub fn from_json(value: &JsonValue) -> Option<Self> {
        match value {
            JsonValue::Number(number) => number
                .parse()
                .ok()
                .filter(|value: &f64| value.is_finite())
                .map(Self::Finite),
            JsonValue::Object(fields) if fields.len() == 1 && fields[0].0 == "number_kind" => {
                match &fields[0].1 {
                    JsonValue::String(kind) if kind == "nan" => Some(Self::NaN),
                    JsonValue::String(kind) if kind == "positive_infinity" => {
                        Some(Self::PositiveInfinity)
                    }
                    JsonValue::String(kind) if kind == "negative_infinity" => {
                        Some(Self::NegativeInfinity)
                    }
                    _ => None,
                }
            }
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NumericMismatch {
    Kind,
    SignedZero,
    OutsideTolerance { delta: f64, limit: f64 },
}

pub fn compare_number(
    expected: EncodedNumber,
    actual: EncodedNumber,
    tolerance: Tolerance,
) -> Result<(), NumericMismatch> {
    match (expected, actual) {
        (EncodedNumber::Finite(expected), EncodedNumber::Finite(actual)) => {
            if expected == 0.0
                && actual == 0.0
                && tolerance.exact_signed_zero
                && expected.to_bits() != actual.to_bits()
            {
                return Err(NumericMismatch::SignedZero);
            }
            let delta = (expected - actual).abs();
            let scale = expected.abs().max(actual.abs());
            let limit = tolerance.absolute + tolerance.relative * scale;
            if delta <= limit {
                Ok(())
            } else {
                Err(NumericMismatch::OutsideTolerance { delta, limit })
            }
        }
        (EncodedNumber::NaN, EncodedNumber::NaN)
        | (EncodedNumber::PositiveInfinity, EncodedNumber::PositiveInfinity)
        | (EncodedNumber::NegativeInfinity, EncodedNumber::NegativeInfinity) => Ok(()),
        _ => Err(NumericMismatch::Kind),
    }
}

pub fn pointer_pattern_matches(pattern: &str, path: &str) -> bool {
    if !pattern.starts_with('/') || !path.starts_with('/') {
        return false;
    }
    let pattern_parts: Vec<_> = pattern[1..].split('/').collect();
    let path_parts: Vec<_> = path[1..].split('/').collect();
    pattern_parts.len() == path_parts.len()
        && pattern_parts
            .iter()
            .zip(path_parts)
            .all(|(expected, actual)| *expected == "*" || *expected == actual)
}

pub fn numeric_class_for_path(
    rules: &[(&str, ToleranceClass)],
    path: &str,
) -> Result<Option<ToleranceClass>, &'static str> {
    let mut matches = rules
        .iter()
        .filter(|(pattern, _)| pointer_pattern_matches(pattern, path))
        .map(|(_, class)| *class);
    let first = matches.next();
    if matches.next().is_some() {
        return Err("ambiguous numeric path classification");
    }
    Ok(first)
}

#[cfg(test)]
mod tests {
    use super::*;

    const VALID: &str = include_str!("../../fixtures/schema/valid-minimal.json");
    const INVALID_SOURCE: &str = include_str!("../../fixtures/schema/invalid-wrong-source.json");
    const INVALID_TOLERANCE: &str =
        include_str!("../../fixtures/schema/invalid-ambiguous-tolerance.json");

    #[test]
    fn schema_positive_and_negative_examples_are_distinguished() {
        let valid = parse_json(VALID.trim_end()).expect("positive example must be JSON");
        assert_eq!(validate_fixture(&valid), Ok(()));

        let invalid_source =
            parse_json(INVALID_SOURCE.trim_end()).expect("negative source example must be JSON");
        assert_eq!(
            validate_fixture(&invalid_source),
            Err(FixtureViolation::WrongSource)
        );

        let invalid_tolerance = parse_json(INVALID_TOLERANCE.trim_end())
            .expect("negative tolerance example must be JSON");
        assert_eq!(
            validate_fixture(&invalid_tolerance),
            Err(FixtureViolation::InvalidComparison)
        );

        let missing_solver_ceilings = VALID.replace(
            "\"weight_policy\":\"diagnostic_only\"",
            "\"weight_policy\":\"residual_feasibility_predictions\"",
        );
        let missing_solver_ceilings = parse_json(missing_solver_ceilings.trim_end()).unwrap();
        assert_eq!(
            validate_fixture(&missing_solver_ceilings),
            Err(FixtureViolation::InvalidComparison)
        );
    }

    #[test]
    fn canonical_fixture_serialization_round_trips_byte_for_byte() {
        let parsed = parse_json(VALID.trim_end()).expect("fixture must parse");
        let first = format!("{}\n", parsed.canonical_json());
        let reparsed = parse_json(first.trim_end()).expect("serialized fixture must parse");
        let second = format!("{}\n", reparsed.canonical_json());
        assert_eq!(first, VALID);
        assert_eq!(first, second);
    }

    #[test]
    fn tolerance_boundary_is_inclusive_and_next_value_fails() {
        let tolerance = Tolerance {
            absolute: 0.25,
            relative: 0.0,
            exact_signed_zero: true,
        };
        let boundary = 1.25_f64;
        let after_boundary = f64::from_bits(boundary.to_bits() + 1);
        assert_eq!(
            compare_number(
                EncodedNumber::Finite(1.0),
                EncodedNumber::Finite(boundary),
                tolerance,
            ),
            Ok(())
        );
        assert!(matches!(
            compare_number(
                EncodedNumber::Finite(1.0),
                EncodedNumber::Finite(after_boundary),
                tolerance,
            ),
            Err(NumericMismatch::OutsideTolerance { .. })
        ));
        assert_eq!(
            compare_number(
                EncodedNumber::Finite(0.0),
                EncodedNumber::Finite(-0.0),
                tolerance,
            ),
            Err(NumericMismatch::SignedZero)
        );
    }

    #[test]
    fn non_finite_values_require_exact_tagged_encoding() {
        assert!(parse_json("NaN").is_err());
        assert!(parse_json("Infinity").is_err());
        let nan = parse_json("{\"number_kind\":\"nan\"}").expect("tagged NaN must parse");
        let positive = parse_json("{\"number_kind\":\"positive_infinity\"}")
            .expect("tagged infinity must parse");
        assert_eq!(EncodedNumber::from_json(&nan), Some(EncodedNumber::NaN));
        assert_eq!(
            EncodedNumber::from_json(&positive),
            Some(EncodedNumber::PositiveInfinity)
        );
        assert_eq!(
            compare_number(
                EncodedNumber::NaN,
                EncodedNumber::PositiveInfinity,
                ToleranceClass::KernelValue.tolerance().unwrap(),
            ),
            Err(NumericMismatch::Kind)
        );
    }

    #[test]
    fn numeric_paths_have_one_explicit_class_and_discrete_fields_default_exact() {
        let rules = [
            (
                "/expected/result/first_derivatives/point_1/*",
                ToleranceClass::FirstDerivative,
            ),
            (
                "/expected/result/mixed_hessian/*",
                ToleranceClass::MixedHessian,
            ),
        ];
        assert_eq!(
            numeric_class_for_path(&rules, "/expected/result/first_derivatives/point_1/dx"),
            Ok(Some(ToleranceClass::FirstDerivative))
        );
        assert_eq!(
            numeric_class_for_path(&rules, "/expected/result/layout/rows"),
            Ok(None)
        );
        let ambiguous = [
            ("/expected/result/*/dx", ToleranceClass::FirstDerivative),
            ("/expected/result/point_1/*", ToleranceClass::KernelValue),
        ];
        assert!(numeric_class_for_path(&ambiguous, "/expected/result/point_1/dx").is_err());
    }

    #[test]
    fn tolerance_classes_are_layered_and_iteration_index_is_exact() {
        assert_ne!(
            ToleranceClass::KernelValue.tolerance(),
            ToleranceClass::AnisotropyF32.tolerance()
        );
        assert_ne!(
            ToleranceClass::FirstDerivative.tolerance(),
            ToleranceClass::MixedHessian.tolerance()
        );
        assert_eq!(ToleranceClass::IterationIndex.tolerance(), None);
    }

    #[test]
    fn matrix_and_vector_shapes_and_numeric_encodings_are_strict() {
        let matrix = parse_json(
            "{\"rows\":2,\"cols\":2,\"order\":\"row_major\",\"data\":[1.0,-0.0,{\"number_kind\":\"nan\"},4.0]}",
        )
        .unwrap();
        let vector = parse_json("{\"length\":2,\"data\":[1.0,2.0]}").unwrap();
        assert!(validate_matrix(&matrix));
        assert!(validate_vector(&vector));

        let wrong_shape =
            parse_json("{\"rows\":2,\"cols\":3,\"order\":\"row_major\",\"data\":[1.0,2.0]}")
                .unwrap();
        let wrong_order =
            parse_json("{\"rows\":1,\"cols\":1,\"order\":\"column_major\",\"data\":[1.0]}")
                .unwrap();
        assert!(!validate_matrix(&wrong_shape));
        assert!(!validate_matrix(&wrong_order));
    }
}
