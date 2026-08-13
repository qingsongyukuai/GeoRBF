#![allow(dead_code)]

use std::{
    hint::black_box,
    sync::{
        mpsc::{self, Receiver, Sender},
        Arc,
    },
    thread,
    thread::JoinHandle,
    time::{Duration, Instant},
};

use georbf::{
    assemble_system_with_layout, constraint_layout, fit_single_surface_linear,
    solve_dense_partial_pivot_lu, ConstraintError, Constraints, DenseMatrix, FunctionalKernel,
    Interface, IsotropicKernel, ModelType, Parameters, Planar, Point, RbfKernel,
    SingleSurfaceLinearError, SingleSurfaceLinearModel, Tangent,
};

pub const BENCHMARK_FORMAT: &str = "georbf-performance-v1";
pub const BENCHMARK_CASE: &str = "single_surface_cubic_dense_v1";
pub const FIXED_MULTI_THREADS: usize = 2;
pub const DEFAULT_SAMPLES: usize = 9;
pub const DEFAULT_WARMUPS: usize = 5;

#[derive(Debug)]
pub struct BenchmarkCase {
    pub constraints: Constraints,
    pub queries: Vec<Point>,
    pub parameters: Parameters,
    pub checksum: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Stage {
    Preprocess,
    Assembly,
    Solve,
    ScalarEvaluation,
    GradientEvaluation,
    EndToEnd,
}

impl Stage {
    pub const ALL: [Self; 6] = [
        Self::Preprocess,
        Self::Assembly,
        Self::Solve,
        Self::ScalarEvaluation,
        Self::GradientEvaluation,
        Self::EndToEnd,
    ];

    pub const fn name(self) -> &'static str {
        match self {
            Self::Preprocess => "preprocess",
            Self::Assembly => "assembly",
            Self::Solve => "solve",
            Self::ScalarEvaluation => "scalar_evaluation",
            Self::GradientEvaluation => "gradient_evaluation",
            Self::EndToEnd => "end_to_end",
        }
    }

    pub fn parse(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|stage| stage.name() == name)
    }

    const fn iterations(self) -> usize {
        match self {
            Self::Preprocess => 64,
            Self::Assembly => 4,
            Self::Solve => 64,
            Self::ScalarEvaluation => 32,
            Self::GradientEvaluation | Self::EndToEnd => 2,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BenchmarkOptions {
    pub samples: usize,
    pub warmups: usize,
    pub stage: Option<Stage>,
    pub threads: Option<usize>,
}

impl Default for BenchmarkOptions {
    fn default() -> Self {
        Self {
            samples: DEFAULT_SAMPLES,
            warmups: DEFAULT_WARMUPS,
            stage: None,
            threads: None,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SampleRow {
    pub stage: Stage,
    pub threads: usize,
    pub sample: usize,
    pub nanoseconds: u128,
    pub checksum: u64,
}

#[derive(Debug)]
pub enum BenchmarkError {
    Constraint(ConstraintError),
    Fit(SingleSurfaceLinearError),
    InvalidCase(&'static str),
}

impl From<ConstraintError> for BenchmarkError {
    fn from(error: ConstraintError) -> Self {
        Self::Constraint(error)
    }
}

impl From<SingleSurfaceLinearError> for BenchmarkError {
    fn from(error: SingleSurfaceLinearError) -> Self {
        Self::Fit(error)
    }
}

pub fn benchmark_case() -> Result<BenchmarkCase, ConstraintError> {
    let mut constraints = Constraints::default();
    for index in 0..96 {
        let x = ((index * 37) % 101) as f64 / 8.0 + (index % 3) as f64 / 1_024.0;
        let y = ((index * 53 + 7) % 103) as f64 / 8.0 + (index % 5) as f64 / 512.0;
        let z = ((index * 71 + 11) % 107) as f64 / 8.0 + (index % 7) as f64 / 256.0;
        let level = x / 8.0 - y / 16.0 + z / 4.0 + ((index * 19) % 17) as f64 / 16_384.0;
        constraints.interfaces.push(Interface::new(x, y, z, level)?);
    }

    for index in 0..16 {
        let x = 0.5 + ((index * 29) % 47) as f64 / 4.0;
        let y = 0.75 + ((index * 31 + 5) % 43) as f64 / 4.0;
        let z = 0.625 + ((index * 17 + 9) % 41) as f64 / 4.0;
        let nx = 0.3125 + (index % 3) as f64 / 128.0;
        let ny = -0.25 + (index % 5) as f64 / 128.0;
        let nz = (1.0 - nx * nx - ny * ny).sqrt();
        constraints
            .planars
            .push(Planar::from_normal(x, y, z, nx, ny, nz)?);
    }

    for index in 0..8 {
        let x = 0.25 + ((index * 19) % 37) as f64 / 4.0;
        let y = 0.5 + ((index * 23 + 3) % 31) as f64 / 4.0;
        let z = 0.75 + ((index * 13 + 7) % 29) as f64 / 4.0;
        let tx = 0.375 + index as f64 / 256.0;
        let ty = 0.625 - index as f64 / 256.0;
        let tz = -0.1875 + index as f64 / 512.0;
        constraints
            .tangents
            .push(Tangent::new(x, y, z, tx, ty, tz)?);
    }

    let mut queries = Vec::with_capacity(4_096);
    for index in 0..4_096 {
        let x = ((index * 43 + 3) % 257) as f64 / 16.0;
        let y = ((index * 89 + 17) % 263) as f64 / 16.0;
        let z = ((index * 131 + 29) % 269) as f64 / 16.0;
        queries.push(Point::new(x, y, z)?);
    }

    let mut parameters = Parameters {
        model_type: ModelType::SingleSurface,
        ..Parameters::default()
    };
    parameters.basis_type = RbfKernel::Cubic;
    parameters.polynomial_order = 1;
    parameters.shape_parameter = 1.0;

    let checksum = dataset_checksum(&constraints, &queries);
    Ok(BenchmarkCase {
        constraints,
        queries,
        parameters,
        checksum,
    })
}

pub fn evaluate_scalars_with_threads(
    model: &SingleSurfaceLinearModel,
    points: &[Point],
    threads: usize,
) -> Result<Vec<f64>, SingleSurfaceLinearError> {
    if threads <= 1 || points.len() < threads {
        return model.evaluate_scalars(points);
    }
    let chunk_size = points.len().div_ceil(threads);
    thread::scope(|scope| {
        let mut handles = Vec::new();
        for chunk in points.chunks(chunk_size) {
            handles.push(scope.spawn(move || model.evaluate_scalars(chunk)));
        }
        let mut values = Vec::with_capacity(points.len());
        for handle in handles {
            values.extend(handle.join().expect("benchmark worker must not panic")?);
        }
        Ok(values)
    })
}

pub fn evaluate_gradients_with_threads(
    model: &SingleSurfaceLinearModel,
    points: &[Point],
    threads: usize,
) -> Result<Vec<[f64; 3]>, SingleSurfaceLinearError> {
    if threads <= 1 || points.len() < threads {
        return model.evaluate_gradients(points);
    }
    let chunk_size = points.len().div_ceil(threads);
    thread::scope(|scope| {
        let mut handles = Vec::new();
        for chunk in points.chunks(chunk_size) {
            handles.push(scope.spawn(move || model.evaluate_gradients(chunk)));
        }
        let mut values = Vec::with_capacity(points.len());
        for handle in handles {
            values.extend(handle.join().expect("benchmark worker must not panic")?);
        }
        Ok(values)
    })
}

pub fn result_checksum(scalars: &[f64], gradients: &[[f64; 3]]) -> u64 {
    let mut hash = FNV_OFFSET;
    mix_u64(&mut hash, scalars.len() as u64);
    for value in scalars {
        mix_u64(&mut hash, value.to_bits());
    }
    mix_u64(&mut hash, gradients.len() as u64);
    for gradient in gradients {
        for value in gradient {
            mix_u64(&mut hash, value.to_bits());
        }
    }
    hash
}

enum EvaluationRequest {
    Scalars {
        model: Arc<SingleSurfaceLinearModel>,
        points: Arc<[Point]>,
        start: usize,
        end: usize,
        repetitions: usize,
    },
    Gradients {
        model: Arc<SingleSurfaceLinearModel>,
        points: Arc<[Point]>,
        start: usize,
        end: usize,
        repetitions: usize,
    },
    Shutdown,
}

enum EvaluationResponse {
    Scalars(usize, Result<Vec<Vec<f64>>, SingleSurfaceLinearError>),
    Gradients(usize, Result<Vec<Vec<[f64; 3]>>, SingleSurfaceLinearError>),
}

struct EvaluationPool {
    senders: Vec<Sender<EvaluationRequest>>,
    responses: Receiver<EvaluationResponse>,
    workers: Vec<JoinHandle<()>>,
}

impl EvaluationPool {
    fn new(worker_count: usize) -> Self {
        let (response_sender, responses) = mpsc::channel();
        let mut senders = Vec::with_capacity(worker_count);
        let mut workers = Vec::with_capacity(worker_count);
        for worker in 0..worker_count {
            let (sender, requests) = mpsc::channel();
            let response_sender = response_sender.clone();
            senders.push(sender);
            workers.push(thread::spawn(move || {
                while let Ok(request) = requests.recv() {
                    let response = match request {
                        EvaluationRequest::Scalars {
                            model,
                            points,
                            start,
                            end,
                            repetitions,
                        } => EvaluationResponse::Scalars(
                            worker,
                            (0..repetitions)
                                .map(|_| model.evaluate_scalars(&points[start..end]))
                                .collect(),
                        ),
                        EvaluationRequest::Gradients {
                            model,
                            points,
                            start,
                            end,
                            repetitions,
                        } => EvaluationResponse::Gradients(
                            worker,
                            (0..repetitions)
                                .map(|_| model.evaluate_gradients(&points[start..end]))
                                .collect(),
                        ),
                        EvaluationRequest::Shutdown => break,
                    };
                    if response_sender.send(response).is_err() {
                        break;
                    }
                }
            }));
        }
        Self {
            senders,
            responses,
            workers,
        }
    }

    fn evaluate_scalars(
        &self,
        model: Arc<SingleSurfaceLinearModel>,
        points: Arc<[Point]>,
    ) -> Result<Vec<f64>, BenchmarkError> {
        Ok(self
            .evaluate_scalars_repeated(model, points, 1)?
            .pop()
            .expect("one scalar repetition was requested"))
    }

    fn evaluate_scalars_repeated(
        &self,
        model: Arc<SingleSurfaceLinearModel>,
        points: Arc<[Point]>,
        repetitions: usize,
    ) -> Result<Vec<Vec<f64>>, BenchmarkError> {
        for (worker, sender) in self.senders.iter().enumerate() {
            let (start, end) = worker_range(points.len(), worker, self.senders.len());
            sender
                .send(EvaluationRequest::Scalars {
                    model: Arc::clone(&model),
                    points: Arc::clone(&points),
                    start,
                    end,
                    repetitions,
                })
                .map_err(|_| BenchmarkError::InvalidCase("scalar worker stopped"))?;
        }
        let mut chunks = (0..self.senders.len())
            .map(|_| None)
            .collect::<Vec<Option<Vec<Vec<f64>>>>>();
        for _ in 0..self.senders.len() {
            match self.responses.recv() {
                Ok(EvaluationResponse::Scalars(worker, result)) => {
                    chunks[worker] = Some(result?);
                }
                Ok(EvaluationResponse::Gradients(_, _)) => {
                    return Err(BenchmarkError::InvalidCase("unexpected gradient response"));
                }
                Err(_) => return Err(BenchmarkError::InvalidCase("scalar worker stopped")),
            }
        }
        let chunks = chunks.into_iter().flatten().collect::<Vec<_>>();
        let mut values = (0..repetitions)
            .map(|_| Vec::with_capacity(points.len()))
            .collect::<Vec<_>>();
        for worker in chunks {
            for (repetition, chunk) in worker.into_iter().enumerate() {
                values[repetition].extend(chunk);
            }
        }
        Ok(values)
    }

    fn evaluate_gradients(
        &self,
        model: Arc<SingleSurfaceLinearModel>,
        points: Arc<[Point]>,
    ) -> Result<Vec<[f64; 3]>, BenchmarkError> {
        Ok(self
            .evaluate_gradients_repeated(model, points, 1)?
            .pop()
            .expect("one gradient repetition was requested"))
    }

    fn evaluate_gradients_repeated(
        &self,
        model: Arc<SingleSurfaceLinearModel>,
        points: Arc<[Point]>,
        repetitions: usize,
    ) -> Result<Vec<Vec<[f64; 3]>>, BenchmarkError> {
        for (worker, sender) in self.senders.iter().enumerate() {
            let (start, end) = worker_range(points.len(), worker, self.senders.len());
            sender
                .send(EvaluationRequest::Gradients {
                    model: Arc::clone(&model),
                    points: Arc::clone(&points),
                    start,
                    end,
                    repetitions,
                })
                .map_err(|_| BenchmarkError::InvalidCase("gradient worker stopped"))?;
        }
        let mut chunks = (0..self.senders.len())
            .map(|_| None)
            .collect::<Vec<Option<Vec<Vec<[f64; 3]>>>>>();
        for _ in 0..self.senders.len() {
            match self.responses.recv() {
                Ok(EvaluationResponse::Gradients(worker, result)) => {
                    chunks[worker] = Some(result?);
                }
                Ok(EvaluationResponse::Scalars(_, _)) => {
                    return Err(BenchmarkError::InvalidCase("unexpected scalar response"));
                }
                Err(_) => return Err(BenchmarkError::InvalidCase("gradient worker stopped")),
            }
        }
        let chunks = chunks.into_iter().flatten().collect::<Vec<_>>();
        let mut values = (0..repetitions)
            .map(|_| Vec::with_capacity(points.len()))
            .collect::<Vec<_>>();
        for worker in chunks {
            for (repetition, chunk) in worker.into_iter().enumerate() {
                values[repetition].extend(chunk);
            }
        }
        Ok(values)
    }
}

impl Drop for EvaluationPool {
    fn drop(&mut self) {
        for sender in &self.senders {
            let _ = sender.send(EvaluationRequest::Shutdown);
        }
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

fn worker_range(length: usize, worker: usize, worker_count: usize) -> (usize, usize) {
    (
        length * worker / worker_count,
        length * (worker + 1) / worker_count,
    )
}

pub fn run_benchmark(options: BenchmarkOptions) -> Result<Vec<SampleRow>, BenchmarkError> {
    if options.samples == 0 {
        return Err(BenchmarkError::InvalidCase("samples must be non-zero"));
    }
    if options
        .threads
        .is_some_and(|threads| threads != 1 && threads != FIXED_MULTI_THREADS)
    {
        return Err(BenchmarkError::InvalidCase(
            "threads must be one or the fixed multi-thread count",
        ));
    }
    let case = benchmark_case()?;
    let prepared = prepare(&case)?;
    let assembled = assemble_system_with_layout(
        prepared.layout.clone(),
        &prepared.constraints,
        &case.parameters,
        FunctionalKernel::from(&prepared.kernel),
    )
    .map_err(|_| BenchmarkError::InvalidCase("fixed assembly failed"))?;
    let rhs = assembled
        .constraints()
        .linear_rhs()
        .ok_or(BenchmarkError::InvalidCase("fixed branch is not linear"))?;
    let model = Arc::new(fit_single_surface_linear(
        &case.constraints,
        &case.parameters,
    )?);
    let query_points: Arc<[Point]> = case.queries.clone().into();
    let evaluation_pool = EvaluationPool::new(FIXED_MULTI_THREADS);

    let mut rows = Vec::with_capacity(Stage::ALL.len() * 2 * options.samples);
    for threads in [1, FIXED_MULTI_THREADS] {
        if options.threads.is_some_and(|selected| selected != threads) {
            continue;
        }
        for stage in Stage::ALL {
            if options.stage.is_some_and(|selected| selected != stage) {
                continue;
            }
            let iterations = stage.iterations();
            for _ in 0..options.warmups {
                black_box(run_stage(
                    stage,
                    threads,
                    &case,
                    &prepared,
                    &assembled,
                    rhs,
                    Arc::clone(&model),
                    Arc::clone(&query_points),
                    &evaluation_pool,
                    iterations,
                )?);
            }
            for sample in 0..options.samples {
                let started = Instant::now();
                let checksum = run_stage(
                    stage,
                    threads,
                    &case,
                    &prepared,
                    &assembled,
                    rhs,
                    Arc::clone(&model),
                    Arc::clone(&query_points),
                    &evaluation_pool,
                    iterations,
                )?;
                let elapsed = started.elapsed();
                rows.push(SampleRow {
                    stage,
                    threads,
                    sample,
                    nanoseconds: per_iteration(elapsed, iterations),
                    checksum: black_box(checksum),
                });
            }
        }
    }
    Ok(rows)
}

struct PreparedCase {
    constraints: Constraints,
    kernel: IsotropicKernel,
    layout: georbf::ConstraintLayout,
    checksum: u64,
}

fn prepare(case: &BenchmarkCase) -> Result<PreparedCase, BenchmarkError> {
    let mut constraints = case.constraints.clone();
    constraints.remove_collocated();
    constraints
        .interface_grouping()
        .ok_or(BenchmarkError::InvalidCase("fixed grouping failed"))?;
    let layout = constraint_layout(ModelType::SingleSurface, &constraints, &case.parameters)
        .map_err(|_| BenchmarkError::InvalidCase("fixed layout failed"))?;
    let kernel = IsotropicKernel::new(RbfKernel::Cubic, 1.0);
    let mut checksum = case.checksum;
    mix_u64(&mut checksum, layout.matrix_size() as u64);
    Ok(PreparedCase {
        constraints,
        kernel,
        layout,
        checksum,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_stage(
    stage: Stage,
    threads: usize,
    case: &BenchmarkCase,
    prepared: &PreparedCase,
    assembled: &georbf::AssembledSystem,
    rhs: &georbf::DenseVector,
    model: Arc<SingleSurfaceLinearModel>,
    query_points: Arc<[Point]>,
    evaluation_pool: &EvaluationPool,
    iterations: usize,
) -> Result<u64, BenchmarkError> {
    if threads > 1 && stage == Stage::ScalarEvaluation {
        let values = evaluation_pool.evaluate_scalars_repeated(
            Arc::clone(&model),
            Arc::clone(&query_points),
            iterations,
        )?;
        let mut checksum = FNV_OFFSET;
        for scalars in values {
            mix_u64(&mut checksum, vector_checksum(&scalars));
        }
        return Ok(checksum);
    }
    if threads > 1 && stage == Stage::GradientEvaluation {
        let values = evaluation_pool.evaluate_gradients_repeated(
            Arc::clone(&model),
            Arc::clone(&query_points),
            iterations,
        )?;
        let mut checksum = FNV_OFFSET;
        for gradients in values {
            mix_u64(&mut checksum, gradient_checksum(&gradients));
        }
        return Ok(checksum);
    }

    let mut checksum = FNV_OFFSET;
    for _ in 0..iterations {
        let value = match stage {
            Stage::Preprocess => prepare(case)?.checksum,
            Stage::Assembly => {
                let value = assemble_system_with_layout(
                    prepared.layout.clone(),
                    &prepared.constraints,
                    &case.parameters,
                    FunctionalKernel::from(&prepared.kernel),
                )
                .map_err(|_| BenchmarkError::InvalidCase("timed assembly failed"))?;
                matrix_checksum(value.interpolation_matrix())
            }
            Stage::Solve => {
                let value = solve_dense_partial_pivot_lu(assembled.interpolation_matrix(), rhs)
                    .map_err(|_| BenchmarkError::InvalidCase("timed solve failed"))?;
                vector_checksum(value.weights().values())
            }
            Stage::ScalarEvaluation => {
                let scalars = if threads == 1 {
                    model.evaluate_scalars(&query_points)?
                } else {
                    evaluation_pool
                        .evaluate_scalars(Arc::clone(&model), Arc::clone(&query_points))?
                };
                vector_checksum(&scalars)
            }
            Stage::GradientEvaluation => {
                let gradients = if threads == 1 {
                    model.evaluate_gradients(&query_points)?
                } else {
                    evaluation_pool
                        .evaluate_gradients(Arc::clone(&model), Arc::clone(&query_points))?
                };
                gradient_checksum(&gradients)
            }
            Stage::EndToEnd => {
                let fitted = Arc::new(fit_single_surface_linear(
                    &case.constraints,
                    &case.parameters,
                )?);
                let scalars = if threads == 1 {
                    fitted.evaluate_scalars(&query_points)?
                } else {
                    evaluation_pool
                        .evaluate_scalars(Arc::clone(&fitted), Arc::clone(&query_points))?
                };
                let gradients = if threads == 1 {
                    fitted.evaluate_gradients(&query_points)?
                } else {
                    evaluation_pool
                        .evaluate_gradients(Arc::clone(&fitted), Arc::clone(&query_points))?
                };
                result_checksum(&scalars, &gradients)
            }
        };
        mix_u64(&mut checksum, value);
    }
    Ok(checksum)
}

fn per_iteration(duration: Duration, iterations: usize) -> u128 {
    duration.as_nanos() / iterations as u128
}

fn dataset_checksum(constraints: &Constraints, queries: &[Point]) -> u64 {
    let mut hash = FNV_OFFSET;
    mix_u64(&mut hash, constraints.interfaces.len() as u64);
    for value in &constraints.interfaces {
        mix_point(&mut hash, value.point());
        mix_u64(&mut hash, value.level().to_bits());
    }
    mix_u64(&mut hash, constraints.planars.len() as u64);
    for value in &constraints.planars {
        mix_point(&mut hash, value.point());
        for component in value.normal() {
            mix_u64(&mut hash, component.to_bits());
        }
    }
    mix_u64(&mut hash, constraints.tangents.len() as u64);
    for value in &constraints.tangents {
        mix_point(&mut hash, value.point());
        for component in value.vector() {
            mix_u64(&mut hash, component.to_bits());
        }
    }
    mix_u64(&mut hash, queries.len() as u64);
    for value in queries {
        mix_point(&mut hash, value);
    }
    hash
}

fn matrix_checksum(matrix: &DenseMatrix) -> u64 {
    let mut hash = FNV_OFFSET;
    mix_u64(&mut hash, matrix.rows() as u64);
    mix_u64(&mut hash, matrix.cols() as u64);
    for value in matrix.data() {
        mix_u64(&mut hash, value.to_bits());
    }
    hash
}

fn vector_checksum(values: &[f64]) -> u64 {
    let mut hash = FNV_OFFSET;
    mix_u64(&mut hash, values.len() as u64);
    for value in values {
        mix_u64(&mut hash, value.to_bits());
    }
    hash
}

fn gradient_checksum(values: &[[f64; 3]]) -> u64 {
    let mut hash = FNV_OFFSET;
    mix_u64(&mut hash, values.len() as u64);
    for value in values {
        for component in value {
            mix_u64(&mut hash, component.to_bits());
        }
    }
    hash
}

fn mix_point(hash: &mut u64, point: &Point) {
    mix_u64(hash, point.x().to_bits());
    mix_u64(hash, point.y().to_bits());
    mix_u64(hash, point.z().to_bits());
    mix_u64(hash, point.c().to_bits());
}

const FNV_OFFSET: u64 = 1_469_598_103_934_665_603;
const FNV_PRIME: u64 = 1_099_511_628_211;

fn mix_u64(hash: &mut u64, value: u64) {
    for byte in value.to_le_bytes() {
        *hash ^= u64::from(byte);
        *hash = hash.wrapping_mul(FNV_PRIME);
    }
}
