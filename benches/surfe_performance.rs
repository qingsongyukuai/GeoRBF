mod support;

use support::{
    benchmark_case, run_benchmark, BenchmarkOptions, Stage, BENCHMARK_CASE, BENCHMARK_FORMAT,
    FIXED_MULTI_THREADS,
};

use georbf::fit_single_surface_linear;

fn main() {
    if std::env::var_os("GEORBF_RUN_PERFORMANCE").is_none() {
        println!("{BENCHMARK_FORMAT} skipped=set-GEORBF_RUN_PERFORMANCE");
        return;
    }

    let arguments = std::env::args().collect::<Vec<_>>();
    let quick = arguments.iter().any(|argument| argument == "--quick");
    let stage = arguments
        .iter()
        .find_map(|argument| argument.strip_prefix("--stage="))
        .map(|name| Stage::parse(name).expect("--stage must name a benchmark stage"));
    let threads = arguments
        .iter()
        .find_map(|argument| argument.strip_prefix("--threads="))
        .map(|value| {
            value
                .parse::<usize>()
                .expect("--threads must be an integer")
        });
    let options = if quick {
        BenchmarkOptions {
            samples: 3,
            warmups: 2,
            stage,
            threads,
        }
    } else {
        BenchmarkOptions {
            stage,
            threads,
            ..BenchmarkOptions::default()
        }
    };
    let case = benchmark_case().expect("fixed GeoRBF performance case must be valid");
    let dataset_checksum = case.checksum;
    let evidence_model = fit_single_surface_linear(&case.constraints, &case.parameters)
        .expect("fixed GeoRBF performance evidence must fit");
    let evidence_indices = [0, 2_048, 4_095];
    let evidence_scalars = evidence_indices.map(|index| {
        evidence_model
            .evaluate_scalar(&case.queries[index])
            .expect("fixed scalar evidence must evaluate")
    });
    let evidence_gradients = evidence_indices.map(|index| {
        evidence_model
            .evaluate_gradient(&case.queries[index])
            .expect("fixed gradient evidence must evaluate")
    });
    let rows = run_benchmark(options).expect("fixed GeoRBF performance benchmark must succeed");
    println!(
        "{BENCHMARK_FORMAT} implementation=georbf case={BENCHMARK_CASE} fixed_multi_threads={FIXED_MULTI_THREADS} samples={} warmups={} dataset_checksum={dataset_checksum:016x}",
        options.samples, options.warmups,
    );
    println!(
        "evidence implementation=georbf scalars={} gradients={}",
        evidence_scalars
            .map(|value| format!("{value:.17e}"))
            .join(","),
        evidence_gradients
            .into_iter()
            .flatten()
            .map(|value| format!("{value:.17e}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    for row in rows {
        println!(
            "sample implementation=georbf case={BENCHMARK_CASE} threads={} stage={} index={} nanoseconds={} checksum={:016x}",
            row.threads,
            row.stage.name(),
            row.sample,
            row.nanoseconds,
            row.checksum
        );
    }
}
