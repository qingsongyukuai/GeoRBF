# T33 性能验收报告

## 结论

GeoRBF 在冻结的同机、同数据、同查询、同线程与同优化级别基准中通过全部 release-blocking 门槛。六个阶段的聚合中位数均不高于 Surfe；两个独立的三轮验收也全通过，最差 GeoRBF/Surfe 比率为 `0.932`（solve）。T32 全局 parity suite、固定 benchmark 数值见证与纯 Rust 边界均未退化。

## 冻结边界与机器

- Surfe：`.cache/surfe-reference@290dbe0ab344f4258a4935f05cad0f153f0f69a4`，工作树干净。
- Eigen gitlink：`36b95962756c1fce8e29b1f8bc45967f30773c00`。
- GeoRBF：本 T33 提交前工作树；零依赖、`#![forbid(unsafe_code)]`。
- CPU：Intel Xeon W-2125，4 physical cores / 8 logical CPUs；固定使用 physical CPU `0,2`。
- 系统：Linux `6.18.33.2-microsoft-standard-WSL2` x86_64。
- Rust：`rustc/cargo 1.85.0`；C++：`g++ 13.3.0`。
- Rust flags：bench profile（`opt-level=3`）与 `-C target-cpu=native`。
- C++ flags：`-std=c++11 -O3 -DNDEBUG -march=native -ffp-contract=off -fopenmp`。Rust 不隐式收缩普通乘加，故 C++ 同样冻结 contraction off；这也使装配矩阵能逐 bit 比较。
- 环境：`LC_ALL=C TZ=UTC OMP_DYNAMIC=FALSE`；固定 1 线程与 2 线程，2 线程分别落在两个 physical cores。

外部 C++ adapter 和 binary 只保存在 `.cache/surfe-oracle/`，不进入生产、正常 Cargo 测试或发布包。最终 adapter/binary SHA-256 分别为 `6fe9ab47b1564cf06f2a2cc0c1aaa792487e1359d7abc24075307e31e9f4dde7` 与 `2d1723f677206a81050048b61b68317fdda9ce834e911d3052d9b5e1b9f8920a`。

## 固定工作负载

`single_surface_cubic_dense_v1` 使用 Cubic ordinary Single Surface、一阶 polynomial、shape parameter 1：

- 96 Interface、16 Planar、8 Tangent；
- 4,096 个固定查询点；
- 所有坐标与 level 扰动由整数和二进制精确分母生成，避免跨语言 libm/FMA 数据漂移；
- Rust/Surfe dataset checksum 均为 `acb535ddc667960e`；
- 输入清洗后系统为 156×156 dense matrix。

每个进程预热 5 次、采样 9 次。比较器执行三轮并交替两端的先后顺序，最后对每组 27 个样本取中位数，从而平均整进程的时间漂移。为降低短样本调度噪声，每个 timed sample 内固定重复：preprocess 64、assembly 4、solve 64、scalar evaluation 32、gradient evaluation 2、end-to-end 2；报告值均除以重复数。scalar 多线程重复被一次派发给持久 worker，等价于 Surfe/OpenMP 的持久 thread team，不把每次创建 OS thread 的成本混入求值吞吐。

六阶段定义如下：

1. preprocess：owned constraint 清洗、interface grouping、layout 与 kernel setup；
2. assembly：使用已验证 layout 构造 dense interpolation matrix 与 RHS；
3. solve：输入有限性验证、partial-pivot LU、triangular solve、solution finite/residual/backward-error evidence；Surfe adapter 同样调用冻结 matrix validation 并计算对应 solution evidence；
4. scalar evaluation：全部 4,096 queries，保持输入顺序；
5. gradient evaluation：全部 4,096 queries，保持输入顺序；冻结 Surfe vector evaluator 本身为串行，共享可变 kernel 不能由 adapter 安全外部并行，GeoRBF 的不可变模型可安全分片；
6. end-to-end：preprocess、assembly、solve、全部 scalar 与 gradient evaluation。

preprocess、assembly、solve 不读取线程参数，因此比较器只使用其 canonical 1-thread median，并在 2-thread 行重复展示；只有 scalar、gradient 与 end-to-end 分别独立验收 1/2-thread 配置。原始输出仍包含全部 12 组并检查每组 checksum 稳定。

## 最终中位数

以下是第二次连续三轮聚合验收；第一次也全 PASS。

| Threads | Stage | Surfe median (ns) | GeoRBF median (ns) | Ratio | Gate |
|---:|---|---:|---:|---:|:---:|
| 1 | preprocess | 38,246 | 20,218 | 0.529 | PASS |
| 1 | assembly | 496,850 | 285,925 | 0.575 | PASS |
| 1 | solve | 434,173 | 377,692 | 0.870 | PASS |
| 1 | scalar evaluation | 5,081,176 | 1,674,073 | 0.329 | PASS |
| 1 | gradient evaluation | 29,864,897 | 7,143,238 | 0.239 | PASS |
| 1 | end-to-end | 21,300,105 | 12,097,362 | 0.568 | PASS |
| 2 | preprocess | 38,246 | 20,218 | 0.529 | PASS |
| 2 | assembly | 496,850 | 285,925 | 0.575 | PASS |
| 2 | solve | 434,173 | 377,692 | 0.870 | PASS |
| 2 | scalar evaluation | 2,769,577 | 1,343,050 | 0.485 | PASS |
| 2 | gradient evaluation | 31,332,577 | 7,011,707 | 0.224 | PASS |
| 2 | end-to-end | 23,105,373 | 12,758,752 | 0.552 | PASS |

最终两个独立的三轮验收连续全 PASS；其第一遍最差比率为 solve `0.932`，第二遍为 solve `0.870`。因此结论不依赖单轮时间漂移或单个异常快样本。

## Profile 证据与优化

机器未提供可用 `perf`，所以使用六阶段 harness 作为确定性 coarse profile，并在每个变更后只重跑相关 stage。最初 quick profile 显示主要成本为 scalar/gradient evaluation 与 end-to-end，generic assembly 重复计算 point-pair radius，naive LU 则有高频边界检查、strided lower-column pass 与逐行 pivot-tail 重载。

最终只保留由测量支持且通过 parity 回归的优化：

- 新增 `assemble_system_with_layout`，把 preprocessing 与 assembly 分离，避免重复 layout work；
- Cubic derivative vector 与 mixed Hessian 对每个 point pair 只计算一次 radius；
- Single Surface/Cubic assembly 按 value/planar/tangent blocks 写入 row-major storage，一个 pair 复用 derivative/Hessian，同时保持每个 matrix cell 的冻结算式与累加顺序；其他 model/kernel 继续走 generic path；
- Cubic scalar evaluation 在进入 query hot loop 前完成 kernel dispatch，内部直接使用冻结 Cubic 公式；
- LU 对 packed storage 直接切片，四行复用 pivot tail，并把 lower-column division 融合进同一次 row traversal；pivot 选择、每个 cell 的 update 次序和 triangular solve 次序不变；
- benchmark 的纯标准库 worker pool 只用于固定多线程测量，不引入 crate、FFI、全局可变状态或生产依赖。

曾测试的 LU 四列手工展开、8-row/2-row microkernel 与单 codegen-unit 配置均未保留：它们在隔离 stage 中退化或无稳定收益。所有失败尝试与修正写入 T33 JOURNAL。

## Parity 与确定性

- dataset checksum exact：`acb535ddc667960e`；
- assembly sample checksum 在 Rust/Surfe 间 exact，比较器不允许不一致；
- 每个 implementation 的全部 stage/sample checksum 必须稳定；
- 固定 query `0, 2048, 4095` 的 scalar 最大绝对差小于 `1.5e-13`，gradient 最大绝对差小于 `2.6e-14`，远小于 T04 `prediction_scalar`（abs `1e-9`, rel `1e-8`）和 `prediction_gradient`（abs `1e-8`, rel `1e-7`）阈值；
- 1/2-thread GeoRBF scalar、gradient 的结果顺序与 bits 相同；
- T32 正式 global-parity fixture、全部 frozen model/kernel/solver tests 和 release parity target 会在最终验证重新执行。

## 复现

在已解析并验证上述 frozen reference 后，先用报告中的 C++ flags、include path 和非可视化 Surfe sources 构建 `.cache/surfe-oracle/t33-surfe-performance`，再运行：

```text
python3 tools/compare_performance.py \
  --surfe-command 'taskset -c 0,2 env OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE LC_ALL=C TZ=UTC .cache/surfe-oracle/t33-surfe-performance' \
  --georbf-command "taskset -c 0,2 env GEORBF_RUN_PERFORMANCE=1 RUSTFLAGS='-C target-cpu=native' cargo bench --bench surfe_performance --"
```

比较器默认执行三轮并交替先运行哪一端，同时验证 header、dataset、样本数量、跨轮/每组 checksum 稳定、assembly exact、三点 prediction parity 与全部 27-sample median gate；任一失败均返回非零退出码。

本报告只证明该冻结机器/工具链/工作负载上的 release-blocking 性能不低于 Surfe，不外推到所有硬件、模型规模或第三方 BLAS 实现。T34 仍负责最终许可证、平台矩阵、发布包和全部 definition-of-done 验收。
