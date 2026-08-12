# Surfe 来源基线与迁移范围

## 冻结身份

- 上游仓库：`https://github.com/MichaelHillier/surfe.git`。
- 唯一参考提交：`290dbe0ab344f4258a4935f05cad0f153f0f69a4`。
- 本文及后续迁移记录中的“Surfe”均指该提交，不指默认分支、标签或更新后的工作树。
- T00 实测 reference：仓库内被忽略的 `.cache/surfe-reference`，处于 detached HEAD；`git rev-parse HEAD` 返回上述完整提交。
- 冻结 `License.txt` 的 SHA-256 为
  `9fd4e80cac11aa3d00278e4d2634ff0c6a169303014155f51873f9c4e6e6441f`。

reference 只用于审计、oracle 和差分验证。不得从 reference 的当前分支名、远端
HEAD 或未提交修改推断冻结身份；每次使用前都必须验证完整提交号。

## Reference 发现协议

从 GeoRBF 仓库根目录按下列唯一顺序选择第一个存在的候选：

1. 非空环境变量 `SURFE_REFERENCE_DIR` 指向的目录；
2. 同级目录 `../surfe`；
3. 被忽略目录 `.cache/surfe-reference`。

选择后必须运行：

```sh
git -C "$resolved_reference" rev-parse HEAD
git -C "$resolved_reference" status --short --branch
```

第一条的输出必须逐字等于
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`。提交不匹配、目录不是 Git
工作树或 reference 不存在时，任何依赖 reference 的 parity、oracle 或性能结论都
只能标为“未验证”，不得回退到网络默认分支，也不得伪造通过。需要 oracle 生成物
时只能使用同样被忽略的 `.cache/surfe-oracle` 或仓库外目录。

`.cache/surfe-reference` 和 `.cache/surfe-oracle` 必须保持在 `.gitignore` 中。
Surfe 源码、submodule 内容、对象、共享库、可执行文件和 oracle 生成物均不得被
GeoRBF 跟踪或打入发布包。

## T00 审阅证据

本任务完整审阅了下列冻结文件；路径后的提交适用于整项，符号用于标明与公开核心
边界的联系：

- `License.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；Surfe 自有源码的
  Crown Copyright/MIT 条款与第三方权利提示。
- `README.md@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；四类约束、建模流水线、
  Eigen/pybind11 依赖以及可选 Qt/VTK 可视化边界。
- `.gitmodules@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；`pybind11` 和
  `eigen-git-mirror` submodule 声明。
- `CMakeLists.txt@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；`math_lib`、
  `surfe_lib`、`surfepy`、OpenMP 与 `GEO_BUILDER` 的构建关系。
- `surfe_lib/surfe_api.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；
  `Surfe_API` 的约束、参数、拟合、查询和评估入口。
- `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`；
  `Surfe_API::ComputeInterpolant`、`EvaluateInterpolantAtPoint(s)`、
  `EvaluateVectorInterpolantAtPoint(s)` 及状态/异常路径。

后续任务必须在各自日志中继续使用准确的
`path@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 和相关符号，不得只写模糊的
“参考了 Surfe”。

## 纳入范围

以下是行为迁移范围，不表示复制文件，也不表示其中每个声明都是可达能力。逐文件、
逐符号覆盖和可达性分别由 T01、T02 固定：

| 来源范围 | 纳入方式 |
| --- | --- |
| `surfe_lib` 非可视化建模核心 | 迁移实际可达的约束语义、参数、核、模型、矩阵、求解、标量场和梯度场。 |
| `math_lib` 被核心调用的数值辅助 | 以纯 Rust 重写并保留调用点可观察行为。 |
| `surfe_lib/surfe_api.{h,cpp}` | 作为公开状态机、输入/输出、模型选择和错误语义的来源；最终 API 使用安全 Rust 表示。 |
| `surfe_pybindings` 到核心的适配调用 | 仅用于证明公共可达性、名称和矩阵形状；不迁移 Python/pybind11 绑定本身。 |
| `test/main.cpp` 与非可视化样例数据 | 仅作使用证据和后续 oracle 数据设计输入；不把 C++ 测试或 VTK 载入路径带入产品。 |
| `README.md`、非可视化输入/输出说明 | 仅作入口和约束语义的辅助证据；源码与 oracle 的实测行为优先。 |

公开核心边界至少包括：添加/设置 interface、planar、tangent、inequality 约束；设置
核、shape、polynomial、smoothing、anisotropy、Greedy 和 restricted range 参数；
`ComputeInterpolant`；单点/批量标量和向量评估；约束读取、bounds/resolution 与
interface 数量。具体可达分支只有在 T01/T02 有定义和调用证据后才算能力。

## 排除范围

以下内容始终排除，不因其在冻结仓库中存在而变成迁移能力：

- `geo_builder/**`、`GEO_BUILDER` 构建分支、Qt UI、Qt resource/moc/uic 代码；
- VTK 数据读取、写出、等值面提取、模型/约束显示和其他纯可视化代码；
- GUI 参数选择、文件浏览器、`CreateGRBFInterpolantFromGUIParameters`；
- `docs/**` 中的图像和 GUI/可视化教程；
- CMake 构建系统、DLL 导出宏、C++ ABI、pybind11/Python 扩展包装；
- Eigen、pybind11 submodule 内容和任何第三方 C/C++ 依赖源码；
- OpenMP 进度输出与其数据竞争等实现细节；最终批量 API 的数值输出属于范围，
  但终端进度显示不是兼容要求；
- 未定义行为、越界、未初始化读取、内存问题和数据竞争本身。

`.csv`/`.vtp` 数据文件不是产品运行时输入格式承诺。后续 fixture 可以依据审阅后的
数值数据重新表达，但不得把 VTK、Qt 或 C++ 解析器带入正常 Cargo 构建。

冻结仓库顶级内容的判定如下；“证据限定”表示可以读取来证明核心行为，但不迁移其
实现或运行时能力：

| 顶级路径 | 判定 |
| --- | --- |
| `surfe_lib/` | 纳入非可视化且实际可达的数学/建模行为；`debug.{h,cpp}` 的 Windows console 和仅诊断输出排除。 |
| `math_lib/` | 纳入被上述核心实际调用的数值行为。 |
| `surfe_pybindings/` | 证据限定；仅审计到核心的公共适配调用。 |
| `test/`、`data/` | 证据限定；不迁移 C++ harness、文件读取器或 VTK 格式运行时。 |
| `README.md`、`docs/` | 文本仅作辅助证据；图片、GUI、可视化和文件导出能力排除。 |
| `geo_builder/` | 全部排除。 |
| `eigen-git-mirror`、`pybind11` | 第三方 submodule 全部排除。 |
| `cmake/`、`CMakeLists.txt`、`setup.py` | 构建/打包证据限定，全部实现排除。 |
| `.github/`、`.travis.yml`、`appveyor.yml` | 上游 CI 证据限定，不迁移其 C++ 构建流程。 |
| `License.txt` | 合规权威来源，不作为产品算法实现。 |

为判定诊断边界，T00 还完整阅读了
`surfe_lib/debug.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的
`open_console_window`；该 Windows console 辅助没有数值输出语义，明确排除。其余
`surfe_lib` 符号的逐项纳入、部分实现、缺陷或排除分类仍由 T01/T02 依据调用证据
完成，不能由顶级目录判定替代。

## 判定规则

- 冻结源码定义、调用点和可运行 oracle 共同决定实际行为；README 或类声明不能单独
  证明能力。
- TODO、空实现、恒定返回、不可达分支和源码缺陷必须由 T02 分类，不得为对称性或
  便利自行补全。
- 后续发现只能归入 `PLAN.md` 已有任务。不得扩大本基线、改变固定顺序或新增顶级
  任务。
- 正常生产构建、测试、文档测试和发布不得要求 reference 或 oracle 存在。
